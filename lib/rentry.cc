#include "rentry.h"
#include "rtree.h"
#include "collection.h"
#include <string.h>
#include <stdlib.h>

// constructor/destructor
RtreeNodeEntry::RtreeNodeEntry(const int a_id, const Hypercube &a_hc) :
        m_id(a_id), m_hc(a_hc) {}

RtreeNodeEntry::~RtreeNodeEntry() {}

RtreeNodeEntry *RtreeNodeEntry::clone() const {
    return new RtreeNodeEntry(m_id, m_hc);
}

// comparison
bool RtreeNodeEntry::operator==(const RtreeNodeEntry &a_entry) const {
    return m_id == a_entry.m_id && m_hc == a_entry.m_hc;
}

bool RtreeNodeEntry::enclose(const RtreeNodeEntry &a_entry) const {
    return m_hc.enclose(a_entry.m_hc);
}

// measures
float RtreeNodeEntry::expansion(const RtreeNodeEntry &a_entry) const {
    Hypercube hc = Hypercube::combine(m_hc, a_entry.m_hc);
    return hc.volume() - m_hc.volume();
}

// update
int RtreeNodeEntry::quadraticSplit(RtreeNodeEntry **a_entry,
                                   const int a_len, const int a_split,
                                   RtreeNodeEntry **a_gp0, int &a_gp0cnt,
                                   RtreeNodeEntry **a_gp1, int &a_gp1cnt) {
    class carrier {
    public:
        const int m_dimen;
        RtreeNodeEntry *m_p;
    public:
        carrier(const int a_dimen, RtreeNodeEntry *a_p) :
                m_dimen(a_dimen), m_p(a_p) {};

        ~carrier() {};

        static int compareL(const void *a_p0, const void *a_p1) {
            carrier *c0 = *(carrier **) a_p0;
            carrier *c1 = *(carrier **) a_p1;
            float val0 = c0->m_p->m_hc.getLower()[c0->m_dimen];
            float val1 = c1->m_p->m_hc.getLower()[c1->m_dimen];
            if (val0 < val1) return -1;
            if (val0 > val1) return +1;
            return 0;
        }

        static int compareU(const void *a_p0, const void *a_p1) {
            carrier *c0 = *(carrier **) a_p0;
            carrier *c1 = *(carrier **) a_p1;
            float val0 = c0->m_p->m_hc.getUpper()[c0->m_dimen];
            float val1 = c1->m_p->m_hc.getUpper()[c1->m_dimen];
            if (val0 < val1) return -1;
            if (val0 > val1) return +1;
            return 0;
        }
    };

    int dimen = m_hc.dimen();
    float totalarea = -1;
    for (int d = 0; d < dimen; d++) {
        Collection::Array listL;
        Collection::Array listU;
        for (int i = 0; i < a_len; i++) {
            carrier *c = new carrier(d, a_entry[i]);
            listL.append(c);
            listU.append(c);
        }
        listL.sort(carrier::compareL);
        listU.sort(carrier::compareU);

        Hypercube seed0 = ((carrier *) listL[a_len - 1])->m_p->m_hc;
        Hypercube seed1 = ((carrier *) listU[0])->m_p->m_hc;

        int cnt0 = 0, cnt1 = 0;
        RtreeNodeEntry **tmp0 = new RtreeNodeEntry *[a_len];
        RtreeNodeEntry **tmp1 = new RtreeNodeEntry *[a_len];
        for (int i = 0; i < a_len; i++) {
            Hypercube test0 = Hypercube::combine(seed0, a_entry[i]->m_hc);
            Hypercube test1 = Hypercube::combine(seed1, a_entry[i]->m_hc);
            if (cnt0 > a_len - a_split)
                tmp1[cnt1++] = a_entry[i];
            else if (cnt1 > a_len - a_split)
                tmp0[cnt0++] = a_entry[i];
            else {
                if ((test0.volume() - seed0.volume()) < (test1.volume() - seed1.volume())) {
                    seed0 = test0;
                    tmp0[cnt0++] = a_entry[i];
                } else {
                    seed1 = test1;
                    tmp1[cnt1++] = a_entry[i];
                }
            }
        }

        float area = seed0.volume() + seed1.volume();
        if (area < totalarea || totalarea == -1) {
            a_gp0cnt = a_gp1cnt = 0;
            for (int i = 0; i < cnt0; i++)
                a_gp0[a_gp0cnt++] = tmp0[i];
            for (int i = 0; i < cnt1; i++)
                a_gp1[a_gp1cnt++] = tmp1[i];
            totalarea = area;
        }
        delete[] tmp0;
        delete[] tmp1;

        for (int i = 0; i < a_len; i++)
            delete (carrier *) listL[i];
        listL.clean();
        listU.clean();
    }
    return 0;
}

int RtreeNodeEntry::goodnessSplit(RtreeNodeEntry **a_entry,
                                  const int a_len, const int a_split,
                                  RtreeNodeEntry **a_gp0, int &a_gp0cnt,
                                  RtreeNodeEntry **a_gp1, int &a_gp1cnt) {
    int dimen = m_hc.dimen();
    int bestdimen = -1;
    int bestsplit = -1;
    float bestmargin = -1;
    float bestoverlap = -1;
    float bestarea = -1;
    for (int d = 0; d < dimen; d++) // ChooseSplitAxis
    {
        sortOnDimen(a_entry, a_len, d);

        float margin = 0;
        for (int i = a_split; i < a_len - a_split; i++) {
            Hypercube hc0(a_entry[0]->m_hc);
            Hypercube hc1(a_entry[i]->m_hc);
            for (int x = 1; x < i; x++)
                hc0 = Hypercube::combine(hc0, a_entry[x]->m_hc);
            for (int y = i + 1; y < a_len; y++)
                hc1 = Hypercube::combine(hc1, a_entry[y]->m_hc);
            // here, goodness is measured by margin
            margin += hc0.perimeter() + hc1.perimeter();
        }
        if (bestdimen == -1 || margin < bestmargin) {
            bestdimen = d;
            bestmargin = margin;
        }
    }

    for (int d = bestdimen; d <= bestdimen; d++)    // ChooseSplitIndex
    {
        if (bestdimen != dimen - 1)
            sortOnDimen(a_entry, a_len, bestdimen);

        for (int i = a_split; i < a_len - a_split; i++) {
            Hypercube hc0(a_entry[0]->m_hc);
            Hypercube hc1(a_entry[i]->m_hc);
            for (int x = 1; x < i; x++)
                hc0 = Hypercube::combine(hc0, a_entry[x]->m_hc);
            for (int y = i + 1; y < a_len; y++)
                hc1 = Hypercube::combine(hc1, a_entry[y]->m_hc);
            // here, goodness is measured by area and overlap
            float area = hc0.volume() + hc1.volume();
            float overlap = Hypercube::intersect(hc0, hc1).volume();    // based on overlap area
            if (bestsplit == -1 || overlap < bestoverlap || (overlap == bestoverlap && area < bestarea)) {
                bestsplit = i;
                bestoverlap = overlap;
                bestarea = area;
            }
        }
    }

    for (int d = bestdimen; d <= bestdimen; d++)    // Split
    {
        for (int i = 0; i < a_len; i++) {
            if (i < bestsplit)
                a_gp0[a_gp0cnt++] = a_entry[i];
            else
                a_gp1[a_gp1cnt++] = a_entry[i];
        }
    }

    return bestsplit;
}

int RtreeNodeEntry::packedSplit(RtreeNodeEntry **a_entry,
                                const int a_len, const int a_bulk,
                                RtreeNodeEntry **a_gp0, int &a_gp0cnt,
                                RtreeNodeEntry **a_gp1, int &a_gp1cnt) {
    int dimen = m_hc.dimen();
    float bestvol = -1;
    int bestcut = -1;
    int bestdim = -1;
    for (int d = 0; d < dimen; d++) {
        sortOnDimen(a_entry, a_len, d);
        int piece = a_len % a_bulk == 0 ? a_len / a_bulk : a_len / a_bulk + 1;
        int splitpt = piece / 2 * a_bulk;

        RtreeNodeEntry *r0;
        RtreeNodeEntry *r1;
        combine(&a_entry[0], splitpt, 0, &r0);
        combine(&a_entry[splitpt], a_len - splitpt, 1, &r1);

        float vol = r0->m_hc.volume() + r1->m_hc.volume();
        if (vol < bestvol || bestcut == -1) {
            bestdim = d;
            bestvol = vol;
            bestcut = splitpt;
        }
        delete r0;
        delete r1;
        /*
        for (int i=a_bulk; i<a_len; i+=a_bulk)
        {
            Hypercube hc0 = a_entry[0]->m_hc;
            Hypercube hc1 = a_entry[i]->m_hc;
            for (int w1=0; w1<i; w1++)
                hc0 = Hypercube::combine(hc0, a_entry[w1]->m_hc);
            for (int w2=i; w2<a_len; w2++)
                hc1 = Hypercube::combine(hc1, a_entry[w2]->m_hc);
            float vol = hc0.volume() + hc1.volume();

            float vol = r0->m_hc.volume() + r1->m_hc.volume();
            if (vol < bestvol || bestcut == -1)
            {
                bestdim = d;
                bestvol = vol;
                bestcut = i;
            }
        }
        */
    }

    if (bestdim != dimen - 1)
        sortOnDimen(a_entry, a_len, bestdim);
    a_gp0cnt = a_gp1cnt = 0;
    for (int i = 0; i < a_len; i++) {
        if (i < bestcut)
            a_gp0[a_gp0cnt++] = a_entry[i];
        else
            a_gp1[a_gp1cnt++] = a_entry[i];
    }
    return 0;
}

int RtreeNodeEntry::pickWorst(RtreeNodeEntry **a_entry,
                              const int a_len, const int a_reinsert,
                              RtreeNodeEntry **a_gp0, int &gpcnt0,
                              RtreeNodeEntry **a_gp1, int &gpcnt1) {
    // logic:
    // 1. sort all entries w.r.t. each dimension
    // 2. apply a sliding window to collect adjacent entries
    // 3. record the smallest cluster and its cut dimen and position
    // 4. put the cluster to a_gp0 and those outside the cluster to a_gp1
    //
    class carrier {
    public:
        RtreeNodeEntry *m_p;
        const int m_d;
    public:
        carrier(RtreeNodeEntry *a_p, const int a_d) :
                m_p(a_p), m_d(a_d) {};

        virtual ~carrier() {};

        static int compare(const void *a_p0, const void *a_p1) {
            carrier *c0 = *(carrier **) a_p0;
            carrier *c1 = *(carrier **) a_p1;
            if (c0->m_p->m_hc.getLower()[c0->m_d] <
                c1->m_p->m_hc.getLower()[c1->m_d])
                return -1;
            if (c0->m_p->m_hc.getLower()[c0->m_d] >
                c1->m_p->m_hc.getLower()[c1->m_d])
                return +1;
            return 0;
        };
    };

    int bestdimen = -1;
    int bestcut = -1;
    float bestvol = 0;
    int dimen = a_entry[0]->m_hc.dimen();
    for (int d = 0; d < dimen; d++) {
        Array a;
        for (int i = 0; i < a_len; i++)
            a.append(new carrier(a_entry[i], d));
        a.sort(carrier::compare);
        for (int i = 0; i < a_reinsert; i++) {
            Hypercube hc = a_entry[i]->m_hc;
            for (int w = i + 1; w < i + (a_len - a_reinsert); w++)
                hc = Hypercube::combine(hc, a_entry[i]->m_hc);

            float vol = hc.volume();
            if (vol < bestvol || bestcut == -1) {
                bestdimen = d;
                bestcut = i;
                bestvol = vol;
            }
        }
        for (int i = 0; i < a_len; i++)
            delete (carrier *) a[i];
        a.clean();
    }

    Array a;
    for (int i = 0; i < a_len; i++)
        a.append(new carrier(a_entry[i]->clone(), bestdimen));
    a.sort(carrier::compare);
    gpcnt0 = gpcnt1 = 0;
    for (int i = 0; i < a_len; i++) {
        carrier *c = (carrier *) a[i];
        if (bestcut <= i && i < bestcut + a_len - a_reinsert)
            a_gp0[gpcnt0++] = c->m_p;
        else
            a_gp1[gpcnt1++] = c->m_p;
        delete c;
    }
    a.clean();
    return 0;

    /* original R*-tree pickWorst
    class carrier
    {
    public:
        RtreeNodeEntry* m_p;
        const float     m_dist;
    public:
        carrier(RtreeNodeEntry* a_p, const float a_dist):
          m_p(a_p), m_dist(a_dist) {};
        ~carrier() {};
        static int compare(const void* a_p0, const void* a_p1)
        {
            carrier* c0 = *(carrier**)a_p0;
            carrier* c1 = *(carrier**)a_p1;
            if (c0->m_dist < c1->m_dist) return -1;
            if (c0->m_dist > c1->m_dist) return +1;
            return 0;
        };
    };

    RtreeNodeEntry* result=0;
    combine(a_entry, a_len, -1, &result);
    Point center(result->m_hc.getCenter());

    Collection::Array a;
    for (int i=0; i<a_len; i++)
        a.append(new carrier(a_entry[i]->clone(), center.distance(a_entry[i]->m_hc.getCenter())));
    a.sort(carrier::compare);

    for (int i=0; i<a_len; i++)
    {
        carrier* c = (carrier*)a[i];
        if (i < (a_len - a_reinsert))
            a_gp0[gpcnt0++] = c->m_p;
        else
            a_gp1[gpcnt1++] = c->m_p;
        delete c;
    }
    a.clean();
    delete result;
    return 0;
    */
}


int RtreeNodeEntry::combine(RtreeNodeEntry **a_entry,
                            const int a_len, const int a_id,
                            RtreeNodeEntry **a_res) {
    int dimen = a_entry[0]->m_hc.dimen();
    float cl[MAXDIMEN], cu[MAXDIMEN];
    for (int d = 0; d < dimen; d++) {
        cl[d] = a_entry[0]->m_hc.getLower()[d];
        cu[d] = a_entry[0]->m_hc.getUpper()[d];
    }
    for (int i = 0; i < a_len; i++) {
        for (int d = 0; d < dimen; d++) {
            cl[d] =
                    cl[d] < a_entry[i]->m_hc.getLower()[d] ?
                    cl[d] : a_entry[i]->m_hc.getLower()[d];
            cu[d] =
                    cu[d] > a_entry[i]->m_hc.getUpper()[d] ?
                    cu[d] : a_entry[i]->m_hc.getUpper()[d];
        }
    }
    Hypercube hc(dimen, cl, cu);
    (*a_res) = new RtreeNodeEntry(a_id, hc);
    return a_len;

    /*
    Hypercube hc = a_entry[0]->m_hc;
    for (int i=1; i<a_len; i++)
        hc = Hypercube::combine(hc, a_entry[i]->m_hc);
    (*a_res) = new RtreeNodeEntry(a_id,hc);
    return a_len;
    */
}

int dimenToSort = 0;

int compareDimen(const void *a0, const void *a1) {
    RtreeNodeEntry *r0 = *(RtreeNodeEntry **) a0;
    RtreeNodeEntry *r1 = *(RtreeNodeEntry **) a1;
    float v0 = r0->m_hc.getLower()[dimenToSort];
    float v1 = r1->m_hc.getLower()[dimenToSort];
    if (v0 < v1) return -1;
    if (v0 > v1) return +1;
    return 0;
}

int RtreeNodeEntry::sortOnDimen(RtreeNodeEntry **a_entry,
                                const int a_len, const int a_dimen) {
    dimenToSort = a_dimen;
    qsort(a_entry, a_len, sizeof(RtreeNodeEntry *), compareDimen);
    return a_len;
}

// memory operations
int RtreeNodeEntry::toMem(char *a_content, int &a_len, const bool a_pt) const {
    int init = a_len;
    memcpy(&a_content[a_len], &m_id, sizeof(m_id));
    a_len += sizeof(m_id);

    for (int i = 0; i < m_hc.dimen(); i++) {
        float l = m_hc.getLower()[i];
        float u = m_hc.getUpper()[i];
        memcpy(&a_content[a_len], &l, sizeof(l));
        a_len += sizeof(l);
        if (!a_pt) {
            memcpy(&a_content[a_len], &u, sizeof(u));
            a_len += sizeof(u);
        }
    }

    return a_len - init;
}

RtreeNodeEntry *RtreeNodeEntry::fromMem(const char *a_content, int &a_len,
                                        const int dimen, const bool a_pt) {
    int id;
    memcpy(&id, &a_content[a_len], sizeof(id));
    a_len += sizeof(id);

    float cl[20];
    float cu[20];
    for (int i = 0; i < dimen; i++) {
        memcpy(&cl[i], &a_content[a_len], sizeof(cl[i]));
        a_len += sizeof(cl[i]);
        cu[i] = cl[i];
        if (!a_pt) {
            memcpy(&cu[i], &a_content[a_len], sizeof(cu[i]));
            a_len += sizeof(cu[i]);
        }
    }

    Hypercube hc(dimen, cl, cu);
    return new RtreeNodeEntry(id, hc);
}

// info
int RtreeNodeEntry::size(const int a_dimen, bool isPoint) {
    if (isPoint)
        return sizeof(int) + Hypercube::size(a_dimen) / 2;
    return sizeof(int) + Hypercube::size(a_dimen);
}
