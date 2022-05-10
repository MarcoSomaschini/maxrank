#include "search.h"
#include "rtree.h"
#include "rentry.h"
#include "rnode.h"
#include "mem.h"
#include "hypercube.h"
#include "psdraw.h"
#include <iostream>
#include <map>
#include <string.h>
#include "stdio.h"

#include <inttypes.h>

#include "virtualRnode.h"

#define MAXPAGEID 9999999
#define ENCLOSEDBY    0      //enclosed by some MBR
#define NOENCLOSURE   1      //enclose some MBR

#define COMPUTE_QUADRANTS
 
using namespace std;

Search::Search()
{
}

Search::~Search()
{}

int Search::window(Rtree& a_rtree, const Point& a_pt,
				   const float* a_l, 
				   Array& a_result, int& a_maxStackSize, Array& a_pageaccessed)
{
    a_maxStackSize = 0;
    Stack s;
    s.push((void*)a_rtree.m_memory.m_rootPageID);
    while (!s.isEmpty())
    {
        a_maxStackSize = s.size() > a_maxStackSize ? s.size() : a_maxStackSize;
        int pageid = (intptr_t)s.pop();
        RtreeNode* node = a_rtree.m_memory.loadPage(pageid);
        
        //my code
        VirtualRNode* n=new VirtualRNode;
        n->copyData(*node);        
        n->copyEntries(*node,node->m_usedSpace);
        
        a_pageaccessed.append((void*)pageid);
        if (n->isLeaf())
        {
            for (int i=0; i<n->m_usedSpace; i++)
            {
                bool covered = true;
                for (int d=0; d<a_rtree.m_dimen && covered; d++)
                    if (n->m_entry[i]->m_hc.mindist(a_pt,d) > a_l[d])
                        covered = false;
                if (covered)
                    a_result.append(n->m_entry[i]->clone());
            }
            /*
            for (int i=0; i<n->m_usedSpace; i++)
				if (n->m_entry[i]->m_hc.mindist(a_pt,0) <= a_l &&
				    n->m_entry[i]->m_hc.mindist(a_pt,1) <= a_h)
                    a_result.append(n->m_entry[i]->clone());
            */
        }
        else
        {
            for (int i=0; i<n->m_usedSpace; i++)
            {
                bool covered = true;
                for (int d=0; d<a_rtree.m_dimen && covered; d++)
                    if (n->m_entry[i]->m_hc.mindist(a_pt,d) > a_l[d])
                        covered = false;
                if (covered)
                    s.push((void*)n->m_entry[i]->m_id);
            }
            /*
            for (int i=0; i<n->m_usedSpace; i++)
				if (n->m_entry[i]->m_hc.mindist(a_pt,0) <= a_l &&
					n->m_entry[i]->m_hc.mindist(a_pt,1) <= a_h)
                    s.push((void*)n->m_entry[i]->m_id);
            */
        }
        delete n;
    }
    return a_result.size();
}


int Search::range(Rtree& a_rtree, const Point& a_pt, const float a_range,
                  Array& a_result, int& maxStackSize, Array& a_pageaccessed)
{
    maxStackSize = 0;
    Stack s;
    s.push((void*)a_rtree.m_memory.m_rootPageID);
    while (!s.isEmpty())
    {
        maxStackSize = s.size() > maxStackSize ? s.size() : maxStackSize;
        int pageid = (intptr_t)s.pop();
        RtreeNode* node = a_rtree.m_memory.loadPage(pageid);

        //my code
        VirtualRNode* n=new VirtualRNode;
        n->copyData(*node);        
        n->copyEntries(*node,node->m_usedSpace);        

        a_pageaccessed.append((void*)pageid);
        if (n->isLeaf())
        {
            for (int i=0; i<n->m_usedSpace; i++)
                if (n->m_entry[i]->m_hc.mindist(a_pt) <= a_range)
                    a_result.append(n->m_entry[i]->clone());
        }
        else
        {
            for (int i=0; i<n->m_usedSpace; i++)
                if (n->m_entry[i]->m_hc.mindist(a_pt) <= a_range)
                    s.push((void*)n->m_entry[i]->m_id);
        }
        delete n;
    }
    return a_result.size();
}

float Search::nn(Rtree& a_rtree, const Point& a_pt, const int k,
                 Array& a_result, int& maxHeapSize, Array& a_pageaccessed)
{
    class carrier
    {
    public:
        const RtreeNodeEntry    *m_entry;
        const float             m_dist;
        const bool              m_isObj;
    public:
        carrier(const RtreeNodeEntry& e, const float a_dist, const bool a_o):
          m_entry(e.clone()), m_dist(a_dist), m_isObj(a_o) {};
        ~carrier() {};
        static int compare(const void* a_p0, const void* a_p1)
        {
            carrier* c0 = *(carrier**)a_p0;
            carrier* c1 = *(carrier**)a_p1;
            if (c0->m_dist < c1->m_dist) return -1;
            if (c0->m_dist > c1->m_dist) return +1;
            if (!c0->m_isObj && c1->m_isObj) return -1;
            if (c0->m_isObj && !c1->m_isObj) return +1;
            return 0;
        };
    };

    maxHeapSize = 0;
    float dist = -1;
    RtreeNodeEntry e(a_rtree.m_memory.m_rootPageID, Hypercube(a_pt.m_dimen,0,0));
    BinHeap h(carrier::compare);
    h.insert(new carrier(e,0,false));

    while (!h.isEmpty())
    {
        maxHeapSize = h.size() > maxHeapSize ? h.size() : maxHeapSize;
        carrier* c = (carrier*)h.removeTop();
        if (a_result.size() < k)
        {
            if (c->m_isObj)
            {
                a_result.append(c->m_entry->clone());
                if (a_result.size() == 0)
                    dist = c->m_dist;
            }
            else
            {
                RtreeNode* n = a_rtree.m_memory.loadPage(c->m_entry->m_id);
                a_pageaccessed.append((void*)c->m_entry->m_id);
                for (int i=0; i<n->m_usedSpace; i++)
                    h.insert(new carrier(
                        *n->m_entry[i],
                        n->m_entry[i]->m_hc.mindist(a_pt),
                        n->m_level == 0));   // is it object?
                delete n;
            }
        }
        delete c->m_entry;
        delete c;
    }
    return dist;
}

void Search::dump(const Rtree& a_rtree,
                  const int a_bottomlevel, const int a_toplevel,
                  const char* epsname)
{
    RtreeNode* root = a_rtree.m_memory.loadPage(a_rtree.m_memory.m_rootPageID);
    RtreeNodeEntry* entry = root->genNodeEntry();
    Hypercube& hc = entry->m_hc;
    PSDraw psdraw(epsname, hc.getLower()[0], hc.getLower()[1], hc.getUpper()[0], hc.getUpper()[1]);
    const int top = a_toplevel < root->m_level ? a_toplevel : root->m_level;
    delete entry;
    delete root;

    Stack s;
    s.push((void*)a_rtree.m_memory.m_rootPageID);
    while (!s.isEmpty())
    {
        int pageid = (intptr_t)s.pop();
        RtreeNode* n = a_rtree.m_memory.loadPage(pageid);
        if (a_bottomlevel <= n->m_level && n->m_level <= a_toplevel)
        {
            for (int i=0; i<n->m_usedSpace; i++)
                psdraw.box(n->m_entry[i]->m_hc, 0, 1, 1-(n->m_level * 1.0 / top));
        }

        if (n->m_level > a_bottomlevel)
        {
            for (int i=0; i<n->m_usedSpace; i++)
                s.push((void*)n->m_entry[i]->m_id);
        }
        delete n;
    }
}
