#include "arentry.h"
#include <string>

// constructor/destructor
ARtreeNodeEntry::ARtreeNodeEntry(const int a_id, const Hypercube& a_hc,
                                 const int a_count):
RtreeNodeEntry(a_id,a_hc),
m_count(a_count)
{}

ARtreeNodeEntry::~ARtreeNodeEntry()
{}

RtreeNodeEntry* ARtreeNodeEntry::clone() const
{
    return new ARtreeNodeEntry(m_id,m_hc,m_count);
}

// comparisons
bool ARtreeNodeEntry::operator==(const RtreeNodeEntry& a_entry) const
{
    ARtreeNodeEntry* a = (ARtreeNodeEntry*)&a_entry;
    return m_count == a->m_count && m_id == a->m_id && m_hc == a->m_hc;
}

// update
int ARtreeNodeEntry::combine(RtreeNodeEntry** a_entry, const int a_len,
                             const int a_id,
                             RtreeNodeEntry** a_res)
{
    RtreeNodeEntry* tmp = 0;
    RtreeNodeEntry::combine(a_entry, a_len, a_id, &tmp);
    int cnt=0;
    for (int i=0; i<a_len; i++)
        cnt += ((ARtreeNodeEntry*)a_entry[i])->m_count;
    (*a_res) = new ARtreeNodeEntry(tmp->m_id, tmp->m_hc, cnt);
    delete tmp;
    return a_len;
}

// static operations
int ARtreeNodeEntry::toMem(char* a_content, int& a_len,
                           const bool a_pt) const
{
    int init=a_len;
    RtreeNodeEntry::toMem(a_content, a_len, a_pt);
    memcpy(&(a_content[a_len]), &m_count, sizeof(int)); a_len += sizeof(int);
    return a_len -init;
}

RtreeNodeEntry* ARtreeNodeEntry::fromMem(const char* a_p, int& a_len,
                                         const int dimen,
                                         const bool a_pt)
{
    RtreeNodeEntry* tmp = RtreeNodeEntry::fromMem(a_p, a_len, dimen, a_pt);
    int cnt;
    memcpy(&cnt, &a_p[a_len], sizeof(int));  a_len += sizeof(int);
    ARtreeNodeEntry* res = new ARtreeNodeEntry(tmp->m_id, tmp->m_hc, cnt);
    delete tmp;
    return res;
}

// info
int ARtreeNodeEntry::size(const int a_dimen, bool isPoint)
{
    return RtreeNodeEntry::size(a_dimen,isPoint) + sizeof(int);
}
