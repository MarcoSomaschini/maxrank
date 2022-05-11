#include "rnode.h"
#include "rtree.h"
#include "rentry.h"
#include <string.h>
#include <math.h>

RtreeNode::RtreeNode(const Rtree *a_rtree,
                     const int a_pageid,
                     const int a_level, const int a_parent) :
        m_rtree(a_rtree),
        m_pageid(a_pageid),
        m_level(a_level),
        m_parent(a_parent),
        m_usedSpace(0) {
    int max = (m_level == 0 ? m_rtree->m_maxLeafChild : m_rtree->m_maxNodeChild) + 1;
    m_entry = new RtreeNodeEntry *[max];
    memset(m_entry, 0, sizeof(RtreeNodeEntry *) * max);
}

RtreeNode::RtreeNode(const int a_pageid,       //my code 
                     const int a_level, const int a_parent) :
        m_pageid(a_pageid),
        m_level(a_level),
        m_parent(a_parent),
        m_usedSpace(0) {
    int max = 2;
    m_entry = new RtreeNodeEntry *[max];
    memset(m_entry, 0, sizeof(RtreeNodeEntry *) * max);
}

RtreeNode::~RtreeNode() {
    for (int i = 0; i < m_usedSpace; i++)
        delete m_entry[i];
    delete[] m_entry;
}

RtreeNode *RtreeNode::clone() const {
    RtreeNode *node = new RtreeNode(m_rtree, m_pageid, m_level, m_parent);
    for (int i = 0; i < m_usedSpace; i++)
        node->m_entry[node->m_usedSpace++] = m_entry[i]->clone();
    return node;
}

// search
bool RtreeNode::isRoot() const {
    return m_parent == -1;
}

bool RtreeNode::isLeaf() const {
    return m_level == 0;
}

RtreeNodeEntry *RtreeNode::genNodeEntry() const {
    RtreeNodeEntry *newentry = 0;
    if (m_usedSpace > 0) {
        m_entry[0]->combine(m_entry, m_usedSpace, m_pageid, &newentry);
    }
    return newentry;
}

// update
int RtreeNode::insert(const RtreeNodeEntry &a_entry) {
    int status = NODE_UNCHANGED;
    m_entry[m_usedSpace++] = a_entry.clone();
    int max = (m_level == 0 ? m_rtree->m_maxLeafChild : m_rtree->m_maxNodeChild) + 1;
    if (m_usedSpace == max) {
        status = NODE_OVERFLOW;
    } else {
        RtreeNodeEntry *oldentry;
        m_entry[0]->combine(m_entry, m_usedSpace - 1, m_pageid, &oldentry);
        RtreeNodeEntry *newentry;
        m_entry[0]->combine(m_entry, m_usedSpace, m_pageid, &newentry);
        status = !(*oldentry == *newentry) ? NODE_CHANGED : status;
        delete oldentry;
        delete newentry;
    }
    return status;
}

/*
int RtreeNode::insert(const RtreeNodeEntry& a_entry, int flag)    //my code
{
    m_entry[0] = a_entry.clone();   
    return 0;
}
//*/

void RtreeNode::quadraticSplit(RtreeNode **a_p0, RtreeNode **a_p1) {
    int split = m_level > 0 ? m_rtree->m_minNodeChild : m_rtree->m_minLeafChild;
    (*a_p0) = new RtreeNode(m_rtree, -1, m_level, m_parent);
    (*a_p1) = new RtreeNode(m_rtree, -1, m_level, m_parent);
    m_entry[0]->quadraticSplit(
            m_entry, m_usedSpace, split,
            (*a_p0)->m_entry, (*a_p0)->m_usedSpace,
            (*a_p1)->m_entry, (*a_p1)->m_usedSpace);
    return;
}

void RtreeNode::goodnessSplit(RtreeNode **a_p0, RtreeNode **a_p1) {
    int split = m_level > 0 ? m_rtree->m_minNodeChild : m_rtree->m_minLeafChild;
    (*a_p0) = new RtreeNode(m_rtree, -1, m_level, m_parent);
    (*a_p1) = new RtreeNode(m_rtree, -1, m_level, m_parent);
    m_entry[0]->goodnessSplit(
            m_entry, m_usedSpace, split,
            (*a_p0)->m_entry, (*a_p0)->m_usedSpace,
            (*a_p1)->m_entry, (*a_p1)->m_usedSpace);
    return;
}

void RtreeNode::pickWorst(const int a_reinsert,
                          RtreeNode **a_remainNode, RtreeNodeEntry **a_residue) {
    (*a_remainNode) = new RtreeNode(m_rtree, m_pageid, m_level, m_parent);
    int cnt = 0;
    m_entry[0]->pickWorst(
            m_entry, m_usedSpace, a_reinsert,
            (*a_remainNode)->m_entry, (*a_remainNode)->m_usedSpace,
            a_residue, cnt);
    return;
}


int RtreeNode::replace(const RtreeNodeEntry &a_entry) {
    int status = NODE_UNCHANGED;
    for (int i = 0; i < m_usedSpace; i++) {
        if (m_entry[i]->m_id == a_entry.m_id) {
            if (!(*m_entry[i] == a_entry))
                status = NODE_CHANGED;
            delete m_entry[i];
            m_entry[i] = a_entry.clone();
            break;
        }
    }
    return status;
}

int RtreeNode::remove(const RtreeNodeEntry &a_entry) {
    // find the target entry
    int status = NODE_UNCHANGED;
    int found = -1;
    for (int i = 0; i < m_usedSpace && found == -1; i++)
        if (m_entry[i]->m_id == a_entry.m_id)
            found = i;

    // remove deleted entry
    if (found != -1) {
        int min = (m_level == 0 ? m_rtree->m_minLeafChild : m_rtree->m_minNodeChild) + 1;
        if (m_usedSpace <= min) {
            m_entry[found] = m_entry[--m_usedSpace];
            status = NODE_UNDERFULL;
        } else {
            RtreeNodeEntry *oldentry;
            m_entry[0]->combine(m_entry, m_usedSpace, m_pageid, &oldentry);
            m_entry[found] = m_entry[--m_usedSpace];

            RtreeNodeEntry *newentry;
            m_entry[0]->combine(m_entry, m_usedSpace, m_pageid, &newentry);

            status = !(*oldentry == *newentry) ? NODE_CHANGED : status;
            delete oldentry;
            delete newentry;
        }
    }
    return status;
}

// memory operations
int RtreeNode::toMem(char *a_content, int &a_len) const {
    memcpy(&a_content[a_len], &m_pageid, sizeof(m_pageid));
    a_len += sizeof(m_pageid);
    memcpy(&a_content[a_len], &m_level, sizeof(m_level));
    a_len += sizeof(m_level);
    memcpy(&a_content[a_len], &m_parent, sizeof(m_parent));
    a_len += sizeof(m_parent);
    memcpy(&a_content[a_len], &m_usedSpace, sizeof(m_usedSpace));
    a_len += sizeof(m_usedSpace);
    for (int i = 0; i < m_usedSpace; i++) {
        int elen = 0;
        m_entry[i]->toMem(&a_content[a_len], elen, m_level == 0 && m_rtree->m_pointOnly);
        a_len += elen;
    }
    return a_len;
}

RtreeNode *RtreeNode::fromMem(const char *a_p, int &a_len, const Rtree *a_rtree,
                              RtreeNodeEntry *(*fromMem)(
                                      const char *a_mem, int &a_len,
                                      const int dimen, const bool a_pt)) {
    int pageid;
    int level;
    int parent;
    int usedSpace;
    int dimen = a_rtree->m_dimen;
    memcpy(&pageid, &a_p[a_len], sizeof(pageid));
    a_len += sizeof(pageid);
    memcpy(&level, &a_p[a_len], sizeof(level));
    a_len += sizeof(level);
    memcpy(&parent, &a_p[a_len], sizeof(parent));
    a_len += sizeof(parent);
    memcpy(&usedSpace, &a_p[a_len], sizeof(usedSpace));
    a_len += sizeof(usedSpace);

    RtreeNode *node = new RtreeNode(a_rtree, pageid, level, parent);

    for (int i = 0; i < usedSpace; i++) {
        RtreeNodeEntry *e = (RtreeNodeEntry *) fromMem(a_p, a_len, dimen, level == 0 && a_rtree->m_pointOnly);
        node->m_entry[node->m_usedSpace++] = e;
    }
    return node;
}

// info
int RtreeNode::size() {
    return 4 * sizeof(int);
}
