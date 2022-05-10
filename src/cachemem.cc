#include "cachemem.h"
#include "rnode.h"
#include "rentry.h"

class CacheBlock
{
public:
    int         m_lastaccess;
    RtreeNode*  m_node;
};

class QueueEntry
{
public:
    int         m_pageid;
    int         m_access;
public:
    QueueEntry(int a_pageid, int a_access):
      m_pageid(a_pageid), m_access(a_access) {};
};

CacheMemory::CacheMemory(Memory& a_mem, const int a_cachesize):
Memory(a_mem.m_pagesize),
m_mem(a_mem),
m_cachesize(a_cachesize)
{}

CacheMemory::~CacheMemory()
{
    flush();
}

RtreeNode* CacheMemory::loadPage(const int a_pageID)
{
    CacheBlock* cb = (CacheBlock*)m_buffer.get(a_pageID);
    if (cb == 0)
    {
        // load a page into cache;
        RtreeNode* n = m_mem.loadPage(a_pageID);
        cb = new CacheBlock;
        cb->m_node = n;
        m_buffer.put(a_pageID, cb);
    }
    // update access q
    cb->m_lastaccess = m_accesstick++;
    m_lruq.enqueue(new QueueEntry(a_pageID,cb->m_lastaccess));
    if (m_lruq.length() > m_cachesize) adjust();
    return cb->m_node;
}

int CacheMemory::allocatePage()
{
    return m_mem.allocatePage();
}

void CacheMemory::writePage(const int a_pageID, const RtreeNode* a_node)
{
    m_mem.writePage(a_pageID, a_node);
    CacheBlock* cb = (CacheBlock*)m_buffer.get(a_pageID);
    if (cb != 0)
    {
        delete cb->m_node;
        cb->m_node = a_node->clone();
    }
    else
    {
        cb = new CacheBlock;
        cb->m_node = a_node->clone();
        m_buffer.put(a_pageID,cb);
    }
    cb->m_lastaccess = m_accesstick++;
    m_lruq.enqueue(new QueueEntry(a_pageID,cb->m_lastaccess));
    if (m_lruq.length() > m_cachesize) adjust();
    return;
}

void CacheMemory::removePage(const int a_pageID)
{
    m_mem.removePage(a_pageID);
    m_buffer.remove(a_pageID);
}

void CacheMemory::flush()
{
    for (HashReader rdr(m_buffer); !rdr.isEnd(); rdr.next())
    {
        CacheBlock* cb = (CacheBlock*)rdr.getVal();
        delete cb->m_node;
        delete cb;
    }
    m_buffer.clean();
    while (!m_lruq.isEmpty())
    {
        QueueEntry* e = (QueueEntry*)m_lruq.dequeue();
        delete e;
    }
}

void CacheMemory::adjust()
{
    while (m_lruq.length() > m_cachesize)
    {
        QueueEntry* e = (QueueEntry*)m_lruq.dequeue();
        CacheBlock* cb = (CacheBlock*)m_buffer.get(e->m_pageid);
        if (cb != 0)
        {
            if (cb->m_lastaccess == e->m_access)
            {
                m_buffer.remove(e->m_pageid);
                delete cb->m_node;
                delete cb;
            }
        }
        delete e;
    }
}
