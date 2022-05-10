#include "mainmem.h"
#include "rnode.h"
#include <inttypes.h>
MainMemory::MainMemory(const int a_pagesize):
Memory(a_pagesize),
m_buffer(1),                   //
m_largest(0)
{
}

MainMemory::~MainMemory()
{
    flush();
}

RtreeNode* MainMemory::loadPage(const int a_pageID)
{
    RtreeNode* n = (RtreeNode*)m_buffer.get(a_pageID);
    return n->clone();
}

int MainMemory::allocatePage()
{
    if (m_avail.isEmpty())
    {
        for (int i=m_largest; i<m_largest+10; i++)
            m_avail.enqueue((void*)i);
        m_largest+=10;
    }
    return (intptr_t)m_avail.dequeue();
}

void MainMemory::writePage(const int a_pageID, const RtreeNode* m_p)
{
    removePage(a_pageID);
    m_buffer.put(a_pageID,m_p->clone());    // the copy of RtreeNode is kept!
}

void MainMemory::removePage(const int a_pageID)
{
    RtreeNode* n = (RtreeNode*)m_buffer.get(a_pageID);
    if (n != 0)
    {
        m_buffer.remove(a_pageID);
        delete n;
    }
}

void MainMemory::flush()
{
    for (HashReader rdr(m_buffer); !rdr.isEnd(); rdr.next())
    {
        RtreeNode* n = (RtreeNode*)rdr.getVal();
        delete n;
    }
    m_buffer.clean();
}
