#ifndef CACHEMEM_DEFINED
#define CACHEMEM_DEFINED

#include "mem.h"
#include "collection.h"

class CacheMemory: public Memory
{
// data members
public:
    Memory& m_mem;          // underlying memory
    Hash    m_buffer;       // buffer of cache Rtree nodes
    int     m_cachesize;    // cache size in terms of # of nodes
    int     m_accesstick;   // logical clock of cache
    Queue   m_lruq;         // a queue to keep track of access (LRU based)
// methods
public:
    CacheMemory(Memory& a_mem, const int a_cachesize);
    virtual ~CacheMemory();
    // search
    virtual RtreeNode* loadPage(    // read a node (page) from cache
        const int a_pageID);
    // update
    virtual int allocatePage();     // obtain ID of an unused page
    virtual void writePage(         // write through a page to memory & cache
        const int a_pageID, const RtreeNode* p);
    virtual void removePage(        // remove a page from cache & memory
        const int a_pageID);
    virtual void flush();
    virtual void adjust();
};

#endif

