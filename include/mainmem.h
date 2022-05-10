/* ----------------------------------------------------------------------------
    This header file includes class MainMemory declaration.
    It stores R-tree data in a main memory.
---------------------------------------------------------------------------- */

#ifndef MAINMEMORY_DEFINED
#define MAINMEMORY_DEFINED

#include "mem.h"
#include "collection.h"

class MainMemory: public Memory
{
// data members
public:
    Hash    m_buffer;       // content of pages
    Queue   m_avail;
    int     m_largest;      // largest page ID available to allocate so far.
// methods
public:
    // constructor/destructor
    MainMemory(const int a_pagesize);
    virtual ~MainMemory();
    // search
    virtual RtreeNode* loadPage(    // read a node (page) from main memory
        const int a_pageID);
    // update
    virtual int allocatePage();     // obtain ID of an unused page
    virtual void writePage(         // write a page to main memory
        const int a_pageID, const RtreeNode* m_p);
    virtual void removePage(        // remove a page from main memory
        const int a_pageID);
    virtual void flush();           // clean the memory
};

#endif  // MAINMEMORY_DEFINED


