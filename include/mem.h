/* ----------------------------------------------------------------------------
    This header file includes abstract class Memory declaration.
    It is only an interface to store R-tree data.

    The user needs to derive a class based on this class for real storage.
---------------------------------------------------------------------------- */

#ifndef MEM_DEFINED
#define MEM_DEFINED

#include "collection.h"

class Rtree;

class RtreeNode;

class Memory {
// data members
public:
    Rtree *m_rtree;        // rtree
    int m_rootPageID;   // root page
    const int m_pagesize;     // page size
// methods
public:
    // constructor/destructor
    Memory(const int a_pagesize) : m_rootPageID(-1), m_pagesize(a_pagesize) {};

    virtual ~Memory() {};

    // search
    virtual RtreeNode *loadPage(    // read a node (page) from memory
            const int a_pageID) = 0;

    // update
    virtual int allocatePage() = 0;   // obtain the ID of an unused page
    virtual void writePage(         // write a page to mem
            const int a_pageID, const RtreeNode *m_p) = 0;

    virtual void removePage(        // remove a page from mem
            const int a_pageID) = 0;

    virtual void flush() = 0;         // clean the memory
};

#endif // MEM_DEFINED
