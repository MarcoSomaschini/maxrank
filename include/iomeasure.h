/* ----------------------------------------------------------------------------
    This header file includes declaration of class IOCount that
    simulate and measure some cached I/Os
---------------------------------------------------------------------------- */
#ifndef IOMEASURE_DEFINED
#define IOMEASURE_DEFINED

#include "collection.h"

class IOMeasure {
public:
    static int lru(                     // measure the I/O based on initially 
            const Array &a_pageaccessed,    // empty LRU cache
            const int a_cachesize);
};

#endif
