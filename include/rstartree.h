/* ----------------------------------------------------------------------------
    This header file declares R*tree class derived from Rtree.
---------------------------------------------------------------------------- */

#ifndef RSTARTREE_DEFINED
#define RSTARTREE_DEFINED

#include "rtree.h"

class Rstartree : public Rtree {
// data members
public:
    const int m_reinsert;      // number of reinserted entry
// methods
public:
    // constructor/destructor
    Rstartree(
            Memory &a_memory,
            const int a_dimen,
            const int a_maxNodeFill, const int a_maxLeafFill,
            const int a_minNodeFill, const int a_minLeafFill,
            const int a_reinsert,
            const bool a_pointOnly = false);

    virtual ~Rstartree();

    // update
    int insert(
            const RtreeNodeEntry &a_entry,          // add an entry to Rtree at
            const int a_level = 0);                   // a specified level (def: 0)
};

#endif // RSTARTREE_DEFINED


