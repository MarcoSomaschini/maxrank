/* ----------------------------------------------------------------------------
    This header file includes class RtreeNodeEntry declaration.
    It is a component in RtreeNode and it is a holder of IDs of child/objects
    and corresponding MBBs.

    The user may need to derive a class based on this class for r-tree
    variants.
---------------------------------------------------------------------------- */

#ifndef RTREENODEENTRY_DEFINED
#define RTREENODEENTRY_DEFINED

#include "hypercube.h"

class RtreeNodeEntry {
// data members
public:
    const int m_id;
    Hypercube m_hc;
// methods
public:
    // constructor/destructor
    RtreeNodeEntry(const int a_id, const Hypercube &a_hc);

    virtual ~RtreeNodeEntry();

    virtual RtreeNodeEntry *clone() const;

    //
    // comparison
    virtual bool operator==(const RtreeNodeEntry &a_entry) const;

    virtual bool enclose(const RtreeNodeEntry &a_entry) const;

    //
    // measures
    virtual float expansion(const RtreeNodeEntry &a_entry) const;

    //
    // update
    virtual int quadraticSplit(     // quadratic cost split SIGMOD84
            RtreeNodeEntry **a_entry,
            const int a_len, const int a_split,
            RtreeNodeEntry **a_gp0, int &a_gp0cnt,
            RtreeNodeEntry **a_gp1, int &a_gp1cnt);

    virtual int goodnessSplit(      // R*tree split algo VLDB90
            RtreeNodeEntry **a_entry,
            const int a_len, const int a_split,
            RtreeNodeEntry **a_gp0, int &a_gp0cnt,
            RtreeNodeEntry **a_gp1, int &a_gp1cnt);

    virtual int packedSplit(        // split algo for TGS
            RtreeNodeEntry **a_entry,
            const int a_len, const int a_bulk,
            RtreeNodeEntry **a_gp0, int &a_gp0cnt,
            RtreeNodeEntry **a_gp1, int &a_gp1cnt);

    virtual int pickWorst(          // pick out entries for R*tree reinsertion
            RtreeNodeEntry **a_entry,
            const int a_len, const int a_reinsert,
            RtreeNodeEntry **a_gp0, int &gpcnt0,
            RtreeNodeEntry **g_gp1, int &gpcnt1);

    virtual int combine(            // merge a number of entries into one
            RtreeNodeEntry **a_entry,
            const int a_len, const int a_id,
            RtreeNodeEntry **a_res);

    //
    // sort
    virtual int sortOnDimen(        // sort a number of entries along dimen
            RtreeNodeEntry **a_entry,
            const int a_len, const int a_dimen);

    //
    // memory operations
    virtual int toMem(
            char *a_content, int &a_len,
            const bool a_pt) const;

    static RtreeNodeEntry *fromMem(     // to convert a bytestring to a
            const char *a_p, int &a_len,    // RtreeNodeEntry
            const int a_dimen,
            const bool a_pt);

    //
    // info
    static int size(const int a_dimen, bool isPoint = false);   // storage size
};

#endif // RTREENODEENTRY_DEFINED
