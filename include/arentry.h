#ifndef ARTREENODEENTRY_DEFINED
#define ARTREENODEENTRY_DEFINED

#include "rentry.h"

class ARtreeNodeEntry: public RtreeNodeEntry
{
// data members
public:
    int m_count;    // count of objects associated with this entry
// methods
public:
    // constructor/destructor
    ARtreeNodeEntry(const int a_id, const Hypercube& a_hc, const int a_count);
    virtual ~ARtreeNodeEntry();
    virtual RtreeNodeEntry* clone() const;
    // comparisons
    virtual bool operator==(const RtreeNodeEntry& a_entry) const;
    // update
    int combine(                    // merge a number of entries into one
        RtreeNodeEntry** a_entry,
        const int a_len, const int a_id,
        RtreeNodeEntry** a_res);
    // static operations
    int toMem(
        char* a_content, int& a_len,
        const bool a_pt) const;
    static RtreeNodeEntry* fromMem(     // to convert a bytestring to a 
        const char* a_p, int& a_len,    // RtreeNodeEntry
        const int dimen,
        const bool a_pt);
    // info
    static int size(const int a_dimen, bool isPoint=false);
};

#endif // ARTREENODEENTRY_DEFINED

