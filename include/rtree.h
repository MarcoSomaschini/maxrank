/* ----------------------------------------------------------------------------
    This header file declares Rtree class.

    The user may need to derive from this class for Rtree variants
---------------------------------------------------------------------------- */

#ifndef RTREE_DEFINED
#define RTREE_DEFINED

class Memory;
class RtreeNode;
class RtreeNodeEntry;
#include "collection.h"

class Rtree
{
// data members
public:
    Memory&     m_memory;       // ref. to mem that stores this Rtree
    const int   m_dimen;        // no. of dimensions
    const int   m_maxNodeChild; // max. number of entries in a node
    const int   m_minNodeChild; // min. number of entries in a node
    const int   m_maxLeafChild; // max. number of entries in a leaf
    const int   m_minLeafChild; // min. number of entries in a leaf
    bool        m_pointOnly;    // indicate if the Rtree stores only points
// methods
public:
    // constructor/destructor
    Rtree(
        Memory& a_memory,
        const int a_dimen,
        const int a_maxNodeFill, const int a_maxLeafFill,
        const int a_minNodeFill, const int a_minLeafFill,
        const bool a_pointOnly=false);
    virtual ~Rtree();
    // update
    int insert(
        const RtreeNodeEntry& a_entry,  // add an entry to Rtree at
        const int a_level=0);           // a specified level (def: 0)
    int remove(                     // del an entry from Rtree
        const RtreeNodeEntry& a_entry);
    // test
    bool integrityTest() const;     // test rtree integrity
    Hash* loadObjects() const;      // test if all objects are stored
    float nodeVolume(               // find the vol. of nodes at level a_level
        const int a_level) const;
    float nodePerimeter(            // find the peri of nodes at level a_level
        const int a_level) const;
    int nodeCount(                  // find the no. of nodes at level a_level
        const int a_level) const;
protected:
    // supporting methods
    RtreeNode* chooseLeaf(
        const RtreeNodeEntry& a_entry,  // find leaf (node) of Rtree
        const int a_level);             // to insert an entry
    RtreeNode* locateLeaf(              // find a leaf that contains
        const RtreeNodeEntry& a_entry); // a_entry
};

#endif // RTREE_DEFINED
