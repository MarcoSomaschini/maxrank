class RtreeNode;

class RtreeNodeEntry;

class VirtualRNode     //this is a class for data-entry node
{
public:
    int m_pageid;       // page ID of the node
    int m_level;        // level in the tree (0: leaf)
    int m_parent;       // page ID of its parent node
    int m_usedSpace;    // number of entries
    RtreeNodeEntry **m_entry;        // an array of entries

    bool m_isLeaf;

    //methods
    VirtualRNode() {};

    ~VirtualRNode();

    int copyData(const RtreeNode &source);

    int copyData(const VirtualRNode *source);

    //RtreeNodeEntry* cloneEntry();
    int copyEntries(const RtreeNode &source, int numbers);

    int insertEntry(const RtreeNodeEntry *source);

    int isLeaf();

    int displayMBR();
};
