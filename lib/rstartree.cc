#include "rstartree.h"
#include "rnode.h"
#include "rentry.h"
#include "collection.h"
#include "mem.h"

Rstartree::Rstartree(Memory& a_memory,
                     const int a_dimen,
                     const int a_maxNodeFill, const int a_maxLeafFill,
                     const int a_minNodeFill, const int a_minLeafFill,
                     const int a_reinsert,
                     const bool a_pointOnly):
Rtree(a_memory, a_dimen, a_maxNodeFill, a_maxLeafFill, a_minNodeFill, a_minLeafFill, a_pointOnly),
m_reinsert(a_reinsert)
{}

Rstartree::~Rstartree()
{}

int Rstartree::insert(const RtreeNodeEntry& a_entry, const int a_level)
{
    class carrier
    {
    public:
        RtreeNodeEntry* m_entry;
        const int       m_level;
    public:
        carrier(RtreeNodeEntry* a_entry, const int a_level):
          m_entry(a_entry), m_level(a_level) {}
        ~carrier() {};
    };

    Set reinserted;
    Queue q;
    q.enqueue(new carrier(a_entry.clone(), a_level));
    while (!q.isEmpty())
    {
        carrier* c = (carrier*)q.dequeue();
        RtreeNodeEntry* entry = c->m_entry;
        int level = c->m_level;

        if (m_memory.m_rootPageID == -1)
        {   // no tree exists
            m_memory.m_rootPageID = m_memory.allocatePage();
            RtreeNode n(this, m_memory.m_rootPageID);
            n.insert(*entry);
            m_memory.writePage(m_memory.m_rootPageID, &n);
        }
        else
        {
            RtreeNode* n = chooseLeaf(*entry, level);
            int status = n->insert(*entry);
            if (level > 0)
            {
                RtreeNode* m = m_memory.loadPage(entry->m_id);
                m->m_parent = n->m_pageid;
                m_memory.writePage(m->m_pageid, m);
                delete m;
            }

            while (status != NODE_UNCHANGED)
            {
                RtreeNode* parent = 0;
                switch (status)
                {
                case NODE_OVERFLOW:
                    {
                        if (reinserted.in((void*)n->m_level) == 0 && !n->isRoot())
                        {
                            reinserted.insert((void*)n->m_level);
                            // reinsertion
                            RtreeNode* newnode;
                            RtreeNodeEntry** residue = new RtreeNodeEntry*[m_reinsert];
                            n->pickWorst(m_reinsert, &newnode, residue);
                            for (int i=0; i<m_reinsert; i++)
                                q.enqueue(new carrier(residue[i],n->m_level));
                            delete[] residue;
                            //
                            m_memory.writePage(newnode->m_pageid, newnode);
                            //
                            parent = m_memory.loadPage(n->m_parent);
                            RtreeNodeEntry* e = (RtreeNodeEntry*)newnode->genNodeEntry();
                            status = parent->replace(*e);
                            delete e;
                            delete n;
                            n = newnode;
                        }
                        else
                        {   // split if this level of Rtree has been reinserted.
                            RtreeNode *n0, *n1;
                            n->goodnessSplit(&n0, &n1);
                            n0->m_pageid = n->m_pageid;
                            n1->m_pageid = m_memory.allocatePage();
                            n->m_usedSpace=0;

                            if (!n->isLeaf())
                            {
                                // update the parent ids of the child nodes
                                for (int i=0; i<n0->m_usedSpace; i++)
                                {
                                    int id = n0->m_entry[i]->m_id;
                                    RtreeNode* child = m_memory.loadPage(id);
                                    child->m_parent = n0->m_pageid;
                                    m_memory.writePage(id, child);
                                    delete child;
                                }
                                for (int i=0; i<n1->m_usedSpace; i++)
                                {
                                    int id = n1->m_entry[i]->m_id;
                                    RtreeNode* child = m_memory.loadPage(id);
                                    child->m_parent = n1->m_pageid;
                                    m_memory.writePage(id, child);
                                    delete child;
                                }
                            }

                            if (n->isRoot())
                            {   
                                // create a new parent (root) node
                                int parentid = m_memory.allocatePage();
                                m_memory.m_rootPageID = parentid;       // set new root
                                parent = new RtreeNode(this, parentid, n->m_level+1, -1);
                                // create node entry e0 to parent node
                                RtreeNodeEntry* e0 = n0->genNodeEntry();
                                n0->m_parent = parentid;
                                parent->insert(*e0);
                                // create node entry e1 to parent node
                                RtreeNodeEntry* e1 = n1->genNodeEntry();
                                n1->m_parent = parentid;
                                parent->insert(*e1);
                                // write them all to the memory
                                m_memory.writePage(n0->m_pageid, n0);
                                m_memory.writePage(n1->m_pageid, n1);
                                status = NODE_UNCHANGED;
                                delete e0;
                                delete e1;
                            }
                            else
                            {
                                // load a parent node and update it
                                parent = m_memory.loadPage(n->m_parent);
                                // create node entry e0 to parent node
                                RtreeNodeEntry* e0 = n0->genNodeEntry();
                                n0->m_parent = parent->m_pageid;
                                int status1 = parent->replace(*e0);
                                // create node entry e1 to parent node
                                RtreeNodeEntry* e1 = n1->genNodeEntry();
                                n1->m_parent = parent->m_pageid;
                                int status2 = parent->insert(*e1);
                                //
                                m_memory.writePage(n0->m_pageid, n0);
                                m_memory.writePage(n1->m_pageid, n1);
                                status = status2;
                                if (status != NODE_OVERFLOW && status1 == NODE_CHANGED)
                                    status = NODE_CHANGED;
                                delete e0;
                                delete e1;
                            }
                            delete n0;
                            delete n1;
                        }
                        break;
                    }
                case NODE_CHANGED:
                    {
                        if (n->isRoot())
                        {
                            status = NODE_UNCHANGED;
                            parent = n->clone();
                            break;  // no parent to update; quit!
                        }
                        else
                        {
                            parent = m_memory.loadPage(n->m_parent);
                            RtreeNodeEntry* e = n->genNodeEntry();
                            status = parent->replace(*e);
                            m_memory.writePage(n->m_pageid, n);
                            delete e;
                        }
                    }
                };
                delete n;
                n = parent;
            }
            m_memory.writePage(n->m_pageid, n);
            delete n;
        }
        delete entry;
        delete c;

    }
    return 0;
}
