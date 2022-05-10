/* ----------------------------------------------------------------------------
    This program:
    1. creates a aRtree based on TGS (aRtree consists of counts in tree nodes)
       for "a point dataset"
    2. performs a number of tests on the created aRtree
       1. print a rtree on eps file
       2. leaf node count, area, perimeter
       3. integrity test
       4. coverage test
       5. search performance

    Suggested arguments:
    > (prog name) -p 4096 -d 2 -f raw_data.txt -i index_file.idx -v
    explanations:
    -p: pagesize, typically, 4096
    -d: data dimensionality, e.g. 2
    -f: a tab delimited file, format:
        id xmin ymin xmax ymax
    -i: index file
    -o: output eps file
    -v: verbose mode on
---------------------------------------------------------------------------- */

#include "rtree.h"
#include "rnode.h"
#include "arentry.h"
#include "filemem.h"
#include "mainmem.h"
#include "hypercube.h"
#include "tgs.h"
#include "param.h"
#include "search.h"
#include <math.h>
#include <fstream>
#include <iostream>

void helpmsg(const char* pgm)
{
    cout << "Suggested arguments:" << endl;
    cout << "> " << pgm << endl;
    cout << "-p 4096 -d 2 -f raw_data.txt -i index_file.idx -v" << endl;
    cout << "explanations:" << endl;
    cout << "-p: pagesize, typically, 4096" << endl;
    cout << "-d: data dimensionality, e.g. 2" << endl;
    cout << "-f: a tab delimited file, format:" << endl;
    cout << "    id xmin ymin xmax ymax" << endl;
    cout << "-i: index file" << endl;
    cout << "-o: output eps file" << endl;
    cout << "-v: verbose mode on" << endl;
}

int main(const int a_argc, const char** a_argv)
{
    if (a_argc == 1)
    {
        helpmsg(a_argv[0]);
        return -1;
    }

    cout << "TGS aRtree" << endl;
    //-------------------------------------------------------------------------
    // initialization
    //-------------------------------------------------------------------------
    int dimen = atol(Param::read(a_argc, a_argv, "-d", ""));
    int pagesize = atol(Param::read(a_argc, a_argv, "-p", ""));
    const char* filename = Param::read(a_argc, a_argv, "-f", "");
    const char* idxname = Param::read(a_argc, a_argv, "-i", "");
    const char* outname = Param::read(a_argc, a_argv, "-o", "out.eps");
    const char* vrbs = Param::read(a_argc, a_argv, "-v", "null");
    bool verbose = strcmp(vrbs,"null") != 0;
    float ml[MAXDIMEN], mu[MAXDIMEN];
    fstream fin1, fin2;

    //-------------------------------------------------------------------------
    // load all objects into an array RtreeNodeEntry
    //-------------------------------------------------------------------------
    cout << "load entries starts " << endl;
    RtreeNodeEntry** p = new RtreeNodeEntry*[1000000];
    fin1.open(filename, ios::in);
    int objcnt = 0;
    while (true)
    {
        int id;
        float c[MAXDIMEN];
        fin1 >> id;
        if (fin1.eof()) break;

        for (int d=0; d<dimen; d++)
        {
            fin1 >> c[d];
            ml[d] = c[d] < ml[d] || objcnt == 0 ? c[d] : ml[d];
            mu[d] = c[d] > mu[d] || objcnt == 0 ? c[d] : mu[d];
        }

        Hypercube hc(dimen, c, c);
        p[objcnt++] = new ARtreeNodeEntry(id, hc, 1);

        if (verbose)
        {
            if (objcnt % 100 == 0)
                cout << ".";
            if (objcnt % 1000 == 0)
                cout << " " << objcnt << " objects loaded." << endl;
        }
    }
    fin1.close();
    cout << objcnt << " objects are loaded." << endl;

    //-------------------------------------------------------------------------
    // create a Rtree based on TGS
    //-------------------------------------------------------------------------
    cout << "bulkloading ...";
    const int maxChildNode =
        (pagesize - RtreeNode::size()) /
        ARtreeNodeEntry::size(dimen);           // no. of entries per node
    const int maxChildLeaf =
        (pagesize - RtreeNode::size()) /
        ARtreeNodeEntry::size(dimen,true);      // no. of entries per leaf node

    //MainMemory mem(pagesize);
    FileMemory mem(pagesize, idxname, ARtreeNodeEntry::fromMem, true);
    Rtree* rtree = 
        TGS::bulkload(mem, dimen, maxChildNode, maxChildLeaf,
        (int)(maxChildNode*0.3), (int)(maxChildLeaf*0.3), p, objcnt, true);
    cout << "[DONE]" << endl;

    //-------------------------------------------------------------------------
    // test 1. draw all MBRs of Rtree
    //-------------------------------------------------------------------------
    cout << "output the Rtree to " << outname << " ...";
    Search::dump(*rtree, 0, 100, outname);
    cout << "[DONE]" << endl;

    //-------------------------------------------------------------------------
    // test 2. find the no. of leaf nodes, their total areas, total perimeters
    //-------------------------------------------------------------------------
    float area = rtree->nodeVolume(0);
    float peri = rtree->nodePerimeter(0);
    int cnt = rtree->nodeCount(0);
    cout << "leaf node cnt,area,perimeter: ";
    cout << cnt << ", " << area << ", " << peri << endl;

    //-------------------------------------------------------------------------
    // test 3. test Rtree integrity
    //-------------------------------------------------------------------------
    cout << "rtree structure integrity test ... ";
    bool ret = rtree->integrityTest();
    if (ret)
        cout << "[PASSED]" << endl;
    else
        cout << "[FAILED]" << endl;

    //-------------------------------------------------------------------------
    // test 4. test if Rtree contains all inserted objects
    //-------------------------------------------------------------------------
    cout << "object test ... ";
    bool result=true;
    fin2.open(filename, ios::in);
    Hash* h = rtree->loadObjects();
    while (true)
    {
        int id;
        float c[MAXDIMEN];
        fin2 >> id;
        if (fin2.eof()) break;

        for (int d=0; d<dimen; d++)
            fin2 >> c[d];
        Hypercube hc(dimen, c, c);
        RtreeNodeEntry e(id, hc);
        RtreeNodeEntry* p = (RtreeNodeEntry*)h->get(id);
        if (p == 0)
        {
            cout << "[FAILED]" << endl;
            cout << "object " << id << " not found." << endl;
            result = false;
        }
        else if (!(e == *p))
        {
            cout << "[FAILED]" << endl;
            cout << "object " << id << " not matched." << endl;
            result = false;
        }
    }
    fin2.close();
    if (result)
        cout << "[PASSED]" << endl;
    for (HashReader rdr(*h); !rdr.isEnd(); rdr.next())
    {
        RtreeNodeEntry* e = (RtreeNodeEntry*)rdr.getVal();
        delete e;
    }
    delete h;

    //-------------------------------------------------------------------------
    // test 5. test if the performance of range queries
    //-------------------------------------------------------------------------
    cout << "query test ... ";
    int pageaccessed=0;
    int maxstack=0;
    for (int i=0; i<1000; i++)
    {
        float c[MAXDIMEN];
        for (int d=0; d<dimen; d++)
            c[d] = (rand() % 100) / 100.0f * (mu[d] - ml[d]) + ml[d];
        Point pt(dimen, c);
        Array res;
        Array page;
        int maxStackSize;
        Search::range(*rtree, pt, 50, res, maxStackSize, page);
        for (int i=0; i<res.size(); i++)
        {
            RtreeNodeEntry* r = (RtreeNodeEntry*)res[i];
            if (verbose)
                cout << r->m_id << " ";
            delete r;
        }
        if (verbose)
            cout << endl;
        pageaccessed += page.size();
        maxstack = maxStackSize > maxstack ? maxStackSize : maxstack;
    }
    cout << "[DONE]" << endl;
    cout << "search cost page, stacksize: ";
    cout << pageaccessed << ", " << maxstack << endl;

    //-------------------------------------------------------------------------
    // clean up
    //-------------------------------------------------------------------------
    delete rtree;
    return 0;
}
