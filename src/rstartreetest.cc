/* ----------------------------------------------------------------------------
    Author: Ken C. K. Lee
    Email:  cklee@cse.psu.edu
    Web:    http://www.cse.psu.edu/~cklee
    Date:   Nov, 2007

    Copyright(c) 2007
    This library is for non-commerical use only.

    This program:
    1. creates a R*tree based on approches described in SIGMOD90 R-tree paper
    2. performs a number of tests on the created Rtree
       1. print a rtree on an eps file
       2. leaf node count, area, perimeter
       3. integrity test
       4. coverage test
       5. search performance
       6. deletion (optional)

    Suggested arguments:
    > (prog name) -p 4096 -d 2 -f raw_data.txt -i index_file.idx -o out.eps -v
    explanations:
    -p: pagesize, typically, 4096
    -d: data dimensionality, e.g. 2
    -f: a tab delimited file, format:
        id xmin ymin xmax ymax
    -i: index file
    -o: output eps file
    -v: verbose mode on
---------------------------------------------------------------------------- */

#include "stdlib.h"
#include "stdio.h"
#include "rtree.h"
#include "rstartree.h"
#include "rnode.h"
#include "rentry.h"
#include "hypercube.h"
#include "param.h"
#include "filemem.h"
#include "mainmem.h"
#include "search.h"
#include <math.h>
#include <fstream>
#include <iostream>
#include "string.h"

void helpmsg(const char* pgm)
{
    cout << "Suggested arguments:" << endl;
    cout << "> " << pgm << endl;
    cout << "-p 4096 -d 2 -f raw_data.txt -i index_file.idx -o out.eps -v" << endl;
    cout << "explanations:" << endl;
    cout << "-p: pagesize, typically, 4096" << endl;
    cout << "-d: data dimensionality, e.g. 2" << endl;
    cout << "-f: a tab delimited file, format:" << endl;
    cout << "    id xmin ymin [zmin] xmax ymax [zmax]" << endl;
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

    cout << "R*tree" << endl;
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
    float ml[5], mu[5];
    fstream fin1, fin2, fin3;

    //-------------------------------------------------------------------------
    // create a R*tree
    //-------------------------------------------------------------------------
    cout << "create a R*tree" << endl;
    const int maxChild =
        (pagesize - RtreeNode::size()) /
        RtreeNodeEntry::size(dimen);        // max. no. of entries per node
    //MainMemory mem(pagesize);
    FileMemory mem(pagesize, idxname, RtreeNodeEntry::fromMem, true);
    Rstartree rtree(mem, dimen, maxChild, maxChild,
        (int)ceil(maxChild * 0.3), (int)ceil(maxChild * 0.3), (int)ceil(maxChild * 0.3), false); 

    //-------------------------------------------------------------------------
    // insert objects to R*tree
    //-------------------------------------------------------------------------
    fin1.open(filename, ios::in);
    cout << "insertion ... ";
    int objcnt = 0;
    while (true)
    {
        int id;
        float cl[5];
        float cu[5];
        fin1 >> id;
        if (fin1.eof()) break;

        for (int d=0; d<dimen; d++)
        {
            fin1 >> cl[d];
            ml[d] = cl[d] < ml[d] || objcnt == 0 ? cl[d] : ml[d];
        }
        for (int d=0; d<dimen; d++)
        {
            fin1 >> cu[d];
            mu[d] = cu[d] > mu[d] || objcnt == 0 ? cu[d] : mu[d];
        }
        Hypercube hc(dimen, cl, cu);
        RtreeNodeEntry e(id, hc);
        rtree.insert(e);


        objcnt++;
        if (verbose)
        {
            if (objcnt % 100 == 0)
                cout << ".";
            if (objcnt % 1000 == 0)
                cout << " " << objcnt << " objects inserted." << endl;
        }
    }
    cout << "[DONE]" << endl;
    fin1.close();

    //-------------------------------------------------------------------------
    // test 1. draw all MBRs of R*tree
    //-------------------------------------------------------------------------
    cout << "output the Rtree to " << outname << " ...";
    Search::dump(rtree, 0, 100, outname);
    cout << "[DONE]" << endl;

    //-------------------------------------------------------------------------
    // test 2. find the no. of leaf nodes, their total areas, total perimeters
    //-------------------------------------------------------------------------
    float area = rtree.nodeVolume(0);
    float peri = rtree.nodePerimeter(0);
    int cnt = rtree.nodeCount(0);
    cout << "leaf node cnt,area,perimeter: ";
    cout << cnt << ", " << area << ", " << peri << endl;

    //-------------------------------------------------------------------------
    // test 3. test R*tree integrity
    //-------------------------------------------------------------------------
    cout << "rtree structure integrity test ... ";
    bool ret = rtree.integrityTest();
    if (ret)
        cout << "[PASSED]" << endl;
    else
        cout << "[FAILED]" << endl;

    //-------------------------------------------------------------------------
    // test 4. test if R*tree contains all inserted objects
    //-------------------------------------------------------------------------
    cout << "object test ... ";
    bool result=true;
    fin2.open(filename, ios::in);
    Hash* h = rtree.loadObjects();
    while (true)
    {
        int id;
        float cl[5];
        float cu[5];
        fin2 >> id;
        if (fin2.eof()) break;

        for (int d=0; d<dimen; d++)
            fin2 >> cl[d];
        for (int d=0; d<dimen; d++)
            fin2 >> cu[d];
        Hypercube hc(dimen, cl, cu);
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
        float c[5];
        for (int d=0; d<dimen; d++)
            c[d] = (rand() % 100) / 100.0f * (mu[d] - ml[d]) + ml[d];
        Point pt(dimen, c);
        Array res;
        Array page;
        int maxStackSize;
        Search::range(rtree, pt, 50, res, maxStackSize, page);
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

    return 0;   // comment it if you want test 6.


    //-------------------------------------------------------------------------
    // test 6. test if R*tree can remove objects
    //-------------------------------------------------------------------------
    cout << "deletion test ... ";
    int objdel=0;
    fin3.open(filename, ios::in);
    while (true)
    {
        int id;
        float cl[5];
        float cu[5];
        fin3 >> id;
        if (fin3.eof()) break;

        for (int d=0; d<dimen; d++)
            fin3 >> cl[d];
        for (int d=0; d<dimen; d++)
            fin3 >> cu[d];
        Hypercube hc(dimen, cl, cu);
        RtreeNodeEntry e(id, hc);
        rtree.remove(e);

        objcnt--;
        objdel++;
        if (objdel % 100 == 0)
            cout << ".";
        if (objdel % 1000 == 0)
            cout << " " << objdel << " objects deleted." << endl;

        bool ret = rtree.integrityTest();
        if (!ret)
            cout << "removal failure at object " << id << endl;
        Hash* h = rtree.loadObjects();
        if (h->size() != objcnt)
            cout << "removal failure at object " << id << endl;
        RtreeNodeEntry* en = (RtreeNodeEntry*)h->get(id);
        if (en != 0)
            cout << "unable to delete object " << id << endl;
        for (HashReader rdr(*h); !rdr.isEnd(); rdr.next())
        {
            RtreeNodeEntry* e = (RtreeNodeEntry*)rdr.getVal();
            delete e;
        }
        delete h;
    }
    fin3.close();
    cout << "[DONE]" << endl;

    return 0;
}
