/* ----------------------------------------------------------------------------
    Author: Ken C. K. Lee
    Email:  cklee@cse.psu.edu
    Web:    http://www.cse.psu.edu/~cklee
    Date:   Nov, 2007

    Copyright(c) 2007
    This program is for non-commerical use only.

    This program:
    1. creates a Rtree based on TGS for "a point dataset"
    2. performs a number of tests on the created Rtree
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
        id x y
    -i: index file
    -o: output eps file
    -v: verbose mode on
---------------------------------------------------------------------------- */

#include "stdlib.h"
#include "stdio.h"
#include "rtree.h"
#include "rnode.h"
#include "rentry.h"
#include "filemem.h"
#include "mainmem.h"
#include "hypercube.h"
#include "tgs.h"
#include "param.h"
#include "search.h"
#include "math.h"
#include "string.h"
#include <fstream>
#include <iostream>
#include <sys/types.h>
#include <sys/timeb.h>

void helpmsg(const char* pgm)
{
    cerr << "Suggested arguments:" << endl;
    cerr << "> " << pgm << endl;
    cerr << "-p 4096 -d 2 -f raw_data.txt -i index_file.idx -v" << endl;
    cerr << "explanations:" << endl;
    cerr << "-p: pagesize, typically, 4096" << endl;
    cerr << "-d: data dimensionality, e.g. 2" << endl;
    cerr << "-f: a tab delimited file, format:" << endl;
    cerr << "    id x y" << endl;
    cerr << "-i: index file" << endl;
    cerr << "-o: output eps file" << endl;
    cerr << "-v: verbose mode on" << endl;
}

int main(const int a_argc, const char** a_argv)
{
    if (a_argc == 1)
    {
        helpmsg(a_argv[0]);
        return -1;
    }

    cerr << "TGS Rtree" << endl;
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
    cerr << "load entries starts ... " << endl;
    RtreeNodeEntry** p = new RtreeNodeEntry*[10000000];
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
        p[objcnt++] = new RtreeNodeEntry(id, hc);

        if (verbose)
        {
            if (objcnt % 100 == 0)
                cerr << ".";
            if (objcnt % 1000 == 0)
                cerr << " " << objcnt << " point objects loaded." << endl;
        }
    }
    fin1.close();
    cerr << objcnt << " objects are loaded." << endl;

    //-------------------------------------------------------------------------
    // create a Rtree based on TGS
    //-------------------------------------------------------------------------
    cerr << "bulkloading ... ";
    const int maxChildNode =
        (pagesize - RtreeNode::size()) /
        RtreeNodeEntry::size(dimen);            // no. of entries per node
    const int maxChildLeaf =
        (pagesize - RtreeNode::size()) /
        RtreeNodeEntry::size(dimen,true);       // no. of entries per leaf node

    //MainMemory mem(pagesize);
    FileMemory mem(pagesize, idxname, RtreeNodeEntry::fromMem, true);
    struct timeb starttime, endtime;
    ftime(&starttime);  // time the algorithm
    Rtree* rtree = 
        TGS::bulkload(mem, dimen, maxChildNode, maxChildLeaf,
        (int)(maxChildNode*0.3), (int)(maxChildLeaf*0.3), p, objcnt,
        true); // true: point dataset
    ftime(&endtime);
    float createtime = 
        ((endtime.time*1000 + endtime.millitm) -
        (starttime.time*1000 + starttime.millitm)) / 1000.0f;
    cerr << "[DONE]" << endl;

    //-------------------------------------------------------------------------
    // test 1. draw all MBRs of Rtree
    //-------------------------------------------------------------------------
    cerr << "output the Rtree to " << outname << " ...";
    Search::dump(*rtree, 0, 100, outname);
    cerr << "[DONE]" << endl;

    //-------------------------------------------------------------------------
    // test 2. find the no. of leaf nodes, their total areas, total perimeters
    //-------------------------------------------------------------------------
    float area = rtree->nodeVolume(0);
    float peri = rtree->nodePerimeter(0);
    int cnt = rtree->nodeCount(0);
    cerr << "leaf node cnt,area,perimeter: ";
    cerr << cnt << ", " << area << ", " << peri << endl;

    //-------------------------------------------------------------------------
    // test 3. test Rtree integrity
    //-------------------------------------------------------------------------
    cerr << "rtree structure integrity test ... ";
    bool ret = rtree->integrityTest();
    if (ret)
        cerr << "[PASSED]" << endl;
    else
        cerr << "[FAILED]" << endl;

    //-------------------------------------------------------------------------
    // test 4. test if Rtree contains all inserted objects
    //-------------------------------------------------------------------------
    cerr << "object test ... ";
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
            cerr << "[FAILED]" << endl;
            cerr << "object " << id << " not found." << endl;
            result = false;
        }
        else if (!(e == *p))
        {
            cerr << "[FAILED]" << endl;
            cerr << "object " << id << " not matched." << endl;
            result = false;
        }
    }
    fin2.close();
    if (result)
        cerr << "[PASSED]" << endl;
    for (HashReader rdr(*h); !rdr.isEnd(); rdr.next())
    {
        RtreeNodeEntry* e = (RtreeNodeEntry*)rdr.getVal();
        delete e;
    }
    delete h;

    //-------------------------------------------------------------------------
    // test 5. test if the performance of range queries
    //-------------------------------------------------------------------------
    cerr << "query test ... ";
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
        int maxStackSize=0;
        Search::range(*rtree, pt, 50, res, maxStackSize, page);
        for (int i=0; i<res.size(); i++)
        {
            RtreeNodeEntry* r = (RtreeNodeEntry*)res[i];
            if (verbose)
                cerr << r->m_id << " ";
            delete r;
        }
        if (verbose)
            cerr << endl;
        pageaccessed += page.size();
        maxstack = maxStackSize > maxstack ? maxStackSize : maxstack;
    }
    cerr << "[DONE]" << endl;
    cerr << "search cost page, stacksize: ";
    cerr << pageaccessed << ", " << maxstack << endl;

    //-------------------------------------------------------------------------
    // report
    //-------------------------------------------------------------------------
    RtreeNode* node = rtree->m_memory.loadPage(rtree->m_memory.m_rootPageID);
    int level = node->m_level;
    int nodecnt = 0;
    for (int i=0; i<=level; i++)
        nodecnt += rtree->nodeCount(i);
    cout << "#nodes:," << nodecnt << ",createtime:," << createtime << endl;


    //-------------------------------------------------------------------------
    // clean up
    //-------------------------------------------------------------------------
    delete rtree;
    return 0;
}
