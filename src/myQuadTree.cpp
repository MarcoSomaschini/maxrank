//#include "QuadTree.h"
#include "Objects.h"
#include "stdio.h"
#include "stdlib.h"
#include "time.h"
#include "math.h"
#include "string.h"
#include <set>
#include <iostream>
#include <fstream>
#include <string>

//headers
#include "collection.h"
#include "rtree.h"
#include "rnode.h"
#include "rentry.h"
#include "hypercube.h"
#include "param.h"
#include "filemem.h"
#include "mainmem.h"
#include "S3.h"
#include "iomeasure.h"
#include "tgs.h"
#include "random.h"
#include "virtualRnode.h"
#include <map>
//end 



#define MAXPTS 1000000   //maximum size of the database
#define MB 1048576

using namespace std;


vector<vector<float> > HalfSpaces;  //format: (coeff_1, coeff_2, ..., coeff_d, offset)
long int NoOfNewlyAddedHalfSpaces;
long int NoOfTotalHalfspacesAdded;
map<long int, long int> entriesInLists;
long int NoOfEntriesInLists;
vector<vector<float> > IncpDB;
int maximumLevel = 0;
long numOfTreeNodes = 0;
long numOfLeafNodes = 0;
long numOfLeavesToIntersect = 0;
long numOfInvalidNodes = 0;
float memRatioOfQuadTree = 0;
long numOfNodeSplits = 0;
int normalizedMax = 1;
int numOfSubdivisions = 0;
float Queryspace[2 * MAXDIM + 1];
int MaxQuadTreeLevels = 6;
int QuadNodeCapacity = 25;
bool verbose;

double timeBuildQuadTree = 0;
double timeNodeIntersection = 0;

long int totalNoOfBitStringsProcessed = 0;
long int totalNoOfPrunedBitStrings = 0;
long int totalNoOfZeroExtentBinStrings = 0;
long int totalNoOfDiscardedCells = 0;

int tao = 0;

bool optWithinNodeIntersection = false;

randomnumber rnd(0);
int CEILING = 1;   //ceiling of the dataspace

//code for computing global immutable region by using R-tree index_file
void helpmsg(const char *pgm) {
    cout << "Suggested arguments:" << endl;
    cout << "> " << pgm << endl;
    cout << "-p 4096 -d 2 -f raw_data.txt -i index_file.idx -o out.eps -v" << endl;
    cout << "explanations:" << endl;
    cout << "-p: pagesize, typically, 4096" << endl;
    cout << "-d: data dimensionality, e.g. 2" << endl;
    cout << "-f: a tab delimited file, format:" << endl;
    cout << "    id xmin ymin [zmin] xmax ymax [zmax]" << endl;
    cout << "-q: query file, format:" << endl;
    cout << "    id xmin ymin xmax ymax" << endl;
    cout << "-k: k for top-k" << endl;
    cout << "-i: index file" << endl;
    cout << "-r: times to repeat" << endl;
    cout << "-m: type of method; AA: AA algorithm, BA: BA algorithm" << endl;
    cout << "-h: maximum height of the quad tree" << endl;
    cout << "-t: parameter \tao for iMaxRank query" << endl;
    cout << "-n: output non-result points after top-k computation" << endl;
    cout << "-o: 0: use naive within-node intersection; 1: use optimized within-node intersection" << endl;
    cout << "-v: verbose mode on" << endl;
}

void myitoa(unsigned long val, char *buf, unsigned radix) {
    char *p;
    char *firstdig;
    char temp;
    unsigned digval;

    p = buf;
    firstdig = p;

    do {
        digval = (unsigned) (val % radix);
        val /= radix;

        if (digval > 9)
            *p++ = (char) (digval - 10 + 'a');
        else
            *p++ = (char) (digval + '0');

    } while (val > 0);

    *p-- = '\0';
    do {
        temp = *p;
        *p = *firstdig;
        *firstdig = temp;
        --p;
        ++firstdig;
    } while (firstdig < p);
}

void
GenString(long int stringLen, long int HammingDistance, long int currentLen, long int start, vector<char> &hammingStr,
          multimap<int, vector<char> > &binString) {
    if (currentLen < 0) return;

    for (long int i = start; i < stringLen; i++) {
        for (long int j = start; j < i; j++)
            hammingStr.push_back('0');
        hammingStr.push_back('1');
        GenString(stringLen, HammingDistance, currentLen - 1, i + 1, hammingStr, binString);
        if (currentLen == 0) {
            for (long int j = i + 1; j < stringLen; j++)
                hammingStr.push_back('0');

            vector<char> tmpHammingStr = hammingStr;
            typedef multimap<int, vector<char> >::value_type VT;
            binString.insert(VT(HammingDistance, tmpHammingStr));

            //cout << "Generated Hamming string: " << string(tmpHammingStr.begin(),tmpHammingStr.end()) << endl;

            for (long int j = i + 1; j < stringLen; j++)
                hammingStr.pop_back();
        }
        hammingStr.pop_back();
        for (long int j = start; j < i; j++)
            hammingStr.pop_back();
    }
}

void GenLenNBinaryString(long int len1, long int HammingDistance, multimap<int, vector<char> > &binString) {
    vector<char> hammingStr;

    if (HammingDistance == 0) {
        for (long int i = 0; i < len1; i++) hammingStr.push_back('0');
        typedef multimap<int, vector<char> >::value_type VT;
        binString.insert(VT(0, hammingStr));
        return;
    }

    if (HammingDistance == 1) {
        for (long int i = 0; i < len1; i++) {
            hammingStr.clear();
            for (long int j = 0; j < i; j++)
                hammingStr.push_back('0');
            hammingStr.push_back('1');
            for (long int j = i + 1; j < len1; j++)
                hammingStr.push_back('0');

            typedef multimap<int, vector<char> >::value_type VT;
            binString.insert(VT(1, hammingStr));
        }
        return;
    }


    GenString(len1, HammingDistance, HammingDistance - 1, 0, hammingStr, binString);

    return;
}


void GenBinaryString(long int len1, long int Max, multimap<int, vector<char> > &binString) {
    //ofstream fOut;
    //fOut.open("CombSubspace.txt", ios::app);

    typedef multimap<int, vector<char> >::value_type VT;

    binString.clear();

    char str[128], *b;
    //int len1=int(log(Max)/log(2));
    for (int i = 0; i < Max; i++) {
        myitoa(i, str, 2);
        int len = len1 - strlen(str);
        b = new char[Max + 1];
        strcpy(b, "");
        for (int j = 1; j <= len; j++) strcat(b, "0");
        strcat(b, str);
        cout << "i= " << i << ", binary string: " << b << endl;
        //fOut << b << endl;

        //count the number of '1's
        vector<char> tmpVec;
        int NoOfOnes = 0;
        for (int j = 0; j < strlen(b); j++) {
            if (b[j] == '1') NoOfOnes++;
            tmpVec.push_back(b[j]);
        }
        binString.insert(VT(NoOfOnes, tmpVec));

        delete b;
    }
    //fOut.close();

    return;
}

bool PointCoveredByMBR(const int &Dimen, const float mbr[], const Point1 &pt) {

    bool isCovered = true;

    for (int i = 0; i < Dimen; i++) {
        if ((pt.coord[i] <= mbr[i]) || (pt.coord[i] >= mbr[Dimen + i])) {
            isCovered = false;
            break;
        }
    }

    return isCovered;
};

int PointVersusHalfSpace(const int &Dimen, const float hs[],
                         const Point1 &pt) {   //position of a point to a halfspace: is the point above, below, or on the halfspace?

    bool isAbove = false;

    float sum = 0;

    for (int i = 0; i < Dimen; i++) {
        sum = sum + hs[i] * pt.coord[i];
    }

    // The absolute difference is considered to avoid approximation errors
    if (fabs(sum - hs[Dimen]) <= 1e-6)
        return ON;
    if (sum < hs[Dimen])
        return ABOVE;
    if (sum > hs[Dimen])
        return BELOW;

    return ON;
}

int MbrVersusHalfSpace(const int &Dimen, const float hs[], const float mbr[],
                       vector<string> &Comb) {   //position of an MBR to a halfspace: is the point above, below, or intersected by the halfspace?

    int numAbove = 0;
    int numBelow = 0;
    int numOn = 0;

    long int numOfVertices = 0;
    numOfVertices = Comb.size();

    for (int i = 0; i < numOfVertices; i++) {
        Point1 pt;

        long int numOfDimen = Comb[i].size();
        for (int j = 0; j < numOfDimen; j++) {
            if (Comb[i][j] == '0')
                pt.coord[j] = mbr[j];
            if (Comb[i][j] == '1')
                pt.coord[j] = mbr[Dimen + j];
        }
        int pos = PointVersusHalfSpace(Dimen, hs, pt);
        if (pos == ABOVE) numAbove++;
        if (pos == BELOW) numBelow++;
        if (pos == ON) numOn++;
    }

    // numOn is summed to account for the halfspaces that lie on one of the leaf sides
    if (numAbove + numOn == numOfVertices) return ABOVE;
    if (numBelow + numOn == numOfVertices) return BELOW;

    return OVERLAPPED;
}

bool MbrIsValid(const int &Dimen, const float hs[], const float mbr[],
                vector<string> &Comb) {   //position of an MBR to a halfspace: is the point above, below, or intersected by the halfspace?

    int numAbove = 0;
    int numBelow = 0;
    int numOn = 0;

    long int numOfVertices = 0;
    numOfVertices = Comb.size();

    for (int i = 0; i < numOfVertices; i++) {
        Point1 pt;

        long int numOfDimen = Comb[i].size();
        for (int j = 0; j < numOfDimen; j++) {
            if (Comb[i][j] == '0')
                pt.coord[j] = mbr[j];
            if (Comb[i][j] == '1')
                pt.coord[j] = mbr[Dimen + j];
        }
        float sum = 0;
        for (int k = 0; k < numOfDimen; k++) sum = sum + pt.coord[k];
        if (sum > hs[Dimen]) numAbove++;
        if (sum < hs[Dimen]) numBelow++;
    }

    if (numAbove == numOfVertices) return false;
    if (numBelow == numOfVertices) return true;

}

int readHalfSpaces(const char *FileName, const int &Dimen) {

    FILE *fp;
    char *token;
    char m_separator[] = " \n\t";
    char buf[512];
    long int numOfHalfSpaces = 0;
    int DimenOfHalfSpace;

    float Min[MAXDIM];
    float Max[MAXDIM];
    std::fill(Min, Min + MAXDIM, 999999);
    std::fill(Max, Max + MAXDIM, -999999);

    fp = fopen(FileName, "r");
    if (fp == NULL) {
        cout << "error in file opening!" << endl;
        exit(0);
    }

    while (fgets(buf, 512, fp) != NULL) {
        token = strtok(buf, m_separator);
        vector<float> hs;                      //a single halfspace
        numOfHalfSpaces++;
        DimenOfHalfSpace = 0;
        while (token != NULL) {
            float tmp = atof(token);
            hs.push_back(tmp);
            if (Min[DimenOfHalfSpace] > tmp) Min[DimenOfHalfSpace] = tmp;
            if (Max[DimenOfHalfSpace] < tmp) Max[DimenOfHalfSpace] = tmp;
            DimenOfHalfSpace++;
            token = strtok(NULL, m_separator);
        }
        HalfSpaces.push_back(hs);
        if (DimenOfHalfSpace != Dimen + 1) {
            cout << "Caution!! The dimensionality of halfspaces is not equal! Halfspace ID = " << endl;
        }
    }
    fclose(fp);

    //normalize halfspaces
    int OutterSz = HalfSpaces.size();
    int InnerSz = HalfSpaces[0].size();
    for (int i = 0; i < InnerSz; i++) {
        float min = Min[i];
        float max = Max[i];
        for (int j = 0; j < OutterSz; j++)
            HalfSpaces[j][i] = (HalfSpaces[j][i] - min) / (max - min);
    }
    //end of normalize halfspaces

    cout << numOfHalfSpaces << " halfspaces have been read!" << endl;

    return 0;
}

int buildHalfSpaces(const char *FileName, const int &Dimen,
                    const Point1 &p) {    //read in the incomparable records, i.e., those not dominate or dominated by query record p

    FILE *fp;
    char *token;
    char m_separator[] = " \n\t";
    char buf[512];
    long int numOfRds = 0, numOfIncpRds = 0;
    int DimenOfRs;

    fp = fopen(FileName, "r");
    if (fp == NULL) {
        cout << "error in opening database file!" << endl;
        exit(0);
    }

    IncpDB.clear();
    while (fgets(buf, 512, fp) != NULL) {
        token = strtok(buf, m_separator);
        vector<float> Rd;                      //a single halfspace
        DimenOfRs = 0;
        while (token != NULL) {
            Rd.push_back(atof(token));
            DimenOfRs++;
            token = strtok(NULL, m_separator);
        }
        if (DimenOfRs != Dimen) {
            cout << "Caution!! The dimensionality of records is not equal! Record ID = " << endl;
        }
        numOfRds++;

        /*
//check whether p is dominated by current record
int dominate = 0;
for (int i=0; i<Dimen; i++)
     if (p.coord[i]>=Rd[i]) dominate++;
        if (dominate >= Dimen) continue;     //current record is dominated by p, so we discard it

//check whether p dominates current record
int isDominatedBy = 0;
for (int i=0; i<Dimen; i++)
     if (p.coord[i]<=Rd[i]) isDominatedBy++;
        if (isDominatedBy >= Dimen) continue;     //current record dominates p, so we discard it too
        //*/

        IncpDB.push_back(Rd);    //found an incomparable record
        numOfIncpRds++;
        vector<float> tmpHS;
        float Rd_d = Rd[Dimen - 1];
        float p_d = p.coord[Dimen - 1];
        for (int i = 0; i < Dimen - 1; i++) {
            tmpHS.push_back((Rd[i] - Rd_d - p.coord[i] + p_d));
        }
        tmpHS.push_back(p_d - Rd_d);
        tmpHS.push_back(numOfRds);            //store the ID of incomparable record in HalfSpace
        HalfSpaces.push_back(tmpHS);          //form the half-space defined by the incomparable record and p
    }
    fclose(fp);

    //cout << numOfIncpRds << " Incomparable records have been found in " << numOfRds << " records." << endl;
    cout << numOfIncpRds << " half-spaces have been constructed from skylines!" << endl;

    return 0;
}

void findTopKfromDB(float *PG[], long int objcnt, int Dimen, int *pos) {
    long int i;

    //ofstream fout;
    //fout.open("testOut.txt",ios::out);

    double point[MAXDIM];

    map<float, long int, greater<float> > Records;
    map<float, long int, greater<float> >::iterator mItr;
    typedef map<float, long int, greater<float> >::value_type mVT;


    for (i = 0; i < objcnt; i++) {
        float score = 0;
        for (int j = 0; j < Dimen; j++)
            score = score + (PG[i + 1][j] + PG[i + 1][Dimen + j]) / 2.0;
        Records.insert(mVT(score, (i + 1)));
    }

    /*
    fstream fin;
    fin.open(dbFile,ios::in);
    while (true)
    {
        fin >> id;
        if (fin.eof()) break;

        float score=0;
        for (i=0;i<Dimen;i++)
        {
             fin >> point[i];
             score=score+point[i];
        }
        Records.insert(mVT(score,id));
    }
    fin.close();
    //*/

    long int counter = 0;
    int idx = 0;
    for (mItr = Records.begin(); mItr != Records.end(); mItr++) {
        //fout << "ID=" << mItr->first << ", score="<< mItr->second << endl;

        counter++;
        switch (counter) {
            case 5:
            case 10:
            case 50:
            case 100:
            case 500:
            case 1000:
                pos[idx] = mItr->second;
                idx++;
        }
    }
    //fout.close();
    return;
}

int main(const int a_argc, const char **a_argv) {
    if (a_argc == 1) {
        helpmsg(a_argv[0]);
        return -1;
    }

    cout << "Rtree" << endl;
    //-------------------------------------------------------------------------
    // initialization
    //-------------------------------------------------------------------------
    int dimen = atoi(Param::read(a_argc, a_argv, "-d", ""));
    int pagesize = atol(Param::read(a_argc, a_argv, "-p", ""));
    int m_K = atoi(Param::read(a_argc, a_argv, "-k", ""));
    int Repeat = atoi(Param::read(a_argc, a_argv, "-r", "1"));
    int Height = atoi(Param::read(a_argc, a_argv, "-h", "5"));
    tao = atoi(Param::read(a_argc, a_argv, "-t", "0"));
    optWithinNodeIntersection = atoi(Param::read(a_argc, a_argv, "-o", ""));
    const char *filename = Param::read(a_argc, a_argv, "-f", "");
    const char *qfilename = Param::read(a_argc, a_argv, "-q", "");
    const char *idxname = Param::read(a_argc, a_argv, "-i", "");
    const char *qhalfoutname = Param::read(a_argc, a_argv, "-h", "");
    const char *nonResultname = Param::read(a_argc, a_argv, "-n", "");
    const char *vrbs = Param::read(a_argc, a_argv, "-v", "null");
    const char *methodName = Param::read(a_argc, a_argv, "-m", "AA");

    verbose = strcmp(vrbs, "null") != 0;
    float ml[5], mu[5];
    fstream fin0, fin1, fin2, fin3;

    ofstream fperfOut;
    fperfOut.open("result/Perf.txt", ios::app);

    //my variables for measuring the performance
    float MemUsage = 0;
    float MemTree = 0;
    float tmpMemUsage = 0;
    float TimeCost = 0;
    long int IOs = 0;
    long int totalNoOfMinCells = 0;
    long int totalMin_K = 0;

    long int TotalMemQuadTree = 0;
    NoOfTotalHalfspacesAdded = 0;
    long int maxNoOfQuadTreeNodes = -1, maxNoOfEntriesInLists = -1, maxNoOfTotalHalfspacesAdded = -1, maxNoOfLeavesToIntersect = -1;

    Array a_PagesRead;
    float AvgIOBufferP1 = 0;
    float AvgIOBufferP2 = 0;
    float AvgIOBufferP5 = 0;
    float AvgIOBufferP10 = 0;
    float AvgIOBufferP20 = 0;
    //end of my variables

    float ds_x0, ds_x1, ds_y0, ds_y1;    //the lower-left and -upper-right corners of the whole data space
    float avg_area, avg_Tmp;

    //initialize variables for QuadTree; note that the dimension of the QuadTree is one dimension less than the original dimension
    MaxQuadTreeLevels = Height;

    for (int d = 0; d < dimen - 1; d++)
        Queryspace[d] = 0;
    for (int d = 0; d < dimen - 1; d++)
        Queryspace[dimen - 1 + d] = normalizedMax;
    numOfSubdivisions = (int) pow(2.0, dimen -
                                       1);   //the number of QuadTree is one dimension less than that of the problem original setting
    //end of QuadTree variables

    char dbFile[1024];
    int posName[6] = {5, 10, 50, 100, 500, 1000};

    //fin0.open(filename,ios::in);
    //while (true)
    //{
    //   fin0 >> dbFile;     //read in the database file
    //   if (fin0.eof()) break;

    //   cout << "Processing file " << dbFile << endl;

    //-------------------------------------------------------------------------
    // bulkload all objects into an array RtreeNodeEntry
    //-------------------------------------------------------------------------
    cout << "begin loading data objects..." << endl;
    char buf[512];
    RtreeNodeEntry **p = new RtreeNodeEntry *[MAXPTS];
    float **PG = new float *[MAXPTS + 1];
    fin1.open(filename, ios::in);
    int objcnt = 0;
    ds_x0 = 999999;
    ds_x1 = -999999;
    ds_y0 = 999999;
    ds_y1 = -999999;
    avg_area = 0;
    while (true) {
        int id;
        float cl[MAXDIMEN], cu[MAXDIMEN];
        fin1 >> id;
        if (fin1.eof()) break;

        PG[objcnt + 1] = new float[2 * dimen];

        strcpy(buf, "");  //for partial hull

        for (int d = 0; d < dimen; d++) {
            fin1 >> cl[d];
            ml[d] = cl[d] < ml[d] || objcnt == 0 ? cl[d] : ml[d];

            PG[objcnt + 1][d] = cl[d];

            //store the string data for partial hull computation in Step 1 onwards  
            //sprintf(buf+strlen(buf),"%lf ",cl[d]+SIDELEN);
            //strcat(buf," ");
        }
        for (int d = 0; d < dimen; d++) {
            fin1 >> cu[d];
            mu[d] = cu[d] > mu[d] || objcnt == 0 ? cu[d] : mu[d];

            PG[objcnt + 1][dimen + d] = cu[d];
        }

        //int strLen=strlen(buf);
        //PGstr[objcnt+1]=new char [strLen+1];
        //strcpy(PGstr[objcnt+1],buf);

        Hypercube hc(dimen, cl, cu);
        p[objcnt++] = new RtreeNodeEntry(id, hc);

        /* 
        cout << "object "<< p[objcnt-1]->m_id << ": ";
        for (int j=0;j<dimen;j++)
            cout << p[objcnt-1]->m_hc.getLower()[j] << " " << p[objcnt-1]->m_hc.getUpper()[j]<< " ";
        cout << endl;
        getchar();    
        //*/

        //computing the data space and the average area of data MBRs
        if (cl[0] < ds_x0) ds_x0 = cl[0];
        if (cl[1] < ds_y0) ds_y0 = cl[1];
        if (cu[0] > ds_x1) ds_x1 = cu[0];
        if (cu[1] > ds_y1) ds_y1 = cu[1];
        avg_Tmp = fabs((cu[0] - cl[0]) * (cu[1] - cl[1]));
        avg_area += avg_Tmp;

        if (objcnt % 100 == 0)
            cout << ".";
        if (objcnt % 1000 == 0)
            cout << " " << objcnt << " objects loaded." << endl;
    }
    fin1.close();
    cout << objcnt << " objects are loaded." << endl;

    //accumulate the memory usage in data set
    tmpMemUsage = sizeof(float) * dimen * objcnt;
    MemTree = tmpMemUsage;
    cout << "Memory: " << MemTree / MB << " MB used in storing objects" << endl;
    //


    //-------------------------------------------------------------------------
    // create a Rtree based on TGS
    //-------------------------------------------------------------------------
    cout << "bulkloading R-tree... ";
    const int maxChild =
            (pagesize - RtreeNode::size()) /
            RtreeNodeEntry::size(dimen);            // no. of entries per node
    //MainMemory mem(pagesize);
    FileMemory mem(pagesize, idxname, RtreeNodeEntry::fromMem, true);
    Rtree *rtree =
            TGS::bulkload(mem, dimen, maxChild, maxChild,
                          (int) (maxChild * 0.3), (int) (maxChild * 0.3), p, objcnt, false);
    cout << "[DONE]" << endl;

    avg_area = avg_area / float(objcnt);   //average area of data MBRs
    cout << "ds_x0:" << ds_x0 << " ds_x1:" << ds_x1 << " ds_y0:" << ds_y0 << " ds_y1:" << ds_y1 << endl;
    cout << "average area: " << avg_area << endl;

    //-------------------------------------------------------------------------
    // test 2. find the no. of leaf nodes, their total areas, total perimeters
    //-------------------------------------------------------------------------
    float area = rtree->nodeVolume(0);
    float peri = rtree->nodePerimeter(0);
    int cnt = rtree->nodeCount(0);
    cout << "leaf node cnt,area,perimeter: ";
    cout << cnt << ", " << area << ", " << peri << endl;
    cout << " Max leaf fanout: " << rtree->m_maxLeafChild << endl;
    cout << " Max node fanout: " << rtree->m_maxNodeChild << endl;

    //accumulate the memory usage in R-tree nodes
    tmpMemUsage = pagesize * (rtree->nodeCount(0));
    MemTree = tmpMemUsage;
    cout << "Memory: " << MemTree / MB << " MB used in R-tree" << endl;

//*
    fstream fp_Queries;
    fp_Queries.open(qfilename, ios::in);
    if (!fp_Queries.is_open()) {
        cout << "Error in opening Queries.txt file!" << endl;
        exit(0);
    }
    //*/

// int pos[6];
// findTopKfromDB(PG, objcnt, dimen, pos);  

    FILE *myout;
    myout = fopen("./myout.txt", "w");
    if (myout == nullptr) {
        printf("My file wasn't correctly opened");
        exit(1);
    }

    //srand(time(NULL));
    float NoOfTotalHSsProcessed = 0;
    for (int round = 0; round < Repeat; round++) {

        /*//initialize all the variables: this initialization part is for the case with parameter t
        totalNoOfMinCells=0;
        totalMin_K=0;
        MemUsage=0;
        tmpMemUsage=0;
        TimeCost=0;
        IOs=0;
        Array a_PagesRead;
        totalNoOfBitStringsProcessed=0;
        totalNoOfPrunedBitStrings=0;
        totalNoOfZeroExtentBinStrings=0;
        totalNoOfDiscardedCells=0;
        /*///end of initialization

        //if (round != 3) continue; //use this condition to focus on t=100 as the default case
        //if (round==5) MaxQuadTreeLevels=9;

        //generate random query workload
        float TmpSum = 0;
        float val[DMAX];
        //Generating Method 1: randomly generate a d-dimensional point as the query
        /*
        for (int j=0; j<dimen; j++)
        {
         val[j] = rnd.frandom();
         //val[j] = float(rand())/RAND_MAX;
         TmpSum += val[j];
        }

        for (int j=0; j<dimen; j++)   //normalization
         val[j] = val[j]/TmpSum;

        Point pt(dimen,val);
        cout << "Query point "<< round+1 <<": ";
        for (int j=0; j<dimen; j++)
        {
        //cout << val[j] << " ";
        cout << pt[j] << " ";
        }
        cout << endl;
        getchar();
        //*/

        //Generating Method 2: randomly pick a record from dataset as the query
        //*
        //long int rdIdx =(long int)(1 + (objcnt-1)*(float(rand())/RAND_MAX));
        long int rdIdx;
        //rdIdx=pos[round];
        //cout << "rdIdx=" << rdIdx << endl;

        fp_Queries >> rdIdx;
        if (fp_Queries.eof()) break;;

        for (int j = 0; j < dimen; j++) val[j] = (PG[rdIdx][j] + PG[rdIdx][dimen + j]) / 2.0;
        Point pt(dimen, val);
        //cout << "Query point "<< round+1 <<" (id="<< rdIdx << ",rank="<< posName[round] << "): ";
        cout << "Query point " << round + 1 << " (id=" << rdIdx << "): ";
        for (int j = 0; j < dimen; j++) cout << pt[j] << " ";
        cout << endl;
        //getchar();
        //*/
        //end of query workload generation


        //-------------------------------------------------------------------------
        // main component of maximum rank query
        //-------------------------------------------------------------------------
        int maxStackSize = 0;
        Array a_DiskPage;
        a_DiskPage.clean();
        vector<long>::iterator Int_vIter;
        multimap<float, int>::iterator DblInt_Iter;

        //process MaxRank query, i.e., compute intervals for current query point
        clock_t st, ed;
        double duration;

        int min_k;
        NoOfTotalHalfspacesAdded = 0;
        NoOfEntriesInLists = 0;
        numOfLeavesToIntersect = 0;
        entriesInLists.clear();
        st = clock();
        if (strcmp(methodName, "AA") == 0) {
            cout << "AA Algorithm invoked." << endl;
            if (dimen == 2) min_k = S3::AA_2D(dimen, *rtree, PG, pt, totalNoOfMinCells, maxStackSize, a_DiskPage);
            if (dimen > 2) min_k = S3::AA_HD(dimen, *rtree, PG, pt, totalNoOfMinCells, maxStackSize, a_DiskPage);
        } else if (strcmp(methodName, "BA") == 0) {
            cout << "BA Algorithm invoked." << endl;
            if (dimen == 2) min_k = S3::BA_2D(dimen, *rtree, PG, pt, totalNoOfMinCells, maxStackSize, a_DiskPage);
            if (dimen > 2) min_k = S3::BA_HD(dimen, *rtree, PG, pt, totalNoOfMinCells, maxStackSize, a_DiskPage);
        } else {
            cout << "Invalid method type! Please specify 'AA' for AA algo, or 'BA' for BA algo." << endl;
            return -1;
        }
        ed = clock();
        TimeCost = TimeCost + float(ed - st) * 1000 / CLOCKS_PER_SEC;
        IOs = IOs + a_DiskPage.size();
        cout << "Minimal cell computation finished! Method=" << methodName << endl;

        fprintf(myout, "%ld,%d\n", rdIdx, min_k);

        totalMin_K = totalMin_K + min_k;

        //accumulate the memory usage
        //tmpMemUsage=sizeof(VirtualRNode)*NonResultEntry.size();
        //MemUsage+=tmpMemUsage;

        //accumulate all the pages accessed during each query execution
        for (int tmp = 0; tmp < a_DiskPage.size(); tmp++)
            a_PagesRead.append((void *) a_DiskPage.get(tmp));


        //compute the size of quad tree
        long int tmpTotalMemQuadTree = 0;
        long int tmpMem1 = 0, tmpNoOfNodes;
        float tmpMemRatio = 0;

        tmpNoOfNodes = numOfTreeNodes - numOfInvalidNodes;
        tmpMem1 = tmpNoOfNodes * 2 + NoOfEntriesInLists + NoOfTotalHalfspacesAdded * (dimen + 1) +
                  numOfLeavesToIntersect * 2;
        tmpTotalMemQuadTree = tmpMem1;
        tmpMemRatio = float(tmpMem1 - numOfLeavesToIntersect * 2) / float(tmpMem1);

        if (TotalMemQuadTree < tmpTotalMemQuadTree) {
            TotalMemQuadTree = tmpTotalMemQuadTree;

            maxNoOfQuadTreeNodes = tmpNoOfNodes;
            maxNoOfEntriesInLists = NoOfEntriesInLists;
            maxNoOfTotalHalfspacesAdded = NoOfTotalHalfspacesAdded;
            maxNoOfLeavesToIntersect = numOfLeavesToIntersect;
        }
        NoOfTotalHSsProcessed += NoOfTotalHalfspacesAdded;
        memRatioOfQuadTree += tmpMemRatio;
        //

        cout << "########################################################" << endl;
        cout << "Result of round " << round + 1 << endl;
        cout << "Average #minCells: " << float(totalNoOfMinCells) / (round + 1) << endl;
        cout << "Average minimum order: " << float(totalMin_K) / (round + 1) << endl;
        cout << "Max. #Quad tree nodes: " << maxNoOfQuadTreeNodes << endl;
        cout << "Max. #elements in quad tree node lists: " << maxNoOfEntriesInLists << endl;
        cout << "Max. #halfspaces inserted: " << maxNoOfTotalHalfspacesAdded << endl;
        cout << "Avg. #halfspaces inserted: " << float(NoOfTotalHSsProcessed) / (round + 1) << endl;
        cout << "Max. #leaf nodes to intersect: " << maxNoOfLeavesToIntersect << endl;
        cout << "Avg. time to build quad tree: " << float(timeBuildQuadTree) / (round + 1) << endl;
        cout << "Avg. time to do node intersection: " << float(timeNodeIntersection) / (round + 1) << endl;

        cout << "Max. memory usage: " << float(TotalMemQuadTree * 4) / MB << endl;
        cout << "Max. memory ratio of quadtree: " << float(memRatioOfQuadTree) / (round + 1) << endl;
        cout << "Average #bin-strings processed: " << float(totalNoOfBitStringsProcessed) / (round + 1) << endl;
        cout << "Average #bin-strings pruned: " << float(totalNoOfPrunedBitStrings) / (round + 1) << endl;
        cout << "Average #zero-extent bin-strings: " << float(totalNoOfZeroExtentBinStrings) / (round + 1) << endl;
        cout << "Average #discarded bin-strings: " << float(totalNoOfDiscardedCells) / (round + 1) << endl;
        cout << "Time used in max rank query: " << float(TimeCost) / (round + 1) << " msec" << endl;
        cout << "Overall I/O cost: " << float(IOs) / (round + 1) << " disk pages" << endl;
        cout << endl;

    }//end of repeat

    fclose(myout);

    /*/measuring IO cost when different size (i.e., 1%,2%,5%,10%) of buffers are given
    int tmpIO1=IOMeasure::lru(a_PagesRead,int(0.01*cnt));
    int tmpIO2=IOMeasure::lru(a_PagesRead,int(0.02*cnt));
    int tmpIO5=IOMeasure::lru(a_PagesRead,int(0.05*cnt));
    int tmpIO10=IOMeasure::lru(a_PagesRead,int(0.1*cnt));
    int tmpIO20=IOMeasure::lru(a_PagesRead,int(0.2*cnt));
    AvgIOBufferP1=tmpIO1;
    AvgIOBufferP2=tmpIO2;
    AvgIOBufferP5=tmpIO5;
    AvgIOBufferP10=tmpIO10;
    AvgIOBufferP20=tmpIO20;
    //*/ //end of measuring buffer effect 

    /*
    cout << "########################################################" << endl;    
    cout << "Average performance summary over " << Repeat << " rounds of executation: " << endl;
    cout << "Average #minCells: " << float(totalNoOfMinCells)/Repeat << endl;
    cout << "Average minimum order: " << float(totalMin_K)/Repeat << endl;
    cout << "Average #bin-strings processed: " << float(totalNoOfBitStringsProcessed)/Repeat << endl;
    cout << "Average #bin-strings pruned: " << float(totalNoOfPrunedBitStrings)/Repeat << endl;
    cout << "Average #zero-extent bin-strings: " << float(totalNoOfZeroExtentBinStrings)/Repeat << endl;
    cout << "Average #discarded bin-strings: " << float(totalNoOfDiscardedCells)/Repeat << endl;
    cout << "Time used in max rank query: " << float(TimeCost)/Repeat << " msec"  << endl;
    cout << "Total memory usage: " << (float(MemUsage)/Repeat+MemTree)/MB << " MB" << endl;    
    cout << "Overall I/O cost: " << float(IOs)/Repeat << " disk pages" << endl;
    cout << "I/O cost when buffer sizes are:   1%      2%      5%     10%     20%"<< endl;
    cout << "                              ";
    cout << AvgIOBufferP1/Repeat << "  "<< AvgIOBufferP2/Repeat << "  "<< AvgIOBufferP5/Repeat << "  "<< AvgIOBufferP10/Repeat << "  "<< AvgIOBufferP20/Repeat << endl;
    //*/


    //output the performance result
    fperfOut << filename << " ";
    //fperfOut << dbFile << " ";
    fperfOut << methodName << " ";
    fperfOut << pagesize << " ";
    //fperfOut << posName[3] << " ";  //set focal point to the top-k record w.r.t ranking function x1+x2+...+xd
    fperfOut << MaxQuadTreeLevels << " ";
    fperfOut << dimen << " ";
    fperfOut << "tao=" << tao << " ";
    fperfOut << float(totalNoOfMinCells) / Repeat << " ";       //#mincells
    fperfOut << float(totalMin_K) / Repeat << " ";              //avg. minimum order (k)
    fperfOut << float(totalNoOfBitStringsProcessed) / Repeat << " ";   //avg. total #processed bin-strings
    fperfOut << float(TimeCost) / Repeat << " ";
    fperfOut << (float(MemUsage) / Repeat + MemTree) / MB << " ";
    fperfOut << float(IOs) / Repeat << endl;
    //fperfOut << AvgIOBufferP1/Repeat << " ";
    //fperfOut << AvgIOBufferP2/Repeat << " ";
    //fperfOut << AvgIOBufferP5/Repeat << " ";
    //fperfOut << AvgIOBufferP10/Repeat << " ";
    //fperfOut << AvgIOBufferP20/Repeat << endl;
    fperfOut.close();

//  }//end of repeat   
    //fperfOut << endl;

    //Free allocated memory
    for (long int i = 1; i <= objcnt; i++)
        delete PG[i];
    delete p;
    delete PG;
    cout << "successfully freed the memory allocated!" << endl;
    //}
    //   fperfOut.close();
    // fin0.close();

    //output memory usage
    fstream finMem;

    finMem.open("/proc/meminfo", ios::in);
    if (finMem.is_open()) {
        while (true) {
            string line;
            finMem >> line;
            if (finMem.eof()) break;
            cout << line << endl;
        }
    } else {
        printf("Couldn't retrieve memory info!\n");
    }
    finMem.close();

    finMem.open("/proc/self/status", ios::in);
    if (finMem.is_open()) {
        while (true) {
            string line;
            finMem >> line;
            if (finMem.eof()) break;
            cout << line << endl;
        }
    } else {
        printf("Couldn't retrieve process info!\n");
    }
    finMem.close();

    return 0;
}


