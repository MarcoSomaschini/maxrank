/* ----------------------------------------------------------------------------
---------------------------------------------------------------------------- */

#include "collection.h"
#include "Query.h"
#include <map>
#include <vector>
#include <set>

#define NUMNONSKY 150

//class PartialHull;

class Rtree;
class Point;
class Hypercube;
class VirtualRNode;
struct Pts;

struct QuadNode;

class S3
{
public:
	//constructor/distructor
	S3();
	~S3();

        static int GetSkylines(
        const int dimen,
        Rtree& a_rtree,             // returns a result and a result size
        std::multimap<long int, VirtualRNode*>& NonResultEntry,  //non-result set, storing rtree node which is simply a data entry
        std::vector<long int>& PrunedNodes, //the heap for storing non-result r-tree nodes and remaining entries in H
        std::set<long>& a_skylines,
        float* PG[],
        int& maxStackSize, Array& a_pageaccessed);

        static int window(
        Rtree& a_rtree, 
        float* PG[],
        const Point& a_pt,
        const float* a_l, 
        std::vector<long int>& a_resultID, 
        int& a_maxStackSize, 
        Array& a_pageaccessed);

        static int incomparableWindows(
        Rtree& a_rtree, 
        float* PG[],
        const Point& a_pt,
        const float* a_l, 
        std::multimap<long int, VirtualRNode*>& NonResultEntry,  //non-result set, storing rtree node which is simply a data entry
        std::vector<long int>& a_resultID, 
        int& a_maxStackSize, 
        Array& a_pageaccessed);

        static int incomparableWindowsHD(
        Rtree& a_rtree, 
        float* PG[],
        Hypercube& Dominee_hc,
        Hypercube& Dominator_hc,
        std::multimap<long int, VirtualRNode*>& NonResultEntry,  //non-result set, storing rtree node which is simply a data entry
        std::vector<long int>& a_resultID, 
        int& a_maxStackSize, 
        Array& a_pageaccessed);

        static int getIncomparableRecords(
        Rtree& a_rtree, 
        float* PG[],
        Hypercube& Dominee_hc,
        Hypercube& Dominator_hc,
        //std::multimap<long int, VirtualRNode*>& NonResultEntry,  //non-result set, storing rtree node which is simply a data entry
        std::vector<long int>& a_resultID,    //the index nodes
        int& a_maxStackSize, 
        Array& a_pageaccessed);
        
        static int AA_2D_old(
        const int dimen, //dimensionality
        Rtree& a_rtree,             
        float* PG[],
        Point& a_pt,
        long int &totalNoOfminCells,  
        int& maxStackSize, Array& a_pageaccessed);

        static int AA_2D(
        const int dimen, //dimensionality
        Rtree& a_rtree,             
        float* PG[],
        Point& a_pt,
        long int &totalNoOfminCells,  
        int& maxStackSize, Array& a_pageaccessed);

        static long int AA_HD(
        const int dimen, //dimensionality
        Rtree& a_rtree,             
        float* PG[],
        Point& a_pt,
        long int &totalNoOfminCells,  
        int& maxStackSize, Array& a_pageaccessed);

        static int BA_2D(    //BA algorithm for 2-d case, i.e., FCA algorithm
        const int dimen, //dimensionality
        Rtree& a_rtree,             
        float* PG[],
        Point& a_pt,
        long int &totalNoOfminCells, 
        int& maxStackSize, Array& a_pageaccessed);

        static int BA_2D_old(    //BA algorithm for 2-d case, i.e., FCA algorithm
        const int dimen, //dimensionality
        Rtree& a_rtree,             
        float* PG[],
        Point& a_pt,
        long int &totalNoOfminCells, 
        int& maxStackSize, Array& a_pageaccessed);
                
        static int BA_HD(    //BA algorithm for general dimensionality
        const int dimen, //dimensionality
        Rtree& a_rtree,             
        float* PG[],
        Point& a_pt,
        long int &totalNoOfminCells,  
        int& maxStackSize, Array& a_pageaccessed);       

        static float minDist(float p1[], float p2[], int dimen);
        static float GetScore(float p[], float weight[], int dimen);    //compute the score

        //check whether a point is dominated by skylines
        static bool IsDominatedBy(
        const int dimen,
        const float pt[], 
        vector<long> a_skylines,
        float* PG[]);

        //check whether point pt is dominated by another point pt0
        static bool IsDominatedBy(
        const int dimen,
        const float pt[], 
        const float pt0[]);
        
        static void OutputForQhalf(
        FILE *fout1, 
        const int dimen, 
        Query& a_query, 
        int K, 
        long int* FacetMin, 
        long int* FacetMax,
        long int& FacetMinPID,
        long int& FacetMaxPID,
        float* PG[]);

        static void OutputNonResultSetVanila(
        FILE *fout1, 
        const int dimen, 
        multimap<long,long> NonResultPtID, 
        float* PG[]);

        static void OutputNonResultSet(
        FILE *fout1, 
        const int dimen, 
        multimap<long,long> NonResultPtID, 
        float* PG[]);
};

