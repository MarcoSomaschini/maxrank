/* ----------------------------------------------------------------------------
    This header file includes class Search declaration.
    It contains search algorithms on Rtree
---------------------------------------------------------------------------- */

#ifndef SEARCH_DEFINED
#define SEARCH_DEFINED

#include "collection.h"


class Rtree;
class Point;
class Hypercube;


class Search
{
public:
	//constructor/distructor
	Search();
	~Search();
	 
	static int window(				// window search at a given point pt
        Rtree& a_rtree,             // returns a result and a result size
        const Point& a_pt, const float* a_l,
        Array& a_result, int& maxStackSize, Array& a_pageaccessed);
    static int range(               // range search at a given point pt
        Rtree& a_rtree,             // returns a result and a result size
        const Point& a_pt, const float a_range,
        Array& a_result, int& maxStackSize, Array& a_pageaccessed);

    static float nn(                // nn search at a given point pt
        Rtree& a_rtree,             // retruns a result and the distance to kNN
        const Point& a_pt, const int a_k,
        Array& a_result, int& maxHeapSize, Array& a_pageaccessed);
    static void dump(               // this draws MBRs of Rtree
        const Rtree& a_rtree,       // in between specified top level and
        const int a_bottomlevel,    // bottom level through PSDraw
        const int a_toplevel,
        const char* epsname);

	//variables
        static int m_idOfMBR;
};

#endif // SEARCH_DEFINED

