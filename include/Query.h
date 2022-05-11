#include "global.h"
#include <map>

#ifndef _Globals
#define DMAX 10	//max dimensionality
#define _Globals
#endif

using namespace std;

class Query {
public:
    Query();

    ~Query();

    int id;
    int k;
    int affected;
    double IR;
    multimap<float, int> score;
    float function[DMAX];

    void initQuery(int i, int n, float *c, int D);
};

