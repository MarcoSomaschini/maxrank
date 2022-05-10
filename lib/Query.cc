#include "Query.h"
#include "stdlib.h"

using namespace std;

Query::Query()
{}

Query::~Query()
{}

void Query::initQuery(int i, int n, float *c, int D)
{
     id = i;
     k = n;
     memcpy(function, c, D*sizeof(float));
     affected=0;
     IR=1e100;
};

