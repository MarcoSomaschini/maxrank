//#define CEILING 1
#define MAXDIM 10
#define LEAFNODE -1
#define NONLEAFNODE 1
#define ABOVE 1
#define BELOW 2
#define ON 3
#define OVERLAPPED 4

#define MAXLOOP 1000
#define MAXNOBINSTRINGTOCHECK 60000
#define ZEROEXTENT 1e-15

#include <set>
#include <vector>

struct Point1 {
    float coord[MAXDIM + 1];
};

struct QuadNode {
    long int NodeID;
    float MBR[2 * MAXDIM + 1];   //minimal bounding box
    int Level;
    int Leaf;
    bool isValid;    //flag to indicate whether current node is valid in the reduced query space, i.e., whether we need to store halfspaces in it or expand it

    //QuadNode *Parent;
    QuadNode **Child;

    //std::set<long> coveredHalfspace;
    long int NoOfcoveredHalfspace;
    std::vector<long> intersectedHalfspace;
};
