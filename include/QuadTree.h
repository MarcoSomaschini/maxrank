#include "Objects.h"
#include "string.h"
#include <string>
#include <vector>
#include <map>
#include <climits>

using namespace std;

class QuadTree{
		
public:
	int IncID;
	int Dimen;
	int MaxLevel;
	int MaxCapacity;
	long int NoOfInValidNodes;
	QuadNode *root;

	float Qspace[2*MAXDIM+1]; 
	vector<string> Comb;   //dimension combinations for generating vertices of MBRs of the subdivisions
	multimap<long int, QuadNode *> minHeap;    //min-heap to explore quadtree leaf nodes so as to perform halfspace intersection
	
	QuadTree(const float mBox[], const int &dim, const int &level, const int &maxLevel, const int &maxCapacity);
	~QuadTree();

        void traversal(QuadNode *node);

	bool readCombinations();
	bool splitNode(QuadNode *node);
	bool checkNodeValidity(const float Qspace[], QuadNode *node, vector<string> &Comb);    //check validity of node with respect to the subdivision; 
        bool insertHalfSpaces(const QuadNode *root);
        bool appendHalfSpace(const long hsID, QuadNode *node);   //append a single halfspace to the Quadtree	
	int  countCoveredHalfSpaces(const Point1 &pt, QuadNode *root, set<long> &ResultHalfSpaces);
        void collectLeafNodes(QuadNode *root, vector<std::pair<long,QuadNode*> > &Leaves, vector<long int> &listOfCoveredHS);
	bool GenHammingHalfSpaces(char *OutFileName, const int Dimen, vector<char>& HammingString, vector<long int>& HalfSpaceIDs, float subDataSpace[]);
        long int naiveInNodeIntersection(vector<std::pair<long,QuadNode*> > &Leaves,vector<set<long int> > &minCellHalfSpaces, vector<vector<char> > &binaryString);   //a naive solution that exhaustively search all combinations of halfspaces (in a leaf node) to intersect
        long int optimizedInNodeIntersection(vector<std::pair<long, QuadNode*> > &Leaves,vector<set<long int> > &minCellHalfSpaces, vector<vector<char> > &binaryString);//optimization of within node intersection
        bool testHalfspacePair(long int HS1, long int IdxHS1,long int HS2, long int IdxHS2,float subDataSpace[],  multimap<int,string>& InValidHammingStr); //test whether two halfspaces are compatible w.r.t Hamming distance 00,01,10,11
	bool genSubdivisions(const float mBox[],vector<vector<float> > &Subdivisions);
};
