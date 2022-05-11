#include "QuadTree.h"
#include "stdlib.h"
#include "stdio.h"
#include "math.h"
#include <string>
#include <cstring>
#include <iostream>
#include <fstream>
#include <ostream>
#include <iterator>
#include <map>

#include <random>

using namespace std;

struct QuadNode;
bool NotfirstTimeInsert;
extern int maximumLevel;
extern long numOfTreeNodes;
extern long numOfLeafNodes;
extern long numOfNodeSplits;
extern int normalizedMax;
extern int numOfSubdivisions;
extern long int totalNoOfBitStringsProcessed;
extern long int totalNoOfPrunedBitStrings;
extern long int totalNoOfZeroExtentBinStrings;
extern long int totalNoOfDiscardedCells;
extern map<long int, long int> entriesInLists;

extern bool verbose;

extern int tao;

extern vector<vector<float> > HalfSpaces;
extern long int NoOfNewlyAddedHalfSpaces;
extern map<long int, long int> RdIDtoHalfplaneID;

bool PointCoveredByMBR(const int &Dimen, const float mbr[], const Point1 &pt);

int PointVersusHalfSpace(const int &Dimen, const float hs[],
                         const Point1 &pt);   //position of a point to halfspace: is the point above, below, or on the halfspace
int MbrVersusHalfSpace(const int &Dimen, const float hs[], const float mbr[],
                       vector<string> &Comb);   //position of an MBR to a halfspace: is the point above, below, or intersected by the halfspace?
bool MbrIsValid(const int &Dimen, const float hs[], const float mbr[],
                vector<string> &Comb);   //position of an MBR to a halfspace: is the point above, below, or intersected by the halfspace?
void myitoa(unsigned long val, char *buf, unsigned radix);

void GenBinaryString(long int len1, long int Max,
                     multimap<int, vector<char> > &binString);  //function to generate combinations
void GenLenNBinaryString(long int len1, long int HammingDistance,
                         multimap<int, vector<char> > &binString); //gen binstring with specified Hamming len

QuadTree::QuadTree(const float mBox[], const int &dim, const int &level, const int &maxLevel, const int &maxCapacity) {

    root = new QuadNode;

    Dimen = dim;
    IncID = 0;
    MaxLevel = maxLevel;
    MaxCapacity = maxCapacity;
    NoOfInValidNodes = 0;

    memcpy(Qspace, mBox, (2 * dim) * sizeof(float));
    memcpy(root->MBR, mBox, (2 * dim) * sizeof(float));
    root->NodeID = IncID++;
    root->Level = level;
    root->Leaf = LEAFNODE;
    root->isValid = true;
    //root->Parent = NULL;
    root->Child = NULL;
    root->NoOfcoveredHalfspace = 0;

    //root->coveredHalfspace = NULL;
    //root->intersectedHalfspace = NULL;
    maximumLevel = level;
    numOfTreeNodes = 1;
    numOfLeafNodes = 1;
    numOfNodeSplits = 0;
    NotfirstTimeInsert = false;
}

QuadTree::~QuadTree() {
    //*
    vector<QuadNode *> H;
    vector<QuadNode *>::iterator sItr;
    QuadNode *n;
    H.push_back(root);
    cout << "Begin freeing QT tree!" << endl;
    while (H.size() > 0) {
        n = H.back();
        H.pop_back();
        if (n->Child != NULL) {
            for (int i = 0; i < numOfSubdivisions; i++)
                if (n->Child[i] != NULL)
                    H.push_back(n->Child[i]);
            delete[] n->Child;
        }
        delete n;
    }
    cout << "Successfully freed QT tree!" << endl;
    //*/
}

bool QuadTree::genSubdivisions(const float mBox[], vector<vector<float> > &Subdivisions) {

    float mid[MAXDIM + 1];
    vector<float> subDiv;

    for (int i = 0; i < Dimen; i++)
        mid[i] = (mBox[i] + mBox[i + Dimen]) / 2.0;

    for (int i = 0; i < numOfSubdivisions; i++) {
        vector<float> subDiv;
        for (int j = 0; j < Dimen; j++) {    //compute the starting point on each dimension
            if (Comb[i][j] == '0') {
                subDiv.push_back(mBox[j]);
            } else if (Comb[i][j] == '1') {
                subDiv.push_back(mid[j]);
            }
        }
        for (int j = 0; j < Dimen; j++) {    //compute the ending point on each dimension
            if (subDiv[j] == mBox[j]) {
                subDiv.push_back(mid[j]);
            } else if (subDiv[j] == mid[j]) {
                subDiv.push_back(mBox[j + Dimen]);
            }
        }
        Subdivisions.push_back(subDiv);
    }
    /*
    cout << "begin printing subdivisions..." << endl;
    for (int i=0;i<Subdivisions.size();i++){
         cout << "[";
         for (int j=0;j<Dimen;j++){
             cout << Subdivisions[i][j] << " ";
         }
         cout << "],";
         cout << "[";
         for (int j=0;j<Dimen;j++){
             cout << Subdivisions[i][Dimen+j] << " ";
         }
         cout << "]" << endl;
    }
        getchar();
    //*/

    return true;
}

bool QuadTree::checkNodeValidity(const float Qspace[], QuadNode *node,
                                 vector<string> &Comb) {    //check validity of node with respect to the subdivision;

    int numOfValidVertices = 0;

    int numOfVertices = Comb.size();

    for (int i = 0; i < numOfVertices; i++) {
        Point1 pt;

        for (int j = 0; j < Dimen; j++) {    //to enumerate each of the 2^d vertices
            if (Comb[i][j] == '0')
                pt.coord[j] = node->MBR[j];
            if (Comb[i][j] == '1')
                pt.coord[j] = node->MBR[Dimen + j];
        }

        float sum = 0;
        for (int k = 0; k < Dimen; k++) {
            sum = sum + pt.coord[k];
        }

        if (sum <= normalizedMax) numOfValidVertices++;               //the normalized maximum value is 1
    }

    if (numOfValidVertices > 1)
        return true;
    else
        return false;
}

bool QuadTree::splitNode(QuadNode *node) {

    if (node->Level == MaxLevel)
        return false;   //if the level of current node is exceeding the maximal level, then stop node-splitting

    if (node->isValid == false) return false;  //current node lies outside q_1+q_2+...q_{d-1}<=1, so it is discarded

    vector<vector<float> > subDiv;
    genSubdivisions(node->MBR, subDiv);    //generate 2^d subdivisions of the MBR of the node that is to split

    node->Child = new QuadNode *[numOfSubdivisions];
    for (int i = 0; i < numOfSubdivisions; i++) {
        node->Child[i] = new QuadNode;

        node->Child[i]->NoOfcoveredHalfspace = 0;
        node->Child[i]->NodeID = IncID++;
        node->Child[i]->Level = node->Level + 1;
        node->Child[i]->Leaf = LEAFNODE;
        memcpy(node->Child[i]->MBR, &(subDiv[i][0]), (2 * Dimen) * sizeof(float));

        if (checkNodeValidity(Qspace, node->Child[i], Comb) ==
            true)  //prune away nodes lying outside q_1+q_2+...q_{d-1}<=1
            node->Child[i]->isValid = true;
        else {
            node->Child[i]->isValid = false;
            NoOfInValidNodes++;
            //continue;
        }

        //node->Child[i]->Parent = node;
        node->Child[i]->Child = NULL;
        node->Leaf = NONLEAFNODE;
    }
    if (maximumLevel < (node->Level + 1)) maximumLevel = (node->Level + 1);
    numOfTreeNodes = numOfTreeNodes + numOfSubdivisions;
    numOfLeafNodes = numOfLeafNodes + numOfSubdivisions -
                     1;    //aggregate the number of leaf nodes whenever a node-split happened

    /*
    cout << "Subspace partitions:" << endl;
    for (int i=0;i<numOfSubdivisions;i++){
         cout << "[";
         for (int j=0;j<Dimen;j++){
             cout << node->Child[i]->MBR[j] << " ";
         }
         cout << "],";
         cout << "[";
         for (int j=0;j<Dimen;j++){
             cout << node->Child[i]->MBR[Dimen+j] << " ";
         }
         cout << "]" << endl;
    }
    //*/

    numOfNodeSplits++;

    return true;
}

bool QuadTree::appendHalfSpace(const long hsID, QuadNode *node) {   //append a single halfspace to Quadtree

    QuadNode *n = node;

    int pos = MbrVersusHalfSpace(Dimen, &(HalfSpaces[hsID][0]), n->MBR, Comb);
    if (pos == BELOW) {    //halfspace lies above current node, then store it into its list
        //(n->coveredHalfspace).insert(hsID);
        (n->NoOfcoveredHalfspace)++;
        return true;
    } else if (pos == OVERLAPPED) {
        if (n->Leaf ==
            LEAFNODE) {   //current node n is a leaf, so we need to insert hs to n's intersectedHS set and then check the capacity of the set
            (n->intersectedHalfspace).push_back(hsID);
            if ((n->intersectedHalfspace).size() >
                MaxCapacity) {   //after insertion, max capacity of the set intersectedHalfspace exceeded, we need to split current node
                if (n->Level < MaxLevel) {
                    splitNode(n);
                    for (vector<long>::iterator itr = (n->intersectedHalfspace).begin();
                         itr != (n->intersectedHalfspace).end(); itr++) {
                        for (int j = 0; j < numOfSubdivisions; j++)
                            if (n->Child[j]->isValid)
                                appendHalfSpace((*itr), n->Child[j]);
                    }
                    (n->intersectedHalfspace).clear();
                }
            }
        } else {
            for (int i = 0; i < numOfSubdivisions; i++)
                if (n->Child[i]->isValid)
                    appendHalfSpace(hsID, n->Child[i]);
        }
    }
    return true;
}

bool QuadTree::insertHalfSpaces(const QuadNode *root) {
    map<long int, long int>::iterator IntInt_mItr;
    set<long int>::iterator sItr;

    if (HalfSpaces.size() == 0) return false;

    long int NoOfExistedHS = 0;
    int idx = HalfSpaces[0].size();

    for (int i = 0; i < HalfSpaces.size(); i++) {        //insert each halfspace
        //
        if (RdIDtoHalfplaneID.size() > 0 && NotfirstTimeInsert) {
            //cout << "id=" << HalfSpaces[i][idx-1] << endl;
            IntInt_mItr = RdIDtoHalfplaneID.find(HalfSpaces[i][idx - 1]);
            if (IntInt_mItr != RdIDtoHalfplaneID.end()) //current halfplane has already been inserted into the Quadtree
            {
                NoOfExistedHS++;
                continue;
            }
        }
        //*/

        for (int j = 0; j < numOfSubdivisions; j++) {                     //to each of the 2^d children of root
            if (root->Child[j]->isValid)
                appendHalfSpace(i, root->Child[j]);
        }
    }
    NotfirstTimeInsert = true;

    /*
    if (NoOfExistedHS == HalfSpaces.size())
        return false;
    else
       cout << HalfSpaces.size() << " halfspaces have been inserted!" << endl;
   //*/
    return true;
}

/*
int QuadTree::countCoveredHalfSpaces(const Point1 &pt, QuadNode *root, set<long> &ResultHalfSpaces){
	
	QuadNode *n = root;
	long int numOfCoveredHalfSpaces = 0;
	int ChildID;

	while (true){
		if (n->Child == NULL) break;

		for (int i=0;i<numOfSubdivisions; i++){
			if (PointCoveredByMBR(Dimen,(n->Child[i])->MBR,pt)) {
				ChildID = i;
			        numOfCoveredHalfSpaces = numOfCoveredHalfSpaces + ((n->Child[i])->coveredHalfspace).size();
				copy(((n->Child[i])->coveredHalfspace).begin(),((n->Child[i])->coveredHalfspace).end(),std::inserter(ResultHalfSpaces,ResultHalfSpaces.begin()));
				break;
			}
		}
	    n = n->Child[ChildID];
	}

	//process the set of intersectedHalfSpace of a leaf node that covers pt
	////Note (according to Prof Mouratidis' comment in email): do not count in the number of halfspaces that overlap with the leaf node
	//for (set<long>::iterator itr = (n->intersectedHalfspace)->begin(); itr != (n->intersectedHalfspace)->end(); itr++){		 
	//	 int position = PointVersusHalfSpace(Dimen, &(HalfSpaces[(*itr)][0]),pt);
	//	 if (position == BELOW || position == ON){
        //     numOfCoveredHalfSpaces++;
        //     ResultHalfSpaces.insert(*itr);
	//	 }
	//}
	///
	
    return numOfCoveredHalfSpaces;
}
//*/

void QuadTree::traversal(QuadNode *node) {
    if (node == NULL) return;

    cout << "Node id: " << node->NodeID;
    cout << ", #coveredHS: " << node->NoOfcoveredHalfspace << endl;

    if (node->Child != NULL) {
        for (int i = 0; i < numOfSubdivisions; i++) {
            if (node->Child[i]->isValid)
                traversal(node->Child[i]);
        }
    }

}

void QuadTree::collectLeafNodes(QuadNode *root, vector<std::pair<long, QuadNode *> > &Leaves,
                                vector<long int> &numOfCoveredHS) {

    if (root == NULL) return;

    numOfCoveredHS.push_back(root->NoOfcoveredHalfspace);

    entriesInLists.insert(std::make_pair(root->NodeID, root->NoOfcoveredHalfspace));  //measuring memory cost

    if (root->Child != NULL) {
        for (int i = 0; i < numOfSubdivisions; i++) {
            if (root->Child[i]->isValid)
                collectLeafNodes(root->Child[i], Leaves, numOfCoveredHS);
        }
    }

    if (root->Child == NULL) {    //current node is a leaf node
        long int NoOfTotalCoveredHalfSpaces = 0;
        for (vector<long int>::iterator itr = numOfCoveredHS.begin(); itr != numOfCoveredHS.end(); itr++) {
            NoOfTotalCoveredHalfSpaces = NoOfTotalCoveredHalfSpaces + (*itr);
        }

        //Leaves.insert(std::pair<long int, QuadNode *>(NoOfTotalCoveredHalfSpaces,root));
        Leaves.push_back(std::make_pair(NoOfTotalCoveredHalfSpaces, root));

        entriesInLists.insert(std::make_pair(root->NodeID, root->intersectedHalfspace.size()));  //measuring memory cost
    }
    if (numOfCoveredHS.size() > 0) numOfCoveredHS.pop_back();
}

//
bool QuadTree::GenHammingHalfSpaces(char *OutFileName, const int Dimen, vector<char> &HammingString,
                                    vector<long int> &HalfSpaceIDs, float subDataSpace[]) {
    int NoOfHyperplanes = 0;
    long int NoOfHalfSpaces = HalfSpaceIDs.size();

    bool interiorPtExists = true;
    float InteriorPt[MAXDIM];

    NoOfHyperplanes = 2 * Dimen + NoOfHalfSpaces;   //total number of halfspaces

    //randomization process to generate an interior point
    long int loops = 0;

    //for (int i=0;i<Dimen;i++)
    //     cout << "[" << subDataSpace[i] << "," << subDataSpace[Dimen+i] << "]"<< endl;

    while (true) {
        loops++;
        //cout << "loop " << loops << endl;
        if (loops >=
            MAXLOOP) //there is a very low probability that an interior point exists for current set of halfspaces
        {
            //cout << "exceed the maximal loop limits, none interior point exists!" << endl;
            interiorPtExists = false;
            break;
        }

        // The original was using rand(), which has a bad randomness
        // Also directly generate numbers between 0 and 1
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> dis(0, 1);

        //Discard query points over y = 1 - x (in 3D), as the last parameter (yield by 1 - x - y) would be negative
        //By doing so eventual mincells in such regions will not be found
        while(true){
            for (int i = 0; i < Dimen; i++)
                InteriorPt[i] = subDataSpace[i] + (subDataSpace[Dimen + i] - subDataSpace[i]) * dis(gen);

            float sum = 0;
            for (int i = 0; i < Dimen - 1; i++)
                sum += InteriorPt[i];

            if (InteriorPt[Dimen - 1] < 1 - sum)
                break;
        }

        /*
        for (int i = 0; i < Dimen; i++)
            InteriorPt[i] = subDataSpace[i] + (subDataSpace[Dimen + i] - subDataSpace[i]) * (float(rand()) / RAND_MAX);
        */

        int index = 0;
        long int count = 0;
        for (vector<long int>::iterator sItr = HalfSpaceIDs.begin(); sItr != HalfSpaceIDs.end(); sItr++) {
            if (HammingString[index] == '0')    //for the case where ax_1+bx_2+... <= d    1
            {
                long int hsID = *sItr;
                float sum = 0;
                for (int i = 0; i < Dimen; i++)
                    sum = sum + HalfSpaces[hsID][i] * InteriorPt[i];
                if (sum <= HalfSpaces[hsID][Dimen]) count++;
            } else if (HammingString[index] == '1')   //for the case where ax_1+bx_2+... > d   0
            {
                long int hsID = *sItr;
                float sum = 0;
                for (int i = 0; i < Dimen; i++)
                    sum = sum + HalfSpaces[hsID][i] * InteriorPt[i];
                if (sum >= HalfSpaces[hsID][Dimen]) count++;
            }
            index++;
        }
        if (count == HalfSpaceIDs.size()) {
            //cout << "We have found an interior point!" << endl;
            //for (int i = 0; i < Dimen; i++) cout << InteriorPt[i] << " ";
            //cout << endl;
            break;
        }
    }
    if (!interiorPtExists) {
        //cout << "Oops, it's unlikely that an interior point exists!!" << endl;
        return false;
    }
    //

    FILE *fout1 = fopen(OutFileName, "w");
    if (fout1 == NULL) {
        cout << "Error in opening file " << OutFileName << endl;
        getchar();
        exit(0);
    }

    fprintf(fout1, "%d 1\n", Dimen);   //the dimensionality
    for (int i = 0; i < Dimen; i++)
        fprintf(fout1, "%f ",
                InteriorPt[i]); //the feasible point found by Monte Carlo process above   *****************
    fprintf(fout1, "\n");
    fprintf(fout1, "%d\n", Dimen + 1);   //dimensionality + 1
    fprintf(fout1, "%d\n", NoOfHyperplanes);  //the total number of hyperplanes for intersection

    //output the 2*dimen bounding facets of the hypercube of the sub-dataspace [x_min,y_min,z_min,...] [x_max,y_max,z_max,...]
    int Cooef = -1;
    for (int i = 1; i <= 2; i++) {
        if (i == 2)
            Cooef = -Cooef;
        for (int j = 1; j <= Dimen; j++) {
            for (int m = 1; m <= Dimen; m++)
                (m == j) ? fprintf(fout1, "%d  ", Cooef) : fprintf(fout1, "0  ");
            if (i == 2)
                fprintf(fout1, "%f\n", -subDataSpace[Dimen + j - 1]);
            else
                fprintf(fout1, "%f\n", subDataSpace[j - 1]);
        }
    }
    ///

    int index = 0;
    for (vector<long int>::iterator sItr = HalfSpaceIDs.begin(); sItr != HalfSpaceIDs.end(); sItr++) {
        if (HammingString[index] == '0')    //for the case where ax_1+bx_2+... <= d  1
        {
            long int hsID = *sItr;
            for (int i = 0; i < Dimen; i++)
                fprintf(fout1, "%f ", HalfSpaces[hsID][i]);
            fprintf(fout1, "%f\n", -HalfSpaces[hsID][Dimen]);   //the offset
        } else if (HammingString[index] == '1')   //for the case where ax_1+bx_2+... > d   0
        {
            long int hsID = *sItr;
            for (int i = 0; i < Dimen; i++)
                fprintf(fout1, "%f ", -HalfSpaces[hsID][i]);
            fprintf(fout1, "%f\n", HalfSpaces[hsID][Dimen]);   //the offset
        }
        index++;
    }
    fclose(fout1);
    //cout << "output data for halfspace intersection finished!" << endl;

    return true;
}
//*/

long int QuadTree::naiveInNodeIntersection(vector<std::pair<long, QuadNode *> > &Leaves,
                                           vector<set<long int> > &minCellHalfSpaces,
                                           vector<vector<char> > &binaryString)   //a naive solution that exhaustively search all combinations of halfspaces (in a leaf node) to intersect
{
    multimap<long int, QuadNode *> nodesToIntersect;   //store the nodes (in ascending order) by using the number of their intersected halfspaces
    //multimap<long int, QuadNode *>::iterator itr, itr1;
    vector<std::pair<long, QuadNode *> >::iterator itr, itr1;

    multimap<int, string>::iterator msItr;

    map<long int, long int>::iterator IntInt_mItr;

    set<string>::iterator ssItr, ssItr1;
    set<long int>::iterator sItr, sItr1;

    vector<string> FilesToRemove;

    long int NoOfBitStringsProcessed = 0;

    if (Leaves.empty()) {
        cout << "There is no leaf nodes to perform intersection!" << endl;
        return -1;
    }

    FILE *fp_tmpIn;

    char Buf[1024];
    char *token;
    char m_seperator[] = " :\n\t";

    char volumeFilename[2048] = "Vol";
    myitoa(Dimen, Buf, 10);
    strcat(volumeFilename, Buf);
    myitoa(rand(), Buf, 10);
    strcat(volumeFilename, Buf);
    strcat(volumeFilename, "D.txt");

    char namePrefix[] = "../tmp/HalfSpaces";
    char nameSuffix[] = ".txt";
    char halfspaceFileName[1024];

    long int minOrder = INT_MAX;

    long int NoOfCoveredHS, NoOfIntersectedHS;
    bool NotFoundAllMinCells = true;

    long int NoOfInvalidLeaves = 0;
    long int NoOfNodesLeft = Leaves.size();

    long int NoOfHalfSpacesInNode = 0;
    long int NoOfZeroExtentBinStrings = 0;
    long int NoOfDiscardedCells = 0;

    for (itr = Leaves.begin(); itr != Leaves.end();)    // && NotFoundAllMinCells;)
    {

        //prune away leaf nodes that lie about hyperplane q_1+q2+...+q_d < 1;
        float queryPlane[MAXDIM];
        for (int i = 0; i < Dimen + 1; i++) queryPlane[i] = 1;
        bool isValid = MbrIsValid(Dimen, queryPlane, (*itr).second->MBR, Comb);
        if (!isValid) {
            if (verbose)
                cout << endl << "Leaf node " << (*itr).second->NodeID << " is pruned!" << endl;
            NoOfInvalidLeaves++;
            ++itr;
            continue;
        }

        NoOfCoveredHS = (*itr).first;

        //modified this line to incorporate iMaxRank query parameter \tao
        if (NoOfCoveredHS - tao >
            minOrder) //terminate searching min-cells, cause no cell with smaller order exists any more
        {
            cout << "Search terminated! No node can contain cells with order less than " << minOrder << endl;
            break;
        }

        NoOfIntersectedHS = ((*itr).second)->intersectedHalfspace.size();
        if (NoOfIntersectedHS > 0)
            nodesToIntersect.insert(std::pair<long int, QuadNode *>(NoOfIntersectedHS, (*itr).second));
        else {
            itr++;
            continue;
        }

        itr1 = itr;
        itr1++;
        //nodesToIntersect.clear();
        while (true) {
            if (itr1 == Leaves.end() || NoOfCoveredHS != (*itr1).first) break;

            NoOfIntersectedHS = 0;
            NoOfIntersectedHS = ((*itr1).second)->intersectedHalfspace.size();
            if (NoOfIntersectedHS > 0)
                nodesToIntersect.insert(std::pair<long int, QuadNode *>(NoOfIntersectedHS, (*itr1).second));
            itr1++;
        }
        if (itr1 == Leaves.end() && nodesToIntersect.size() == 0) break;
        itr = itr1;

        //perform in-node halfspace intersection for nodes sorting in ascending order according to the number of intersected halfspaces
        for (multimap<long, QuadNode *>::iterator itr1 = nodesToIntersect.begin();
             itr1 != nodesToIntersect.end(); itr1++) {
            if (verbose)
                cout << endl << "Processing Leaf node " << (*itr1).second->NodeID << endl;

            NoOfHalfSpacesInNode = itr1->first;

            NoOfNodesLeft--;
            if (NoOfHalfSpacesInNode < 1) {
                if (verbose)
                    cout << "The node do not contain any halfspaces, skipping..." << endl;
                continue;
            }
            if (NoOfNodesLeft % 1000 == 0)
                cout << "Number of nodes left to process :" << NoOfNodesLeft << endl;
            //intersect the halfspaces in each leaf node
            long int NoOfCombinations = (int) pow(2.0, NoOfHalfSpacesInNode);
            multimap<int, vector<char> > binString;
            long int HammingDistance = 0;
            long int LoopCounter = 0;
            bool stopIncrHammingDist = false;
            long int HSidx1, HSidx2;
            size_t substrPos1, substrPos2;
            while ((HammingDistance <= NoOfHalfSpacesInNode) && !stopIncrHammingDist) {
                //LoopCounter++;

                if (minOrder < INT_MAX && (HammingDistance + NoOfCoveredHS - tao) > minOrder)
                    break;

                //if (LoopCounter >= MAXNOBINSTRINGTOCHECK) break;

                binString.clear();
                GenLenNBinaryString(NoOfHalfSpacesInNode, HammingDistance, binString);  //generate all the combinations

                //getchar();
                multimap<int, vector<char> >::iterator itrHamming;
                for (itrHamming = binString.begin(); itrHamming != binString.end(); itrHamming++) {
                    LoopCounter++;
                    if (LoopCounter >= MAXNOBINSTRINGTOCHECK) {
                        stopIncrHammingDist = true;
                        break;
                    }

                    NoOfBitStringsProcessed++;
                    /* Too much blob
                    if (verbose)
                        cout << "Hamming Dist=" << HammingDistance << ", Mininmal Order=" << minOrder << ", #coveredHS="
                             << NoOfCoveredHS << ", #intersectedHS=" << NoOfHalfSpacesInNode << endl;
                    if (verbose)
                        cout << "Testing Hamming string: "
                             << string(itrHamming->second.begin(), itrHamming->second.end()) << endl;*/

                    if (HammingDistance != itrHamming->first) {
                        if (verbose)
                            cout << "Intersecting halfspaces with hamming distance = " << itrHamming->first << endl;
                        HammingDistance = itrHamming->first;
                    }

                    //generate halfspaces for intersection based on hamming distance
                    char tmpBuf[64];
                    myitoa(itr1->second->NodeID, tmpBuf, 10);
                    strcpy(halfspaceFileName, namePrefix);
                    strcat(halfspaceFileName, tmpBuf);
                    strcat(halfspaceFileName, "_");
                    //strcat(halfspaceFileName,itrHamming->second);
                    //string tmpStr=string(itrHamming->second.begin(),itrHamming->second.end());
                    myitoa(rand(), tmpBuf, 10);
                    strcat(halfspaceFileName, tmpBuf);
                    strcat(halfspaceFileName, nameSuffix);
                    //cout << "Hamming string will be written to " << halfspaceFileName << endl;
                    bool nonZeroExtent = GenHammingHalfSpaces(halfspaceFileName, Dimen, itrHamming->second,
                                                              itr1->second->intersectedHalfspace, itr1->second->MBR);
                    if (nonZeroExtent == false) {
                        //cout << "Discard Hamming binstring for zero-extent! " << endl;
                        NoOfZeroExtentBinStrings++;
                        continue;
                    }

                    //perform intersection for all the halfspaces inside the set 'intersected halfspaces' of a leaf node
                    char sys_string[4096] = "qhalf Fp < ";
                    strcat(sys_string, halfspaceFileName);

                    strcat(sys_string, " | qconvex FA > ");
                    strcat(sys_string, volumeFilename);
                    system(sys_string);
                    if (verbose)
                        cout << "Command: " << sys_string << " performed..." << endl;

                    //open the file to read the volumes of the intersection of the halfspaces
                    ifstream fp_in(volumeFilename, ios::in);
                    string fileLine, volText;
                    size_t substrPos;
                    if (verbose)
                        cout << "Read and output Vol.txt:" << endl;
                    while (getline(fp_in, fileLine)) {
                        substrPos = fileLine.find("volume:");
                        if (substrPos != std::string::npos)
                            volText = fileLine.substr(substrPos + 7);
                    }
                    fp_in.close();

                    float volume = atof(volText.c_str());
                    if (verbose)
                        cout << "Volume=" << volume << endl;

                    string tmpString = halfspaceFileName;
                    FilesToRemove.push_back(tmpString);  //collect the files for removal later
                    FilesToRemove.push_back(volumeFilename);  //collect the files for removal later


                    if (volume > ZEROEXTENT) //ZEROEXTENT: 1e-10
                    {
                        //store the min-cell found
                        if (minOrder >= (HammingDistance + NoOfCoveredHS - tao)) {
                            if (minOrder == INT_MAX)
                                minOrder = HammingDistance + NoOfCoveredHS;
                            if (minOrder > (HammingDistance + NoOfCoveredHS - tao)) {
                                if (verbose)
                                    cout << "Found a cell with lower overall order!" << endl;
                                minOrder = HammingDistance + NoOfCoveredHS;
                                minCellHalfSpaces.clear();
                                binaryString.clear();
                            }
                            set<long int> tmpSet;
                            std::copy((itr1->second->intersectedHalfspace).begin(),
                                      (itr1->second->intersectedHalfspace).end(),
                                      std::inserter(tmpSet, tmpSet.begin()));
                            minCellHalfSpaces.push_back(tmpSet);
                            vector<char> tmpVec = itrHamming->second;
                            binaryString.push_back(tmpVec);
                            if (verbose)
                                cout << "order=" << HammingDistance << ", #coveredHS=" << NoOfCoveredHS << endl;
                            if (verbose)
                                cout << "Node: " << itr1->second->NodeID << ", found a min-cell with binstring:"
                                     << string(itrHamming->second.begin(), itrHamming->second.end()) << endl;
                        } else {
                            NoOfDiscardedCells++;
                            stopIncrHammingDist = true;
                            string tmpString = halfspaceFileName;
                            FilesToRemove.push_back(tmpString);  //collect the files for removal
                            break;
                        }
                    } else {
                        NoOfZeroExtentBinStrings++;
                        if (verbose)
                            cout << "Node: " << itr1->second->NodeID << ", HammingDist: " << HammingDistance
                                 << " is empty!" << endl;
                    }
                    if (verbose)
                        cout << "#min-cells found so far:" << minCellHalfSpaces.size() << endl;
                }
                if (!stopIncrHammingDist) HammingDistance++;
                if (minOrder < INT_MAX && minOrder < (HammingDistance + NoOfCoveredHS - tao)) {
                    if (verbose)
                        cout << "There aren't any other mincells within this node, proceding..." << endl;
                    break;
                }
                else {
                    if (verbose)
                        cout << "Next, process bit-string with Hamming dist: " << HammingDistance << endl;
                    //getchar();
                }
            }

            if (FilesToRemove.size() >= 100) {
                vector<string>::iterator vItr = FilesToRemove.begin();
                while (!FilesToRemove.empty()) {
                    remove((*vItr).c_str());  //delete the halfspace data file
                    FilesToRemove.erase(vItr);
                    vItr = FilesToRemove.begin();
                }
            }
        }
        //if (minCellHalfSpaces.size()>0) NotFoundAllMinCells=false;
    }
    vector<string>::iterator vItr = FilesToRemove.begin();
    while (!FilesToRemove.empty()) {
        remove((*vItr).c_str());  //delete the halfspace data file
        FilesToRemove.erase(vItr);
        vItr = FilesToRemove.begin();
    }

    //cout << "Number of invalid leaf nodes = " << NoOfInvalidLeaves << endl;
    if (verbose) cout << "#min-Cells found: " << minCellHalfSpaces.size() << ", minOrder=" << minOrder << endl;

    totalNoOfBitStringsProcessed = totalNoOfBitStringsProcessed + NoOfBitStringsProcessed;
    totalNoOfZeroExtentBinStrings = totalNoOfZeroExtentBinStrings + NoOfZeroExtentBinStrings;
    totalNoOfDiscardedCells = totalNoOfDiscardedCells + NoOfDiscardedCells;

    if (minOrder == INT_MAX) return NoOfCoveredHS;
    return minOrder;
}

long int QuadTree::optimizedInNodeIntersection(vector<std::pair<long, QuadNode *> > &Leaves,
                                               vector<set<long int> > &minCellHalfSpaces,
                                               vector<vector<char> > &binaryString) //optimization of within node intersection
{
    multimap<long int, QuadNode *> nodesToIntersect;   //store the nodes (in ascending order) by using the number of their intersected halfspaces
    //multimap<long int, QuadNode *>::iterator itr, itr1;
    vector<std::pair<long, QuadNode *> >::iterator itr, itr1;

    multimap<int, string>::iterator msItr;

    map<long int, long int>::iterator IntInt_mItr;

    set<string>::iterator ssItr, ssItr1;
    set<long int>::iterator sItr, sItr1;

    vector<string> FilesToRemove;

    if (Leaves.empty()) {
        cout << "There is no leaf nodes to perform intersection!" << endl;
        return -1;
    }

    FILE *fp_tmpIn;

    char Buf[1024];
    char *token;
    char m_seperator[] = " :\n\t";

    char volumeFilename[2048] = "Vol";
    myitoa(Dimen, Buf, 10);
    strcat(volumeFilename, Buf);
    myitoa(rand(), Buf, 10);
    strcat(volumeFilename, Buf);
    strcat(volumeFilename, "D.txt");

    char namePrefix[] = "./tmp/HalfSpaces";
    char nameSuffix[] = ".txt";
    char halfspaceFileName[1024];

    long int minOrder = INT_MAX;

    long int NoOfCoveredHS, NoOfIntersectedHS;
    bool NotFoundAllMinCells = true;

    long int NoOfInvalidLeaves = 0;
    long int NoOfNodesLeft = Leaves.size();

    long int NoOfBitStringsProcessed = 0;
    long int NoOfHalfSpacesInNode;
    long int NoOfPrunedBitStrings = 0;
    long int NoOfZeroExtentBinStrings = 0;
    long int NoOfDiscardedCells = 0;

    for (itr = Leaves.begin(); itr != Leaves.end();)    // && NotFoundAllMinCells;)
    {

        //prune away leaf nodes that lie about hyperplane q_1+q2+...+q_d < 1;
        float queryPlane[MAXDIM];
        for (int i = 0; i < Dimen + 1; i++) queryPlane[i] = 1;
        bool isValid = MbrIsValid(Dimen, queryPlane, (*itr).second->MBR, Comb);
        if (!isValid) {
            //cout << "Leaf node " << (*itr).second->NodeID << " is pruned!" << endl;
            NoOfInvalidLeaves++;
            ++itr;
            continue;
        }

        NoOfCoveredHS = (*itr).first;
        if (NoOfCoveredHS > minOrder)
            break; //terminate searching min-cells, cause no cell with smaller order exists any more

        NoOfIntersectedHS = 0;
        NoOfIntersectedHS = ((*itr).second)->intersectedHalfspace.size();
        if (NoOfIntersectedHS > 0)
            nodesToIntersect.insert(std::pair<long int, QuadNode *>(NoOfIntersectedHS, (*itr).second));
        else {
            itr++;
            continue;
        }
        itr1 = itr;
        itr1++;
        nodesToIntersect.clear();
        while (true) {
            if (itr1 == Leaves.end() || NoOfCoveredHS != (*itr1).first) break;

            NoOfIntersectedHS = 0;
            NoOfIntersectedHS = ((*itr1).second)->intersectedHalfspace.size();
            if (NoOfIntersectedHS > 0)
                nodesToIntersect.insert(std::pair<long int, QuadNode *>(NoOfIntersectedHS, (*itr1).second));
            itr1++;
        }
        if (itr1 == Leaves.end()) break;
        itr = itr1;

        //cout << "Processing node " << itr1->second->NodeID << endl;

        //perform in-node halfspace intersection for nodes sorting in ascending order according to the number of intersected halfspaces
        //cout << "Size of nodesToIntersect: " << nodesToIntersect.size() << endl;
        for (multimap<long, QuadNode *>::iterator itr1 = nodesToIntersect.begin();
             itr1 != nodesToIntersect.end(); itr1++) {


            NoOfHalfSpacesInNode = itr1->first;

            NoOfNodesLeft--;
            if (NoOfHalfSpacesInNode <= 1) continue;
            //cout << "Number of nodes left to process :" << NoOfNodesLeft << endl;
            //cout <<"examining node " << itr1->second->NodeID << ", #inter.HS:" << itr1->second->intersectedHalfspace->size() << endl;

            //test compatibility of Hamming distance for each pair of halfspaces
            multimap<int, string> InValidHammingStr;
            long int idx1 = 0, idx2;
            for (vector<long>::iterator sItr = (itr1->second->intersectedHalfspace).begin();
                 sItr != (itr1->second->intersectedHalfspace).end(); sItr++) {
                long int hs1 = (*sItr);
                vector<long>::iterator sItr1;
                sItr1 = sItr;
                sItr1++;
                idx2 = idx1 + 1;
                for (; sItr1 != (itr1->second->intersectedHalfspace).end(); sItr1++) {
                    long int hs2 = (*sItr1);
                    testHalfspacePair(hs1, idx1, hs2, idx2, itr1->second->MBR, InValidHammingStr);
                    idx2++;
                }
                idx1++;
            }
            /*    test output of the incompatible Halfspace pair
            size_t substrPos1,substrPos2;
            for (msItr=InValidHammingStr.begin();msItr!=InValidHammingStr.end();msItr++)
            {
                 cout << "Hamming Distance: " << msItr->first << ", String=" << msItr->second << endl;
                 substrPos1=msItr->second.find("|");
                 cout << msItr->second.substr(0,substrPos1) << " ";
                 substrPos2=msItr->second.find("|",substrPos1+1);
                 cout << msItr->second.substr(substrPos1+1,substrPos2-substrPos1-1) << " ";
                 cout << msItr->second.substr(substrPos2+1) << endl;
            }
            //*/
            //end of testing compatibility

            //intersect the halfspaces in each leaf node
            long int NoOfCombinations = (int) pow(2.0, NoOfHalfSpacesInNode);
            multimap<int, vector<char> > binString;
            long int HammingDistance = 0;
            long int LoopCounter = 0;
            bool stopIncrHammingDist = false;
            long int HSidx1, HSidx2;
            size_t substrPos1, substrPos2;
            while ((HammingDistance <= NoOfHalfSpacesInNode) && !stopIncrHammingDist) {
                //if (LoopCounter >= MAXNOBINSTRINGTOCHECK) break;
                if (minOrder < INT_MAX && (HammingDistance + NoOfCoveredHS) > minOrder)
                    break;

                binString.clear();
                GenLenNBinaryString(NoOfHalfSpacesInNode, HammingDistance, binString);  //generate all the combinations
                multimap<int, vector<char> >::iterator itrHamming;
                for (itrHamming = binString.begin(); itrHamming != binString.end(); itrHamming++) {
                    LoopCounter++;
                    if (LoopCounter >= MAXNOBINSTRINGTOCHECK) {
                        cout << "Maximum loop limit exceeds!" << endl;
                        stopIncrHammingDist = true;
                        break;
                    }

                    NoOfBitStringsProcessed++;
                    if (verbose)
                        cout << "Testing Hamming string: "
                             << string(itrHamming->second.begin(), itrHamming->second.end()) << endl;
                    //Optimzed part Task 4: prune away Hamming String that contain incompatible halfspace pairs
                    bool isValid = true;
                    for (msItr = InValidHammingStr.begin(); msItr != InValidHammingStr.end(); msItr++) {
                        //cout << "Hamming Distance: " << msItr->first << ", String=" << msItr->second << endl;
                        substrPos1 = msItr->second.find("|");
                        HSidx1 = atoi((msItr->second.substr(0, substrPos1)).c_str());
                        substrPos2 = msItr->second.find("|", substrPos1 + 1);
                        HSidx2 = atoi((msItr->second.substr(substrPos1 + 1, substrPos2 - substrPos1 - 1)).c_str());
                        string tmpString = msItr->second.substr(substrPos2 + 1);
                        if (itrHamming->second[HSidx1] == tmpString[0] && itrHamming->second[HSidx2] == tmpString[1]) {
                            //cout << "String " << msItr->second << " matched pattern in " << msItr->second << endl;
                            isValid = false;
                            break;
                        }
                    }
                    if (!isValid) {
                        NoOfPrunedBitStrings++;
                        if (verbose)
                            cout << "Current testing Hamming string is invalid!" << endl;
                        if (verbose)
                            cout << "#pruned bit-strings: " << NoOfPrunedBitStrings << endl;

                        continue;
                    }
                    //end of Optimized Task 4

                    if (verbose)
                        cout << "Mininmal Order=" << minOrder << endl;

                    //prune current Hamming string if there exist two halfspaces that do not compatible
                    if (HammingDistance != itrHamming->first) {
                        if (verbose)
                            cout << "Intersecting halfspaces with hamming distance = " << itrHamming->first << endl;
                        HammingDistance = itrHamming->first;
                    }

                    //generate halfspaces for intersection based on hamming distance
                    char tmpBuf[64];
                    myitoa(itr1->second->NodeID, tmpBuf, 10);
                    strcpy(halfspaceFileName, namePrefix);
                    strcat(halfspaceFileName, tmpBuf);
                    strcat(halfspaceFileName, "_");
                    //strcat(halfspaceFileName,itrHamming->second);
                    //string tmpStr=string(itrHamming->second.begin(),itrHamming->second.end());
                    myitoa(rand(), tmpBuf, 10);
                    strcat(halfspaceFileName, tmpBuf);
                    strcat(halfspaceFileName, nameSuffix);
                    //cout << "Hamming string will be written to " << halfspaceFileName << endl;

                    bool nonZeroExtent = GenHammingHalfSpaces(halfspaceFileName, Dimen, itrHamming->second,
                                                              itr1->second->intersectedHalfspace, itr1->second->MBR);
                    if (nonZeroExtent == false) {
                        //cout << "Discard Hamming binstring " << itrHamming->second << ", for zero-extent! " << endl;
                        NoOfZeroExtentBinStrings++;
                        continue;
                    }

                    //perform intersection for all the halfspaces inside the set 'intersected halfspaces' of a leaf node
                    char sys_string[4096] = "./qhalf Fp < ";
                    strcat(sys_string, halfspaceFileName);
                    strcat(sys_string, " | ./qconvex FA > ");
                    strcat(sys_string, volumeFilename);
                    system(sys_string);
                    if (verbose)
                        cout << "Command: " << sys_string << " performed..." << endl;

                    //open the file to read the volumes of the intersection of the halfspaces
                    ifstream fp_in(volumeFilename, ios::in);
                    string fileLine, volText;
                    size_t substrPos;
                    if (verbose)
                        cout << "Read and output cell volume:" << endl;
                    while (getline(fp_in, fileLine)) {
                        substrPos = fileLine.find("volume:");
                        if (substrPos != std::string::npos)
                            volText = fileLine.substr(substrPos + 7);
                    }
                    fp_in.close();

                    float volume = atof(volText.c_str());
                    if (verbose)
                        cout << "Volume=" << volume << endl;

                    string tmpString = halfspaceFileName;
                    FilesToRemove.push_back(tmpString);  //collect the files for removal later

                    if (volume > ZEROEXTENT)  //ZEROEXTENT: 1e-10
                    {
                        //store the min-cell found
                        if (minOrder >= (HammingDistance + NoOfCoveredHS)) {
                            minOrder = HammingDistance + NoOfCoveredHS;
                            set<long int> tmpSet;
                            std::copy((itr1->second->intersectedHalfspace).begin(),
                                      (itr1->second->intersectedHalfspace).end(),
                                      std::inserter(tmpSet, tmpSet.begin()));
                            minCellHalfSpaces.push_back(tmpSet);
                            vector<char> tmpVec = itrHamming->second;
                            binaryString.push_back(tmpVec);
                            if (verbose)
                                cout << "Node: " << itr1->second->NodeID << ", found a min-cell with binstring:"
                                     << string(itrHamming->second.begin(), itrHamming->second.end()) << endl;
                        } else {
                            NoOfDiscardedCells++;
                            stopIncrHammingDist = true;
                            string tmpString = halfspaceFileName;
                            FilesToRemove.push_back(tmpString);//collect the files for removal later
                            break;
                        }
                    } else {
                        NoOfZeroExtentBinStrings++;
                        if (verbose)
                            cout << "Node: " << itr1->second->NodeID << ", HammingDist: " << HammingDistance
                                 << " is empty!" << endl;
                    }
                    if (verbose)
                        cout << "#min-cells found so far:" << minCellHalfSpaces.size() << endl;
                }
                if (!stopIncrHammingDist) HammingDistance++;
                if (minOrder < INT_MAX && minOrder < (HammingDistance + NoOfCoveredHS))
                    break;
                else {
                    if (verbose)
                        cout << "Next, process bit-string with Hamming dist: " << HammingDistance << endl;
                    //getchar();
                }
            }

            if (FilesToRemove.size() >= 100) {
                vector<string>::iterator vItr = FilesToRemove.begin();
                while (!FilesToRemove.empty()) {
                    remove((*vItr).c_str());  //delete the halfspace data file
                    FilesToRemove.erase(vItr);
                    vItr = FilesToRemove.begin();
                }
            }
        }
        //if (minCellHalfSpaces.size()>0) NotFoundAllMinCells=false;
    }
    vector<string>::iterator vItr = FilesToRemove.begin();
    while (!FilesToRemove.empty()) {
        remove((*vItr).c_str());  //delete the halfspace data file
        FilesToRemove.erase(vItr);
        vItr = FilesToRemove.begin();
    }

    if (verbose) cout << "Number of invalid leaf nodes = " << NoOfInvalidLeaves << endl;
    if (verbose) cout << "#total bit-strings processed: " << NoOfBitStringsProcessed << endl;
    if (verbose) cout << "#pruned bit-strings: " << NoOfPrunedBitStrings << endl;
    if (verbose)
        cout << "#min-Cells found: " << minCellHalfSpaces.size() << ", minOrder=" << minOrder + NoOfCoveredHS << endl;

    totalNoOfBitStringsProcessed = totalNoOfBitStringsProcessed + NoOfBitStringsProcessed;
    totalNoOfPrunedBitStrings = totalNoOfPrunedBitStrings + NoOfPrunedBitStrings;
    totalNoOfZeroExtentBinStrings = totalNoOfZeroExtentBinStrings + NoOfZeroExtentBinStrings;
    totalNoOfDiscardedCells = totalNoOfDiscardedCells + NoOfDiscardedCells;

    if (minOrder == INT_MAX) return NoOfCoveredHS;
    return minOrder + NoOfCoveredHS;
}

bool QuadTree::testHalfspacePair(long int HS1, long int IdxHS1, long int HS2, long int IdxHS2, float subDataSpace[],
                                 multimap<int, string> &InValidHammingStr) //test whether two halfspaces are compatible w.r.t Hamming distance 00,01,10,11
{

    map<long int, long int>::iterator IntInt_mItr;
    typedef map<long int, long int>::value_type IntIntVT;

    typedef multimap<int, string>::value_type msVT;

    char HammingStr[4][3] = {"00", "01", "10", "11"};

    bool interiorPtExists;
    float InteriorPt[MAXDIM];

    //randomization process to generate an interior point
    long int loops = 0;
    long int count = 0;
    int HammingDistance;

    for (int i = 0; i < 4; i++) {
        loops = 0;
        interiorPtExists = true;
        while (true) {
            loops++;
            if (loops >= MAXLOOP + 200) //low probability that an interior point exists for current halfspaces
            {
                //cout << "exceed the maximal loop limits, none interior point exists!" << endl;
                interiorPtExists = false;
                break;
            }
            for (int j = 0; j < Dimen; j++)
                InteriorPt[j] =
                        subDataSpace[j] + (subDataSpace[Dimen + j] - subDataSpace[j]) * (float(rand()) / RAND_MAX);

            count = 0;
            if (strcmp(HammingStr[i], "00") == 0)    //for case: ax_1+bx_2+... > d, lx_1+mx_2+... > t
            {
                HammingDistance = 0;
                float sum = 0;
                for (int j = 0; j < Dimen; j++)
                    sum = sum + HalfSpaces[HS1][j] * InteriorPt[j];
                if (sum > HalfSpaces[HS1][Dimen]) count++;
                sum = 0;
                for (int j = 0; j < Dimen; j++)
                    sum = sum + HalfSpaces[HS2][j] * InteriorPt[j];
                if (sum > HalfSpaces[HS2][Dimen]) count++;
            } else if (strcmp(HammingStr[i], "01") == 0)   //for case: ax_1+bx_2+... > d, lx_1+mx_2+... <= t
            {
                HammingDistance = 1;
                float sum = 0;
                for (int j = 0; j < Dimen; j++)
                    sum = sum + HalfSpaces[HS1][j] * InteriorPt[j];
                if (sum > HalfSpaces[HS1][Dimen]) count++;
                sum = 0;
                for (int j = 0; j < Dimen; j++)
                    sum = sum + HalfSpaces[HS2][j] * InteriorPt[j];
                if (sum <= HalfSpaces[HS2][Dimen]) count++;
            } else if (strcmp(HammingStr[i], "10") == 0)   //for case: ax_1+bx_2+... <= d, lx_1+mx_2+... > t
            {
                HammingDistance = 1;
                float sum = 0;
                for (int j = 0; j < Dimen; j++)
                    sum = sum + HalfSpaces[HS1][j] * InteriorPt[j];
                if (sum <= HalfSpaces[HS1][Dimen]) count++;
                sum = 0;
                for (int j = 0; j < Dimen; j++)
                    sum = sum + HalfSpaces[HS2][j] * InteriorPt[j];
                if (sum > HalfSpaces[HS2][Dimen]) count++;
            } else if (strcmp(HammingStr[i], "11") == 0)   //for case: ax_1+bx_2+... <= d, lx_1+mx_2+... <= t
            {
                HammingDistance = 2;
                float sum = 0;
                for (int j = 0; j < Dimen; j++)
                    sum = sum + HalfSpaces[HS1][j] * InteriorPt[j];
                if (sum <= HalfSpaces[HS1][Dimen]) count++;
                sum = 0;
                for (int j = 0; j < Dimen; j++)
                    sum = sum + HalfSpaces[HS2][j] * InteriorPt[j];
                if (sum <= HalfSpaces[HS2][Dimen]) count++;
            }
            if (count == 2) {
                //cout << "We have found an interior point!" << endl;
                break;
            }
        }
        if (!interiorPtExists) {
            //long int Val1=HS1,Val2=HS2;        //form incompatible pairs by IDs of the halfspaces
            long int Val1 = IdxHS1, Val2 = IdxHS2;    //form incompatible pairs by indices of the halfspaces in set 'intersectedHalfSpace'

            if (Val1 > Val2) {
                long int tmpInt;
                tmpInt = Val1;
                Val1 = Val2;
                Val2 = tmpInt;
            }
            char m_Buf[1024], m_Buf1[256];
            myitoa(Val1, m_Buf1, 10);
            strcpy(m_Buf, m_Buf1);
            strcat(m_Buf, "|");
            myitoa(Val2, m_Buf1, 10);
            strcat(m_Buf, m_Buf1);
            strcat(m_Buf, "|");
            strcat(m_Buf, HammingStr[i]);
            string tmpString = m_Buf;

            InValidHammingStr.insert(msVT(HammingDistance, tmpString));

            //cout << "Hamming Distance: " << HammingDistance << ", String=" << tmpString << endl;

        }
    }
    return true;
}

bool QuadTree::readCombinations() {

    FILE *fp;
    char *token;
    char m_separator[] = " \n\t";
    char buf[512];
    string FileName[] = {"Comb2D.txt", "Comb3D.txt", "Comb4D.txt", "Comb5D.txt", "Comb6D.txt", "Comb7D.txt",
                         "Comb8D.txt", "Comb9D.txt"};
    string strComb;

    fp = fopen(FileName[Dimen - 2].c_str(), "r");
    if (fp == NULL) {
        cout << "error in fileopen!" << endl;
        exit(0);
    }

    fgets(buf, 512, fp);
    if (atoi(buf) != Dimen) {
        cout << "Error! Dimensions are not equal!" << endl;
        exit(0);
    }
    while (fgets(buf, 512, fp) != NULL) {
        token = strtok(buf, m_separator);
        //while (token != NULL){
        //	token = strtok(NULL,m_separator);
        //}
        string strComb = token;
        Comb.push_back(strComb);
    }
    fclose(fp);

    if (verbose)
        cout << "QuadTree Dimen=" << Dimen << ", reading combination finished!" << endl;

    return true;
}
