#include "stdio.h"
#include "S3.h"
#include "rtree.h"
#include "rentry.h"
#include "rnode.h"
#include "filemem.h"
#include "mainmem.h"
#include "mem.h"
#include "tgs.h"
#include "hypercube.h"
#include "psdraw.h"
#include "virtualRnode.h"

#include <iostream>
#include <map>
#include <string.h>
#include <climits>
#include <algorithm>

//headers from partial hull
#include <float.h>
#include <math.h>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <locale.h>
#include <string.h>
#include <ctype.h>
#include <set>
//#include "PartialHull.h"
//end of headers from partial hull

//headers from QuadTree
#include "QuadTree.h"
//end of QuadTree headers

#define MAXPAGEID 4999999
//#ifndef _Globals
#define DMAX 10    //max dimensionality
//#define SIDELEN 2
//#define _Globals
//#endif

#define OPTDIM 3

extern int CEILING;
extern int MaxQuadTreeLevels;
extern int QuadNodeCapacity;
extern int maximumLevel;
extern long numOfTreeNodes;
extern long numOfLeafNodes;
extern long numOfNodeSplits;
extern bool optWithinNodeIntersection;
extern float Queryspace[2 * MAXDIM + 1];
extern vector<vector<float> > HalfSpaces;  //format: (coeff_1, coeff_2, ..., coeff_d, offset)
extern long int NoOfNewlyAddedHalfSpaces;
extern long int NoOfTotalHalfspacesAdded;
extern map<long int, long int> entriesInLists;
extern long int NoOfEntriesInLists;
extern long int numOfLeavesToIntersect;
extern double timeBuildQuadTree;
extern double timeNodeIntersection;
extern bool verbose;

vector<std::pair<long, QuadNode *> > Leaves;
map<long int, long int> RdIDtoHalfplaneID;

using namespace std;


S3::S3() {
}

S3::~S3() {}

/*
int S3::GetTopK(
        const int dimen, //dimensionality
        Rtree& a_rtree,             // returns a result and a result size
        Query& a_query,
        std::multimap<long int, VirtualRNode*>& NonResultEntry,  //non-result set, storing rtree node which is simply a data entry
        std::vector<long int>& PrunedNodes, //the heap for storing non-result r-tree nodes and remaining entries in H
        const int K,  //the k for top-k query
        int& maxStackSize, Array& a_pageaccessed, long *nonSky10)
{

   multimap<long int, VirtualRNode*> DataEntryQ;  //queue storing rtree node which is simply a data entry
   multimap<long int, VirtualRNode*>::iterator IntVRN_Iter;
   typedef multimap<long int, VirtualRNode*>::value_type IntVRN_Pair;

   multimap<float, long int, std::greater<float> > H;  //the max-heap H
   multimap<float, long int>::iterator DblInt_Iter;
   typedef multimap<float, long int>::value_type DblInt_Pair;
   typedef multimap<float, float>::value_type DblDbl_Pair;


   float topscore=-1e100;

   //these variables are for Step 1 optimization for solution FP
   set<float> Nonscores;
   set<float>::iterator sDbl_Iter;
   float bottomNonscore=-1e100;
   int NonSkyPtr=0;

   float m_stack=-1;
   float m_SizeOfDataEntQ=-1;
   float score=0;
   float pt[DMAX];
   bool isAnObject;

   a_query.score.clear();
   //NonResultSet.clear();
   NonResultEntry.clear();
   PrunedNodes.clear();

   H.insert(DblInt_Pair(0,a_rtree.m_memory.m_rootPageID));
   while (H.size()!=0)
   {
        isAnObject=false;

        if (m_stack<=H.size())
        {
            m_stack=H.size();
            m_SizeOfDataEntQ=DataEntryQ.size();
        }

        //maxStackSize = s.size() > maxStackSize ? s.size() : maxStackSize;
    	DblInt_Iter=H.begin();
    	long int pageid = DblInt_Iter->second;
        float ScoreTmp= DblInt_Iter->first;

    	H.erase(DblInt_Iter);

    	VirtualRNode* VirNode=new VirtualRNode;
        RtreeNodeEntry* e0;            // create node entry e0 for node n, so that its MBR can be obtained
    	if (pageid>=MAXPAGEID)     //current element in H is a data entry node (for distinction,its pageid is equal to m_id+MAXPAGED)
    	{
            isAnObject=true;

            IntVRN_Iter=DataEntryQ.find(pageid);
            if (IntVRN_Iter==DataEntryQ.end())
            {
            	cout << "Error! there is no node " << pageid << " in DataEntryQ!" << endl;
            }
            else
            {
            	VirNode->copyData(IntVRN_Iter->second);

            	if (VirNode->m_usedSpace>1)
                {
            	   cout << "There is a fatal error in retrieving data-entry node from queue!" << endl;
                   return -1;
                }
            	else
                   e0 = VirNode->m_entry[0]->clone();  //retrieve the sole entry from the data-entry node

                delete IntVRN_Iter->second;         //free the virtualRnode; Note, this operation is VERY important
                DataEntryQ.erase(IntVRN_Iter);
            }
            pageid=pageid-MAXPAGEID;
    	}
    	else
    	{
            RtreeNode* node=a_rtree.m_memory.loadPage(pageid);

            //my code
            VirNode->copyData(*node);
            VirNode->copyEntries(*node,node->m_usedSpace);
            e0 = node->genNodeEntry();     //compute the enclosing MBR for current index node

            delete node;
        }
        //cout << "examining page: " << pageid << " is an object "<< isAnObject << ", Heap size: " << H.size() << " stack " << m_stack << endl;
        //cout << "topscore= " << topscore << endl;

        if (isAnObject)  //if current entry is a data object, then compute its score by using its values
        {
            for (int i=0;i<dimen;i++)
                 pt[i]=VirNode->m_entry[0]->m_hc.getUpper()[i];
            score=GetScore(pt,a_query.function,dimen);
        }
        else   //if current entry is NOT a data object, then compute score by using the top-right corner of its MBR
        {
            for (int i=0; i<dimen; i++)
                 pt[i]=e0->m_hc.getUpper()[i];
            score=GetScore(pt,a_query.function,dimen);
        }

        //current node is below topscore, so we explore the node and put all its data points to non-result set.
        if (score<=topscore)
        {
            if (VirNode->isLeaf())
            {
                if (VirNode->m_usedSpace>1)
                {
                    //a_pageaccessed.append((void*)pageid);   //*******************
                    PrunedNodes.push_back(pageid);  //keep the non-result rtree node for computation later
                }
                else
                {
                         //long int pid=VirNode->m_entry[0]->m_id+MAXPAGEID;
                         long int pid=pageid+MAXPAGEID;

                         VirtualRNode* vNode=new VirtualRNode;
                         vNode->copyData(VirNode);
                         NonResultEntry.insert(IntVRN_Pair(pid,vNode));

                         //code below for Step 1 optimization for solution FP
                         if (NonSkyPtr<NUMNONSKY)
                         {
                             nonSky10[NonSkyPtr++]=pageid;
                             Nonscores.insert(score);
                             sDbl_Iter=Nonscores.begin();
                             bottomNonscore=(*sDbl_Iter);
                         }
                         else
                         {
                            if (score>=bottomNonscore)
                            {
                                 nonSky10[NonSkyPtr++]=pageid;
                                 Nonscores.insert(score);

                                 sDbl_Iter=Nonscores.begin();
                                 bottomNonscore=(*sDbl_Iter);
                             }
                         }
                }
            }
            else //current entry is an index node, we kep it for later use
            {
                //a_pageaccessed.append((void*)pageid);  //**********************
                PrunedNodes.push_back(pageid);  //keep the non-result rtree node for skyline computation later
            }
            delete VirNode;
            delete e0;

            continue;
        }//end compute score

        //check whether current entry is leaf node or index node
        if (VirNode->isLeaf())
        {
        	    if (VirNode->m_usedSpace>1)   //current node is a leaf node, so all its data entries must be put into priority queue H
        	    {
        	    	a_pageaccessed.append((void*)pageid);
                        for (int i=0; i<VirNode->m_usedSpace; i++)
                        {
                             long int NegPageid=VirNode->m_entry[i]->m_id+MAXPAGEID;
                	     VirtualRNode* node=new VirtualRNode; //for each data entries of n, insert it into data-entry queue

                	     RtreeNodeEntry* Nentry=VirNode->m_entry[i]->clone();
                	     node->insertEntry(Nentry);
                	     DataEntryQ.insert(IntVRN_Pair(NegPageid,node));
                	     delete Nentry;

                             //compute the score of current entry w.r.t query q
                             for (int j=0;j<dimen;j++)
                                  pt[j]=VirNode->m_entry[i]->m_hc.getUpper()[j];
                             score=GetScore(pt,a_query.function,dimen);

                             H.insert(DblInt_Pair(score,NegPageid));
                         }
        	    }
        	    else    //current node is a data entry node
        	    {
                        multimap<float,long int, std::greater<float> > entryHeap;
                        multimap<float,long int>::iterator entryIter;
                        entryHeap.clear();

        		for (int i=0; i<VirNode->m_usedSpace; i++)
                        {
                             //compute the score of current entry w.r.t query q
                             for (int j=0;j<dimen;j++)
                                  pt[j]=VirNode->m_entry[i]->m_hc.getUpper()[j]-SIDELEN;
                             score=GetScore(pt,a_query.function,dimen);

                             entryHeap.insert(DblInt_Pair(score,i));
                        }
                        int processedEntries=0;
                        while (processedEntries<VirNode->m_usedSpace)
        		{
                             entryIter=entryHeap.begin();
                             int Pos=entryIter->second;
                             score=entryIter->first; //point (x,y) is expanded to polygon (x-SIDELEN,y-SIDELEN)(x+SIDELEN,y+SIDELEN)

                             //cout << "found a top-k result " << VirNode->m_entry[Pos]->m_id << " with score "<< score << endl;

                             //insert current point into the top-k list of query q
                             a_query.score.insert(DblInt_Pair(score,VirNode->m_entry[Pos]->m_id));
                             if (a_query.score.size()>K)
                             {
                                 int pid=(a_query.score.begin())->second;

 			         a_query.score.erase(a_query.score.begin());
				 topscore = (a_query.score.begin())->first; //update topscore
                                 //NonResultSet.push_back(pid);        //for storing id of pruned points only

                                 VirtualRNode* vNode=new VirtualRNode;
                                 vNode->copyData(VirNode);
                                 NonResultEntry.insert(IntVRN_Pair(pid+MAXPAGEID,vNode));

                                 //code below for Step 1 optimization for solution FP
                                 if (NonSkyPtr<NUMNONSKY)
                                 {
                                     nonSky10[NonSkyPtr++]=pid;
                                     Nonscores.insert(score);
                                     sDbl_Iter=Nonscores.begin();
                                     bottomNonscore=(*sDbl_Iter);
                                 }
                                 else
                                 {
                                     if (score>=bottomNonscore)
                                     {
                                         nonSky10[NonSkyPtr++]=pid;
                                         Nonscores.insert(score);

                                         sDbl_Iter=Nonscores.begin();
                                         bottomNonscore=(*sDbl_Iter);
                                     }
                                 }
                             }
                             else if (a_query.score.size() == K)
				      topscore = (a_query.score.begin())->first; //set topscore

                             entryHeap.erase(entryIter);
                             processedEntries++;
                        }
                    }
        }
        else      //current entry is an index, so all its entries are inserted into the heap
        {
                  //cout << "expand index node "<< pageid << endl;
             	  a_pageaccessed.append((void*)pageid);
                  for (int i=0; i<VirNode->m_usedSpace; i++)
                  {
                       //compute the score of current entry w.r.t query q
                       for (int j=0;j<dimen;j++)
                            pt[j]=VirNode->m_entry[i]->m_hc.getUpper()[j];
                       score=GetScore(pt,a_query.function,dimen);
                       H.insert(DblInt_Pair(score,VirNode->m_entry[i]->m_id));
                  }
        }//end check whether it's leaf or index node

        delete VirNode;
        delete e0;
   }
   //cout << "#points remained in working heap: " << DataEntryQ.size() << endl;
   //cout << "size of heap H: " << H.size() << endl;
   //cout << "#pruned index nodes: " << PrunedNodes.size() << endl;
   //cout << "Top-k computation finished!"<< endl;

   return 1;
}
//*/

/*
int S3::Step1(
        const int dimen, //dimensionality
        Query& a_query,
        //std::multimap<long int, VirtualRNode*>& NonResultEntry,  //non-result set, storing rtree node which is simply a data entry
        std::multimap<long,long>& NonResultPtID,
        char* PGstr[],
        float* PG[],
        long *nonSky10,
        long* HullIDtoPGID[],
        long& IDMapPtr,
        PartialHull &myHull)
{
   long int i, j, k, m, n, pos;
   long int sizeOfdominated=0;
   long int numOfPtsInsideHull=0;

   char buf[512];   //the buffer to store the data point for partial hull
   char *pSite;
   set<long> FirstDpoints;
   set<long>::iterator Int_sIter;

   multimap<long,long>::iterator IntInt_Iter;
   multimap<long,long>::reverse_iterator IntInt_rIter;
   typedef multimap<long, long>::value_type IntInt_Pair;

   multimap<double,long> RankNonResult[DMAX];
   multimap<double,long>::reverse_iterator DblInt_rIter;
   multimap<double,long>::iterator DblInt_Iter;
   typedef multimap<double, long>::value_type DblInt_Pair;

   std::multimap<long int, VirtualRNode*> TmpVec;  //temp vector
   std::multimap<long int, VirtualRNode*>::iterator IntVRN_Iter;
   typedef multimap<long int, VirtualRNode*>::value_type IntVRN_Pair;

   float pt[DMAX];
   float kthPt[DMAX];
   int kthID=(a_query.score.begin())->second;    //id of the k-th data point

   //check the status of myHull
   if (myHull.isInitialized==false) {cout << "Please initialize the partial hull with myHull.isInitialized(dim)!"<< endl; return -1;}

   //
   //map id of the k-th result to the id of the first partial hull vertex (which is 0)
   //HullIDtoPGID[IDMapPtr]=new long [2];
   //HullIDtoPGID[IDMapPtr][0]=0;   //the hull id of the k-th result
   //HullIDtoPGID[IDMapPtr][1]=kthID;   //the original id of the k-th result
   //IDMapPtr++;
   ///

   //retrieve the k-th result
   //strcpy(buf,"");
   //for (i=0;i<dimen;i++)
   //{
        //kthPt[i]=PG[kthID][i]+SIDELEN;
        //sprintf(buf+strlen(buf),"%lf ",kthPt[i]);
        //strcat(buf," ");
   //}
   ///
   //sprintf(buf, "\n");

   strcpy(buf,PGstr[kthID]);

   //use the k-th result as the first point for building the initial convex hull_infinity
   myHull.getInitSimplex(buf);

   //
   //prune away the non-result points that are dominated by the k-th result
   //while (NonResultEntry.size()!=0)
   //{
     //   IntVRN_Iter=NonResultEntry.begin();

     //   long Pos=IntVRN_Iter->first-MAXPAGEID;
     //   for (i=0;i<dimen;i++) pt[i]=PG[Pos][i]+SIDELEN;

     //   bool dominated=IsDominatedBy(dimen, pt, kthPt);
     //   if (dominated)
     //   {
     //       sizeOfdominated++;
     //       delete IntVRN_Iter->second;
     //   }
     //   else
     //       TmpVec.insert(IntVRN_Pair(IntVRN_Iter->first, IntVRN_Iter->second));

     //   NonResultEntry.erase(IntVRN_Iter);
   //}
   //NonResultEntry=TmpVec;

   //while (TmpVec.size()!=0)
   //{
    //    IntVRN_Iter=TmpVec.begin();
    //    NonResultEntry.insert(IntVRN_Pair(IntVRN_Iter->first, IntVRN_Iter->second));
    //    TmpVec.erase(IntVRN_Iter);
   //}
   //cout << "number of non-result points dominated by the k-th result: " << sizeOfdominated <<endl;
   //cout << "after pruning, the size of Nonresult: " << NonResultEntry.size() << endl;
   ///

 if (Opt4InitSimp==1)
 {
   //find the points with the maximum value along each dimension
   for (IntInt_Iter=NonResultPtID.begin();IntInt_Iter!=NonResultPtID.end();IntInt_Iter++)
   {
        long Pos=IntInt_Iter->first;
        for (i=0;i<dimen;i++)
        {
             pt[i]=PG[Pos][i]+SIDELEN;
             RankNonResult[i].insert(DblInt_Pair(pt[i],Pos));
        }
   }
   for (i=0;i<dimen;i++)
   {
        //printf("size of the %d-th map: %d\n",i,RankSkyline[i].size());
        DblInt_rIter=RankNonResult[i].rbegin();
        int szSet=FirstDpoints.size();
        while (DblInt_rIter!=RankNonResult[i].rend())
        {
            FirstDpoints.insert(DblInt_rIter->second);
            if (szSet==FirstDpoints.size())
                DblInt_rIter++;
            else
                break;
        }
   }

   //use the first d data points to compute the initial convex hull
   for (Int_sIter=FirstDpoints.begin();Int_sIter!=FirstDpoints.end();Int_sIter++)
   {
        //
        //printf("FirstDpoints: %ld\n",(*Int_sIter));
        //strcpy(buf,"");
        //for (i=0;i<dimen;i++)
        //{
        //     sprintf(buf+strlen(buf),"%lf ",PG[(*Int_sIter)][i]+SIDELEN);
        //     strcat(buf," ");
        //}
        ///
        //sprintf(buf+strlen(buf),"\n");

        strcpy(buf,PGstr[(*Int_sIter)]);

        //use the point for building initial simplex
        myHull.getInitSimplex(buf);

        //
        //build the mapping from id of partial hull vertex to original id of the points
        //HullIDtoPGID[IDMapPtr]=new long [2];
        //HullIDtoPGID[IDMapPtr][0]=myHull.numInitSitesRead-1;   //the id of partial hull vertex
        //HullIDtoPGID[IDMapPtr][1]=(*Int_sIter);   //the original id of the point
        //IDMapPtr++;
        ///
   }
   if (myHull.numInitSitesRead<dimen+1) cout << "Caution! there are less than d+1 points for building the initial hull!" << endl;

   //cout << "initial simplex got, beging building pHull..." << endl;

   //the initial simplex has been built; preceed to compute the facets
   //adjacent to the k-th result (i.e., maintain the partial convex hull)
   for (IntInt_Iter=NonResultPtID.begin();IntInt_Iter!=NonResultPtID.end();IntInt_Iter++) //this version is the original access order of points
   {
          long Pos=IntInt_Iter->first;        //this version is the original access order of points
          Int_sIter=FirstDpoints.find(Pos);
          if (Int_sIter!=FirstDpoints.end()) continue;
          //
          //strcpy(buf,"");
          //for (i=0;i<dimen;i++)
          //{
            //   sprintf(buf+strlen(buf),"%lf ",PG[Pos][i]+SIDELEN);
            //   strcat(buf," ");
          //}
          ///
          //sprintf(buf+strlen(buf),"\n");
          //printf("%s\n",buf);

          strcpy(buf,PGstr[Pos]);

          //build the mapping from id of partial hull vertex to original id of the points
          //HullIDtoPGID[IDMapPtr]=new long [2];

          if (myHull.cdim<myHull.rdim)
          {
              //cout << "feeding " << buf << " to getInitSimplex..." << endl;
              myHull.getInitSimplex(buf);
              //HullIDtoPGID[IDMapPtr][0]=myHull.numInitSitesRead-1;   //the id of partial hull vertex
          }
          else
          {
              myHull.buildPartialHull(buf);
              //HullIDtoPGID[IDMapPtr][0]=myHull.pnum-2;   //the id of partial hull vertex
          }
          if (myHull.isProned==true) numOfPtsInsideHull++;

          //
          //build the mapping from id of partial hull vertex to original id of the points
          //HullIDtoPGID[IDMapPtr][1]=Pos;   //the original id of the point
          //IDMapPtr++;
          ///
   }
   ///
 }
 else
 {
   for (IntInt_rIter=NonResultPtID.rbegin();IntInt_rIter!=NonResultPtID.rend();IntInt_rIter++) //the original access order of points
   {
          long Pos=IntInt_rIter->first;        //this version is the original access order of points
          /
          //strcpy(buf,"");
          //for (i=0;i<dimen;i++)
          //{
          //     sprintf(buf+strlen(buf),"%lf ",PG[Pos][i]+SIDELEN);
          //     strcat(buf," ");
         // }
          ///
          //sprintf(buf+strlen(buf),"\n");
          //printf("%s\n",buf);

          strcpy(buf,PGstr[Pos]);

          //build the mapping from id of partial hull vertex to original id of the points
          //HullIDtoPGID[IDMapPtr]=new long [2];

          if (myHull.cdim<myHull.rdim)
          {
              myHull.getInitSimplex(buf);
              //HullIDtoPGID[IDMapPtr][0]=myHull.numInitSitesRead-1;   //the id of partial hull vertex
          }
          else
          {
              myHull.buildPartialHull(buf);
              //HullIDtoPGID[IDMapPtr][0]=myHull.pnum-2;   //the id of partial hull vertex
          }
          if (myHull.isProned==true) numOfPtsInsideHull++;

          //
          //build the mapping from id of partial hull vertex to original id of the points
          //HullIDtoPGID[IDMapPtr][1]=Pos;   //the original id of the point
          //IDMapPtr++;
          ///
   }
 }

   //cout << "number of non-results: " << NonResultPtID.size() << endl;
   cout << "number of non-results proned: " << numOfPtsInsideHull << endl;
   //cout << "The initial convex hull contains " << myHull.numInitSitesRead << " points." << endl;
   //cout << "Step 1 finished!" << endl;

    return 1;
}
//*/

/*
int S3::Step2(
        const int dimen, //dimensionality
        Rtree& a_rtree,             // returns a result and a result size
        Query& a_query,
        std::vector<long int>& PrunedNodes, //the heap for storing non-result r-tree nodes and remaining entries in H
        char* PGstr[],
        float* PG[],
        long *nonSky10,
        long* HullIDtoPGID[],
        long& IDMapPtr,
        long& NumOfPtsInHull,
        long& NumOfNodesInHull,
        int& maxStackSize, Array& a_pageaccessed,
        PartialHull &myHull)
{
   multimap<float, long int> H0;  //the min-heap H0
   multimap<float, long int>::iterator DblInt_Iter;
   typedef multimap<float, long int>::value_type DblInt_Pair;
   typedef multimap<float, float>::value_type DblDbl_Pair;

   multimap<long int, VirtualRNode*> DataEntryQ;  //queue storing rtree node which is simply a data entry
   multimap<long int, VirtualRNode*>::iterator IntVRN_Iter;
   typedef multimap<long int, VirtualRNode*>::value_type IntVRN_Pair;

   vector<long int>::iterator Int_vIter;

   float mindist;
   float m_stack=-1;
   float pt[DMAX];
   float nonPt[DMAX];
   float ORIGIN[DMAX];
   float CornerPoint[DMAX];
   float kthPt[DMAX];
   int kthID=(a_query.score.begin())->second;    //id of the k-th data point
   bool isAnObject;

   int MBRisInsideHull;
   char buf[512];

   for (int i=0;i<dimen;i++)
   {
        kthPt[i]=PG[kthID][i]+SIDELEN;
        ORIGIN[i]=CEILING+20;
   };

   //put the non-result points and pruned index entries to a min-heap sorted on mindist
   RtreeNodeEntry* e0;
   for (Int_vIter=PrunedNodes.begin();Int_vIter!=PrunedNodes.end();Int_vIter++)
   {
        RtreeNode* Rnode=a_rtree.m_memory.loadPage(*Int_vIter);
        e0 = Rnode->genNodeEntry();     //compute the enclosing MBR for current index node
        for (int i=0;i<dimen;i++) pt[i]=e0->m_hc.getUpper()[i];
        mindist=minDist(pt,ORIGIN,dimen);
        H0.insert(DblInt_Pair(mindist,*Int_vIter));
        delete Rnode;
        delete e0;
   }
   //cout << "size of H0: " << H0.size() << endl;

   //begin explore the min-heap to refine the initial facets adjacent to the k-th result
   while (H0.size()!=0)
   {
        isAnObject=false;

        if (m_stack<=H0.size())
        {
            m_stack=H0.size();
        }

        //maxStackSize = s.size() > maxStackSize ? s.size() : maxStackSize;
    	DblInt_Iter=H0.begin();
    	long int pageid = DblInt_Iter->second;
        float ScoreTmp= DblInt_Iter->first;

    	H0.erase(DblInt_Iter);

    	VirtualRNode* VirNode=new VirtualRNode;
        RtreeNodeEntry* e0;            // create node entry e0 for node n, so that its MBR can be obtained
    	if (pageid>=MAXPAGEID)     //current element in H is a data entry node (for distinction,its pageid is equal to m_id+MAXPAGED)
    	{
            isAnObject=true;

            IntVRN_Iter=DataEntryQ.find(pageid);
            if (IntVRN_Iter==DataEntryQ.end())
            {
            	cout << "Error! there is no node " << pageid << " in DataEntryQ!" << endl;
            }
            else
            {
            	VirNode->copyData(IntVRN_Iter->second);

            	if (VirNode->m_usedSpace>1)
                {
            	   cout << "There is a fatal error in retrieving data-entry node from queue!" << endl;
                   return -1;
                }
            	else
                   e0 = VirNode->m_entry[0]->clone();  //retrieve the sole entry from the data-entry node

                delete IntVRN_Iter->second;         //free the virtualRnode; Note, this operation is VERY important
                DataEntryQ.erase(IntVRN_Iter);
            }
            pageid=pageid-MAXPAGEID;
    	}
    	else
    	{
            RtreeNode* node=a_rtree.m_memory.loadPage(pageid);

            //my code
            VirNode->copyData(*node);
            VirNode->copyEntries(*node,node->m_usedSpace);
            e0 = node->genNodeEntry();     //compute the enclosing MBR for current index node

            delete node;
        }
        //cout << "exam page " << pageid << s", Heap size: " << H0.size() << " mindist:" << DistTmp << endl;

        for (int i=0;i<dimen;i++)
             pt[i]=VirNode->m_entry[0]->m_hc.getUpper()[i];

        //throw away the points/MBRs inside the hull; if the dimensionality is small, then don't apply this optimization
        if (dimen>OPTDIM)
        {
            bool dominated=false;
            for (int i=0;i<NUMNONSKY;i++)
            {
                 for (int j=0;j<dimen;j++) nonPt[j]=PG[nonSky10[i]][j]+SIDELEN;
                 dominated=IsDominatedBy(dimen, pt, nonPt);
                 if (dominated) break;
            }
            if (dominated)
            {
                delete VirNode;
                delete e0;
                continue;
            }
        }
        if (isAnObject==true)
        {
            strcpy(buf,PGstr[(VirNode->m_entry[0]->m_id)]);
        }
        else
        {
            //begin updating the partial hull by examining MBR and data points
            strcpy(buf,"");
            for (int i=0;i<dimen;i++)
            {
                 //pt[i]=VirNode->m_entry[0]->m_hc.getUpper()[i];

                 sprintf(buf+strlen(buf),"%lf ",pt[i]);
                 strcat(buf," ");
            }
        }
        MBRisInsideHull=myHull.isMBRVisible(buf);
        //
        //cout << "MBR "<< VirNode->m_entry[0]->m_id << ": "<< buf;
        //if (MBRisInsideHull==IS_VISIBLE)
            //cout << " " << "visible" << endl;
        //else
            //cout << " " << " NOT visible" << endl;
        ///

        if (MBRisInsideHull==IS_VISIBLE)
        {
           //check whether current entry is leaf node or index node
            if (VirNode->isLeaf())
            {
        	    if (VirNode->m_usedSpace>1)   //current node is a leaf node, so all its data entries must be examined
        	    {
        	    	a_pageaccessed.append((void*)pageid);
                        for (int i=0; i<VirNode->m_usedSpace; i++)
                        {
                             //
                             //strcpy(buf,"");
                             //for (int j=0;j<dimen;j++)
                             //{
                               //   pt[j]=(VirNode->m_entry[i]->m_hc.getUpper()[j])-SIDELEN;

                                 // sprintf(buf+strlen(buf),"%lf ",pt[j]);
                                 // strcat(buf," ");
                             //}
                             ///

                             //throw away the points/MBRs inside the hull
                             if (dimen>OPTDIM)
                             {
                                 bool dominated=false;
                                 for (int j=0;j<dimen;j++) pt[j]=(VirNode->m_entry[i]->m_hc.getUpper()[j])-SIDELEN;
                                 for (int j=0;j<NUMNONSKY;j++)
                                 {
                                      for (int k=0;k<dimen;k++) nonPt[k]=PG[nonSky10[j]][k]+SIDELEN;
                                      dominated=IsDominatedBy(dimen, pt, nonPt);
                                      if (dominated) break;
                                 }
                                 if (dominated) continue;
                             }

                             strcpy(buf,PGstr[(VirNode->m_entry[i]->m_id)]);

                             //feed the point to partial hull
                             myHull.buildPartialHull(buf);
                             if (myHull.isProned==true) NumOfPtsInHull++;

                             //
                             //build the mapping from id of partial hull vertex to original id of the points
                             //HullIDtoPGID[IDMapPtr]=new long [2];
                             //HullIDtoPGID[IDMapPtr][0]=myHull.pnum-2;   //the id of partial hull vertex
                             //HullIDtoPGID[IDMapPtr][1]=VirNode->m_entry[i]->m_id;   //the original id of the point
                             //IDMapPtr++;
                             ///
                        }
        	    }
        	    else    //if current entry is a data point, then update the initial facets accordingly
                    {
                             //
                             //strcpy(buf,"");
                             //for (int j=0;j<dimen;j++)
                             //{
                                  ////pt[j]=pt[j]-SIDELEN;   //obtain the real values of a data entry
                               //   pt[j]=VirNode->m_entry[0]->m_hc.getUpper()[j]-SIDELEN;

                                 // sprintf(buf+strlen(buf),"%lf ",pt[j]);
                                  //strcat(buf," ");
                             //}
                             ///

                             strcpy(buf,PGstr[(VirNode->m_entry[0]->m_id)]);

                             //feed the point to partial hull
                             myHull.buildPartialHull(buf);
                             if (myHull.isProned==true) NumOfPtsInHull++;

                             //
                             //build the mapping from id of partial hull vertex to original id of the points
                             //HullIDtoPGID[IDMapPtr]=new long [2];
                             //HullIDtoPGID[IDMapPtr][0]=myHull.pnum-2;   //the id of partial hull vertex
                             //HullIDtoPGID[IDMapPtr][1]=VirNode->m_entry[0]->m_id;   //the original id of the point
                             //IDMapPtr++;
                             ///
                    }
            }
            else   //current entry is an index, so all its entries are inserted into the min-heap
            {
             	  a_pageaccessed.append((void*)pageid);
                  for (int i=0; i<VirNode->m_usedSpace; i++)
                  {
                       //check whether current index node is inside the partial hull or not
                       strcpy(buf,"");
                       for (int j=0;j<dimen;j++)
                       {
                            pt[j]=VirNode->m_entry[i]->m_hc.getUpper()[j];   //obtain the top-right corner point of current MBR

                            //sprintf(buf+strlen(buf),"%lf ",pt[j]);
                            //strcat(buf," ");
                       }
                       //
                       //int isInsideHull=myHull.isMBRVisible(buf);

                       //if (isInsideHull!=IS_VISIBLE)
                       //{
                         //  NumOfNodesInHull++;
                           //continue;
                       //}
                       ///

                       //
                       //throw away the points/MBRs inside the hull
                       //if (dimen>OPTDIM)
                       //{
                         //  bool dominated=false;
                           //for (int j=0;j<NUMNONSKY;j++)
                           //{
                                //for (int k=0;k<dimen;k++) nonPt[k]=PG[nonSky10[j]][k]+SIDELEN;
                                //dominated=IsDominatedBy(dimen, pt, nonPt);
                                //if (dominated) break;
                           //}
                           //if (dominated) continue;
                       //}
                       ///

                       mindist=minDist(pt,ORIGIN,dimen);
                       H0.insert(DblInt_Pair(mindist,VirNode->m_entry[i]->m_id));
                  }
            }//end check whether it's leaf or index node
        }
        else
            NumOfNodesInHull++;

        delete VirNode;
        delete e0;
   }

   return 1;
}
//*/

/*
int S3::OptimizedStep1(
        const int dimen, //dimensionality
        Query& a_query,
        //std::multimap<long int, VirtualRNode*>& NonResultEntry,  //non-result set, storing rtree node which is simply a data entry
        std::multimap<long,long>& NonResultPtID,
        char* PGstr[],
        float* PG[],
        long *nonSky10,
        long* HullIDtoPGID[],
        long& IDMapPtr,
        PartialHull &myHull)
{
   long int i, j, k, m, n, pos;
   long int sizeOfdominated=0;
   long int numOfPtsInsideHull=0;

   //data structure for R-tree
   multimap<double, long int, std::greater<float> > H0;  //the max-heap H0
   //


   char buf[512];   //the buffer to store the data point for partial hull
   char *pSite;
   set<long> FirstDpoints;
   set<long>::iterator Int_sIter;

   multimap<long,long>::iterator IntInt_Iter;
   multimap<long,long>::reverse_iterator IntInt_rIter;
   typedef multimap<long, long>::value_type IntInt_Pair;

   multimap<double,long> RankNonResult[DMAX];
   multimap<double,long>::reverse_iterator DblInt_rIter;
   multimap<double,long>::iterator DblInt_Iter;
   typedef multimap<double, long>::value_type DblInt_Pair;

   std::multimap<long int, VirtualRNode*> TmpVec;  //temp vector
   std::multimap<long int, VirtualRNode*>::iterator IntVRN_Iter;
   typedef multimap<long int, VirtualRNode*>::value_type IntVRN_Pair;

   float pt[DMAX];
   float kthPt[DMAX];
   int kthID=(a_query.score.begin())->second;    //id of the k-th data point

   cout << "performing Step 1 pruning..." << endl;

   //check the status of myHull
   if (myHull.isInitialized==false) {cout << "Please initialize the partial hull with myHull.isInitialized(dim)!"<< endl; return -1;}

   //
   //map id of the k-th result to the id of the first partial hull vertex (which is 0)
   //HullIDtoPGID[IDMapPtr]=new long [2];
   //HullIDtoPGID[IDMapPtr][0]=0;   //the hull id of the k-th result
   //HullIDtoPGID[IDMapPtr][1]=kthID;   //the original id of the k-th result
   //IDMapPtr++;
   ///

   strcpy(buf,PGstr[kthID]);

   //use the k-th result as the first point for building the initial convex hull_infinity
   myHull.getInitSimplex(buf);

   //find the points with the maximum value along each dimension
   //    //this piece of code uses the first d points with the maximal values along each dimension
   //for (IntInt_Iter=NonResultPtID.begin();IntInt_Iter!=NonResultPtID.end();IntInt_Iter++)
   //{
     //   long Pos=IntInt_Iter->first;
        //for (i=0;i<dimen;i++)
        //{
          //   pt[i]=PG[Pos][i]+SIDELEN;
           //  RankNonResult[i].insert(DblInt_Pair(pt[i],Pos));
        //}
   //}
   ///

   // //this piece of code uses the first d points that have the highest score in non-result points
   //for (IntInt_Iter=NonResultPtID.begin();IntInt_Iter!=NonResultPtID.end();IntInt_Iter++)
   for (i=0;i<NUMNONSKY;i++)
   {
        long Pos=nonSky10[i];
        for (j=0;j<dimen;j++)
        {
             pt[j]=PG[Pos][j]+SIDELEN;
             RankNonResult[j].insert(DblInt_Pair(pt[j],Pos));
        }
   }
   ///

   for (i=0;i<dimen;i++)
   {
        //printf("size of the %d-th map: %d\n",i,RankSkyline[i].size());
        DblInt_rIter=RankNonResult[i].rbegin();
        int szSet=FirstDpoints.size();
        while (DblInt_rIter!=RankNonResult[i].rend())
        {
            FirstDpoints.insert(DblInt_rIter->second);
            if (szSet==FirstDpoints.size())
                DblInt_rIter++;
            else
                break;
        }
   }

   //use the first d data points to compute the initial convex hull
   for (Int_sIter=FirstDpoints.begin();Int_sIter!=FirstDpoints.end();Int_sIter++)
   {
        //
        //printf("FirstDpoints: %ld\n",(*Int_sIter));
        //strcpy(buf,"");
        //for (i=0;i<dimen;i++)
        //{
          //   sprintf(buf+strlen(buf),"%lf ",PG[(*Int_sIter)][i]+SIDELEN);
                       //strcat(buf," ");
        //}
        ///
        //sprintf(buf+strlen(buf),"\n");

        strcpy(buf,PGstr[(*Int_sIter)]);

        //use the point for building initial simplex
        myHull.getInitSimplex(buf);

        //
        //build the mapping from id of partial hull vertex to original id of the points
        //HullIDtoPGID[IDMapPtr]=new long [2];
        //HullIDtoPGID[IDMapPtr][0]=myHull.numInitSitesRead-1;   //the id of partial hull vertex
        //HullIDtoPGID[IDMapPtr][1]=(*Int_sIter);   //the original id of the point
        //IDMapPtr++;
        ///
   }
   if (myHull.numInitSitesRead<dimen+1) cout << "Caution! there are less than d+1 points for building the initial hull!" << endl;


    //-------------------------------------------------------------------------
    // create a Rtree based on TGS
    //-------------------------------------------------------------------------
    long NumOfNonResultPt=0;
    float cl[MAXDIMEN], cu[MAXDIMEN];
    RtreeNodeEntry** p = new RtreeNodeEntry*[NonResultPtID.size()];
    for (IntInt_Iter=NonResultPtID.begin();IntInt_Iter!=NonResultPtID.end();IntInt_Iter++)
    {
        long id=IntInt_Iter->first;
        Int_sIter=FirstDpoints.find(id);
        if (Int_sIter!=FirstDpoints.end()) continue;

        for (int d=0; d<dimen; d++)
        {
             cl[d]=PG[id][d];
             cu[d]=PG[id][dimen+d];
        }
        Hypercube hc(dimen, cl, cu);
        p[NumOfNonResultPt++] = new RtreeNodeEntry(id, hc);
    }

    cout << "bulkloading R-tree... " << endl;;
    const int maxChild =
        (2048 - RtreeNode::size()) /
        RtreeNodeEntry::size(dimen);            // no. of entries per node
    //MainMemory mem(pagesize);
    FileMemory mem(2048, "tmpIdx", RtreeNodeEntry::fromMem, true);
    Rtree* rtreeNon =
        TGS::bulkload(mem, dimen, maxChild, maxChild,
        (int)(maxChild*0.3), (int)(maxChild*0.3), p, NonResultPtID.size()-FirstDpoints.size(), false);

    cout << "leaf node cnt: "<< rtreeNon->nodeCount(0) << endl;
    cout << " Max leaf fanout: " << rtreeNon->m_maxLeafChild << endl;
    cout << " Max node fanout: " << rtreeNon->m_maxNodeChild << endl;


   //begin explore the R-tree of non-result points, in order to construct the initial facets
   H0.insert(DblInt_Pair(0,rtreeNon->m_memory.m_rootPageID));

   //begin explore the min-heap to refine the initial facets adjacent to the k-th result
   float mindist;
   int MBRisInsideHull;
   bool isAnObject;
   int NumOfNonResultInHull=0;
   float ORIGIN[DMAX];

   for (int i=0;i<dimen;i++)
        ORIGIN[i]=CEILING+20;

   while (H0.size()!=0)
   {
        isAnObject=false;

    	DblInt_Iter=H0.begin();
    	long int pageid = DblInt_Iter->second;
        float ScoreTmp= DblInt_Iter->first;

    	H0.erase(DblInt_Iter);

    	VirtualRNode* VirNode=new VirtualRNode;
        RtreeNodeEntry* e0;            // create node entry e0 for node n, so that its MBR can be obtained
        RtreeNode* node=rtreeNon->m_memory.loadPage(pageid);
        VirNode->copyData(*node);
        VirNode->copyEntries(*node,node->m_usedSpace);
        e0 = node->genNodeEntry();     //compute the enclosing MBR for current index node
        delete node;

        //cout << "exam page " << pageid << s", Heap size: " << H0.size() << " mindist:" << DistTmp << endl;

        for (int i=0;i<dimen;i++)
             pt[i]=VirNode->m_entry[0]->m_hc.getUpper()[i];

        if (isAnObject==true)
        {
            strcpy(buf,PGstr[(VirNode->m_entry[0]->m_id)]);
        }
        else
        {
            //begin updating the partial hull by examining MBR and data points
            strcpy(buf,"");
            for (int i=0;i<dimen;i++)
            {
                 //pt[i]=VirNode->m_entry[0]->m_hc.getUpper()[i];

                 sprintf(buf+strlen(buf),"%lf ",pt[i]);
                 strcat(buf," ");
            }
        }
        MBRisInsideHull=myHull.isMBRVisible(buf);

        //
        //cout << "MBR "<< VirNode->m_entry[0]->m_id << ": "<< buf;
        //if (MBRisInsideHull==IS_VISIBLE)
          //  cout << " " << "visible" << endl;
        //else
          //  cout << " " << " NOT visible" << endl;
        ///

        if (MBRisInsideHull==IS_VISIBLE)
        {
           //check whether current entry is leaf node or index node
            if (VirNode->isLeaf())
            {
        	    if (VirNode->m_usedSpace>1)   //current node is a leaf node, so all its data entries must be examined
        	    {
                        for (int i=0; i<VirNode->m_usedSpace; i++)
                        {
                             //
                             //strcpy(buf,"");
                             //for (int j=0;j<dimen;j++)
                             //{
                               //   pt[j]=(VirNode->m_entry[i]->m_hc.getUpper()[j])-SIDELEN;

                                  //sprintf(buf+strlen(buf),"%lf ",pt[j]);
                                  //strcat(buf," ");
                             //}
                             ///

                             strcpy(buf,PGstr[(VirNode->m_entry[i]->m_id)]);

                             //feed the point to partial hull
                             myHull.buildPartialHull(buf);
                             if (myHull.isProned==true) NumOfNonResultInHull++;
                        }
        	    }
        	    else    //if current entry is a data point, then update the initial facets accordingly
                    {
                             //
                             //strcpy(buf,"");
                             //for (int j=0;j<dimen;j++)
                             //{
                                  //pt[j]=pt[j]-SIDELEN;   //obtain the real values of a data entry
                               //   pt[j]=VirNode->m_entry[0]->m_hc.getUpper()[j]-SIDELEN;

                                  //sprintf(buf+strlen(buf),"%lf ",pt[j]);
                                  //strcat(buf," ");
                             //}
                             ///

                             strcpy(buf,PGstr[(VirNode->m_entry[0]->m_id)]);

                             //feed the point to partial hull
                             myHull.buildPartialHull(buf);
                             if (myHull.isProned==true) NumOfNonResultInHull++;
                    }
            }
            else   //current entry is an index, so all its entries are inserted into the min-heap
            {
                  for (int i=0; i<VirNode->m_usedSpace; i++)
                  {
                       //check whether current index node is inside the partial hull or not
                       strcpy(buf,"");
                       for (int j=0;j<dimen;j++)
                       {
                            pt[j]=VirNode->m_entry[i]->m_hc.getUpper()[j];   //obtain the top-right corner point of current MBR

                            //sprintf(buf+strlen(buf),"%lf ",pt[j]);
                            //strcat(buf," ");
                       }

                       //
                       //int isInsideHull=myHull.isMBRVisible(buf);
                       //if (isInsideHull!=IS_VISIBLE)
                       //{
                         //  NumOfNodesInHull++;
                           //continue;
                       //}
                       ///


                       mindist=minDist(pt,ORIGIN,dimen);
                       H0.insert(DblInt_Pair(mindist,VirNode->m_entry[i]->m_id));
                  }
            }//end check whether it's leaf or index node
        }

        delete VirNode;
        delete e0;
   }


   //
   //the initial simplex has been built; preceed to compute the facets
   //adjacent to the k-th result (i.e., maintain the partial convex hull)
//   for (IntInt_Iter=NonResultPtID.begin();IntInt_Iter!=NonResultPtID.end();IntInt_Iter++) //this version is the original access order of points
  // {
    //      long Pos=IntInt_Iter->first;        //this version is the original access order of points
      //    Int_sIter=FirstDpoints.find(Pos);
      //    if (Int_sIter!=FirstDpoints.end()) continue;
      //    //sprintf(buf+strlen(buf),"\n");
      //    //printf("%s\n",buf);

      //    strcpy(buf,PGstr[Pos]);

      //    //build the mapping from id of partial hull vertex to original id of the points
      //    //HullIDtoPGID[IDMapPtr]=new long [2];

//          if (myHull.cdim<myHull.rdim)
  //        {
      //        //cout << "feeding " << buf << " to getInitSimplex..." << endl;
    //          myHull.getInitSimplex(buf);
        //      //HullIDtoPGID[IDMapPtr][0]=myHull.numInitSitesRead-1;   //the id of partial hull vertex
     //     }
     //     else
     //     {
      //        myHull.buildPartialHull(buf);
      //        //HullIDtoPGID[IDMapPtr][0]=myHull.pnum-2;   //the id of partial hull vertex
        //  }
          //if (myHull.isProned==true) numOfPtsInsideHull++;


          //build the mapping from id of partial hull vertex to original id of the points
          ////HullIDtoPGID[IDMapPtr][1]=Pos;   //the original id of the point
          ////IDMapPtr++;
          ////
   //}
   ///

   //cout << "number of non-results: " << NonResultPtID.size() << endl;
   cout << "number of non-results proned: " << NumOfNonResultInHull << endl;
   //cout << "The initial convex hull contains " << myHull.numInitSitesRead << " points." << endl;
   //cout << "Step 1 finished!" << endl;

    return 1;
}
//*/

int S3::window(
        Rtree &a_rtree,
        float *PG[],
        const Point &a_pt,
        const float *a_l,
        std::vector<long int> &a_resultID,
        int &a_maxStackSize,
        Array &a_pageaccessed) {
    multimap<long int, VirtualRNode *> DataEntryQ;  //queue storing rtree node which is simply a data entry
    multimap<long int, VirtualRNode *>::iterator IntVRN_Iter;
    typedef multimap<long int, VirtualRNode *>::value_type IntVRN_Pair;

    multimap<float, long int>::iterator DblInt_Iter;
    typedef multimap<float, long int>::value_type DblInt_Pair;
    typedef multimap<float, float>::value_type DblDbl_Pair;

    a_maxStackSize = 0;
    vector<long int> H;
    bool isAnObject;
    float pt[DMAX];

    H.push_back(a_rtree.m_memory.m_rootPageID);
    while (!H.empty()) {
        a_maxStackSize = H.size() > a_maxStackSize ? H.size() : a_maxStackSize;
        long int pageid = H.back();
        H.pop_back();

        isAnObject = false;

        VirtualRNode *VirNode = new VirtualRNode;
        RtreeNodeEntry *e0;// create node entry e0 for node n, so that its MBR can be obtained
        if (pageid >= MAXPAGEID) //current element in H is a data entry node (its pageid equals to m_id+MAXPAGED)
        {
            isAnObject = true;

            IntVRN_Iter = DataEntryQ.find(pageid);
            if (IntVRN_Iter == DataEntryQ.end()) {
                cout << "Error! there is no node " << pageid << " in DataEntryQ!" << endl;
            } else {
                VirNode->copyData(IntVRN_Iter->second);

                if (VirNode->m_usedSpace > 1) {
                    cout << "There is a fatal error in retrieving data-entry node from queue!" << endl;
                    return -1;
                } else
                    e0 = VirNode->m_entry[0]->clone();  //retrieve the sole entry from the data-entry node

                delete IntVRN_Iter->second; //free the virtualRnode; Note, this operation is VERY important
                DataEntryQ.erase(IntVRN_Iter);
            }
            pageid = pageid - MAXPAGEID;
        } else {
            RtreeNode *node = a_rtree.m_memory.loadPage(pageid);

            //my code
            VirNode->copyData(*node);
            VirNode->copyEntries(*node, node->m_usedSpace);
            e0 = node->genNodeEntry(); //compute the enclosing MBR for current index node

            delete node;
        }
        //cout << "examining page with id=" << pageid << endl;

        if (VirNode->isLeaf()) {
            if (VirNode->m_usedSpace >
                1)   //current node is a leaf node, so all its data entries must be put into priority queue H
            {
                a_pageaccessed.append((void *) pageid);
                for (int i = 0; i < VirNode->m_usedSpace; i++) {
                    long int NegPageid = VirNode->m_entry[i]->m_id + MAXPAGEID;
                    VirtualRNode *node = new VirtualRNode; //for each data entries of n, insert it into data-entry queue

                    RtreeNodeEntry *Nentry = VirNode->m_entry[i]->clone();
                    node->insertEntry(Nentry);
                    DataEntryQ.insert(IntVRN_Pair(NegPageid, node));
                    delete Nentry;

                    H.push_back(NegPageid);
                }
            } else    //current node is a data entry node
            {
                bool covered = true;
                long int id = VirNode->m_entry[0]->m_id;
                if (id != pageid) cout << "error!! page ids are not equal!" << endl;

                float rd[2] = {(PG[id][0] + PG[id][2]) / 2.0, (PG[id][1] + PG[id][3]) / 2.0};

                for (int d = 0; d < a_rtree.m_dimen && covered; d++)
                    if (fabs(rd[d] - a_pt[d]) > a_l[d]) covered = false;
                if (covered)
                    a_resultID.push_back(id);
            }
        } else//current entry is an index, so insert its entries into the heap
        {
            a_pageaccessed.append((void *) pageid);
            for (int i = 0; i < VirNode->m_usedSpace; i++) {
                bool covered = true;
                for (int d = 0; d < a_rtree.m_dimen && covered; d++)
                    if (fabs(VirNode->m_entry[i]->m_hc.mindist(a_pt, d)) > a_l[d])
                        covered = false;
                if (covered)
                    H.push_back(VirNode->m_entry[i]->m_id);
            }
        }//end check whether it's leaf or index node

        delete VirNode;
        delete e0;
    }

    return a_resultID.size();
}
//*/

int S3::incomparableWindows(
        Rtree &a_rtree,
        float *PG[],
        const Point &a_pt,
        const float *a_l,
        std::multimap<long int, VirtualRNode *> &NonResultEntry,  //non-result set, storing rtree node which is simply a data entry
        std::vector<long int> &a_resultID,    //the index nodes
        int &a_maxStackSize,
        Array &a_pageaccessed) {
    multimap<long int, VirtualRNode *> DataEntryQ;  //queue storing rtree node which is simply a data entry
    multimap<long int, VirtualRNode *>::iterator IntVRN_Iter;
    typedef multimap<long int, VirtualRNode *>::value_type IntVRN_Pair;

    multimap<float, long int>::iterator DblInt_Iter;
    typedef multimap<float, long int>::value_type DblInt_Pair;
    typedef multimap<float, float>::value_type DblDbl_Pair;

    a_maxStackSize = 0;
    vector<long int> H;
    bool isAnObject;
    float pt[DMAX];

    //generate a hypercube for current search quadrant w.r.t query point q
    float cl[DMAX], cu[DMAX];
    cl[0] = fabs(a_pt[0] - a_l[0]);
    cu[0] = a_pt[0] + a_l[0];
    cl[1] = fabs(a_pt[1] - a_l[1]);
    cu[1] = a_pt[1] + a_l[1];
    Hypercube q_hc(a_rtree.m_dimen, cl, cu);

    H.push_back(a_rtree.m_memory.m_rootPageID);
    long int round = 1;
    while (!H.empty()) {
        a_maxStackSize = H.size() > a_maxStackSize ? H.size() : a_maxStackSize;
        long int pageid = H.back();
        H.pop_back();

        isAnObject = false;

        VirtualRNode *VirNode = new VirtualRNode;
        RtreeNodeEntry *e0;// create node entry e0 for node n, so that its MBR can be obtained
        if (pageid >= MAXPAGEID) //current element in H is a data entry node (its pageid equals to m_id+MAXPAGED)
        {
            isAnObject = true;

            IntVRN_Iter = DataEntryQ.find(pageid);
            if (IntVRN_Iter == DataEntryQ.end()) {
                cout << "Error! there is no node " << pageid << " in DataEntryQ!" << endl;
            } else {
                VirNode->copyData(IntVRN_Iter->second);

                if (VirNode->m_usedSpace > 1) {
                    cout << "There is a fatal error in retrieving data-entry node from queue!" << endl;
                    return -1;
                } else
                    e0 = VirNode->m_entry[0]->clone();  //retrieve the sole entry from the data-entry node

                delete IntVRN_Iter->second; //free the virtualRnode; Note, this operation is VERY important
                DataEntryQ.erase(IntVRN_Iter);
            }
            pageid = pageid - MAXPAGEID;
        } else {
            RtreeNode *node = a_rtree.m_memory.loadPage(pageid);

            //my code
            VirNode->copyData(*node);
            VirNode->copyEntries(*node, node->m_usedSpace);
            e0 = node->genNodeEntry(); //compute the enclosing MBR for current index node

            delete node;
        }
        //cout << "examining page with id=" << pageid << endl;

        if (VirNode->isLeaf()) {
            if (VirNode->m_usedSpace >
                1)   //current node is a leaf node, so all its data entries must be put into priority queue H
            {
                bool covered = false;
                bool intersected = false;
                if (q_hc.enclose(e0->m_hc) == true)
                    covered = true;
                else if (q_hc.isIntersected(q_hc, e0->m_hc) == true)
                    intersected = true;

                if (intersected) {
                    a_pageaccessed.append((void *) pageid);
                    for (int i = 0; i < VirNode->m_usedSpace; i++) {
                        long int NegPageid = VirNode->m_entry[i]->m_id + MAXPAGEID;
                        VirtualRNode *node = new VirtualRNode; //for each data entries of n, insert it into data-entry queue

                        RtreeNodeEntry *Nentry = VirNode->m_entry[i]->clone();
                        node->insertEntry(Nentry);
                        DataEntryQ.insert(IntVRN_Pair(NegPageid, node));
                        delete Nentry;

                        //cout << "Pushing record " << NegPageid-MAXPAGEID << " into search heap." << endl;
                        H.push_back(NegPageid);
                    }
                }
                if (covered) {
                    a_resultID.push_back(pageid);  //keep the id of index node
                    //cout << "Pushing leaf " << pageid << " into result index set." << endl;
                }
            } else    //current node is a data entry node
            {
                bool covered = false;
                long int id = VirNode->m_entry[0]->m_id;
                if (id != pageid) cout << "error!! page ids are not equal!" << endl;

                float rd[2] = {(PG[id][0] + PG[id][2]) / 2.0, (PG[id][1] + PG[id][3]) / 2.0};

                //cout << "Record "<< id <<":"<< rd[0] << ", " << rd[1] << endl;
                //cout << "q:" << a_pt[0] << ", " << a_pt[1] << endl;
                //cout << "len:" << a_l[0]<< ", " << a_l[1] << endl;

                //for (int d=0; d<a_rtree.m_dimen && covered; d++)
                //     if (fabs(rd[d]-a_pt[d]) > a_l[d]) covered = false;

                Point tmpPt(a_rtree.m_dimen, rd);
                if (q_hc.enclose(tmpPt) == true) covered = true;
                if (covered) {
                    VirtualRNode *vNode = new VirtualRNode;
                    vNode->copyData(VirNode);
                    NonResultEntry.insert(IntVRN_Pair(id + MAXPAGEID, vNode));
                    //cout << "pushing record " << id << " into the result set." << endl;
                }
            }
        } else//current entry is an index, so insert its entries into the heap
        {
            bool covered = false;
            bool intersected = false;

            if (q_hc.enclose(e0->m_hc) == true)
                covered = true;
            else if (q_hc.isIntersected(q_hc, e0->m_hc) == true)
                intersected = true;

            if (intersected) {
                a_pageaccessed.append((void *) pageid);
                for (int i = 0; i < VirNode->m_usedSpace; i++) {
                    bool covered = true;
                    for (int d = 0; d < a_rtree.m_dimen && covered; d++)
                        if (fabs(VirNode->m_entry[i]->m_hc.mindist(a_pt, d)) > a_l[d])
                            covered = false;
                    if (covered) {
                        H.push_back(VirNode->m_entry[i]->m_id);
                        //cout << "pushing entry " << VirNode->m_entry[i]->m_id << " into the search heap." << endl;
                    }
                }
            }
            if (covered) {
                a_resultID.push_back(pageid);  //keep the id of index node
                //cout << "pushing node " << pageid << " into the result index set." << endl;
            }
        }//end check whether it's leaf or index node
        round++;

        delete VirNode;
        delete e0;
    }

    return a_resultID.size();
}
//*/

int S3::incomparableWindowsHD(
        Rtree &a_rtree,
        float *PG[],
        Hypercube &Dominee_hc,
        Hypercube &Dominator_hc,
        std::multimap<long int, VirtualRNode *> &NonResultEntry,  //non-result set, storing rtree node which is simply a data entry
        std::vector<long int> &a_resultID,    //the index nodes
        int &a_maxStackSize,
        Array &a_pageaccessed) {
    multimap<long int, VirtualRNode *> DataEntryQ;  //queue storing rtree node which is simply a data entry
    multimap<long int, VirtualRNode *>::iterator IntVRN_Iter;
    typedef multimap<long int, VirtualRNode *>::value_type IntVRN_Pair;

    multimap<float, long int>::iterator DblInt_Iter;
    typedef multimap<float, long int>::value_type DblInt_Pair;
    typedef multimap<float, float>::value_type DblDbl_Pair;

    a_maxStackSize = 0;
    vector<long int> H;
    bool isAnObject;
    float pt[DMAX];

    int dimen = a_rtree.m_dimen;

    long int NoOfDominators = 0;

    H.push_back(a_rtree.m_memory.m_rootPageID);
    long int round = 1;
    while (!H.empty()) {
        a_maxStackSize = H.size() > a_maxStackSize ? H.size() : a_maxStackSize;
        long int pageid = H.back();
        H.pop_back();

        isAnObject = false;

        VirtualRNode *VirNode = new VirtualRNode;
        RtreeNodeEntry *e0;// create node entry e0 for node n, so that its MBR can be obtained
        if (pageid >= MAXPAGEID) //current element in H is a data entry node (its pageid equals to m_id+MAXPAGED)
        {
            isAnObject = true;

            IntVRN_Iter = DataEntryQ.find(pageid);
            if (IntVRN_Iter == DataEntryQ.end()) {
                cout << "Error! there is no node " << pageid << " in DataEntryQ!" << endl;
            } else {
                VirNode->copyData(IntVRN_Iter->second);

                if (VirNode->m_usedSpace > 1) {
                    cout << "There is a fatal error in retrieving data-entry node from queue!" << endl;
                    return -1;
                } else
                    e0 = VirNode->m_entry[0]->clone();  //retrieve the sole entry from the data-entry node

                delete IntVRN_Iter->second; //free the virtualRnode; Note, this operation is VERY important
                DataEntryQ.erase(IntVRN_Iter);
            }
            pageid = pageid - MAXPAGEID;
        } else {
            RtreeNode *node = a_rtree.m_memory.loadPage(pageid);

            //my code
            VirNode->copyData(*node);
            VirNode->copyEntries(*node, node->m_usedSpace);
            e0 = node->genNodeEntry(); //compute the enclosing MBR for current index node

            delete node;
        }
        //cout << "examining page with id=" << pageid << endl;

        if (VirNode->isLeaf()) {
            if (VirNode->m_usedSpace >
                1)//current node is a leaf node, so all its data entries must be put into priority queue H
            {
                bool intersected = false;
                bool inDomineeWindow = false;
                bool inDominatorWindow = false;
                if (Dominator_hc.enclose(e0->m_hc) == true) {
                    inDominatorWindow = true;
                    NoOfDominators = NoOfDominators + VirNode->m_usedSpace;
                } else if (Dominator_hc.isIntersected(Dominator_hc, e0->m_hc) == true) {
                    intersected = true;
                }
                if (Dominee_hc.enclose(e0->m_hc) == true) {
                    inDomineeWindow = true;
                } else if (Dominee_hc.isIntersected(Dominee_hc, e0->m_hc) == true) {
                    intersected = true;
                }

                if (intersected) {
                    a_pageaccessed.append((void *) pageid);
                    for (int i = 0; i < VirNode->m_usedSpace; i++) {
                        long int NegPageid = VirNode->m_entry[i]->m_id + MAXPAGEID;
                        VirtualRNode *node = new VirtualRNode; //for each data entries of n, insert it into data-entry queue

                        RtreeNodeEntry *Nentry = VirNode->m_entry[i]->clone();
                        node->insertEntry(Nentry);
                        DataEntryQ.insert(IntVRN_Pair(NegPageid, node));
                        delete Nentry;

                        //cout << "Pushing record " << NegPageid-MAXPAGEID << " into search heap." << endl;
                        H.push_back(NegPageid);
                    }
                } else if (!inDomineeWindow && !inDominatorWindow) {
                    a_resultID.push_back(pageid);  //keep the id of index node
                    //cout << "Pushing leaf " << pageid << " into result index set." << endl;
                }
            } else    //current node is a data entry node
            {
                long int id = VirNode->m_entry[0]->m_id;
                if (id != pageid) cout << "error!! page ids are not equal!" << endl;

                float rd[DMAX];
                for (int d = 0; d < dimen; d++)
                    rd[d] = (PG[id][d] + PG[id][dimen + d]) / 2.0;

                Point tmpPt(dimen, rd);
                if (Dominator_hc.enclose(tmpPt) == true)   //current point lies in the dominator window
                {
                    NoOfDominators++;
                } else {
                    if (Dominee_hc.enclose(tmpPt) == false)  //current point lies in some incomparable window
                    {
                        VirtualRNode *vNode = new VirtualRNode;
                        vNode->copyData(VirNode);
                        NonResultEntry.insert(IntVRN_Pair(id + MAXPAGEID, vNode));
                        //cout << "pushing record " << id << " into the result set." << endl;
                    }
                }
            }
        } else//current entry is an index, so insert its entries into the heap
        {
            bool intersected = false;
            bool inDomineeWindow = false;
            bool inDominatorWindow = false;

            if (Dominee_hc.enclose(e0->m_hc) == true)
                inDomineeWindow = true;
            else if (Dominee_hc.isIntersected(Dominee_hc, e0->m_hc) == true)
                intersected = true;
            if (Dominator_hc.enclose(e0->m_hc) == true)
                inDominatorWindow = true;
            else if (Dominator_hc.isIntersected(Dominator_hc, e0->m_hc) == true)
                intersected = true;

            if (intersected) {
                a_pageaccessed.append((void *) pageid);
                for (int i = 0; i < VirNode->m_usedSpace; i++) {
                    if (Dominee_hc.enclose(VirNode->m_entry[i]->m_hc) == false) {
                        H.push_back(VirNode->m_entry[i]->m_id);
                    }
                }
            } else if (!inDomineeWindow && !inDominatorWindow) {
                a_resultID.push_back(pageid);  //keep the id of index node
            }
        }//end check whether it's leaf or index node
        round++;

        delete VirNode;
        delete e0;
    }

    return NoOfDominators;
}
//*/

int S3::getIncomparableRecords(
        Rtree &a_rtree,
        float *PG[],
        Hypercube &Dominee_hc,
        Hypercube &Dominator_hc,
        //std::multimap<long int, VirtualRNode*>& NonResultEntry,  //non-result set, storing rtree node which is simply a data entry
        std::vector<long int> &a_resultID,    //IDs of incomparable records
        int &a_maxStackSize,
        Array &a_pageaccessed) {
    multimap<long int, VirtualRNode *> DataEntryQ;  //queue storing rtree node which is simply a data entry
    multimap<long int, VirtualRNode *>::iterator IntVRN_Iter;
    typedef multimap<long int, VirtualRNode *>::value_type IntVRN_Pair;

    multimap<float, long int>::iterator DblInt_Iter;
    typedef multimap<float, long int>::value_type DblInt_Pair;
    typedef multimap<float, float>::value_type DblDbl_Pair;

    a_maxStackSize = 0;
    vector<long int> H;
    bool isAnObject;
    float pt[DMAX];

    int dimen = a_rtree.m_dimen;

    long int NoOfDominators = 0;

    H.push_back(a_rtree.m_memory.m_rootPageID);
    long int round = 1;
    while (!H.empty()) {
        a_maxStackSize = H.size() > a_maxStackSize ? H.size() : a_maxStackSize;
        long int pageid = H.back();
        H.pop_back();

        isAnObject = false;

        VirtualRNode *VirNode = new VirtualRNode;
        RtreeNodeEntry *e0;// create node entry e0 for node n, so that its MBR can be obtained
        if (pageid >= MAXPAGEID) //current element in H is a data entry node (its pageid equals to m_id+MAXPAGED)
        {
            isAnObject = true;

            IntVRN_Iter = DataEntryQ.find(pageid);
            if (IntVRN_Iter == DataEntryQ.end()) {
                cout << "Error! there is no node " << pageid << " in DataEntryQ!" << endl;
            } else {
                VirNode->copyData(IntVRN_Iter->second);

                if (VirNode->m_usedSpace > 1) {
                    cout << "There is a fatal error in retrieving data-entry node from queue!" << endl;
                    return -1;
                } else
                    e0 = VirNode->m_entry[0]->clone();  //retrieve the sole entry from the data-entry node

                delete IntVRN_Iter->second; //free the virtualRnode; Note, this operation is VERY important
                DataEntryQ.erase(IntVRN_Iter);
            }
            pageid = pageid - MAXPAGEID;
        } else {
            RtreeNode *node = a_rtree.m_memory.loadPage(pageid);

            //my code
            VirNode->copyData(*node);
            VirNode->copyEntries(*node, node->m_usedSpace);
            e0 = node->genNodeEntry(); //compute the enclosing MBR for current index node

            delete node;
        }
        //cout << "examining page with id=" << pageid << endl;

        if (VirNode->isLeaf()) {
            if (VirNode->m_usedSpace >
                1)//current node is a leaf node, so all its data entries must be put into priority queue H
            {
                bool intersected = false;
                bool inDomineeWindow = false;
                bool inDominatorWindow = false;
                if (Dominator_hc.enclose(e0->m_hc) == true) {
                    inDominatorWindow = true;
                    NoOfDominators = NoOfDominators + VirNode->m_usedSpace;
                } else if (Dominator_hc.isIntersected(Dominator_hc, e0->m_hc) == true) {
                    intersected = true;
                }
                if (Dominee_hc.enclose(e0->m_hc) == true) {
                    inDomineeWindow = true;
                } else if (Dominee_hc.isIntersected(Dominee_hc, e0->m_hc) == true) {
                    intersected = true;
                }

                if (intersected || (!inDomineeWindow && !inDominatorWindow)) {
                    a_pageaccessed.append((void *) pageid);
                    for (int i = 0; i < VirNode->m_usedSpace; i++) {
                        long int NegPageid = VirNode->m_entry[i]->m_id + MAXPAGEID;
                        VirtualRNode *node = new VirtualRNode; //for each data entries of n, insert it into data-entry queue

                        RtreeNodeEntry *Nentry = VirNode->m_entry[i]->clone();
                        node->insertEntry(Nentry);
                        DataEntryQ.insert(IntVRN_Pair(NegPageid, node));
                        delete Nentry;

                        //cout << "Pushing record " << NegPageid-MAXPAGEID << " into search heap." << endl;
                        H.push_back(NegPageid);
                    }
                }
            } else    //current node is a data entry node
            {
                long int id = VirNode->m_entry[0]->m_id;
                if (id != pageid) cout << "error!! page ids are not equal!" << endl;

                float rd[DMAX];
                for (int d = 0; d < dimen; d++)
                    rd[d] = (PG[id][d] + PG[id][dimen + d]) / 2.0;

                Point tmpPt(dimen, rd);
                if (Dominator_hc.enclose(tmpPt) == true)   //current point lies in the dominator window
                {
                    NoOfDominators++;
                } else {
                    if (Dominee_hc.enclose(tmpPt) == false)  //current point lies in some incomparable window
                    {
                        a_resultID.push_back(id);
                        //cout << "pushing record " << id << " into the result set." << endl;
                    }
                }
            }
        } else//current entry is an index, so insert its entries into the heap
        {
            bool intersected = false;
            bool inDomineeWindow = false;
            bool inDominatorWindow = false;

            if (Dominee_hc.enclose(e0->m_hc) == true)
                inDomineeWindow = true;
            else if (Dominee_hc.isIntersected(Dominee_hc, e0->m_hc) == true)
                intersected = true;
            if (Dominator_hc.enclose(e0->m_hc) == true)
                inDominatorWindow = true;
            else if (Dominator_hc.isIntersected(Dominator_hc, e0->m_hc) == true)
                intersected = true;

            if (intersected || (!inDomineeWindow && !inDominatorWindow)) {
                a_pageaccessed.append((void *) pageid);
                for (int i = 0; i < VirNode->m_usedSpace; i++) {
                    if (Dominee_hc.enclose(VirNode->m_entry[i]->m_hc) == false) {
                        H.push_back(VirNode->m_entry[i]->m_id);
                    }
                }
            }
        }//end check whether it's leaf or index node
        round++;

        delete VirNode;
        delete e0;
    }

    return NoOfDominators;
}
//*/


typedef struct vdp {      //value-direction pair
    std::set<long int> HalflineL;
    std::set<long int> HalflineR;
} vdp;

int S3::AA_2D_old(      //Advanced Algorithm for 2-d
        const int dimen, //dimensionality
        Rtree &a_rtree,
        float *PG[],
        Point &a_pt,
        long int &totalNoOfminCells,
        int &maxStackSize, Array &a_pageaccessed) {
    multimap<float, vdp *> Cells;
    multimap<float, vdp *>::iterator c_Itr, c_Itr1, c_Itr2;
    typedef multimap<float, vdp *>::value_type c_VT;
    multimap<long int, float>::iterator IntDbl_Itr, IntDbl_Itr1;
    typedef multimap<long int, float>::value_type IntDbl_m_VT;

    multimap<float, float>::iterator DblDbl_Itr, DblDbl_Itr1;
    typedef multimap<float, float>::value_type DblDbl_VT;

    map<float, set<float> *> ResultCells;
    map<float, set<float> *>::iterator DblS_mItr, DblS_mItr1;
    typedef map<float, set<float> *>::value_type m_VT;
    set<float>::iterator Dbl_sItr, Dbl_sItr1;

    std::set<long int> Singular;
    std::pair<set<long int>::iterator, bool> rset;

    std::vector<long int> ResultID[3];
    vector<long int>::iterator vecItr, vecItr1;

    std::vector<long int> skyline;

    std::set<long int>::iterator sItr, sItr1;

    std::multimap<long int, VirtualRNode *> NonResultEntry;   //data nodes that are kept for incremental skyline computation

    std::map<float, long int> Intersections;
    std::map<long int, float> RdAndCrossPt;
    std::map<float, float>::iterator mDblDblItr;
    typedef map<float, float>::value_type DblDbl_Pair;
    std::map<long int, float>::iterator mDblItr;
    typedef map<long int, float>::value_type IntDbl_Pair;

    map<float, long int>::iterator mItr, mItr1;
    typedef map<float, long int>::value_type DblInt_Pair;

    float dataspace[4][2] = {{0, 0},
                             {1, 0},
                             {1, 1},
                             {0, 1}};
    float TMPcenter[DMAX];
    float a_l[DMAX];    //storing the ranges (w.r.t point pt) along each dimension of a quadrant

    long int min_k = INT_MAX;
    long int NoOfDominators = 0;
    long int NoOfDominees = 0;

    long int CntLeft = 0, CntRight = 0;
    bool FoundMink = false;
    float theMinCell[2];
    long int NoOfMinCells = 0;

    char direction;

    NonResultEntry.clear();
    for (int i = 0; i < 4; i++)  //for each of the four quadrants (defined by a_pt), perform a window query
    {
        TMPcenter[0] = float(dataspace[i][0] + a_pt[0]) / 2.0;
        TMPcenter[1] = float(dataspace[i][1] + a_pt[1]) / 2.0;
        a_l[0] = fabs(TMPcenter[0] - a_pt[0]);
        a_l[1] = fabs(TMPcenter[1] - a_pt[1]);

        //cout << "center:" << TMPcenter[0] << ", " << TMPcenter[1] << endl;
        //cout << "len:" << a_l[0] << ", " << a_l[1] << endl;

        Point center(dimen, TMPcenter);

        if (i == 0) {
            window(a_rtree, PG, center, a_l, ResultID[i], maxStackSize, a_pageaccessed);
            NoOfDominees = ResultID[i].size();
            continue;
        }

        if (i == 2) {
            window(a_rtree, PG, center, a_l, ResultID[i], maxStackSize, a_pageaccessed);
            NoOfDominators = ResultID[i].size();
            continue;
        }

        if (i != 0 && i != 2)   //process quadrants that contain incomparable records
        {
            incomparableWindows(a_rtree, PG, center, a_l, NonResultEntry, ResultID[1], maxStackSize, a_pageaccessed);
        }
    }
    cout << "NoOfDominators=" << NoOfDominators << ", NoOfDominees=" << NoOfDominees << endl;

    //incrementally update skylines and then compute the cells with minimal order
    std::set<long int> tmpSkyline;
    tmpSkyline.clear();

    //add in the two boundary half-lines, i.e., 0 and 1
    vdp *vdPair = new vdp;
    vdPair->HalflineL.clear();
    vdPair->HalflineR.clear();
    Cells.insert(c_VT(0, vdPair));
    vdPair = new vdp;
    vdPair->HalflineL.clear();
    vdPair->HalflineR.clear();
    Cells.insert(c_VT(1, vdPair));

    CntLeft = 0;
    CntRight = 0;
    while (!FoundMink) {
        //update the skyline by calling GetSkyline() with tmpSkyline containing the skyline points to delete
        //cout << "Number of NonResultEntry: "<< NonResultEntry.size() << ", Number of index nodes:"<< ResultID[1].size() << endl;
        GetSkylines(dimen, a_rtree, NonResultEntry, ResultID[1], tmpSkyline, PG, maxStackSize, a_pageaccessed);
        cout << "Number of skylines: " << tmpSkyline.size() << endl;

        //compute intersection point, i.e., the value-direction pair, for each newly updated skyline points, and insert it into Cells
        for (sItr = tmpSkyline.begin(); sItr != tmpSkyline.end(); sItr++) {
            long int id = (*sItr);
            float rd[2] = {(PG[id][0] + PG[id][2]) / 2.0, (PG[id][1] + PG[id][3]) / 2.0};

            int tmpDim = 0;
            for (int d = 0; d < dimen; d++)
                if (fabs(rd[d] - a_pt[d]) <= 1e-6) tmpDim++;
            if (tmpDim == dimen) {
                //cout << "CAUTION: there exist two exactly the same points!" << endl;  //i.e.,query point is a randomly picked record
                //getchar();
                continue;
            }

            cout << "processing skyline " << id << ": (" << rd[0] << "," << rd[1] << ")" << endl;

            float denominator = a_pt[1] - a_pt[1] + rd[0] - rd[1];
            float crosspoint;
            crosspoint = float(a_pt[1] - rd[1]) / denominator;
            cout << "denominator: " << denominator << ", crosspoint:" << crosspoint << " ";

            if (crosspoint < 0 || crosspoint > 1) {
                //cout << endl;
                continue;
            }

            if (rd[1] > a_pt[1]) cout << " <--" << endl;
            else cout << " -->" << endl;

            //determine the direction

            if (rd[1] > a_pt[1]) //the halfline lies at the left of the intersection point
                direction = 'L';   //half-line pointing to the left of the intersection point
            else            //the halfline lies at the right of the intersection point
                direction = 'R';   //half-line pointing to the right of the intersection point

            c_Itr = Cells.find(crosspoint);
            if (c_Itr != Cells.end()) {
                if (direction == 'L') {
                    rset = c_Itr->second->HalflineL.insert(id);
                    if (rset.second == true) CntLeft++;  //id is a newly inserted crosspoint
                } else if (direction == 'R') {
                    rset = c_Itr->second->HalflineR.insert(id);
                    if (rset.second == true) CntRight++;  //id is a newly inserted crosspoint
                }
            } else {
                vdp *vdPair = new vdp;
                if (direction == 'L') {
                    vdPair->HalflineL.insert(id);
                    CntLeft++;
                } else if (direction == 'R') {
                    vdPair->HalflineR.insert(id);
                    CntRight++;
                }
                Cells.insert(c_VT(crosspoint, vdPair));
            }
        }
        //cout << "Number of total interim cells: " << Cells.size()-1 << endl;
        //cout << "CntLeft=" << CntLeft << ", CntRight="<< CntRight << endl;
        //getchar();

        //Determine the cell(s) with the minimal order
        long int minOrder = INT_MAX, order;
        multimap<float, float> minCells;
        vector<float> prevLeftBoundary;
        long int currentCntL = CntLeft, currentCntR = 0;
        c_Itr = Cells.begin();
        c_Itr1 = c_Itr;
        c_Itr1++;
        prevLeftBoundary.push_back(0);
        for (; c_Itr1 != Cells.end();) {
            if (c_Itr1->second->HalflineR.size() > 0) {
                long int tmpNoOfSingulars = 0;
                for (sItr = c_Itr1->second->HalflineR.begin(); sItr != c_Itr1->second->HalflineR.end(); sItr++) {
                    sItr1 = Singular.find((*sItr));
                    if (sItr1 != Singular.end()) tmpNoOfSingulars++;
                }
                if (tmpNoOfSingulars == c_Itr1->second->HalflineR.size())
                    prevLeftBoundary.push_back(c_Itr1->first);
            }

            if (c_Itr->first == 0)  //the first cell
            {
                order = currentCntL + c_Itr->second->HalflineR.size();
                currentCntL = currentCntL - c_Itr1->second->HalflineL.size();
                currentCntR = currentCntR + c_Itr->second->HalflineR.size();
            } else if (c_Itr1->first == 1) //the last cell
            {
                order = CntRight + c_Itr1->second->HalflineL.size();
            } else {
                currentCntR = currentCntR + c_Itr->second->HalflineR.size();
                order = currentCntL + currentCntR;
                currentCntL = currentCntL - c_Itr1->second->HalflineL.size();
            }
            cout << "order: " << order << ", minOrder: " << minOrder << endl;
            if (order <= minOrder) {
                if (order < minOrder) minCells.clear();
                minOrder = order;

                //determine the two boundaries of current cell
                float tmpLeftBoundary = prevLeftBoundary.back();
                float tmpRightBoundary = 1;
                c_Itr2 = c_Itr1;
                while (c_Itr2 != Cells.end()) {
                    if (c_Itr2->second->HalflineL.size() > 0) {
                        long int tmpNoOfSingulars = 0;
                        for (sItr = c_Itr2->second->HalflineL.begin();
                             sItr != c_Itr2->second->HalflineL.end(); sItr++) {
                            sItr1 = Singular.find((*sItr));
                            if (sItr1 != Singular.end()) tmpNoOfSingulars++;
                        }
                        if (tmpNoOfSingulars == c_Itr2->second->HalflineL.size()) {
                            tmpRightBoundary = c_Itr2->first;
                            break;
                        }
                    }
                    c_Itr2++;
                }
                minCells.insert(DblDbl_VT(tmpLeftBoundary, tmpRightBoundary));
                cout << "tmpLeftBound: " << tmpLeftBoundary << ", tmpRightBound: " << tmpRightBoundary << endl;
            }
            c_Itr++;
            c_Itr1 = c_Itr;
            c_Itr1++;
        }
        //cout << "minOrder=" << minOrder << ", min_k=" << min_k << endl;
        //cout << "#minCells=" << minCells.size() << ", size of Singular="<< Singular.size() << endl;
        //getchar();

        //check the min_k with the current minimal order; if min_k is less than the currrent minimal order, stops;
        //otherwise, if the halflines of the minimal-ordered cell(s) are not singular, then expand the halflines
        if (minOrder > min_k) return min_k + NoOfDominators;

        tmpSkyline.clear();
        for (DblDbl_Itr = minCells.begin(); DblDbl_Itr != minCells.end(); DblDbl_Itr++) {
            //cout << "Left: " << DblDbl_Itr->first << ", Right: " << DblDbl_Itr->second << endl;
            //insert the bounding halflines of current minimal-ordered cell to expand
            c_Itr = Cells.find(DblDbl_Itr->first);
            c_Itr1 = Cells.find(DblDbl_Itr->second);
            if (c_Itr != Cells.end()) {
                long int NoOfSingulars = 0;
                long int NoOfHalfLines = 0;
                std::pair<set<long int>::iterator, bool> iRes;
                for (int i = 0; i < 2; i++) {
                    theMinCell[i] = c_Itr->first;

                    //cout << "processing min-cell with intersection=" << c_Itr->first << ", left sz=" << c_Itr->second->HalflineL.size() << ", right sz="<< c_Itr->second->HalflineR.size() << endl;
                    for (sItr = c_Itr->second->HalflineL.begin(); sItr != c_Itr->second->HalflineL.end(); sItr++) {
                        //cout << "processing left halfline "<< *sItr << endl;
                        iRes = Singular.insert(*sItr);
                        if (iRes.second == true)   //current halfline with id *sItr is augmented, so remove it from SL
                            tmpSkyline.insert(*sItr);
                        else
                            NoOfSingulars++;
                    }
                    for (sItr = c_Itr->second->HalflineR.begin(); sItr != c_Itr->second->HalflineR.end(); sItr++) {
                        //cout << "processing right halfline "<< *sItr << endl;
                        iRes = Singular.insert(*sItr);
                        if (iRes.second == true)   //current halfline with id *sItr is augmented, so remove it from SL
                            tmpSkyline.insert(*sItr);
                        else
                            NoOfSingulars++;
                    }
                    NoOfHalfLines = NoOfHalfLines + c_Itr->second->HalflineL.size() + c_Itr->second->HalflineR.size();
                    c_Itr = c_Itr1;
                }
                if (NoOfSingulars == NoOfHalfLines) {
                    //current cell has only singular halflines, so we fix its order and compare it with min_k
                    if (minOrder > min_k) {
                        FoundMink = true;
                        break;
                    } else {
                        min_k = minOrder;

                        DblS_mItr = ResultCells.find(theMinCell[0]);
                        if (DblS_mItr != ResultCells.end()) {
                            Dbl_sItr = (DblS_mItr->second)->find(theMinCell[1]);
                            if (Dbl_sItr != (DblS_mItr->second)->end()) {
                                FoundMink = true;
                                break;
                            } else {
                                (DblS_mItr->second)->insert(theMinCell[1]);
                                NoOfMinCells++;
                            }
                        } else {
                            set<float> *cellRbound = new set<float>;
                            cellRbound->insert(theMinCell[1]);
                            ResultCells.insert(m_VT(theMinCell[0], cellRbound));
                            NoOfMinCells++;
                        }
                    }
                }
            }
        }
        cout << "#skylines to remove: " << tmpSkyline.size() << endl;
        //for (sItr=tmpSkyline.begin();sItr!=tmpSkyline.end();sItr++)
        //     cout << *sItr << " ";
        //cout << endl;
    }

    cout << "minimal k= " << min_k << ", NoOfDominators=" << NoOfDominators << ", NoOfDominees=" << NoOfDominees
         << endl;
    cout << "#Cell with order " << min_k << ": " << NoOfMinCells << endl;

    totalNoOfminCells = totalNoOfminCells + NoOfMinCells;

    return min_k + NoOfDominators;
}

int S3::AA_2D(      //Advanced Algorithm for 2-d
        const int dimen, //dimensionality
        Rtree &a_rtree,
        float *PG[],
        Point &a_pt,
        long int &totalNoOfminCells,
        int &maxStackSize, Array &a_pageaccessed) {
    multimap<float, vdp *> Cells;
    multimap<float, vdp *>::iterator c_Itr, c_Itr1, c_Itr2;
    typedef multimap<float, vdp *>::value_type c_VT;
    multimap<long int, float>::iterator IntDbl_Itr, IntDbl_Itr1;
    typedef multimap<long int, float>::value_type IntDbl_m_VT;

    multimap<float, float>::iterator DblDbl_Itr, DblDbl_Itr1;
    typedef multimap<float, float>::value_type DblDbl_VT;

    set<float>::iterator Dbl_sItr, Dbl_sItr1;

    std::set<long int> Singular;
    std::pair<set<long int>::iterator, bool> rset;

    std::vector<long int> ResultID[3];
    vector<long int>::iterator vecItr, vecItr1;

    std::vector<long int> skyline;

    std::set<long int>::iterator sItr, sItr1;

    std::multimap<long int, VirtualRNode *> NonResultEntry;   //data nodes that are kept for incremental skyline computation

    std::map<float, long int> Intersections;
    std::map<long int, float> RdAndCrossPt;
    std::map<float, float>::iterator mDblDblItr;
    typedef map<float, float>::value_type DblDbl_Pair;
    std::map<long int, float>::iterator mDblItr;
    typedef map<long int, float>::value_type IntDbl_Pair;

    map<float, long int>::iterator mItr, mItr1;
    typedef map<float, long int>::value_type DblInt_Pair;

    float dataspace[4][2] = {{0, 0},
                             {1, 0},
                             {1, 1},
                             {0, 1}};
    float TMPcenter[DMAX];
    float a_l[DMAX];    //storing the ranges (w.r.t point pt) along each dimension of a quadrant

    long int min_k = INT_MAX;
    long int NoOfDominators = 0;
    long int NoOfDominees = 0;

    long int CntLeft = 0, CntRight = 0;
    bool FoundMink = false;
    float theMinCell[2];
    long int NoOfMinCells = 0;

    char direction;

    NonResultEntry.clear();
    for (int i = 0; i < 4; i++)  //for each of the four quadrants (defined by a_pt), perform a window query
    {
        TMPcenter[0] = float(dataspace[i][0] + a_pt[0]) / 2.0;
        TMPcenter[1] = float(dataspace[i][1] + a_pt[1]) / 2.0;
        a_l[0] = fabs(TMPcenter[0] - a_pt[0]);
        a_l[1] = fabs(TMPcenter[1] - a_pt[1]);

        //cout << "center:" << TMPcenter[0] << ", " << TMPcenter[1] << endl;
        //cout << "len:" << a_l[0] << ", " << a_l[1] << endl;

        Point center(dimen, TMPcenter);

        if (i == 0) {
            window(a_rtree, PG, center, a_l, ResultID[i], maxStackSize, a_pageaccessed);
            NoOfDominees = ResultID[i].size();
            continue;
        }

        if (i == 2) {
            window(a_rtree, PG, center, a_l, ResultID[i], maxStackSize, a_pageaccessed);
            NoOfDominators = ResultID[i].size();
            continue;
        }

        if (i != 0 && i != 2)   //process quadrants that contain incomparable records
        {
            incomparableWindows(a_rtree, PG, center, a_l, NonResultEntry, ResultID[1], maxStackSize, a_pageaccessed);
        }
    }
    cout << "NoOfDominators=" << NoOfDominators << ", NoOfDominees=" << NoOfDominees << endl;

    //incrementally update skylines and then compute the cells with minimal order
    std::set<long int> tmpSkyline;
    tmpSkyline.clear();

    //add in the two boundary half-lines, i.e., 0 and 1
    vdp *vdPair = new vdp;
    vdPair->HalflineL.clear();
    vdPair->HalflineR.clear();
    Cells.insert(c_VT(0, vdPair));
    vdPair = new vdp;
    vdPair->HalflineL.clear();
    vdPair->HalflineR.clear();
    Cells.insert(c_VT(1, vdPair));

    CntLeft = 0;
    CntRight = 0;
    bool firstTimeScan = true;
    while (!FoundMink) {
        //update the skyline by calling GetSkyline() with tmpSkyline containing the skyline points to delete
        //cout << "Number of NonResultEntry: "<< NonResultEntry.size() << ", Number of index nodes:"<< ResultID[1].size() << endl;
        GetSkylines(dimen, a_rtree, NonResultEntry, ResultID[1], tmpSkyline, PG, maxStackSize, a_pageaccessed);
        //cout << "Number of skylines: "<< tmpSkyline.size() << endl;

        //compute intersection point, i.e., the value-direction pair, for each newly updated skyline points, and insert it into Cells
        for (sItr = tmpSkyline.begin(); sItr != tmpSkyline.end(); sItr++) {
            long int id = (*sItr);
            float rd[2] = {(PG[id][0] + PG[id][2]) / 2.0, (PG[id][1] + PG[id][3]) / 2.0};

            int tmpDim = 0;
            for (int d = 0; d < dimen; d++)
                if (fabs(rd[d] - a_pt[d]) <= 1e-6) tmpDim++;
            if (tmpDim == dimen) {
                //cout << "CAUTION: there exist two exactly the same points!" << endl;  //i.e.,query point is a randomly picked record
                //getchar();
                continue;
            }

            //cout << "processing skyline " << id <<": (" << rd[0] << ","<< rd[1] <<")"<< endl;

            float denominator = a_pt[1] - a_pt[1] + rd[0] - rd[1];
            float crosspoint;
            crosspoint = float(a_pt[1] - rd[1]) / denominator;
            //cout << "denominator: " << denominator << ", crosspoint:" << crosspoint << " ";

            //if (crosspoint < 0 || crosspoint > 1)
            //{
            //cout << endl;
            //   continue;
            //}

            //if (rd[1]>a_pt[1]) cout << " <--" << endl;
            //else cout << " -->" << endl;

            //determine the direction

            if (rd[1] > a_pt[1]) //the halfline lies at the left of the intersection point
                direction = 'L';   //half-line pointing to the left of the intersection point
            else            //the halfline lies at the right of the intersection point
                direction = 'R';   //half-line pointing to the right of the intersection point

            c_Itr = Cells.find(crosspoint);
            if (c_Itr != Cells.end()) {
                if (direction == 'L') {
                    rset = c_Itr->second->HalflineL.insert(id);
                    if (rset.second == true) CntLeft++;  //id is a newly inserted crosspoint
                } else if (direction == 'R') {
                    rset = c_Itr->second->HalflineR.insert(id);
                    if (rset.second == true) CntRight++;  //id is a newly inserted crosspoint
                }
            } else {
                vdp *vdPair = new vdp;
                if (direction == 'L') {
                    vdPair->HalflineL.insert(id);
                    CntLeft++;
                } else if (direction == 'R') {
                    vdPair->HalflineR.insert(id);
                    CntRight++;
                }
                Cells.insert(c_VT(crosspoint, vdPair));
            }
        }
        //cout << "Number of total interim cells: " << Cells.size()-1 << endl;
        //cout << "CntLeft=" << CntLeft << ", CntRight="<< CntRight << endl;

        //Determine the cell(s) with the minimal order
        long int minOrder = INT_MAX, order;
        multimap<float, float> minCells;
        vector<float> prevLeftBoundary;
        long int currentCntL = CntLeft, currentCntR = 0;
        c_Itr = Cells.begin();
        c_Itr1 = c_Itr;
        c_Itr1++;
        prevLeftBoundary.push_back(0);
        for (; c_Itr1 != Cells.end();) {
            if (c_Itr1->second->HalflineR.size() > 0) {
                if (!firstTimeScan) {
                    long int tmpNoOfSingulars = 0;
                    for (sItr = c_Itr1->second->HalflineR.begin(); sItr != c_Itr1->second->HalflineR.end(); sItr++) {
                        sItr1 = Singular.find((*sItr));
                        if (sItr1 != Singular.end()) {
                            //cout << "Singular@Left: " << *sItr << endl;
                            tmpNoOfSingulars++;
                        }
                    }
                    if (tmpNoOfSingulars < c_Itr1->second->HalflineR.size())
                        prevLeftBoundary.push_back(c_Itr1->first);
                } else prevLeftBoundary.push_back(c_Itr1->first);
            }

            if (c_Itr->first == 0)  //the first cell
            {
                order = currentCntL + c_Itr->second->HalflineR.size();
                currentCntL = currentCntL - c_Itr1->second->HalflineL.size();
                currentCntR = currentCntR + c_Itr->second->HalflineR.size();
            } else if (c_Itr1->first == 1) //the last cell
            {
                order = CntRight + c_Itr1->second->HalflineL.size();
            } else {
                currentCntR = currentCntR + c_Itr->second->HalflineR.size();
                order = currentCntL + currentCntR;
                currentCntL = currentCntL - c_Itr1->second->HalflineL.size();
            }
            //cout << "order: " << order << ", minOrder: " << minOrder << endl;
            if (order <= minOrder) {
                if (order < minOrder) minCells.clear();
                minOrder = order;

                //determine the two boundaries of current cell
                float tmpLeftBoundary = prevLeftBoundary.back();
                float tmpRightBoundary = 1;
                c_Itr2 = c_Itr1;
                while (c_Itr2 != Cells.end()) {
                    if (c_Itr2->second->HalflineL.size() > 0) {
                        if (!firstTimeScan) {
                            long int tmpNoOfSingulars = 0;
                            for (sItr = c_Itr2->second->HalflineL.begin();
                                 sItr != c_Itr2->second->HalflineL.end(); sItr++) {
                                sItr1 = Singular.find((*sItr));
                                if (sItr1 != Singular.end()) {
                                    //cout << "Singular@Right: " << *sItr << endl;
                                    tmpNoOfSingulars++;
                                }
                            }
                            if (tmpNoOfSingulars < c_Itr2->second->HalflineL.size()) {
                                tmpRightBoundary = c_Itr2->first;
                                break;
                            }
                        } else {
                            tmpRightBoundary = c_Itr2->first;
                            break;
                        }
                    }
                    c_Itr2++;
                }

                minCells.insert(DblDbl_VT(tmpLeftBoundary, tmpRightBoundary));
                //cout << "tmpLeftBound: " << tmpLeftBoundary << ", tmpRightBound: " << tmpRightBoundary << endl;
            }
            c_Itr++;
            c_Itr1 = c_Itr;
            c_Itr1++;
        }
        firstTimeScan = false;
        //cout << "minOrder=" << minOrder << ", min_k=" << min_k << endl;
        //cout << "#minCells=" << minCells.size() << ", size of Singular="<< Singular.size() << endl;
        //getchar();

        //check the min_k with the current minimal order; if min_k is less than the currrent minimal order, stops;
        //otherwise, if the halflines of the minimal-ordered cell(s) are not singular, then expand the halflines
        if (minOrder > min_k) return min_k + NoOfDominators;

        tmpSkyline.clear();
        for (DblDbl_Itr = minCells.begin(); DblDbl_Itr != minCells.end(); DblDbl_Itr++) {
            //cout << "Left: " << DblDbl_Itr->first << ", Right: " << DblDbl_Itr->second << endl;
            //insert the bounding halflines of current minimal-ordered cell to expand
            c_Itr = Cells.find(DblDbl_Itr->first);
            c_Itr1 = Cells.find(DblDbl_Itr->second);
            if (c_Itr != Cells.end()) {
                long int NoOfSingulars = 0;
                long int NoOfHalfLines = 0;
                std::pair<set<long int>::iterator, bool> iRes;
                for (int i = 0; i < 2; i++) {
                    theMinCell[i] = c_Itr->first;

                    //cout << "processing min-cell with intersection=" << c_Itr->first << ", left sz=" << c_Itr->second->HalflineL.size() << ", right sz="<< c_Itr->second->HalflineR.size() << endl;
                    for (sItr = c_Itr->second->HalflineL.begin(); sItr != c_Itr->second->HalflineL.end(); sItr++) {
                        //cout << "processing the left halfline "<< *sItr << endl;
                        iRes = Singular.insert(*sItr);
                        if (iRes.second == true)   //current halfline with id *sItr is augmented, so remove it from SL
                            tmpSkyline.insert(*sItr);
                        else
                            NoOfSingulars++;
                    }
                    for (sItr = c_Itr->second->HalflineR.begin(); sItr != c_Itr->second->HalflineR.end(); sItr++) {
                        //cout << "processing right halfline "<< *sItr << endl;
                        iRes = Singular.insert(*sItr);
                        if (iRes.second == true)   //current halfline with id *sItr is augmented, so remove it from SL
                            tmpSkyline.insert(*sItr);
                        else
                            NoOfSingulars++;
                    }
                    NoOfHalfLines = NoOfHalfLines + c_Itr->second->HalflineL.size() + c_Itr->second->HalflineR.size();
                    c_Itr = c_Itr1;
                }
                if (NoOfSingulars == NoOfHalfLines) {
                    //current cell has only singular halflines, so we fix its order and compare it with min_k
                    if (minOrder > min_k) {
                        FoundMink = true;
                        break;
                    } else {
                        min_k = minOrder;
                        NoOfMinCells++;
                    }
                }
            }
        }
        //cout << "#skylines to remove: "<< tmpSkyline.size() << endl;
        //for (sItr=tmpSkyline.begin();sItr!=tmpSkyline.end();sItr++)
        //     cout << *sItr << " ";
        //cout << endl;
    }

    cout << "minimal k= " << min_k << ", NoOfDominators=" << NoOfDominators << ", NoOfDominees=" << NoOfDominees
         << endl;
    cout << "#Cell with order " << min_k << ": " << NoOfMinCells << endl;

    totalNoOfminCells = totalNoOfminCells + NoOfMinCells;

    return min_k + NoOfDominators;
}

long int S3::AA_HD(   //Advanced Algorithm for general dimensionality
        const int dimen, //dimensionality
        Rtree &a_rtree,
        float *PG[],
        Point &a_pt,
        long int &totalNoOfminCells,
        int &maxStackSize, Array &a_pageaccessed) {
    multimap<long int, float>::iterator IntDbl_Itr, IntDbl_Itr1;
    typedef multimap<long int, float>::value_type IntDbl_m_VT;

    map<float, set<float> *> ResultCells;
    map<float, set<float> *>::iterator DblS_mItr, DblS_mItr1;
    typedef map<float, set<float> *>::value_type m_VT;
    set<float>::iterator Dbl_sItr, Dbl_sItr1;

    std::set<long int> Singular;
    std::pair<set<long int>::iterator, bool> rset;

    std::vector<long int> ResultID;
    vector<long int>::iterator vecItr, vecItr1;

    std::vector<long int> skyline;

    std::set<long int>::iterator sItr, sItr1;

    std::multimap<long int, VirtualRNode *> NonResultEntry;   //data nodes that are kept for incremental skyline computation
    std::multimap<long int, VirtualRNode *>::iterator IntVRN_Iter;

    std::map<float, long int> Intersections;
    std::map<long int, float> RdAndCrossPt;
    std::map<float, float>::iterator mDblDblItr;
    typedef map<float, float>::value_type DblDbl_Pair;
    std::map<long int, float>::iterator mDblItr;
    typedef map<long int, float>::value_type IntDbl_Pair;

    map<float, long int>::iterator mItr, mItr1;
    typedef map<float, long int>::value_type DblInt_Pair;

    map<long int, long int>::iterator IntInt_mItr, IntInt_mItr1;
    typedef map<long int, long int>::value_type IntInt_mPair;

    long int min_k = INT_MAX;
    long int NoOfDominators = 0;
    long int NoOfDominees = 0;

    clock_t st, ed;
    double duration = 0;

    bool FoundMink = false;
    long int NoOfSkylinesToRemove = 0;
    long int NoOfMinCells = 0;

    HalfSpaces.clear();
    NonResultEntry.clear();
    RdIDtoHalfplaneID.clear();

    //compute the number of Dominators w.r.t query point q
    float cl[DMAX], cu[DMAX];
    for (int d = 0; d < dimen; d++) {
        cl[d] = 0;
        cu[d] = a_pt[d];
    }
    Hypercube Dominee_hc(dimen, cl, cu);
    for (int d = 0; d < dimen; d++) {
        cl[d] = a_pt[d];
        cu[d] = 1;
    }
    Hypercube Dominator_hc(dimen, cl, cu);

    //collect the records and index nodes that belong to the two incomparable windows for any d-dimensional space
    NoOfDominators = incomparableWindowsHD(a_rtree, PG, Dominee_hc, Dominator_hc, NonResultEntry, ResultID,
                                           maxStackSize, a_pageaccessed);
    //if (verbose)
    cout << "NoOfDominators=" << NoOfDominators << endl;

    //incrementally update skylines and then compute the cells with minimal order
    set<long int> tmpSkyline;

    //build a QuadTree
    QuadTree QT(Queryspace, dimen - 1, 0, MaxQuadTreeLevels, QuadNodeCapacity);
    QT.readCombinations();
    QT.splitNode(QT.root);

    while (!FoundMink) {
        //update the skyline by calling GetSkyline() with tmpSkyline containing the skyline points to delete
        NoOfNewlyAddedHalfSpaces = 0;
        if (verbose)
            cout << "Number of NonResultEntry: " << NonResultEntry.size() << ", Number of index nodes:"
                 << ResultID.size() << endl;
        GetSkylines(dimen, a_rtree, NonResultEntry, ResultID, tmpSkyline, PG, maxStackSize, a_pageaccessed);
        //if (verbose)
        cout << "Number of skylines: " << tmpSkyline.size() << endl;
        long int NoOfExistedSkylines = 0;
        st = clock();
        for (sItr = tmpSkyline.begin(); sItr != tmpSkyline.end(); sItr++) {
            long int id = (*sItr);

            float rd[DMAX];
            int tmpDim = 0;
            for (int d = 0; d < dimen; d++) {
                rd[d] = (PG[id][d] + PG[id][dimen + d]) / 2.0;
                if (fabs(rd[d] - a_pt[d]) <= 1e-6) tmpDim++;
            }
            if (tmpDim == dimen) {
                //cout << "CAUTION: there exist two exactly the same points!" << endl;  //i.e.,query point is a randomly picked record
                //getchar();
                continue;
            }

            if (RdIDtoHalfplaneID.size() > 0) {
                IntInt_mItr = RdIDtoHalfplaneID.find(id);
                if (IntInt_mItr != RdIDtoHalfplaneID.end())    //current skyline point has been added to the Quadtree
                {
                    NoOfExistedSkylines++;
                    continue;
                }
            }

            vector<float> tmpHS;
            float rd_d = rd[dimen - 1];
            float p_d = a_pt[dimen - 1];
            for (int d = 0; d < dimen - 1; d++)
                tmpHS.push_back((rd[d] - rd_d - a_pt[d] + p_d));

            tmpHS.push_back(p_d - rd_d);
            tmpHS.push_back(id);            //store the ID of incomparable record in HalfSpace
            HalfSpaces.push_back(tmpHS);          //form the half-space defined by the incomparable record and p
            RdIDtoHalfplaneID.insert(IntInt_mPair(id, HalfSpaces.size()));
            NoOfNewlyAddedHalfSpaces++;
            NoOfTotalHalfspacesAdded++;
        }
        if (verbose) cout << "No. of existed skylines: " << NoOfExistedSkylines << endl;
        if (verbose) cout << "No. of newly added skylines: " << NoOfNewlyAddedHalfSpaces << endl;
        if (verbose) cout << "No. of halfspaces: " << HalfSpaces.size() << endl;
        //getchar();

        QT.insertHalfSpaces(QT.root); //for each skyline, form a halfspace and insert it into the QuadTree
        //if (verbose)
        {
            cout << "QuadTree height: " << maximumLevel << endl;
            cout << "Number of tree nodes: " << numOfTreeNodes << endl;
            cout << "Number of leaf nodes: " << numOfLeafNodes << endl;
            cout << "Number of invalid leaf nodes: " << QT.NoOfInValidNodes << endl;
            cout << "Total number of node splits: " << numOfNodeSplits << endl;
        }
        //getchar();

        ed = clock();
        duration = float(ed - st) * 1000 / CLOCKS_PER_SEC;
        timeBuildQuadTree = timeBuildQuadTree + duration;

        /*
        //code to test how many leaf nodes cover a specified query point
        Point1 pt;
        pt.coord[0] = 0.14;
        pt.coord[1] = 0.27;

        set<long> ResultHalfSpaces;
        int numOfCoveredHalfSpaces = QT->countCoveredHalfSpaces(pt, QT->root,ResultHalfSpaces);
        cout << "There are " << numOfCoveredHalfSpaces << " halfspaces covered the query point:"<< endl;
        for (set<long>::iterator itr = ResultHalfSpaces.begin(); itr != ResultHalfSpaces.end(); itr++)
	     cout << (*itr) << "  ";
        getchar();
        //*/

        st = clock();
        entriesInLists.clear();  //measure memory usage for the quad tree
        vector<long int> numOfCoveredHS;
        Leaves.clear();
        QT.collectLeafNodes(QT.root, Leaves, numOfCoveredHS);
        //if (verbose)
        cout << "Number of leaf nodes to intersect: " << Leaves.size() << endl;
        //clock_t st=clock();
        std::sort(Leaves.begin(), Leaves.end());
        //clock_t ed=clock();
        //cout << "All leaves are sorted! time used in sorting: "<< float(ed-st)*1000/CLOCKS_PER_SEC << " ms." << endl;
        //getchar();

        map<long int, long int>::iterator tmpItr1;
        long int tmpCntr = 0;
        for (tmpItr1 = entriesInLists.begin(); tmpItr1 != entriesInLists.end(); tmpItr1++) {
            tmpCntr = tmpCntr + tmpItr1->second;
        }
        NoOfEntriesInLists = tmpCntr;
        numOfLeavesToIntersect = Leaves.size();


        /*
        //multimap<long int, QuadNode *>::iterator mQItr;
        vector<std::pair<long,QuadNode*> >::iterator mQItr;
        long int sZ=0;
        multimap<long int, long int> tmpMap;
        multimap<long int, long int>::iterator llItr;
        for (mQItr=Leaves.begin();mQItr!=Leaves.end();mQItr++)
        {
                 sZ=((*mQItr).second)->intersectedHalfspace.size();
                 //cout << "Node id: " << mQItr->second->NodeID << ", #inter.HS: " << sZ << endl;
                 tmpMap.insert(std::pair<long int, long int>(sZ,(*mQItr).second->NodeID));
        }
        for (llItr=tmpMap.begin();llItr!=tmpMap.end();llItr++)
             cout << "#inter.HS: " << llItr->first << ", node ID: " << llItr->second << endl;
        getchar();
        //*/

        long int m_k;
        vector<set<long int> > minCellHalfSpaces;
        vector<vector<char> > binaryString;
        if (optWithinNodeIntersection) {
            //use the optimized version of within-node intersection (i.e., the optimized Task 4)
            m_k = QT.optimizedInNodeIntersection(Leaves, minCellHalfSpaces, binaryString);
            if (verbose) cout << "Optimized task 4 testing is finished! " << endl;
        } else {
            //use the naive version of within-node intersection (i.e., the naive Task 4)
            m_k = QT.naiveInNodeIntersection(Leaves, minCellHalfSpaces, binaryString);
            if (verbose) cout << "Naive task 4 testing is finished! " << endl;

        }
        ed = clock();
        duration = float(ed - st) * 1000 / CLOCKS_PER_SEC;
        timeNodeIntersection = timeNodeIntersection + duration;

        if (verbose)
            cout << "m_k=" << m_k << ", min_k=" << min_k << ", #tmpMinCells: " << minCellHalfSpaces.size()
                 << ", #dominators= " << NoOfDominators << endl;
        //getchar();

        //if (m_k>min_k) return min_k+NoOfDominators;

        //augment the bounding halfspaces of current cell if they are not augmented yet
        std::pair<set<long int>::iterator, bool> iRes;
        long int NoOfSingulars;
        long int tmpCellCnt = minCellHalfSpaces.size();
        tmpSkyline.clear();
        NoOfMinCells = 0;
        if (verbose)
            cout << "#singulars=" << Singular.size() << ", #temp minCells=" << tmpCellCnt << endl;
        for (int i = 0; i < tmpCellCnt; i++)  //examine each interim minCells
        {
            NoOfSingulars = 0;
            for (sItr = minCellHalfSpaces[i].begin(); sItr != minCellHalfSpaces[i].end(); sItr++) {
                //check the bounding halfspaces of current mincell to see whether they are augmented yet
                long int Rid = HalfSpaces[(*sItr)][dimen];
                //cout << "Rid=" << Rid << endl;
                iRes = Singular.insert(Rid);
                if (iRes.second ==
                    true)   //current halfline with id *sItr is augmented, so remove it from the skyline set
                    tmpSkyline.insert(Rid);
                else
                    NoOfSingulars++;
            }
            if (NoOfSingulars == minCellHalfSpaces[i].size()) {
                if (m_k <= min_k)  //found the min_k and corresponding cell with minimal order
                {
                    NoOfMinCells++;
                    min_k = m_k;
                } else {          //skip current cell whose order is greater than min_k
                    continue;
                }
            }
            //cout << "#singulars: " << NoOfSingulars << endl;
            //if (NoOfMinCells > 0)
            //    cout << "Found " << NoOfMinCells << " mincells with m_k=" << m_k << ", min_k="<< min_k << endl;
        }

        if (verbose)
            cout << "#skylines to remove: " << tmpSkyline.size() << endl;
        if (tmpSkyline.size() == 0)
            FoundMink = true;
        else
            NoOfSkylinesToRemove = tmpSkyline.size();

        /*
        for (sItr=tmpSkyline.begin();sItr!=tmpSkyline.end();sItr++)
             cout << *sItr << " ";
        cout << endl;
        //*/
        //getchar();
    }

    if (verbose)
        cout << "#minCells with order " << min_k << ": " << NoOfMinCells << ", #Dominators=" << NoOfDominators << endl;
    //getchar();

    totalNoOfminCells = totalNoOfminCells + NoOfMinCells;
    cout << "Return to main..." << endl;

    //free memory
    long int tmpFreeCount = 0;
    for (IntVRN_Iter = NonResultEntry.begin(); IntVRN_Iter != NonResultEntry.end(); IntVRN_Iter++) {
        tmpFreeCount++;
        delete IntVRN_Iter->second;
    }
    cout << "Freed " << tmpFreeCount << " nodes!" << endl;

    if (min_k == INT_MAX) return NoOfDominators;
    return min_k + NoOfDominators;
}

int S3::BA_2D(    //BA algorithm for 2-d case, i.e., FCA algorithm
        const int dimen, //dimensionality
        Rtree &a_rtree,
        float *PG[],
        Point &a_pt,
        long int &totalNoOfminCells,
        int &maxStackSize, Array &a_pageaccessed) {
    vector<long int> ResultID[4];
    vector<long int>::iterator vecItr;

    map<float, set<long int> > Intersections;
    map<float, set<long int> >::iterator mItr, mItr1;
    typedef map<float, set<long int> >::value_type VT;
    set<long int>::iterator sItr;

    map<float, long int> CellsWithOrder;
    map<float, long int>::iterator mDblItr;
    typedef map<float, long int>::value_type VT1;

    float dataspace[4][2] = {{0, 0},
                             {1, 0},
                             {1, 1},
                             {0, 1}};
    float TMPcenter[DMAX];
    float a_l[DMAX];    //storing the ranges (w.r.t point pt) along each dimension of a quadrant

    long int NoOfDominators = 0;
    long int NoOfDominees = 0;

    long int min_k = INT_MAX;
    long int initOrder = 0;

    for (int i = 0; i < 4; i++)  //for each of the four quadrants (defined by a_pt), perform a window query
    {
        TMPcenter[0] = float(dataspace[i][0] + a_pt[0]) / 2.0;
        TMPcenter[1] = float(dataspace[i][1] + a_pt[1]) / 2.0;
        a_l[0] = fabs(TMPcenter[0] - a_pt[0]);
        a_l[1] = fabs(TMPcenter[1] - a_pt[1]);

        Point center(dimen, TMPcenter);
        window(a_rtree, PG, center, a_l, ResultID[i], maxStackSize, a_pageaccessed);

        cout << "i= " << i << ", #results of window query: " << ResultID[i].size() << endl;
        //cout << "Center[0]=" << TMPcenter[0] << ", Center[1]=" << TMPcenter[1] << endl;
        //cout << "length[0]=" << a_l[0] << ", length[1]=" << a_l[1] << endl;

        if (i == 0) //skip the dominees of a_pt
        {
            NoOfDominees = ResultID[0].size();
            continue;
        }

        if (i == 2)  //skip the dominators of a_pt
        {
            NoOfDominators = ResultID[2].size();
            continue;
        }

        if (i != 2)   //process quadrants that contain incomparable records
        {
            for (vecItr = ResultID[i].begin(); vecItr != ResultID[i].end(); vecItr++)  //find the intersection points
            {
                long int id = (*vecItr);

                float rd[2] = {(PG[id][0] + PG[id][2]) / 2.0, (PG[id][1] + PG[id][3]) / 2.0};

                int tmpDim = 0;
                for (int d = 0; d < dimen; d++)
                    if (fabs(rd[d] - a_pt[d]) <= 1e-6) tmpDim++;
                if (tmpDim == dimen) {
                    //query point is a randomly picked record from dataset
                    //cout << "CAUTION: there exist two exactly the same points!" << endl;
                    //getchar();
                    continue;
                }

                if (rd[1] > a_pt[1]) initOrder++;  //obtain the initial order of current query point

                float denominator = rd[0] - rd[1] - a_pt[0] + a_pt[1];
                float crosspoint;
                if (fabs(denominator) >= 1e-6)
                    crosspoint = float(a_pt[1] - rd[1]) / denominator;
                else
                    continue;
                if (crosspoint < 0 || crosspoint > 1) continue;

                //float score = rd[0]*crosspoint+rd[1]*(1-crosspoint);

                mItr = Intersections.find(crosspoint);
                if (mItr == Intersections.end()) {
                    set<long int> tmpSet;
                    tmpSet.insert(id);
                    Intersections.insert(VT(crosspoint, tmpSet));
                } else {
                    mItr->second.insert(id);
                }
            }
        }
    }
    set<long int> tmpSet;
    Intersections.insert(VT(1, tmpSet));
    initOrder++;
    //cout << "Initial order: " << initOrder << ", #intersections: "<< Intersections.size() << endl;
    //getchar();

    CellsWithOrder.insert(VT1(0, initOrder));  //the first cell from the left most intersection point 0
    min_k = initOrder;
    for (mItr = Intersections.begin(); mItr != Intersections.end(); mItr++) {
        for (sItr = mItr->second.begin(); sItr != mItr->second.end(); sItr++) {
            long int id = (*sItr);
            float rd[2] = {(PG[id][0] + PG[id][2]) / 2.0, (PG[id][1] + PG[id][3]) / 2.0};

            if (rd[0] > a_pt[0])
                initOrder++;
            else if (rd[0] <= a_pt[0])
                initOrder--;
        }
        CellsWithOrder.insert(VT1(mItr->first, initOrder));
        if (min_k > initOrder) min_k = initOrder;
        //cout << "intersection: " << mItr->first <<", order: "<< initOrder <<", min_k: "<< min_k << endl;
        //getchar();
    }

    long int NoOfMinCells = 0;
    for (mDblItr = CellsWithOrder.begin(); mDblItr != CellsWithOrder.end(); mDblItr++)
        if (min_k == mDblItr->second) NoOfMinCells++;

    cout << "minimal k= " << min_k << ", #minCells= " << NoOfMinCells << ", NoOfDominators=" << NoOfDominators << endl;

    totalNoOfminCells = totalNoOfminCells + NoOfMinCells;

    return min_k + NoOfDominators;
}

int S3::BA_2D_old(    //my version of BA algorithm for 2-d case
        const int dimen, //dimensionality
        Rtree &a_rtree,
        float *PG[],
        Point &a_pt,
        long int &totalNoOfminCells,
        int &maxStackSize, Array &a_pageaccessed) {
    std::vector<long int> ResultID[4];
    vector<long int>::iterator vecItr;

    std::map<float, long int> Intersections;
    std::map<long int, float> RdAndCrossPt;
    std::map<float, float>::iterator mDblDblItr;
    typedef map<float, float>::value_type DblDbl_Pair;
    std::map<long int, float>::iterator mDblItr;
    typedef map<long int, float>::value_type IntDbl_Pair;

    multimap<float, long int> LeftFixedIntervals, RightFixedIntervals;
    multimap<float, long int>::iterator mmDblIntItr, mmDblIntItr1;
    typedef multimap<float, long int>::value_type mDblInt_Pair;

    map<float, long int>::iterator mItr, mItr1;
    typedef map<float, long int>::value_type DblInt_Pair;
    map<float, set<long int> *>::iterator mVecItr, mVecItr1;
    typedef map<float, set<long int> *>::value_type DblVec_Pair;

    float dataspace[4][2] = {{0, 0},
                             {1, 0},
                             {1, 1},
                             {0, 1}};
    float TMPcenter[DMAX];
    float a_l[DMAX];    //storing the ranges (w.r.t point pt) along each dimension of a quadrant

    long int NoOfDominators = 0;
    long int NoOfDominees = 0;

    for (int i = 0; i < 4; i++)  //for each of the four quadrants (defined by a_pt), perform a window query
    {
        TMPcenter[0] = float(dataspace[i][0] + a_pt[0]) / 2.0;
        TMPcenter[1] = float(dataspace[i][1] + a_pt[1]) / 2.0;
        a_l[0] = fabs(TMPcenter[0] - a_pt[0]);
        a_l[1] = fabs(TMPcenter[1] - a_pt[1]);

        Point center(dimen, TMPcenter);
        window(a_rtree, PG, center, a_l, ResultID[i], maxStackSize, a_pageaccessed);

        cout << "i= " << i << ", #results of window query: " << ResultID[i].size() << endl;
        //cout << "Center[0]=" << TMPcenter[0] << ", Center[1]=" << TMPcenter[1] << endl;
        //cout << "length[0]=" << a_l[0] << ", length[1]=" << a_l[1] << endl;

        if (i == 0) //skip the dominees of a_pt
        {
            NoOfDominees = ResultID[0].size();
            continue;
        }

        if (i == 2)  //skip the dominators of a_pt
        {
            NoOfDominators = ResultID[2].size();
            continue;
        }

        if (i != 2)   //process quadrants that contain incomparable records
        {
            for (vecItr = ResultID[i].begin(); vecItr != ResultID[i].end(); vecItr++)  //find the intersection points
            {
                long int id = (*vecItr);

                float rd[2] = {(PG[id][0] + PG[id][2]) / 2.0, (PG[id][1] + PG[id][3]) / 2.0};

                int tmpDim = 0;
                for (int d = 0; d < dimen; d++)
                    if (fabs(rd[d] - a_pt[d]) <= 1e-6) tmpDim++;
                if (tmpDim == dimen) {
                    //query point is a randomly picked record from dataset
                    //cout << "CAUTION: there exist two exactly the same points!" << endl;
                    //getchar();
                    continue;
                }

                float denominator = rd[0] - rd[1] - a_pt[0] + a_pt[1];
                float crosspoint;
                //if (fabs(denominator)>=1e-6)
                //if (denominator>1e-6)
                crosspoint = float(a_pt[1] - rd[1]) / denominator;
                //else
                //    continue;
                if (crosspoint < 0 || crosspoint > 1) continue;

                float score = rd[0] * crosspoint + rd[1] * (1 - crosspoint);
                RdAndCrossPt.insert(IntDbl_Pair(id, crosspoint));

                Intersections.insert(DblInt_Pair(crosspoint, 0));

                /*/
                 mVecItr = Intervals.find(crosspoint);
                 if (mVecItr != Intervals.end())
                    (mVecItr->second)->insert(id);  //keep the points that are responsable for an intersection point with respect to a_pt
                 else
                 {
                     set<long int>* vec = new set<long int>;
                     vec->insert(id);
                     Intervals.insert(DblVec_Pair(crosspoint,vec));
                 }
                 //*/
            }
        }
    }
    Intersections.insert(DblInt_Pair(1, 0));

    /*/
    mVecItr = Intervals.find(1);
    if (mVecItr == Intervals.end())
    {
        set<long int>* vec = new set<long int>;
        Intervals.insert(DblVec_Pair(1,vec));
    }
    //*/

    //compute the number of records whose score is greater/smaller than S(a_pt) within each intervals
    cout << "begin examining " << Intersections.size() << " intersection points..." << endl;
    long int NumOfProcessed = 0;
    for (int i = 0; i < 4; i++) {
        if (i == 0 || i == 2) continue;
        cout << "i=" << i << endl;

        for (vecItr = ResultID[i].begin(); vecItr != ResultID[i].end(); vecItr++) {
            long int id = *vecItr;
            float rd[2] = {(PG[id][0] + PG[id][2]) / 2.0, (PG[id][1] + PG[id][3]) / 2.0};

            int tmpDim = 0;
            for (int d = 0; d < dimen; d++)
                if (fabs(rd[d] - a_pt[d]) <= 1e-6) tmpDim++;
            if (tmpDim == dimen) {
                //query point is a randomly picked record from dataset
                //cout << "CAUTION: there exist two exactly the same points!" << endl;
                //getchar();
                continue;
            }

            NumOfProcessed++;
            if ((NumOfProcessed) % 2000 == 0)
                cout << NumOfProcessed << " result records processed..." << endl;

            mDblItr = RdAndCrossPt.find(id);
            if (mDblItr != RdAndCrossPt.end()) {
                float tmpQ1 = mDblItr->second;
                float score = rd[0] * tmpQ1 + rd[1] * (1 - tmpQ1);

                if (rd[1] < a_pt[1])   //current record will score higher than a_pt in interval (crosspoint,1]
                {
                    //optimized version begin
                    RightFixedIntervals.insert(mDblInt_Pair(tmpQ1, id));
                    //optimized version end

                    /*
                 mItr = Intersections.find(tmpQ1);
                 if (mItr != Intersections.end())
                 {
                     mItr1 = mItr;
                     mItr1++;
                     for (; mItr1 != Intersections.end(); mItr1++)
                          mItr1->second++;

                     //
                     //mVecItr = Intervals.find(tmpQ1);
                     //if (mVecItr != Intervals.end())
                     //{
                     //    mVecItr1 = mVecItr;
                     //    mVecItr1++;
                     //    for (; mVecItr1!=Intervals.end();mVecItr1++)
                     //        (mVecItr1->second)->insert(id);  //keep the points that are responsable for
                     //}
                     ///
                 }
                 //*/
                } else if (rd[1] >
                           a_pt[1])//(rd[0] < a_pt[0])  //current record will score higher than a_pt in interval [0,crosspoint)
                {
                    LeftFixedIntervals.insert(mDblInt_Pair(tmpQ1, id));

                    /*/
                  mItr = Intersections.find(tmpQ1);
                  if (mItr != Intersections.end())
                  {
                      mItr++;
                      for (mItr1 = Intersections.begin(); mItr1!=mItr; mItr1++)
                           mItr1->second++;

                     //
                     //mVecItr = Intervals.find(tmpQ1);
                     //if (mVecItr != Intervals.end())
                     //{
                     //    mVecItr++;
                     //    for (mVecItr1=Intervals.begin();mVecItr1!=mVecItr;mVecItr1++)
                     //        (mVecItr1->second)->insert(id);  //keep the points that are responsable for
                     //}
                     ///
                  }
                  //*/
                }
            }
        }
    }
    cout << "finish examining intersection points..." << endl;
    cout << "Size of Left=" << LeftFixedIntervals.size() << ", Size of Right=" << RightFixedIntervals.size() << endl;

    /*/
            cout << "Intersections (#rds above p): " << endl;
            for (mItr1 = Intersections.begin(); mItr1!=Intersections.end(); mItr1++)
            {
                 cout << "("<< mItr1->first << ", " << mItr1->second << ") ";
            }
            cout << "Intersections (#rds below p): " << endl;
            for (mItr1 = InterBelow.begin(); mItr1!=InterBelow.end(); mItr1++)
            {
                 cout << "("<< mItr1->first << ", " << mItr1->second << ") ";
            }
            //*/

    /*/
          cout << "Number of intervals: " << Intervals.size() << endl;
         for (mVecItr=Intervals.begin();mVecItr!=Intervals.end();mVecItr++)
         {
              cout << "crosspoint=" << mVecItr->first << "(" << (mVecItr->second)->size() << "): ";
              //for (set<long int>::iterator vItr=(mVecItr->second)->begin();vItr!=(mVecItr->second)->end();vItr++)
              //     cout << *vItr << " ";
              //cout << endl;
         }
         //*/

    //optimized version to compute the min_k begin
    long int min_k = 99999999, max_k = -9999999;
    long int SizeLeft = LeftFixedIntervals.size();
    long int SizeRight = RightFixedIntervals.size();
    long int index = 1;
    for (mmDblIntItr = LeftFixedIntervals.begin(); mmDblIntItr != LeftFixedIntervals.end(); mmDblIntItr++)
        mmDblIntItr->second = index++;

    index = 1;
    for (mmDblIntItr = RightFixedIntervals.begin(); mmDblIntItr != RightFixedIntervals.end(); mmDblIntItr++)
        mmDblIntItr->second = index++;

    for (mmDblIntItr = LeftFixedIntervals.begin(); mmDblIntItr != LeftFixedIntervals.end(); mmDblIntItr++) {
        long int tmpRank = mmDblIntItr->second;
        float Rbound = mmDblIntItr->first;

        if ((tmpRank % 2000) == 0)
            cout << "Examining left, " << tmpRank << " records processed..." << endl;

        mmDblIntItr1 = RightFixedIntervals.lower_bound(Rbound);
        if (mmDblIntItr1 != RightFixedIntervals.end()) {
            tmpRank = tmpRank + (SizeRight - mmDblIntItr1->second);
        }
        if (min_k > tmpRank) min_k = tmpRank;
        if (max_k < tmpRank) max_k = tmpRank;
    }

    for (mmDblIntItr = RightFixedIntervals.begin(); mmDblIntItr != RightFixedIntervals.end(); mmDblIntItr++) {
        float Lbound = mmDblIntItr->first;

        long int tmpRank = mmDblIntItr->second; //(SizeRight - mmDblIntItr->second);

        if ((tmpRank % 2000) == 0)
            cout << "Examining right, distance=" << tmpRank << endl;

        mmDblIntItr1 = LeftFixedIntervals.lower_bound(Lbound);
        if (mmDblIntItr1 != LeftFixedIntervals.end()) {
            tmpRank = tmpRank + (SizeLeft - mmDblIntItr1->second);
        }

        if (min_k > tmpRank) min_k = tmpRank;
        if (max_k < tmpRank) max_k = tmpRank;
    }
    if (RightFixedIntervals.size() == 0 && LeftFixedIntervals.size() > 0) {
        min_k = 0;
        max_k = LeftFixedIntervals.size();
    }
    if (RightFixedIntervals.size() > 0 && LeftFixedIntervals.size() == 0) {
        min_k = 0;
        max_k = RightFixedIntervals.size();
    }
    //optimized version to compute the min_k end

    /*/count in the number of dominators
    long int min_k = 99999999, max_k=-9999999;
    for (mItr1 = Intersections.begin(); mItr1 != Intersections.end(); mItr1++)
    {
         if (mItr1->second < min_k) min_k = mItr1->second;
         if (mItr1->second > max_k) max_k = mItr1->second;
    }
    //*/
    cout << "minimal k= " << min_k + NoOfDominators << ", maximum k=" << max_k + NoOfDominators << " NoOfDominators="
         << NoOfDominators << ", NoOfDominees=" << NoOfDominees << endl;

    //totalNoOfminCells=totalNoOfminCells+;

    return min_k + NoOfDominators;
}

int S3::BA_HD(    //BA algorithm for general dimensionality
        const int dimen, //dimensionality
        Rtree &a_rtree,
        float *PG[],
        Point &a_pt,
        long int &totalNoOfminCells,
        int &maxStackSize, Array &a_pageaccessed) {
    multimap<long int, float>::iterator IntDbl_Itr, IntDbl_Itr1;
    typedef multimap<long int, float>::value_type IntDbl_m_VT;

    vector<set<long int> > minCellHalfSpaces;
    vector<vector<char> > binaryString;

    map<float, set<float> *> ResultCells;
    map<float, set<float> *>::iterator DblS_mItr, DblS_mItr1;
    typedef map<float, set<float> *>::value_type m_VT;
    set<float>::iterator Dbl_sItr, Dbl_sItr1;

    std::set<long int> Singular;
    std::pair<set<long int>::iterator, bool> rset;

    std::vector<long int> ResultID;
    vector<long int>::iterator vecItr, vecItr1;

    std::set<long int>::iterator sItr, sItr1;

    std::map<float, float>::iterator mDblDblItr;
    typedef map<float, float>::value_type DblDbl_Pair;
    std::map<long int, float>::iterator mDblItr;
    typedef map<long int, float>::value_type IntDbl_Pair;

    map<float, long int>::iterator mItr, mItr1;
    typedef map<float, long int>::value_type DblInt_Pair;

    map<long int, long int>::iterator IntInt_mItr, IntInt_mItr1;
    typedef map<long int, long int>::value_type IntInt_mPair;

    long int min_k = INT_MAX;
    long int NoOfDominators = 0;
    long int NoOfDominees = 0;

    bool FoundMink = false;
    long int NoOfMinCells = 0;

    clock_t st, ed;
    double duration;

    HalfSpaces.clear();
    RdIDtoHalfplaneID.clear();

    //build the dominator and dominee window to eliminate records falling inside the two windows
    float cl[DMAX], cu[DMAX];
    for (int d = 0; d < dimen; d++) {
        cl[d] = 0;
        cu[d] = a_pt[d];
    }
    Hypercube Dominee_hc(dimen, cl, cu);
    for (int d = 0; d < dimen; d++) {
        cl[d] = a_pt[d];
        cu[d] = 1;
    }
    Hypercube Dominator_hc(dimen, cl, cu);

    //collect all the incomparable records
    NoOfDominators = getIncomparableRecords(a_rtree, PG, Dominee_hc, Dominator_hc, ResultID, maxStackSize,
                                            a_pageaccessed);
    cout << "NoOfDominators: " << NoOfDominators << ", #Incomparable Records: " << ResultID.size() << endl;

//return -1; //this line is for measuring IO only

    //build a QuadTree
    QuadTree QT(Queryspace, dimen - 1, 0, MaxQuadTreeLevels, QuadNodeCapacity);
    QT.readCombinations();
    QT.splitNode(QT.root);

    //insert into QuadTree all halfspaces where each of which corresponds to an incomparable rd
    cout << "Begin to insert halfspaces into the QuadTree..." << endl;
    long int counter = 0;
    st = clock();
    for (vecItr = ResultID.begin(); vecItr != ResultID.end(); vecItr++) {
        if ((counter % 5000) == 0)
            cout << counter << " halfspaces have been inserted..." << endl;

        long int id = (*vecItr);

        float rd[DMAX];
        int tmpDim = 0;
        for (int d = 0; d < dimen; d++) {
            rd[d] = (PG[id][d] + PG[id][dimen + d]) / 2.0;
            if (fabs(rd[d] - a_pt[d]) <= 1e-6) tmpDim++;
        }
        if (tmpDim == dimen)  //skip the focal point
        {
            //cout << "CAUTION: there are two records exactly the same!" << endl;
            continue;
        }

        vector<float> tmpHS;
        float rd_d = rd[dimen - 1];
        float p_d = a_pt[dimen - 1];
        for (int d = 0; d < dimen - 1; d++)
            tmpHS.push_back((rd[d] - rd_d - a_pt[d] + p_d));

        tmpHS.push_back(p_d - rd_d);
        tmpHS.push_back(id);         //store the ID of incomparable record in HalfSpace
        HalfSpaces.push_back(tmpHS); //form a half-space defined by a incomparable rd and p
        counter++;
        NoOfTotalHalfspacesAdded++;
    }

    QT.insertHalfSpaces(QT.root);
    cout << "QuadTree height: " << maximumLevel << endl;
    cout << "Number of tree nodes: " << numOfTreeNodes << endl;
    cout << "Number of leaf nodes: " << numOfLeafNodes << endl;
    cout << "Total number of node splits: " << numOfNodeSplits << endl;

    ed = clock();
    duration = float(ed - st) * 1000 / CLOCKS_PER_SEC;
    timeBuildQuadTree = timeBuildQuadTree + duration;

    st = clock();
    vector<long int> numOfCoveredHS;
    Leaves.clear();
    QT.collectLeafNodes(QT.root, Leaves, numOfCoveredHS);
    cout << "Number of leaf nodes to intersect: " << Leaves.size() << endl;
    //clock_t st=clock();
    std::sort(Leaves.begin(), Leaves.end());
    //clock_t ed=clock();
    //cout <<"All " << Leaves.size() << " leaves are sorted! Time used in sorting: " << float(ed-st)*1000/CLOCKS_PER_SEC <<" ms." << endl;

    numOfLeavesToIntersect = Leaves.size();
    // count the entries in lists of each quadtree node
    map<long int, long int>::iterator tmpItr1;
    long int tmpCntr = 0;
    for (tmpItr1 = entriesInLists.begin(); tmpItr1 != entriesInLists.end(); tmpItr1++) {
        tmpCntr = tmpCntr + tmpItr1->second;
    }
    NoOfEntriesInLists = tmpCntr;
    //

//return -1; //this line is for measuring IO only

    //output the memory usage info
    /*
    fstream finMem;
    finMem.open("/proc/meminfo",ios::in);
    while (true)
    {
       string line;
       finMem >> line;
       if (finMem.eof())break;
       cout << line << endl;
    }
    finMem.close();
    finMem.open("/proc/self/status",ios::in);
    while (true)
    {
       string line;
       finMem >> line;
       if (finMem.eof())break;
       cout << line << endl;
    }
    finMem.close();
    //*/

    //QT.traversal(QT->root);

    /*
        vector<std::pair<long,QuadNode*> >::iterator mQItr;
        long int sZ=0;
        multimap<long int, long int> tmpMap;
        multimap<long int, long int>::iterator llItr;
        long int tmpCnter=0;
        for (mQItr=Leaves.begin();mQItr!=Leaves.end();mQItr++)
        {
             sZ=((*mQItr).second)->intersectedHalfspace.size();
             cout << "Node id: " << (*mQItr).second->NodeID << ", #inter.HS: " << sZ << ", #coveredHS: "<< (*mQItr).first << endl;
             tmpMap.insert(std::pair<long int, long int>(sZ,(*mQItr).second->NodeID));
             tmpCnter++;
             //if (tmpCnter>=10) break;
        }
        for (llItr=tmpMap.begin();llItr!=tmpMap.end();llItr++)
             cout << "#inter.HS: " << llItr->first << ", node ID: " << llItr->second << endl;
        //getchar();
        //*/



    //Within-node halfspace intersection
    if (optWithinNodeIntersection) {
        min_k = QT.optimizedInNodeIntersection(Leaves, minCellHalfSpaces, binaryString);
        cout << "Optimzed within-node intersection finished!" << endl;
    } else {
        min_k = QT.naiveInNodeIntersection(Leaves, minCellHalfSpaces, binaryString);
        cout << "Naive halfspace intersection finished!" << endl;
    }

    ed = clock();
    duration = float(ed - st) * 1000 / CLOCKS_PER_SEC;
    timeNodeIntersection = timeNodeIntersection + duration;

    totalNoOfminCells = totalNoOfminCells + minCellHalfSpaces.size();
    Leaves.clear();

    cout << "Return to main..." << endl;

    if (min_k == INT_MAX) return NoOfDominators;
    return min_k + NoOfDominators;
}

int S3::GetSkylines(
        const int dimen,
        Rtree &a_rtree,             // returns a result and a result size
        std::multimap<long int, VirtualRNode *> &NonResultEntry,  //non-result set, storing rtree node which is simply a data entry
        std::vector<long int> &PrunedNodes,  //the heap for storing non-result r-tree nodes and remaining entries in H
        std::set<long> &a_skylines,
        float *PG[],
        int &maxStackSize, Array &a_pageaccessed) {
    multimap<float, long int> H0;  //the min-heap H0
    multimap<float, long int>::iterator DblInt_Iter;
    typedef multimap<float, long int>::value_type DblInt_Pair;
    typedef multimap<float, float>::value_type DblDbl_Pair;

    multimap<long int, VirtualRNode *> DataEntryQ;  //queue storing rtree node which is simply a data entry
    multimap<long int, VirtualRNode *>::iterator IntVRN_Iter;
    typedef multimap<long int, VirtualRNode *>::value_type IntVRN_Pair;

    vector<long int> skylines;
    vector<long int>::iterator Int_vIter;
    set<long int>::iterator Int_sIter;

    float mindist;

    float m_stack = -1;
    float pt[DMAX];
    bool isAnObject;
    float ORIGIN[DMAX];
    int m_SizeOfDataEntQ;

    for (int i = 0; i < dimen; i++) ORIGIN[i] = CEILING;

    //put the non-result points and pruned index entries to a min-heap sorted on mindist
    RtreeNodeEntry *e0;
    for (Int_vIter = PrunedNodes.begin(); Int_vIter != PrunedNodes.end(); Int_vIter++) {
        RtreeNode *Rnode = a_rtree.m_memory.loadPage(*Int_vIter);
        e0 = Rnode->genNodeEntry();     //compute the enclosing MBR for current index node
        for (int i = 0; i < dimen; i++) pt[i] = e0->m_hc.getUpper()[i];
        mindist = minDist(pt, ORIGIN, dimen);
        H0.insert(DblInt_Pair(mindist, *Int_vIter));
        delete Rnode;
        delete e0;
    }
    PrunedNodes.clear();
    if (verbose)
        cout << "size of H0: " << H0.size() << endl;

    vector<long int> nodeIdToDelete;
    for (IntVRN_Iter = NonResultEntry.begin(); IntVRN_Iter != NonResultEntry.end(); IntVRN_Iter++) {
        //each time GetSkylines is invoked, parameter a_skylines contains the skylines to remove,
        //so we remove those skylines and update the skyline set
        if (a_skylines.size() > 0) {
            Int_sIter = a_skylines.find(IntVRN_Iter->first - MAXPAGEID);
            if (Int_sIter != a_skylines.end()) {
                nodeIdToDelete.push_back(IntVRN_Iter->first);
                continue;
            }
        }

        for (int i = 0; i < dimen; i++) pt[i] = (IntVRN_Iter->second)->m_entry[0]->m_hc.getUpper()[i];
        mindist = minDist(pt, ORIGIN, dimen);
        H0.insert(DblInt_Pair(mindist, IntVRN_Iter->first));
    }
    //each time GetSkylines is invoked, parameter a_skylines contains the skylines to remove,
    //so we remove those skylines and update the skyline set
    for (Int_vIter = nodeIdToDelete.begin(); Int_vIter != nodeIdToDelete.end(); Int_vIter++) {
        IntVRN_Iter = NonResultEntry.find(*Int_vIter);
        delete IntVRN_Iter->second;
        NonResultEntry.erase(IntVRN_Iter);
    }
    if (verbose)
        cout << "size of H0: " << H0.size() << endl;

    //begin exploring the min-heap H0 in order to find skylines
    while (H0.size() != 0) {
        isAnObject = false;

        if (m_stack <= H0.size()) {
            m_stack = H0.size();
            m_SizeOfDataEntQ = NonResultEntry.size();
        }

        //maxStackSize = s.size() > maxStackSize ? s.size() : maxStackSize;
        DblInt_Iter = H0.begin();
        long int pageid = DblInt_Iter->second;
        float DistTmp = DblInt_Iter->first;

        H0.erase(DblInt_Iter);

        VirtualRNode *VirNode = new VirtualRNode;
        RtreeNodeEntry *e0;            // create node entry e0 for node n, so that its MBR can be obtained
        if (pageid >=
            MAXPAGEID)     //current element in H is a data entry node (for distinction,its pageid is equal to m_id+MAXPAGED)
        {
            isAnObject = true;

            IntVRN_Iter = NonResultEntry.find(pageid);
            if (IntVRN_Iter == NonResultEntry.end()) {
                cout << "Error! there is no node " << pageid << " in NonResultEntry!" << endl;
            } else {
                VirNode->copyData(IntVRN_Iter->second);

                if (VirNode->m_usedSpace > 1) {
                    cout << "There is a fatal error in retrieving data-entry node from queue!" << endl;
                    return -1;
                } else
                    e0 = VirNode->m_entry[0]->clone();  //retrieve the sole entry from the data-entry node

                //comment the two lines out, because we need to maintain all data points for skyline update later
                //delete IntVRN_Iter->second;         //free the virtualRnode; Note, this operation is VERY important
                //NonResultEntry.erase(IntVRN_Iter);
            }
            pageid = pageid - MAXPAGEID;
        } else {
            RtreeNode *node = a_rtree.m_memory.loadPage(pageid);

            //my code
            VirNode->copyData(*node);
            VirNode->copyEntries(*node, node->m_usedSpace);
            e0 = node->genNodeEntry();     //compute the enclosing MBR for current index node

            delete node;
        }
        //cout << "exam page " << pageid << " isObject "<< isAnObject << ", Heap size: " << H0.size() << " mindist:" << DistTmp << " #skylines:" << a_skylines.size() << endl;

        //dominance check for current node
        if (isAnObject)
            for (int i = 0; i < dimen; i++) pt[i] = VirNode->m_entry[0]->m_hc.getLower()[i] + SIDELEN;
        else
            for (int i = 0; i < dimen; i++) pt[i] = e0->m_hc.getUpper()[i];
        bool dominated = IsDominatedBy(dimen, pt, skylines, PG);

        if (!dominated) {
            //cout << "object " << pageid << " is not dominated!" << endl;

            //check whether current entry is leaf node or index node
            if (VirNode->isLeaf()) {
                if (VirNode->m_usedSpace >
                    1) //current node is a leaf node, so all its data entries must be put into priority queue H
                {
                    a_pageaccessed.append((void *) pageid);
                    for (int i = 0; i < VirNode->m_usedSpace; i++) {
                        for (int j = 0; j < dimen; j++)
                            pt[j] = VirNode->m_entry[i]->m_hc.getUpper()[j];

                        bool dominated = IsDominatedBy(dimen, pt, skylines, PG);

                        long int NegPageid = VirNode->m_entry[i]->m_id + MAXPAGEID;
                        VirtualRNode *node = new VirtualRNode; //for each data entries of n, insert it into data-entry queue

                        RtreeNodeEntry *Nentry = VirNode->m_entry[i]->clone();
                        node->insertEntry(Nentry);
                        NonResultEntry.insert(IntVRN_Pair(NegPageid, node));
                        delete Nentry;

                        //if current data entry is dominated, then we only keep it in the nonresultentry for later skyline update; otherwise, we put it into the search heap as well
                        if (dominated) continue;

                        //compute the mindist of current entry
                        mindist = minDist(pt, ORIGIN, dimen);
                        H0.insert(DblInt_Pair(mindist, NegPageid));
                    }
                } else    //if current node is a data node, then put it to the skyline set
                {
                    skylines.push_back(pageid);

                    /*  comment out the code below, because the data entry has already been kept in NonResultEntry
               	         VirtualRNode* node=new VirtualRNode; //for each entries of n, insert it into data-entry queue
                	 RtreeNodeEntry* Nentry=VirNode->m_entry[i]->clone();
                	 node->insertEntry(Nentry);
                	 NonResultEntry.insert(IntVRN_Pair(NegPageid,node));
                	 delete Nentry;
                         //*/

                    //cout << "adding data point " << pageid << " to the skyline set!" << endl;
                }
            } else   //current entry is an index, so all its entries are inserted into the min-heap
            {
                a_pageaccessed.append((void *) pageid);
                for (int i = 0; i < VirNode->m_usedSpace; i++) {
                    //compute the score of current entry w.r.t query q
                    for (int j = 0; j < dimen; j++)
                        pt[j] = VirNode->m_entry[i]->m_hc.getUpper()[j];

                    bool dominated = IsDominatedBy(dimen, pt, skylines, PG);
                    if (!dominated) {
                        mindist = minDist(pt, ORIGIN, dimen);
                        H0.insert(DblInt_Pair(mindist, VirNode->m_entry[i]->m_id));
                    } else {
                        PrunedNodes.push_back(VirNode->m_entry[i]->m_id);
                    }
                }
            }//end check whether it's leaf or index node
        } else {//current object is dominated; so we put it in either prunednodes or nonresultentry for later skyline update
            if (VirNode->isLeaf()) {
                if (VirNode->m_usedSpace > 1)  //leaf node
                {
                    PrunedNodes.push_back(pageid);
                }
            } else {   //index node
                PrunedNodes.push_back(pageid);
            }
        }

        delete VirNode;
        delete e0;
    }

    a_skylines.clear();
    for (Int_vIter = skylines.begin(); Int_vIter != skylines.end(); Int_vIter++) {
        a_skylines.insert(*Int_vIter);
        //cout << *Int_vIter << " ";
    }
    if (verbose)
        cout << endl;
    if (verbose)
        cout << "Number of skylines: " << a_skylines.size() << endl;

    return 1;
}

float S3::GetScore(float p[], float weight[], int dimen)    //compute the score
{
    double score = 0;

    for (int i = 0; i < dimen; i++)
        score += p[i] * weight[i];

    return score;
}

float S3::minDist(float p1[], float p2[], int dimen) {
    float mindist = 0;

    for (int i = 0; i < dimen; i++) {
        float dist = p1[i] - p2[i];
        mindist += (dist * dist);
    }
    return (float) sqrt(mindist);
}

//check whether a MBR is dominated by a skyline
bool S3::IsDominatedBy(
        const int dimen,
        const float pt[],
        vector<long> a_skylines,
        float *PG[]) {
    vector<long>::iterator Int_vIter;

    if (a_skylines.size() == 0) return false;

    for (Int_vIter = a_skylines.begin(); Int_vIter != a_skylines.end(); Int_vIter++) {
        long pid = *Int_vIter;
        float s[DMAX];

        bool dominated = true;
        for (int i = 0; i < dimen; i++) {
            if (PG[pid][i] + SIDELEN < pt[i]) {
                dominated = false;
                break;
            }
        }
        if (dominated)
            return dominated;
    }

    return false;
}

//check whether a point pt is dominated by another point pt0
bool S3::IsDominatedBy(
        const int dimen,
        const float pt[],
        const float pt0[]) {
    for (int i = 0; i < dimen; i++)
        if (pt0[i] < pt[i]) return false;

    return true;
}

void S3::OutputForQhalf(
        FILE *fout1,
        const int dimen,
        Query &a_query,
        int K,
        long int *FacetMin,
        long int *FacetMax,
        long int &FacetMinPID,
        long int &FacetMaxPID,
        float *PG[]) {
    int i, j, m;
    int NoOfHyperplanes = 0;
    float pFirst[DMAX], pSecond[DMAX];

    vector<long>::iterator Int_vIter;
    multimap<float, int>::reverse_iterator rpIter;
    multimap<float, int>::iterator pIter;

    NoOfHyperplanes = 2 * dimen + (K - 1) + 2;    //the constant '2' is the number of facets adjacent to theTopk

    fprintf(fout1, "%d 1\n", dimen);   //the dimensionality
    for (i = 0; i < dimen; i++)
        fprintf(fout1, "%f ", a_query.function[i]); //the feasible point is the query point itself
    fprintf(fout1, "\n");
    fprintf(fout1, "%d\n", dimen + 1);   //the dimensionality
    fprintf(fout1, "%d\n", NoOfHyperplanes);  //the total number of hyperplanes for intersection

    //output the 2*dimen bounding facets of the hypercube of the dataspace
    double Offset = 0;
    int Cooef = -1;
    for (i = 1; i <= 2; i++) {
        if (i == 2) {
            Offset = -CEILING + 20;
            Cooef = -Cooef;
        }
        for (j = 1; j <= dimen; j++) {
            for (m = 1; m <= dimen; m++)
                (m == j) ? fprintf(fout1, "%d  ", Cooef) : fprintf(fout1, "0  ");
            fprintf(fout1, "%f\n", Offset);
        }
    }

    //output the hyperplanes formed by consecutive pairs of top-k result points
    rpIter = a_query.score.rbegin();
    for (int i = 0; i < dimen; i++) pFirst[i] = PG[rpIter->second][i] + SIDELEN;
    //pFirst=S[rpIter->second];
    rpIter++;
    for (; rpIter != a_query.score.rend(); rpIter++) {
        for (int i = 0; i < dimen; i++) pSecond[i] = PG[rpIter->second][i] + SIDELEN;

        for (int i = 0; i < dimen; i++)
            fprintf(fout1, "%f ", -(pFirst[i] - pSecond[i]));  //for qhalf
        fprintf(fout1, "0\n");

        for (int i = 0; i < dimen; i++) pFirst[i] = pSecond[i];
    }
    pIter = a_query.score.begin();  //obtain the k-th point in the top-k result set
    for (int i = 0; i < dimen; i++) pFirst[i] = PG[pIter->second][i] + SIDELEN;


    //output the two hyperplanes formed by theTopk and the two facets adjacent to theTopk
    for (int i = 0; i < dimen; i++) fprintf(fout1, "%f ", -(pFirst[i] - FacetMin[i]));
    fprintf(fout1, "0\n");

    for (int i = 0; i < dimen; i++) fprintf(fout1, "%f ", -(pFirst[i] - FacetMax[i]));
    fprintf(fout1, "0\n");

    cout << "output data for halfspace intersection finished!" << endl;

    return;
}

void S3::OutputNonResultSetVanila(FILE *fout1, const int dimen, multimap<long, long> NonResultPtID, float *PG[]) {
    int i;
    multimap<long, long>::iterator IntInt_Iter;

    //fprintf(fout1, "%d\n", dimen);   //the dimensionality
    //fprintf(fout1, "%d\n", NonResultPtID.size());   //the number of skylines

    //output the the skyline points
    for (IntInt_Iter = NonResultPtID.begin(); IntInt_Iter != NonResultPtID.end(); IntInt_Iter++) {
        long Pos = IntInt_Iter->first;
        for (i = 0; i < dimen; i++)
            fprintf(fout1, "%f ", PG[Pos][i] + SIDELEN);
        fprintf(fout1, "\n");
    }
    //cout << "Nonresult set output finished!" << endl;
}

void S3::OutputNonResultSet(FILE *fout1, const int dimen, multimap<long, long> NonResultPtID, float *PG[]) {
    int i;
    vector<long>::iterator Int_vIter;

    set<long> FirstDpoints;
    set<long>::iterator Int_sIter;

    multimap<long, long>::iterator IntInt_Iter;

    multimap<double, long> RankNonResult[DMAX];
    multimap<double, long>::reverse_iterator DblInt_rIter;
    typedef multimap<double, long>::value_type DblInt_Pair;

    float pt[DMAX];

    for (IntInt_Iter = NonResultPtID.begin(); IntInt_Iter != NonResultPtID.end(); IntInt_Iter++) {
        //long Pos=IntVRN_Iter->first-MAXPAGEID;
        long Pos = IntInt_Iter->first;
        for (i = 0; i < dimen; i++) {
            pt[i] = PG[Pos][i] + SIDELEN;
            RankNonResult[i].insert(DblInt_Pair(pt[i], Pos));
        }
    }
    for (i = 0; i < dimen; i++) {
        //printf("size of the %d-th map: %d\n",i,RankSkyline[i].size());
        DblInt_rIter = RankNonResult[i].rbegin();
        int szSet = FirstDpoints.size();
        while (DblInt_rIter != RankNonResult[i].rend()) {
            FirstDpoints.insert(DblInt_rIter->second);
            if (szSet == FirstDpoints.size())
                DblInt_rIter++;
            else
                break;
        }
        //printf("size of FirstDpoints in %d-th dim: %d\n",i,FirstDpoints.size());
    }
    //output the first d data points who have the maximum value along any dimension
    for (Int_sIter = FirstDpoints.begin(); Int_sIter != FirstDpoints.end(); Int_sIter++) {
        //printf("FirstDpoints: %ld\n",(*Int_sIter));
        for (i = 0; i < dimen; i++)
            fprintf(fout1, "%lf  ", PG[(*Int_sIter)][i] + SIDELEN);
        fprintf(fout1, "\n");
    }

    //output the rest points to file
    for (IntInt_Iter = NonResultPtID.begin(); IntInt_Iter != NonResultPtID.end(); IntInt_Iter++) {
        long Pos = IntInt_Iter->first;
        Int_sIter = FirstDpoints.find(Pos);
        if (Int_sIter != FirstDpoints.end()) continue;
        for (i = 0; i < dimen; i++)
            fprintf(fout1, "%lf  ", PG[Pos][i] + SIDELEN);
        fprintf(fout1, "\n");
    }
}
