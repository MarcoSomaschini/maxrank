#include "stdio.h"
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <math.h>
#include <time.h>
#include <fstream>
#include <iostream>

#define MAXDIM 20     
using namespace std;

int main(int argc, char** argv)
{
    fstream  fin;
    ofstream fout;

    int i,j;
    long LineCnt=0;
    double point[MAXDIM];
    double min[MAXDIM],max[MAXDIM];

    long D;
    int ExpandRatio;
   
    if (argc != 5) 
    {
	printf("Syntax: %s [infile] [outfile] [dimensionality] [expand ratio]\n", argv[0]);
        return -1;
    }

    D=atoi(argv[3]);
    ExpandRatio=atoi(argv[4]);
     
    for (i=0;i<D;i++)
    {
         min[i]=99999999;
         max[i]=-99999999;
    }

    fin.open(argv[1], ios::in);
    while (true)
    {   
         fin >> point[0];
         if (fin.eof()) break;                

         LineCnt++;
         if ((LineCnt%5000)==0) cout << LineCnt << " lines have been read..." << endl;

         if (point[0]<min[0]) min[0]=point[0];
         if (point[0]>max[0]) max[0]=point[0];

         for (i=1;i<D;i++)
         {
              fin >> point[i];
              if (point[i]<min[i]) min[i]=point[i];
              if (point[i]>max[i]) max[i]=point[i];
         }
    }
    fin.close();
    
    cout << "min   max:"<< endl;
    for (i=0;i<D;i++)
    {
         cout << min[i] << " ";
         cout << max[i] << " "<< endl;
    }

    //normalization
    long zeroLines=0;
    LineCnt = 0;
    fin.open(argv[1], ios::in);
    fout.open(argv[2], ios::out);
    while (true)
    {   
         fin >> point[0];
         if (fin.eof()) break;         

         LineCnt++;
         if ((LineCnt%5000)==0) cout << LineCnt << " lines have been normalized..." << endl;

         double tmpSum;
         point[0]=((point[0]-min[0])/(max[0]-min[0]))*ExpandRatio;
         tmpSum=point[0];
         for (i=1;i<D;i++)
         {
              fin >> point[i];
              point[i]=((point[i]-min[i])/(max[i]-min[i]))*ExpandRatio;
              tmpSum+=point[i];
         }

         if (tmpSum>0)
         {
             for (i=0;i<D;i++)
                  fout << point[i] << " ";
             fout << endl;         
         }
         else
             zeroLines++;  
    } 
    fin.close();
    fout.close();
    cout << "Normalization finished!  " << zeroLines << " lines of all zero values are discarded!" << endl;

    return 0;
}
