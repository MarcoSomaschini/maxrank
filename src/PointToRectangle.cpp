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
    double point[MAXDIM];
    int xTimes[MAXDIM];
    long D;
    float Len;
    int Prefix,StartID;
    long int ID=0;
    int Addup=0;
    long LineCnt=0;
    int mulFactor=1;
   
    if (argc != 9) 
    {
	printf("Syntax: %s [infile] [outfile] [dimensionality] [SideLength] [prefix] [StartID] [mult. factor] [add up]\n", argv[0]);
        return -1;
    }

    D=atoi(argv[3]);
    Len=atof(argv[4]);
    Prefix=atoi(argv[5]);    
    StartID=atoi(argv[6]);        
    mulFactor=atoi(argv[7]);
    Addup=atoi(argv[8]);        

    xTimes[0]=1;
    xTimes[1]=1;
    xTimes[2]=1000;
    xTimes[3]=1000;
    xTimes[4]=1000;
 
    fin.open(argv[1], ios::in);
    fout.open(argv[2], ios::out);
    while (true)
    {   
         fin >> point[0];
         if (fin.eof()) break;         

         LineCnt++;
         if ((LineCnt%5000)==0) cout << LineCnt << " lines have been processed..." << endl;

         for (i=1;i<D;i++)
              fin >> point[i];

         if (Prefix>=0) fout << Prefix << " ";
         if (StartID>=0) fout << StartID++ << " ";          
         for (i=0;i<D;i++)         
             fout << (point[i]*mulFactor+Addup)-Len/2.0 << " ";
         
         for (i=0;i<D;i++)         
             fout << (point[i]*mulFactor+Addup)+Len/2.0 << " ";
         fout << endl;         
    }
    fin.close();
    fout.close();
    cout << "Expansion from point to rectangle is finished!" << endl;

    return 0;
}
