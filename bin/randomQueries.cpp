#include "stdlib.h"
#include "stdio.h"
#include "string.h"
#include "math.h"
#include <climits>
#include <iostream>

using namespace std;

int main(int argc, char** argv)
{
    if (argc<4)
    {
        cout << "Error! Input parameters are not correct!" << endl;
        cout << "Usage: " << argv[0] << ": [#queries] [Upperlimit] [output file]" << endl;
        return -1;
    }
    
    int n=atoi(argv[1]);
    long int UL=atoi(argv[2]);
    
    FILE *fp_out;
    
    fp_out=fopen(argv[3],"w+");
    if (fp_out==NULL)
    {
       cout << "There is an error in openning file!" << endl;
       return -1;
    }
    
    srand(time(NULL));
    for (long int i=0;i<n;i++)
    {
        long int idx=(long int)(1 + (UL-1)*(float(rand())/RAND_MAX));
        fprintf(fp_out,"%ld\n",idx);
    }
    fclose(fp_out);
    
    cout << n << " random queries generated!" << endl;
  
    return 1; 
}
