#include "lineseg.h"
#include "point.h"
#include <iostream>

using namespace std;

int main(int argc, const char** argv)
{
    float c0[2], c1[2];
    c0[0] = c0[1] = 50;

    float ca[2], cb[2];
    ca[0] = ca[1] = 25;

    for (c1[0]=0; c1[0]<=100; c1[0]+=50)
    {
        for (c1[1]=0; c1[1]<=100; c1[1]+=50)
        {
            Point p0(2,c0);
            Point p1(2,c1);
            LineSeg l0(p0,p1);

            for (cb[0]=0; cb[0]<=100; cb[0]+=50)
            {
                for (cb[1]=0; cb[1]<=100; cb[1]+=50)
                {
                    Point pa(2,ca);
                    Point pb(2,cb);

                    LineSeg l1(pa,pb);

                    Point pt(2);
                    bool result = LineSeg::intersect(l0, l1, pt);

                    cout << "[(" << c0[0] << "," << c0[1] << "),(" << c1[0] << "," << c1[1] << ")] ";
                    cout << "[(" << ca[0] << "," << ca[1] << "),(" << cb[0] << "," << cb[1] << ")] ";
                    cout << "==> (" << pt[0] << "," << pt[1] << ")";
                    if (result == true)
                        cout << " ... valid";
                    cout << endl;
                }
            }
        }
    }

    return 0;
}