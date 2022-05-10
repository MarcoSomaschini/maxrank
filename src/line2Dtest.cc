#include "line2d.h"
#include "point.h"
#include <iostream>

using namespace std;

int main(int argc, const char** argv)
{
    float c0[2], c1[2];
    c0[0] = c0[1] = 50;

    float ca[2];
    ca[0] = 20; ca[1] = 30;
    Point origin(2,ca);

    for (c1[0]=0; c1[0]<=100; c1[0]+=50)
    {
        for (c1[1]=0; c1[1]<=100; c1[1]+=50)
        {
            Point p0(2,c0);
            Point p1(2,c1);
            Line2D line = Line2D::line(p0, p1);

            float angle=0;
            float distance=0;
            line.normal(origin, angle, distance);

            cout << "at (" << ca[0] << "," << ca[1] << ") ";
            cout << "(" << c0[0] << "," << c0[1] << "),(" << c1[0] << "," << c1[1] << ") ";
            cout << "ang: " << angle << " dis: " << distance << endl;
        }
    }

    return 0;
}
