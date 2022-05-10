#include "point.h"
#include <string.h>
#include <math.h>

// constructor/destructor
Point::Point(const int a_dimen, const float* a_coor):
m_dimen(a_dimen)
{
    if (a_coor != 0)
        memcpy(m_coor,a_coor,sizeof(float)*m_dimen);
    else
        memset(m_coor,0,sizeof(float)*m_dimen);
}

Point::Point(const Point& a_pt):
m_dimen(a_pt.m_dimen)
{
    memcpy(m_coor,a_pt.m_coor,sizeof(float)*m_dimen);
}

Point::~Point()
{}

// operation
Point& Point::operator=(const Point& a_pt)
{
    if (this != &a_pt)
        memcpy(m_coor, a_pt.m_coor, sizeof(float)*m_dimen);
    return *this;
}

Point Point::midpoint(const Point& a_p0, const Point& a_p1)
{
	float m[MAXDIMEN];
	for (int i=0; i<a_p0.m_dimen; i++)
		m[i] = (a_p0[i] + a_p1[i])/2;
	return Point(a_p0.m_dimen, m);
}


// search
float Point::operator[](const int a_i) const
{
    return m_coor[a_i];
}

// update
void Point::set(const int a_i, const float a_c)
{
    m_coor[a_i] = a_c;
}

// measures
float Point::distance(const Point& a_pt) const
{
    float dist=0;
    for (int d=0; d<m_dimen; d++)
    {
        float diff = (a_pt.m_coor[d] - m_coor[d]);
        dist += diff*diff;
    }
    return sqrt(dist);
}

float Point::angle(const Point& a_pt) const
{
    // this function considers 2d only
    double xdiff = m_coor[0] - a_pt.m_coor[0];
    double ydiff = m_coor[1] - a_pt.m_coor[1];
    bool xpos = xdiff >= 0;
    bool ypos = ydiff >= 0;
    xdiff = xdiff > 0 ? xdiff : -xdiff;
    ydiff = ydiff > 0 ? ydiff : -ydiff;

    double radian = atan(ydiff/xdiff);
	float degree = (float)(radian * 180 / PI);
    if (xpos)
    {
        if (ypos)                   // quadrant I (+ve x-axis and +ve y-axis)
            return degree;
        else                        // quadrant IV (-ve x-axis and +ve y-axis)
            return 360 - degree;
    }
    else
    {
        if (ypos)
            return 180-degree;		// quadrant III (-ve x-axis and -ve y-axis)
        else
            return degree + 180;	// quadrant IV (-ve x-axis and -ve y-axis)
    }


    /*
    if (xdiff > 0)
    {
        if (ydiff > 0)
            return atan(ydiff/xdiff);
        return atan(ydiff/xdiff)+2*PI;
    }
    else
    {
        if (ydiff > 0)
            return atan(ydiff/xdiff) + PI;
        return atan(ydiff/xdiff) + PI;
    }
    */

}

// comparison
bool Point::operator==(const Point& a_pt) const
{
    for (int d=0; d<m_dimen; d++)
        if (m_coor[d] != a_pt.m_coor[d])
            return false;
    return true;
}

// info
int Point::size(const int a_dimen)
{
    return sizeof(float)*a_dimen;
}
