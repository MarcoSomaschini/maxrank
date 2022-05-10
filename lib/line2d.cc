#include "line2d.h"

Line2D::Line2D(const float a_slope, const float a_intercept):
m_slope(a_slope), m_intercept(a_intercept)
{}

Line2D::Line2D(const Line2D& a_line):
m_slope(a_line.m_slope), m_intercept(a_line.m_intercept)
{}

Line2D::~Line2D()
{}

//-----------------------------------------------------------------------------
// properties
//-----------------------------------------------------------------------------
float Line2D::slope() const
{
	return m_slope;
}

//-----------------------------------------------------------------------------
// properties
//-----------------------------------------------------------------------------
float Line2D::intercept() const
{
	return m_intercept;
}

//-----------------------------------------------------------------------------
// assignment
//-----------------------------------------------------------------------------
Line2D& Line2D::operator=(const Line2D& a_line)
{
	m_slope = a_line.m_slope;
	m_intercept = a_line.m_intercept;
	return *this;
}

//-----------------------------------------------------------------------------
// manipulation
//-----------------------------------------------------------------------------
Line2D Line2D::rotate(const Point& a_pt) const
{
	float slope = INFTY;
	float i     = INFTY;
	if (m_slope == 0)
	{
		slope = INFTY;
		i     = a_pt[0];
	}
	else if (m_slope == INFTY)
	{
		slope = 0;
		i = a_pt[1];
	}
	else
	{
		slope = -1/m_slope;
		i = a_pt[1] - slope * a_pt[0];
	}
	return Line2D(slope,i);
}

bool Line2D::intersection(const Line2D& a_l, Point& a_pt) const
{
	if (m_slope == a_l.m_slope)
		return false;

	if (m_slope == 0)						// A: y = c_y_a;
	{
		a_pt.set(1, m_intercept);			// y = c_y_a
		if (a_l.m_slope == INFTY)			// B: x = c_x_b;
			a_pt.set(0,a_l.m_intercept);	// x = c_x_b
		else								// B: y = m_b * x + c_y_b
			a_pt.set(0,(m_intercept - a_l.m_intercept)/a_l.m_slope);	// x = (c_y_a - c_y_b)/m_b
	}
	else if (m_slope == INFTY)				// A: x = c_x_a;
	{
		a_pt.set(0,m_intercept);			// x = c_x_a
		if (a_l.m_slope == 0)				// B: x = c_y_b;
			a_pt.set(1,a_l.m_intercept);	// y = c_y_b
		else							    // B: y = m_b * x + c_y_b;
			a_pt.set(1,a_l.m_slope * m_intercept + a_l.m_intercept);	// y = m_b * c_x_a + c_y_b
	}
	else								    // A: y = m_a * x + c_y_a;
	{
		if (a_l.m_slope == 0)			    // B: y = c_y_b
		{
			a_pt.set(0,(a_l.m_intercept - m_intercept)/m_slope);		// x = (c_y_b - c_y_a)/m_a
			a_pt.set(1,a_l.m_intercept);								// y = c_y_b
		}
		else if (a_l.m_slope == INFTY)	    // B: x = c_x_b
		{
			a_pt.set(0,a_l.m_intercept);								// x = c_x_b
			a_pt.set(1,m_slope * a_l.m_intercept + m_intercept);		// y = m_a * c_x_b + c_y_a
		}
		else							    // B: y = m_b * x + c_y_b;
		{
			a_pt.set(0,(m_intercept - a_l.m_intercept) /
				(a_l.m_slope - m_slope));								// x = (c_y_a - c_y_b) / (m_b - m_a);
			a_pt.set(1,m_slope * a_pt[0] + m_intercept);				// y = m_a * x + c_y_a;
		}
	}
	return true;
}

//-----------------------------------------------------------------------------
// determine the normal with respect to pt
//-----------------------------------------------------------------------------
void Line2D::normal(Point& a_pt, float& a_angle, float& a_dist) const
{
    // formulate a percendicular line
    float slope, intercept;
    if (m_slope == INFTY)
    {
        slope = 0;
        intercept = a_pt[1];    // y-intercept
    }
    else if (m_slope == 0)
    {
        slope = INFTY;
        intercept = a_pt[0];    // x-intercept
    }
    else
    {
        slope = -1/m_slope;
        intercept = a_pt[1] - slope * a_pt[0];
    }

    // find an intersection
    Line2D l(slope, intercept);
    Point nn(2);
    Line2D::intersection(l, nn);

    // find an angle and a distance w.r.t a_pt
    a_angle = nn.angle(a_pt);
    a_dist = nn.distance(a_pt);
}

//-----------------------------------------------------------------------------
// generate new line object
//-----------------------------------------------------------------------------
Line2D Line2D::line(const Point& a_p0, const Point& a_p1)
{
	float deltax = a_p0[0] - a_p1[0];
	float deltay = a_p0[1] - a_p1[1];
	float slope;
	float inter;
    if (deltay == 0)
    {
        if (deltax == 0)
        {
            // identical point?
            slope = INFTY;
            inter = INFTY;
        }
        else
        {
            slope = 0;
            inter = a_p0[1];
        }
    }
    else
    {
        if (deltax == 0)
        {
            slope = INFTY;
            inter = a_p0[0];
        }
        else
        {
		    slope = deltay / deltax;
		    inter = a_p0[1] - slope * a_p0[0];
        }
    }
	return Line2D(slope,inter);

    /*
	float deltax = a_p0[0] - a_p1[0];
	float deltay = a_p0[1] - a_p1[1];
	float slope = INFTY;
	float y     = INFTY;
	if (deltax != 0)
	{
		slope = deltay / deltax;
		y = a_p0[1] - slope * a_p0[0];
	}
	return Line2D(slope,y);
    */
}

const int Line2D::CROSS = 0;
const int Line2D::OVERLAP = 1;
const int Line2D::PARALLEL = 2;

//-----------------------------------------------------------------------------
// determine the relationship between two lines, that can be
// 0. CROSS
// 1. OVERLAP
// 2. PARALLEL
//-----------------------------------------------------------------------------
int Line2D::relationship(const Line2D& a_l1, const Line2D& a_l2)
{
    if (a_l1.m_slope == a_l2.m_slope)
    {
        if (a_l1.m_intercept == a_l2.m_intercept)
            return OVERLAP;
        else
            return PARALLEL;
    }
    return CROSS;
}
