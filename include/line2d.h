/* ----------------------------------------------------------------------------
    This header file includes class Line2D declaration.
---------------------------------------------------------------------------- */

#ifndef LINE2D_DEFINED
#define LINE2D_DEFINED

#include "point.h"

#define INFTY 100000000

class Line2D {
protected:
    float m_slope;        // slope
    float m_intercept;    // y-intercept; if slope = INFTY, x-intercept
public:
    Line2D(const float a_slope, const float a_intercept);

    Line2D(const Line2D &a_line);

    virtual ~Line2D();

    //
    // properties
    float slope() const;

    float intercept() const;

    void normal(            // determine normal angle/distance with respect to pt
            Point &a_pt,
            float &a_angle, float &a_dist) const;

    //
    // assignment
    Line2D &operator=(const Line2D &a_line);

    //
    // manipulation
    Line2D rotate(const Point &a_pt) const;

    bool intersection(const Line2D &a_l, Point &pt) const;

    //
    // generate new line object
    static Line2D line(const Point &a_p1, const Point &a_p2);

    static const int CROSS, OVERLAP, PARALLEL;

    static int relationship(const Line2D &a_l1, const Line2D &a_l2);
};

#endif    // LINE2D_DEFINED


