/*---------------------------------------------------------------------------- 
    This header file includes class Point declaration, that represents
    d-dimensional points
---------------------------------------------------------------------------- */

#ifndef POINT_DEFINED
#define POINT_DEFINED

#define PI 3.1415926535897932384626433832795f

#define MAXDIMEN    16   // this is used to allocate for max dimensionality

class Point {
// data members
public:
    const int m_dimen;            // dimensionality
protected:
    float m_coor[MAXDIMEN];   // dimensional values
// methods
public:
    // constructor/destructor
    Point(const int a_dimen, const float *a_coor = 0);

    Point(const Point &a_pt);

    virtual ~Point();

    //
    // operation
    Point &operator=(const Point &a_pt);        // assignment
    //
    // search
    float operator[](const int a_i) const;      // return d-dimensional value
    static Point midpoint(const Point &a_p0, const Point &a_p1);

    //
    // update
    void set(const int a_i, const float a_c);   // update
    //
    // measures
    float distance(const Point &a_pt) const;    // distance to point a_pt
    float angle(const Point &a_pt) const;       // polar angle (radian) with
    // respect to point a_pt
    //
    // comparison
    bool operator==(const Point &a_pt) const;   // if it equals a point
    //
    // info
    static int size(const int a_dimen);         // storage size (bytes)
};

#endif // POINT_DEFINED

