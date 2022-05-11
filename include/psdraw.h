/* ----------------------------------------------------------------------------
    This header file includes class PSDraw declaration.
    This class is used to draw points, lines, polygons in an EPS file
---------------------------------------------------------------------------- */

#ifndef PSDRAW_DEFINED
#define PSDRAW_DEFINED

class Point;

class Hypercube;

#include <fstream>

using namespace std;

class PSDraw {
// data members
protected:
    fstream m_f;
    const float m_cwminx;               // clipping window
    const float m_cwmaxx;
    const float m_cwminy;
    const float m_cwmaxy;
    const int m_minx;                 // bounding box in EPS
    const int m_miny;
    const int m_maxx;
    const int m_maxy;
    const float m_xscale;
    const float m_yscale;
// methods
public:
    // constructor/destructor
    PSDraw(                             // create an EPS file with a specified
            const char *a_filename,         // clipping window
            float a_cwminx, float a_cwminy,
            float a_cwmaxx, float a_cwmaxy,
            int a_minx = 50, int a_miny = 50,   // bounding box on a page
            int a_maxx = 550, int a_maxy = 550);// default: 50, 50, 550, 550
    virtual ~PSDraw();

    void point(                         // point
            float a_x, float a_y,
            float a_lineGrayScale = 0);

    void line(                          // line
            float a_xstart, float a_ystart,
            float a_xend, float a_yend,
            float a_lineGrayScale = 0);

    void circle(                        // circle
            float a_xcenter, float a_ycenter,
            float a_radius,
            float a_lineGrayScale = 0);

    void box(                           // rectangular box
            float a_xmin, float a_ymin,
            float a_xmax, float a_ymax,
            float a_lineGrayScale = 0);

    void polygon(                       // polygon with n points
            float *a_x, float *a_y, const int a_n,
            float a_lineGrayScale = 0);

    void point(                         // point on selected x-,y- dimen
            const Point &pt,
            const int x, const int y,
            float a_lineGrayScale = 0);

    void box(                           // rect. box on selected x-,y- dimen
            const Hypercube &hc,
            const int x, const int y,
            float a_lineGrayScale = 0);
};

#endif // PSDRAW_DEFINED

