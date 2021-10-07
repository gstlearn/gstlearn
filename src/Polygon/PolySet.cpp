/******************************************************************************/
/* COPYRIGHT ARMINES, ALL RIGHTS RESERVED                                     */
/*                                                                            */
/* THE CONTENT OF THIS WORK CONTAINS CONFIDENTIAL AND PROPRIETARY             */
/* INFORMATION OF ARMINES. ANY DUPLICATION, MODIFICATION,                     */
/* DISTRIBUTION, OR DISCLOSURE IN ANY FORM, IN WHOLE, OR IN PART, IS STRICTLY */
/* PROHIBITED WITHOUT THE PRIOR EXPRESS WRITTEN PERMISSION OF ARMINES         */
/*                                                                            */
/* TAG_SOURCE_CG                                                              */
/******************************************************************************/
#include "Polygon/PolySet.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/Utilities.hpp"
#include "geoslib_f.h"

PolySet::PolySet()
  : _x(0)
  , _y(0)
  , _zmin(TEST)
  , _zmax(TEST)
{
}

PolySet::PolySet(const VectorDouble& x,
                 const VectorDouble& y,
                 double zmin,
                 double zmax)
  : _x(0)
  , _y(0)
  , _zmin(TEST)
  , _zmax(TEST)
{
  init(x,y,zmin,zmax);
}

PolySet::PolySet(const PolySet& r)
    : _x(r._x),
      _y(r._y),
      _zmin(r._zmin),
      _zmax(r._zmax)
{
}

PolySet& PolySet::operator=(const PolySet& r)
{
  if (this != &r)
  {
    _x = r._x;
    _y = r._y;
    _zmin = r._zmin;
    _zmax = r._zmax;
  }
  return *this;
}

PolySet::~PolySet()
{
}

void PolySet::init(const VectorDouble& x,
                   const VectorDouble& y,
                   double zmin,
                   double zmax)
{
  int nvert = static_cast<int> (x.size());

  /* Check if the polygon must be closed */

  bool flag_close = false;
  if (ABS(x[0] - x[nvert-1]) > EPSILON5 ||
      ABS(y[0] - y[nvert-1]) > EPSILON5) flag_close = true;
  int nvert1 = nvert + flag_close;

  /* Load the new PolySet */

  _x.resize(nvert1,0);
  _y.resize(nvert1,0);

  /* Copy the arrays */

  for (int i=0; i<nvert; i++)
  {
    _x[i] = x[i];
    _y[i] = y[i];
  }
  if (flag_close)
  {
    _x[nvert] = x[0];
    _y[nvert] = y[0];
  }
  _zmin  = zmin;
  _zmax  = zmax;
}

String PolySet::toString(int level) const
{
  std::stringstream sstr;

  int nvert = static_cast<int> (_x.size());
  sstr << "Number of vertices = " << nvert << std::endl;

  if (! IFFFF(level) && level > 0)
  {
    VectorDouble tab = VectorDouble(2 * nvert);
    for (int i = 0; i < nvert; i++)
    {
      tab[i] = _x[i];
      tab[i + nvert] = _y[i];
    }
    sstr << toMatrix("Polygon Vertex Coordinates", VectorString(), VectorString(),
                     true, 2, nvert, tab);

    if (!FFFF(_zmin) || !FFFF(_zmax))
      sstr << toInterval(_zmin, _zmax);
  }
  return sstr.str();
}

void PolySet::getExtension(double *xmin,
                           double *xmax,
                           double *ymin,
                           double *ymax) const
{
  *xmin = ut_vector_min(_x);
  *ymin = ut_vector_min(_y);
  *xmax = ut_vector_max(_x);
  *ymax = ut_vector_max(_y);
}

double PolySet::getSurface() const
{
  int np = getNVertices();
  double x0 = _x[0];
  double y0 = _y[0];
  double surface = 0.;
  for (int i=1; i<np-2; i++)
  {
    double x1 = _x[i] - x0;
    double y1 = _y[i] - y0;
    double x2 = _x[i + 1] - x0;
    double y2 = _y[i + 1] - y0;
    surface += 0.5 * ((x1 * y2) - (x2 * y1));
  }
  surface = ABS(surface);
  return(surface);
}

