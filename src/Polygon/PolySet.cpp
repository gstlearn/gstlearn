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
#include "Basic/ASerializable.hpp"
#include "Basic/Utilities.hpp"
#include "geoslib_f.h"

PolySet::PolySet()
  : AStringable(),
    ASerializable(),
    _x(0)
  , _y(0)
  , _zmin(TEST)
  , _zmax(TEST)
{
}

PolySet::PolySet(const VectorDouble& x,
                 const VectorDouble& y,
                 double zmin,
                 double zmax)
  : AStringable(),
    ASerializable(),
    _x(0)
  , _y(0)
  , _zmin(TEST)
  , _zmax(TEST)
{
  init(x,y,zmin,zmax);
}

PolySet::PolySet(const PolySet& r)
    : AStringable(r),
      ASerializable(r),
      _x(r._x),
      _y(r._y),
      _zmin(r._zmin),
      _zmax(r._zmax)
{
}

PolySet& PolySet::operator=(const PolySet& r)
{
  if (this != &r)
  {
    AStringable::operator=(r);
    ASerializable::operator=(r);
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

  /* Load the new PolySet */

  _x.resize(nvert,0);
  _y.resize(nvert,0);

  /* Copy the arrays */

  for (int i=0; i<nvert; i++)
  {
    _x[i] = x[i];
    _y[i] = y[i];
  }

  // Check if the PolySet must be closed

  closePolySet();

  _zmin  = zmin;
  _zmax  = zmax;
}

String PolySet::toString(const AStringFormat* strfmt) const
{
  std::stringstream sstr;

  int nvert = static_cast<int> (_x.size());
  sstr << "Number of vertices = " << nvert << std::endl;

  AStringFormat sf;
  if (strfmt != nullptr) sf = *strfmt;
  if (sf.getLevel() > 0)
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

int PolySet::_serialize(std::ostream& os, bool /*verbose*/) const
{
  if (getNVertices() <= 0) return 0;
  bool ret = _recordWrite<int>(os, "Number of Vertices", getNVertices());

  for (int i = 0; i < getNVertices(); i++)
  {
    ret = ret && _recordWrite<double>(os, "", getX(i));
    ret = ret && _recordWrite<double>(os, "", getY(i));
    ret = ret && _commentWrite(os, "");
  }
  return ret ? 0 : 1;
}

int PolySet::_deserialize(std::istream& is, bool /*verbose*/)
{
  int nvert;
  double zmin = TEST;
  double zmax = TEST;

  bool ret = _recordRead<int>(is, "Number of Vertices", nvert);
  VectorDouble x(nvert);
  VectorDouble y(nvert);

  /* Loop on the Vertices */

  for (int i = 0; i < nvert; i++)
  {
    ret = ret && _recordRead<double>(is, "X-Coordinate of a Polyset", x[i]);
    ret = ret && _recordRead<double>(is, "Y-Coordinate of a Polyset", y[i]);
  }
  if (! ret) return 1;

  /* Add the polyset */

  init(x,y,zmin,zmax);

  return 0;
}

int PolySet::dumpToNF(const String& neutralFilename, bool verbose) const
{
  std::ofstream os;
  int ret = 1;
  if (_fileOpenWrite(neutralFilename, "PolySet", os, verbose))
  {
    ret = _serialize(os, verbose);
    if (ret && verbose) messerr("Problem writing in the Neutral File.");
    os.close();
  }
  return ret;
}

PolySet* PolySet::create()
{
  return new PolySet();
}

PolySet* PolySet::createFromNF(const String& neutralFilename, bool verbose)
{
  PolySet* polyset = nullptr;
  std::ifstream is;
  if (_fileOpenRead(neutralFilename, "PolySet", is, verbose))
  {
    polyset = new PolySet();
    if (polyset->_deserialize(is, verbose))
    {
      if (verbose) messerr("Problem reading the Neutral File.");
      delete polyset;
      polyset = nullptr;
    }
    is.close();
  }
  return polyset;
}

/**
 * Check if the PolySet must be closed
 */
void PolySet::closePolySet()
{
  int nvert = static_cast<int> (_x.size());

  bool flag_close = false;
  if (ABS(_x[0] - _x[nvert-1]) > EPSILON5 ||
      ABS(_y[0] - _y[nvert-1]) > EPSILON5) flag_close = true;
  if (! flag_close) return;

  // Duplicate the first point at the end of the PolySet

  int nvert1 = nvert + 1;
  _x.resize(nvert1);
  _y.resize(nvert1);

  _x[nvert] = _x[0];
  _y[nvert] = _y[0];
}
