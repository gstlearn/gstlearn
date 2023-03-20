/******************************************************************************/
/* gstlearn C++ Library                                                       */
/*                                                                            */
/* Authors: <authors>                                                         */
/*                                                                            */
/* License: BSD 3 Clause                                                      */
/******************************************************************************/
#include "Polygon/PolySet.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/ASerializable.hpp"
#include "Basic/Utilities.hpp"
#include "Basic/PolyLine2D.hpp"

PolySet::PolySet(const VectorDouble& x,
                 const VectorDouble& y,
                 double zmin,
                 double zmax)
    : PolyLine2D(x, y),
      _zmin(TEST),
      _zmax(TEST)
{
  init(x,y,zmin,zmax);
}

PolySet::PolySet(const PolySet& r)
    : PolyLine2D(r),
      _zmin(r._zmin),
      _zmax(r._zmax)
{
}

PolySet& PolySet::operator=(const PolySet& r)
{
  if (this != &r)
  {
    PolyLine2D::operator=(r);
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
  PolyLine2D::init(x,y);

  _zmin  = zmin;
  _zmax  = zmax;
}

String PolySet::toString(const AStringFormat* strfmt) const
{
  std::stringstream sstr;

  sstr << PolyLine2D::toString(strfmt);

  if (!FFFF(_zmin) || !FFFF(_zmax)) sstr << toInterval(_zmin, _zmax);

  return sstr.str();
}

void PolySet::getExtension(double *xmin,
                           double *xmax,
                           double *ymin,
                           double *ymax) const
{
  *xmin = getXmin();
  *ymin = getYmin();
  *xmax = getXmax();
  *ymax = getYmax();
}

double PolySet::getSurface() const
{
  int np = getNPoints();
  double x0 = getX(0);
  double y0 = getY(0);
  double surface = 0.;
  for (int i=1; i<np-2; i++)
  {
    double x1 = getX(i) - x0;
    double y1 = getY(i) - y0;
    double x2 = getX(i + 1) - x0;
    double y2 = getY(i + 1) - y0;
    surface += 0.5 * ((x1 * y2) - (x2 * y1));
  }

  // Check if the PolySet is closed

  if (! _isClosed())
  {
    double x1 = getX(np-1) - x0;
    double y1 = getY(np-1) - y0;
    double x2 = getX(0) - x0;
    double y2 = getY(0) - y0;
    surface += 0.5 * ((x1 * y2) - (x2 * y1));
  }

  surface = ABS(surface);
  return(surface);
}

bool PolySet::_serialize(std::ostream& os, bool verbose) const
{
  if (getNPoints() <= 0) return false;
  bool ret = true;
  ret = ret && _recordWrite<double>(os, "Z-Minimum", _zmin);
  ret = ret && _recordWrite<double>(os, "Z-Maximum", _zmax);
  ret = ret && PolyLine2D::_serialize(os, verbose);
  return ret;
}

bool PolySet::_deserialize(std::istream& is, bool verbose)
{
  _zmin = TEST;
  _zmax = TEST;
  bool ret = true;
  ret = ret && _recordRead<double>(is, "Z-Minimum", _zmin);
  ret = ret && _recordRead<double>(is, "Z-Maximum", _zmax);
  ret = ret && PolyLine2D::_deserialize(is, verbose);
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
  polyset = new PolySet();
  bool success = false;
  if (polyset->_fileOpenRead(neutralFilename, is, verbose))
  {
    success = polyset->deserialize(is, verbose);
  }

  if (! success)
  {
    delete polyset;
    polyset = nullptr;
  }
  return polyset;
}

bool PolySet::_isClosed() const
{
  int nvert = getNPoints();
  if (ABS(getX(0) - getX(nvert-1)) > EPSILON5 ||
      ABS(getY(0) - getY(nvert-1)) > EPSILON5) return false;
  return true;
}

/**
 * Close the PolySet if necessary
 */
void PolySet::closePolySet()
{
  if (!_isClosed()) addPoint(getX(0), getY(0));
}

/****************************************************************************/
/*!
 **  Check if one point belongs to a 2-D polyset
 **
 ** \return  True if the point belongs to the polygon; False otherwise
 **
 ** \param[in]  xx  array of point coordinates of the point along X
 ** \param[in]  yy  array of point coordinates of the point along Y
 **
 *****************************************************************************/
bool PolySet::inside(double xx, double yy)
{
  double dx, dy, xinter;
  int inter, j, sel;

  inter = 0;
  int np = getNPoints();

  /* Loop on the polygon vertices */

  for (j = 0; j < np - 1; j++)
  {

    /* Horizontal segment */

    dy = getY(j + 1) - getY(j);
    if (dy == 0 && yy == getY(j))
    {
      if (getX(j + 1) > getX(j) && xx > getX(j) && xx < getX(j + 1))
      {
        inter = 1;
        continue;
      }
      if (getX(j + 1) < getX(j) && xx < getX(j) && xx > getX(j + 1))
      {
        inter = 1;
        continue;
      }
    }

    /* One vertex below and one vertex above: point distinct from segment */

    if (dy != 0 && ((getY(j) > yy && getY(j + 1) < yy)
        || (getY(j) < yy && getY(j + 1) > yy)))
    {
      dx = getX(j + 1) - getX(j);
      xinter = (dx * yy + dy * getX(j) - dx * getY(j)) / dy;
      if (xinter > xx) inter++;
    }

    /* One vertex below and one vertex above: point belongs to segment */

    if (dy != 0 && ((getY(j) > yy && getY(j + 1) < yy)
        || (getY(j) < yy && getY(j + 1) > yy)))
    {
      dx = getX(j + 1) - getX(j);
      xinter = (dx * yy + dy * getX(j) - dx * getY(j)) / dy;
      if (xinter == xx)
      {
        inter = 1;
        continue;
      }
    }

    /* Point is in contact with the highest vertex */

    if (yy == getY(j) && getY(j) > getY(j + 1) && xx < getX(j)) inter++;
    if (yy == getY(j + 1) && getY(j + 1) > getY(j) && xx < getX(j + 1)) inter++;

    /* Point coincides with a vertex */

    if (xx == getX(j) && yy == getY(j))
    {
      inter = 1;
      continue;
    }
  }
  sel = ((inter % 2) != 0);
  return (sel);
}

/****************************************************************************/
/*!
 **  Check if one point belongs to a vertical interval of a (limited) polyset
 **
 ** \return  True if the point belongs to the polygon; False otherwise
 **
 ** \param[in]  zz   array of point coordinates of the point along Z or TEST
 **
 *****************************************************************************/
bool PolySet::inside3D(double zz)
{
  if (FFFF(zz)) return true;
  if (!FFFF(_zmin) && zz < _zmin) return false;
  if (!FFFF(_zmax) && zz > _zmax) return false;
  return true;
}

