/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES Paris / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://gstlearn.org                                              */
/* License: BSD 3-clause                                                      */
/*                                                                            */
/******************************************************************************/
#include "Polygon/PolyElem.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/Utilities.hpp"
#include "Basic/PolyLine2D.hpp"

PolyElem::PolyElem(const VectorDouble& x,
                 const VectorDouble& y,
                 double zmin,
                 double zmax)
    : PolyLine2D(x, y),
      _zmin(TEST),
      _zmax(TEST)
{
  init(x,y,zmin,zmax);
}

PolyElem::PolyElem(const PolyElem& r)
    : PolyLine2D(r),
      _zmin(r._zmin),
      _zmax(r._zmax)
{
}

PolyElem& PolyElem::operator=(const PolyElem& r)
{
  if (this != &r)
  {
    PolyLine2D::operator=(r);
    _zmin = r._zmin;
    _zmax = r._zmax;
  }
  return *this;
}

PolyElem::~PolyElem()
{
}

void PolyElem::init(const VectorDouble& x,
                   const VectorDouble& y,
                   double zmin,
                   double zmax)
{
  PolyLine2D::init(x,y);

  _zmin  = zmin;
  _zmax  = zmax;
}

String PolyElem::toString(const AStringFormat* strfmt) const
{
  std::stringstream sstr;

  sstr << PolyLine2D::toString(strfmt);

  if (!FFFF(_zmin) || !FFFF(_zmax)) sstr << toInterval(_zmin, _zmax);

  return sstr.str();
}

void PolyElem::getExtension(double *xmin,
                           double *xmax,
                           double *ymin,
                           double *ymax) const
{
  *xmin = getXmin();
  *ymin = getYmin();
  *xmax = getXmax();
  *ymax = getYmax();
}

double PolyElem::getSurface() const
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

  // Check if the PolyElem is closed

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

bool PolyElem::_serialize(std::ostream& os, bool verbose) const
{
  if (getNPoints() <= 0) return false;
  bool ret = true;
  ret = ret && _recordWrite<double>(os, "Z-Minimum", _zmin);
  ret = ret && _recordWrite<double>(os, "Z-Maximum", _zmax);
  ret = ret && PolyLine2D::_serialize(os, verbose);
  return ret;
}

bool PolyElem::_deserialize(std::istream& is, bool verbose)
{
  _zmin = TEST;
  _zmax = TEST;
  bool ret = true;
  ret = ret && _recordRead<double>(is, "Z-Minimum", _zmin);
  ret = ret && _recordRead<double>(is, "Z-Maximum", _zmax);
  ret = ret && PolyLine2D::_deserialize(is, verbose);
  return ret;
}

PolyElem* PolyElem::create()
{
  return new PolyElem();
}

PolyElem* PolyElem::createFromNF(const String& neutralFilename, bool verbose)
{
  PolyElem* polyelem = nullptr;
  std::ifstream is;
  polyelem = new PolyElem();
  bool success = false;
  if (polyelem->_fileOpenRead(neutralFilename, is, verbose))
  {
    success = polyelem->deserialize(is, verbose);
  }

  if (! success)
  {
    delete polyelem;
    polyelem = nullptr;
  }
  return polyelem;
}

bool PolyElem::_isClosed() const
{
  int nvert = getNPoints();
  return (ABS(getX(0) - getX(nvert - 1)) <= EPSILON5 &&
          ABS(getY(0) - getY(nvert - 1)) <= EPSILON5);
}

/**
 * Close the PolyElem if necessary
 */
void PolyElem::closePolyElem()
{
  if (!_isClosed()) addPoint(getX(0), getY(0));
}

/****************************************************************************/
/*!
 **  Check if one point belongs to a 2-D polyelem
 **
 ** \return  True if the point belongs to the polygon; False otherwise
 **
 ** \param[in]  coor  Vector giving the coordinates of the target point
 **
 *****************************************************************************/
bool PolyElem::inside(const VectorDouble& coor)
{
  double dx, dy, xj0, xj1, yj0, yj1, xinter;

  int inter = 0;
  int np = getNPoints();
  double xx = coor[0];
  double yy = coor[1];

  /* Loop on the polygon vertices */

  for (int j = 0; j < np - 1; j++)
  {
    xj0 = getX(j);
    xj1 = getX(j+1);
    yj0 = getY(j);
    yj1 = getY(j+1);

    dx = xj1 - xj0;
    dy = yj1 - yj0;

    /* Horizontal segment */

    if (dy == 0 && yy == yj0)
    {
      if (xj1 > xj0 && xx > xj0 && xx < xj1)
      {
        inter = 1;
        continue;
      }
      if (xj1 < xj0 && xx < xj0 && xx > xj1)
      {
        inter = 1;
        continue;
      }
    }

    /* One vertex below and one vertex above */

    if (dy != 0 && ( (yj0 > yy && yj1 < yy) || (yj0 < yy && yj1 > yy) ))
    {
      xinter = (dx * yy + dy * xj0 - dx * yj0) / dy;

      /* Point distinct from segment */

      if (xinter > xx) inter++;

      /* Point belongs to segment */

      if (xinter == xx)
      {
        inter = 1;
        continue;
      }
    }

    /* Point is in contact with the highest vertex */

    if (yy == yj0 && yj0 > yj1 && xx < xj0) inter++;
    if (yy == yj1 && yj1 > yj0 && xx < xj1) inter++;

    /* Point coincides with a vertex */

    if (xx == xj0 && yy == yj0)
    {
      inter = 1;
      continue;
    }
  }
  return ((inter % 2) != 0);
}

/****************************************************************************/
/*!
 **  Check if one point belongs to a vertical interval of a (limited) polyelem
 **
 ** \return  True if the point belongs to the polygon; False otherwise
 **
 ** \param[in]  zz   array of point coordinates of the point along Z or TEST
 **
 *****************************************************************************/
bool PolyElem::inside3D(double zz) const
{
  if (FFFF(zz)) return true;
  if (!FFFF(_zmin) && zz < _zmin) return false;
  if (!FFFF(_zmax) && zz > _zmax) return false;
  return true;
}

PolyElem PolyElem::reduceComplexity(double distmin) const
{
  int np = getNPoints();
  double dmin2 = distmin * distmin;
  PolyElem newpolyelem;

  /* Loop on the polygon vertices */

  double xcur = getX(0);
  double ycur = getY(0);
  newpolyelem.addPoint(xcur, ycur);

  int ecr = 1;
  while (ecr < np)
  {
    double xnext = getX(ecr);
    double ynext = getY(ecr);
    double dx = xnext - xcur;
    double dy = ynext - ycur;
    double dist2 = (dx * dx + dy * dy);
    if (dist2 >= dmin2)
    {
      // This point belongs to the new PolyElem
      newpolyelem.addPoint(xnext, ynext);
      xcur = xnext;
      ycur = ynext;
    }
    ecr++;
  }
  return newpolyelem;
}
