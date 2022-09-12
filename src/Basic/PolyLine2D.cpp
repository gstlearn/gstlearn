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
#include "geoslib_f_private.h"

#include "Basic/PolyLine2D.hpp"

PolyLine2D::PolyLine2D(const VectorDouble& x,
               const VectorDouble& y)
    : AStringable(),
      ASerializable(),
      _x(x),
      _y(y)
{
  // Check that both vectors have the same dimension
  if (x.size() != y.size())
  {
    _x.clear();
    _y.clear();
  }
}

PolyLine2D::PolyLine2D(const PolyLine2D &m)
    : AStringable(m),
      ASerializable(m),
      _x(m._x),
      _y(m._y)
{

}

PolyLine2D& PolyLine2D::operator=(const PolyLine2D &m)
{
  if (this != &m)
  {
    AStringable::operator=(m);
    ASerializable::operator=(m);
    _x = m._x;
    _y = m._y;
  }
  return *this;
}

PolyLine2D::~PolyLine2D()
{

}

PolyLine2D* PolyLine2D::create(const VectorDouble& x, const VectorDouble& y)
{
  return new PolyLine2D(x, y);
}

String PolyLine2D::toString(const AStringFormat* /*strfmt*/) const
{
  std::stringstream sstr;

  int npoints = getNPoints();

  VectorDouble tab = VectorDouble(2 * npoints);
  for (int i = 0; i < npoints; i++)
  {
    tab[i] = _x[i];
    tab[i + npoints] = _y[i];
  }
  sstr << toMatrix("Line Vertex Coordinates", VectorString(), VectorString(),
                   true, 2, npoints, tab);
  return sstr.str();
}

PolyLine2D* PolyLine2D::createFromNF(const String& neutralFilename, bool verbose)
{
  PolyLine2D* line2D = nullptr;
  std::ifstream is;
  line2D = new PolyLine2D();
  bool success = false;
  if (line2D->_fileOpenRead(neutralFilename, is, verbose))
  {
    success =  line2D->deserialize(is, verbose);
  }
  if (! success)
  {
    delete line2D;
    line2D = nullptr;
  }
  return line2D;
}

/**
 * Serialization (by Point rather than by Coordinate)
 * This is maintained for all classes using this interface for serialization
 * @param os Output Stream
 * @param verbose Verbose flag
 * @return
 */
bool PolyLine2D::_serialize(std::ostream& os, bool /*verbose*/) const
{
  if (getNPoints() <= 0) return false;
  bool ret = true;
  ret = ret && _recordWrite<int>(os, "Number of Points", (int) _x.size());

  VectorDouble buffer(2);
  for (int i = 0; i < (int) _x.size(); i++)
  {
    buffer[0] = _x[i];
    buffer[1] = _y[i];
    ret = ret && _recordWriteVec<double>(os, "", buffer);
  }
  return ret;
}

/**
 * Deserialization (by sample)
 * @param is Input stream
 * @param verbose Verbose flag
 * @return
 */
bool PolyLine2D::_deserialize(std::istream& is, bool /*verbose*/)
{
  int np = 0;
  bool ret = true;
  VectorDouble buffer(2);
  ret = ret && _recordRead<int>(is, "Number of Points", np);
  if (np < 0)
  {
    messerr("Something wrong in your file: the Number of Points read is equal to %d",np);
    return false;
  }
  _x.resize(np);
  _y.resize(np);
  for (int i = 0; i < np; i++)
  {
    ret = ret && _recordReadVec<double>(is, "", buffer, 2);
    _x[i] = buffer[0];
    _y[i] = buffer[1];
  }
  return ret;
}

void PolyLine2D::init(const VectorDouble& x, const VectorDouble& y)
{
  int nvert = static_cast<int> (x.size());

  /* Load the new Line */

  _x.resize(nvert,0);
  _y.resize(nvert,0);

  /* Copy the arrays */

  for (int i=0; i<nvert; i++)
  {
    _x[i] = x[i];
    _y[i] = y[i];
  }
}

void PolyLine2D::addPoint(double x, double y)
{
  int n = getNPoints();
  _x.resize(n+1);
  _y.resize(n+1);
  _x[n] = x;
  _y[n] = y;
}
