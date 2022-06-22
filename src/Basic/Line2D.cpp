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
#include "geoslib_f.h"
#include "geoslib_f_private.h"
#include "Basic/Line2D.hpp"

Line2D::Line2D(const VectorDouble& x,
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

Line2D::Line2D(const Line2D &m)
    : AStringable(m),
      ASerializable(m),
      _x(m._x),
      _y(m._y)
{

}

Line2D& Line2D::operator=(const Line2D &m)
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

Line2D::~Line2D()
{

}

Line2D* Line2D::create(const VectorDouble& x, const VectorDouble& y)
{
  return new Line2D(x, y);
}

String Line2D::toString(const AStringFormat* strfmt) const
{
  std::stringstream sstr;

  int npoints = getNPoints();
  sstr << "Number of Points = " << npoints << std::endl;

  AStringFormat sf;
  if (strfmt != nullptr) sf = *strfmt;
  if (sf.getLevel() > 0)
  {
    VectorDouble tab = VectorDouble(2 * npoints);
    for (int i = 0; i < npoints; i++)
    {
      tab[i] = _x[i];
      tab[i + npoints] = _y[i];
    }
    sstr << toMatrix("PolyLine Vertex Coordinates", VectorString(), VectorString(),
                     true, 2, npoints, tab);
  }
  return sstr.str();
}

int Line2D::dumpToNF(const String& neutralFilename, bool verbose) const
{
  std::ofstream os;
  int ret = 1;
  if (_fileOpenWrite(neutralFilename, "Line2D", os, verbose))
  {
    ret = _serialize(os, verbose);
    if (ret && verbose) messerr("Problem writing in the Neutral File.");
    os.close();
  }
  return ret;
}

Line2D* Line2D::createFromNF(const String& neutralFilename, bool verbose)
{
  Line2D* line2D = nullptr;
  std::ifstream is;
  if (_fileOpenRead(neutralFilename, "Line2D", is, verbose))
  {
    line2D = new Line2D();
    if (line2D->_deserialize(is, verbose))
    {
      if (verbose) messerr("Problem reading the Neutral File.");
      delete line2D;
      line2D = nullptr;
    }
    is.close();
  }
  return line2D;
}

int Line2D::_serialize(std::ostream& os, bool /*verbose*/) const
{
  if (getNPoints() <= 0) return 0;
  bool ret = true;
  ret = ret && _recordWriteVec<double>(os, "X-Coordinates of Line2D", _x);
  ret = ret && _recordWriteVec<double>(os, "Y-Coordinates of Line2D", _y);
  return ret ? 0 : 1;
}

int Line2D::_deserialize(std::istream& is, bool /*verbose*/)
{
  bool ret = true;
  ret = ret && _recordReadVec<double>(is, "X-Coordinates of Line2D", _x);
  ret = ret && _recordReadVec<double>(is, "Y-Coordinates of Line2D", _y);
  if (! ret) return 1;
  return 0;
}

void Line2D::init(const VectorDouble& x, const VectorDouble& y)
{
  int nvert = static_cast<int> (x.size());

  /* Load the new PolyLine */

  _x.resize(nvert,0);
  _y.resize(nvert,0);

  /* Copy the arrays */

  for (int i=0; i<nvert; i++)
  {
    _x[i] = x[i];
    _y[i] = y[i];
  }
}
