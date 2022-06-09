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
#include "Faults/PolyLine.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/ASerializable.hpp"
#include "Basic/Utilities.hpp"
#include "geoslib_f.h"

PolyLine::PolyLine()
  : AStringable(),
    ASerializable(),
    _x(0)
  , _y(0)
{
}

PolyLine::PolyLine(const VectorDouble& x, const VectorDouble& y)
    : AStringable(),
      ASerializable(),
      _x(0),
      _y(0)
{
  init(x, y);
}

PolyLine::PolyLine(const PolyLine& r)
    : AStringable(r),
      ASerializable(r),
      _x(r._x),
      _y(r._y)
{
}

PolyLine& PolyLine::operator=(const PolyLine& r)
{
  if (this != &r)
  {
    AStringable::operator=(r);
    ASerializable::operator=(r);
    _x = r._x;
    _y = r._y;
  }
  return *this;
}

PolyLine::~PolyLine()
{
}

PolyLine* PolyLine::create(const VectorDouble& x, const VectorDouble& y)
{
  return new PolyLine(x, y);
}

void PolyLine::init(const VectorDouble& x, const VectorDouble& y)
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

String PolyLine::toString(const AStringFormat* strfmt) const
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
    sstr << toMatrix("PolyLine Vertex Coordinates", VectorString(), VectorString(),
                     true, 2, nvert, tab);
  }
  return sstr.str();
}

int PolyLine::_serialize(std::ostream& os, bool /*verbose*/) const
{
  if (getNVertices() <= 0) return 0;
  bool ret = true;
  ret = ret && _recordWriteVec<double>(os, "X-Coordinates of PolyLine", _x);
  ret = ret && _recordWriteVec<double>(os, "Y-Coordinates of PolyLine", _y);
  return ret ? 0 : 1;
}

int PolyLine::_deserialize(std::istream& is, bool /*verbose*/)
{
  bool ret = true;
  ret = ret && _recordReadVec<double>(is, "X-Coordinates of PolyLine", _x);
  ret = ret && _recordReadVec<double>(is, "Y-Coordinates of PolyLine", _y);
  if (! ret) return 1;
  return 0;
}

int PolyLine::dumpToNF(const String& neutralFilename, bool verbose) const
{
  std::ofstream os;
  int ret = 1;
  if (_fileOpenWrite(neutralFilename, "PolyLine", os, verbose))
  {
    ret = _serialize(os, verbose);
    if (ret && verbose) messerr("Problem writing in the Neutral File.");
    os.close();
  }
  return ret;
}

PolyLine* PolyLine::createFromNF(const String& neutralFilename, bool verbose)
{
  PolyLine* polyLine = nullptr;
  std::ifstream is;
  if (_fileOpenRead(neutralFilename, "PolyLine", is, verbose))
  {
    polyLine = new PolyLine();
    if (polyLine->_deserialize(is, verbose))
    {
      if (verbose) messerr("Problem reading the Neutral File.");
      delete polyLine;
      polyLine = nullptr;
    }
    is.close();
  }
  return polyLine;
}

