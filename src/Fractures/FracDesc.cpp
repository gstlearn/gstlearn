/******************************************************************************/
/* gstlearn C++ Library                                                       */
/*                                                                            */
/* Authors: <authors>                                                         */
/*                                                                            */
/* License: BSD 3 Clause                                                      */
/******************************************************************************/
#include "Fractures/FracDesc.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/Utilities.hpp"

#include <math.h>

FracDesc::FracDesc()
  : AStringable(),
    _family(0),
    _orient(0.),
    _x(),
    _y()
{
}

FracDesc::FracDesc(const FracDesc& r)
    : AStringable(r),
      _family(r._family),
      _orient(r._orient),
      _x(r._x),
      _y(r._y)
{
}

FracDesc& FracDesc::operator=(const FracDesc& r)
{
  if (this != &r)
  {
    AStringable::operator=(r);
    _family = r._family;
    _orient = r._orient;
    _x = r._x;
    _y = r._y;
  }
  return *this;
}

FracDesc::~FracDesc()
{
}

String FracDesc::toString(const AStringFormat* strfmt) const
{
  std::stringstream sstr;
  if (getNPoint() <= 0) return sstr.str();

  sstr << "Fracture: family=" << _family+1 << " : " <<
      getNPoint()-1 << " segment(s) starting at level #" << getYYF(0) << std::endl;

  AStringFormat sf;
  if (strfmt != nullptr) sf = *strfmt;
  if (sf.getLevel() > 1)
  {
    sstr << "     X           Y" << std::endl;
    for (int j = 0; j < getNPoint(); j++)
      sstr << " " << getXXF(j) << " " << getYYF(j) << std::endl;
  }
  return sstr.str();
}

void FracDesc::addPoint(double x, double y)
{
  int np = getNPoint();
  _x.resize(np+1);
  _y.resize(np+1);
  _x[np] = x;
  _y[np] = y;
}

/****************************************************************************/
/*!
 **  Calculate the fracture extension
 **
 ** \return  The fracture extension
 **
 ** \param[in]  cote         Selected layer or TEST for all layers
 ** \param[in]  dcote        Tolerance on the layer elevation
 **
 *****************************************************************************/
double FracDesc::fractureExtension(double cote, double dcote)
{
  double dist = 0.;
  for (int i = 0; i < getNPoint() - 1; i++)
  {
    double distx = getXXF(i+1) - getXXF(i);
    double disty = getYYF(i+1) - getYYF(i);
    if (!FFFF(cote) && (getYYF(i) < cote - dcote || getYYF(i+1) < cote - dcote))
      continue;
    dist += sqrt(distx * distx + disty * disty);
  }
  return (dist);
}

