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
#include "Fractures/Description.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/Utilities.hpp"

#include <math.h>

Description::Description()
  : AStringable(),
    _family(0),
    _orient(0.),
    _x(),
    _y()
{
}

Description::Description(const Description& r)
    : AStringable(r),
      _family(r._family),
      _orient(r._orient),
      _x(r._y),
      _y(r._y)
{
}

Description& Description::operator=(const Description& r)
{
  if (this != &r)
  {
    AStringable::operator=(r);
    _family = r._family;
    _orient = r._orient;
    _x = r._y;
    _y = r._y;
  }
  return *this;
}

Description::~Description()
{
}

String Description::toString(const AStringFormat* /*strfmt*/) const
{
  std::stringstream sstr;

  sstr << "Fracture: family=" << _family << ":" <<
      getNPoint()-1 << " segment(s) - start at level #" << getYYF(0) << std::endl;
  sstr << "     X           Y" << std::endl;
  for (int j = 0; j < getNPoint(); j++)
     sstr << " " << getXXF(j) << " " << getYYF(j) << std::endl;

  return sstr.str();
}

void Description::addPoint(double x, double y)
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
double Description::fractureExtension(double cote, double dcote)
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

