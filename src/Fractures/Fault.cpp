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
#include "Fractures/Fault.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/Utilities.hpp"

#include <math.h>

Fault::Fault(double coord, double orient)
  : AStringable(),
    _coord(coord),
    _orient(orient),
    _thetal(),
    _thetar(),
    _rangel(),
    _ranger()
{
}

Fault::Fault(const Fault& r)
    : AStringable(r),
      _coord(r._coord),
      _orient(r._orient),
      _thetal(r._thetal),
      _thetar(r._thetar),
      _rangel(r._rangel),
      _ranger(r._ranger)
{
}

Fault& Fault::operator=(const Fault& r)
{
  if (this != &r)
  {
    AStringable::operator=(r);
    _coord = r._coord;
    _orient = r._orient;
    _thetal = r._thetal;
    _thetar = r._thetar;
    _rangel = r._rangel;
    _ranger = r._ranger;
  }
  return *this;
}

Fault::~Fault()
{
}

String Fault::toString(const AStringFormat* /*strfmt*/) const
{
  std::stringstream sstr;
  sstr << "Location of the Fault           = " << _coord << std::endl;
  sstr << "Fault orientation               = " << _orient << "(deg)" << std::endl;

  int number = (int) _thetal.size();
  for (int j = 0; j < number; j++)
  {
    sstr << toTitle(2, "Family #%d/%d", j + 1, number);
    sstr << "Intensity maximum value (left)  = " << _thetal[j] << std::endl;
    sstr << "Intensity range (left)          = " << _rangel[j] << std::endl;
    sstr << "Intensity maximum value (right) = " << _thetar[j] << std::endl;
    sstr << "Intensity range (right)         = " << _ranger[j] << std::endl;
  }
  return sstr.str();
}

/****************************************************************************/
/*!
 **  Calculate the abscissae of a fault at a given elevation
 **
 ** \return The fault abscissae
 **
 ** \param[in]  cote   Ordinate of the fracture starting point
 **
 *****************************************************************************/
double Fault::faultAbscissae(double cote) const
{
  return (_coord + cote * tan(ut_deg2rad(_orient)));
}

void Fault::addFaultPerFamily(double thetal,
                              double thetar,
                              double rangel,
                              double ranger)
{
  int nfam  = getNFamilies();
  _thetal.resize(nfam + 1);
  _thetar.resize(nfam + 1);
  _rangel.resize(nfam + 1);
  _ranger.resize(nfam + 1);

  _thetal[nfam] = thetal;
  _thetar[nfam] = thetar;
  _rangel[nfam] = rangel;
  _ranger[nfam] = ranger;
}
