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
#include "Fractures/FracFault.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/Utilities.hpp"

#include <math.h>

FracFault::FracFault(double coord, double orient)
  : AStringable(),
    ASerializable(),
    _coord(coord),
    _orient(orient),
    _thetal(),
    _thetar(),
    _rangel(),
    _ranger()
{
}

FracFault::FracFault(const FracFault& r)
    : AStringable(r),
      ASerializable(r),
      _coord(r._coord),
      _orient(r._orient),
      _thetal(r._thetal),
      _thetar(r._thetar),
      _rangel(r._rangel),
      _ranger(r._ranger)
{
}

FracFault& FracFault::operator=(const FracFault& r)
{
  if (this != &r)
  {
    AStringable::operator=(r);
    ASerializable::operator=(r);
    _coord = r._coord;
    _orient = r._orient;
    _thetal = r._thetal;
    _thetar = r._thetar;
    _rangel = r._rangel;
    _ranger = r._ranger;
  }
  return *this;
}

FracFault::~FracFault()
{
}

String FracFault::toString(const AStringFormat* /*strfmt*/) const
{
  std::stringstream sstr;

  sstr << "Location of the Fault           = " << _coord << std::endl;
  sstr << "Fault orientation               = " << _orient << " (degree)" << std::endl;

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
double FracFault::faultAbscissae(double cote) const
{
  return (_coord + cote * tan(ut_deg2rad(_orient)));
}

void FracFault::addFaultPerFamily(double thetal,
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

bool FracFault::_deserialize(std::istream& is, bool /*verbose*/)
{
  bool ret = true;
  int nfam = getNFamilies();
  ret = ret && _recordRead<double>(is, "Abscissa of the first Fault point", _coord);
  ret = ret && _recordRead<double>(is, "Fault orientation", _orient);
  ret = ret && _recordRead<int>   (is, "Number of Families", nfam);
  ret = ret && _recordReadVec<double>(is, "Maximum Density on the left", _thetal, nfam);
  ret = ret && _recordReadVec<double>(is, "Maximum Density on the right", _thetar, nfam);
  ret = ret && _recordReadVec<double>(is, "Decrease Range on the left", _rangel, nfam);
  ret = ret && _recordReadVec<double>(is, "Decrease Range on the right", _ranger, nfam);
  return ret;
}

bool FracFault::_serialize(std::ostream& os, bool /*verbose*/) const
{
  bool ret = true;
  ret = ret && _recordWrite<double>(os, "Abscissa of the first Fault point", _coord);
  ret = ret && _recordWrite<double>(os, "Fault orientation", _orient);
  ret = ret && _recordWrite<int>   (os, "Number of Families", getNFamilies());
  ret = ret && _recordWriteVec<double>(os, "Maximum Density on the left", _thetal);
  ret = ret && _recordWriteVec<double>(os, "Maximum Density on the right", _thetar);
  ret = ret && _recordWriteVec<double>(os, "Decrease Range on the left", _rangel);
  ret = ret && _recordWriteVec<double>(os, "Decrease Range on the right", _ranger);
  return ret;
}
