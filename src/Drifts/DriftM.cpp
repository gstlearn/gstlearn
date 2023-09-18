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
#include "Enum/EDrift.hpp"

#include "Drifts/DriftM.hpp"
#include "Drifts/ADriftElem.hpp"
#include "Db/Db.hpp"

DriftM::DriftM(const VectorInt &powers,
               double coeff0,
               const VectorDouble &coeffs,
               const CovContext &ctxt)
    : ADriftElem(EDrift::UC, 0, ctxt),
      _monomialPower(powers),
      _coeff0(coeff0),
      _monomialCoeffs(coeffs)
{
  if (! _monomialPower.empty() &&_monomialCoeffs.empty())
    _monomialCoeffs.resize((int) _monomialPower.size(), 1.);
}

DriftM::DriftM(const DriftM &r)
    : ADriftElem(r)
{
}

DriftM& DriftM::operator=(const DriftM &r)
{
  if (this != &r)
  {
    ADriftElem::operator =(r);
  }
  return *this;
}

DriftM::~DriftM()
{
}

double DriftM::eval(const Db* db, int iech) const
{
  double value = _coeff0;
  for (int idim = 0, ndim = _monomialPower.size(); idim < ndim; idim++)
  {
    double locoor = db->getCoordinate(iech,idim);
    double locpow = _monomialPower[idim];
    double coeff  = _monomialCoeffs[idim];
    value *= coeff * pow(locoor, locpow);
  }
  return value;
}

String DriftM::getDriftName() const
{
  std::stringstream sstr;
  if (_monomialPower.empty())
    return "Universality Condition";
  else
  {
    sstr << "Drift ";
    for (int idim = 0, ndim = _monomialPower.size(); idim < ndim; idim++)
    {
      double locpow = _monomialPower[idim];
      if (locpow <= 0) continue;
      double coeff  = _monomialCoeffs[idim];
      if (coeff != 1.)
        sstr << "(" << coeff << "*)";
      sstr << "x_" << idim+1 ;
      if (locpow <= 1) continue;
      sstr << "^" << locpow;
    }
  }
  return sstr.str();
}

int DriftM::getOrderIRF() const
{
  int irf = -1;
  for (int idim = 0, ndim = _monomialPower.size(); idim < ndim; idim++)
  {
    double locpow = _monomialPower[idim] - 1;
    if (locpow > irf) irf = locpow;
  }
  return irf;
}

int  DriftM::getNDim() const
{
  return (int) _monomialPower.size();
}
