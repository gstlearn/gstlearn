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
               const CovContext &ctxt)
    : ADriftElem(EDrift::MONO, ctxt),
      _monomialPower(powers)
{
}

DriftM::DriftM(const DriftM &r)
    : ADriftElem(r),
      _monomialPower(r._monomialPower)
{
}

DriftM& DriftM::operator=(const DriftM &r)
{
  if (this != &r)
  {
    ADriftElem::operator =(r);
    _monomialPower = r._monomialPower;
  }
  return *this;
}

DriftM::~DriftM()
{
}

double DriftM::eval(const Db* db, int iech) const
{
  double value = 1.;
  for (int idim = 0, ndim = _monomialPower.size(); idim < ndim; idim++)
  {
    double locoor = db->getCoordinate(iech,idim);
    double locpow = _monomialPower[idim];
    value *= pow(locoor, locpow);
  }
  return value;
}

String DriftM::getDriftName() const
{
  std::stringstream sstr;
  if (_monomialPower.empty())
    sstr << "Universality_Condition";
  else
  {
    sstr << "Drift:";
    bool flag_first = true;
    for (int idim = 0, ndim = _monomialPower.size(); idim < ndim; idim++)
    {
      double locpow = _monomialPower[idim];
      if (locpow > 0)
      {
        if (!flag_first) sstr << "*";
        sstr << "x" << idim+1;
        if (locpow > 1) sstr << "^" << locpow;
        flag_first = false;
      }
    }
  }
  return sstr.str();
}

int DriftM::getOrderIRF() const
{
  int irf = -1;
  for (int idim = 0, ndim = _monomialPower.size(); idim < ndim; idim++)
  {
    double locpow = _monomialPower[idim];
    if (locpow > irf) irf = locpow;
  }
  return irf;
}

int  DriftM::getOrderIRFIdim(int idim) const
{
  if (idim < getNDim()) return -1;
  return _monomialPower[idim];
}

int  DriftM::getNDim() const
{
  return (int) _monomialPower.size();
}
