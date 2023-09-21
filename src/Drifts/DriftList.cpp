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
#include "Drifts/DriftList.hpp"
#include "Space/ASpace.hpp"
#include "Space/SpaceRN.hpp"
#include "Basic/AException.hpp"
#include "Basic/Utilities.hpp"
#include "Drifts/ADriftElem.hpp"
#include "Drifts/DriftFactory.hpp"
#include "Drifts/DriftF.hpp"
#include "Drifts/DriftM.hpp"
#include "Db/Db.hpp"

DriftList::DriftList(const ASpace* space)
    : ADrift(space),
      _flagLinked(false),
      _driftCoef(),
      _drifts(),
      _filtered()
{
}

DriftList::DriftList(const DriftList &r)
    : ADrift(r),
      _flagLinked(r._flagLinked),
      _driftCoef(r._driftCoef),
      _drifts(),
      _filtered(r._filtered)
{
  for (auto e: r._drifts)
  {
    _drifts.push_back(dynamic_cast<ADriftElem*>(e->clone()));
  }
}

DriftList& DriftList::operator=(const DriftList &r)
{
  if (this != &r)
  {
    ADrift::operator=(r);
    _flagLinked = r._flagLinked;
    _driftCoef  = r._driftCoef;
    for (auto e: r._drifts)
    {
      _drifts.push_back(dynamic_cast<ADriftElem*>(e->clone()));
    }
    _filtered = r._filtered;
  }
  return *this;
}

DriftList::~DriftList()
{
  delAllDrifts();
}

bool DriftList::isConsistent(const ASpace* /*space*/) const
{
  return true;
}

String DriftList::toString(const AStringFormat* /*strfmt*/) const
{
  std::stringstream sstr;
  for (int i = 0; i < (int) getDriftNumber(); i++)
  {
    sstr << _drifts[i]->toString();
    if (_filtered[i])
      sstr << " (This component is filtered)";
    sstr << std::endl;
  }

//   Display the coefficients
  for (int ivar = 0; ivar < getNVariables(); ivar++)
    for (int ib = 0; ib < getDriftEquationNumber(); ib++)
    {
      sstr << "ivar = " << ivar << " ib = " << ib << " : ";
      for (int il = 0; il < getDriftNumber(); il++)
         sstr << " " << getDriftCoef(ivar, il, ib);
      sstr << std::endl;
    }
  return sstr.str();
}

double DriftList::eval(const Db* /*db*/, int /*iech1*/) const
{
  double drift = 0;
  return drift;
}

void DriftList::addDrift(const ADriftElem* drift)
{
  if (drift == nullptr) return;
  _drifts.push_back(dynamic_cast<ADriftElem*>(drift->clone()));
  _filtered.push_back(false);
  updateDriftList();
}

void DriftList::delDrift(unsigned int i)
{
  if (_drifts.empty()) return;
  if (! _isDriftIndexValid(i)) return;
  _drifts.erase(_drifts.begin() + i);
  _filtered.erase(_filtered.begin() + i);
  updateDriftList();
}

void DriftList::delAllDrifts()
{
  if (! _drifts.empty())
    for (auto e: _drifts)
    {
      delete e;
    }
  _drifts.clear();
  _filtered.clear();
  _driftCoef.clear();
}

bool DriftList::isFiltered(int i) const
{
  if (! _isDriftIndexValid(i)) return false;
  return _filtered[i];
}

void DriftList::setFiltered(int i, bool filter)
{
  if (! _isDriftIndexValid(i)) return;
  _filtered[i] = filter;
}

const ADriftElem* DriftList::getDrift(int il) const
{
  if (! _isDriftIndexValid(il)) return nullptr;
  return _drifts[il];
}

ADriftElem* DriftList::getDrift(int il)
{
  if (! _isDriftIndexValid(il)) return nullptr;
  return _drifts[il];
}

int DriftList::getRankFex(int il) const
{
  if (! _isDriftIndexValid(il)) return 0;
  return _drifts[il]->getRankFex();
}

String DriftList::getDriftName(int il) const
{
  if (! _isDriftIndexValid(il)) return String();
  return _drifts[il]->getDriftName();
}

int DriftList::getDriftEquationNumber() const
{
  int ndrift = getDriftNumber();
  int nvar = getNVariables();
  int ndriftEquationNumber = (_flagLinked) ? ndrift : ndrift * nvar;
  return ndriftEquationNumber;
}

/**
 * Check that the set of drift functions is valid
 * @return
 */
bool DriftList::isValid() const
{
  int ndrift = getDriftNumber();

  // Check that the same drift function has not been called twice
  for (int il=0; il<ndrift; il++)
  {
    const String str_il = _drifts[il]->getDriftName();

    for (int jl=0; jl<il; jl++)
    {
      const String str_jl = _drifts[jl]->getDriftName();

      if (str_il.compare(str_jl) == 0)
      {
        messerr("Set of drift functions is invalid: %d and %d are similar",il+1,jl+1);
        return false;
      }
    }
  }
  return true;
}

int DriftList::getNVariables() const
{
  if (getDriftNumber() > 0)
    return _drifts[0]->getNVariables();
  return 0;
}

bool DriftList::_isDriftIndexValid(int i) const
{
  if (i < 0 || i >= getDriftNumber())
  {
    mesArg("Drift Rank",i,getDriftNumber());
    return false;
  }
  return true;
}

bool DriftList::_isDriftEquationValid(int ib) const
{
  if (ib < 0 || ib >= getDriftEquationNumber())
  {
    mesArg("Drift Equation",ib,getDriftEquationNumber());
    return false;
  }
  return true;
}

void DriftList::updateDriftList()
{
  int nvar = getNVariables();
  int nfeq = getDriftEquationNumber();
  int nbfl = getDriftNumber();

  /* Copy the coefficients from the old to the new structure */

  _driftCoef.resize(nvar * nfeq * nbfl);
  if (_flagLinked)
  {
    for (int ivar = 0; ivar < nvar; ivar++)
      for (int ib = 0; ib < nfeq; ib++)
        for (int il = 0; il < nbfl; il++)
        {
          setDriftCoef(ivar, il, ib, (ib == il));
        }
  }
  else
  {
    for (int ivar = 0; ivar < nvar; ivar++)
      for (int jvar = 0; jvar < nvar; jvar++)
        for (int jl = 0; jl < nbfl; jl++)
          for (int il = 0; il < nbfl; il++)
          {
            int ib = jvar + nvar * jl;
            setDriftCoef(ivar, il, ib, (ivar == jvar && il == jl));
          }
  }

  // Resize the 'filtered' array (if necessary)

  if (nbfl != (int) _filtered.size())
    _filtered.resize(nbfl, false);
}

VectorDouble DriftList::getDriftByColumn(const Db* db, int ib, bool useSel) const
{
  if (! _isDriftIndexValid(ib)) return VectorDouble();

  int nech = db->getSampleNumber(useSel);
  VectorDouble vec(nech);

  int ecr = 0;
  for (int iech=0; iech<db->getSampleNumber(); iech++)
  {
    if (useSel && ! db->isActive(iech)) continue;
    vec[ecr++] = _drifts[ib]->eval(db, iech);
  }
  return vec;
}

VectorDouble DriftList::getDriftBySample(const Db* db, int iech) const
{
  int ndrift = getDriftNumber();
  VectorDouble vec(ndrift);

  for (int ib=0; ib<ndrift; ib++)
  {
    vec[ib] = _drifts[ib]->eval(db, iech);
  }
  return vec;
}

VectorDouble DriftList::getDriftCoefByPart(int ivar, int ib) const
{
  int number = getDriftNumber();
  VectorDouble coef(number,0.);
  for (int il = 0; il < number; il++)
    coef[il] = getDriftCoef(ivar, il, ib);
  return coef;
}

void DriftList::setDriftCoefByPart(int ivar, int ib, const VectorDouble& coef)
{
  int number = getDriftNumber();
  if (number != (int) coef.size())
  {
    messerr("The dimension of 'vec' (%d) is not equal to the number of drift functions (%d)",
            (int) coef.size(), number);
    return;
  }
  for (int il = 0; il < number; il++)
    setDriftCoef(ivar, il, ib, coef[il]);
}

double DriftList::getDrift(const Db* db, int ib, int iech) const
{
  if (! _isDriftIndexValid(ib)) return TEST;
  return _drifts[ib]->eval(db, iech);
}

VectorVectorDouble DriftList::getDrifts(const Db* db, bool useSel) const
{
  VectorVectorDouble vec;
  int ndrift = getDriftNumber();
  for (int ib=0; ib<ndrift; ib++)
  {
    vec.push_back(getDriftByColumn(db, ib, useSel));
  }
  return vec;
}

/**
 * Evaluate the Linear combination of drift terms at each sample of a Db
 * @param db     Target Db
 * @param coeffs Vector of coefficients (must have dimension of number of drift elements)
 * @param useSel True if the selection must be taken into account
 * @return
 */
VectorDouble DriftList::evalDrifts(const Db* db,
                                   const VectorDouble& coeffs,
                                   bool useSel) const
{
  VectorDouble vec;
  int ndrift = getDriftNumber();
  int ncoeff = (int) coeffs.size();
  int nech   = db->getSampleNumber();
  if (ncoeff != ndrift)
  {
    messerr("'coeffs' dimension (%d) should match number of drift functions (%d)",
        ncoeff, ndrift);
    return vec;
  }

  for (int iech=0; iech<nech; iech++)
  {
    if (useSel && ! db->isActive(iech)) continue;
    double value = 0.;
    for (int ib=0; ib<ndrift; ib++)
    {
      double drift = getDrift(db, ib, iech);
      if (FFFF(drift))
      {
        value = TEST;
        break;
      }
      value += coeffs[ib] * getDrift(db, ib, iech);
    }

    vec.push_back(value);
  }
  return vec;
}

/**
 * @return Maximum IRF-order (-1 for order-2 stationarity)
 */
int DriftList::getDriftMaxIRFOrder(void) const
{
  int max_order = 0;
  for (int il = 0; il < getDriftNumber(); il++)
  {
    const ADriftElem* drft = _drifts[il];
    int order = drft->getOrderIRF();
    if (order > max_order) max_order = order;
  }
  return (max_order);
}

/**
 * Check if a given drift type is defined among the drift functions
 * @param type0 Target drift type (EDrift.hpp)
 * @return
 */
bool DriftList::isDriftDefined(const VectorInt &powers, int rank_fex) const
{
  for (int il = 0; il < getDriftNumber(); il++)
  {
    if (_drifts[il]->isDriftExternal())
    {
      if (_drifts[il]->getRankFex() == rank_fex) return true;
    }
    else
    {
      if (_drifts[il]->getPowers() == powers) return true;
    }
  }
  return false;
}

/**
 * Check if at least one drift function exists whose type is different
 * from the target type
 * @param powers Vector of exponent for monomials of a polynomial drift
 * @param rank_fex Rank of the variable for external Drift
 * @return
 */
bool DriftList::isDriftDifferentDefined(const VectorInt &powers, int rank_fex) const
{
  for (int il = 0; il < getDriftNumber(); il++)
  {
    if (_drifts[il]->isDriftExternal())
    {
      if (_drifts[il]->getRankFex() != rank_fex) return true;
    }
    else
    {
      if (_drifts[il]->getPowers() != powers) return true;
    }
  }
  return false;
}

void DriftList::copyCovContext(const CovContext& ctxt)
{
  int number = (int) _drifts.size();
  for (int i = 0; i < number; i++)
    _drifts[i]->copyCovContext(ctxt);
}

bool DriftList::hasExternalDrift() const
{
  for (int il = 0; il < getDriftNumber(); il++)
  {
    if (getDrift(il)->isDriftExternal()) return true;
  }
  return false;
}
