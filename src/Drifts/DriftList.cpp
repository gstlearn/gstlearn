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
#include <Drifts/ADrift.hpp>
#include "Drifts/DriftList.hpp"
#include "Space/ASpace.hpp"
#include "Space/SpaceRN.hpp"
#include "Basic/AException.hpp"
#include "Basic/Utilities.hpp"
#include "Drifts/DriftFactory.hpp"
#include "Drifts/DriftF.hpp"
#include "Drifts/DriftM.hpp"
#include "Db/Db.hpp"

DriftList::DriftList(const CovContext &ctxt)
    : AStringable(),
      _flagLinked(false),
      _driftCL(),
      _drifts(),
      _filtered(),
      _ctxt(ctxt)
{
}

DriftList::DriftList(const DriftList &r)
    : AStringable(r),
      _flagLinked(r._flagLinked),
      _driftCL(r._driftCL),
      _drifts(),
      _filtered(r._filtered),
      _ctxt(r._ctxt)
{
  for (auto e: r._drifts)
  {
    _drifts.push_back(dynamic_cast<ADrift*>(e->clone()));
  }
}

DriftList& DriftList::operator=(const DriftList &r)
{
  if (this != &r)
  {
    AStringable::operator=(r);
    _flagLinked = r._flagLinked;
    _driftCL  = r._driftCL;
    for (auto e: r._drifts)
    {
      _drifts.push_back(dynamic_cast<ADrift*>(e->clone()));
    }
    _filtered = r._filtered;
    _ctxt = r._ctxt;
  }
  return *this;
}

DriftList::~DriftList()
{
  delAllDrifts();
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
  return sstr.str();
}

void DriftList::addDrift(const ADrift* drift)
{
  if (drift == nullptr) return;
  _drifts.push_back(dynamic_cast<ADrift*>(drift->clone()));
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
  _driftCL.clear();
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

const ADrift* DriftList::getDrift(int il) const
{
  if (! _isDriftIndexValid(il)) return nullptr;
  return _drifts[il];
}

ADrift* DriftList::getDrift(int il)
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

  _driftCL.resize(nvar * nfeq * nbfl);
  if (_flagLinked)
  {
    for (int ivar = 0; ivar < nvar; ivar++)
      for (int ib = 0; ib < nfeq; ib++)
        for (int il = 0; il < nbfl; il++)
        {
          setDriftCL(ivar, il, ib, (ib == il));
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
            setDriftCL(ivar, il, ib, (ivar == jvar && il == jl));
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

VectorDouble DriftList::getDriftCLByPart(int ivar, int ib) const
{
  int nbfl = getDriftNumber();
  VectorDouble coef(nbfl,0.);
  for (int il = 0; il < nbfl; il++)
    coef[il] = getDriftCL(ivar, il, ib);
  return coef;
}

void DriftList::setDriftCLByPart(int ivar, int ib, const VectorDouble& coef)
{
  int number = getDriftNumber();
  if (number != (int) coef.size())
  {
    messerr("The dimension of 'vec' (%d) is not equal to the number of drift functions (%d)",
            (int) coef.size(), number);
    return;
  }
  for (int il = 0; il < number; il++)
    setDriftCL(ivar, il, ib, coef[il]);
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

double DriftList::evalDriftCoef(const Db *db, int iech, const VectorDouble &coeffs) const
{
  int ndrift = getDriftNumber();
  int ncoeff = (int) coeffs.size();
  if (ndrift != ncoeff)
  {
    messerr("Dimension of 'coeffs' (%d) should match number of drift functions (%d)",
        ncoeff, ndrift);
    return TEST;
  }

  double value = 0.;
  for (int ib = 0; ib < ndrift; ib++)
  {
    double drift = getDrift(db, ib, iech);
    if (FFFF(drift)) return TEST;
    value += coeffs[ib] * drift;
  }
  return value;
}

/**
 * Evaluate the Linear combination of drift terms at each sample of a Db
 * @param db     Target Db
 * @param coeffs Vector of coefficients (must have dimension of number of drift elements)
 * @param useSel True if the selection must be taken into account
 * @return
 */
VectorDouble DriftList::evalDriftCoefVec(const Db *db,
                                         const VectorDouble &coeffs,
                                         bool useSel) const
{
  VectorDouble vec;
  int ndrift = getDriftNumber();
  int ncoeff = (int) coeffs.size();
  if (ncoeff != ndrift)
  {
    messerr("'coeffs' dimension (%d) should match number of drift functions (%d)",
        ncoeff, ndrift);
    return vec;
  }

  for (int iech=0, nech=db->getSampleNumber(); iech<nech; iech++)
  {
    if (useSel && ! db->isActive(iech)) continue;
    double value = evalDriftCoef(db, iech, coeffs);
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
    const ADrift* drft = _drifts[il];
    int order = drft->getOrderIRF();
    if (order > max_order) max_order = order;
  }
  return (max_order);
}

/**
 * Check if a given drift type is defined among the drift functions
 * @param powers Vector of exponents for monomials
 * @param rank_fex Rank of the variable for external dift
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

bool DriftList::hasExternalDrift() const
{
  for (int il = 0; il < getDriftNumber(); il++)
  {
    if (getDrift(il)->isDriftExternal()) return true;
  }
  return false;
}

/****************************************************************************/
/*!
 **  Establish the drift rectangular matrix for a given Db
 **
 ** \return Returned matrix as a VD (Dimension = nech * nvar * nfeq * nvar)
 **
 ** \param[in]  db     Db structure
 ** \param[in]  member Member of the Kriging System (ECalcMember)
 **
 *****************************************************************************/
MatrixRectangular DriftList::evalDriftMat(const Db *db, const ECalcMember &member)
{
  int nvar = getNVariables();
  int nbfl = getDriftNumber();
  int nfeq = getDriftEquationNumber();
  int nech = db->getSampleNumber(true);
  int nrows = nech * nvar;
  int ncols = (isFlagLinked()) ? nfeq : nvar * nbfl;
  MatrixRectangular drfmat(nrows, ncols);

  /* Loop on the variables */

  int irow = 0;
  for (int ivar = 0; ivar < nvar; ivar++)
  {

    /* Loop on the samples */

    for (int iech = 0; iech < nech; iech++)
    {
      if (!db->isActive(iech)) continue;
      VectorDouble drftab = evalDriftVec(db, iech, member);

      /* Loop on the drift functions */

      int icol = 0;
      if (isFlagLinked())
      {
        for (int ib = 0; ib < nfeq; ib++)
        {
          drfmat.setValue(irow, icol, evalDriftValue(ivar, ib, drftab));
          icol++;
        }
      }
      else
      {
        for (int jvar = 0; jvar < nvar; jvar++)
          for (int jl = 0; jl < nbfl; jl++)
          {
            int jb = jvar + nvar * jl;
            drfmat.setValue(irow, icol, evalDriftValue(ivar, jb, drftab));
            icol++;
          }
      }
      irow++;
    }
  }
  return 0;
}

VectorDouble DriftList::evalDriftVec(const Db *db,
                                     int iech,
                                     const ECalcMember &member) const
{
  int ndrift = getDriftNumber();
  VectorDouble drftab(ndrift);
  for (int il = 0; il < ndrift; il++)
     drftab[il] = evalDrift(db, iech, il, member);
  return drftab;
}

void DriftList::evalDriftVecInPlace(const Db *db,
                                    int iech,
                                    const ECalcMember &member,
                                    VectorDouble &drftab) const
{
  int ndrift = getDriftNumber();
  for (int il = 0; il < ndrift; il++)
     drftab[il] = evalDrift(db, iech, il, member);
}

double DriftList::evalDriftValue(int ivar, int ib, const VectorDouble &drftab) const
{
  double value = 0.;
  for (int il = 0, nbfl = getDriftNumber(); il < nbfl; il++)
    value += drftab[il] * getDriftCL(ivar, il, ib);
  return value;
}

/**
 * Evaluate a given drift function for a given sample
 * @param db     Db structure
 * @param iech   Rank of the target sample
 * @param il     Rank of the drift function
 * @param member Member type (used to check filtering)
 * @return
 */
double DriftList::evalDrift(const Db *db,
                            int iech,
                            int il,
                            const ECalcMember &member) const
{
  if (member != ECalcMember::LHS && isFiltered(il))
    return 0.;
  else
  {
    if (! _isDriftIndexValid(il)) return TEST;
    return _drifts[il]->eval(db, iech);
  }
  return TEST;
}
