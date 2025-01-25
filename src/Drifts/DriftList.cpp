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
#include "Drifts/ADrift.hpp"
#include "Drifts/DriftList.hpp"
#include "Basic/Utilities.hpp"
#include "Basic/VectorHelper.hpp"
#include "Drifts/DriftFactory.hpp"
#include "Db/Db.hpp"

DriftList::DriftList(const CovContext &ctxt)
    : AStringable(),
      _flagLinked(false),
      _flagCombined(false),
      _driftCL(),
      _drifts(),
      _betaHat(),
      _filtered(),
      _ctxt(ctxt)
{
}

DriftList::DriftList(const DriftList &r)
    : AStringable(r),
      _flagLinked(r._flagLinked),
      _flagCombined(r._flagCombined),
      _driftCL(r._driftCL),
      _drifts(),
      _betaHat(r._betaHat),
      _filtered(r._filtered),
      _ctxt(r._ctxt)
{
  for (const auto& e: r._drifts)
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
    _flagCombined = r._flagCombined;
    _driftCL  = r._driftCL;
    for (const auto& e: r._drifts)
    {
      _drifts.push_back(dynamic_cast<ADrift*>(e->clone()));
    }
    _betaHat = r._betaHat;
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
  for (int i = 0, nbfl = getDriftNumber(); i < nbfl; i++)
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
  _betaHat.push_back(0.);
  resetDriftList();
}

void DriftList::delDrift(unsigned int i)
{
  if (_drifts.empty()) return;
  if (! _isDriftIndexValid(i)) return;
  _drifts.erase(_drifts.begin() + i);
  _filtered.erase(_filtered.begin() + i);
  _betaHat.erase(_betaHat.begin() + i);
  resetDriftList();
}

void DriftList::delAllDrifts()
{
  if (! _drifts.empty())
    for (const auto& e: _drifts)
    {
      delete e;
    }
  _drifts.clear();
  _filtered.clear();
  _betaHat.clear();
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
  int nbfl = getDriftNumber();
  int nvar = getNVariables();
  int ndriftEquationNumber = (_flagLinked) ? nbfl : nbfl * nvar;
  return ndriftEquationNumber;
}

/**
 * Check that the set of drift functions is valid
 * @return
 */
bool DriftList::isValid() const
{
  int nbfl = getDriftNumber();

  // Check that the same drift function has not been called twice
  for (int il=0; il<nbfl; il++)
  {
    const String str_il = _drifts[il]->getDriftName();

    for (int jl=0; jl<il; jl++)
    {
      const String str_jl = _drifts[jl]->getDriftName();

      if (str_il == str_jl)
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
  return checkArg("Drift Rank", i, getDriftNumber());
}

/**
 * Returns if the rank of the Drift Equation valid or not
 *
 * @param ib Rank of the drift equation
 */
bool DriftList::_isDriftEquationValid(int ib) const
{
  return checkArg("Drift Equation", ib, getDriftEquationNumber());
}

void DriftList::resetDriftList()
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
          _setDriftCL(ivar, il, ib, (ib == il));
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
            _setDriftCL(ivar, il, ib, (ivar == jvar && il == jl));
          }
  }

  // Resize the 'filtered' array (if necessary)

  if (nbfl != (int) _filtered.size())
    _filtered.resize(nbfl, false);
}

/**
 * Set the values of the internal array 'driftCL'
 * This feature is used when the drift equation does not coincide with a drift function
 * (when working with the gradient of the target variable, for example)
 *
 * @param ivar Rank of the variable (_nVar)
 * @param ib Rank of the drift equation (_driftEquationNumber)
 * @param coef Vector of coefficients
 */
void DriftList::setDriftCLByPart(int ivar, int ib, const VectorDouble& coef)
{
  int nbfl = getDriftNumber();
  if (nbfl != (int) coef.size())
  {
    messerr("The dimension of 'vec' (%d) is not equal to the number of drift functions (%d)",
            (int) coef.size(), nbfl);
    return;
  }
  for (int il = 0; il < nbfl; il++)
    _setDriftCL(ivar, il, ib, coef[il]);

  _flagCombined = true;
}

/**
 * Check that the drift 'ib' is defined for at least one variable one sample
 *
 * @param db Data file used for reading the drift contents (depends on its type)
 * @param ib Rank of the drift function
 * @param nech Number of samples to be checked
 * @param nbgh Vector of sample indices within the data base
 * @param loctype Locator to be checked
 */
bool DriftList::isDriftSampleDefined(const Db *db,
                                     int ib,
                                     int nech,
                                     const VectorInt &nbgh,
                                     const ELoc &loctype) const
{
  int nbfl = getDriftNumber();
  int nvar = db->getNLoc(loctype);

  if (_flagCombined)
  {
    for (int il = 0; il < nbfl; il++)
      for (int ivar = 0; ivar < nvar; ivar++)
      {
        if (isZero(_getDriftCL(ivar, il, ib))) continue;
        for (int iech = 0; iech < nech; iech++)
          if (!FFFF(db->getLocVariable(loctype, nbgh[iech], ivar))) return true;
      }
  }
  else
  {
    for (int ivar = 0; ivar < nvar; ivar++)
      for (int iech = 0; iech < nech; iech++)
        if (!FFFF(db->getLocVariable(loctype, nbgh[iech], ivar))) return true;
  }
  return false;
}

double DriftList::getDrift(const Db* db, int ib, int iech) const
{
  if (! _isDriftIndexValid(ib)) return TEST;
  return _drifts[ib]->eval(db, iech);
}

VectorVectorDouble DriftList::getDrifts(const Db* db, bool useSel) const
{
  VectorVectorDouble vecvec;
  int nbfl = getDriftNumber();
  int nech = db->getNSample(useSel);
  VectorDouble vec(nech);

  for (int ib=0; ib<nbfl; ib++)
  {
    int ecr = 0;
    for (int iech=0; iech<db->getNSample(); iech++)
    {
      if (useSel && ! db->isActive(iech)) continue;
      vec[ecr++] = _drifts[ib]->eval(db, iech);
    }
    vecvec.push_back(vec);
  }
  return vecvec;
}

double DriftList::evalDriftCoef(const Db *db, int iech, const VectorDouble &coeffs) const
{
  int nbfl = getDriftNumber();
  int ncoeff = (int) coeffs.size();
  if (nbfl != ncoeff)
  {
    messerr("Dimension of 'coeffs' (%d) should match number of drift functions (%d)", ncoeff, nbfl);
    return TEST;
  }

  double value = 0.;
  for (int ib = 0; ib < nbfl; ib++)
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
VectorDouble DriftList::evalDriftCoefs(const Db *db,
                                       const VectorDouble &coeffs,
                                       bool useSel) const
{
  VectorDouble vec;
  int nbfl = getDriftNumber();
  int ncoeff = (int) coeffs.size();
  if (ncoeff != nbfl)
  {
    messerr("'coeffs' dimension (%d) should match number of drift functions (%d)",
        ncoeff, nbfl);
    return vec;
  }

  for (int iech=0, nech=db->getNSample(); iech<nech; iech++)
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
  for (int il = 0, nbfl = getDriftNumber(); il < nbfl; il++)
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
  for (int il = 0, nbfl = getDriftNumber(); il < nbfl; il++)
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
  for (int il = 0, nbfl = getDriftNumber(); il < nbfl; il++)
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
  for (int il = 0, nbfl = getDriftNumber(); il < nbfl; il++)
  {
    if (getDrift(il)->isDriftExternal()) return true;
  }
  return false;
}

int DriftList::getExternalDriftNumber() const
{
  int nfex = 0;
  for (int il = 0; il < getDriftNumber(); il++)
  {
    if (getDrift(il)->isDriftExternal()) nfex++;
  }
  return nfex;
}

VectorInt DriftList::_getActiveVariables(int ivar0) const
{
  int nvar = getNVariables();

  VectorInt ivars;
  if (ivar0 >= 0)
  {
    if (!checkArg("Argument 'ivar0'", ivar0, nvar)) return VectorInt();
    ivars.push_back(ivar0);
  }
  else
  {
    ivars = VH::sequence(nvar);
  }
  return ivars;
}

/****************************************************************************/
/*!
 **  Establish the drift rectangular matrix for a given Db
 **
 ** \return Returned matrix (Dimension/ nrows = nvar * nech; ncols = nfeq * nvar)
 **
 ** \param[in]  db     Db structure
 ** \param[in]  ivar0  Rank of the variable (-1 for all variables)
 ** \param[in]  nbgh   Vector of indices of active samples in db (optional)
 ** \param[in]  member Member of the Kriging System (ECalcMember)
 **
 *****************************************************************************/
MatrixRectangular DriftList::evalDriftMat(const Db* db,
                                          int ivar0,
                                          const VectorInt& nbgh,
                                          const ECalcMember& member)
{
  MatrixRectangular drfmat;
  VectorInt ivars = _getActiveVariables(ivar0);
  if (ivars.empty()) return drfmat;

  // Create the sets of Vector of valid sample indices per variable (not masked and defined)
  VectorVectorInt index = db->getSampleRanks(ivars, nbgh, true, true, true);

  return evalDriftMatByRanks(db, index, ivar0, member);
}

MatrixRectangular DriftList::evalDriftMatByRanks(const Db* db,
                                                 const VectorVectorInt& sampleRanks,
                                                 int ivar0,
                                                 const ECalcMember& member)
{
  MatrixRectangular drfmat;
  VectorInt ivars = _getActiveVariables(ivar0);
  if (ivars.empty()) return drfmat;

  int nvar  = getNVariables();
  int nbfl  = getDriftNumber();
  int nfeq  = getDriftEquationNumber();
  int ncols = (isFlagLinked()) ? nfeq : nvar * nbfl;

  // Creating the matrix
  int neq = VH::count(sampleRanks);
  if (neq <= 0)
  {
    messerr("The returned matrix does not have any valid sample for any valid variable");
    return drfmat;
  }
  drfmat.resize(neq, ncols);

  /* Loop on the variables */

  int irow = 0;
  for (int ivar = 0, nvars = (int)ivars.size(); ivar < nvars; ivar++)
  {
    int ivar1 = ivars[ivar];

    /* Loop on the samples */

    int nechs = (int)sampleRanks[ivar].size();
    for (int jech = 0; jech < nechs; jech++)
    {
      int iech = sampleRanks[ivar][jech];

      /* Loop on the drift functions */

      if (isFlagLinked())
      {
        VectorDouble drftab = evalDriftBySample(db, iech, member);
        for (int ib = 0; ib < nfeq; ib++)
        {
          drfmat.setValue(irow, ib, drftab[ib]); // TODO: to be generalized
        }
      }
      else
      {
        int icol = 0;
        for (int jvar = 0; jvar < nvar; jvar++)
          for (int jl = 0; jl < nbfl; jl++)
          {
            int jb = jl + jvar * nbfl;
            drfmat.setValue(irow, icol, evalDriftValue(db, iech, ivar1, jb, member));
            icol++;
          }
      }
      irow++;
    }
  }
  return drfmat;
}

/****************************************************************************/
/*!
 **  Establish the drift rectangular matrix for a given Db
 **
 ** \return Returned matrix
 ** (Dimension/ nrows = nvar * nech; ncols = nfeq * nvar)
 **
 ** \param[in]  db     Db structure
 ** \param[in]  ivar0  Rank of the variable (-1 for all variables)
 ** \param[in]  iech2  Index of active samples in db
 ** \param[in]  member Member of the Kriging System (ECalcMember)
 **
 *****************************************************************************/
MatrixRectangular DriftList::evalDriftMatByTarget(const Db* db,
                                                   int ivar0,
                                                   int iech2,
                                                   const ECalcMember& member)
{
  MatrixRectangular drfmat;
  int nvar        = getNVariables();
  int nbfl        = getDriftNumber();
  int nfeq        = getDriftEquationNumber();
  int ncols       = (isFlagLinked()) ? nfeq : nvar * nbfl;
  VectorInt ivars = _getActiveVariables(ivar0);
  if (ivars.empty()) return drfmat;

  // Create the sets of Vector of valid sample indices per variable
  // (not masked and defined)
  VectorVectorInt index =
    db->getSampleRanks(ivars, {iech2}, true, false, false);

  // Creating the matrix
  int neq = VH::count(index);
  if (neq <= 0)
  {
    messerr("The returned matrix does not have any valid sample for any valid "
            "variable");
    return drfmat;
  }
  drfmat.resize(neq, ncols);

  /* Loop on the variables */

  int irow = 0;
  for (int ivar = 0, nvars = (int)ivars.size(); ivar < nvars; ivar++)
  {
    int ivar1 = ivars[ivar];

    /* Loop on the samples */

    int nechs = (int)index[ivar].size();
    for (int jech = 0; jech < nechs; jech++)
    {
      int iech = index[ivar][jech];

      /* Loop on the drift functions */

      if (isFlagLinked())
      {
        VectorDouble drftab = evalDriftBySample(db, iech, member);
        for (int ib = 0; ib < nfeq; ib++)
        {
          drfmat.setValue(irow, ib, drftab[ib]); // TODO: to be generalized
        }
      }
      else
      {
        int icol = 0;
        for (int jvar = 0; jvar < nvar; jvar++)
          for (int jl = 0; jl < nbfl; jl++)
          {
            int jb = jl + jvar * nbfl;
            drfmat.setValue(irow, icol,
                            evalDriftValue(db, iech, ivar1, jb, member));
            icol++;
          }
      }
      irow++;
    }
  }
  return drfmat;
}

VectorDouble DriftList::evalDriftBySample(const Db *db,
                                          int iech,
                                          const ECalcMember &member) const
{
  int nbfl = getDriftNumber();
  VectorDouble drftab(nbfl);
  evalDriftBySampleInPlace(db, iech, member, drftab);
  return drftab;
}

void DriftList::evalDriftBySampleInPlace(const Db *db,
                                         int iech,
                                         const ECalcMember &member,
                                         VectorDouble &drftab) const
{
  int nbfl = getDriftNumber();
  if (nbfl != (int) drftab.size()) drftab.resize(nbfl);
  for (int il = 0; il < nbfl; il++)
  {
    if (member != ECalcMember::LHS && isFiltered(il))
      drftab[il] = 0.;
    else
      drftab[il] = _drifts[il]->eval(db, iech);
  }
}

/**
 * Returns the Drift value for a given variable and a given drift function
 *
 * @param db Db structure
 * @param iech Rank of the sample
 * @param ivar Rank of the variable
 * @param ib   Rank of the Drift function
 * @param member ECalcMember characteristics
 */
double DriftList::evalDriftValue(const Db *db,
                                 int iech,
                                 int ivar,
                                 int ib,
                                 const ECalcMember &member) const
{
  int nbfl = getDriftNumber();
  double value = 0.;
  if (_flagCombined)
  {
    for (int il = 0; il < nbfl; il++)
    {
      double local = evalDrift(db, iech, il, member);
      if (FFFF(local)) return TEST;
      value += local * _getDriftCL(ivar, il, ib);
    }
  }
  else
  {
    int il = ib;
    if (! _flagLinked) il = ib - ivar * nbfl;
    if (il < 0 || il >= nbfl) return 0.;
    value = evalDrift(db, iech, il, member);
  }
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
  if (member != ECalcMember::LHS && isFiltered(il)) return 0.;
  if (!_isDriftIndexValid(il)) return TEST;
  return _drifts[il]->eval(db, iech);
  return TEST;
}
