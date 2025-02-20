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
#include "Basic/VectorNumT.hpp"
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
      _ctxt(ctxt),
      _mean(VectorDouble(ctxt.getNVar(), 0.))
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
      _ctxt(r._ctxt),
      _mean(r._mean)
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
    _mean = r._mean;
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

void DriftList::copyCovContext(const CovContext& ctxt)
{
  _ctxt = ctxt;
  _update();
}

void DriftList::_update()
{
  if ((int)_mean.size() != _ctxt.getNVar())
  _mean =  VectorDouble(_ctxt.getNVar(), 0.);
}

String DriftList::toString(const AStringFormat* /*strfmt*/) const
{
  std::stringstream sstr;
  if (getNDrift() <= 0)
    sstr << toVector("Known Mean(s)", getMeans());
  // TODO: could be added but changes all non-regression files
  //    sstr << "(Note: Simple Kriging will be used)" << std::endl;
  for (int i = 0, nbfl = getNDrift(); i < nbfl; i++)
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

void DriftList::delDrift(unsigned int rank)
{
  if (_drifts.empty()) return;
  if (! _isDriftIndexValid(rank)) return;
  _drifts.erase(_drifts.begin() + rank);
  _filtered.erase(_filtered.begin() + rank);
  _betaHat.erase(_betaHat.begin() + rank);
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
  _mean.resize(0);

}

bool DriftList::isDriftFiltered(int i) const
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

int DriftList::getNDriftEquation() const
{
  int nbfl = getNDrift();
  int nvar = getNVar();
  int ndriftEquationNumber = (_flagLinked) ? nbfl : nbfl * nvar;
  return ndriftEquationNumber;
}

/**
 * Check that the set of drift functions is valid
 * @return
 */
bool DriftList::isValid() const
{
  int nbfl = getNDrift();

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
  return checkArg("Drift Rank", i, getNDrift());
}

/**
 * Returns if the rank of the Drift Equation valid or not
 *
 * @param ib Rank of the drift equation
 */
bool DriftList::_isDriftEquationValid(int ib) const
{
  return checkArg("Drift Equation", ib, getNDriftEquation());
}

void DriftList::resetDriftList()
{
  int nvar = getNVar();
  int nfeq = getNDriftEquation();
  int nbfl = getNDrift();

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
  int nbfl = getNDrift();
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
  int nbfl = getNDrift();
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

double DriftList::computeDrift(const Db* db, int ib, int iech) const
{
  if (! _isDriftIndexValid(ib)) return TEST;
  return _drifts[ib]->eval(db, iech);
}

VectorVectorDouble DriftList::getDrifts(const Db* db, bool useSel) const
{
  VectorVectorDouble vecvec;
  int nbfl = getNDrift();
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
  int nbfl = getNDrift();
  int ncoeff = (int) coeffs.size();
  if (nbfl != ncoeff)
  {
    messerr("Dimension of 'coeffs' (%d) should match number of drift functions (%d)", ncoeff, nbfl);
    return TEST;
  }

  double value = 0.;
  for (int ib = 0; ib < nbfl; ib++)
  {
    double drift = computeDrift(db, ib, iech);
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
  int nbfl = getNDrift();
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
int DriftList::getDriftMaxIRFOrder() const
{
  int max_order = 0;
  for (int il = 0, nbfl = getNDrift(); il < nbfl; il++)
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
  for (int il = 0, nbfl = getNDrift(); il < nbfl; il++)
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
  for (int il = 0, nbfl = getNDrift(); il < nbfl; il++)
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
  for (int il = 0, nbfl = getNDrift(); il < nbfl; il++)
  {
    if (getDrift(il)->isDriftExternal()) return true;
  }
  return false;
}

int DriftList::getNExtDrift() const
{
  int nfex = 0;
  for (int il = 0; il < getNDrift(); il++)
  {
    if (getDrift(il)->isDriftExternal()) nfex++;
  }
  return nfex;
}

VectorInt DriftList::_getActiveVariables(int ivar0) const
{
  int nvar = getNVar();

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
                                          const ECalcMember& member) const
{
  MatrixRectangular mat;
  VectorInt ivars = _getActiveVariables(ivar0);
  if (ivars.empty()) return mat;

  // Create the sets of Vector of valid sample indices per variable (not masked and defined)
  VectorVectorInt index = db->getSampleRanks(ivars, nbgh, true, true, true);

  int error = evalDriftMatByRanks(mat, db, index, ivar0, member);
  return (error == 0) ? mat : MatrixRectangular();
}

/**
 * @brief Calculate the Drift matrix
 *
 * @param mat Drift matrix (possibly resized)
 * @param db Data Db
 * @param sampleRanks Vector of sample ranks in 'db'
 * @param ivar0 Rank of the variable (-1 for all)
 * @param member CalcMember
 *
 * @return int Error returned code
 */
int DriftList::evalDriftMatByRanks(MatrixRectangular& mat,
                                   const Db* db,
                                   const VectorVectorInt& sampleRanks,
                                   int ivar0,
                                   const ECalcMember& member) const
{
  VectorInt ivars = _getActiveVariables(ivar0);
  if (ivars.empty()) return 1;

  // Creating the matrix
  int neq = VH::count(sampleRanks);
  if (neq <= 0)
  {
    messerr("The returned matrix has no valid sample and no valid variable");
    return 1;
  }

  int nvar  = getNVar();
  int nbfl  = getNDrift();
  int nfeq  = getNDriftEquation();
  int ncols = (isFlagLinked()) ? nfeq : nvar * nbfl;
  if (ncols <= 0) return 0;
  mat.resize(neq, ncols);
  mat.fill(0.);

  for (int ivar = 0, irow = 0, nvars = (int)ivars.size(); ivar < nvars; ivar++)
  {
    int ivar1 = ivars[ivar];

    /* Loop on the samples */

    int nechs = (int)sampleRanks[ivar1].size();
    for (int jech = 0; jech < nechs; jech++, irow++)
    {
      int iech = sampleRanks[ivar][jech];
      {
        for (int ib = 0; ib < nfeq; ib++)
        {
          double value = evalDriftValue(db, iech, ivar, ib, member);
          mat.setValue(irow, ib, value);
        }
      }
    }
  }
  return 0;
}

/****************************************************************************/
/*!
 **  Establish the drift rectangular matrix for a given Db
 **
 ** \return Returned matrix
 ** (Dimension/ nrows = nvar * nech; ncols = nfeq * nvar)
 **
 ** \param[in]  mat     Drift matrix (possibly resized)
 ** \param[in]  db     Db structure
 ** \param[in]  iech2  Index of active samples in db
 ** \param[in]  krigopt KrigOpt structure
 **
 *****************************************************************************/
int DriftList::evalDriftMatByTarget(MatrixRectangular& mat,
                                    const Db* db,
                                    int iech2,
                                    const KrigOpt& krigopt) const
{
  VectorInt ivars = VH::sequence(getNVar());
  if (ivars.empty()) return 1;

  // Create the sets of Vector of valid sample indices per variable
  // (not masked and defined)
  VectorVectorInt index = db->getSampleRanks(ivars, {iech2}, true, false, false);

  // Creating the matrix
  int neq = VH::count(index);
  if (neq <= 0)
  {
    messerr("The returned matrix has no valid sample and no valid variable");
    return 1;
  }

  int nvar  = getNVar();
  int nbfl  = getNDrift();
  int nfeq  = getNDriftEquation();
  int ncols = (isFlagLinked()) ? nfeq : nvar * nbfl;
  if (ncols <= 0) return 0;
  mat.resize(neq, ncols);
  mat.fill(0.);

  for (int ivar = 0; ivar < nvar; ivar++)
    for (int ib = 0; ib < nfeq; ib++)
    {
      double value = evalDriftValue(db, iech2, ivar, ib, ECalcMember::RHS);
      if (FFFF(value)) return 1;
      mat.setValue(ivar, ib, value);
    }

  // In case of combined R.H.S., modify the output matrix
  if (krigopt.isMatLC()) mat = mat.compressMatLC(*krigopt.getMatLC(), true);
  return 0;
}

VectorDouble
DriftList::evalDriftBySample(const Db* db, int iech, const ECalcMember& member) const
{
  int nbfl = getNDrift();
  VectorDouble drftab(nbfl);
  evalDriftBySampleInPlace(db, iech, member, drftab);
  return drftab;
}

void DriftList::evalDriftBySampleInPlace(const Db* db,
                                         int iech,
                                         const ECalcMember& member,
                                         VectorDouble& drftab) const
{
  int nbfl = getNDrift();
  if (nbfl != (int)drftab.size()) drftab.resize(nbfl);
  for (int il = 0; il < nbfl; il++)
  {
    if (member != ECalcMember::LHS && isDriftFiltered(il))
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
double DriftList::evalDriftValue(
  const Db* db, int iech, int ivar, int ib, const ECalcMember& member) const
{
  int nbfl     = getNDrift();
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
    if (!_flagLinked) il = ib - ivar * nbfl;
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
double
DriftList::evalDrift(const Db* db, int iech, int il, const ECalcMember& member) const
{
  if (member != ECalcMember::LHS && isDriftFiltered(il)) return 0.;
  if (!_isDriftIndexValid(il)) return TEST;
  return _drifts[il]->eval(db, iech);
  return TEST;
}

const DriftList* DriftList::createReduce(const VectorInt& validVars) const
{
  int ecr = 0;
  int lec = 0;
  VectorBool valids(getNVar(), false);
  int nvar = (int)validVars.size();
  VectorDouble mean(nvar, 0);

  for (int ivar = 0; ivar < nvar; ivar++) valids[validVars[ivar]] = true;
  for (int ivar = 0; ivar < getNVar(); ivar++)
  {
    if (valids[ivar])
    {
      mean[ecr++] = _mean[lec];
    }
    lec++;
  }
  DriftList* driftlist = new DriftList(_ctxt);
  driftlist->setMeans(mean);
  return driftlist;
}

void DriftList::setMeans(const VectorDouble& mean)
{
  if (_mean.size() == mean.size()) _mean = mean;
}

double DriftList::getMean(int ivar) const
{
  if (ivar < 0 || ivar >= (int)_mean.size()) my_throw("Invalid argument in _getMean");
  return _mean[ivar];
}

/**
 * Define the Mean for one variable
 * @param mean Value for the mean
 * @param ivar Rank of the variable (starting from 0)
 */
void DriftList::setMean(const double mean, int ivar)
{
  if (ivar < 0 || ivar >= (int)_mean.size()) my_throw("Invalid argument in _setMean");
  _mean[ivar] = mean;
}

/****************************************************************************/
/*!
 **  Evaluate the drift with a given sample and a given variable
 **  The value is scaled by 'coeffs'
 **
 ** \param[in]  db      Db structure
 ** \param[in]  iech    Rank of the sample
 ** \param[in]  ivar    Rank of the variable
 ** \param[in]  coeffs  Vector of coefficients
 **
 *****************************************************************************/
double DriftList::evalDriftVarCoef(const Db* db,
                                   int iech,
                                   int ivar,
                                   const VectorDouble& coeffs) const
{

  double drift = 0.;
  for (int ib = 0, nfeq = getNDriftEquation(); ib < nfeq; ib++)
    drift += evalDriftValue(db, iech, ivar, ib, ECalcMember::LHS) * coeffs[ib];
  return drift;
}

/**
 * A vector of the drift evaluation (for all samples)
 * @param db     Db structure
 * @param coeffs Vector of drift coefficients
 * @param useSel When TRUE, only non masked samples are returned
 * @return The vector of values
 *
 * @remark When no drift is defined, a vector is returned filled with the variable
 * mean
 */

VectorDouble
DriftList::evalDriftVarCoefs(const Db* db, const VectorDouble& coeffs, bool useSel) const
{
  VectorDouble vec;
  vec = evalDriftCoefs(db, coeffs, useSel);
  return vec;
}