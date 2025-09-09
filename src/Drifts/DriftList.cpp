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
#include "Basic/Utilities.hpp"
#include "Basic/VectorHelper.hpp"
#include "Basic/VectorNumT.hpp"
#include "Db/Db.hpp"
#include "Drifts/ADrift.hpp"
#include "Drifts/DriftFactory.hpp"

namespace gstlrn
{
DriftList::DriftList(const CovContext& ctxt)
  : AStringable()
  , _flagLinked(false)
  , _flagCombined(false)
  , _driftCL()
  , _drifts()
  , _betaHat()
  , _filtered()
  , _ctxt(ctxt)
  , _mean(VectorDouble(ctxt.getNVar(), 0.))
{
}

DriftList::DriftList(const DriftList& r)
  : AStringable(r)
  , _flagLinked(r._flagLinked)
  , _flagCombined(r._flagCombined)
  , _driftCL(r._driftCL)
  , _drifts()
  , _betaHat(r._betaHat)
  , _filtered(r._filtered)
  , _ctxt(r._ctxt)
  , _mean(r._mean)
{
  for (const auto& e: r._drifts)
  {
    _drifts.push_back(dynamic_cast<ADrift*>(e->clone()));
  }
}

DriftList& DriftList::operator=(const DriftList& r)
{
  if (this != &r)
  {
    AStringable::operator=(r);
    _flagLinked   = r._flagLinked;
    _flagCombined = r._flagCombined;
    _driftCL      = r._driftCL;
    _mean         = r._mean;
    for (const auto& e: r._drifts)
    {
      _drifts.push_back(dynamic_cast<ADrift*>(e->clone()));
    }
    _betaHat  = r._betaHat;
    _filtered = r._filtered;
    _ctxt     = r._ctxt;
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
  if (static_cast<Id>(_mean.size()) != _ctxt.getNVar())
    _mean = VectorDouble(_ctxt.getNVar(), 0.);
}

String DriftList::toString(const AStringFormat* /*strfmt*/) const
{
  std::stringstream sstr;
  if (getNDrift() <= 0)
    sstr << toVector("Known Mean(s)", getMeans());
  // TODO: could be added but changes all non-regression files
  //    sstr << "(Note: Simple Kriging will be used)" << std::endl;
  for (Id i = 0, nbfl = getNDrift(); i < nbfl; i++)
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

void DriftList::delDrift(size_t rank)
{
  if (_drifts.empty()) return;
  if (!_isDriftIndexValid(static_cast<Id>(rank))) return;
  _drifts.erase(_drifts.begin() + rank);
  _filtered.erase(_filtered.begin() + rank);
  _betaHat.erase(_betaHat.begin() + rank);
  resetDriftList();
}

void DriftList::delAllDrifts()
{
  if (!_drifts.empty())
    for (const auto& e: _drifts)
    {
      delete e;
    }
  _drifts.clear();
  _filtered.clear();
  _betaHat.clear();
  _driftCL.clear();
  _mean.resize(0);
  _update();
}

bool DriftList::isDriftFiltered(Id i) const
{
  if (!_isDriftIndexValid(i)) return false;
  return _filtered[i];
}

void DriftList::setFiltered(Id i, bool filter)
{
  if (!_isDriftIndexValid(i)) return;
  _filtered[i] = filter;
}

const ADrift* DriftList::getDrift(Id il) const
{
  if (!_isDriftIndexValid(il)) return nullptr;
  return _drifts[il];
}

Id DriftList::getRankFex(Id il) const
{
  if (!_isDriftIndexValid(il)) return 0;
  return _drifts[il]->getRankFex();
}

String DriftList::getDriftName(Id il) const
{
  if (!_isDriftIndexValid(il)) return String();
  return _drifts[il]->getDriftName();
}

Id DriftList::getNDriftEquation() const
{
  auto nbfl               = getNDrift();
  auto nvar               = getNVar();
  Id ndriftEquationNumber = (_flagLinked) ? nbfl : nbfl * nvar;
  return ndriftEquationNumber;
}

/**
 * Check that the set of drift functions is valid
 * @return
 */
bool DriftList::isValid() const
{
  auto nbfl = getNDrift();

  // Check that the same drift function has not been called twice
  for (Id il = 0; il < nbfl; il++)
  {
    const String str_il = _drifts[il]->getDriftName();

    for (Id jl = 0; jl < il; jl++)
    {
      const String str_jl = _drifts[jl]->getDriftName();

      if (str_il == str_jl)
      {
        messerr("Set of drift functions is invalid: %d and %d are similar", il + 1, jl + 1);
        return false;
      }
    }
  }
  return true;
}

bool DriftList::_isDriftIndexValid(Id i) const
{
  return checkArg("Drift Rank", i, getNDrift());
}

/**
 * Returns if the rank of the Drift Equation valid or not
 *
 * @param ib Rank of the drift equation
 */
bool DriftList::_isDriftEquationValid(Id ib) const
{
  return checkArg("Drift Equation", ib, getNDriftEquation());
}

void DriftList::resetDriftList()
{
  auto nvar = getNVar();
  auto nfeq = getNDriftEquation();
  auto nbfl = getNDrift();

  /* Copy the coefficients from the old to the new structure */

  _driftCL.resize(nvar * nfeq * nbfl);
  if (_flagLinked)
  {
    for (Id ivar = 0; ivar < nvar; ivar++)
      for (Id ib = 0; ib < nfeq; ib++)
        for (Id il = 0; il < nbfl; il++)
        {
          _setDriftCL(ivar, il, ib, (ib == il));
        }
  }
  else
  {
    for (Id ivar = 0; ivar < nvar; ivar++)
      for (Id jvar = 0; jvar < nvar; jvar++)
        for (Id jl = 0; jl < nbfl; jl++)
          for (Id il = 0; il < nbfl; il++)
          {
            Id ib = jvar + nvar * jl;
            _setDriftCL(ivar, il, ib, (ivar == jvar && il == jl));
          }
  }

  // Resize the 'filtered' array (if necessary)

  if (nbfl != static_cast<Id>(_filtered.size()))
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
void DriftList::setDriftCLByPart(Id ivar, Id ib, const VectorDouble& coef)
{
  auto nbfl = getNDrift();
  if (nbfl != static_cast<Id>(coef.size()))
  {
    messerr("The dimension of 'vec' (%d) is not equal to the number of drift functions (%d)",
            static_cast<Id>(coef.size()), nbfl);
    return;
  }
  for (Id il = 0; il < nbfl; il++)
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
bool DriftList::isDriftSampleDefined(const Db* db,
                                     Id ib,
                                     Id nech,
                                     const VectorInt& nbgh,
                                     const ELoc& loctype) const
{
  auto nbfl = getNDrift();
  Id nvar   = db->getNLoc(loctype);

  if (_flagCombined)
  {
    for (Id il = 0; il < nbfl; il++)
      for (Id ivar = 0; ivar < nvar; ivar++)
      {
        if (isZero(_getDriftCL(ivar, il, ib))) continue;
        for (Id iech = 0; iech < nech; iech++)
          if (!FFFF(db->getLocVariable(loctype, nbgh[iech], ivar))) return true;
      }
  }
  else
  {
    for (Id ivar = 0; ivar < nvar; ivar++)
      for (Id iech = 0; iech < nech; iech++)
        if (!FFFF(db->getLocVariable(loctype, nbgh[iech], ivar))) return true;
  }
  return false;
}

double DriftList::computeDrift(const Db* db, Id ib, Id iech) const
{
  if (!_isDriftIndexValid(ib)) return TEST;
  return _drifts[ib]->eval(db, iech);
}

VectorVectorDouble DriftList::getDrifts(const Db* db, bool useSel) const
{
  VectorVectorDouble vecvec;
  auto nbfl = getNDrift();
  Id nech   = db->getNSample(useSel);
  VectorDouble vec(nech);

  for (Id ib = 0; ib < nbfl; ib++)
  {
    Id ecr = 0;
    for (Id iech = 0; iech < db->getNSample(); iech++)
    {
      if (useSel && !db->isActive(iech)) continue;
      vec[ecr++] = _drifts[ib]->eval(db, iech);
    }
    vecvec.push_back(vec);
  }
  return vecvec;
}

double DriftList::evalDriftCoef(const Db* db, Id iech, const VectorDouble& coeffs) const
{
  auto nbfl = getNDrift();
  Id ncoeff = static_cast<Id>(coeffs.size());
  if (nbfl != ncoeff)
  {
    messerr("Dimension of 'coeffs' (%d) should match number of drift functions (%d)", ncoeff, nbfl);
    return TEST;
  }

  double value = 0.;
  for (Id ib = 0; ib < nbfl; ib++)
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
VectorDouble DriftList::evalDriftCoefs(const Db* db,
                                       const VectorDouble& coeffs,
                                       bool useSel) const
{
  VectorDouble vec;
  auto nbfl = getNDrift();
  Id ncoeff = static_cast<Id>(coeffs.size());
  if (ncoeff != nbfl)
  {
    messerr("'coeffs' dimension (%d) should match number of drift functions (%d)",
            ncoeff, nbfl);
    return vec;
  }

  for (Id iech = 0, nech = db->getNSample(); iech < nech; iech++)
  {
    if (useSel && !db->isActive(iech)) continue;
    double value = evalDriftCoef(db, iech, coeffs);
    vec.push_back(value);
  }
  return vec;
}

/**
 * @return Maximum IRF-order (-1 for order-2 stationarity)
 */
Id DriftList::getDriftMaxIRFOrder() const
{
  Id max_order = 0;
  for (Id il = 0, nbfl = getNDrift(); il < nbfl; il++)
  {
    const ADrift* drft = _drifts[il];
    Id order           = drft->getOrderIRF();
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
bool DriftList::isDriftDefined(const VectorInt& powers, Id rank_fex) const
{
  for (Id il = 0, nbfl = getNDrift(); il < nbfl; il++)
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
bool DriftList::isDriftDifferentDefined(const VectorInt& powers, Id rank_fex) const
{
  for (Id il = 0, nbfl = getNDrift(); il < nbfl; il++)
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
  for (Id il = 0, nbfl = getNDrift(); il < nbfl; il++)
  {
    if (getDrift(il)->isDriftExternal()) return true;
  }
  return false;
}

Id DriftList::getNExtDrift() const
{
  Id nfex = 0;
  for (Id il = 0; il < getNDrift(); il++)
  {
    if (getDrift(il)->isDriftExternal()) nfex++;
  }
  return nfex;
}

VectorInt DriftList::_getActiveVariables(Id ivar0) const
{
  auto nvar = getNVar();

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
 ** \param[in]  nbgh   Vector of indices of active samples in db (optional)
 ** \param[in]  member Member of the Kriging System (ECalcMember)
 **
 *****************************************************************************/
MatrixDense DriftList::evalDriftMat(const Db* db,
                                    const VectorInt& nbgh,
                                    const ECalcMember& member) const
{
  MatrixDense mat;
  VectorInt ivars = _getActiveVariables(-1);
  if (ivars.empty()) return mat;

  // Create the sets of Vector of valid sample indices per variable (not masked and defined)
  VectorVectorInt index = db->getSampleRanks(ivars, nbgh, true, true, true);

  Id error = evalDriftMatByRanksInPlace(mat, db, index, member);
  return (error == 0) ? mat : MatrixDense();
}

Id DriftList::evalDriftMatInPlace(MatrixDense& mat,
                                  const Db* db,
                                  const VectorInt& nbgh,
                                  const ECalcMember& member) const
{
  VectorInt ivars = _getActiveVariables(-1);
  if (ivars.empty()) return 1;
  VectorVectorInt index = db->getSampleRanks(ivars, nbgh, true, true, true);

  return evalDriftMatByRanksInPlace(mat, db, index, member);
}

/**
 * @brief Calculate the Drift matrix
 *
 * @param mat Drift matrix (possibly resized)
 * @param db Data Db
 * @param sampleRanks Vector of sample ranks in 'db'
 * @param member CalcMember
 *
 * @return Id Error returned code
 */
Id DriftList::evalDriftMatByRanksInPlace(MatrixDense& mat,
                                         const Db* db,
                                         const VectorVectorInt& sampleRanks,
                                         const ECalcMember& member) const
{
  VectorVectorInt sampleRanksLoc = sampleRanks;
  if (sampleRanksLoc.empty())
  {
    sampleRanksLoc = db->getSampleRanks();
  }

  // Creating the matrix
  Id neq = VH::count(sampleRanksLoc);
  if (neq <= 0)
  {
    messerr("The returned matrix has no valid sample and no valid variable");
    return 1;
  }

  auto nvar = getNVar();
  auto nbfl = getNDrift();
  auto nfeq = (isFlagLinked()) ? nbfl : nvar * nbfl;
  if (nfeq <= 0) return 0;
  mat.resize(neq, nfeq);
  mat.fill(0.);

  for (Id ivar = 0, irow = 0; ivar < nvar; ivar++)
  {

    /* Loop on the samples */

    Id nechs = static_cast<Id>(sampleRanksLoc[ivar].size());
    for (Id jech = 0; jech < nechs; jech++, irow++)
    {
      Id iech = sampleRanksLoc[ivar][jech];
      {
        for (Id ib = 0; ib < nfeq; ib++)
        {
          double value = evalDriftValue(db, iech, ivar, ib, member);
          mat.setValue(irow, ib, value);
        }
      }
    }
  }
  return 0;
}

/**
 * @brief Returns the Matrix of the Drift elements
 *
 * @param db
 * @param sampleRanks
 * @param member
 * @return MatrixDense
 */
MatrixDense DriftList::evalDriftMatByRanks(const Db* db,
                                           const VectorVectorInt& sampleRanks,
                                           const ECalcMember& member) const
{
  MatrixDense mat;

  Id error = evalDriftMatByRanksInPlace(mat, db, sampleRanks, member);
  if (error) mat.resize(0, 0);
  return mat;
}

VectorDouble DriftList::evalMeanVecByRanks(const Db* db,
                                           const VectorVectorInt& sampleRanks) const
{
  VectorVectorInt sampleRanksLoc = sampleRanks;
  if (sampleRanksLoc.empty())
  {
    sampleRanksLoc = db->getSampleRanks();
  }

  // Creating the matrix
  Id neq = VH::count(sampleRanksLoc);
  if (neq <= 0)
  {
    messerr("The returned matrix has no valid sample and no valid variable");
    return 1;
  }
  auto nvar = getNVar();

  VectorDouble vec(neq, 0.);

  for (Id ivar = 0, ecr = 0, irow = 0; ivar < nvar; ivar++)
  {
    Id nechs = static_cast<Id>(sampleRanksLoc[ivar].size());
    for (Id jech = 0; jech < nechs; jech++, irow++, ecr++)
    {
      vec[ecr] = getMean(ivar);
    }
  }
  return vec;
}

thread_local VectorInt viech2(1);
thread_local VectorInt ivars;

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
Id DriftList::evalDriftMatByTargetInPlace(MatrixDense& mat,
                                          const Db* db,
                                          Id iech2,
                                          const KrigOpt& krigopt) const
{
  VH::sequenceInPlace(getNVar(), ivars);
  if (ivars.empty()) return 1;

  // Create the sets of Vector of valid sample indices per variable
  // (not masked and defined)
  viech2[0]                    = iech2;
  const VectorVectorInt& index = db->getSampleRanks(ivars, viech2, true, false, false);

  // Creating the matrix
  Id neq = VH::count(index);
  if (neq <= 0)
  {
    messerr("The returned matrix has no valid sample and no valid variable");
    return 1;
  }

  auto nvar = getNVar();
  auto nbfl = getNDrift();
  auto nfeq = getNDriftEquation();
  Id ncols  = (isFlagLinked()) ? nfeq : nvar * nbfl;
  if (ncols <= 0) return 0;
  mat.resize(neq, ncols);
  mat.fill(0.);

  for (Id ivar = 0; ivar < nvar; ivar++)
    for (Id ib = 0; ib < nfeq; ib++)
    {
      double value = evalDriftValue(db, iech2, ivar, ib, ECalcMember::RHS);
      if (FFFF(value)) return 1;
      mat.setValue(ivar, ib, value);
    }

  // In case of combined R.H.S., modify the output matrix
  if (krigopt.hasMatLC()) mat = mat.compressMatLC(*krigopt.getMatLC(), true);
  return 0;
}

VectorDouble DriftList::evalDriftBySample(const Db* db,
                                          Id iech,
                                          const ECalcMember& member) const
{
  auto nbfl = getNDrift();
  VectorDouble drftab(nbfl);
  evalDriftBySampleInPlace(db, iech, member, drftab);
  return drftab;
}

void DriftList::evalDriftBySampleInPlace(const Db* db,
                                         Id iech,
                                         const ECalcMember& member,
                                         VectorDouble& drftab) const
{
  auto nbfl = getNDrift();
  if (nbfl != static_cast<Id>(drftab.size())) drftab.resize(nbfl);
  for (Id il = 0; il < nbfl; il++)
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
double DriftList::evalDriftValue(const Db* db,
                                 Id iech,
                                 Id ivar,
                                 Id ib,
                                 const ECalcMember& member) const
{
  auto nbfl    = getNDrift();
  double value = 0.;
  if (_flagCombined)
  {
    for (Id il = 0; il < nbfl; il++)
    {
      double local = evalDrift(db, iech, il, member);
      if (FFFF(local)) return TEST;
      value += local * _getDriftCL(ivar, il, ib);
    }
  }
  else
  {
    Id il = ib;
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
DriftList::evalDrift(const Db* db, Id iech, Id il, const ECalcMember& member) const
{
  if (member != ECalcMember::LHS && isDriftFiltered(il)) return 0.;
  if (!_isDriftIndexValid(il)) return TEST;
  return _drifts[il]->eval(db, iech);
}

const DriftList* DriftList::createReduce(const VectorInt& validVars) const
{
  Id ecr = 0;
  Id lec = 0;
  VectorBool valids(getNVar(), false);
  Id nvar = static_cast<Id>(validVars.size());
  VectorDouble mean(nvar, 0);

  for (Id ivar = 0; ivar < nvar; ivar++) valids[validVars[ivar]] = true;
  for (Id ivar = 0; ivar < getNVar(); ivar++)
  {
    if (valids[ivar])
    {
      mean[ecr++] = _mean[lec];
    }
    lec++;
  }
  auto* driftlist = new DriftList(_ctxt);
  driftlist->setMeans(mean);
  return driftlist;
}

void DriftList::setMeans(const VectorDouble& mean)
{
  if (_mean.size() == mean.size()) _mean = mean;
}

double DriftList::getMean(Id ivar) const
{
  if (ivar < 0 || ivar >= static_cast<Id>(_mean.size()))
  {
    messerr("Invalid argument in DriftList::getMean");
    return TEST;
  }
  return _mean[ivar];
}

/**
 * Define the Mean for one variable
 * @param mean Value for the mean
 * @param ivar Rank of the variable (starting from 0)
 */
void DriftList::setMean(const double mean, Id ivar)
{
  if (ivar < 0 || ivar >= static_cast<Id>(_mean.size()))
  {
    messerr("Invalid argument in DriftList::setMean - nothing changed");
    return;
  }
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
                                   Id iech,
                                   Id ivar,
                                   const VectorDouble& coeffs) const
{

  double drift = 0.;
  for (Id ib = 0, nfeq = getNDriftEquation(); ib < nfeq; ib++)
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

VectorDouble DriftList::evalDriftVarCoefs(const Db* db,
                                          const VectorDouble& coeffs,
                                          bool useSel) const
{
  VectorDouble vec;
  vec = evalDriftCoefs(db, coeffs, useSel);
  return vec;
}

} // namespace gstlrn