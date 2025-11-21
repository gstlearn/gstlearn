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
#include "Estimation/Vecchia.hpp"

#include "Basic/VectorHelper.hpp"
#include "Db/Db.hpp"
#include "Estimation/ALikelihood.hpp"
#include "LinearOp/CholeskyDense.hpp"
#include "LinearOp/CholeskySparse.hpp"
#include "Matrix/MatrixSparse.hpp"
#include "Matrix/MatrixSymmetric.hpp"
#include "Matrix/MatrixT.hpp"
#include "Model/ModelGeneric.hpp"
#include "Stats/Classical.hpp"
#include "Tree/Ball.hpp"
#include "geoslib_define.h"

namespace gstlrn
{
Vecchia::Vecchia(ModelGeneric* model,
                 Id nb_vecchia,
                 const Db* db1,
                 const Db* db2,
                 bool reml,
                 bool verbose)
  : ALikelihood(model, db1, reml)
  , _nbVecchia(nb_vecchia)
  , _db1(db1)
  , _db2(db2)
  , _Ranks()
  , _matCov()
  , _vectCov()
  , _work()
  , _DFull()
  , _LFull()
  , _NumberAbs1(0)
  , _NumberRel1(0)
  , _NumberRel2(0)
  , _cumulRanks1()
  , _cumulRanks2()
  , _varRanks1()
  , _varRanks2()
  , _varInverse1()
  , _varInverse2()
{
  setAuthorizedAnalyticalGradients(false);
  _chol = new CholeskyDense();
  _init(verbose);
}

Vecchia::Vecchia(const Vecchia& r)
  : ALikelihood(r)
  , _nbVecchia(r._nbVecchia)
  , _db1(r._db1)
  , _db2(r._db2)
  , _Ranks(r._Ranks)
  , _matCov(r._matCov)
  , _vectCov(r._vectCov)
  , _work(r._work)
  , _DFull(r._DFull)
  , _LFull(r._LFull)
  , _NumberAbs1(r._NumberAbs1)
  , _NumberRel1(r._NumberRel1)
  , _NumberRel2(r._NumberRel2)
  , _cumulRanks1(r._cumulRanks1)
  , _cumulRanks2(r._cumulRanks2)
  , _varRanks1(r._varRanks1)
  , _varRanks2(r._varRanks2)
  , _varInverse1(r._varInverse1)
  , _varInverse2(r._varInverse2)
{
  _chol = new CholeskyDense(*r._chol);
}
Vecchia& Vecchia::operator=(const Vecchia& r)
{
  if (this != &r)
  {
    ALikelihood::operator=(r);
    _nbVecchia   = r._nbVecchia;
    _db1         = r._db1;
    _db2         = r._db2;
    _Y           = r._Y;
    _DFull       = r._DFull;
    _LFull       = r._LFull;
    _NumberAbs1  = r._NumberAbs1;
    _NumberRel1  = r._NumberRel1;
    _NumberRel2  = r._NumberRel2;
    _cumulRanks1 = r._cumulRanks1;
    _cumulRanks2 = r._cumulRanks2;
    _varRanks1   = r._varRanks1;
    _varRanks2   = r._varRanks2;
    _varInverse1 = r._varInverse1;
    _varInverse2 = r._varInverse2;

    delete _chol;
    _chol = new CholeskyDense(*r._chol);
  }
  return *this;
}

Vecchia::~Vecchia()
{
  delete _chol;
}

void Vecchia::_init(bool verbose)
{
  _Ranks = findNN(_db1, _db2, _nbVecchia + 1, false, verbose);
}

Id Vecchia::_getCase() const
{
  Id icase = 0;
  if (_db1 != nullptr) icase = 1;
  if (_db2 != nullptr) icase = 2;
  if (_db1 != nullptr && _db2 != nullptr) icase = 2;
  return icase;
}

bool Vecchia::_identifyDbAndAbsoluteRank(const MatrixT<Id>& Ranks,
                                         Id irow,
                                         Id icol,
                                         Id* icaseDb,
                                         Id* ipAbs) const
{
  Id rank = Ranks(irow, icol);
  if (IFFFF(rank)) return false;
  if (rank < _NumberAbs1)
  {
    // Information belongs to the first Db
    *icaseDb = 1;
    *ipAbs   = rank;
  }
  else
  {
    *icaseDb = 2;
    *ipAbs   = rank - _NumberAbs1;
  }
  return true;
}

/**
 * @brief Returns the address in the covariance matrix for a given sample and variable
 *
 * @param icaseDb Rank of the current Db (1 or 2)
 * @param ipAbs Rank of the sample
 * @param ivar  Rank of the variable
 * @return Rank within the matrix
 */
Id Vecchia::_getAddressInMatrix(Id icaseDb, Id ipAbs, Id ivar) const
{
  Id irel = 0;
  if (icaseDb == 1)
  {
    if (ivar > 0) irel = _cumulRanks1[ivar];
    irel += _varInverse1[ivar][ipAbs];
  }
  else
  {
    irel = _NumberRel1;
    if (ivar > 0) irel = _cumulRanks2[ivar];
    irel += _varInverse2[ivar][ipAbs];
  }
  return irel;
}

Id Vecchia::_buildNeighborhood(const MatrixT<Id>& Ranks,
                               Id ndim,
                               Id icase1,
                               Id isample,
                               Id ivar,
                               Id nb_vecchia,
                               std::vector<std::array<Id, 4>>& neighDescr,
                               bool verbose) const
{
  // Loop on the ranks of the neighboring samples
  Id icaseDb = 0;
  Id ipAbs   = 0;
  Id ipAbs0  = 0;

  Id nitems = 0;
  for (Id jp = 0; jp < nb_vecchia; jp++)
  {
    // Identify the target Db and absolute sample rank; discard missing information
    if (!_identifyDbAndAbsoluteRank(Ranks, isample, jp + 1, &icaseDb, &ipAbs)) continue;

    for (Id jvar = 0; jvar <= ivar; jvar++)
    {
      // Identify the Db and absolute sample rank

      if (!_isVariableDefined(icaseDb, ipAbs, jvar)) continue;
      neighDescr[nitems][0] = icaseDb;
      neighDescr[nitems][1] = jvar;
      neighDescr[nitems][2] = _getAddressInMatrix(icaseDb, ipAbs, jvar);
      neighDescr[nitems][3] = ipAbs;
      nitems++;
    }
  }

  // For high variable rank, ensure that the same sample for previous variables is selected
  if (ivar > 0)
  {
    // Identify the target Db and absolute sample rank; discard missing information
    (void)_identifyDbAndAbsoluteRank(Ranks, isample, 0, &icaseDb, &ipAbs0);

    for (Id jvar = 0; jvar < ivar; jvar++)
    {
      if (!_isVariableDefined(icaseDb, ipAbs0, jvar)) continue;
      neighDescr[nitems][0] = icaseDb;
      neighDescr[nitems][1] = jvar;
      neighDescr[nitems][2] = _getAddressInMatrix(icaseDb, ipAbs0, jvar);
      neighDescr[nitems][3] = ipAbs0;
      nitems++;
    }
  }

  // Optional printout
  if (verbose)
  {
    VectorDouble coor(ndim);
    (void)_identifyDbAndAbsoluteRank(Ranks, isample, 0, &icase1, &ipAbs);
    if (icase1 == 1)
      _db1->getCoordinatesInPlace(coor, ipAbs);
    else
      _db2->getCoordinatesInPlace(coor, ipAbs);

    message("Row Number %4d (Db %d Var %d)", isample, icase1, ivar);
    message(" - Abs Rank = %4d", ipAbs);
    print_vector(" - Coors:", coor, true, false);
    if (nitems > 0) message(" Db | Var | Col Number | Abs Rank |           Coors\n");
    for (Id item = 0; item < nitems; item++)
    {
      auto icase2 = neighDescr[item][0];
      auto iabs2  = neighDescr[item][3];
      if (icase2 == 1)
        _db1->getCoordinatesInPlace(coor, iabs2);
      else
        _db2->getCoordinatesInPlace(coor, iabs2);
      message(" %2d |  %2d |    %4d    |   %4d   |", neighDescr[item][0], neighDescr[item][1],
              neighDescr[item][2], neighDescr[item][3]);
      print_vector("", coor, true, false);
    }
  }
  return nitems;
}

void Vecchia::_buildLHS(Id nitems,
                        const std::vector<std::array<Id, 4>>& neighDescr,
                        MatrixSymmetric& matCov)
{
  SpacePoint p1;
  SpacePoint p2;

  for (Id i1 = 0; i1 < nitems; i1++)
  {
    Id icase1 = neighDescr[i1][0];
    Id ivar1  = neighDescr[i1][1];
    Id iabs1  = neighDescr[i1][3];
    if (icase1 == 1)
      _db1->getSampleAsSPInPlace(p1, iabs1);
    else
      _db2->getSampleAsSPInPlace(p1, iabs1);

    for (Id i2 = 0; i2 <= i1; i2++)
    {
      Id icase2 = neighDescr[i2][0];
      Id ivar2  = neighDescr[i2][1];
      Id iabs2  = neighDescr[i2][3];
      if (icase2 == 1)
        _db1->getSampleAsSPInPlace(p2, iabs2);
      else
        _db2->getSampleAsSPInPlace(p2, iabs2);

      _model->getCov()->updateCovByPoints(icase1, iabs1, icase2, iabs2);

      double value = _model->evalCov(p1, p2, ivar1, ivar2);
      matCov.setValue(i1, i2, value);
    }

    // Update the Diagonal due to the presence of Variance of Measurement Error
    if (_db1->hasLocVariable(ELoc::V) && icase1 == 1)
    {
      Id icolVerr = _db1->getColIdxByLocator(ELoc::V, ivar1);

      double verr = 0.;

      if (icolVerr >= 0)
        verr = _db1->getValueByColIdx(iabs1, icolVerr);

      // Update the Covariance matrix
      if (verr > 0) matCov.updValue(i1, i1, EOperator::ADD, verr);
    }
  }
}

void Vecchia::_buildRHS(Id icase2,
                        Id iabs2,
                        Id ivar2,
                        Id nitems,
                        const std::vector<std::array<Id, 4>>& neighDescr,
                        MatrixDense& vectCov)
{
  SpacePoint p1;
  SpacePoint p2;

  if (icase2 == 1)
    _db1->getSampleAsSPInPlace(p2, iabs2);
  else
    _db2->getSampleAsSPInPlace(p2, iabs2);

  for (Id i1 = 0; i1 < nitems; i1++)
  {
    Id icase1 = neighDescr[i1][0];
    Id ivar1  = neighDescr[i1][1];
    Id iabs1  = neighDescr[i1][3];
    if (icase1 == 1)
      _db1->getSampleAsSPInPlace(p1, iabs1);
    else
      _db2->getSampleAsSPInPlace(p1, iabs1);

    _model->getCov()->updateCovByPoints(icase1, iabs1, icase2, iabs2);
    double value = _model->evalCov(p1, p2, ivar1, ivar2);
    vectCov.setValue(i1, 0, value);
  }
}

/**
 * @brief Store in the member _Y the list of data information
 * This information corresponds to the list of data defined values over the active samples
 * for all the variables
 */
void Vecchia::_loadDataFlattened()
{
  auto icase = _getCase();
  Id ecr     = 0;
  if (icase == 1)
  {
    _Y.resize(_NumberRel1);
    Id nvar = static_cast<Id>(_cumulRanks1.size());
    for (Id ivar = 0; ivar < nvar; ivar++)
      for (Id iech = 0, nech = static_cast<Id>(_varRanks1[ivar].size()); iech < nech; iech++)
        _Y[ecr++] = _db1->getLocVariable(ELoc::Z, _varRanks1[ivar][iech], ivar);
  }
  else
  {
    _Y.resize(_NumberRel2);
    Id nvar = static_cast<Id>(_cumulRanks2.size());
    for (Id ivar = 0; ivar < nvar; ivar++)
      for (Id iech = 0, nech = static_cast<Id>(_varRanks2[ivar].size()); iech < nech; iech++)
        _Y[ecr++] = _db2->getLocVariable(ELoc::Z, _varRanks2[ivar][iech], ivar);
  }
}

bool Vecchia::_isVariableDefined(Id icaseDb, Id ipAbs, Id ivar) const
{
  Id nvar = 0;
  const Db* db;
  if (icaseDb == 1)
    db = _db1;
  else
    db = _db2;
  nvar = _db->getNLoc(ELoc::Z);
  if (nvar <= 0) return true;

  double value = db->getLocVariable(ELoc::Z, ipAbs, ivar);
  return !FFFF(value);
}

void Vecchia::_convertAbsToRel(Id nech,
                               const VectorVectorInt& varRanks,
                               VectorVectorInt& varInverse)
{
  Id nvar = static_cast<Id>(varRanks.size());
  varInverse.resize(nvar);
  for (Id ivar = 0; ivar < nvar; ivar++)
  {
    varInverse[ivar].resize(nech, ITEST);
    for (Id irel = 0, n = static_cast<Id>(varRanks[ivar].size()); irel < n; irel++)
    {
      Id ipAbs                = varRanks[ivar][irel];
      varInverse[ivar][ipAbs] = irel;
    }
  }
}

double Vecchia::_buildC00(Id icaseDb,
                          Id ipAbs,
                          Id ivar)
{
  SpacePoint p1;
  if (icaseDb == 1)
    _db1->getSampleAsSPInPlace(p1, ipAbs);
  else
    _db2->getSampleAsSPInPlace(p1, ipAbs);

  if (_model->isNoStat())
  {
    _model->getCov()->updateCovByPoints(icaseDb, ipAbs, icaseDb, ipAbs);
  }
  double var0 = _model->getCov()->evalCov(p1, p1, ivar, ivar);

  if (icaseDb == 1)
  {
    if (_db1->hasLocVariable(ELoc::V))
    {
      Id icolVerr = _db1->getColIdxByLocator(ELoc::V, ivar);
      if (icolVerr >= 0)
        var0 += _db1->getValueByColIdx(ipAbs, icolVerr);
    }
  }
  return var0;
}

/**
 * @brief Construct the Vecchia approximation starting from 'Ranks'
 *
 * @param Ranks MatrixT<Id> which the ranks of the sample indices for each target
 * @param verbose Verbose flag
 * @return Id Error returned code
 *
 * @note The dimension of 'Rank' is:
 * - ncols = Dimension of the Neighborhood + 1
 * - nrows = Number of samples (dbin and dbout [optional])
 *
 * @note For each row, the first element of 'Ranks' is the sample number
 * - if smaller than N_dbin, it refers to the sample absolute rank in 'dbin'
 * - if larger, its value (after subtracting N_dbin) gives the sample absolute
 *   rank in 'dbout'
 */
Id Vecchia::computeLower(const MatrixT<Id>& Ranks, bool verbose)
{
  // Handle non stationarities if any
  // It has to be made before covariance computations since the model may have changed.
  _model->manage(_db1, _db2);

  // Preliminary check
  Id ndim;
  if (!haveSameNDim(_db1, _db2, _model, &ndim)) return 1;

  Id nvar;
  if (!haveCompatibleNVar(_db1, _db2, _model, &nvar)) return 1;

  if (verbose) mestitle(1, "Processing Neighborhood for constructing Lower Matrix for Vecchia");

  // 'nsample' corresponds to the count of active samples over 1 or 2 Dbs,
  // whatever the variable contents.
  Id nsample    = static_cast<Id>(Ranks.getNRows());
  Id nb_vecchia = static_cast<Id>(Ranks.getNCols()) - 1;

  _NumberRel1 = 0;
  if (_db1 != nullptr)
  {
    _NumberRel1 = _db1->getListOfSampleIndicesInPlace(nvar, _cumulRanks1, _varRanks1, true);
    _convertAbsToRel(_db1->getNSample(false), _varRanks1, _varInverse1);
  }
  _NumberRel2 = 0;
  if (_db2 != nullptr)
  {
    _NumberRel2 = _db2->getListOfSampleIndicesInPlace(nvar, _cumulRanks2, _varRanks2, true);
    _convertAbsToRel(_db2->getNSample(false), _varRanks2, _varInverse2);
  }
  _NumberAbs1 = (_db1 != nullptr) ? _db1->getNSample() : 0;

  // Resizing
  Id ntot = _NumberRel1 + _NumberRel2;
  _DFull.resize(ntot);
  _LFull = MatrixSparse(ntot, ntot, nb_vecchia + 1);

  // Loop on the samples
  Id nmax = nvar * (nb_vecchia + 1); // Multivariate neighborhood + Collocation
  std::vector<std::array<Id, 4>> neighDescr(nmax);
  Id icaseDb = 0; // Indication of the current Db (1 or 2)
  Id ipAbs   = 0; // Absolute sample rank in the current Db
  for (Id ivar = 0; ivar < nvar; ivar++)
  {
    // Loop over the samples referenced in 'Ranks' (over 1 or 2 Dbs)
    for (Id isample = 0; isample < nsample; isample++)
    {
      // Identify the Db and absolute sample rank
      (void)_identifyDbAndAbsoluteRank(Ranks, isample, 0, &icaseDb, &ipAbs);

      // Check if the sample corresponds to a valid (defined) variable value
      if (!_isVariableDefined(icaseDb, ipAbs, ivar)) continue;

      // Return the rank of the sample within the covariance matrix
      auto irel1 = _getAddressInMatrix(icaseDb, ipAbs, ivar);

      // Build the list of neighboring information
      Id nitems = _buildNeighborhood(Ranks, ndim, icaseDb, isample, ivar, nb_vecchia, neighDescr, verbose);

      // Fill the full matrix

      _LFull.setValue(irel1, irel1, 1.);
      double varK = _buildC00(icaseDb, ipAbs, ivar);

      if (nitems <= 0)
      {
        // Case with no previous information available
        _DFull[irel1] = 1. / varK;
      }
      else
      {
        // Constitute the local Kriging matrix
        _matCov.resize(nitems, nitems);
        _vectCov.resize(nitems, 1);
        _buildLHS(nitems, neighDescr, _matCov);
        _buildRHS(icaseDb, ipAbs, ivar, nitems, neighDescr, _vectCov);

        // Solve the local system
        _chol->setMatrix(_matCov);
        constvect vect = _vectCov.getViewOnColumn(0);
        _work.resize(vect.size());
        _chol->solve(vect, _work);

        // Patch the global matrix
        for (Id i = 0; i < nitems; i++)
        {
          Id irel2 = neighDescr[i][2];
          _LFull.setValue(irel2, irel1, -_work[i]);
        }
        _DFull[irel1] = 1. / (varK - VH::innerProductCV(_work, vect));
      }
    }
  }

  // Final calculations
  _LFull.transposeInPlace();

  // Optional printout
  if (verbose)
  {
    mestitle(1, "Lower Vecchia Matrix");
    _LFull.display();
    print_vector("Diagonal of Vecchia Matrix", _DFull, true, true);
  }
  return 0;
}

// Calculate LdY = Ldat %*% Y
VectorDouble Vecchia::calculateLdY(const VectorDouble& Y) const
{
  auto nd = getND();
  auto nt = getNT();
  VectorDouble LdY(nd);

  for (Id id = 0; id < nd; id++)
  {
    double value = 0.;
    for (Id jd = 0; jd < nd; jd++)
      value += getLFull(id + nt, jd + nt) * Y[jd];
    LdY[id] = value;
  }
  return LdY;
}

// Calculate FtLdY = Ft %*% Ldat %*% Y
VectorDouble Vecchia::calculateFtLdY(const VectorDouble& LdY) const
{
  auto nd = getND();
  auto nt = getNT();
  VectorDouble FtLdY(nt);
  for (Id it = 0; it < nt; it++)
  {
    double value = 0.;
    for (Id id = 0; id < nd; id++)
      value += getLFull(id + nt, it) * LdY[id];
    FtLdY[it] = value;
  }
  return FtLdY;
}

MatrixSparse* Vecchia::calculateW(const VectorDouble& D_dd) const
{
  auto nd = getND();
  auto nt = getNT();

  // Extract sub-part of 'Diagonal' vector
  VectorDouble D_tt(nt);
  VH::extractInPlace(getDFull(), D_tt, 0);

  // Extracting information from 'LFull'
  VectorInt indT(nt + nd, -1);
  for (Id it = 0; it < nt; it++) indT[it] = it;
  VectorInt indD(nt + nd, -1);
  for (Id id = 0; id < nd; id++) indD[id + nt] = id;
  MatrixSparse* Ltt = getLFull().extractSubmatrixByRanks(indT, indT);
  Ltt->forceDimension(nt, nt);
  MatrixSparse* Ldt = getLFull().extractSubmatrixByRanks(indD, indT);
  Ldt->forceDimension(nd, nt);

  /*! Product 't(A)' %*% 'M' %*% 'A' or 'A' %*% 'M' %*% 't(A)' */
  MatrixSparse* mat1 = prodNormMatVec(Ltt, D_tt, true);
  MatrixSparse* mat2 = prodNormMatVec(Ldt, D_dd, true);
  mat1->forceDimension(nt, nt);
  mat2->forceDimension(nt, nt);

  // Cleaning
  indT.clear();
  indD.clear();
  delete Ltt;
  delete Ldt;

  MatrixSparse* W = MatrixSparse::addMatMat(mat1, mat2);
  delete mat1;
  delete mat2;
  return W;
}

VectorDouble Vecchia::computeAndGetY()
{
  _loadDataFlattened();
  return _Y;
}
Id krigingVecchia(Db* dbin,
                  Db* dbout,
                  ModelGeneric* model,
                  Id nb_vecchia,
                  bool verbose,
                  const NamingConvention& namconv)
{
  Vecchia V(model, nb_vecchia, dbout, dbin);

  MatrixT<Id> Ranks = findNN(dbout, dbin, nb_vecchia + 1, false, verbose);
  if (V.computeLower(Ranks, verbose)) return 1;

  // Extract sub-part of 'Diagonal' vector
  const auto& DFull = V.getDFull();
  auto nd           = V.getND();
  auto nt           = V.getNT();
  Id nvar           = model->getNVar();
  VectorDouble D_dd(nd);
  VH::extractInPlace(DFull, D_dd, nt);

  // Calculate LdY
  const VectorDouble& Y = V.computeAndGetY();
  VectorDouble LdY      = V.calculateLdY(Y);
  LdY.multiply(D_dd);

  // Calculate FtLdY
  VectorDouble FtLdY = V.calculateFtLdY(LdY);

  // Calculating 'W'
  MatrixSparse* W = V.calculateW(D_dd);

  // Compute the Cholesky decomposition of 'W'
  CholeskySparse cholW(*W);
  if (!cholW.isReady())
  {
    messerr("Cholesky decomposition of Covariance matrix failed");
    return 1;
  }

  // Perform the estimation
  VectorDouble result = cholW.solveX(FtLdY);
  for (Id i = 0; i < nt; i++) result[i] = -result[i];

  // Saving the results
  Id iptr = dbout->addColumns(result, String(), ELoc::UNKNOWN, 0, true);
  namconv.setNamesAndLocators(dbout, iptr, "estim", nvar);

  return 0;
}

void Vecchia::productVecchia(constvect Y, vect res) const
{
  _LdY.resize(_LFull.getNRows());
  _LFull.prodMatVecInPlaceC(Y, _LdY, false);
  _LdY.multiply(_DFull);
  _LFull.prodMatVecInPlaceC(_LdY, res, true);
}

void Vecchia::productMatVecchia(const MatrixDense& X, MatrixDense& resmat) const
{
  auto nrows = X.getNRows();
  auto ncols = X.getNCols();
  resmat.resize(nrows, ncols);

  // Loop on the columns
  for (Id icol = 0; icol < ncols; icol++)
  {
    constvect colin = X.getViewOnColumn(icol);
    vect colout     = resmat.getViewOnColumnModify(icol);
    productVecchia(colin, colout);
  }
}

/**
 * Compute the log-likelihood (based on Vecchia approximation for covMat)
 *
 * @param db  Db structure where variable are loaded from
 * @param model ModelGeneric structure used for the calculation
 * @param nb_vecchia Number of neighbors to consider in the Vecchia approximation
 * @param flagPrint Print the output results
 * @param verbose Verbose flag
 *
 * @remarks The calculation considers all the active samples.
 * @remarks It can work in multivariate case with or without drift conditions (linked or not)
 * @remarks The algorithm is stopped (with a message) in the heterotopic case
 */
double logLikelihoodVecchia(const Db* db,
                            ModelGeneric* model,
                            Id nb_vecchia,
                            bool flagPrint,
                            bool verbose)
{
  Vecchia* vec  = Vecchia::createForOptim(model, db, nb_vecchia, false, verbose);
  double result = -vec->computeCost(flagPrint, verbose);
  delete vec;
  return result;
}

Vecchia* Vecchia::createForOptim(ModelGeneric* model,
                                 const Db* db,
                                 Id nb_vecchia,
                                 bool reml,
                                 bool verbose)
{

  auto* vec = new Vecchia(model, nb_vecchia, db, nullptr, reml, verbose);

  vec->_initLikelihood(verbose);
  return vec;
}

void Vecchia::_computeCm1X()
{
  productMatVecchia(_X, _Cm1X);
}

void Vecchia::_computeCm1Yc()
{
  _Cm1Yc.resize(_Yc.size());
  productVecchia(_Yc, _Cm1Yc);
}

double Vecchia::_computeLogDet() const
{
  return -VH::cumulLog(getDFull());
}

void Vecchia::_updateModel(bool verbose)
{
  computeLower(_Ranks, verbose);
}
} // namespace gstlrn
