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
  Vecchia::Vecchia(
    ModelGeneric* model,
    Id nb_vecchia,
    const Db* dbin,
    const Db* dbout,
    bool reml,
    bool verbose)
    : ALikelihood(model, (dbin && dbout) ? dbin : dbout, reml)
    , _nbVecchia(nb_vecchia)
    , _dbout(dbout)
    , _dbin(dbin)
    , _Ranks()
    , _matCov()
    , _vectCov()
    , _work()
    , _DFull()
    , _LFull()
    , _NumberShift(0)
    , _NumberRelOut(0)
    , _NumberRelIn(0)
    , _cumulRanksOut()
    , _cumulRanksIn()
    , _varRanksOut()
    , _varRanksIn()
    , _varInverseOut()
    , _varInverseIn()
  {
    setAuthorizedAnalyticalGradients(false);
    _chol = new CholeskyDense();
    _init(verbose);

    // The likelihood is initialized here, since the calculation of the Drift coefficients
    // may be used in Vecchia for centering the input variable and decentering the results
    initLikelihood(verbose);
  }

  Vecchia::Vecchia(const Vecchia& r)
    : ALikelihood(r)
    , _nbVecchia(r._nbVecchia)
    , _dbout(r._dbout)
    , _dbin(r._dbin)
    , _Ranks(r._Ranks)
    , _matCov(r._matCov)
    , _vectCov(r._vectCov)
    , _work(r._work)
    , _LdY(r._LdY)
    , _DFull(r._DFull)
    , _LFull(r._LFull)
    , _Qmat(r._Qmat)
    , _chol(nullptr)
    , _NumberShift(r._NumberShift)
    , _NumberRelOut(r._NumberRelOut)
    , _NumberRelIn(r._NumberRelIn)
    , _cumulRanksOut(r._cumulRanksOut)
    , _cumulRanksIn(r._cumulRanksIn)
    , _varRanksOut(r._varRanksOut)
    , _varRanksIn(r._varRanksIn)
    , _varInverseOut(r._varInverseOut)
    , _varInverseIn(r._varInverseIn)
  {
    _chol = new CholeskyDense(*r._chol);
  }

  Vecchia& Vecchia::operator=(const Vecchia& r)
  {
    if (this != &r)
    {
      ALikelihood::operator=(r);
      _nbVecchia = r._nbVecchia;
      _dbout = r._dbout;
      _dbin = r._dbin;
      _Ranks = r._Ranks;
      _matCov = r._matCov;
      _vectCov = r._vectCov;
      _work = r._work;
      _LdY = r._LdY;
      _DFull = r._DFull;
      _LFull = r._LFull;
      _Qmat = r._Qmat;
      _NumberShift = r._NumberShift;
      _NumberRelOut = r._NumberRelOut;
      _NumberRelIn = r._NumberRelIn;
      _cumulRanksOut = r._cumulRanksOut;
      _cumulRanksIn = r._cumulRanksIn;
      _varRanksOut = r._varRanksOut;
      _varRanksIn = r._varRanksIn;
      _varInverseOut = r._varInverseOut;
      _varInverseIn = r._varInverseIn;

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
    _Ranks = findNN(_dbout, _dbin, _nbVecchia + 1, false, verbose);
  }

  Id Vecchia::_getCase() const
  {
    Id icase = 0;
    if (_dbout != nullptr) icase = 1;
    if (_dbin != nullptr) icase = 2;
    if (_dbout != nullptr && _dbin != nullptr) icase = 2;
    return icase;
  }

  MatrixSparse& Vecchia::getQMat() const
  {
    _Qmat.prodNormMatVecInPlace(&_LFull, _DFull, true);
    return _Qmat;
  }

  bool Vecchia::_identifyDbAndAbsoluteRank(
    const MatrixT<Id>& Ranks,
    Id irow,
    Id icol,
    Id* icaseDb,
    Id* ipAbs) const
  {
    Id rank = Ranks(irow, icol);
    if (isNA(rank)) return false;
    if (rank < _NumberShift)
    {
      // Information belongs to the first Db
      *icaseDb = 1;
      *ipAbs = rank;
    }
    else
    {
      *icaseDb = 2;
      *ipAbs = rank - _NumberShift;
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
      if (ivar > 0) irel = _cumulRanksOut[ivar];
      irel += _varInverseOut[ivar][ipAbs];
    }
    else
    {
      irel = _NumberRelOut;
      if (ivar > 0) irel = _cumulRanksIn[ivar];
      irel += _varInverseIn[ivar][ipAbs];
    }
    return irel;
  }

  Id Vecchia::_buildNeighborhood(
    const MatrixT<Id>& Ranks,
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
    Id ipAbs = 0;
    Id ipAbs0 = 0;

    Id nitems = 0;
    for (Id jp = 0; jp < nb_vecchia; jp++)
    {
      // Identify the target Db and absolute sample rank; discard missing information
      if (!_identifyDbAndAbsoluteRank(Ranks, isample, jp + 1, &icaseDb, &ipAbs))
        continue;

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
        _dbout->getCoordinatesInPlace(coor, ipAbs);
      else
        _dbin->getCoordinatesInPlace(coor, ipAbs);

      message("Row Number %4d (Db %d Var %d)", isample, icase1, ivar);
      message(" - Abs Rank = %4d", ipAbs);
      printVector(coor, " - Coors:", true, false);
      if (nitems > 0)
        message(" Db | Var | Col Number | Abs Rank |           Coors\n");
      for (Id item = 0; item < nitems; item++)
      {
        auto icase2 = neighDescr[item][0];
        auto iabs2 = neighDescr[item][3];
        if (icase2 == 1)
          _dbout->getCoordinatesInPlace(coor, iabs2);
        else
          _dbin->getCoordinatesInPlace(coor, iabs2);
        message(
          " %2d |  %2d |    %4d    |   %4d   |", neighDescr[item][0],
          neighDescr[item][1], neighDescr[item][2], neighDescr[item][3]);
        printVector(coor, "", true, false);
      }
    }
    return nitems;
  }

  void Vecchia::_buildLHS(
    Id nitems,
    const std::vector<std::array<Id, 4>>& neighDescr,
    MatrixSymmetric& matCov)
  {
    SpacePoint p1(_model->getSpace());
    SpacePoint p2(_model->getSpace());

    for (Id i1 = 0; i1 < nitems; i1++)
    {
      Id icase1 = neighDescr[i1][0];
      Id ivar1 = neighDescr[i1][1];
      Id iabs1 = neighDescr[i1][3];
      if (icase1 == 1)
        _dbout->getSampleAsSPInPlace(p1, iabs1);
      else
        _dbin->getSampleAsSPInPlace(p1, iabs1);

      for (Id i2 = 0; i2 <= i1; i2++)
      {
        Id icase2 = neighDescr[i2][0];
        Id ivar2 = neighDescr[i2][1];
        Id iabs2 = neighDescr[i2][3];
        if (icase2 == 1)
          _dbout->getSampleAsSPInPlace(p2, iabs2);
        else
          _dbin->getSampleAsSPInPlace(p2, iabs2);

        _model->getCov()->updateCovByPoints(icase1, iabs1, icase2, iabs2);

        double value = _model->evalCov(p1, p2, ivar1, ivar2);
        matCov.setValue(i1, i2, value);
      }

      // Update the Diagonal due to the presence of Variance of Measurement Error
      if (_dbout->hasLocVariable(ELoc::V) && icase1 == 1)
      {
        Id icolVerr = _dbout->getColIdxByLocator(ELoc::V, ivar1);

        double verr = 0.;

        if (icolVerr >= 0) verr = _dbout->getValueByColIdx(iabs1, icolVerr);

        // Update the Covariance matrix
        if (verr > 0) matCov.updValue(i1, i1, EOperator::ADD, verr);
      }
    }
  }

  void Vecchia::_buildRHS(
    Id icase2,
    Id iabs2,
    Id ivar2,
    Id nitems,
    const std::vector<std::array<Id, 4>>& neighDescr,
    MatrixDense& vectCov)
  {
    SpacePoint p1(_model->getSpace());
    SpacePoint p2(_model->getSpace());

    if (icase2 == 1)
      _dbout->getSampleAsSPInPlace(p2, iabs2);
    else
      _dbin->getSampleAsSPInPlace(p2, iabs2);

    for (Id i1 = 0; i1 < nitems; i1++)
    {
      Id icase1 = neighDescr[i1][0];
      Id ivar1 = neighDescr[i1][1];
      Id iabs1 = neighDescr[i1][3];
      if (icase1 == 1)
        _dbout->getSampleAsSPInPlace(p1, iabs1);
      else
        _dbin->getSampleAsSPInPlace(p1, iabs1);

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
    Id ecr = 0;
    Id nech;
    if (icase == 1)
    {
      _Y.resize(_NumberRelOut);
      Id nvar = static_cast<Id>(_cumulRanksOut.size());
      for (Id ivar = 0; ivar < nvar; ivar++)
      {
        nech = static_cast<Id>(_varRanksOut[ivar].size());
        for (Id iech = 0; iech < nech; iech++)
          _Y[ecr++] =
            _dbout->getLocVariable(ELoc::Z, _varRanksOut[ivar][iech], ivar);
      }
    }
    else
    {
      _Y.resize(_NumberRelIn);
      Id nvar = static_cast<Id>(_cumulRanksIn.size());
      for (Id ivar = 0; ivar < nvar; ivar++)
      {
        nech = static_cast<Id>(_varRanksIn[ivar].size());
        for (Id iech = 0; iech < nech; iech++)
          _Y[ecr++] =
            _dbin->getLocVariable(ELoc::Z, _varRanksIn[ivar][iech], ivar);
      }
    }
  }

  /**
   * @brief Suppress a drift from the Data vector
   */
  void Vecchia::centerDataInPlace()
  {
    if (!_model->hasDrift())
    {
      // The Model has no Drift, simply remove the constant mean from the data
      auto nvar = static_cast<Id>(_cumulRanksIn.size());
      for (Id ivar = 0, ecr = 0; ivar < nvar; ivar++)
      {
        double meanValue = _model->getMean(ivar);
        auto nech = static_cast<Id>(_varRanksIn[ivar].size());
        for (Id iech = 0; iech < nech; iech++) _Y[ecr++] -= meanValue;
      }
    }
    else
    {
      // Evaluate the Drift matrix at the Data points ('dbin')
      auto X = _model->evalDriftMat(_dbin);
      // Center the data by the optimal drift: Y = Y - beta * X
      VH::subtractInPlace(AMatrix::product(X, _model->getBetaHat()), _Y, _Y);
    }
  }

  /**
   * @brief Add a drift to the Result vector
   *
   * @param tab: VectorDouble to be modified
   */
  void Vecchia::uncenterResultsInPlace(VectorDouble& tab)
  {
    if (!_model->hasDrift())
    {
      auto nt = getNT();
      auto nvar = _model->getNVar();
      for (Id ivar = 0, ecr = 0; ivar < nvar; ivar++)
      {
        Id nech = nt / nvar;
        double meanValue = _model->getMean(ivar);
        for (Id iech = 0; iech < nech; iech++) tab[ecr++] += meanValue;
      }
    }
    else
    {
      // Evaluate the Drift matrix at the Target points ('dbout')
      auto X = _model->evalDriftMat(_dbout);
      // Uncenter the results by the optimal drift: Y = Y + beta * X
      VH::addInPlace(AMatrix::product(X, _model->getBetaHat()), tab);
    }
  }

  bool Vecchia::_isVariableDefined(Id icaseDb, Id ipAbs, Id ivar) const
  {
    const Db* db;
    if (icaseDb == 1)
      db = _dbout;
    else
      db = _dbin;
    Id nvar = db->getNLoc(ELoc::Z);
    if (nvar <= 0) return true;
    double value = db->getLocVariable(ELoc::Z, ipAbs, ivar);
    return !FFFF(value);
  }

  void Vecchia::_convertAbsToRel(
    Id nech,
    const VectorVectorInt& varRanks,
    VectorVectorInt& varInverse)
  {
    Id nvar = static_cast<Id>(varRanks.size());
    varInverse.resize(nvar);
    for (Id ivar = 0; ivar < nvar; ivar++)
    {
      varInverse[ivar].resize(nech, ITEST);
      for (Id irel = 0, n = static_cast<Id>(varRanks[ivar].size()); irel < n;
           irel++)
      {
        Id ipAbs = varRanks[ivar][irel];
        varInverse[ivar][ipAbs] = irel;
      }
    }
  }

  double Vecchia::_buildC00(Id icaseDb, Id ipAbs, Id ivar)
  {
    SpacePoint p1(_model->getSpace());
    if (icaseDb == 1)
      _dbout->getSampleAsSPInPlace(p1, ipAbs);
    else
      _dbin->getSampleAsSPInPlace(p1, ipAbs);

    if (_model->isNoStat())
    {
      _model->getCov()->updateCovByPoints(icaseDb, ipAbs, icaseDb, ipAbs);
    }
    double var0 = _model->getCov()->evalCov(p1, p1, ivar, ivar);

    if (icaseDb == 1)
    {
      if (_dbout->hasLocVariable(ELoc::V))
      {
        Id icolVerr = _dbout->getColIdxByLocator(ELoc::V, ivar);
        if (icolVerr >= 0) var0 += _dbout->getValueByColIdx(ipAbs, icolVerr);
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
    _model->manage(_dbout, _dbin);

    // Preliminary check
    Id ndim;
    if (!haveSameNDim(_dbout, _dbin, _model, &ndim)) return 1;

    Id nvar;
    if (!haveCompatibleNVar(_dbout, _dbin, _model, &nvar)) return 1;

    if (verbose)
      mestitle(
        1, "Processing Neighborhood for constructing Lower Matrix for Vecchia");

    // 'nsample' corresponds to the count of active samples over 1 or 2 Dbs,
    // whatever the variable contents.
    Id nsample = static_cast<Id>(Ranks.getNRows());
    Id nb_vecchia = static_cast<Id>(Ranks.getNCols()) - 1;

    _NumberRelOut = 0;
    if (_dbout != nullptr)
    {
      _NumberRelOut = _dbout->getListOfSampleIndicesInPlace(
        nvar, _cumulRanksOut, _varRanksOut, true, false);
      _convertAbsToRel(_dbout->getNSample(false), _varRanksOut, _varInverseOut);
    }
    _NumberRelIn = 0;
    if (_dbin != nullptr)
    {
      _NumberRelIn = _dbin->getListOfSampleIndicesInPlace(
        nvar, _cumulRanksIn, _varRanksIn, true, true);
      _convertAbsToRel(_dbin->getNSample(false), _varRanksIn, _varInverseIn);
    }
    _NumberShift = (_dbout != nullptr) ? _dbout->getNSample() : 0;

    // Resizing
    Id ntot = _NumberRelOut + _NumberRelIn;
    _DFull.resize(ntot);
    _LFull = MatrixSparse(ntot, ntot, nb_vecchia + 1);

    // Loop on the samples
    Id nmax =
      nvar * (nb_vecchia + 1); // Multivariate neighborhood + Collocation
    std::vector<std::array<Id, 4>> neighDescr(nmax);
    Id icaseDb = 0; // Indication of the current Db (1 or 2)
    Id ipAbs = 0; // Absolute sample rank in the current Db
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
        Id nitems = _buildNeighborhood(
          Ranks, ndim, icaseDb, isample, ivar, nb_vecchia, neighDescr, verbose);

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
          constvect vec = _vectCov.getViewOnColumn(0);
          _work.resize(vec.size());
          _chol->solve(vec, _work);

          // Patch the global matrix
          for (Id i = 0; i < nitems; i++)
          {
            Id irel2 = neighDescr[i][2];
            _LFull.setValue(irel2, irel1, -_work[i]);
          }
          _DFull[irel1] = 1. / (varK - VH::innerProductCV(_work, vec));
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
      printVector(_DFull, "Diagonal of Vecchia Matrix", true, true);
    }
    return 0;
  }

  // Calculate LdY = Ldat %*% Y
  VectorDouble Vecchia::calculateLdY(const VectorDouble& Y) const
  {
    auto nd = getND();
    auto nt = getNT();
    VectorDouble LdY(nd);
    double value;

    for (Id id = 0; id < nd; id++)
    {
      value = 0.;
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
    double value;

    for (Id it = 0; it < nt; it++)
    {
      value = 0.;
      for (Id id = 0; id < nd; id++) value += getLFull(id + nt, it) * LdY[id];
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
    // Load the Data from the dbin and dbout into the vector _Y
    _loadDataFlattened();

    // Centering the Data by the Mean or Drift
    centerDataInPlace();

    return _Y;
  }

  void Vecchia::productVecchia(constvect Y, vect res) const
  {
    _LdY.resize(_LFull.getNRows());
    AMatrix::productInPlace(_LdY.asVect(), _LFull, Y, false);
    _LdY *= _DFull;
    AMatrix::productInPlace(res, _LFull, _LdY.asConstVect(), true);
  }

  void
    Vecchia::productMatVecchia(const MatrixDense& X, MatrixDense& resmat) const
  {
    auto nrows = X.getNRows();
    auto ncols = X.getNCols();
    resmat.resize(nrows, ncols);

    // Loop on the columns
    for (Id icol = 0; icol < ncols; icol++)
    {
      constvect colin = X.getViewOnColumn(icol);
      vect colout = resmat.getViewOnColumnModify(icol);
      productVecchia(colin, colout);
    }
  }

  Vecchia* Vecchia::createForOptim(
    ModelGeneric* model,
    const Db* db,
    Id nb_vecchia,
    bool reml,
    bool verbose)
  {
    // The 'db' is passed as the 'dbout' in the Vecchia constructor, since this is the compulsory Data Base
    // However it is the one used for the likelihood calculation.
    auto* vec = new Vecchia(model, nb_vecchia, nullptr, db, reml, verbose);
    vec->_initLikelihoodForOptim(verbose);
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

  void Vecchia::_solveQ(constvect inv, vect outv) const
  {
    productVecchia(inv, outv);
  }

  /**
   * Perform the estimation by Kriging (based on Vecchia approximation for covMat)
   *
   * @param dbin  Db structure where variable are loaded from
   * @param dbout Db structure where results are stored
   * @param model ModelGeneric structure used for the calculation
   * @param nb_vecchia Number of neighbors to consider in the Vecchia approximation
   * @param verbose Verbose flag
   * @param namconv NamingConvention structure used to store the results in 'dbout'
   */
  Id krigingVecchia(
    Db* dbin,
    Db* dbout,
    ModelGeneric* model,
    Id nb_vecchia,
    bool verbose,
    const NamingConvention& namconv)
  {
    // Instantiate a small Vecchia based on dbin alone
    // This is meant to calculate the optimal drift coefficients (which will be stored in model)
    // It will be used to center the data and uncenter the results
    Vecchia Vdat(model, nb_vecchia, nullptr, dbin, false, verbose);
    Vdat.updateModel(verbose);
    if (model->hasDrift()) Vdat.calculateBeta(verbose);

    Vecchia V(model, nb_vecchia, dbin, dbout);
    MatrixT<Id> Ranks = findNN(dbout, dbin, nb_vecchia + 1, false, verbose);
    if (V.computeLower(Ranks, verbose)) return 1;

    // Extract sub-part of 'Diagonal' vector
    const auto& DFull = V.getDFull();
    auto nd = V.getND();
    auto nt = V.getNT();
    Id nvar = model->getNVar();
    VectorDouble D_dd(nd);
    VH::extractInPlace(DFull, D_dd, nt);

    // Calculate LdY
    const VectorDouble& Y = V.computeAndGetY();
    auto LdY = V.calculateLdY(Y);
    LdY *= D_dd;

    // Calculate FtLdY
    auto FtLdY = V.calculateFtLdY(LdY);

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

    // UnCentering the result by the Mean
    V.uncenterResultsInPlace(result);

    // Saving the results
    Id iptr = dbout->addColumns(result, String(), ELoc::UNDEFINED, 0, true);
    namconv.setNamesAndLocators(
      dbin, VectorString(), ELoc::Z, nvar, dbout, iptr, "estim");

    return 0;
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
  double logLikelihoodVecchia(
    const Db* db,
    ModelGeneric* model,
    Id nb_vecchia,
    bool flagPrint,
    bool verbose)
  {
    // The 'db' is passed as the 'dbout' in the Vecchia constructor, since this is the compulsory Data Base
    // However it is the one used for the likelihood calculation.
    auto* vec = new Vecchia(model, nb_vecchia, nullptr, db, false, verbose);
    double result = -vec->computeCost(flagPrint, verbose);
    delete vec;
    return result;
  }

} // namespace gstlrn
