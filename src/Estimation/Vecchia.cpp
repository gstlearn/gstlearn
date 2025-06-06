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
#include "Estimation/ALikelihood.hpp"
#include "LinearOp/CholeskySparse.hpp"
#include "Tree/Ball.hpp"
#include "Db/Db.hpp"
#include "LinearOp/CholeskyDense.hpp"
#include "Matrix/MatrixSparse.hpp"
#include "Matrix/MatrixSymmetric.hpp"
#include "Matrix/MatrixT.hpp"
#include "Model/ModelGeneric.hpp"
#include "geoslib_define.h"

Vecchia::Vecchia(const ModelGeneric* model,
                 int nb_neigh,
                 const Db* db1,
                 const Db* db2)
  : ALikelihood(model, db1)
  ,_nbNeigh(nb_neigh)
  , _db1(db1)
  , _db2(db2)
  , _DFull()
  , _LFull()
  , _dbTemp(nullptr)
  , _dbOnePoint(nullptr)
  , _Dmat()
{
  _nt   = (db1 == nullptr) ? 0 : db1->getNSample();
  _nd   = (db2 == nullptr) ? 0 : db2->getNSample();
  _chol = new CholeskyDense();
}

Vecchia::Vecchia(const Vecchia& r)
  : ALikelihood(r)
  , _nbNeigh(r._nbNeigh)
  , _db1(r._db1)
  , _db2(r._db2)
  , _DFull(r._DFull)
  , _LFull(r._LFull)
  , _dbTemp(nullptr)
  , _dbOnePoint(nullptr)
  , _Dmat(r._Dmat)
{
  if (r._dbTemp != nullptr) _dbTemp = r._dbTemp->clone();
  if (r._dbOnePoint != nullptr) _dbOnePoint = r._dbOnePoint->clone();
  _chol = new CholeskyDense(*r._chol);
}
Vecchia& Vecchia::operator=(const Vecchia& r)
{
  if (this != &r)
  {
    ALikelihood::operator=(r);
    _nbNeigh = r._nbNeigh;
    _db1     = r._db1;
    _db2     = r._db2;
    _DFull   = r._DFull;
    _LFull   = r._LFull;
    _Dmat    = r._Dmat;

    delete _dbTemp;
    delete _dbOnePoint;
    delete _chol;

    if (r._dbTemp != nullptr) _dbTemp = r._dbTemp->clone();
    if (r._dbOnePoint != nullptr) _dbOnePoint = r._dbOnePoint->clone();
    _chol = new CholeskyDense(*r._chol);
  }
  return *this;
}

Vecchia::~Vecchia()
{
  delete _dbTemp;
  delete _dbOnePoint;
  delete _chol;
}

void Vecchia::_init(bool verbose)
{
  _Ranks = findNN(_db, nullptr, _nbNeigh + 1, false, verbose);
}
/**
 * @brief Construct the Vecchia approximation starting from 'Ranks'
 *
 * @param Ranks MatrixT<int> which the ranks of the sample indices for each target
 * @param verbose Verbose flag
 * @return int Error returned code
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
int Vecchia::computeLower(const MatrixT<int>& Ranks, bool verbose)
{
  int ndim     = _model->getNDim();
  int ntot     = (int)Ranks.getNRows();
  int nb_neigh = (int)Ranks.getNCols() - 1;
  int Ndb1     = _db1->getNSample();
  double varK  = _model->eval0();
  double value;

  // Resizing
  _DFull.resize(ntot);
  //if (_LFull.empty())
    _LFull = MatrixSparse(ntot, ntot, nb_neigh + 1);
  //if (_Dmat.empty())
    _Dmat = MatrixSparse(ntot, ntot);

  // Creating empty Dbs
  if (_dbTemp == nullptr)
    _dbTemp = Db::createEmpty(nb_neigh, ndim, 0, 0, 0, false, true);
  if (_dbOnePoint == nullptr)
    _dbOnePoint = Db::createEmpty(1, ndim, 0);

  // Loop on the samples
  for (int ind = 0; ind < ntot; ind++)
  {
    int icur = 0;
    int iabs = Ranks(ind, 0);
    for (int idim = 0; idim < ndim; idim++)
    {
      if (iabs < Ndb1)
        value = _db1->getCoordinate(iabs, idim);
      else
        value = _db2->getCoordinate(iabs - Ndb1, idim);
      _dbOnePoint->setCoordinate(0, idim, value);
    }

    for (int jp = nb_neigh; jp >= 1; jp--)
    {
      int ip = Ranks(ind, jp);
      if (IFFFF(ip))
      {
        _dbTemp->setLocVariable(ELoc::SEL, icur, 0, 0.);
      }
      else
      {
        _dbTemp->setLocVariable(ELoc::SEL, icur, 0, 1.);
        for (int idim = 0; idim < ndim; idim++)
        {
          if (ip < Ndb1)
            value = _db1->getCoordinate(ip, idim);
          else
            value = _db2->getCoordinate(ip - Ndb1, idim);
          _dbTemp->setCoordinate(icur, idim, value);
        }
      }
      icur++;
    }

    if (_dbTemp->getNSample(true) <= 0)
    {
      _LFull.setValue(ind, ind, 1.);
      _DFull[ind] = 1. / varK;
    }
    else
    {

      _matCov.resize(nb_neigh, nb_neigh);
      _vectCov.resize(nb_neigh, 1);
      _model->evalCovMatSymInPlace(_matCov, _dbTemp);
      _chol->setMatrix(&_matCov);
      _model->evalCovMatInPlace(_vectCov, _dbTemp, _dbOnePoint);
      constvect vect = _vectCov.getViewOnColumn(0);
      _work.resize(vect.size());
      _chol->solve(vect, _work);

      icur = 0;
      _LFull.setValue(ind, ind, 1.);
      for (int jp = nb_neigh; jp >= 1; jp--)
      {
        int ip = Ranks(ind, jp);
        if (IFFFF(ip)) continue;
        _LFull.setValue(ip, ind, -_work[icur]);
        icur++;
      }
      _DFull[ind] = 1. / (varK - VH::innerProduct(_work, vect));
    }
  }
  _Dmat.setDiagonal(_DFull);
  _LFull.transposeInPlace();

  // Optional printout
  if (verbose)
  {
    message("Matrix L\n");
    _LFull.display();
    VH::dump("Diagonal D", _DFull);
  }
  return 0;
}

// Calculate LdY = Ldat %*% Y
VectorDouble Vecchia::calculateLdY(const VectorDouble& Y) const
{
  int nd = getND();
  int nt = getNT();
  VectorDouble LdY(nd);

  for (int id = 0; id < nd; id++)
  {
    double value = 0.;
    for (int jd = 0; jd < nd; jd++)
      value += getLFull(id + nt, jd + nt) * Y[jd];
    LdY[id] = value;
  }
  return LdY;
}

// Calculate FtLdY = Ft %*% Ldat %*% Y
VectorDouble Vecchia::calculateFtLdY(const VectorDouble& LdY) const
{
  int nd = getND();
  int nt = getNT();
  VectorDouble FtLdY(nt);
  for (int it = 0; it < nt; it++)
  {
    double value = 0.;
    for (int id = 0; id < nd; id++)
      value += getLFull(id + nt, it) * LdY[id];
    FtLdY[it] = value;
  }
  return FtLdY;
}

MatrixSparse* Vecchia::calculateW(const VectorDouble& D_dd) const
{
  int nd = getND();
  int nt = getNT();

  // Extract sub-part of 'Diagonal' vector
  VectorDouble D_tt(nt);
  VH::extractInPlace(getDFull(), D_tt, 0);

  // Extracting information from 'LFull'
  VectorInt indT = VectorInt(nt + nd, -1);
  for (int it = 0; it < nt; it++) indT[it] = it;
  VectorInt indD = VectorInt(nt + nd, -1);
  for (int id = 0; id < nd; id++) indD[id + nt] = id;
  MatrixSparse* Ltt = getLFull().extractSubmatrixByRanks(indT, indT);
  MatrixSparse* Ldt = getLFull().extractSubmatrixByRanks(indD, indT);

  /*! Product 't(A)' %*% 'M' %*% 'A' or 'A' %*% 'M' %*% 't(A)' */
  MatrixSparse* mat1 = prodNormMat(Ltt, D_tt, true);
  MatrixSparse* mat2 = prodNormMat(Ldt, D_dd, true);
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

int krigingVecchia(Db* dbin,
                   Db* dbout,
                   ModelGeneric* model,
                   int nb_neigh,
                   bool verbose,
                   const NamingConvention& namconv)
{
  Vecchia V = Vecchia(model, nb_neigh, dbout, dbin);

  MatrixT<int> Ranks = findNN(dbout, dbin, nb_neigh + 1, false, verbose);
  if (V.computeLower(Ranks, verbose)) return 1;

  // Extract sub-part of 'Diagonal' vector
  VectorDouble DFull = V.getDFull();
  int nd             = V.getND();
  int nt             = V.getNT();
  VectorDouble D_dd(nd);
  VH::extractInPlace(DFull, D_dd, nt);

  // Calculate LdY
  VectorDouble Y   = dbin->getColumnByLocator(ELoc::Z, 0);
  VectorDouble LdY = V.calculateLdY(Y);
  VH::multiplyInPlace(LdY, D_dd);

  // Calculate FtLdY
  VectorDouble FtLdY = V.calculateFtLdY(LdY);

  // Calculating 'W'
  MatrixSparse* W = V.calculateW(D_dd);

  // Compute the Cholesky decomposition of 'W'
  CholeskySparse cholW(W);
  if (!cholW.isReady())
  {
    messerr("Cholesky decomposition of Covariance matrix failed");
    delete W;
    return 1;
  }

  // Perform the estimation
  VectorDouble result = cholW.solveX(FtLdY);
  for (int i = 0; i < nt; i++) result[i] = -result[i];

  // Saving the results
  int iptr = dbout->addColumns(result, String(), ELoc::UNKNOWN, 0, true);
  namconv.setNamesAndLocators(dbout, iptr, "estim", 1);

  delete W;
  return 0;
}

void Vecchia::productVecchia(constvect Y, vect res) const
{
  _LdY.resize(_LFull.getNRows());
  _LFull.prodMatVecInPlace(Y, _LdY, false);
  VH::multiplyInPlace(_LdY, _DFull);
  _LFull.prodMatVecInPlace(_LdY, res, true);
}

void Vecchia::productMatVecchia(const MatrixDense& X, MatrixDense& resmat) const
{
  int nrows = X.getNRows();
  int ncols = X.getNCols();
  resmat.resize(nrows, ncols);

  // Loop on the columns
  for (int icol = 0; icol < ncols; icol++)
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
 * @param nb_neigh Number of neighbors to consider in the Vecchia approximation
 * @param verbose Verbose flag
 *
 * @remarks The calculation considers all the active samples.
 * @remarks It can work in multivariate case with or without drift conditions (linked or not)
 * @remarks The algorithm is stopped (with a message) in the heterotopic case
 */
double logLikelihoodVecchia(const Db* db,
                            ModelGeneric* model,
                            int nb_neigh,
                            bool verbose)
{
  Vecchia* vec = Vecchia::createForOptim(model, db, nb_neigh);
  double result = vec->computeCost(verbose);
  delete vec;
  return result;
}

Vecchia* Vecchia::createForOptim(ModelGeneric* model,
                                 const Db* db,
                                 int nb_neigh)
{
  auto* vec = new Vecchia(model, nb_neigh, db, nullptr);
  vec->init();
  return vec;
}

void Vecchia::_computeCm1X()
{
  productMatVecchia(_X, _Cm1X);
}

void Vecchia::_computeCm1Z()
{
  _Cm1Y.resize(_Y.size());
  productVecchia(_Y, _Cm1Y);
}

double Vecchia::_computeLogDet() const 
{
  return -VH::cumulLog(getDFull());
}

void Vecchia::_updateModel(bool verbose)
{
  computeLower(_Ranks, verbose);
 
}