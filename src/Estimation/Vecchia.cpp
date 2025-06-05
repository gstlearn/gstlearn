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
#include "LinearOp/CholeskySparse.hpp"
#include "Tree/Ball.hpp"
#include "Db/Db.hpp"
#include "LinearOp/CholeskyDense.hpp"
#include "Matrix/MatrixSparse.hpp"
#include "Matrix/MatrixSymmetric.hpp"
#include "Matrix/MatrixT.hpp"
#include "Matrix/MatrixFactory.hpp"
#include "Model/ModelGeneric.hpp"
#include "geoslib_define.h"

Vecchia::Vecchia(const ModelGeneric* model,
                 const Db* db1,
                 const Db* db2)
  : _db1(db1)
  , _db2(db2)
  , _model(model)
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

Vecchia::~Vecchia()
{
  delete _dbTemp;
  delete _dbOnePoint;
  delete _chol;
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
  if (_LFull.empty())
    _LFull = MatrixSparse(ntot, ntot, nb_neigh + 1);
  if (_Dmat.empty())
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
  Vecchia V = Vecchia(model, dbout, dbin);

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
  Vecchia* vec = Vecchia::createForOptim(model, db, nb_neigh, verbose);
  double result = vec->computeCost(verbose);
  delete vec;
  return result;
}

Vecchia* Vecchia::createForOptim(const ModelGeneric* model,
                                 const Db* db,
                                 int nb_neigh,
                                 bool verbose)
{
  auto* vec = new Vecchia(model, db, nullptr);

  int nvar = db->getNLoc(ELoc::Z);
  if (nvar < 1)
  {
    messerr("The 'db' should have at least one variable defined");
    return nullptr;
  }

  // Establish the vector of multivariate data
  int nDrift = model->getNDriftEquation();
  if (nDrift > 0)
    vec->_Y = db->getColumnsByLocator(ELoc::Z, true, true);
  else
    vec->_Y = db->getColumnsByLocator(ELoc::Z, true, true, model->getMeans());

  int size = (int)vec->_Y.size();
  if (verbose)
  {
    message("Likelihood calculation:\n");
    message("- Number of active samples     = %d\n", db->getNSample(true));
    message("- Number of variables          = %d\n", nvar);
    message("- Length of Information Vector = %d\n", size);
    if (nDrift > 0)
      message("- Number of drift conditions = %d\n", model->getNDriftEquation());
    else
      VH::dump("Constant Mean(s)", model->getMeans());
  }

  // Calculate the Vecchia approximation
  vec->_Ranks = findNN(db, nullptr, nb_neigh + 1, false);

  // If Drift function is present, evaluate the optimal Drift coefficients
  if (nDrift > 0)
  {
    // Extract the matrix of drifts at samples X
    vec->_X = model->evalDriftMat(db);

    vec->_beta.resize(nDrift);
  }
  return vec;
}

double Vecchia::computeCost(bool verbose)
{
  computeLower(_Ranks, verbose);
  int nDrift = _X.getNCols();

  if (nDrift > 0)
  {
    // Calculate t(L-1) %*% D-1 %*% L-1 applied to X (L and D from Vecchia)
    productMatVecchia(_X, _Cm1X);

    // Calculate XtCm1X = Xt * Cm1 * X
    MatrixSymmetric* XtCm1X =
      MatrixFactory::prodMatMat<MatrixSymmetric>(&_X, &_Cm1X, true, false);

    // Construct ZtCm1X = Zt * Cm1 * X and perform its Cholesky decomposition
    VectorDouble ZtCm1X = _Cm1X.prodVecMat(_Y);
    CholeskyDense XtCm1XChol(XtCm1X);
    if (!XtCm1XChol.isReady())
    {
      messerr("Cholesky decomposition of XtCm1X matrix failed");
      delete XtCm1X;
      return TEST;
    }

  // Calculate beta = (XtCm1X)-1 * ZtCm1X
    if (XtCm1XChol.solve(ZtCm1X, _beta))
    {
      messerr("Error when calculating Likelihood");
      delete XtCm1X;
      return TEST;
    }
    // model->setBetaHat(beta);
    delete XtCm1X;

    if (verbose)
    {
      VH::dump("Optimal Drift coefficients = ", _beta);
    }

    // Center the data by the optimal drift: Y = Y - beta * X
    VH::subtractInPlace(_Y, _X.prodMatVec(_beta));
  }

  // Calculate t(L-1) %*% D-1 %*% L-1 applied to Y (L and D from Vecchia)
  _Cm1Z.resize(_Y.size());
  productVecchia(_Y, _Cm1Z);

// Calculate the log-determinant
  double logdet = -VH::cumulLog(getDFull());

  // Calculate quad = Zt * Cm1Z
  double quad = VH::innerProduct(_Y, _Cm1Z);

  // Derive the log-likelihood
  int size = (int)_Y.size();
  double loglike = -0.5 * (logdet + quad + size * log(2. * GV_PI));

  // Optional printout
  if (verbose)
  {
    message("Log-Determinant = %lf\n", logdet);
    message("Quadratic term  = %lf\n", quad);
    message("Log-likelihood  = %lf\n", loglike);
  }
  return loglike;
}