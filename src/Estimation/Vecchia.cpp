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

Vecchia::Vecchia(const ModelGeneric* model,
                 const Db* db1,
                 const Db* db2)
  : _db1(db1)
  , _db2(db2)
  , _model(model)
  , _DFull()
  , _LFull()
  , _Dmat()
{
  _nt = (db1 == nullptr) ? 0 : db1->getNSample();
  _nd = (db2 == nullptr) ? 0 : db2->getNSample();
}

Vecchia::~Vecchia()
{
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
  int ntot     = (int) Ranks.getNRows();
  int nb_neigh = (int) Ranks.getNCols() - 1;
  int Ndb1     = _db1->getNSample();
  double varK  = _model->eval0();
  double value;

  // Resizing
  _DFull.resize(ntot);
  _LFull = MatrixSparse(ntot, ntot, nb_neigh+1);
  _Dmat = MatrixSparse(ntot, ntot);

  // Creating empty Dbs
  Db* Dbtemp = Db::createEmpty(nb_neigh, ndim, 0, 0, 0, false, true);
  Db* DbOnePoint = Db::createEmpty(1, ndim, 0);

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
      DbOnePoint->setCoordinate(0, idim, value);
    }

    for (int jp = nb_neigh; jp >= 1; jp--)
    {
      int ip = Ranks(ind, jp);
      if (IFFFF(ip))
      {
        Dbtemp->setLocVariable(ELoc::SEL, icur, 0, 0.);
      }
      else
      {
        Dbtemp->setLocVariable(ELoc::SEL, icur, 0, 1.);
        for (int idim = 0; idim < ndim; idim++)
        {
          if (ip < Ndb1)
            value = _db1->getCoordinate(ip, idim);
          else
            value = _db2->getCoordinate(ip - Ndb1, idim);
          Dbtemp->setCoordinate(icur, idim, value);
        }
      }
      icur++;
    }

    if (Dbtemp->getNSample(true) <= 0)
    {
      _LFull.setValue(ind, ind, 1.);
      _DFull[ind] = 1. / varK;
    }
    else
    {
      MatrixSymmetric mat = _model->evalCovMatSym(Dbtemp);
      CholeskyDense chol         = CholeskyDense(&mat);
      MatrixDense crossmat = _model->evalCovMat(Dbtemp, DbOnePoint);
      VectorDouble vect = crossmat.getColumn(0);
      VectorDouble res  = chol.solveX(vect);

      icur = 0;
      _LFull.setValue(ind, ind, 1.);
      for (int jp = nb_neigh; jp >= 1; jp--)
      {
        int ip = Ranks(ind, jp);
        if (IFFFF(ip)) continue;
        _LFull.setValue(ip, ind, -res[icur]);
        icur++;
      }
      _DFull[ind] = 1. / (varK - VH::innerProduct(res, vect));
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
static VectorDouble _calculateLdY(const Vecchia& V, const VectorDouble& Y)
{
  int nd = V.getND();
  int nt = V.getNT();
  VectorDouble LdY(nd);

  for (int id = 0; id < nd; id++)
  {
    double value = 0.;
    for (int jd = 0; jd < nd; jd++)
      value += V.getLFull(id + nt, jd + nt) * Y[jd];
    LdY[id] = value;
  }
  return LdY;
}

// Calculate FtLdY = Ft %*% Ldat %*% Y
static VectorDouble _calculateFtLdY(const Vecchia& V, const VectorDouble& LdY)
{
  int nd = V.getND();
  int nt = V.getNT();
  VectorDouble FtLdY(nt);
  for (int it = 0; it < nt; it++)
  {
    double value = 0.;
    for (int id = 0; id < nd; id++)
      value += V.getLFull(id + nt, it) * LdY[id];
    FtLdY[it] = value;
  }
  return FtLdY;
}

static MatrixSparse* _calculateW(const Vecchia& V, const VectorDouble& D_dd)
{
  int nd = V.getND();
  int nt = V.getNT();

  // Extract sub-part of 'Diagonal' vector
  VectorDouble D_tt(nt);
  VH::extractInPlace(V.getDFull(), D_tt, 0);

  // Extracting information from 'LFull'
  VectorInt indT = VectorInt(nt + nd, -1);
  for (int it = 0; it < nt; it++) indT[it] = it;
  VectorInt indD = VectorInt(nt + nd, -1);
  for (int id = 0; id < nd; id++) indD[id + nt] = id;
  MatrixSparse* Ltt = V.getLFull().extractSubmatrixByRanks(indT, indT);
  MatrixSparse* Ldt = V.getLFull().extractSubmatrixByRanks(indD, indT);

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
  VectorDouble LdY = _calculateLdY(V, Y);
  VH::multiplyInPlace(LdY, D_dd);

  // Calculate FtLdY
  VectorDouble FtLdY = _calculateFtLdY(V, LdY);

  // Calculating 'W'
  MatrixSparse* W = _calculateW(V, D_dd);

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

static VectorDouble _productVecchia(const Vecchia& V, const VectorDouble& Y)
{
  const MatrixSparse& L = V.getLFull();
  const VectorDouble& D = V.getDFull();

  VectorDouble LdY = L.prodMatVec(Y, false);
  VectorDouble DLdY = VH::multiply(D, LdY);
  VectorDouble LDLdY = L.prodMatVec(DLdY, true);
  return LDLdY;
}

static MatrixDense _productMatVecchia(const Vecchia& V, const MatrixDense& X)
{
  int nrows = X.getNRows();
  int ncols = X.getNCols();
  MatrixDense resmat(nrows, ncols);
  VectorDouble colin(nrows);
  VectorDouble colout(nrows);

  // Loop on the columns
  for (int icol = 0; icol < ncols; icol++)
  {
    VectorDouble colin = X.getColumn(icol);
    colout = _productVecchia(V, colin);
    resmat.setColumn(icol, colout);
  }
  return resmat;
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
double logLikelihoodVecchia(Db* db,
                            ModelGeneric* model,
                            int nb_neigh,
                            bool verbose)
{
  int nvar = db->getNLoc(ELoc::Z);
  if (nvar < 1)
  {
    messerr("The 'db' should have at least one variable defined");
    return TEST;
  }
 
  // Establish the vector of multivariate data
  VectorDouble Y;
  int nDrift = model->getNDriftEquation();
  if (nDrift > 0)
    Y = db->getColumnsByLocator(ELoc::Z, true, true);
  else
    Y = db->getColumnsByLocator(ELoc::Z, true, true, model->getMeans());

  int size = (int)Y.size();
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
  Vecchia V = Vecchia(model, db);
  MatrixT<int> Ranks = findNN(db, nullptr, nb_neigh + 1, false);

  if (V.computeLower(Ranks)) return 1;
  // If Drift functions are present, evaluate the optimal Drift coefficients
  if (nDrift > 0)
  {
    // Extract the matrix of drifts at samples X
    MatrixDense X = model->evalDriftMat(db);

    // Calculate t(L-1) %*% D-1 %*% L-1 applied to X (L and D from Vecchia)
    MatrixDense Cm1X = _productMatVecchia(V, X);

    // Calculate XtCm1X = Xt * Cm1 * X
    MatrixSymmetric* XtCm1X =
      MatrixFactory::prodMatMat<MatrixSymmetric>(&X, &Cm1X, true, false);

    // Construct ZtCm1X = Zt * Cm1 * X and perform its Cholesky decomposition
    VectorDouble ZtCm1X = Cm1X.prodVecMat(Y);
    CholeskyDense XtCm1XChol(XtCm1X);
    if (!XtCm1XChol.isReady())
    {
      messerr("Cholesky decomposition of XtCm1X matrix failed");
      delete XtCm1X;
      return TEST;
    }

    // Calculate beta = (XtCm1X)-1 * ZtCm1X
    VectorDouble beta(nDrift);
    if (XtCm1XChol.solve(ZtCm1X, beta))
    {
      messerr("Error when calculating Likelihood");
      delete XtCm1X;
      return TEST;
    }
    model->setBetaHat(beta);
    delete XtCm1X;

    if (verbose)
    {
      VH::dump("Optimal Drift coefficients = ", beta);
    }

    // Center the data by the optimal drift: Y = Y - beta * X
    VH::subtractInPlace(Y, X.prodMatVec(beta));
  }

  // Calculate t(L-1) %*% D-1 %*% L-1 applied to Y (L and D from Vecchia)
  VectorDouble Cm1Z = _productVecchia(V, Y);

  // Calculate the log-determinant
  double logdet = -VH::cumulLog(V.getDFull());

  // Calculate quad = Zt * Cm1Z
  double quad = VH::innerProduct(Y, Cm1Z);

  // Derive the log-likelihood
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
