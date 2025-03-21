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
#include "Matrix/MatrixSquareSymmetric.hpp"
#include "Matrix/MatrixT.hpp"
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
  int ntot     = Ranks.getNRows();
  int nb_neigh = Ranks.getNCols() - 1;
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
      MatrixSquareSymmetric mat = _model->evalCovMatSym(Dbtemp);
      CholeskyDense chol         = CholeskyDense(&mat);
      MatrixRectangular crossmat = _model->evalCovMat(Dbtemp, DbOnePoint);
      VectorDouble vect = crossmat.getColumn(0);
      VectorDouble res  = chol.solveX(vect);

      icur = 0;
      _LFull.setValue(ind, ind, 1.);
      for (int jp = nb_neigh; jp >= 1; jp--)
      {
        int ip = Ranks(ind, jp);
        if (IFFFF(ip)) continue;
        _LFull.setValue(ind, ip, -res[icur]);
        icur++;
      }
      _DFull[ind] = 1. / (varK - VH::innerProduct(res, vect));
    }
  }
  _Dmat.setDiagonal(_DFull);

  // Optional printout
  if (verbose)
  {
    message("Matrix L\n");
    _LFull.display();
    VH::dump("Diagonal D", _DFull);
  }
  return 0;
}

int krigingVecchia(Db* dbin,
                   Db* dbout,
                   ModelGeneric* model,
                   int nb_neigh,
                   bool verbose,
                   const NamingConvention& namconv) 
{
  MatrixT<int> Ranks = findNN(dbout, dbin, nb_neigh+1, false, verbose);

  Vecchia V = Vecchia(model, dbout, dbin);

  if (V.computeLower(Ranks, verbose)) return 1;

  // Main dimensions
  int nd         = dbin->getNSample();
  int nt         = dbout->getNSample();
  VectorDouble Y = dbin->getColumnByLocator(ELoc::Z, 0);

  // Extract sub-part of 'Diagonal' vector
  VectorDouble DFull = V.getDFull();
  VectorDouble D_tt(nt);
  VectorDouble D_dd(nd);
  VectorDouble D_dd3(nd);

  VH::extractInPlace(DFull, D_tt, 0);
  VH::extractInPlace(DFull, D_dd3, nt);
  VH::extractInPlace(DFull, D_dd, nt);
  VH::transformVD(D_dd3, 3);

  // Calculate LdY = Ldat %*% Y
  VectorDouble LdY(nd);
  for (int id = 0; id < nd; id++)
  {
    double value = 0.;
    for (int jd = 0; jd < nd; jd++)
      value += V.getLFull(id + nt, jd + nt) * Y[jd];
    LdY[id] = value;
  }
  VH::multiplyInPlace(LdY, D_dd);

  // Calculate FtLdY = Ft %*% Ldat %*% Y
  VectorDouble FtLdY(nt);
  for (int it = 0; it < nt; it++)
  {
    double value = 0.;
    for (int id = 0; id < nd; id++)
      value += V.getLFull(id + nt, it) * LdY[id];
    FtLdY[it] = value;
  }
  LdY.clear();

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
  MatrixSparse* W = MatrixSparse::addMatMat(mat1, mat2);
  CholeskySparse cholW(W);
  VectorDouble result = cholW.solveX(FtLdY);
  for (int i = 0; i < nt; i++) result[i] = -result[i];

  // Cleaning stage
  indT.clear();
  indD.clear();
  FtLdY.clear();
  delete Ltt;
  delete Ldt;
  delete mat1;
  delete mat2;
  delete W;

  int iptr = dbout->addColumns(result, String(), ELoc::UNKNOWN, 0, true);
  namconv.setNamesAndLocators(dbout, iptr, "estim", 1);

  return 0;
}
