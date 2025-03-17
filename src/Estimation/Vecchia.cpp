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
#include "LinearOp/CholeskyDense.hpp"
#include "Matrix/MatrixSparse.hpp"
#include "Model/Model.hpp"

Vecchia::Vecchia(const Model* model,
                 const Db* dbin,
                 const Db* dbout)
  : _dbin(dbin)
  , _dbout(dbout)
  , _model(model)
  , _m()
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
int Vecchia::computeLower(const MatrixT<int>& Ranks)
{
  int ndim     = _model->getNDim();
  int ntot     = Ranks.getNRows();
  int nb_neigh = Ranks.getNCols();
  int Ndbin    = _dbin->getNSample();
  VectorDouble D(ntot);
  double value;

  // Resizing
  _m = MatrixSparse(ntot, ntot, nb_neigh+1);
  _Dmat = MatrixSparse(ntot, ntot);

  // Creating empty Dbs
  Db* Dbtemp = Db::createEmpty(nb_neigh, ndim, 0, 0, 0, false, true);
  Dbtemp->setLocVariable(ELoc::SEL, 0, 0, 0.);
  Db* DbOnePoint = Db::createEmpty(1, ndim, 0);

  // Loop on the samples
  for (int ind = 0; ind < ntot; ind++)
  {
    int icur = 1;
    int iabs = Ranks(ind, 0);
    for (int idim = 0; idim < ndim; idim++)
    {
      if (iabs < Ndbin)
        value = _dbin->getCoordinate(iabs, idim);
      else
        value = _dbout->getCoordinate(iabs - Ndbin, idim);
      DbOnePoint->setCoordinate(0, idim, value);
    }

    for (int jp = nb_neigh-1; jp >= 1; jp--)
    {
      int ip = Ranks(ind,jp);
      if (IFFFF(ip))
      {
        Dbtemp->setLocVariable(ELoc::SEL, icur, 0, 0.);
        continue;
      }
      for (int idim = 0; idim < ndim; idim++)
      {
        if (ip < Ndbin)
          value = _dbin->getCoordinate(ip, idim);
        else
          value = _dbout->getCoordinate(ip - Ndbin, idim);
        Dbtemp->setCoordinate(icur, idim, value);
      }
      Dbtemp->setLocVariable(ELoc::SEL, icur, 0, 1.);
      icur++;
    }

    MatrixSquareSymmetric mat = _model->evalCovMatSym(Dbtemp);
    CholeskyDense chol = CholeskyDense(&mat);
    MatrixRectangular crossmat = _model->evalCovMat(Dbtemp, DbOnePoint);
    VectorDouble vect = crossmat.getColumn(0);
    VectorDouble res = chol.solveX(vect);

    icur = 0;
    _m.setValue(ind, ind, 1.);
    for (int jp = nb_neigh-1; jp >= 1; jp--)
    {
      int ip = Ranks(ind,jp);
      if (IFFFF(ip)) continue;
      _m.setValue(ind, ip, -res[icur]);
      icur++;
    }
    D[ind] = 1. / (1. - VH::innerProduct(res, vect));
  }
  _Dmat.setDiagonal(D);

  return 0;
}