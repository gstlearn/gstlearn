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
#pragma once

#include "Basic/VectorNumT.hpp"
#include "LinearOp/CholeskyDense.hpp"
#include "Matrix/MatrixDense.hpp"
#include "Matrix/MatrixSymmetric.hpp"
#include "Mesh/AMesh.hpp"
#include "Model/AModelOptimNew.hpp"
#include "gstlearn_export.hpp"

#include "Matrix/MatrixT.hpp"
#include "Matrix/MatrixSparse.hpp"
#include "Basic/NamingConvention.hpp"

class Db;
class ModelGeneric;

class GSTLEARN_EXPORT Vecchia: public AModelOptimNew
{
public:
  Vecchia(const ModelGeneric* model,
          const Db* db1,
          const Db* db2 = nullptr);
  Vecchia(const Vecchia& r)            = delete;
  Vecchia& operator=(const Vecchia& r) = delete;
  virtual ~Vecchia();

  static Vecchia* createForOptim(ModelGeneric* model,
                                 const Db* db1,
                                 int nb_neigh = 30,
                                 bool verbose = false);

  double computeCost(bool verbose = false) override;
  int computeLower(const MatrixT<int>& Ranks, bool verbose = false);
  const MatrixSparse& getLFull() const { return _LFull; }
  const VectorDouble& getDFull() const { return _DFull; }
  const MatrixSparse& getDMat() const { return _Dmat; }

  double getLFull(int i, int j) const { return _LFull.getValue(i, j); }
  int getND() const { return _nd; }
  int getNT() const { return _nt; }

  void productMatVecchia(const MatrixDense& X, MatrixDense& resmat) const;
  void productVecchia(constvect Y, vect res) const;
  VectorDouble calculateLdY(const VectorDouble& Y) const;
  VectorDouble calculateFtLdY(const VectorDouble& LdY) const;
  MatrixSparse* calculateW(const VectorDouble& D_dd) const;


private:
  // Following members are copies of pointers (not to be deleted)
  const Db* _db1;
  const Db* _db2;

  int _nt;
  int _nd;
  VectorDouble _Y; // Vector of multivariate data
  MatrixT<int> _Ranks; // Matrix of ranks for the Vecchia approximation
  MatrixDense _X; // Matrix of drifts
  VectorDouble _beta;
  MatrixDense _Cm1X;
  VectorDouble _Cm1Z;
  MatrixSymmetric _matCov; 
  MatrixDense _vectCov;
  VectorDouble _work; // Work vector for calculations
  mutable VectorDouble _LdY;
  mutable VectorDouble _DFull;
  mutable MatrixSparse _LFull;
  mutable Db* _dbTemp; // Temporary Db for calculations
  mutable Db* _dbOnePoint; // Temporary Db for one point calculations
  mutable CholeskyDense* _chol; // Cholesky decomposition of the covariance matrix
  // Local calculation results (to be deleted later)
  mutable MatrixSparse _Dmat;
};

GSTLEARN_EXPORT int krigingVecchia(Db* dbin,
                                   Db* dbout,
                                   ModelGeneric* model,
                                   int nb_neigh                    = 5,
                                   bool verbose                    = false,
                                   const NamingConvention& namconv = NamingConvention("Vecchia"));
GSTLEARN_EXPORT double logLikelihoodVecchia(const Db* db,
                                            ModelGeneric* model,
                                            int nb_neigh = 5,
                                            bool verbose = false);
