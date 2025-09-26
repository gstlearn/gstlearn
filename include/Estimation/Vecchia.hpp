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

#include "Basic/NamingConvention.hpp"
#include "Basic/VectorNumT.hpp"
#include "Estimation/ALikelihood.hpp"
#include "LinearOp/CholeskyDense.hpp"
#include "Matrix/MatrixDense.hpp"
#include "Matrix/MatrixSparse.hpp"
#include "Matrix/MatrixSymmetric.hpp"
#include "Matrix/MatrixT.hpp"
#include "Mesh/AMesh.hpp"
#include "gstlearn_export.hpp"

namespace gstlrn

{
class Db;
class ModelGeneric;

class GSTLEARN_EXPORT Vecchia: public ALikelihood
{
public:
  Vecchia(ModelGeneric* model,
          Id nb_neigh,
          const Db* db1,
          const Db* db2 = nullptr,
          bool reml     = false);
  Vecchia(const Vecchia& r);
  Vecchia& operator=(const Vecchia& r);
  virtual ~Vecchia();

public:
  static Vecchia* createForOptim(ModelGeneric* model,
                                 const Db* db1,
                                 Id nb_neigh = 30,
                                 bool reml   = false);

  Id computeLower(const MatrixT<Id>& Ranks, bool verbose = false);
  const MatrixSparse& getLFull() const { return _LFull; }
  const VectorDouble& getDFull() const { return _DFull; }
  const VectorDouble& getY() const { return _Y; }

  double getLFull(Id i, Id j) const { return _LFull.getValue(i, j); }
  Id getND() const { return _Ntot2; }
  Id getNT() const { return _Ntot1; }
  Id getNonZeros() const { return _LFull.getNonZeros(); }

  void productMatVecchia(const MatrixDense& X, MatrixDense& resmat) const;
  void productVecchia(constvect Y, vect res) const;
  VectorDouble calculateLdY(const VectorDouble& Y) const;
  VectorDouble calculateFtLdY(const VectorDouble& LdY) const;
  MatrixSparse* calculateW(const VectorDouble& D_dd) const;
  VectorDouble computeAndGetY();
private:
  void _init(bool verbose = false) override;
  void _updateModel(bool verbose = false) override;
  void _computeCm1X() override;
  void _computeCm1Y() override;
  double _computeLogDet() const override;
  Id _buildNeighborhood(const MatrixT<Id>& Ranks,
                        Id isample,
                        Id ivar,
                        Id nb_neigh,
                        std::vector<std::array<Id, 4>>& neighDescr) const;
  void _buildLHS(Id nitems,
                 const std::vector<std::array<Id, 4>>& neighDescr,
                 MatrixSymmetric& matCov);
  void _buildRHS(Id icase2,
                 Id iabs2,
                 Id ivar2,
                 Id nitems,
                 const std::vector<std::array<Id, 4>>& neighDescr,
                 MatrixDense& vectCov);
  void _loadDataFlattened();
  Id _getAddressInMatrix(Id ip, Id ivar) const;
  Id _getAddressAbsolute(Id ip) const;
  Id _getSampleCase(Id ip) const;
  Id _getCase() const;

private:
  // Following members are copies of pointers (not to be deleted)
  Id _nbNeigh;
  const Db* _db1;
  const Db* _db2;

  MatrixT<Id> _Ranks; // Matrix of ranks for the Vecchia approximation
  MatrixSymmetric _matCov;
  MatrixDense _vectCov;
  VectorDouble _work; // Work vector for calculations
  mutable VectorDouble _LdY;
  mutable VectorDouble _DFull;
  mutable MatrixSparse _LFull;
  mutable CholeskyDense* _chol; // Cholesky decomposition of the covariance matrix
  // Local calculation results (to be deleted later)
  mutable Id _Ndb1; // Number of samples in Db1 (used for shift calculations)
  mutable Id _Ntot1;
  mutable Id _Ntot2;
  mutable VectorInt _cumulRanks1;
  mutable VectorInt _cumulRanks2;
  mutable VectorVectorInt _varRanks1;
  mutable VectorVectorInt _varRanks2;
};

GSTLEARN_EXPORT Id krigingVecchia(Db* dbin,
                                  Db* dbout,
                                  ModelGeneric* model,
                                  Id nb_neigh                     = 5,
                                  bool verbose                    = false,
                                  const NamingConvention& namconv = NamingConvention("Vecchia"));
GSTLEARN_EXPORT double logLikelihoodVecchia(const Db* db,
                                            ModelGeneric* model,
                                            Id nb_neigh  = 5,
                                            bool verbose = false);
} // namespace gstlrn