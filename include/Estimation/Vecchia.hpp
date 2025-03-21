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

#include "gstlearn_export.hpp"

#include "Matrix/MatrixT.hpp"
#include "Matrix/MatrixSparse.hpp"
#include "Basic/NamingConvention.hpp"

class Db;
class ModelGeneric;

class GSTLEARN_EXPORT Vecchia
{
public:
  Vecchia(const ModelGeneric* model, const Db* db1, const Db* db2 = nullptr);
  Vecchia(const Vecchia &r) = delete;
  Vecchia& operator=(const Vecchia &r) = delete;
  virtual ~Vecchia();

  int computeLower(const MatrixT<int>& Ranks, bool verbose = false);
  const MatrixSparse& getLFull() const { return _LFull; }
  const VectorDouble& getDFull() const { return _DFull; }
  const MatrixSparse& getDMat() const { return _Dmat; }

  double getLFull(int i, int j) const { return _LFull.getValue(i, j); }

private:
  // Following members are copies of pointers (not to be deleted)
  const Db* _db1;
  const Db* _db2;
  const ModelGeneric* _model;

  mutable VectorDouble _DFull;
  mutable MatrixSparse _LFull;

  // Local calculation results (to be deleted later)
  mutable MatrixSparse _Dmat;
};

GSTLEARN_EXPORT int krigingVecchia(Db* dbin,
                                   Db* dbout,
                                   ModelGeneric* model,
                                   int nb_neigh = 5,
                                   bool verbose = false,
                                   const NamingConvention& namconv = NamingConvention("Vecchia"));

