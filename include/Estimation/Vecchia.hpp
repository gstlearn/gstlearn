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

class Db;
class Model;

class GSTLEARN_EXPORT Vecchia
{
public:
  Vecchia(const Db* dbin, const Db* dbout, const Model* model);
  Vecchia(const Vecchia &r) = delete;
  Vecchia& operator=(const Vecchia &r) = delete;
  virtual ~Vecchia();

  int computeLower(const MatrixT<int>& Ranks);
  const MatrixSparse& getM() const { return _m; }
  const MatrixSparse& getDMat() const { return _Dmat; }

private:
  // Following members are copies of pointers (not to be deleted)
  const Db* _dbin;
  const Db* _dbout;
  const Model* _model;

  mutable MatrixSparse _m;
  mutable MatrixSparse _Dmat;
};
