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

#include "LinearOp/ProjMulti.hpp"
#include "Matrix/MatrixSparse.hpp"

namespace gstlrn
{
class ProjMatrix;
class AMesh;
class Db;

class GSTLEARN_EXPORT ProjMultiMatrix: public ProjMulti
{
public:
  ProjMultiMatrix(const std::vector<std::vector<const ProjMatrix*>>& proj,
                  bool toClean = false,
                  bool silent  = false);
  virtual ~ProjMultiMatrix();
  static std::vector<std::vector<const ProjMatrix*>> create(std::vector<const ProjMatrix*>& vectproj,
                                                            Id nvariable);
  static ProjMultiMatrix* createFromDbAndMeshes(const Db* db,
                                                const std::vector<const AMesh*>& meshes,
                                                Id ncov,
                                                Id nvar,
                                                bool checkOnZVariable = true,
                                                bool verbose          = false);

  const MatrixSparse* getProj() const { return &_Proj; }
#ifndef SWIG

protected:
  Id _addPoint2mesh(const constvect inv, vect outv) const override;
  Id _addMesh2point(const constvect inv, vect outv) const override;
#endif

private:
  MatrixSparse _Proj;
  void _clear() override;

private:
  bool _toClean;
};
} // namespace gstlrn