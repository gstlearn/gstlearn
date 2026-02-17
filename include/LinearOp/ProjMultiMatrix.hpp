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

#include "LinearOp/ProjMatrix.hpp"
#include "LinearOp/ProjMulti.hpp"
#include "Matrix/MatrixSparse.hpp"

namespace gstlrn
{
class AMesh;
class Db;
class VectorMeshes;

class GSTLEARN_EXPORT ProjMultiMatrix: public ProjMulti
{
  using ProjsStore = std::vector<std::vector<std::optional<ProjMatrix>>>;

public:
  ProjMultiMatrix(const std::vector<std::vector<const ProjMatrix*>>& proj,
                  bool silent = false);
#ifndef SWIG
  ProjMultiMatrix(ProjsStore&& proj,
                  bool silent = false);
#endif // SWIG

  ~ProjMultiMatrix() override = default;
  static std::vector<std::vector<const ProjMatrix*>> create(std::vector<const ProjMatrix*>& vectproj,
                                                            Id nvariable);
  static ProjMultiMatrix* createFromDbAndMeshes(const Db* db,
                                                const VectorMeshes& meshes,
                                                Id ncov,
                                                Id nvar,
                                                bool checkOnZVariable = true,
                                                bool verbose          = false);

  const MatrixSparse* getProj() const { return &_Proj; }
#ifndef SWIG

protected:
  void init();

  Id _addPoint2mesh(const constvect inv, vect outv) const override;
  Id _addMesh2point(const constvect inv, vect outv) const override;
  ProjsStore _projsStore;
#endif

private:
  MatrixSparse _Proj;
};
} // namespace gstlrn
