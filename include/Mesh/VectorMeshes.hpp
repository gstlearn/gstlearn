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

#include "Basic/AStringable.hpp"
#include "gstlearn_export.hpp"

namespace gstlrn
{
/**
 * This corresponds to a list of AMeshes
 */

class AMesh;

class GSTLEARN_EXPORT VectorMeshes: public AStringable
{
public:
  VectorMeshes(const std::vector<const AMesh*>& meshes = {});
  VectorMeshes(Id nmesh);
  VectorMeshes(const VectorMeshes& r);
  VectorMeshes& operator=(const VectorMeshes& r);
  virtual ~VectorMeshes();

  /// Interface to AStringable
  String toString(const AStringFormat* strfmt = nullptr) const override;

  Id size() const { return static_cast<Id>(_meshes.size()); }
  bool empty() const { return _meshes.empty(); }
  bool isTurbo() const;
  bool allDefined() const;
  void replace(Id ind, const AMesh* mesh);

  // Operator overload
  const AMesh* operator()(Id ind) const
  {
    return _meshes[ind];
  }

private:
  std::vector<const AMesh*> _meshes;
};

} // namespace gstlrn
