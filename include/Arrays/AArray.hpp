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

#include "Basic/AStringable.hpp"
#include "Basic/VectorHelper.hpp"

namespace gstlrn
{
class GSTLEARN_EXPORT AArray: public AStringable
{
public:
  AArray(const VectorInt& ndims = VectorInt());
  AArray(const AArray& r);
  AArray& operator=(const AArray& r);
  virtual ~AArray();

  /// Interface for AStringable
  String toString(const AStringFormat* strfmt = nullptr) const override;

  void init(const VectorInt& ndims);
  Id indiceToRank(const VectorInt& indice) const;
  VectorInt rankToIndice(Id rank) const;
  void rankToIndice(Id rank, VectorInt& indices) const;

  Id getNDim() const { return static_cast<Id>(_ndims.size()); }
  Id getNPixels() const { return VH::product(_ndims); }
  const VectorInt& getNDims() const { return _ndims; }
  VectorInt getNDimsExt(Id ndimMax) const;
  Id getNDims(Id idim) const;

protected:
  bool _isValidIndice(const VectorInt& indice) const;

private:
  VectorInt _ndims;
};
} // namespace gstlrn