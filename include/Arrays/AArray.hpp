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

#include "Basic/VectorHelper.hpp"
#include "Basic/AStringable.hpp"

class GSTLEARN_EXPORT AArray : public AStringable
{
public:
  AArray(const VectorInt& ndims = VectorInt());
  AArray(const AArray &r);
  AArray& operator=(const AArray &r);
  virtual ~AArray();

  /// Interface for AStringable
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  void init(const VectorInt& ndims);
  int  indiceToRank(const VectorInt& indice) const;
  VectorInt rankToIndice(int rank) const;
  void rankToIndice(int rank, VectorInt& indices) const;

  int getNDim() const { return (int) _ndims.size(); }
  int getNPixels() const { return  VH::product(_ndims); }
  const VectorInt& getNDims() const { return _ndims; }
  VectorInt getNDimsExt(int ndimMax) const;
  int getNDims(int idim) const;

protected:
  bool _isValidIndice(const VectorInt& indice) const;

private:
  VectorInt _ndims;
};
