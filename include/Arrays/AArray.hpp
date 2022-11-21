/******************************************************************************/
/* COPYRIGHT ARMINES, ALL RIGHTS RESERVED                                     */
/*                                                                            */
/* THE CONTENT OF THIS WORK CONTAINS CONFIDENTIAL AND PROPRIETARY             */
/* INFORMATION OF ARMINES. ANY DUPLICATION, MODIFICATION,                     */
/* DISTRIBUTION, OR DISCLOSURE IN ANY FORM, IN WHOLE, OR IN PART, IS STRICTLY */
/* PROHIBITED WITHOUT THE PRIOR EXPRESS WRITTEN PERMISSION OF ARMINES         */
/*                                                                            */
/* TAG_SOURCE_CG                                                              */
/******************************************************************************/
#pragma once

#include "gstlearn_export.hpp"

#include "Basic/VectorHelper.hpp"
#include "Basic/AStringable.hpp"

class GSTLEARN_EXPORT AArray : public AStringable
{
public:
  AArray(const VectorInt& ndims = VectorInt());
  AArray(const AArray &m);
  AArray& operator=(const AArray &m);
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
