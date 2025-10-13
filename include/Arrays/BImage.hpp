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

#include "Arrays/AArray.hpp"
#include "gstlearn_export.hpp"

namespace gstlrn
{
class GSTLEARN_EXPORT BImage: public AArray
{
public:
  BImage(const VectorInt& ndims = VectorInt());
  BImage(const BImage& r);
  BImage& operator=(const BImage& r);
  virtual ~BImage();

  /// Interface for AStringable
  String toString(const AStringFormat* strfmt = nullptr) const override;

  void init(const VectorInt& ndims);
  const VectorUChar& getValues() const { return _values; }
  void setValues(const VectorUChar& values) { _values = values; }

  unsigned char getBImage(Id i, Id j, Id k) const { return _values[_divide(i, j, k)]; }
  unsigned char getOffset(Id i, Id j, Id k) const;
  unsigned char getMaskoff(Id i, Id j, Id k) const;

  unsigned char getValue(Id i) const { return _values[i]; }
  void setValue(Id i, unsigned char c) { _values[i] = c; }

  bool getValue(Id i, Id j, Id k) const;
  void setMaskoff(Id i, Id j, Id k);
  void setOffset(Id i, Id j, Id k);

  Id getAllocSize() const;
  bool isInside(Id i, Id j, Id k) const;
  Id getAddress(Id i, Id j, Id k) const;

private:
  void _update();
  Id _divide(Id i, Id j, Id k) const { return getAddress(i, j, k) / 8; }
  Id _residu(Id i, Id j, Id k) const { return getAddress(i, j, k) % 8; }

private:
  VectorUChar _values;
};
} // namespace gstlrn