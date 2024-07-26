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
#include "Arrays/AArray.hpp"

class GSTLEARN_EXPORT BImage : public AArray
{
public:
  BImage(const VectorInt& ndims = VectorInt());
  BImage(const BImage &m);
  BImage& operator=(const BImage &m);
  virtual ~BImage();

  /// Interface for AStringable
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  void init(const VectorInt& ndims);
  const VectorUChar& getValues() const { return _values; }
  void setValues(const VectorUChar& values) { _values = values; }

  unsigned char getBImage (int i, int j, int k) const { return _values[_divide(i,j,k)]; }
  unsigned char getOffset (int i, int j, int k) const;
  unsigned char getMaskoff(int i, int j, int k) const;

  unsigned char getValue(int i) const { return _values[i]; }
  void setValue(int i, unsigned char c) { _values[i] = c; }

  bool getValue(int i, int j, int k) const;
  void setMaskoff(int i, int j, int k);
  void setOffset(int i, int j, int k);

  int getAllocSize() const;
  bool isInside(int i, int j, int k) const;
  int getAddress(int i, int j, int k) const;

private:
  void _update();
  int _divide(int i, int j, int k) const { return getAddress(i,j,k) / 8; }
  int _residu(int i, int j, int k) const { return getAddress(i,j,k) % 8; }

private:
  VectorUChar _values;
};
