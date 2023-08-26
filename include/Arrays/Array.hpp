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
#include "Basic/AStringable.hpp"

class GSTLEARN_EXPORT Array : public AArray
{
public:
  Array(const VectorInt& ndims = VectorInt());
  Array(const Array &m);
  Array& operator=(const Array &m);
  virtual ~Array();

  void init(const VectorInt& ndims);
  double getValue(const VectorInt& indice) const;
  void setValue(const VectorInt& indice, double value);
  const VectorDouble& getValues() const { return _values; }
  void setValues(const VectorDouble& values) { _values = values; }

private:
  void _update();

private:
  VectorDouble _values;
};
