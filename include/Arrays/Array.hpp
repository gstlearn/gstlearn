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
#include "Arrays/AArray.hpp"
#include "Basic/Vector.hpp"
#include "Basic/AStringable.hpp"

class GSTLEARN_EXPORT Array : public AArray
{
public:
  Array(const VectorInt& ndims = VectorInt());
  Array(const Array &m);
  Array& operator=(const Array &m);
  virtual ~Array();

  /// Interface for AStringable
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

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
