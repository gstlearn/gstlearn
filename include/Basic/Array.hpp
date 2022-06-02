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
#include "Basic/Vector.hpp"
#include "Basic/AStringable.hpp"

class GSTLEARN_EXPORT Array : public AStringable
{
public:
  Array(const VectorInt& ndims = VectorInt());
  Array(const Array &m);
  Array& operator=(const Array &m);
  virtual ~Array();

  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  void init(const VectorInt& ndims);
  int  indiceToRank(const VectorInt& indice) const;
  double getValue(const VectorInt& indice) const;
  void setValue(const VectorInt& indice, double value);

  const VectorInt& getNdims() const { return _ndims; }
  const VectorDouble& getValues() const { return _values; }

private:
  void _update();
  bool _isValidIndice(const VectorInt& indice) const;

private:
  VectorInt _ndims;
  VectorDouble _values;
};
