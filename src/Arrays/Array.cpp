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
#include "Arrays/Array.hpp"
#include "Basic/VectorNumT.hpp"

namespace gstlrn
{
Array::Array(const VectorInt& ndims)
    : AArray(ndims),
      _values()
{
  _update();
}

Array::Array(const Array &r)
    : AArray(r),
      _values(r._values)
{

}

Array& Array::operator=(const Array &r)
{
  if (this != &r)
  {
    AArray::operator=(r);
    _values = r._values;
  }
  return *this;
}

Array::~Array()
{

}

void Array::init(const VectorInt& ndims)
{
  AArray::init(ndims);
  _update();
}

void Array::_update()
{
  auto total = getNPixels();
  _values.resize(total,0.);
}

double Array::getValue(const VectorInt& indice) const
{
  if (! _isValidIndice(indice)) return TEST;
  Id iad = indiceToRank(indice);
  return _values[iad];
}

void Array::setValue(const VectorInt& indice, double value)
{
  if (! _isValidIndice(indice)) return;
  Id iad = indiceToRank(indice);
  _values[iad] = value;
}
}