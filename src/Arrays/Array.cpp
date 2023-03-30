/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://github.com/gstlearn                                       */
/* License: BSD 3 clause                                                      */
/*                                                                            */
/******************************************************************************/
#include "Arrays/Array.hpp"
#include "Basic/VectorNumT.hpp"

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
  int total = getNPixels();
  _values.resize(total,0.);
}

double Array::getValue(const VectorInt& indice) const
{
  if (! _isValidIndice(indice)) return TEST;
  int iad = indiceToRank(indice);
  return _values[iad];
}

void Array::setValue(const VectorInt& indice, double value)
{
  if (! _isValidIndice(indice)) return;
  int iad = indiceToRank(indice);
  _values[iad] = value;
}
