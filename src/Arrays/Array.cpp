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
#include "Arrays/Array.hpp"
#include "Basic/Vector.hpp"
#include "geoslib_f.h"

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

String Array::toString(const AStringFormat* /*strfmt*/) const
{
  std::stringstream sstr;
  if (getNDim() <= 0) return sstr.str();

  sstr << "Array dimension = " << getNDim() << std::endl;

  for (int idim = 0; idim < getNDim(); idim++)
  {
    sstr << "- Dimension #" << idim+1 << " : " << getNDims(idim) << std::endl;
  }
  return sstr.str();
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
