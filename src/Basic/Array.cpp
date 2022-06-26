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
#include "Basic/Array.hpp"
#include "Basic/Vector.hpp"
#include "geoslib_f.h"

Array::Array(const VectorInt& ndims)
    : AStringable(),
      _ndims(ndims),
      _values()
{
  _update();
}

Array::Array(const Array &r)
    : AStringable(r),
      _ndims(r._ndims),
      _values(r._values)
{

}

Array& Array::operator=(const Array &r)
{
  if (this != &r)
  {
    AStringable::operator=(r);
    _ndims = r._ndims;
    _values = r._values;
  }
  return *this;
}

Array::~Array()
{

}

void Array::init(const VectorInt& ndims)
{
  _ndims = ndims;
  _update();
}

void Array::_update()
{
  int total = ut_vector_prod(_ndims);
  _values.resize(total,0.);
}

String Array::toString(const AStringFormat* /*strfmt*/) const
{
  std::stringstream sstr;
  if (_ndims.size() <= 0) return sstr.str();

  sstr << "Array dimension = " << (int) _ndims.size() << std::endl;

  for (int idim = 0; idim < (int) _ndims.size(); idim++)
  {
    sstr << "- Dimension #" << idim+1 << " : " << _ndims[idim] << std::endl;
  }
  return sstr.str();
}

int Array::indiceToRank(const VectorInt& indice) const
{
  if (! _isValidIndice(indice)) return ITEST;
  int ndim = (int) _ndims.size();
  int ival = indice[ndim-1];
  if (ival < 0 || ival >= _ndims[ndim-1])
    return(-1);
  for (int idim=ndim-2; idim>=0; idim--)
  {
    if (indice[idim] < 0 || indice[idim] >= _ndims[idim])
      return(-1);
    ival = ival * _ndims[idim] + indice[idim];
  }
  return ival;
}

void Array::rankToIndice(int rank, VectorInt& indices) const
{
  int ndim = (int) _ndims.size();

  int nval = 1;
  for (int idim=0; idim<ndim; idim++) nval *= _ndims[idim];

  for (int idim=ndim-1; idim>=0; idim--)
  {
    nval /= _ndims[idim];
    indices[idim] = rank / nval;
    rank -= indices[idim] * nval;
  }
}

VectorInt Array::rankToIndice(int rank) const
{
  int ndim = (int) _ndims.size();
  VectorInt indices(ndim);
  rankToIndice(rank,indices);
  return indices;
}

bool Array::_isValidIndice(const VectorInt& indice) const
{
  int ndim = (int) _ndims.size();
  if ((int) indice.size() != ndim)
  {
    messerr("Argument 'indice' does not have the correct dimension (%d)",(int) indice.size());
    messerr("It should match the Array dimension (%d)", ndim);
    return false;
  }

  for (int idim = 0; idim < ndim; idim++)
  {
    if (indice[idim] < 0 || indice[idim] >= _ndims[idim])
    {
      mesArg("Element of 'indice'",indice[idim],_ndims[idim]);
      return false;
    }
  }
  return true;
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
