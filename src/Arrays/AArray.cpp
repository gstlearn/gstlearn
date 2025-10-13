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
#include "Arrays/AArray.hpp"
#include "Basic/VectorNumT.hpp"

namespace gstlrn
{ 
AArray::AArray(const VectorInt& ndims)
    : AStringable(),
      _ndims(ndims)
{
}

AArray::AArray(const AArray &r)
    : AStringable(r),
      _ndims(r._ndims)
{
}

AArray& AArray::operator=(const AArray &r)
{
  if (this != &r)
  {
    AStringable::operator=(r);
    _ndims = r._ndims;
  }
  return *this;
}

AArray::~AArray()
{

}

void AArray::init(const VectorInt& ndims)
{
  _ndims = ndims;
}

String AArray::toString(const AStringFormat* /*strfmt*/) const
{
  std::stringstream sstr;
  if (getNDim() <= 0) return sstr.str();

  sstr << "Array dimension = " << getNDim() << std::endl;

  for (Id idim = 0; idim < getNDim(); idim++)
  {
    sstr << "- Dimension #" << idim+1 << " : " << getNDims(idim) << std::endl;
  }
  return sstr.str();
}


Id AArray::indiceToRank(const VectorInt& indice) const
{
  if (! _isValidIndice(indice)) return ITEST;
  auto ndim = getNDim();
  Id ival = indice[ndim-1];
  if (ival < 0 || ival >= _ndims[ndim-1])
    return(-1);
  for (Id idim=ndim-2; idim>=0; idim--)
  {
    if (indice[idim] < 0 || indice[idim] >= _ndims[idim])
      return(-1);
    ival = ival * _ndims[idim] + indice[idim];
  }
  return ival;
}

void AArray::rankToIndice(Id rank, VectorInt& indices) const
{
  auto ndim = getNDim();
  if (static_cast<Id>(indices.size()) != ndim)
  {
    messerr("Argument indices should have the correct size. Output argument 'indices' not modified");
    return;
  }

  Id nval = 1;
  for (Id idim=0; idim<ndim; idim++) nval *= _ndims[idim];

  for (Id idim=ndim-1; idim>=0; idim--)
  {
    nval /= _ndims[idim];
    indices[idim] = rank / nval;
    rank -= indices[idim] * nval;
  }
}

VectorInt AArray::rankToIndice(Id rank) const
{
  auto ndim = getNDim();
  VectorInt indices(ndim);
  rankToIndice(rank,indices);
  return indices;
}

Id AArray::getNDims(Id idim) const
{
  if (idim < getNDim())
    return _ndims[idim];
  return 1;
}

VectorInt AArray::getNDimsExt(Id ndimMax) const
{
  VectorInt ndims = _ndims;
  ndims.resize(ndimMax, 1);
  return ndims;
}

bool AArray::_isValidIndice(const VectorInt& indice) const
{
  auto ndim = getNDim();
  if (static_cast<Id>(indice.size()) != ndim)
  {
    messerr("Argument 'indice' does not have the correct dimension (%d)",static_cast<Id>(indice.size()));
    messerr("It should match the AArray dimension (%d)", ndim);
    return false;
  }

  for (Id idim = 0; idim < ndim; idim++)
  {
    if (!checkArg("Element of 'indice'", indice[idim], _ndims[idim]))
      return false;
  }
  return true;
}
}