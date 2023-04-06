/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://github.com/gstlearn                                       */
/* License: BSD 3 clauses                                                     */
/*                                                                            */
/******************************************************************************/
#include "Arrays/BImageStringFormat.hpp"

BImageStringFormat::BImageStringFormat(char zero,
                                       char one,
                                       const VectorInt &indMin,
                                       const VectorInt indMax)
    : AStringFormat(1),
      _indMin(indMin),
      _indMax(indMax),
      _charZero(zero),
      _charOne(one)
{
}

BImageStringFormat::BImageStringFormat(const BImageStringFormat& r)
    : AStringFormat(r),
      _indMin(r._indMin),
      _indMax(r._indMax),
      _charZero(r._charZero),
      _charOne(r._charOne)
{
}

BImageStringFormat& BImageStringFormat::operator=(const BImageStringFormat& r)
{
  if (this != &r)
  {
    AStringFormat::operator=(r);
    _indMin = r._indMin;
    _indMax = r._indMax;
    _charZero = r._charZero;
    _charOne = r._charOne;
  }
  return *this;
}

BImageStringFormat::~BImageStringFormat()
{
}

int BImageStringFormat::getIndMin(int idim) const
{
  if (idim < (int) _indMin.size())
    return _indMin[idim];
  else
    return 0;
}

int BImageStringFormat::getIndMax(int idim) const
{
  if (idim < (int) _indMax.size())
    return _indMax[idim];
  else
    return ITEST;
}
