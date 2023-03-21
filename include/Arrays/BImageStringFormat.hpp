/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Authors: <authors>                                                         */
/* Website: <website>                                                         */
/* License: BSD 3 clause                                                      */
/*                                                                            */
/******************************************************************************/
#pragma once

#include "gstlearn_export.hpp"
#include "Basic/AStringFormat.hpp"
#include "Basic/VectorNumT.hpp"

#include "geoslib_define.h"

class GSTLEARN_EXPORT BImageStringFormat: public AStringFormat
{
public:
  BImageStringFormat(char zero = '0',
                     char one = '1',
                     const VectorInt &indMin = VectorInt(),
                     const VectorInt indMax = VectorInt());
  BImageStringFormat(const BImageStringFormat& r);
  BImageStringFormat& operator=(const BImageStringFormat& r);
  virtual ~BImageStringFormat();

  char getCharOne() const { return _charOne; }
  char getCharZero() const { return _charZero; }
  const VectorInt getIndMax() const { return _indMax; }
  int getIndMin(int idim) const;
  const VectorInt& getIndMin() const { return _indMin; }
  int getIndMax(int idim) const;

private:
  VectorInt _indMin;
  VectorInt _indMax;
  char _charZero;
  char _charOne;
};
