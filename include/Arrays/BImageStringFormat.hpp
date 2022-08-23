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
#include "Basic/AStringFormat.hpp"
#include "Basic/Vector.hpp"

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
