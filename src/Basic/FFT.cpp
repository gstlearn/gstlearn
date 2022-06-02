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
#include "geoslib_old_f.h"
#include "Basic/AException.hpp"
#include "Basic/Vector.hpp"

#include <math.h>

/****************************************************************************/
/*!
 **  Calculate the FFT in a space of dimension N
 **
 ** \return  Error return code
 **
 *****************************************************************************/
int FFTn(int ndim,
         const VectorInt& dims,
         VectorDouble& Re,
         VectorDouble& Im,
         int iSign,
         double scaling)
{
  int n = (int) Re.size();
  Im.resize(n,0.);
  return fftn(ndim, dims.data(), Re.data(), Im.data(), iSign, scaling);
}
