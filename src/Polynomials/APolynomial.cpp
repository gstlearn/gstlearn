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
#include "math.h"
#include "Polynomials/APolynomial.hpp"
#include "Basic/AException.hpp"
#include "Basic/Vector.hpp"
#include "csparse_d.h"

VectorDouble APolynomial::evalOp(cs* Op, const VectorDouble& in) const
{
  VectorDouble result(in.size());
  evalOp(Op,in,result);
  return result;
}

APolynomial::APolynomial(VectorDouble coeffs)
{
  init(coeffs);
}

void APolynomial::init(VectorDouble coeffs)
{
  _coeffs = coeffs;
}
