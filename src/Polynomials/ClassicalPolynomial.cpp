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

#include "Polynomials/ClassicalPolynomial.hpp"

#include "Basic/Vector.hpp"

#include "csparse_d.h"
#include "csparse_f.h"

ClassicalPolynomial::ClassicalPolynomial()
{
  // TODO Auto-generated constructor stub

}


ClassicalPolynomial::ClassicalPolynomial(const VectorDouble& coeffs)
{
  init(coeffs);
}

ClassicalPolynomial::~ClassicalPolynomial()
{
  // TODO Auto-generated destructor stub
}

double ClassicalPolynomial::eval(double x) const
{
  double result=_coeffs.back();
  for (int i = _coeffs.size()-2; i >= 0; i--)
  {
    result *= x;
    result +=_coeffs[i];
  }
  return result;
}


void ClassicalPolynomial::evalOp(cs* Op, const VectorDouble& in, VectorDouble& out) const
{
  int n = in.size();
  VectorDouble work(n);

  for(int i = 0; i < n ;i++)
  {
    out[i] = _coeffs.back() * in[i];
  }

  for(int j=_coeffs.size()-2;j>=0;j--)
  {
    cs_vecmult(Op,out.data(),work.data());
    for (int i = 0; i<n ; i++)
    {
        out[i] = _coeffs[j] * in[i] + work[i];
    }
  }
}
