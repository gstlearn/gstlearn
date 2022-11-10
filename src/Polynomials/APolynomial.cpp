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
#include "Polynomials/APolynomial.hpp"
#include "Basic/AException.hpp"

#include "csparse_d.h"
#include "math.h"

#include <string>
#include <algorithm>
#include <sstream>
#include <iterator>
#include <iostream>

APolynomial::APolynomial()
    : AStringable()
{
}

APolynomial::APolynomial(VectorDouble coeffs)
  : AStringable()
{
  init(coeffs);
}

APolynomial::APolynomial(const APolynomial &m)
    : AStringable(m),
      _coeffs(m._coeffs)
{
}

APolynomial::~APolynomial()
{

}

APolynomial & APolynomial::operator=(const APolynomial& p)
{
  if (this !=& p)
  {
    AStringable::operator=(p);
    _coeffs=p._coeffs;
  }
  return *this;
}

VectorDouble APolynomial::evalOp(cs* Op, const VectorDouble& in) const
{
  VectorDouble result(in.size());
  evalOp(Op,in,result);
  return result;
}

String APolynomial::toString(const AStringFormat* /*strfmt*/) const
{
  String str;
  if (_coeffs.size() <= 0) return str;

  std::ostringstream oss;
  str += "Polynomials of degree " + std::to_string(_coeffs.size()-1) + "\n";
  if (!_coeffs.empty())
  {
    std::copy(_coeffs.begin(), _coeffs.end(),
    std::ostream_iterator<double>(oss, " "));
  }
  str +=oss.str() + "\n";
  return str;
}
void APolynomial::init(VectorDouble coeffs)
{
  _coeffs = coeffs;
}
