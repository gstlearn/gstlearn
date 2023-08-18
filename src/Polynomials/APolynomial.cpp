/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://gstlearn.org                                              */
/* License: BSD 3 clauses                                                     */
/*                                                                            */
/******************************************************************************/
#include "Polynomials/APolynomial.hpp"
#include "Basic/AException.hpp"
#include "Matrix/LinkMatrixSparse.hpp"

#include <string>
#include <algorithm>
#include <sstream>
#include <iterator>
#include <iostream>
#include <math.h>

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
#ifndef SWIG
VectorDouble APolynomial::evalOp(cs* Op, const VectorDouble& in) const
{
  VectorDouble result(in.size());
  evalOp(Op,in,result);
  return result;
}
#endif
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
