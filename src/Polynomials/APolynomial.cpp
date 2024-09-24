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
#include "Polynomials/APolynomial.hpp"
#include "Basic/VectorNumT.hpp"

#include "Matrix/MatrixSparse.hpp"
#include <string>
#include <algorithm>
#include <sstream>
#include <iterator>
#include <math.h>

APolynomial::APolynomial()
    : AStringable()
{
}

APolynomial::APolynomial(const VectorDouble& coeffs)
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
VectorDouble APolynomial::evalOp(MatrixSparse* Op, const constvect& in) const
{
  VectorDouble result(in.size());
  vect results(result);
  evalOp(Op,in,results);
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
void APolynomial::init(const VectorDouble& coeffs)
{
  _coeffs = coeffs;
}
