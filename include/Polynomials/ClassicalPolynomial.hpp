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

#include "Polynomials/APolynomial.hpp"
#include "Basic/Vector.hpp"

class ClassicalPolynomial : public APolynomial
{
public:
  ClassicalPolynomial();
  ClassicalPolynomial(const VectorDouble&);
  virtual ~ClassicalPolynomial();

  double eval(double x) const override;
  void evalOp(cs* Op, const VectorDouble& in, VectorDouble& out) const override;
};
