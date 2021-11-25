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
#include "Polynomials/APolynomial.hpp"
#include "Basic/IClonable.hpp"
#include "Basic/Vector.hpp"

class ShiftOpCs;
class GSTLEARN_EXPORT ClassicalPolynomial : public APolynomial
{
public:
  ClassicalPolynomial();
  ClassicalPolynomial(const VectorDouble&);
  virtual ~ClassicalPolynomial();
  virtual IClonable* clone() const override { return new ClassicalPolynomial(*this); }
  double eval(double x) const override;
  void evalDerivOp(ShiftOpCs* shiftOp,const VectorDouble& in,
                   VectorDouble& out,int iapex,int igparam)const;
  void evalOpTraining(cs* Op, const VectorDouble& in,VectorVectorDouble& store,VectorDouble& work) const override;
  void evalDerivOpOptim(ShiftOpCs* shiftOp,VectorDouble& temp1,VectorDouble& temp2,
                       VectorDouble& out,const VectorVectorDouble workpoly,int iapex,int igparam)const;
  void evalOpCumul(cs* Op, const VectorDouble& in, VectorDouble& out) const ;
  void evalOp(cs* Op, const VectorDouble& in, VectorDouble& out) const override;
};
