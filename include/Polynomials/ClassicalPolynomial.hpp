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

#include "Polynomials/APolynomial.hpp"
#include "Basic/ICloneable.hpp"
#include "Basic/VectorNumT.hpp"

class ShiftOpCs;
class ALinearOp;

class GSTLEARN_EXPORT ClassicalPolynomial : public APolynomial
{
public:
  ClassicalPolynomial();
  ClassicalPolynomial(const VectorDouble&);
  virtual ~ClassicalPolynomial();

  /// ICloneable interface
  IMPLEMENT_CLONING(ClassicalPolynomial)

  double eval(double x) const override;
  void evalDerivOp(ShiftOpCs* shiftOp,const VectorDouble& inv,
                   VectorDouble& outv,int iapex,int igparam)const;
  void evalOpTraining(cs* Op, const VectorDouble& inv,VectorVectorDouble& store,VectorDouble& work) const override;
  void evalDerivOpOptim(ShiftOpCs* shiftOp,VectorDouble& temp1,VectorDouble& temp2,
                       VectorDouble& outv,const VectorVectorDouble workpoly,int iapex,int igparam)const;
  void evalOpCumul(cs* Op, const VectorDouble& inv, VectorDouble& outv) const ;
  void evalOp(const ALinearOpMulti* /*Op*/,
              const VectorVectorDouble& /*inv*/,
              VectorVectorDouble& /*outv*/) const override { }
  void evalOp(cs* Op, const VectorDouble& inv, VectorDouble& outv) const override;
};
