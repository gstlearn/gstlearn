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
#pragma once

#include "gstlearn_export.hpp"

#include "Polynomials/APolynomial.hpp"
#include "Basic/ICloneable.hpp"
#include "Basic/VectorNumT.hpp"
#include <Eigen/src/Core/Matrix.h>

#ifndef SWIG
#  include <Eigen/Core>
#  include <Eigen/Dense>
#endif

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
  static void evalDerivOp(ShiftOpCs* shiftOp,
                          const VectorDouble& inv,
                          VectorDouble& outv,
                          int iapex,
                          int igparam);
  static void evalDerivOpOptim(ShiftOpCs* shiftOp,
                               VectorDouble& temp1,
                               VectorDouble& temp2,
                               VectorDouble& outv,
                               const VectorVectorDouble& workpoly,
                               int iapex,
                               int igparam);
#ifndef SWIG
  void evalDerivOp(ShiftOpCs* shiftOp,const Eigen::VectorXd& inv,
                   Eigen::VectorXd& outv,int iapex,int igparam)const;
 
  void evalDerivOpOptim(ShiftOpCs* shiftOp,Eigen::VectorXd& temp1,Eigen::VectorXd& temp2,
                       Eigen::VectorXd& outv,const std::vector<Eigen::VectorXd>& workpoly,int iapex,int igparam)const;
  // void evalOp(const ALinearOpMulti* /*Op*/,
  //            const std::vector<Eigen::VectorXd>& /*inv*/,
  //           std::vector<Eigen::VectorXd>& /*outv*/) const override { } 
              

  void evalOpTraining(MatrixSparse *Op,
                      const Eigen::VectorXd &inv,
                      std::vector<Eigen::VectorXd> &store,
                      Eigen::VectorXd &work) const override;
  void evalOpCumul(MatrixSparse *Op,
                   const Eigen::VectorXd &inv,
                   Eigen::VectorXd &outv) const;
  void evalOp(MatrixSparse* Op, const Eigen::VectorXd& inv, Eigen::VectorXd& outv) const override;
  void addEvalOp(ALinearOp* Op,const Eigen::VectorXd& inv, Eigen::VectorXd& outv) const override;
#endif
  
#ifndef SWIG
  
private:
  mutable Eigen::VectorXd _work;
  mutable Eigen::VectorXd _work2;
#endif
};
