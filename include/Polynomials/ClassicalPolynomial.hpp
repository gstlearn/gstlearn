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

#include "LinearOp/ALinearOp.hpp"

#include "Polynomials/APolynomial.hpp"
#include "Basic/ICloneable.hpp"
#include "Basic/VectorNumT.hpp"

#include <vector>

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
  // void evalDerivOp(ShiftOpCs* shiftOp,
  //                         const constvect& inv,
  //                         vect& outv,
  //                         int iapex,
  //                         int igparam);
  // static void evalDerivOpOptim(ShiftOpCs* shiftOp,
  //                              vect& temp1,
  //                              vect& temp2,
  //                              vect& outv,
  //                              const VectorVectorDouble& workpoly,
  //                              int iapex,
  //                              int igparam);
#ifndef SWIG
  // void evalDerivOp(ShiftOpCs* shiftOp,const std::vector<double>& inv,
  //                  std::vector<double>& outv,int iapex,int igparam)const;
  
  // void evalDerivOpOptim(ShiftOpCs* shiftOp,Eigen::VectorXd& temp1,Eigen::VectorXd& temp2,
  //                      Eigen::VectorXd& outv,const std::vector<Eigen::VectorXd>& workpoly,int iapex,int igparam)const;
  // void evalOp(const ALinearOpMulti* /*Op*/,
  //            const std::vector<Eigen::VectorXd>& /*inv*/,
  //           std::vector<Eigen::VectorXd>& /*outv*/) const override { }

  void evalOpTraining(MatrixSparse* Op,
                      const constvect inv,
                      std::vector<std::vector<double>>& store,
                      std::vector<double>& work) const override;
  void evalOpCumul(MatrixSparse* Op, const constvect inv, vect outv) const;
  void evalOp(MatrixSparse* Op, const constvect inv, vect outv) const override;
  void addEvalOp(ALinearOp* Op, const constvect inv, vect outv) const override;

  double evalOpByRank(MatrixSparse* Op, int rank) const;
#endif
  
#ifndef SWIG
  
private:
  mutable std::vector<double> _work;
  mutable std::vector<double> _work2;
#endif
};
