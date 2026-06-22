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

#include "Basic/ICloneable.hpp"
#include "Basic/VectorNumT.hpp"
#include "Polynomials/APolynomial.hpp"

class AShiftOp;

namespace gstlrn
{
  class ALinearOp;

  class GSTLEARN_EXPORT ClassicalPolynomial: public APolynomial
  {
  public:
    ClassicalPolynomial();
    ClassicalPolynomial(const VectorDouble&);
    virtual ~ClassicalPolynomial();

    /// ICloneable interface
    IMPLEMENT_CLONING(ClassicalPolynomial)

    double eval(double x) const override;
    // void evalDerivOp(ShiftOpMatrix* shiftOp,
    //                         const constvect& inv,
    //                         vect& outv,
    //                         Id iapex,
    //                         Id igparam);
    // static void evalDerivOpOptim(ShiftOpMatrix* shiftOp,
    //                              vect& temp1,
    //                              vect& temp2,
    //                              vect& outv,
    //                              const VectorVectorDouble& workpoly,
    //                              Id iapex,
    //                              Id igparam);
#ifndef SWIG
    // void evalDerivOp(ShiftOpMatrix* shiftOp,const VectorDouble& inv,
    //                  VectorDouble& outv,Id iapex,Id igparam)const;

    // void evalDerivOpOptim(ShiftOpMatrix* shiftOp,Eigen::VectorXd& temp1,Eigen::VectorXd& temp2,
    //                      Eigen::VectorXd& outv,const std::vector<Eigen::VectorXd>& workpoly,Id iapex,Id igparam)const;
    // void evalOp(const ALinearOpMulti* /*Op*/,
    //            const std::vector<Eigen::VectorXd>& /*inv*/,
    //           std::vector<Eigen::VectorXd>& /*outv*/) const override { }

    void evalOpTraining(
      MatrixSparse* Op,
      const constvect inv,
      VectorVectorDouble& store,
      VectorDouble& work) const override;
    void evalOpCumul(MatrixSparse* Op, const constvect inv, vect outv) const;
    void
      evalOp(MatrixSparse* Op, const constvect inv, vect outv) const override;
    double evalOpByRank(MatrixSparse* S, Id rank) const override;
#endif

#ifndef SWIG

    void _addEvalOp(const ALinearOp* Op, const constvect inv, vect outv)
      const override;
    Id setWorkArrays(vect work1, vect work2);

  protected:
    Id initWorkArrays(vect& work1, vect& work2, size_t size) const;

  private:
    mutable VectorDouble _w1, _w2; // local work arrays
    mutable vect _work1, _work2; // shared work arrays
#endif
  };
} // namespace gstlrn
