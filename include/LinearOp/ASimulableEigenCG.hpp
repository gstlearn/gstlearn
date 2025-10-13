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

#include "LinearOp/ASimulable.hpp"

namespace gstlrn
{
/**
 * @brief This class extends ASimulable to make it working with
 *        Eigen conjugate gradient algorithm
 *
 *        This class is similar to ALinearOpEigenCG.
 * 
 * @tparam TLinOP Concrete class that inherits from ASimulableEigenCG
 */
template<typename TLinOP>
class ASimulableEigenCG : public Eigen::EigenBase<TLinOP>, // No Export because it's a template
                          public ASimulable
{
public:
  virtual ~ASimulableEigenCG() {};

#ifndef SWIG
public:
  // Required typedefs, constants, and method:
  typedef double Scalar;
  typedef double RealScalar;
  typedef Id StorageIndex;
  enum
  {
    ColsAtCompileTime    = Eigen::Dynamic,
    MaxColsAtCompileTime = Eigen::Dynamic,
    IsRowMajor           = false
  };
  
  Eigen::Index rows() const { return getSize(); }
  Eigen::Index cols() const { return getSize(); }

  template<typename Rhs>
  Eigen::Product<TLinOP,Rhs,Eigen::AliasFreeProduct> operator*(const Eigen::MatrixBase<Rhs>& x) const {
    return Eigen::Product<TLinOP,Rhs,Eigen::AliasFreeProduct>(*(dynamic_cast<const TLinOP*>(this)), x.derived());
  }
#endif
};

}