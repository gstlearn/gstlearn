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

#include "Basic/VectorNumT.hpp"
#include "Matrix/VectorEigen.hpp"
#include "LinearOp/LogStats.hpp"

#define DECLARE_LINEAROP_EIGEN_CG_INTERFACE                                    \
  void evalInverse(const VectorDouble& inv, VectorDouble& outv) const override;\
  void evalInverse(const VectorEigen& inv, VectorEigen& outv) const override;  \
  void evalDirect(const VectorDouble& inv, VectorDouble& outv) const override; \
  void evalDirect(const VectorEigen& inv, VectorEigen& outv) const override;   \
  void setX0(const VectorDouble& x0) override;                                 \
  void mustShowStats(bool status) override;                                    \
  const LogStats& getLogStats() const override;

#define IMPLEMENT_LINEAROP_EIGEN_CG_INTERFACE(TLinOP)                          \
  void TLinOP::evalInverse(const VectorDouble& inv, VectorDouble& outv) const  \
  {                                                                            \
    ALinearOpEigenCG<TLinOP>::evalInverse(inv, outv);                          \
  }                                                                            \
  void TLinOP::evalInverse(const VectorEigen& inv, VectorEigen& outv) const    \
  {                                                                            \
    ALinearOpEigenCG<TLinOP>::evalInverse(inv, outv);                          \
  }                                                                            \
  void TLinOP::evalDirect(const VectorDouble& inv, VectorDouble& outv) const   \
  {                                                                            \
    ALinearOpEigenCG<TLinOP>::evalDirect(inv, outv);                           \
  }                                                                            \
  void TLinOP::evalDirect(const VectorEigen& inv, VectorEigen& outv) const     \
  {                                                                            \
    ALinearOpEigenCG<TLinOP>::evalDirect(inv, outv);                           \
  }                                                                            \
  void TLinOP::setX0(const VectorDouble& x0)                                   \
  {                                                                            \
    ALinearOpEigenCG<TLinOP>::setX0(x0);                                       \
  }                                                                            \
  void TLinOP::mustShowStats(bool status)                                      \
  {                                                                            \
    ALinearOpEigenCG<TLinOP>::mustShowStats(status);                           \
  }                                                                            \
  const LogStats& TLinOP::getLogStats() const                                  \
  {                                                                            \
    return ALinearOpEigenCG<TLinOP>::getLogStats();                            \
  }

class GSTLEARN_EXPORT ILinearOpEigenCG
{
public:
  virtual ~ILinearOpEigenCG() {}
  virtual int getSize() const = 0;
  
  virtual void evalInverse(const VectorDouble& inv, VectorDouble& outv) const = 0;
  virtual void evalInverse(const VectorEigen& inv, VectorEigen& outv) const = 0;

  virtual void evalDirect(const VectorDouble& inv, VectorDouble& outv) const = 0;
  virtual void evalDirect(const VectorEigen& inv, VectorEigen& outv) const = 0;

  virtual void setX0(const VectorDouble& x0) = 0;
  virtual void mustShowStats(bool status) = 0;

  virtual const LogStats& getLogStats() const = 0;
};
