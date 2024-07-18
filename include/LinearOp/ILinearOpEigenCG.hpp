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

class GSTLEARN_EXPORT ILinearOpEigenCG
{
public:
  virtual int getSize() const = 0;
  
  virtual void evalInverse(const VectorDouble& inv, VectorDouble& outv) const = 0;
  virtual void evalInverse(const VectorEigen& inv, VectorEigen& outv) const =  0;

  virtual void evalDirect(const VectorDouble& inv, VectorDouble& outv) const = 0;
  virtual void evalDirect(const VectorEigen& inv, VectorEigen& outv) const = 0;

  virtual void setX0(const VectorDouble& x0) = 0;
  virtual void mustShowStats(bool status) = 0;

  virtual const LogStats& getLogStats() const = 0;
};
