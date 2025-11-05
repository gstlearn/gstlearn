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

#include "Basic/ListParams.hpp"
#include "Basic/VectorNumT.hpp"
#include "Covariances/ACov.hpp"
#include "geoslib_define.h"
#include "gstlearn_export.hpp"
namespace gstlrn
{

class ListParams;

class GSTLEARN_EXPORT ATransform
{
public:
  ATransform() {};
  ATransform(const ATransform& r)            = delete;
  ATransform& operator=(const ATransform& r) = delete;

  virtual double transform(double h) const
  {
    DECLARE_UNUSED(h);
    return TEST;
  }
  virtual String getName() const = 0;
  virtual double inverseTransform(double y) const;
#ifndef SWIG
  void transformVec(constvect in, vect out) const;
  void condExpVec(constvect mu, constvect sigma,  vect out, Id power = 1) const;
  void inverseTransformVec(constvect in, vect out) const;
#endif
  VectorDouble transformVec(const VectorDouble& in) const;
  VectorDouble inverseTransformVec(const VectorDouble& in) const;
  virtual double evalJacobian(double x) const;
  virtual double condExp(double mu, double sigma, Id power = 1) const;
  VectorDouble condExpVec(const VectorDouble& mu, const VectorDouble& sigma, Id power = 1) const;
  double evalLogJacobianVec(constvect in) const;
  virtual void updateTransform() {}
  virtual void initParams(){};
  virtual void appendParams(ListParams& listParams)
  {
    DECLARE_UNUSED(listParams);
  }
  void setNMonteCarlo(Id n) { _nMonteCarlo = n; }
  Id getNMonteCarlo() const { return _nMonteCarlo; }

  virtual ~ATransform() = default;
  private:
  Id _nMonteCarlo = 100;
};
} // namespace gstlrn