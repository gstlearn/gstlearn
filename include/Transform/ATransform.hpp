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

#include "geoslib_define.h"

#include "Basic/AStringable.hpp"
#include "Basic/ICloneable.hpp"
#include "Basic/ListParams.hpp"
#include "Basic/VectorNumT.hpp"
#include "Covariances/ACov.hpp"

namespace gstlrn
{

class ListParams;

class GSTLEARN_EXPORT ATransform: public ICloneable, public AStringable
{
public:
  ATransform() {};
  ATransform(const ATransform& r)            = default;
  ATransform& operator=(const ATransform& r) = default;
  virtual ~ATransform()                      = default;

  /**
   * @brief Gaussian to Raw data transformation
   * @param h
   * @return double
   */
  virtual double transform(double h) const
  {
    DECLARE_UNUSED(h);
    return TEST;
  }
  virtual bool hasParameters() const = 0;
  virtual void _printParams(std::stringstream& sstr, const AStringFormat* strfmt) const
  {
    DECLARE_UNUSED(sstr, strfmt);
  }
  /**
   * @brief Raw to Gaussian data transformation
   * @param y
   * @return double
   */
  virtual double inverseTransform(double y) const;
  String toString(const AStringFormat* strfmt = nullptr) const override;

  virtual String getName() const = 0;
#ifndef SWIG
  void transformVec(constvect in, vect out) const;
  void condExpVec(constvect mu, constvect sigma, vect out, Id power = 1) const;
  void inverseTransformVec(constvect in, vect out) const;
#endif
  VectorDouble transformVec(const VectorDouble& in) const;
  VectorDouble inverseTransformVec(const VectorDouble& in) const;
  VectorDouble condExpVec(const VectorDouble& mu, const VectorDouble& sigma, Id power = 1) const;
  virtual double evalJacobian(double x) const;
  virtual double condExp(double mu, double sigma, Id power = 1) const;

  virtual void updateTransform() {}
  virtual void initParams(double min = 0., double max = INF) { DECLARE_UNUSED(min, max) };
  virtual void appendParams(ListParams& listParams)
  {
    DECLARE_UNUSED(listParams);
  }
  virtual VectorDouble getParams() const { return VectorDouble(); }
  double evalLogJacobianVec(constvect in) const;

  void setNMonteCarlo(Id n) { _nMonteCarlo = n; }
  Id getNMonteCarlo() const { return _nMonteCarlo; }

private:
  Id _nMonteCarlo = 100;
};
} // namespace gstlrn