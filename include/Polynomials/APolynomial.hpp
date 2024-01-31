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

#include "Basic/AStringable.hpp"
#include "Basic/ICloneable.hpp"
#include "Basic/VectorNumT.hpp"
#include "geoslib_define.h"

#include <functional>

class MatrixSparse;
class ALinearOpMulti;
class cs;

class GSTLEARN_EXPORT APolynomial: public AStringable, public ICloneable
{
public:
  APolynomial();
  APolynomial(VectorDouble coeffs);
  APolynomial(const APolynomial& p);
  APolynomial & operator=(const APolynomial& p);
  virtual ~APolynomial();

  /// Interface for AStringable
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  void init(VectorDouble coeffs);
  virtual double eval(double x) const = 0;
#ifndef SWIG
  virtual void evalOp(MatrixSparse* Op,
                      const VectorDouble& inv,
                      VectorDouble& outv) const { DECLARE_UNUSED(Op,inv,outv); };
  virtual void evalOpTraining(MatrixSparse* Op,
                      const VectorDouble& inv,
                      VectorVectorDouble& outv,
                      VectorDouble& work) const { DECLARE_UNUSED(Op,inv,outv,work); };
  VectorDouble evalOp(MatrixSparse* Op, const VectorDouble& inv) const;
#endif
  VectorDouble getCoeffs() const { return _coeffs; }
  void setCoeffs(const VectorDouble coeffs) {_coeffs = coeffs;}

  int getDegree() const { return static_cast<int>(_coeffs.size());}
  virtual void evalOp(const ALinearOpMulti* Op,const VectorVectorDouble& inv, VectorVectorDouble& outv) const = 0;
  virtual int fit(std::function<double(double)> f,
                  double from = 0.,
                  double to = 1.,
                  double tol = EPSILON5)
  {
    DECLARE_UNUSED(f,from,to,tol);
    return 1;
  }

protected:
  VectorDouble _coeffs;
};
