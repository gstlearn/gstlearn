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

#include "Basic/AStringable.hpp"
#include "Basic/ICloneable.hpp"
#include "Basic/VectorNumT.hpp"
#include "geoslib_define.h"
#include "Matrix/csparse_d.h"

#include <functional>

class ALinearOpMulti;
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

  virtual void evalOp(cs* Op,
                      const VectorDouble& inv,
                      VectorDouble& outv) const { DECLARE_UNUSED(Op,inv,outv); };
  virtual void evalOpTraining(cs* Op,
                      const VectorDouble& inv,
                      VectorVectorDouble& outv,
                      VectorDouble& work) const { DECLARE_UNUSED(Op,inv,outv,work); };
  VectorDouble evalOp(cs* Op, const VectorDouble& inv) const;
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
