/******************************************************************************/
/* COPYRIGHT ARMINES, ALL RIGHTS RESERVED                                     */
/*                                                                            */
/* THE CONTENT OF THIS WORK CONTAINS CONFIDENTIAL AND PROPRIETARY             */
/* INFORMATION OF ARMINES. ANY DUPLICATION, MODIFICATION,                     */
/* DISTRIBUTION, OR DISCLOSURE IN ANY FORM, IN WHOLE, OR IN PART, IS STRICTLY */
/* PROHIBITED WITHOUT THE PRIOR EXPRESS WRITTEN PERMISSION OF ARMINES         */
/*                                                                            */
/* TAG_SOURCE_CG                                                              */
/******************************************************************************/
#pragma once

#include "gstlearn_export.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/IClonable.hpp"
#include "Basic/Vector.hpp"
#include "geoslib_define.h"
#include "csparse_d.h"

#include <functional>

class ALinearOpMulti;
class GSTLEARN_EXPORT APolynomial: public AStringable, public IClonable
{
public:
  APolynomial();
  APolynomial(VectorDouble coeffs);
  APolynomial(const APolynomial& p);
  APolynomial & operator=(const APolynomial& p);
  virtual ~APolynomial();

  void init(VectorDouble coeffs);
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;
  virtual double eval(double x) const = 0;

  virtual void evalOp(cs* /*Op*/,
                      const VectorDouble& /*in*/,
                      VectorDouble& /*out*/) const {};
  virtual void evalOpTraining(cs* /*Op*/,
                      const VectorDouble& /*in*/,
                      VectorVectorDouble& /*out*/,
                      VectorDouble& /*work*/) const {};
  VectorDouble evalOp(cs* /*Op*/, const VectorDouble& /*in*/) const;
  VectorDouble getCoeffs() const { return _coeffs; }
  void setCoeffs(const VectorDouble coeffs){_coeffs = coeffs;}

  int getDegree() const { return static_cast<int>(_coeffs.size());}
  virtual void evalOp(const ALinearOpMulti* Op,const VectorVectorDouble& in, VectorVectorDouble& out) const =0;
  virtual int fit(std::function<double(double)> /*f*/,
                  double /*from*/ = 0.,
                  double /*to*/ = 1.,
                  double /*tol*/ = EPSILON5)
  {
    return 1;
  }

protected:
  VectorDouble _coeffs;
};
