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
#include "geoslib_define.h"

#include "Anamorphosis/AnamDiscrete.hpp"

class GSTLEARN_EXPORT AnamDiscreteIR: public AnamDiscrete
{
public:
  AnamDiscreteIR();
  AnamDiscreteIR(const AnamDiscreteIR &m);
  AnamDiscreteIR& operator= (const AnamDiscreteIR &m);
  virtual ~AnamDiscreteIR();

  int fit(const VectorDouble& tab, int verbose=0);
  void calculateMeanAndVariance() override;
  VectorDouble z2f(int nfact, const VectorInt& ifacs, double z) const override;
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  double getRCoef() const { return _rCoef; }
  void   setRCoef(double rcoef) { _rCoef = rcoef; }

private:
  int _stats_residuals(int verbose,
                       int nech,
                       const VectorDouble& tab,
                       int *nsorted,
                       double *mean,
                       double *residuals,
                       double *T,
                       double *Q);
  double _getResidual(int iclass, double z) const;

private:
  double _rCoef;
};
