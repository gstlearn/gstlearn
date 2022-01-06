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

#include "Anamorphosis/AnamContinuous.hpp"

class GSTLEARN_EXPORT AnamEmpirical: public AnamContinuous
{
private:
  int    _nDisc;
  double _sigma2e;
  VectorDouble _tDisc;

public:
  AnamEmpirical(int ndisc = 100,
                 double sigma2e = TEST);
  AnamEmpirical(const AnamEmpirical &m);
  AnamEmpirical& operator= (const AnamEmpirical &m);
  virtual ~AnamEmpirical();

  void    calculateMeanAndVariance() override;
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  int    getNDisc() const { return _nDisc; }
  double getSigma2e() const { return _sigma2e; }
  const  VectorDouble& getTDisc() const { return _tDisc; }
  void   setSigma2e(double sigma2e) { _sigma2e = sigma2e; }

  void   setNDisc(int ndisc);
  void   setTDisc(const VectorDouble& tdisc);

  double  RawToGaussianValue(double zz) const override;
  double  GaussianToRawValue(double yy) const override;
  int     fit(const VectorDouble& tab);

public:
  bool _isTDiscIndexValid(int i) const;
};
