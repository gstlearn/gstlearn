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
#include "Anamorphosis/EAnam.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/ASerializable.hpp"

class GSTLEARN_EXPORT AnamDiscreteIR: public AnamDiscrete
{
public:
  AnamDiscreteIR(double rcoef = 0.);
  AnamDiscreteIR(const AnamDiscreteIR &m);
  AnamDiscreteIR& operator= (const AnamDiscreteIR &m);
  virtual ~AnamDiscreteIR();

  /// Interface AStringable
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  /// ASerializable Interface
  int dumpToNF(const String& neutralFilename, bool verbose = false) const;
  static AnamDiscreteIR* createFromNF(const String& neutralFilename, bool verbose = false);

  /// AAnam Interface
  const EAnam&  getType() const override { return EAnam:: DISCRETE_IR; }

  /// AnamDiscrete Interface
  void calculateMeanAndVariance() override;
  VectorDouble z2f(int nfact, const VectorInt& ifacs, double z) const override;

  AnamDiscreteIR* create(double rcoef = 0.);
  int fit(const VectorDouble& tab, int verbose=0);
  double getRCoef() const { return _rCoef; }
  void   setRCoef(double rcoef) { _rCoef = rcoef; }

protected:
  virtual int _deserialize(FILE* file, bool verbose = false) override;
  virtual int _serialize(FILE* file, bool verbose = false) const override;

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
