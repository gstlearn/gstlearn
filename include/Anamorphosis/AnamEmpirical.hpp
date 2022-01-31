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

#include "Basic/ASerializable.hpp"
#include "Anamorphosis/AnamContinuous.hpp"
#include "Anamorphosis/EAnam.hpp"

class GSTLEARN_EXPORT AnamEmpirical: public AnamContinuous
{
public:
  AnamEmpirical(int ndisc = 100, double sigma2e = TEST);
  AnamEmpirical(const AnamEmpirical &m);
  AnamEmpirical& operator= (const AnamEmpirical &m);
  virtual ~AnamEmpirical();

  /// ASerializable Interface
  int dumpToNF(const String& neutralFilename, bool verbose = false) const;
  static AnamEmpirical* createFromNF(const String& neutralFilename, bool verbose = false);

  /// AAnam Interface
  const EAnam&  getType() const override { return EAnam:: EMPIRICAL; }

  /// AnamContinuous Interface
  void    calculateMeanAndVariance() override;
  double  RawToGaussianValue(double zz) const override;
  double  GaussianToRawValue(double yy) const override;

  /// AStringable Interface
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  AnamEmpirical* create(int ndisc = 100, double sigma2e = TEST);
  int    getNDisc() const { return _nDisc; }
  double getSigma2e() const { return _sigma2e; }
  const  VectorDouble& getTDisc() const { return _tDisc; }
  void   setSigma2e(double sigma2e) { _sigma2e = sigma2e; }

  void   setNDisc(int ndisc);
  void   setTDisc(const VectorDouble& tdisc);
  int    fit(const VectorDouble& tab);
  bool   isTDiscIndexValid(int i) const;

protected:
  /// ASerializable Interface
  virtual int _deserialize(FILE* file, bool verbose = false) override;
  virtual int _serialize(FILE* file, bool verbose = false) const override;

private:
  int    _nDisc;
  double _sigma2e;
  VectorDouble _tDisc;
};
