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

#include "Stats/PCA.hpp"

class GSTLEARN_EXPORT AnamDiscreteDD: public AnamDiscrete
{
public:
  AnamDiscreteDD(double mu = 1., double scoef = 0.);
  AnamDiscreteDD(const AnamDiscreteDD &m);
  AnamDiscreteDD& operator= (const AnamDiscreteDD &m);
  virtual ~AnamDiscreteDD();

  /// ASerializable Interface
  int dumpToNF(const String& neutralFilename, bool verbose = false) const;
  static AnamDiscreteDD* createFromNF(const String& neutralFilename, bool verbose = false);
  static AnamDiscreteDD* createFromNF2(const String& neutralFilename, bool verbose = false);

  /// AAnam Interface
  const EAnam&  getType() const override { return EAnam:: DISCRETE_DD; }

  /// AnamDiscrete Interface
  void calculateMeanAndVariance() override;
  VectorDouble z2f(int nfact, const VectorInt& ifacs, double z) const override;

  /// AStringable Interface
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  VectorDouble factors_exp(int verbose);
  VectorDouble factors_maf(int verbose);
  VectorDouble factors_mod();
  VectorDouble chi2I(const VectorDouble& chi, int mode);

  AnamDiscreteDD* create(double mu = 1., double scoef = 0.);
  int  fit(const VectorDouble& tab, int verbose=0);

  PCA& getMAF() { return _maf; }
  double getMu() const { return _mu; }
  double getSCoef() const { return _sCoef; }
  const VectorDouble& getI2Chi() const { return _i2Chi; }
  VectorDouble getPcaZ2F() const { return _maf.getZ2F(); }
  VectorDouble getPcaF2Z() const { return _maf.getF2Z(); }

  void setMu(double mu) { _mu = mu; }
  void setSCoef(double scoef) { _sCoef = scoef; }
  void setPcaZ2F(VectorDouble pcaz2f) { _maf.setPcaZ2F(pcaz2f); }
  void setPcaF2Z(VectorDouble pcaf2z) { _maf.setPcaF2Z(pcaf2z); }
  void setI2Chi(const VectorDouble& i2Chi) { _i2Chi = i2Chi; }

protected:
  virtual int _deserialize(FILE* file, bool verbose = false);
  virtual int _serialize(FILE* file, bool verbose = false) const override;

  virtual int _deserialize2(std::istream& is, bool verbose) override;

private:
  int _stats(int nech, const VectorDouble& tab);
  VectorDouble _generator(const VectorDouble& vecc,
                          const VectorDouble& veca,
                          const VectorDouble& vecb,
                          VectorDouble& eigvec,
                          VectorDouble& eigval);
  void _lambda_to_mul();

private:
  double _mu;
  double _sCoef;
  PCA    _maf;
  VectorDouble _i2Chi;
};
