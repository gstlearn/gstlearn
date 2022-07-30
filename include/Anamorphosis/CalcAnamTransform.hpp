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

#include "Calculators/ACalcDbVarCreator.hpp"
#include "Anamorphosis/AAnam.hpp"

class Db;

class GSTLEARN_EXPORT CalcAnamTransform: public ACalcDbVarCreator
{
public:
  CalcAnamTransform(AAnam* anam);
  CalcAnamTransform(const CalcAnamTransform &r) = delete;
  CalcAnamTransform& operator=(const CalcAnamTransform &r) = delete;
  virtual ~CalcAnamTransform();

  void setIatt(int iatt) { _iatt = iatt; }
  void setAnam(AAnam *anam) { _anam = anam; }
  const AAnam* getAnam() const { return _anam; }
  void setFlagNormalScore(bool flagNormalScore) { _flagNormalScore = flagNormalScore; }
  void setFlagZToY(bool flagZToY) { _flagZToY = flagZToY; }
  void setIfacs(const VectorInt &ifacs) { _ifacs = ifacs; }

private:
  virtual bool _check() override;
  virtual bool _preprocess() override;
  virtual bool _run() override;
  virtual bool _postprocess() override;
  virtual void _rollback() override;

  bool _ZToYByHermite();
  bool _YToZByHermite();
  bool _ZToYByNormalScore();
  bool _ZToFactors();

  int _getNfact() const { return (int) _ifacs.size(); }
  bool _toFactors() const { return (! _ifacs.empty()); }

private:
  int _iatt;
  int _number;
  bool _flagZToY;
  bool _flagNormalScore;
  VectorInt _ifacs;
  AAnam* _anam;
};

GSTLEARN_EXPORT int RawToGaussian(Db* db,
                                  AAnam* anam,
                                  const ELoc &locatorType = ELoc::Z,
                                  const NamingConvention &namconv = NamingConvention(
                                      "Y"));
GSTLEARN_EXPORT int RawToGaussian(Db *db,
                                  AAnam *anam,
                                  const String &name,
                                  const NamingConvention &namconv = NamingConvention(
                                      "Y"));
GSTLEARN_EXPORT int NormalScore(Db *db,
                                const NamingConvention &namconv = NamingConvention(
                                    "Gaussian"));
GSTLEARN_EXPORT int GaussianToRaw(Db *db,
                                  AAnam *anam,
                                  const ELoc &locatorType = ELoc::Z,
                                  const NamingConvention &namconv = NamingConvention(
                                      "Z"));

GSTLEARN_EXPORT int GaussianToRaw(Db *db,
                                  AAnam *anam,
                                  const String& name,
                                  const NamingConvention &namconv = NamingConvention(
                                      "Z"));

GSTLEARN_EXPORT int RawToFactor(Db *db,
                                AAnam* anam,
                                const VectorInt &ifacs,
                                const NamingConvention &namconv = NamingConvention(
                                    "Factor"));
GSTLEARN_EXPORT int RawToFactor(Db *db,
                                AAnam* anam,
                                int nfactor,
                                const NamingConvention &namconv = NamingConvention(
                                    "Factor"));

