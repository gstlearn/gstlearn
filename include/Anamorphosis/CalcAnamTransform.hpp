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
#include "Stats/Selectivity.hpp"

class Db;

class GSTLEARN_EXPORT CalcAnamTransform: public ACalcDbVarCreator
{
public:
  CalcAnamTransform(AAnam* anam);
  CalcAnamTransform(const CalcAnamTransform &r) = delete;
  CalcAnamTransform& operator=(const CalcAnamTransform &r) = delete;
  virtual ~CalcAnamTransform();

  void setAnam(AAnam *anam) { _anam = anam; }
  const AAnam* getAnam() const { return _anam; }
  void setFlagVars(bool flagVars) { _flagVars = flagVars; }
  void setFlagToFactors(bool flagToFactors) { _flagToFactors = flagToFactors; }
  void setFlagFromFactors(bool flagFromFactors) { _flagFromFactors = flagFromFactors; }
  void setFlagZToY(bool flagZToY) { _flagZToY = flagZToY; }
  void setFlagNormalScore(bool flagNormalScore) { _flagNormalScore = flagNormalScore; }
  void setIfacs(const VectorInt &ifacs) { _ifacs = ifacs; }
  void setIptrEst(const VectorInt& iptrEst) { _iptrEst = iptrEst; }
  void setIptrStd(const VectorInt& iptrStd) { _iptrStd = iptrStd; }
  void setSelectivity(Selectivity *selectivity) { _selectivity = selectivity; }
  void setFlagOk(bool flagOk) { _flagOK = flagOk; }
  void setNbsimu(int nbsimu) { _nbsimu = nbsimu; }
  void setProba(double proba) { _proba = proba; }
  void setVerbose(bool verbose) { _verbose = verbose; }
  void setFlagCondExp(bool flagCondExp) { _flagCondExp = flagCondExp; }
  void setFlagUniCond(bool flagUniCond) { _flagUniCond = flagUniCond; }
  void setCvv(double cvv) { _cvv = cvv; }

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
  bool _FactorsToSelectivity();

  int _getNfact() const { return (int) _ifacs.size(); }
  int _getNSel() const { return _selectivity->getVariableNumber(); }

  bool _hasAnam(const EAnam& anamType = EAnam::UNKNOWN) const;
  bool _hasInputVarDefined(int mode = 0) const;
  bool _hasSelectivity() const;
  bool _hasVariableNumber(bool equal1 = false) const;

private:
  int _iattVar;
  int _iattFac;
  int _iattSel;
  bool _flagVars;
  bool _flagToFactors;
  bool _flagFromFactors;
  bool _flagCondExp;
  bool _flagUniCond;
  bool _flagZToY;
  bool _flagNormalScore;
  VectorInt _ifacs;
  VectorInt _iptrEst;
  VectorInt _iptrStd;
  int _nbsimu;
  bool _flagOK;
  bool _verbose;
  double _proba;
  double _cvv;
  AAnam* _anam;
  Selectivity* _selectivity;
};

GSTLEARN_EXPORT int RawToGaussianByLocator(Db* db,
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
GSTLEARN_EXPORT int GaussianToRawByLocator(Db *db,
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

GSTLEARN_EXPORT int FactorToSelectivity(Db *db,
                                        AAnam *anam,
                                        Selectivity *selectivity,
                                        const VectorString &names_est,
                                        const VectorString &names_std,
                                        const NamingConvention &namconv = NamingConvention(
                                            "QT"));
GSTLEARN_EXPORT int ConditionalExpectation(Db *db,
                                           AAnam *anam,
                                           Selectivity *selectivity,
                                           const String &name_est,
                                           const String &name_std,
                                           bool flag_OK = false,
                                           double proba = TEST,
                                           int nbsimu = 0,
                                           bool verbose = false,
                                           const NamingConvention &namconv = NamingConvention(
                                               "CE"));
GSTLEARN_EXPORT int UniformConditioning(Db *db,
                                        AAnam *anam,
                                        Selectivity *selectivity,
                                        const String &names_est,
                                        const String &names_std,
                                        double cvv,
                                        bool verbose = false,
                                        const NamingConvention &namconv = NamingConvention(
                                            "UC"));
