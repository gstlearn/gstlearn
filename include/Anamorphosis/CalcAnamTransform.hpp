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
  void setFlagZToY(bool flagZToY) { _flagZToY = flagZToY; }
  void setFlagNormalScore(bool flagNormalScore) { _flagNormalScore = flagNormalScore; }
  void setIfacs(const VectorInt &ifacs) { _ifacs = ifacs; }
  void setIptrEst(const VectorInt& iptrEst) { _iptrEst = iptrEst; }
  void setIptrStd(const VectorInt& iptrStd) { _iptrStd = iptrStd; }
  void setSelectivity(Selectivity *selectivity) { _selectivity = selectivity; }
  void setFlagOk(bool flagOk) { _flagOK = flagOk; }
  void setNbsimu(int nbsimu) { _nbsimu = nbsimu; }
  void setProba(double proba) { _proba = proba; }
  void setFlagDisjKrig(bool flagDisjKrig) { _flagDisjKrig = flagDisjKrig; }
  void setFlagCondExp(bool flagCondExp) { _flagCondExp = flagCondExp; }
  void setFlagUniCond(bool flagUniCond) { _flagUniCond = flagUniCond; }

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

  bool _hasAnam(const EAnam& anamType = EAnam::fromKey("UNKNOWN")) const;
  bool _hasInputVarDefined(int mode = 0) const;
  bool _hasSelectivity() const;
  bool _hasVariableNumber(bool equal1 = false) const;

private:
  int _iattVar;
  int _iattFac;
  int _iattSel;
  bool _flagVars;
  bool _flagToFactors;
  bool _flagDisjKrig;
  bool _flagCondExp;
  bool _flagUniCond;
  bool _flagZToY;
  bool _flagNormalScore;
  VectorInt _ifacs;
  VectorInt _iptrEst;
  VectorInt _iptrStd;
  int _nbsimu;
  bool _flagOK;
  double _proba;
  AAnam* _anam;
  Selectivity* _selectivity;
};

// TODO : rename functions with a lower case at the beginning
GSTLEARN_EXPORT int DisjunctiveKriging(Db *db,
                                       AAnam *anam,
                                       Selectivity *selectivity,
                                       const VectorString &name_est,
                                       const VectorString &name_std,
                                       const NamingConvention &namconv = NamingConvention(
                                           "DK"));
GSTLEARN_EXPORT int ConditionalExpectation(Db *db,
                                           AAnam *anam,
                                           Selectivity *selectivity = nullptr,
                                           const String &name_est = "",
                                           const String &name_std = "",
                                           bool flag_OK = false,
                                           double proba = TEST,
                                           int nbsimu = 0,
                                           const NamingConvention &namconv = NamingConvention(
                                               "CE"));
GSTLEARN_EXPORT int UniformConditioning(Db *db,
                                        AAnam *anam,
                                        Selectivity *selectivity,
                                        const String &name_est,
                                        const String &name_varz,
                                        const NamingConvention &namconv = NamingConvention(
                                            "UC"));
