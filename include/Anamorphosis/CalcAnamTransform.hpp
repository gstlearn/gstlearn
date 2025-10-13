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

#include "geoslib_define.h"
#include "gstlearn_export.hpp"

#include "Anamorphosis/AAnam.hpp"
#include "Anamorphosis/AnamHermite.hpp"
#include "Calculators/ACalcDbVarCreator.hpp"
#include "Stats/Selectivity.hpp"

namespace gstlrn
{
class Db;

class GSTLEARN_EXPORT CalcAnamTransform: public ACalcDbVarCreator
{
public:
  CalcAnamTransform(AAnam* anam);
  CalcAnamTransform(const CalcAnamTransform& r)            = delete;
  CalcAnamTransform& operator=(const CalcAnamTransform& r) = delete;
  virtual ~CalcAnamTransform();

  void setAnam(AAnam* anam) { _anam = anam; }
  const AAnam* getAnam() const { return _anam; }
  void setFlagVars(bool flagVars) { _flagVars = flagVars; }
  void setFlagToFactors(bool flagToFactors) { _flagToFactors = flagToFactors; }
  void setFlagZToY(bool flagZToY) { _flagZToY = flagZToY; }
  void setFlagNormalScore(bool flagNormalScore) { _flagNormalScore = flagNormalScore; }
  void setIfacs(const VectorInt& ifacs) { _ifacs = ifacs; }
  void setIptrEst(const VectorInt& iptrEst) { _iptrEst = iptrEst; }
  void setIptrStd(const VectorInt& iptrStd) { _iptrStd = iptrStd; }
  void setSelectivity(Selectivity* selectivity) { _selectivity = selectivity; }
  void setFlagOk(bool flagOk) { _flagOK = flagOk; }
  void setNbsimu(Id nbsimu) { _nbsimu = nbsimu; }
  void setProba(double proba) { _proba = proba; }
  void setFlagDisjKrig(bool flagDisjKrig) { _flagDisjKrig = flagDisjKrig; }
  void setFlagCondExp(bool flagCondExp) { _flagCondExp = flagCondExp; }
  void setFlagUniCond(bool flagUniCond) { _flagUniCond = flagUniCond; }

private:
  bool _check() override;
  bool _preprocess() override;
  bool _run() override;
  bool _postprocess() override;
  void _rollback() override;

  bool _ZToYByHermite();
  bool _YToZByHermite();
  bool _ZToYByNormalScore();
  bool _ZToFactors();
  bool _FactorsToSelectivity();

  Id _getNfact() const { return static_cast<Id>(_ifacs.size()); }
  Id _getNSel() const { return _selectivity->getNVar(); }

  bool _hasAnam(const EAnam& anamType = EAnam::fromKey("UNKNOWN")) const;
  bool _hasInputVarDefined(Id mode = 0) const;
  bool _hasSelectivity() const;
  bool _hasVariableNumber(bool equal1 = false) const;

  static void _correctForOK(Db* db,
                            Id iech,
                            Id col_est,
                            Id col_std,
                            bool flag_OK,
                            double* krigest,
                            double* krigstd);
  static void _getVectorsForCE(Db* db,
                               Id col_est,
                               Id col_std,
                               bool flag_OK,
                               VectorDouble& krigest,
                               VectorDouble& krigstd);
  static Id _conditionalExpectation(Db* db,
                                     AAnam* anam,
                                     const Selectivity* selectivity,
                                     Id iptr0,
                                     Id col_est,
                                     Id col_std,
                                     bool flag_OK,
                                     double proba,
                                     Id nbsimu);
  static Id _uniformConditioning(Db* db,
                                  AnamHermite* anam,
                                  Selectivity* selectivity,
                                  Id iptr0,
                                  Id col_est,
                                  Id col_var);
  static Id _ceZ(Db* db,
                  const AnamHermite* anam,
                  const Selectivity* selectivity,
                  Id iptr0,
                  Id col_est,
                  Id col_std,
                  Id nbsimu,
                  bool flag_OK);
  static Id _ceT(Id mode,
                  Db* db,
                  const Selectivity* selectivity,
                  Id iptr0,
                  Id col_est,
                  Id col_std,
                  const VectorDouble& ycuts,
                  Id nbsimu,
                  bool flag_OK);
  static Id _ceQ(Db* db,
                  const AnamHermite* anam,
                  const Selectivity* selectivity,
                  Id iptr0,
                  Id col_est,
                  Id col_std,
                  const VectorDouble& ycuts,
                  Id nbsimu,
                  bool flag_OK);
  static Id _ceB(Db* db,
                  const Selectivity* selectivity,
                  Id iptr0,
                  const VectorDouble& ycuts);
  static Id _ceM(Db* db, const Selectivity* selectivity, Id iptr0);
  static Id _ceQuant(Db* db,
                      const AnamHermite* anam,
                      const Selectivity* selectivity,
                      Id iptr0,
                      Id col_est,
                      Id col_std,
                      double proba,
                      bool flag_OK);

private:
  Id _iattVar;
  Id _iattFac;
  Id _iattSel;
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
  Id _nbsimu;
  bool _flagOK;
  double _proba;
  AAnam* _anam;
  Selectivity* _selectivity;
};

// TODO : rename functions with a lower case at the beginning
GSTLEARN_EXPORT Id DisjunctiveKriging(Db* db,
                                       AAnam* anam,
                                       Selectivity* selectivity,
                                       const VectorString& name_est,
                                       const VectorString& name_std,
                                       const NamingConvention& namconv = NamingConvention(
                                         "DK"));
GSTLEARN_EXPORT Id ConditionalExpectation(Db* db,
                                           AAnam* anam,
                                           Selectivity* selectivity        = nullptr,
                                           const String& name_est          = "",
                                           const String& name_std          = "",
                                           bool flag_OK                    = false,
                                           double proba                    = TEST,
                                           Id nbsimu                      = 0,
                                           const NamingConvention& namconv = NamingConvention(
                                             "CE"));
GSTLEARN_EXPORT Id UniformConditioning(Db* db,
                                        AAnam* anam,
                                        Selectivity* selectivity,
                                        const String& name_est,
                                        const String& name_varz,
                                        const NamingConvention& namconv = NamingConvention(
                                          "UC"));

GSTLEARN_EXPORT Id anamPointToBlock(AAnam* anam,
                                     Id verbose,
                                     double cvv,
                                     double coeff,
                                     double mu);
} // namespace gstlrn