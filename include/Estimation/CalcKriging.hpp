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

#include "Calculators/ACalcInterpolator.hpp"

class Db;
class DbGrid;
class KrigingSystem;

class GSTLEARN_EXPORT CalcKriging: public ACalcInterpolator
{
public:
  CalcKriging(bool flag_est = true, bool flag_std = true, bool flag_varZ = false);
  CalcKriging(const CalcKriging &r) = delete;
  CalcKriging& operator=(const CalcKriging &r) = delete;
  virtual ~CalcKriging();

  void setCalcul(const EKrigOpt &calcul) { _calcul = calcul; }
  void setMatCl(const VectorVectorDouble &matCl) { _matCL = matCl; }
  void setNdisc(const VectorInt &ndisc) { _ndisc = ndisc; }
  void setRankColCok(const VectorInt &rankColCok) { _rankColCok = rankColCok; }
  void setFlagDgm(bool flagDgm) { _flagDGM = flagDgm; }
  void setRCoeff(double rCoeff) { _rCoeff = rCoeff; }
  void setPriorCov(const VectorDouble &priorCov) { _priorCov = priorCov; }
  void setPriorMean(const VectorDouble &priorMean) { _priorMean = priorMean; }
  void setFlagBayes(bool flagBayes) { _flagBayes = flagBayes; }
  void setFlagProf(bool flagProf) { _flagProf = flagProf; }
  void setIechSingleTarget(int iechSingleTarget) { _iechSingleTarget = iechSingleTarget; }
  void setFlagPerCell(bool flagPerCell) { _flagPerCell = flagPerCell; }
  void setAnam(AAnam *anam) { _anam = anam; }
  void setFlagGam(bool flagGam) { _flagGam = flagGam; }
  void setFlagXvalidEst(int flagXvalidEst) { _flagXvalidEst = flagXvalidEst; }
  void setFlagXvalidStd(int flagXvalidStd) { _flagXvalidStd = flagXvalidStd; }
  void setFlagXvalidVarZ(int flagXvalidVarZ) { _flagXvalidVarZ = flagXvalidVarZ; }
  void setFlagXvalid(bool flagXvalid) { _flagXvalid = flagXvalid; }
  void setFlagKfold(bool flag_kfold) { _flagKfold = flag_kfold; }
  void setFlagNeighOnly(bool flagNeighOnly) { _flagNeighOnly = flagNeighOnly; }

  Krigtest_Res getKtest() const { return _ktest; }

private:
  virtual bool _check() override;
  virtual bool _preprocess() override;
  virtual bool _run() override;
  virtual bool _postprocess() override;
  virtual void _rollback() override;
  int _getNVar() const override;

  void _storeResultsForExport(const KrigingSystem& ksys);

private:
  bool _flagEst;
  bool _flagStd;
  bool _flagVarZ;

  EKrigOpt  _calcul;
  VectorInt _ndisc;
  VectorInt _rankColCok;
  VectorVectorDouble _matCL;

  bool _flagDGM;
  double _rCoeff;

  bool _flagBayes;
  VectorDouble _priorMean;
  VectorDouble _priorCov;

  bool _flagProf;

  int _iechSingleTarget;

  bool _flagPerCell;

  bool _flagGam;
  AAnam* _anam;

  bool _flagXvalid;
  bool _flagKfold;
  int  _flagXvalidEst;
  int  _flagXvalidStd;
  int  _flagXvalidVarZ;

  bool _flagNeighOnly;
  int  _nbNeigh;

  int _iptrEst;
  int _iptrStd;
  int _iptrVarZ;
  int _iptrNeigh;

  Krigtest_Res _ktest;
};

GSTLEARN_EXPORT int kriging(Db *dbin,
                            Db *dbout,
                            Model *model,
                            ANeighParam *neighparam,
                            const EKrigOpt &calcul = EKrigOpt::PONCTUAL,
                            bool flag_est = true,
                            bool flag_std = true,
                            bool flag_varz = false,
                            VectorInt ndisc = VectorInt(),
                            VectorInt rank_colcok = VectorInt(),
                            VectorVectorDouble matCL = VectorVectorDouble(),
                            const NamingConvention& namconv = NamingConvention("Kriging"));
GSTLEARN_EXPORT int krigcell(Db *dbin,
                             Db *dbout,
                             Model *model,
                             ANeighParam *neighparam,
                             bool flag_est,
                             bool flag_std,
                             VectorInt ndisc,
                             VectorInt rank_colcok = VectorInt(),
                             const NamingConvention& namconv = NamingConvention("KrigCell"));
GSTLEARN_EXPORT int krigdgm(Db *dbin,
                            DbGrid *dbout,
                            Model *model,
                            ANeighParam *neighparam,
                            bool flag_est = true,
                            bool flag_std = true,
                            bool flag_varz = false,
                            double rval = 1.,
                            const NamingConvention& namconv = NamingConvention("DGM"));
GSTLEARN_EXPORT int kribayes(Db *dbin,
                             Db *dbout,
                             Model *model,
                             ANeighParam *neighparam,
                             const VectorDouble& dmean = VectorDouble(),
                             const VectorDouble& dcov = VectorDouble(),
                             bool flag_est = true,
                             bool flag_std = true,
                             const NamingConvention& namconv = NamingConvention("Bayes"));
GSTLEARN_EXPORT int krigprof(Db *dbin,
                             Db *dbout,
                             Model *model,
                             ANeighParam *neighparam,
                             bool flag_est = true,
                             bool flag_std = true,
                             const NamingConvention& namconv = NamingConvention("KrigProf"));
GSTLEARN_EXPORT int kriggam(Db *dbin,
                            Db *dbout,
                            Model *model,
                            ANeighParam *neighparam,
                            AAnam *anam,
                            const NamingConvention& namconv = NamingConvention("KrigGam"));
GSTLEARN_EXPORT Krigtest_Res krigtest(Db *dbin,
                                      Db *dbout,
                                      Model *model,
                                      ANeighParam *neighparam,
                                      int iech0,
                                      const EKrigOpt &calcul = EKrigOpt::PONCTUAL,
                                      VectorInt ndisc = VectorInt(),
                                      bool forceDebug = true);
GSTLEARN_EXPORT int xvalid(Db *db,
                           Model *model,
                           ANeighParam *neighparam,
                           bool flag_kfold = false,
                           int flag_xvalid_est = 1,
                           int flag_xvalid_std = 1,
                           int flag_xvalid_varz = 0,
                           VectorInt rank_colcok = VectorInt(),
                           const NamingConvention& namconv = NamingConvention("Xvalid"));
GSTLEARN_EXPORT int test_neigh(Db *dbin,
                               Db *dbout,
                               Model *model,
                               ANeighParam *neighparam,
                               const NamingConvention& namconv = NamingConvention("Neigh"));
