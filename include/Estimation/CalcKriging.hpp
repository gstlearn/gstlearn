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

#include "Enum/EKrigOpt.hpp"

#include "Calculators/ACalcInterpolator.hpp"
#include "Matrix/MatrixRectangular.hpp"
#include "Matrix/MatrixSquareSymmetric.hpp"
#include "Anamorphosis/AAnam.hpp"

class Db;
class DbGrid;
class KrigingSystem;

class GSTLEARN_EXPORT Krigtest_Res
{
public:
  int ndim; // Space dimension
  int nvar; // Number of variables
  int nech; // Number of Neighboring samples
  int neq;  // Number of Equations in the Kriging/CoKriging system
  int nrhs; // Number of R.H.S. vectors (= nvar)
  VectorInt nbgh;    // Ranks of the neighboring samples
  VectorVectorDouble xyz;    // Coordinates of the neighboring samples [ndim][nech]
  VectorDouble data;         // Usable values at neighboring samples [neq]
  MatrixSquareSymmetric lhs; // L.H.S. of the Kriging system (neq * neq)
  MatrixRectangular rhs;     // R.H.S. of the Kriging system (neq * nvar)
  MatrixRectangular wgt;     // Vector of weights [nvar][nech]
  MatrixSquareGeneral var;   // Matrix of Target-Target Variance (nvar * nvar)
  MatrixRectangular zam;     // Vector of pre-calculations

  /// Has a specific implementation in the Target language
  DECLARE_TOTL;
};

// TODO : Create KrigingParam which inherits from InterpolatorParam
class GSTLEARN_EXPORT CalcKriging: public ACalcInterpolator
{
public:
  CalcKriging(bool flag_est = true, bool flag_std = true, bool flag_varZ = false);
  CalcKriging(const CalcKriging &r) = delete;
  CalcKriging& operator=(const CalcKriging &r) = delete;
  virtual ~CalcKriging();

  void setCalcul(const EKrigOpt &calcul);
  void setMatLC(const MatrixRectangular* matLC) { _matLC = matLC; }
  void setNdisc(const VectorInt& ndiscs) { _ndiscs = ndiscs; }
  void setRankColCok(const VectorInt& rankColCok) { _rankColCok = rankColCok; }
  void setFlagDgm(bool flagDgm) { _flagDGM = flagDgm; }
  void setPriorCov(const MatrixSquareSymmetric &priorCov) { _priorCov = priorCov; }
  void setPriorMean(const VectorDouble &priorMean) { _priorMean = priorMean; }
  void setFlagBayes(bool flagBayes) { _flagBayes = flagBayes; }
  void setFlagProf(bool flagProf) { _flagProf = flagProf; }
  void setIechSingleTarget(int iechSingleTarget) { _iechSingleTarget = iechSingleTarget; }
  void setVerboseSingleTarget(bool verbose) { _verboseSingleTarget = verbose; }
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
  VectorInt _ndiscs;
  VectorInt _rankColCok;
  const MatrixRectangular* _matLC;

  bool _flagDGM;
  VectorString _nameCoord;

  bool _flagBayes;
  VectorDouble _priorMean;
  MatrixSquareSymmetric _priorCov;

  bool _flagProf;

  int  _iechSingleTarget;
  bool _verboseSingleTarget;

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
                            std::shared_ptr<ANeigh>& neigh,
                            const EKrigOpt &calcul = EKrigOpt::fromKey("POINT"),
                            bool flag_est = true,
                            bool flag_std = true,
                            bool flag_varz = false,
                            const VectorInt& ndiscs = VectorInt(),
                            const VectorInt& rank_colcok = VectorInt(),
                            const MatrixRectangular* matLC = nullptr,
                            const NamingConvention& namconv = NamingConvention("Kriging"));
GSTLEARN_EXPORT int krigcell(Db *dbin,
                             Db *dbout,
                             Model *model,
                             std::shared_ptr<ANeigh>& neigh,
                             bool flag_est = true,
                             bool flag_std = true,
                             const VectorInt& ndiscs = VectorInt(),
                             const VectorInt& rank_colcok = VectorInt(),
                             const NamingConvention& namconv = NamingConvention("KrigCell"));
GSTLEARN_EXPORT int kribayes(Db *dbin,
                             Db *dbout,
                             Model *model,
                             std::shared_ptr<ANeigh>& neigh,
                             const VectorDouble& prior_mean = VectorDouble(),
                             const MatrixSquareSymmetric& prior_cov = MatrixSquareSymmetric(),
                             bool flag_est = true,
                             bool flag_std = true,
                             const NamingConvention& namconv = NamingConvention("Bayes"));
GSTLEARN_EXPORT int krigprof(Db *dbin,
                             Db *dbout,
                             Model *model,
                             std::shared_ptr<ANeigh>& neigh,
                             bool flag_est = true,
                             bool flag_std = true,
                             const NamingConvention& namconv = NamingConvention("KrigProf"));
GSTLEARN_EXPORT int kriggam(Db *dbin,
                            Db *dbout,
                            Model *model,
                            std::shared_ptr<ANeigh>& neigh,
                            AAnam *anam,
                            const NamingConvention& namconv = NamingConvention("KrigGam"));
GSTLEARN_EXPORT Krigtest_Res krigtest(Db *dbin,
                                      Db *dbout,
                                      Model *model,
                                      std::shared_ptr<ANeigh>& neigh,
                                      int iech0 = 0,
                                      const EKrigOpt &calcul = EKrigOpt::fromKey("POINT"),
                                      const VectorInt& ndiscs = VectorInt(),
                                      bool flagPerCell = false,
                                      bool verbose = true);
GSTLEARN_EXPORT int xvalid(Db *db,
                           Model *model,
                           std::shared_ptr<ANeigh>& neigh,
                           bool flag_kfold = false,
                           int flag_xvalid_est = 1,
                           int flag_xvalid_std = 1,
                           int flag_xvalid_varz = 0,
                           const VectorInt& rank_colcok = VectorInt(),
                           const NamingConvention& namconv = NamingConvention("Xvalid"));
GSTLEARN_EXPORT int test_neigh(Db *dbin,
                               Db *dbout,
                               Model *model,
                               std::shared_ptr<ANeigh>& neigh,
                               const NamingConvention& namconv = NamingConvention("Neigh"));
