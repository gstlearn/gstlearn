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

#include "Anamorphosis/AAnam.hpp"
#include "Calculators/ACalcInterpolator.hpp"
#include "Matrix/MatrixDense.hpp"
#include "Matrix/MatrixSymmetric.hpp"

namespace gstlrn
{
class KrigingSystem;
class Db;
class DbGrid;

class GSTLEARN_EXPORT Krigtest_Res
{
public:
  Id ndim;                // Space dimension
  Id nvar;                // Number of variables
  Id nech;                // Number of Neighboring samples
  Id CSize;               // Number of drift equations in the Drift part
  Id DSize;               // Number of Equations of the Covariance part
  Id nrhs;                // Number of R.H.S. vectors
  VectorInt nbgh;         // Ranks of the neighboring samples
  VectorVectorDouble xyz; // Coordinates of the neighboring samples [ndim][nech]
  VectorDouble data;      // Usable values at neighboring samples [neq]
  MatrixSymmetric lhs;    // L.H.S. Covariance part (neq * neq)
  MatrixDense lhsF;       // L.H.S. Drift part
  MatrixDense rhs;        // R.H.S. Covariance part (neq * nrhs)
  MatrixDense rhsF;       // R.H.S. Drift part  (nbfl * nrhs)
  MatrixDense wgt;        // Vector of weights (neq * nrhs)
  MatrixDense mu;         // Vector of Lagrange parameters (nbfl * nrhs)
  MatrixSquare var;       // Matrix of Target-Target Variance (nvar * nvar)

  /// Has a specific implementation in the Target language
  DECLARE_TOTL;
};

// TODO : Create KrigingParam which inherits from InterpolatorParam
class GSTLEARN_EXPORT CalcKriging: public ACalcInterpolator
{
public:
  CalcKriging(bool flag_est = true, bool flag_std = true, bool flag_varZ = false);
  CalcKriging(const CalcKriging& r)            = delete;
  CalcKriging& operator=(const CalcKriging& r) = delete;
  virtual ~CalcKriging();

  void setPriorCov(const MatrixSymmetric& priorCov) { _priorCov = priorCov; }
  void setPriorMean(const VectorDouble& priorMean) { _priorMean = priorMean; }
  void setFlagBayes(bool flagBayes) { _flagBayes = flagBayes; }
  void setIechSingleTarget(Id iechSingleTarget) { _iechSingleTarget = iechSingleTarget; }
  void setVerboseSingleTarget(bool verbose) { _verboseSingleTarget = verbose; }
  void setAnam(AAnam* anam) { _anam = anam; }
  void setFlagGam(bool flagGam) { _flagGam = flagGam; }
  void setFlagXvalidEst(Id flagXvalidEst) { _flagXvalidEst = flagXvalidEst; }
  void setFlagXvalidStd(Id flagXvalidStd) { _flagXvalidStd = flagXvalidStd; }
  void setFlagXvalidVarZ(Id flagXvalidVarZ) { _flagXvalidVarZ = flagXvalidVarZ; }
  void setFlagXvalid(bool flagXvalid) { _flagXvalid = flagXvalid; }
  void setFlagKfold(bool flag_kfold) { _flagKfold = flag_kfold; }
  void setFlagNeighOnly(bool flagNeighOnly) { _flagNeighOnly = flagNeighOnly; }

  Krigtest_Res getKtest() const { return _ktest; }

private:
  bool _check() override;
  bool _preprocess() override;
  bool _run() override;
  bool _postprocess() override;
  void _rollback() override;

  void _storeResultsForExport(const KrigingSystem& ksys);

private:
  bool _flagEst;
  bool _flagStd;
  bool _flagVarZ;

  VectorString _nameCoord;

  bool _flagBayes;
  VectorDouble _priorMean;
  MatrixSymmetric _priorCov;

  Id _iechSingleTarget;
  bool _verboseSingleTarget;

  bool _flagGam;
  AAnam* _anam;

  bool _flagXvalid;
  bool _flagKfold;
  Id _flagXvalidEst;
  Id _flagXvalidStd;
  Id _flagXvalidVarZ;

  bool _flagNeighOnly;
  Id _nbNeigh;

  Id _iptrEst;
  Id _iptrStd;
  Id _iptrVarZ;
  Id _iptrNeigh;

  Krigtest_Res _ktest;
};

GSTLEARN_EXPORT Id kriging(Db* dbin,
                           Db* dbout,
                           ModelGeneric* model,
                           ANeigh* neigh,
                           bool flag_est                   = true,
                           bool flag_std                   = true,
                           bool flag_varz                  = false,
                           const KrigOpt& krigopt          = KrigOpt(),
                           const NamingConvention& namconv = NamingConvention("Kriging"));
GSTLEARN_EXPORT Id krigcell(Db* dbin,
                            Db* dbout,
                            ModelGeneric* model,
                            ANeigh* neigh,
                            bool flag_est                   = true,
                            bool flag_std                   = true,
                            const KrigOpt& krigopt          = KrigOpt(),
                            const NamingConvention& namconv = NamingConvention("KrigCell"));
GSTLEARN_EXPORT Id kribayes(Db* dbin,
                            Db* dbout,
                            ModelGeneric* model,
                            ANeigh* neigh,
                            const VectorDouble& prior_mean   = VectorDouble(),
                            const MatrixSymmetric& prior_cov = MatrixSymmetric(),
                            bool flag_est                    = true,
                            bool flag_std                    = true,
                            const NamingConvention& namconv  = NamingConvention("Bayes"));
GSTLEARN_EXPORT Id kriggam(Db* dbin,
                           Db* dbout,
                           ModelGeneric* model,
                           ANeigh* neigh,
                           AAnam* anam,
                           const NamingConvention& namconv = NamingConvention("KrigGam"));
GSTLEARN_EXPORT Krigtest_Res krigtest(Db* dbin,
                                      Db* dbout,
                                      ModelGeneric* model,
                                      ANeigh* neigh,
                                      Id iech0               = 0,
                                      const KrigOpt& krigopt = KrigOpt(),
                                      bool verbose           = true);
GSTLEARN_EXPORT Id xvalid(Db* db,
                          ModelGeneric* model,
                          ANeigh* neigh,
                          bool flag_kfold                 = false,
                          Id flag_xvalid_est              = 1,
                          Id flag_xvalid_std              = 1,
                          Id flag_xvalid_varz             = 0,
                          const KrigOpt& krigopt          = KrigOpt(),
                          const NamingConvention& namconv = NamingConvention("Xvalid"));
GSTLEARN_EXPORT Id test_neigh(Db* dbin,
                              Db* dbout,
                              ModelGeneric* model,
                              ANeigh* neigh,
                              const NamingConvention& namconv = NamingConvention("Neigh"));
} // namespace gstlrn