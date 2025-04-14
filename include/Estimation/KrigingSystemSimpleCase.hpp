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

#include "Basic/VectorNumT.hpp"
#include "gstlearn_export.hpp"

#include "Estimation/KrigingAlgebraSimpleCase.hpp"
#include "Estimation/KrigOpt.hpp"
#include "Model/ModelGeneric.hpp"
#include "Space/SpaceRN.hpp"
#include "Space/SpacePoint.hpp"
#include "Neigh/ANeigh.hpp"
#include "Matrix/MatrixSquareGeneral.hpp"
#include "Matrix/MatrixSymmetric.hpp"
#include "Matrix/MatrixDense.hpp"
#include "LinearOp/CholeskyDense.hpp"
#include "Enum/EKrigOpt.hpp"

class Db;
class DbGrid;
class Model;
class ModelGeneric;
class ANeigh;
class CovCalcMode;
class ECalcMember;
class NeighImage;
class AAnam;
class ACov;
class KrigingAlgebra;
class KrigOpt;

class GSTLEARN_EXPORT KrigingSystemSimpleCase
{
public:
  KrigingSystemSimpleCase(Db* dbin,
                          Db* dbout,
                          const ModelGeneric* model,
                          ANeigh* neigh);
  KrigingSystemSimpleCase(const KrigingSystemSimpleCase& m)            = delete;
  KrigingSystemSimpleCase& operator=(const KrigingSystemSimpleCase& m) = delete;
  virtual ~KrigingSystemSimpleCase();

  int setKrigOptCalcul(const EKrigOpt& calcul);
  int setKrigOptDataWeights(int iptrWeights, bool flagSet = true);
  int setKrigOptFlagGlobal(bool flag_global);
  int setKrigOptFlagLTerm(bool flag_lterm);

  // The subsequent methods do not require isReady() validation
  int updKrigOptEstim(int iptrEst, int iptrStd, int iptrVarZ, bool forceNoDual = false);
  bool isReady();
  void updateLHS(KrigingAlgebraSimpleCase& algebra, ModelGeneric& model);
  int estimate(int iechout,
               SpacePoint& pin,
               SpacePoint& pout,
               VectorDouble& tabwork,
               KrigingAlgebraSimpleCase& algebra,
               ModelGeneric& model,
               ANeigh* neigh = nullptr);
  void conclusion();

  KrigingAlgebraSimpleCase& getAlgebra() { return _algebra; }

  // Methods used to return algebraic internal information
  int getNDim() const { return (_model != nullptr) ? _model->getNDim() : 0; }
  int getNVar() const { return (_model != nullptr) ? _model->getNVar() : 0; }

  VectorVectorDouble getSampleCoordinates(KrigingAlgebraSimpleCase& algebra, int iechout) const;
  static MatrixDense getWeights(KrigingAlgebraSimpleCase& algebra);
  static MatrixDense getMu(KrigingAlgebraSimpleCase& algebra);
  double getLTerm() const { return _algebra.getLTerm(); }
  ModelGeneric* getModel() const { return _model; }

private:
  int _getNVar() const;
  int _getNbfl() const;
  int _getNeq(int nech) const;
  int _getNFeq() const;

  void _resetMemoryGeneral();

  static void _dumpOptions();
  void _rhsDump(KrigingAlgebraSimpleCase& algebra) const;
  static void _wgtDump(KrigingAlgebraSimpleCase& algebra);
  void _estimateCalcul(int status, int iechout, KrigingAlgebraSimpleCase& algebra) const;
  void _simulateCalcul(int status);
  void _neighCalcul(int status, const VectorDouble& tab, int iechout);
  void _estimateVarZ(int status, int iechout, KrigingAlgebraSimpleCase& algebra) const;
  void _estimateStdv(int status, int iechout, KrigingAlgebraSimpleCase& algebra) const;
  void _estimateEstim(int status, KrigingAlgebraSimpleCase& algebra, int iechout) const;
  void _dumpKrigingResults(int status, int iechout, KrigingAlgebraSimpleCase* algebra) const;
  bool _isCorrect();
  bool _preparNoStat();

  void _setInternalShortCutVariablesGeneral();
  void _setInternalShortCutVariablesModel();

private:
  Db* _dbin;
  Db* _dbout;
  ModelGeneric* _model;
  ANeigh* _neigh;
  bool _isReady;

  // Pointers used when plugging KrigingAlgebra (not to be deleted)
  // Note that 'algebra' is mutable not to destroy constness when calling getLambda.
  mutable KrigingAlgebraSimpleCase _algebra;
  mutable KrigOpt _krigopt;
  VectorDouble _means;            // Means of the variables (used to center variables)
  VectorDouble _meansTarget;      // Means for target (possible using matLC)

  /// UID for storage
  int _iptrEst;
  int _iptrStd;
  int _iptrVarZ;
  bool _flagEst;
  bool _flagStd;
  bool _flagVarZ;
  bool _flagDataChanged;

  /// Option for Weights at Data locations
  int _iptrWeights;
  bool _flagWeights;
  bool _flagSet;

  /// Option for asking for Z * A-1 * Z
  bool _flagLTerm;

  /// Option for Neighboring test
  bool _neighUnique;
  int _iptrNeigh;

  /// Local variables
  int _ndim;
  int _nvar;
  int _nfeq;

  /// Working arrays
  VectorInt _dbinUidToBeDeleted;
  VectorInt _dboutUidToBeDeleted;

  /// Some local flags defined in order to speed up the process
  bool _flagVerr;
  bool _flagNoStat;
};
