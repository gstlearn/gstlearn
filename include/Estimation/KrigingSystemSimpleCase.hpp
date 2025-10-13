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

#include "Enum/EKrigOpt.hpp"
#include "Estimation/KrigOpt.hpp"
#include "Estimation/KrigingAlgebraSimpleCase.hpp"
#include "LinearOp/CholeskyDense.hpp"
#include "Matrix/MatrixDense.hpp"
#include "Model/ModelGeneric.hpp"
#include "Neigh/ANeigh.hpp"
#include "Space/SpacePoint.hpp"
#include "Space/SpaceRN.hpp"

namespace gstlrn
{

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

  Id setKrigOptCalcul(const EKrigOpt& calcul);
  Id setKrigOptDataWeights(Id iptrWeights, bool flagSet = true);
  Id setKrigOptFlagGlobal(bool flag_global);
  Id setKrigOptFlagLTerm(bool flag_lterm);

  // The subsequent methods do not require isReady() validation
  Id updKrigOptEstim(Id iptrEst, Id iptrStd, Id iptrVarZ, bool forceNoDual = false);
  bool isReady(VectorDouble& tabwork);
  void updateLHS(KrigingAlgebraSimpleCase& algebra,
                 ModelGeneric& model,
                 VectorDouble& tabwork);
  Id estimate(Id iechout,
              SpacePoint& pin,
              SpacePoint& pout,
              VectorDouble& tabwork,
              KrigingAlgebraSimpleCase& algebra,
              ModelGeneric& model,
              ANeigh* neigh = nullptr);
  void conclusion();

  KrigingAlgebraSimpleCase& getAlgebra() { return _algebra; }

  // Methods used to return algebraic internal information
  Id getNDim() const { return (_model != nullptr) ? static_cast<Id>(_model->getNDim()) : 0; }
  Id getNVar() const { return (_model != nullptr) ? _model->getNVar() : 0; }

  VectorVectorDouble getSampleCoordinates(KrigingAlgebraSimpleCase& algebra, Id iechout) const;
  static MatrixDense getWeights(KrigingAlgebraSimpleCase& algebra);
  static MatrixDense getMu(KrigingAlgebraSimpleCase& algebra);
  double getLTerm() const { return _algebra.getLTerm(); }
  ModelGeneric* getModel() const { return _model; }

private:
  Id _getNVar() const;
  Id _getNbfl() const;
  Id _getNeq(Id nech) const;
  Id _getNFeq() const;

  void _resetMemoryGeneral();

  static void _dumpOptions();
  void _rhsDump(KrigingAlgebraSimpleCase& algebra) const;
  static void _wgtDump(KrigingAlgebraSimpleCase& algebra);
  void _estimateCalcul(Id status, Id iechout, KrigingAlgebraSimpleCase& algebra) const;
  void _simulateCalcul(Id status);
  void _neighCalcul(Id status, const VectorDouble& tab, Id iechout);
  void _estimateVarZ(Id status, Id iechout, KrigingAlgebraSimpleCase& algebra) const;
  void _estimateStdv(Id status, Id iechout, KrigingAlgebraSimpleCase& algebra) const;
  void _estimateEstim(Id status, KrigingAlgebraSimpleCase& algebra, Id iechout) const;
  void _dumpKrigingResults(Id status, Id iechout, KrigingAlgebraSimpleCase* algebra) const;
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
  VectorDouble _means;       // Means of the variables (used to center variables)
  VectorDouble _meansTarget; // Means for target (possible using matLC)

  /// UID for storage
  Id _iptrEst;
  Id _iptrStd;
  Id _iptrVarZ;
  bool _flagEst;
  bool _flagStd;
  bool _flagVarZ;
  bool _flagDataChanged;

  /// Option for Weights at Data locations
  Id _iptrWeights;
  bool _flagWeights;
  bool _flagSet;

  /// Option for asking for Z * A-1 * Z
  bool _flagLTerm;

  /// Option for Neighboring test
  bool _neighUnique;
  Id _iptrNeigh;

  /// Local variables
  Id _ndim;
  Id _nvar;
  Id _nfeq;

  /// Working arrays
  VectorInt _dbinUidToBeDeleted;
  VectorInt _dboutUidToBeDeleted;

  /// Some local flags defined in order to speed up the process
  bool _flagVerr;
  bool _flagNoStat;
};
} // namespace gstlrn