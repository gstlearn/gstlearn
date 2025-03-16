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

#include "Estimation/KrigingAlgebra.hpp"
#include "Estimation/KrigOpt.hpp"
#include "Model/ModelGeneric.hpp"
#include "Space/SpaceRN.hpp"
#include "Space/SpacePoint.hpp"
#include "Neigh/ANeigh.hpp"
#include "Matrix/MatrixSquareGeneral.hpp"
#include "Matrix/MatrixSquareSymmetric.hpp"
#include "Matrix/MatrixRectangular.hpp"
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
  KrigingSystemSimpleCase(const KrigingSystemSimpleCase &m) = delete;
  KrigingSystemSimpleCase& operator=(const KrigingSystemSimpleCase &m) = delete;
  virtual ~KrigingSystemSimpleCase();

  int resetData();
  int setKrigOptCalcul(const EKrigOpt& calcul);
  int setKrigOptDataWeights(int iptrWeights, bool flagSet = true);
  int setKrigOptFlagGlobal(bool flag_global);
  int setKrigOptFlagLTerm(bool flag_lterm);

  // The subsequent methods do not require isReady() validation
  int  updKrigOptEstim(int iptrEst, int iptrStd, int iptrVarZ, bool forceNoDual = false);
  int  updKrigOptNeighOnly(int iptrNeigh);
  bool isReady();
  int  estimate(int iech_out);
  void conclusion();

  // Methods used to return algebraic internal information
  int getNDim() const { return (_model != nullptr) ? _model->getNDim() : 0; }
  int getNVar() const { return (_model != nullptr) ? _model->getNVar() : 0; }
  int getNech() const { return (int)_nbgh.size(); }
  int getCovSize() const { return (!_Sigma.empty()) ? _Sigma.getNRows() : 0; }
  int getDriftSize() const { return (!_X.empty()) ? _X.getNCols() : 0; }
  int getNrhs() const { return (!_Sigma0.empty()) ? _Sigma0.getNCols() : 0; }

  VectorInt             getSampleNbgh() const { return _nbgh; }
  VectorVectorDouble    getSampleCoordinates() const;
  VectorDouble          getSampleData() const { return _Z; };
  MatrixSquareSymmetric getLHS() const { return _Sigma; }
  MatrixRectangular     getLHSF() const { return _Sigma0; }
  MatrixRectangular     getRHS() const { return _Sigma0; }
  MatrixRectangular     getRHSF() const { return _X0; }
  MatrixSquareGeneral   getVariance() const { return _Sigma00; }
  MatrixRectangular     getWeights() const;
  MatrixRectangular     getMu() const;
  double                getLTerm() const { return _algebra.getLTerm(); }

private:
  int    _getNVar() const;
  int    _getNbfl() const;
  int    _getNeq() const;
  int    _getNFeq() const;

  void _resetMemoryGeneral();
  bool _isAuthorized() const;

  static void _dumpOptions();
  void _rhsDump();
  void _wgtDump();
  void _estimateCalcul(int status);
  void _simulateCalcul(int status);
  void _neighCalcul(int status, const VectorDouble& tab);
  void _estimateVarZ(int status);
  void _estimateStdv(int status);
  void _estimateEstim(int status);
  void _dumpKrigingResults(int status);
  bool _isCorrect();
  bool _preparNoStat();


  void   _setInternalShortCutVariablesGeneral();
  void   _setInternalShortCutVariablesModel();
  int    _setInternalShortCutVariablesNeigh();

private:
  Db* _dbin;
  Db* _dbout;
  ModelGeneric* _model;
  ANeigh*       _neigh;
  bool          _isReady;

  // Pointers used when plugging KrigingAlgebra (not to be deleted)
  // Note that 'algebra' is mutable not to destroy constness when calling getLambda.
  mutable KrigingAlgebra _algebra;
  mutable KrigOpt        _krigopt;
  VectorVectorInt        _sampleRanks; // Vector of vector of sample indices
  MatrixSquareSymmetric  _Sigma00; // Covariance part for variance
  MatrixSquareSymmetric  _Sigma;   // Covariance part for LHS
  MatrixRectangular      _X;       // Drift part for LHS
  MatrixRectangular      _Sigma0;  // Covariance part for RHS
  MatrixRectangular      _X0;      // Drift par for RHS
  VectorDouble           _Z;       // Vector of Data
  VectorDouble _means;            // Means of the variables (used to center variables)
  VectorDouble _meansTarget;      // Means for target (possible using matLC)

  /// UID for storage
  int  _iptrEst;
  int  _iptrStd;
  int  _iptrVarZ;
  bool _flagEst;
  bool _flagStd;
  bool _flagVarZ;
  bool _flagDataChanged;

  /// Option for Calculation
  EKrigOpt _calcul;

  /// Option for Weights at Data locations
  int  _iptrWeights;
  bool _flagWeights;
  bool _flagSet;


  /// Option for asking for Z * A-1 * Z
  bool _flagLTerm;

  /// Option for Neighboring test
  bool _flagNeighOnly;
  int _iptrNeigh;

  /// Local variables
  int _iechOut;
  int _ndim;
  int _nvar;
  int _nech;
  int _nfeq;
  int _neq;

  /// Working arrays
  mutable VectorInt    _nbgh;
  mutable VectorInt    _dbinUidToBeDeleted;
  mutable VectorInt    _dboutUidToBeDeleted;

  /// Some Space Point allocated once for all
  mutable ASpaceSharedPtr    _space;

  /// Some local flags defined in order to speed up the process
  mutable bool _flagVerr;
  mutable bool _flagNoStat;
};
