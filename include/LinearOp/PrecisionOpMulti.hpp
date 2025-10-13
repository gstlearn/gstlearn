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

#include "Basic/AStringFormat.hpp"
#include "LinearOp/ASimulable.hpp"
#include "Matrix/MatrixDense.hpp"
#include "Matrix/MatrixSymmetric.hpp"
#include "gstlearn_export.hpp"

#include "Basic/AStringable.hpp"
#include "Basic/VectorNumT.hpp"
#include "Basic/VectorT.hpp"
#include "LinearOp/CholeskyDense.hpp"
#include "LinearOp/PrecisionOp.hpp"
#include "Model/Model.hpp"
#include <vector>

#define IND(i, j, nvar) j* nvar + i - (j * (j + 1)) / 2

namespace gstlrn
{
class CholeskyDense;
class Model;
class AStringable;
class AStringFormat;
class ASimulable;

// This class is dedicated to the multivariate Model.
// It creates a vector of precision operators (matrix-free).
class GSTLEARN_EXPORT PrecisionOpMulti: public AStringable, public ASimulable
{
public:
  PrecisionOpMulti(Model* model               = nullptr,
                   const VectorMeshes& meshes = VectorMeshes(),
                   bool stencil               = false,
                   bool buildOp               = true);
  PrecisionOpMulti(const PrecisionOpMulti& m)            = delete;
  PrecisionOpMulti& operator=(const PrecisionOpMulti& m) = delete;
  virtual ~PrecisionOpMulti();

  /// AStringable Interface
  String toString(const AStringFormat* strfmt = nullptr) const override;

  Id getSize() const override;

  double computeLogDet(Id nMC = 1) const override;
  std::pair<double, double> rangeEigenValQ() const;

protected:
#ifndef SWIG
  Id _addToDest(const constvect vecin, vect vecout) const override;
  Id _addSimulateToDest(const constvect vecin, vect vecout) const override;
#endif

  void buildQop(bool stencil = false);
  Id size(Id imesh) const;
  Id _getNCov() const;
  Id _getCovInd(Id i) const { return _covList[i]; }
  Id _getNVar() const;
  Id _getNMesh() const;

private:
  bool _checkReady() const;
  virtual void _buildQop(bool stencil = false);
  bool _isValidModel(Model* model);
  bool _isValidMeshes(const std::vector<const AMesh*>& meshes);
  bool _isNoStat(Id istruct) const { return _isNoStatForVariance[istruct]; }
  bool _matchModelAndMeshes() const;

  Id _buildGlobalMatricesStationary(Id icov);
  Id _buildLocalMatricesNoStat(Id icov);
  Id _buildMatrices();
  void _popsClear();
  void _computeSize();

protected:
  std::vector<PrecisionOp*> _pops;
  VectorBool _isNoStatForVariance;
  std::vector<MatrixSymmetric> _sills;
  std::vector<std::vector<MatrixSymmetric>> _localSills; // Local Sills for non-stationary covariances
  std::vector<VectorVectorDouble> _invCholSillsNoStat;
  std::vector<VectorVectorDouble> _cholSillsNoStat;
  std::vector<CholeskyDense> _invCholSillsStat; // Stationary Sills
  std::vector<CholeskyDense> _cholSillsStat;    // Cholesky of the Sills

  Model* _model;                     // Not to be deleted. TODO : make it const
  std::vector<const AMesh*> _meshes; // Not to be deleted
  Id _size;

private:
  bool _isValid;
  VectorInt _covList;
  VectorInt _nmeshList;
  bool _allStat;
  bool _ready;

  mutable VectorVectorDouble _works;
  mutable VectorDouble _workTot;
};
} // namespace gstlrn
