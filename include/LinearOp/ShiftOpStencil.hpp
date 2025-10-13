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
#include "LinearOp/AShiftOp.hpp"
#include "geoslib_define.h"

namespace gstlrn
{
class CovAniso;
class MeshETurbo;
class AMesh;

/**
 * @brief This is an implementation of ShiftOp dedicated to case where:
 * - the target is a regular grid
 * - the meshing is elaborated as a TurboMeshing
 * - the covariance is stationary
 *
 * The different members are:
 * _relativeShifts For each vector, gives the vector of shifts, with respect
 *                 to the target node (in relative indices)
 * _absoluteShifts Vector of shifts to calculate where the weights should apply
 *                 calculated on the global target grid.
 *                 This can only be used if the grid has no selection
 * _weights        Vector of weights (only significative ones are kept)
 * _isInside       Vector telling if each node of the grid is located on its edge
 *                 and should be bypassed for matrix calculations, or not
 */
class GSTLEARN_EXPORT ShiftOpStencil: public AShiftOp
{
public:
  ShiftOpStencil(const MeshETurbo* mesh = nullptr,
                 const CovAniso* cova   = nullptr,
                 bool verbose           = false);
  ShiftOpStencil(const ShiftOpStencil& shift);
  ShiftOpStencil& operator=(const ShiftOpStencil& shift);
  virtual ~ShiftOpStencil();
  /// ICloneable interface
  IMPLEMENT_CLONING(ShiftOpStencil)

  void normalizeLambdaBySills(const AMesh* mesh) override;
  void multiplyByValueAndAddDiagonal(double v1 = 1., double v2 = 0.) const override;
  void resetModif() const override;

  double getLambda(Id iapex) const override;
  double logDetLambda() const override;

#ifndef SWIG
  Id _addToDest(const constvect inv, vect outv) const override;
#endif

private:
  double _getMaxEigenValue() const override;
  Id _buildInternal(const MeshETurbo* mesh, const CovAniso* cova, bool verbose);
  void _printStencil() const;
  Id _getNWeights() const { return static_cast<Id>(_weights.size()); }

private:
  VectorVectorInt _relativeShifts;
  VectorInt _absoluteShifts;
  VectorDouble _weights;
  mutable VectorDouble _weightsSimu;
  VectorBool _isInside;
  double _lambdaVal;
  bool _useLambdaSingleVal;
  mutable bool _useModifiedShift;
  const MeshETurbo* _mesh; // not to be deleted
};
} // namespace gstlrn