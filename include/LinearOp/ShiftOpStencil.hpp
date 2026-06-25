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
  class Grid;

  /**
   * @brief Helper class for ShiftOpStencil, it stores the weights and the
   * relative shifts needed to apply them.
   */
  class GSTLEARN_EXPORT ShiftStencil
  {
  public:
    ShiftStencil() = default;
    ShiftStencil(const MeshETurbo* mesh, const CovAniso* cova, bool verbose);

    const VectorVectorInt& relativeShifts() const { return _relativeShifts; }

    const VectorInt& absoluteShifts() const { return _absoluteShifts; }

    Id getNWeights() const { return static_cast<Id>(_weights.size()); }

    const VectorDouble& weights() const;

    double getLambda() const { return _lambdaVal; }

    double& getLambda() { return _lambdaVal; }

    void multiplyByValueAndAddDiagonal(double v1 = 1., double v2 = 0.) const;
    void resetModif() const;

    void print() const;

  private:
    VectorVectorInt _relativeShifts;
    VectorInt _absoluteShifts;
    VectorDouble _weights;
    mutable VectorDouble _weightsSimu;
    double _lambdaVal = 0.;
    mutable bool _useModifiedShift = false;
  };

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
   * _weights        Vector of weights (only significant ones are kept)
   * _isInside       Vector telling if each node of the grid is located on its edge
   *                 and should be bypassed for matrix calculations, or not
   */
  class GSTLEARN_EXPORT ShiftOpStencil: public AShiftOp
  {
  public:
    ShiftOpStencil(
      const MeshETurbo* mesh = nullptr,
      const CovAniso* cova = nullptr,
      bool verbose = false);
    ShiftOpStencil(const ShiftOpStencil& shift) = default;
    ShiftOpStencil& operator=(const ShiftOpStencil& shift) = default;
    ~ShiftOpStencil() override = default;
    /// ICloneable interface
    IMPLEMENT_CLONING(ShiftOpStencil)

    void normalizeLambdaBySills(const AMesh* mesh) override;
    void multiplyByValueAndAddDiagonal(double v1 = 1., double v2 = 0.)
      const override;
    void resetModif() const override;

    double getLambda(Id iapex) const override;
    double logDetLambda() const override;

#ifndef SWIG
    Id _addToDest(const constvect inv, vect outv) const override;
#endif

  private:
    double _getMaxEigenValue() const override;
    Id _buildInternal(
      const MeshETurbo* mesh,
      const CovAniso* cova,
      bool verbose);

  private:
    ShiftStencil _stencil;
    VectorBool _isInside;
    bool _useLambdaSingleVal;
    const MeshETurbo* _mesh; // not to be deleted
  };
} // namespace gstlrn
