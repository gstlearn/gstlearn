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

#include "LinearOp/AShiftOp.hpp"
#include "geoslib_define.h"

class CovAniso;
class MeshETurbo;

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
    void normalizeLambdaBySills(const AMesh* mesh) override;
    void multiplyByValueAndAddDiagonal(double v1 = 1., double v2 = 0.) override
    {
      DECLARE_UNUSED(v1);
      DECLARE_UNUSED(v2);
    };

    double getMaxEigenValue() const override;
#ifndef SWIG
    void addProdLambda(const constvect x, vect y, const EPowerPT& power) const override;
#endif 

#ifndef SWIG
  int _addToDest(const constvect inv, vect outv) const override;
#endif

private:
  int _buildInternal(const MeshETurbo* mesh, const CovAniso* cova, bool verbose);
  void _printStencil() const;
  int _getNWeights() const { return (int) _weights.size(); }

private:
  VectorVectorInt _relativeShifts;
  VectorInt _absoluteShifts;
  VectorDouble _weights;
  VectorBool _isInside; 

  const MeshETurbo* _mesh; // not to be deleted
};
