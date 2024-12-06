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

class CovAniso;
class MeshETurbo;

class GSTLEARN_EXPORT ShiftOpStencil: public AShiftOp
{
public:
  ShiftOpStencil(const MeshETurbo* mesh = nullptr,
                 const CovAniso* cova   = nullptr,
                 bool verbose           = true);
  ShiftOpStencil(const ShiftOpStencil& shift);
  ShiftOpStencil& operator=(const ShiftOpStencil& shift);
  virtual ~ShiftOpStencil();
  void normalizeLambdaBySills(const AMesh* mesh) override;

  double getMaxEigenValue() const override;
#ifndef SWIG
  void prodLambda(const constvect x, vect y, const EPowerPT& power) const override;
#endif

#ifndef SWIG
  int _addToDest(const constvect inv, vect outv) const override;
#endif

private:
  void _buildInternal(const MeshETurbo* mesh, const CovAniso* cova, bool verbose);
  void _printStencil() const;

private:
  int _ndim;
  VectorVectorInt _relativeRanks;
  VectorDouble _weights;
  VectorInt _locations;
  VectorBool _onEdge;
};
