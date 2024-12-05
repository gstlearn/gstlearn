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
class AMesh;


class GSTLEARN_EXPORT ShiftOpStencil: public AShiftOp
{
  public:
    ShiftOpStencil(const AMesh* amesh = nullptr, const CovAniso* cova = nullptr, bool verbose = false);
    ShiftOpStencil(const ShiftOpStencil& shift);
    ShiftOpStencil& operator=(const ShiftOpStencil& shift);
    virtual ~ShiftOpStencil();
    void normalizeLambdaBySills(const AMesh* mesh) override;
    void multiplyByValueAndAddDiagonal(double v1 = 1.,double v2 = 0.) override;

    double getMaxEigenValue() const override;
#ifndef SWIG
    void addProdLambda(const constvect x, vect y, const EPowerPT& power) const override;
#endif 

#ifndef SWIG
    int _addToDest(const constvect inv, vect outv) const override;
#endif

};

