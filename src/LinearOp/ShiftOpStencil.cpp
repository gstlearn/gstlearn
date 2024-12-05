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
#include "LinearOp/ShiftOpStencil.hpp"

#include "LinearOp/AShiftOp.hpp"
#include "Mesh/AMesh.hpp"
#include "Covariances/CovAniso.hpp"
#include "geoslib_define.h"



ShiftOpStencil::ShiftOpStencil(const AMesh* amesh, const CovAniso* cova, bool verbose)
{
    DECLARE_UNUSED(amesh,cova,verbose)
}


ShiftOpStencil::ShiftOpStencil(const ShiftOpStencil& shift)
  : AShiftOp(shift)
{
 DECLARE_UNUSED(shift)
}

ShiftOpStencil& ShiftOpStencil::operator=(const ShiftOpStencil &shift)
{
  if (this != &shift)
  {
  
  }
  return *this;
}

ShiftOpStencil::~ShiftOpStencil()
{
 
}

int ShiftOpStencil::_addToDest(const constvect inv, vect outv) const
{
  DECLARE_UNUSED(inv,outv)
  return 0;
}

void ShiftOpStencil::normalizeLambdaBySills(const AMesh* mesh) 
{
    DECLARE_UNUSED(mesh)
}

void ShiftOpStencil::addProdLambda(const constvect x, vect y, const EPowerPT& power) const
{
    DECLARE_UNUSED(x,y,power)
    messerr("ShiftOpStencil::addProdLambda not implemented");
}

double ShiftOpStencil::getMaxEigenValue() const
{
    return TEST;
}
void ShiftOpStencil::multiplyByValueAndAddDiagonal(double v1, double v2)
{
    DECLARE_UNUSED(v1,v2)
} 