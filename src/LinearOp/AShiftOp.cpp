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
#include "LinearOp/AShiftOp.hpp"

#include <math.h>

AShiftOp::AShiftOp()
: _Lambda()
, _napices(0)
{
}


AShiftOp::AShiftOp(const AShiftOp& shift)
: _Lambda()
, _napices(0)
{
    DECLARE_UNUSED(shift)
}


AShiftOp& AShiftOp::operator=(const AShiftOp &shift)
{
  if (this != &shift)
  {
    _napices = shift._napices;
    _Lambda = shift._Lambda;
  }
  return *this;
}

AShiftOp::~AShiftOp()
{
}


void AShiftOp::prodLambda(const VectorDouble& x,
                          vect y,
                           const EPowerPT& power) const
{
  constvect xv(x.data(),x.size());
  prodLambda(xv,y,power);
  
}

void AShiftOp::prodLambda(const constvect x,
                           VectorDouble& y,
                           const EPowerPT& power) const
{
    vect yv(y.data(),y.size());
    prodLambda(x,yv,power);
}

void AShiftOp::prodLambda(const VectorDouble& x,
                           VectorDouble& y,
                           const EPowerPT& power) const
{
  constvect xv(x.data(),x.size());
  vect yv(y.data(),y.size());
  prodLambda(xv,yv,power);
}