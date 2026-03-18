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
#include "LinearOp/ASimulableMatrix.hpp"
#include "geoslib_define.h"
#include "Matrix/MatrixSparse.hpp"

namespace gstlrn
{

ASimulableMatrix::ASimulableMatrix()
:
 _chol(nullptr)
{
}

ASimulableMatrix::~ASimulableMatrix()
{
  delete _chol;
}

Id ASimulableMatrix::_addToDest(const constvect whitenoise, vect outv) const
{
  return getQMat().addToDest(whitenoise, outv);
}

const CholeskySparse& ASimulableMatrix::getChol() const
{
    if (_chol == nullptr)
    {
      _chol = new CholeskySparse(getQMat());
    }
    return *_chol;
}

Id ASimulableMatrix::_addSimulateToDest(const constvect whitenoise, vect outv) const
{
  return getChol().addSimulateToDest(whitenoise, outv);
}

double ASimulableMatrix::computeLogDet(Id nMC) const
{
  DECLARE_UNUSED(nMC)
  const auto& chol = getChol();
  return chol.computeLogDeterminant();
}
} // namespace gstlrn