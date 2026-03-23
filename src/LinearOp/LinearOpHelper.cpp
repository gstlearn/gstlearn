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

#include "LinearOp/LinearOpHelper.hpp"
#include "Basic/VectorHelper.hpp"

namespace gstlrn
{

  double LinearOpHelper::powerIteration(ALinearOp* op, Id niter)
  {

    Id size = op->getSize();
    auto b = VH::simulateGaussian(size);
    VectorDouble bnext(size);
    double norm2 = VH::innerProductCV(b, b);

    for (Id i = 0; i < niter; ++i)
    {
      op->evalDirect(b, bnext); // b_next = A * b
      norm2 = VH::innerProductCV(bnext, bnext);
      double norm = std::sqrt(norm2);
      if (i == niter - 1)
      {
        break;
      }
      for (Id j = 0; j < size; ++j)
      {
        b[j] = bnext[j] / norm; // Normalize
      }
    }

    return VH::innerProductCV(b, bnext);
  }

} // namespace gstlrn
