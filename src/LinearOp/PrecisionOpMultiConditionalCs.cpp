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
#include "LinearOp/PrecisionOpMultiConditionalCs.hpp"
#include "LinearOp/PrecisionOpCs.hpp"

#include "Basic/Law.hpp"
#include "Basic/VectorHelper.hpp"
#include "Matrix/MatrixSquareSymmetric.hpp"
#include "Polynomials/Chebychev.hpp"

#include "csparse_f.h"

#include <functional>

#include <math.h>

PrecisionOpMultiConditionalCs::PrecisionOpMultiConditionalCs()
    : _qChol()
{

}

PrecisionOpMultiConditionalCs::~PrecisionOpMultiConditionalCs()
{
}

int PrecisionOpMultiConditionalCs::push_back(PrecisionOp* pmatElem, IProjMatrix* projDataElem)
{
  PrecisionOpCs* pmatElemCs = dynamic_cast<PrecisionOpCs*>(pmatElem);
  if (pmatElemCs == nullptr)
  {
    messerr("The first argument of 'push_back' should be a pointer to PrecisionOpCs");
    return 1;
  }
  return PrecisionOpMultiConditional::push_back(pmatElem, projDataElem);
}

int PrecisionOpMultiConditionalCs::_buildQmult()
{
  cs *Qmult;
  int number = sizes();


}
