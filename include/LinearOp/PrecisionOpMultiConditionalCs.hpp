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

#include "gstlearn_export.hpp"

#include "LinearOp/Cholesky.hpp"
#include "LinearOp/PrecisionOpMultiConditional.hpp"

#include <vector>

class PrecisionOp;
class IProjMatrix;

/**
 * Class to store objects for SPDE
 */
class GSTLEARN_EXPORT PrecisionOpMultiConditionalCs : public PrecisionOpMultiConditional {

public:
  PrecisionOpMultiConditionalCs();
  PrecisionOpMultiConditionalCs(const PrecisionOpMultiConditionalCs &m)= delete;
  PrecisionOpMultiConditionalCs& operator= (const PrecisionOpMultiConditionalCs &m)= delete;
  virtual ~PrecisionOpMultiConditionalCs();

  void push_back(PrecisionOp* pmatElem, IProjMatrix* projDataElem) override;

private:
  cs* _Q;
  Cholesky _qChol;
};
