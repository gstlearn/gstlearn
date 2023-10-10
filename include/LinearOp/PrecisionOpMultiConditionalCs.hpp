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

  /// Interface to PrecisionOpMultiConditional
  void makeReady() override;
  int push_back(PrecisionOp* pmatElem, IProjMatrix* projDataElem) override;

  /// Interface to ALinearOp
  void evalInverse(const VectorVectorDouble &vecin,
                   VectorVectorDouble &vecout) const override;

  void mustShowStats(bool status) const { getLogStats().mustShowStats(status); }

private:
  cs* _buildQmult() const;
  ProjMatrix* _buildAmult() const;
  int _buildQpAtA();

private:
  Cholesky _qChol;
};
