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

#include "LinearOp/PrecisionOpMultiConditional.hpp"

#ifndef SWIG
  #include <Eigen/src/Core/Matrix.h>
#endif

class PrecisionOp;
class IProjMatrix;

/**
 * Class to store objects for SPDE
 */
class GSTLEARN_EXPORT PrecisionOpMultiConditionalCs : public PrecisionOpMultiConditional 
{
public:
  PrecisionOpMultiConditionalCs();
  PrecisionOpMultiConditionalCs(const PrecisionOpMultiConditionalCs &m)= delete;
  PrecisionOpMultiConditionalCs& operator= (const PrecisionOpMultiConditionalCs &m)= delete;
  virtual ~PrecisionOpMultiConditionalCs();

  /// Interface to PrecisionOpMultiConditional
  void makeReady() override;
  int push_back(PrecisionOp* pmatElem, IProjMatrix* projDataElem) override;
  double computeLogDetOp(int nbsimu = 1, int seed = 123) const override;

  /// Interface to ALinearOp
  void evalInverse(const std::vector<Eigen::VectorXd> &vecin,
                   std::vector<Eigen::VectorXd> &vecout) const override;

  void mustShowStats(bool status) const { getLogStats().mustShowStats(status); }

private:
  MatrixSparse* _buildQmult() const;
  ProjMatrix*   _buildAmult() const;
  int _buildQpAtA();

  MatrixSparse* _Q;
};
