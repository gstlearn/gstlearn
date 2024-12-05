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

#include "LinearOp/PrecisionOpMultiConditional.hpp"
#include <vector>

class PrecisionOp;
class IProj;
class CholeskySparse;

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
  int push_back(PrecisionOp* pmatElem, IProj* projDataElem) override;
  double computeLogDetOp(int nbsimu = 1) const override;

  /// Interface to ALinearOp
  void evalInverse(const std::vector<std::vector<double>> &vecin,
                   std::vector<std::vector<double>> &vecout) const override;

  void mustShowStats(bool status) const { getLogStats().mustShowStats(status); }
private :
  void _clear();
private:
  MatrixSparse* _buildQmult() const;
  ProjMatrix*   _buildAmult() const;
  int _buildQpAtA();

  MatrixSparse* _Q;
  mutable CholeskySparse* _chol;
};
