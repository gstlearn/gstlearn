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

#include "Model/AModelOptim.hpp"

class Model;
class Db;

/**
 * \brief
 * Class which, starting from an experimental variogram, enables fitting the
 * various parameters of a Covariance part of a Model
 */
class GSTLEARN_EXPORT ModelOptimLikelihood: public AModelOptim
{
public:
  ModelOptimLikelihood(Model* model);
  ModelOptimLikelihood(const ModelOptimLikelihood& m);
  ModelOptimLikelihood& operator=(const ModelOptimLikelihood& m);
  virtual ~ModelOptimLikelihood();

  int loadEnvironment(Db* db, bool flagSPDE = false, bool verbose = false);

private:
  bool _checkConsistency();

private:
  // Use SPDE approach (TRUE); use the Covariance Matrix approach (FALSE)
  bool _flagSPDE;
  Db* _db;
};
