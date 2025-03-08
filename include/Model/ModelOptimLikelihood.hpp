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

  int fit(Db* db, bool flagSPDE = false, bool verbose = false);
  int loadEnvironment(Db* db, bool flagSPDE = false, bool verbose = false);

#ifndef SWIG
  static double evalCost(unsigned int nparams,
                         const double* current,
                         double* grad,
                         void* my_func_data);
#endif

private:
  typedef struct
  {
    // If TRUE: use SPDE approach
    // otherwise: use the Covariance Matrix approach
    bool _flagSPDE;

    // Pointer to the Vario structure
    Db* _db;
  } Db_Part;

  typedef struct
  {
    // Part of the structure dedicated to the Model
    Model_Part& _modelPart;

    // Part relative to the Experimental variograms
    ModelOptimLikelihood::Db_Part& _dbPart;

  } AlgorithmLikelihood;

  void _copyDbPart(const Db_Part& dbPart);
  bool _checkConsistency();

private:
  // Part relative to the Db
  Db_Part _dbPart;
};
