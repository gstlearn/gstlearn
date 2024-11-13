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

#include "Basic/VectorNumT.hpp"
#include "Model/ModelOptim.hpp"

class Model;
class Db;

/**
 * \brief
 * Class which, starting from an experimental variogram, enables fitting the
 * various parameters of a Covariance part of a Model
 */
class GSTLEARN_EXPORT ModelOptimLikelihood: public ModelOptim
{
public:
  ModelOptimLikelihood();
  ModelOptimLikelihood(const ModelOptimLikelihood& m);
  ModelOptimLikelihood& operator=(const ModelOptimLikelihood& m);
  virtual ~ModelOptimLikelihood();

  typedef struct
  {
    // Pointer to the Vario structure
    const Db* _db;

  } Db_Part;

  typedef struct
  {
    // Part of the structure dedicated to the Model
    Model_Part& _modelPart;

    // Part relative to the Experimental variograms
    Db_Part& _dbPart;

  } AlgorithmLikelihood;

  int fit(const Db* db, Model* model, bool verbose = false);

  static double evalCost(unsigned int nparams,
                         const double* current,
                         double* grad,
                         void* my_func_data);

private:
  void _copyDbPart(const Db_Part& dbPart);
  bool _checkConsistency();

private:
  // Part relative to the Db
  Db_Part _dbPart;
};
