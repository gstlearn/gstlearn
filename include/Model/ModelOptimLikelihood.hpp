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
    // If TRUE: use SPDE approach
    // otherwise: use the Covariance Matrix approach
    bool _flagSPDE;

    // Pointer to the Vario structure
    Db* _db;

  } Db_Part;

  int fit(Db* db, Model* model, bool flagSPDE = false, bool verbose = false);

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
