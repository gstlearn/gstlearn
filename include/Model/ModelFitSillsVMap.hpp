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
#include "Model/AModelFitSills.hpp"
#include "Model/ModelOptimParam.hpp"

class ModelGeneric;
class DbGrid;
class Constraints;
class MatrixDense;
class MatrixSymmetric;

/**
 * \brief
 * Class which, starting from an experimental variogram, enables fitting the
 * sills of all Covariance parts of a Model
 */
class GSTLEARN_EXPORT ModelFitSillsVMap: public AModelFitSills
{
public:
  ModelFitSillsVMap(const DbGrid* dbmap,
                    ModelCovList* model,
                    Constraints* constraints   = nullptr,
                    const ModelOptimParam& mop = ModelOptimParam());
  ModelFitSillsVMap(const ModelFitSillsVMap& m);
  ModelFitSillsVMap& operator=(const ModelFitSillsVMap& m);
  virtual ~ModelFitSillsVMap();

  IMPLEMENT_CLONING(ModelFitSillsVMap)

  int fitSills(bool verbose = false, bool trace = false) override;

  static ModelFitSillsVMap* createForOptim(const DbGrid* dbmap,
                                           ModelGeneric* model,
                                           Constraints* constraints   = nullptr,
                                           const ModelOptimParam& mop = ModelOptimParam());

private:
  int _prepare();
  int  _getDimensions();
  void _computeVMap();
  void _updateFromModel();

private:
  const DbGrid* _dbmap;
  VectorInt _indg1;
  VectorInt _indg2;
  int _nech;
};
