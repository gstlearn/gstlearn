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
#include "Model/AModelOptimSills.hpp"
#include "Model/Option_AutoFit.hpp"
#include "Model/Option_VarioFit.hpp"

class Model;
class DbGrid;
class Constraints;
class MatrixRectangular;
class MatrixSquareSYmmetric;

/**
 * \brief
 * Class which, starting from an experimental variogram, enables fitting the
 * sills of all Covariance parts of a Model
 */
class GSTLEARN_EXPORT ModelOptimSillsVMap: public AModelOptimSills
{
public:
  ModelOptimSillsVMap(Model* model,
                      Constraints* constraints      = nullptr,
                      const Option_AutoFit& mauto   = Option_AutoFit(),
                      const Option_VarioFit& optvar = Option_VarioFit());
  ModelOptimSillsVMap(const ModelOptimSillsVMap& m);
  ModelOptimSillsVMap& operator=(const ModelOptimSillsVMap& m);
  virtual ~ModelOptimSillsVMap();

  int fit(DbGrid* dbmap);
  int loadEnvironment(DbGrid* dbmap);

private:
  int  _getDimensions();
  void _computeVMap();
  void _updateFromModel();

private:
  DbGrid* _dbmap;
  VectorInt _indg1;
  VectorInt _indg2;
  int _nech;
};
