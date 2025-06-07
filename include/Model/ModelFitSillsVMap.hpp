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
#include "Model/Option_AutoFit.hpp"
#include "Model/Option_VarioFit.hpp"

class Model;
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
                    Constraints* constraints      = nullptr,
                    const Option_AutoFit& mauto   = Option_AutoFit(),
                    const Option_VarioFit& optvar = Option_VarioFit());
  ModelFitSillsVMap(const ModelFitSillsVMap& m);
  ModelFitSillsVMap& operator=(const ModelFitSillsVMap& m);
  virtual ~ModelFitSillsVMap();

  int fit(bool verbose = false);
  void updateFromModel();

private:
  int  _prepare();
  int  _getDimensions();
  void _computeVMap();

private:
  const DbGrid* _dbmap;
  VectorInt _indg1;
  VectorInt _indg2;
  int _nech;
  CovCalcMode _calcmode;
};
