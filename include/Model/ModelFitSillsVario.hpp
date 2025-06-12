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
class Vario;
class Constraints;
class MatrixDense;
class MatrixSymmetric;

/**
 * \brief
 * Class which, starting from an experimental variogram, enables fitting the
 * sills of all Covariance parts of a Model
 */
class GSTLEARN_EXPORT ModelFitSillsVario: public AModelFitSills
{
public:
  ModelFitSillsVario(Vario* vario,
                     ModelCovList* model,
                     Constraints* constraints   = nullptr,
                     const ModelOptimParam& mop = ModelOptimParam());
  ModelFitSillsVario(const ModelFitSillsVario& m);
  ModelFitSillsVario& operator=(const ModelFitSillsVario& m);
  virtual ~ModelFitSillsVario();

  IMPLEMENT_CLONING(ModelFitSillsVario)

  int fitSills(bool verbose = false) override;

  static ModelFitSillsVario* createForOptim(Vario* vario,
                                            ModelGeneric* model,
                                            Constraints* constraints   = nullptr,
                                            const ModelOptimParam& mop = ModelOptimParam());

private:
  int _prepare();
  int _getDimensions();
  void _computeGg();
  void _compressArray(const VectorDouble& tabin, VectorDouble& tabout);
  void _prepareGoulard();
  void _updateFromModel();

private:
  Vario* _vario;
};
