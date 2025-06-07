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
class Vario;
class Constraints;
class MatrixDense;
class MatrixSymmetric;

/**
 * \brief
 * Class which, starting from an experimental variogram, enables fitting the
 * sills of all Covariance parts of a Model
 */
class GSTLEARN_EXPORT AModelFitSillsVario: public AModelFitSills
{
public:
  AModelFitSillsVario(Vario* vario,
                      Model* model,
                      Constraints* constraints      = nullptr,
                      const Option_AutoFit& mauto   = Option_AutoFit(),
                      const Option_VarioFit& optvar = Option_VarioFit());
  AModelFitSillsVario(const AModelFitSillsVario& m);
  AModelFitSillsVario& operator=(const AModelFitSillsVario& m);
  virtual ~AModelFitSillsVario();

  int fit(bool verbose = false);
  void updateFromModel();

private:
  int _prepare();
  int _getDimensions();
  void _computeGg();
  void _compressArray(const VectorDouble& tabin, VectorDouble& tabout);
  void _prepareGoulard();

private:
  Vario* _vario;
  CovCalcMode _calcmode;
};
