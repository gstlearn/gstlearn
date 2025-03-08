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

#include "Enum/EConsElem.hpp"
#include "gstlearn_export.hpp"

#include "Basic/VectorNumT.hpp"
#include "Covariances/CovCalcMode.hpp"
#include "Matrix/MatrixSquareSymmetric.hpp"
#include "Model/Option_AutoFit.hpp"
#include "Model/Option_VarioFit.hpp"
#include <vector>

class Model;
class Constraints;

typedef struct
{
  int _icov;
  EConsElem _type;
  int _rank;
  double _scale;
} OneParam;

typedef struct
{
  // Pointer to the Model structure
  Model* _model;

  // Model fitting options
  Option_VarioFit _optvar;

  // Model parametrization
  std::vector<OneParam> _params;

  // Model parametrization
  VectorDouble _tabval;
  VectorDouble _tablow;
  VectorDouble _tabupp;

  // Verbosity flag
  bool _verbose;
  int _niter;
  CovCalcMode _calcmode;
} Model_Part;

/**
 * \brief
 * Class which, starting from an experimental variogram, enables fitting the
 * various parameters of a Covariance part of a Model
 */
class GSTLEARN_EXPORT AModelOptim
{
public:
  AModelOptim(Model* model,
              Constraints* constraints      = nullptr,
              const Option_AutoFit& mauto   = Option_AutoFit(),
              const Option_VarioFit& optvar = Option_VarioFit());
  AModelOptim(const AModelOptim& m);
  AModelOptim& operator=(const AModelOptim& m);
  virtual ~AModelOptim();



protected:
  int _buildModelParamList();
  int _getNParam() const { return (int) _modelPart._params.size(); }

  static void _patchModel(Model_Part& modelPart, const double* current);
  static void _printResult(const String& title, const Model_Part& modelPart, double result);
  void _setSill(int icov, int ivar, int jvar, double value) const;

  void _performOptimization(double (*optim_func)(unsigned n,
                                                 const double* x,
                                                 double* gradient,
                                                 void* func_data),
                            void* f_data,
                            double distmax_def = TEST,
                            const MatrixSquareSymmetric& vars_def = MatrixSquareSymmetric());

private:
  void _updateModelParamList(double distmax_def = TEST,
                             const MatrixSquareSymmetric& vars_def = MatrixSquareSymmetric());
  void _dumpParamList() const;
  static void _dumpOneModelParam(const OneParam& param, double value);
  void _addOneModelParam(int icov,
                         const EConsElem& type,
                         int rank,
                         double lbound = TEST,
                         double ubound = TEST);
  void _copyModelPart(const Model_Part& modelPart);

protected:
  // Part of the structure dedicated to the Model
  Model_Part _modelPart;

  Constraints* _constraints;
  Option_AutoFit _mauto;

};
