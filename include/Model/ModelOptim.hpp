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
#include <vector>

class Model;
class MatrixSquareSymmetric;

/**
 * \brief
 * Class which, starting from an experimental variogram, enables fitting the
 * various parameters of a Covariance part of a Model
 */
class GSTLEARN_EXPORT ModelOptim
{
public:
  ModelOptim();
  ModelOptim(const ModelOptim& m);
  ModelOptim& operator=(const ModelOptim& m);
  virtual ~ModelOptim();

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

    // Model parametrization
    std::vector<OneParam> _params;

    // Model parametrization
    VectorDouble _tabval;
    VectorDouble _tablow;
    VectorDouble _tabupp;

    // Verbosity flag
    bool _verbose;
    int  _niter;
  } Model_Part;

protected:
  int _buildModelParamList();
  int _getParamNumber() const { return (int) _modelPart._params.size(); }
  void updateModelParamList(double hmax, const MatrixSquareSymmetric& vars);
  void dumpParamList() const;
  static void _patchModel(Model_Part& modelPart, const double* current);
  static void _printResult(const String& title, const Model_Part& modelPart, double result);
  void _setSill(int icov, int ivar, int jvar, double value) const;

private:
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
};
