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
  } Model_Part;

protected:
  int _buildModelParamList();
  int _getParamNumber() const { return (int) _modelPart._params.size(); }

  static void _patchModel(Model_Part& modelPart, const double* current);
  void _copyModelPart(const Model_Part& modelPart);

  void _addOneModelParam(int icov, const EConsElem& type, int rank, double lbound, double ubound);
  
protected:
  // Part of the structure dedicated to the Model
  Model_Part _modelPart;
};
