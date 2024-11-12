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
#include "Space/SpacePoint.hpp"

class Model;
class Vario;

typedef struct
{
  int _ivar;
  int _jvar;
  double _weight;
  double _gg;
  SpacePoint _P;
} OneLag;

typedef struct
{
  int _icov;
  int _type;
  int _rank;
} OneParam;

typedef struct
{
  // Pointer to the Model structure (updated by algorithm()
  Model* _model;

  // Experimental quantities
  std::vector<OneLag> _lags;

  // Model parametrization
  std::vector<OneParam> _params;
} Algorithm;

/**
 * \brief
 * Class which, starting from an experimental variogram, enables fitting the
 * various parameters of a Covariance part of a Model
 */
class GSTLEARN_EXPORT Model_Optim
{
public:
  Model_Optim(const Vario* vario = nullptr);
  Model_Optim(const Model_Optim& m);
  Model_Optim& operator=(const Model_Optim& m);
  virtual ~Model_Optim();

  int fit(Model* model, int wmode = 2);

  static double evalCost(unsigned int nparams,
                         const double* current,
                         double* grad,
                         void* my_func_data);

  Model* getModel();

private:
  int _buildExperimental();
  int _buildModelParamList();
  bool _isLagCorrect(int idir, int k);
  double _getC00(int idir, int ivar, int jvar);
  VectorDouble _computeWeight();
  VectorDouble _computeWeightPerDirection();
  int _getTotalLagsPerDirection() const;

  static void _patchModel(Algorithm* algo, const double* current);

  OneLag
  _createOneLag(int ndim, int idir, int ivar, int jvar, double gg, double dist);
  void
  _addOneModelParam(int icov, int type, int rank, double lbound, double ubound);

private:
  const Vario* _vario;
  int _wmode;

  // Model parametrization
  VectorDouble _tabval;
  VectorDouble _tablow;
  VectorDouble _tabupp;

  // Minimization algorithm
  Algorithm _algorithm;
};
