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
#include "Model/ModelOptim.hpp"

class Model;
class Vario;

/**
 * \brief
 * Class which, starting from an experimental variogram, enables fitting the
 * various parameters of a Covariance part of a Model
 */
class GSTLEARN_EXPORT ModelOptimVario: public ModelOptim
{
public:
  ModelOptimVario();
  ModelOptimVario(const ModelOptimVario& m);
  ModelOptimVario& operator=(const ModelOptimVario& m);
  virtual ~ModelOptimVario();

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
    // Pointer to the Vario structure
    Vario* _vario;

    // Parametrization
    int _wmode;

    // Experimental quantities
    std::vector<OneLag> _lags;

  } Vario_Part;

  int fit(Vario* vario, Model* model, int wmode = 2, bool verbose = false);

  static double evalCost(unsigned int nparams,
                         const double* current,
                         double* grad,
                         void* my_func_data);

private:
  int _buildExperimental();
  int _getTotalLagsPerDirection() const;
  int _getParamNumber() const { return (int)_modelPart._params.size(); }

  VectorDouble _computeWeightPerDirection();

  void _copyVarioPart(const Vario_Part& varioPart);
  bool _checkConsistency();
  bool _isLagCorrect(int idir, int k) const;

  OneLag _createOneLag(int ndim, int idir, int ivar, int jvar, double gg, double dist) const;

protected:
  void _computeWt();
  double _getC00(int idir, int ivar, int jvar) const;

protected:
  // Part relative to the Experimental variograms
  Vario_Part _varioPart;

  VectorDouble _wt;
};
