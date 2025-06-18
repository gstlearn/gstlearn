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

#include "Space/SpacePoint.hpp"
#include "Estimation/AModelOptim.hpp"

class ModelGeneric;
class Vario;

/**
 * \brief
 * Class which, starting from an experimental variogram, enables fitting the
 * various parameters of a Covariance part of a Model
 */
class GSTLEARN_EXPORT ModelOptimVario: public AModelOptim
{
public:
  ModelOptimVario(ModelGeneric* model,
                  const Constraints* constraints   = nullptr,
                  const ModelOptimParam& mop = ModelOptimParam());
  ModelOptimVario(const ModelOptimVario& m);
  ModelOptimVario& operator=(const ModelOptimVario& m);
  virtual ~ModelOptimVario();

  double computeCost(bool verbose = false) override;
  double computeDerivatives(std::vector<double>& params);

  static ModelOptimVario* createForOptim(ModelGeneric* model,
                                         const Vario* vario,
                                         const Constraints* constraints   = nullptr,
                                         const ModelOptimParam& mop = ModelOptimParam());

protected:
  struct OneLag
  {
    int _ivar;
    int _jvar;
    double _weight;
    double _gg;
    SpacePoint _P;
  };

private:
  void _updateGradients() override;
  int  _buildExperimental();
  bool _checkConsistency();
  OneLag _createOneLag(int ndim, int idir, int ivar, int jvar, double gg, double dist) const;

protected:
  // Model fitting options
  ModelOptimParam _mop;

  // Set of constraints
  const Constraints* _constraints;

  // Calculation option
  CovCalcMode _calcmode;

  // Part relative to the Experimental variograms
  const Vario* _vario;
  std::vector<OneLag> _lags;
};
