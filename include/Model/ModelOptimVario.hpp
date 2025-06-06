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
#include "Model/AModelOptim.hpp"
#include "Model/AModelOptimNew.hpp"
#include "Model/ModelOptimSillsVario.hpp"
#include "Model/Option_AutoFit.hpp"
#include "Model/Option_VarioFit.hpp"

class ModelGeneric;
class Vario;

/**
 * \brief
 * Class which, starting from an experimental variogram, enables fitting the
 * various parameters of a Covariance part of a Model
 */
class GSTLEARN_EXPORT ModelOptimVario: public AModelOptimNew, public AModelOptim
{
public:
  ModelOptimVario(ModelGeneric* model,
                  Constraints* constraints      = nullptr,
                  const Option_AutoFit& mauto   = Option_AutoFit(),
                  const Option_VarioFit& optvar = Option_VarioFit());
  ModelOptimVario(const ModelOptimVario& m);
  ModelOptimVario& operator=(const ModelOptimVario& m);
  virtual ~ModelOptimVario();

  int fit(Vario* vario, bool flagGoulard = true, int wmode = 2, bool verbose = false);

  int loadEnvironment(Vario* vario,
                      bool flagGoulard = true,
                      int wmode        = 2,
                      bool verbose     = false);
  double computeCost(bool verbose = false) override;

  static ModelOptimVario* createForOptim(ModelGeneric* model,
                                         Vario* vario,
                                         Constraints* constraints      = nullptr,
                                         const Option_AutoFit& mauto   = Option_AutoFit(),
                                         const Option_VarioFit& optvar = Option_VarioFit());

#ifndef SWIG
  static double evalCost(unsigned int nparams,
                         const double* current,
                         double* grad,
                         void* my_func_data);
#endif

protected:
  struct OneLag
  {
    int _ivar;
    int _jvar;
    double _weight;
    double _gg;
    SpacePoint _P;
  };
  struct Vario_Part
  {
    const Vario* _vario;
    int _wmode;
    std::vector<OneLag> _lags;
  };
  struct AlgorithmVario
  {
    Model_Part& _modelPart;
    Vario_Part& _varioPart;

    // Part relative to Sill fitting procedure
    ModelOptimSillsVario& _goulardPart;
  };

private:
  int  _buildExperimental();
  void _copyVarioPart(const Vario_Part& varioPart);
  bool _checkConsistency();
  OneLag _createOneLag(int ndim, int idir, int ivar, int jvar, double gg, double dist) const;

protected:
  // Model fitting options
  Option_VarioFit _optvar;

  // Model fitting parameters
  Option_AutoFit _mauto;

  // Set of constraints
  Constraints* _constraints;
  CovCalcMode _calcmode;

  // Part relative to the Experimental variograms
  Vario_Part _varioPart;

  // Only used for Goulard Option
  ModelOptimSillsVario _goulardPart;
};
