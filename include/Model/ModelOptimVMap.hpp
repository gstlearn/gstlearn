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

#include "Model/AModelOptimNew.hpp"
#include "gstlearn_export.hpp"

#include "Basic/VectorNumT.hpp"
#include "Model/AModelOptim.hpp"
#include "Model/ModelOptimSillsVMap.hpp"
#include "Model/Option_AutoFit.hpp"
#include "Model/Option_VarioFit.hpp"

class Model;
class DbGrid;
class Constraints;

/**
 * \brief
 * Class which, starting from an experimental variogram, enables fitting the
 * various parameters of a Covariance part of a Model
 */
class GSTLEARN_EXPORT ModelOptimVMap: public AModelOptimNew, public AModelOptim
{
public:
  ModelOptimVMap(ModelGeneric* model,
                 Constraints* constraints      = nullptr,
                 const Option_AutoFit& mauto   = Option_AutoFit(),
                 const Option_VarioFit& optvar = Option_VarioFit());
  ModelOptimVMap(const ModelOptimVMap& m);
  ModelOptimVMap& operator=(const ModelOptimVMap& m);
  virtual ~ModelOptimVMap();

  int fit(const DbGrid* dbmap, bool flagGoulard = true, bool verbose = false);

  int loadEnvironment(const DbGrid* dbmap, bool flagGoulard = true, bool verbose = false);

  double computeCost(bool verbose = false) override;

  static ModelOptimVMap* createForOptim(ModelGeneric* model,
                                        const DbGrid* dbmap,
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
  struct VMap_Part
  {
    const DbGrid* _dbmap;
    VectorInt _indg1;
    VectorInt _indg2;
    int _npadir;
  };

  struct AlgorithmVMap
  {
    Model_Part& _modelPart;
    VMap_Part& _vmapPart;
    ModelOptimSillsVMap& _goulardPart;
  };

private:
  void _copyVMapPart(const VMap_Part& vmapPart);
  bool _checkConsistency();
  int  _getDimensions();
  void _allocateInternalArrays();
  void _computeFromVMap();

protected:
  // Model fitting options
  Option_VarioFit _optvar;

  // Model fitting parameters
  Option_AutoFit _mauto;

  // Set of constraints
  Constraints* _constraints;
  CovCalcMode _calcmode;

  // Part relative to the Experimental VMap
  VMap_Part _vmapPart;

  // Only used for Goulard Option
  ModelOptimSillsVMap _goulardPart;

  // Following members are simply there to accelerate the computation
  int _ndim;
  int _nvar;
  int _nech;
};
