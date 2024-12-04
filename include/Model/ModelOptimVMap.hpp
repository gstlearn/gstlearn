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
class GSTLEARN_EXPORT ModelOptimVMap: public AModelOptim
{
public:
  ModelOptimVMap(Model* model,
                 Constraints* constraints      = nullptr,
                 const Option_AutoFit& mauto   = Option_AutoFit(),
                 const Option_VarioFit& optvar = Option_VarioFit());
  ModelOptimVMap(const ModelOptimVMap& m);
  ModelOptimVMap& operator=(const ModelOptimVMap& m);
  virtual ~ModelOptimVMap();

  int fit(const DbGrid* dbmap, bool flagGoulard = true, bool verbose = false);

  int loadEnvironment(const DbGrid* dbmap, bool flagGoulard = true, bool verbose = false);

#ifndef SWIG
  static double evalCost(unsigned int nparams,
                         const double* current,
                         double* grad,
                         void* my_func_data);
#endif

private:
  typedef struct
  {
    const DbGrid* _dbmap;
    VectorInt _indg1;
    VectorInt _indg2;
    int _npadir;
  } VMap_Part;

  typedef struct
  {
    // Part of the structure dedicated to the Model
    Model_Part& _modelPart;

    // Part relative to the Experimental variograms
    VMap_Part& _vmapPart;

    // Only used for Goulard Option
    ModelOptimSillsVMap& _goulardPart;
  } AlgorithmVMap;

  void _copyVMapPart(const VMap_Part& vmapPart);
  bool _checkConsistency();
  int  _getDimensions();
  void _allocateInternalArrays();
  void _computeFromVMap();

protected:
  // Part relative to the Experimental VMap
  VMap_Part _vmapPart;

  // Only used for Goulard Option
  ModelOptimSillsVMap _optGoulard;
};
