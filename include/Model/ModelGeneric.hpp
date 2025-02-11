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

#include "Basic/ICloneable.hpp"
#include "geoslib_define.h"
#include "gstlearn_export.hpp"

#include "Covariances/ACov.hpp"
#include "Covariances/CovContext.hpp"
#include "Drifts/DriftList.hpp"

class Model;
class Db;
class DbGrid;
class CovCalcMode;
/**
 * \brief
 * Class containing the Model Information describing the formal Spatial (or Temporal) Characteristics
 * of the (set of) random variable(s) under study.
 *
 * The Model is essentially a container with two main contents:
 * - the **covariance** part: see ACov.hpp for more information
 * - the **drift** part: see DriftList.hpp for more information
 *
 * The additional member **CovContext** only serves in carrying the following information:
 * - the number of variables: if more than 1, the Model becomes multivariate
 * - the field extension: this information is needed to get a *stationary* version to any covariance
 * - the experimental mean vector and the variance-covariance matrix (used to calibrate the Model)
 */
class GSTLEARN_EXPORT ModelGeneric : public ICloneable
{
public:
  ModelGeneric(const CovContext& ctxt = CovContext());
  ModelGeneric(const ModelGeneric &m) = delete;
  ModelGeneric& operator= (const ModelGeneric &m) = delete;
  virtual ~ModelGeneric();

  //getters for member pointers
  const ACov*       getCov()             const { return  _cova;     }
  const CovContext* getContext()         const { return &_ctxt;     }
  const DriftList*  getDriftList()       const { return  _driftList;}

public:
  ACov* _getCovModify() { return _cova; }
  CovContext* _getContextModify() { return &_ctxt; }
  DriftList* _getDriftListModify() { return _driftList; }

public:

  // Forwarding the methods from _cova
  FORWARD_METHOD(getCov, evalCovMatBiPointInPlace)
  FORWARD_METHOD(getCov, eval0CovMatBiPointInPlace)
  FORWARD_METHOD(getCov, evalCovMatOptim)
  FORWARD_METHOD(getCov, evalCovMatSym)
  FORWARD_METHOD(getCov, evalCovMatSymOptim)
  FORWARD_METHOD(getCov, evalCovMatOptimByTarget)
  FORWARD_METHOD(getCov, evalCovMatSymOptimByRanks)
  FORWARD_METHOD(getCov, evalCovMatSymByRanks)
  FORWARD_METHOD(getCov, eval0Mat)
  FORWARD_METHOD(getCov, eval0MatByTarget)
  FORWARD_METHOD(getCov, evalCovMat)
  FORWARD_METHOD(getCov, evalCovMatV)
  FORWARD_METHOD(getCov, evalCovMatByTarget)
  FORWARD_METHOD(getCov, evalCovMatSparse)
  FORWARD_METHOD(getCov, eval0)
  FORWARD_METHOD(getCov, eval)
  FORWARD_METHOD(getCov, evalNvarIpas)
  FORWARD_METHOD(getCov, evalMat)
  FORWARD_METHOD(getCov, evalNvarIpasIncr)
  FORWARD_METHOD(getCov, evalIvarNpas)
  FORWARD_METHOD(getCov, evalIvarIpas)
  FORWARD_METHOD(getCov, evalCvv)
  FORWARD_METHOD(getCov, evalCvvShift)
  FORWARD_METHOD(getCov, evalCvvM)
  FORWARD_METHOD(getCov, evalCxv)
  FORWARD_METHOD(getCov, evalCxvM)
  FORWARD_METHOD(getCov, evalPointToDb)
  FORWARD_METHOD(getCov, evalPointToDbAsSP)
  FORWARD_METHOD(getCov, evalAverageDbToDb,TEST)
  FORWARD_METHOD(getCov, evalAverageIncrToIncr,TEST)
  FORWARD_METHOD(getCov, evalAveragePointToDb,TEST)
  FORWARD_METHOD(getCov, samplingDensityVariance, TEST)
  FORWARD_METHOD(getCov, specificVolume, TEST)
  FORWARD_METHOD(getCov, coefficientOfVariation, TEST)
  FORWARD_METHOD(getCov, specificVolumeFromCoV, TEST)
  FORWARD_METHOD(getCov, extensionVariance, TEST)
  FORWARD_METHOD(getCov, calculateStDev, TEST)
  FORWARD_METHOD(getCov, evaluateMatInPlace)
  FORWARD_METHOD(getCov, evaluateOneGeneric, TEST)
  FORWARD_METHOD(getCov, evaluateOneIncr, TEST)
  FORWARD_METHOD(getCov, buildVmapOnDbGrid)
  FORWARD_METHOD(getCov, sample)
  FORWARD_METHOD(getCov, sampleUnitary)
  FORWARD_METHOD(getCov, envelop)
  FORWARD_METHOD(getCov, gofToVario, TEST)

  FORWARD_METHOD_NON_CONST(_getCovModify, setContext)

  // Forwarding the methods from _driftList
  
  FORWARD_METHOD(getDriftList, getDrift)
  FORWARD_METHOD(getDriftList, computeDrift, TEST)
  FORWARD_METHOD(getDriftList, evalDriftValue, TEST)
  FORWARD_METHOD(getDriftList, evalDriftMat)
  FORWARD_METHOD(getDriftList, evalDriftMatByRanks)
  FORWARD_METHOD(getDriftList, evalDriftMatByTarget)
  FORWARD_METHOD(getDriftList, getNDrift)
  FORWARD_METHOD(getDriftList, getNDriftEquation)
  FORWARD_METHOD(getDriftList, getNExtDrift)
  FORWARD_METHOD(getDriftList, isFlagLinked)
  FORWARD_METHOD(getDriftList, getDriftMaxIRFOrder,-1)
  FORWARD_METHOD(getDriftList, getRankFex)
  FORWARD_METHOD(getDriftList, isDriftSampleDefined)
  FORWARD_METHOD(getDriftList, isDriftFiltered)
  FORWARD_METHOD(getDriftList, isDriftDefined)
  FORWARD_METHOD(getDriftList, isDriftDifferentDefined)
  FORWARD_METHOD(getDriftList, getDrifts)
  FORWARD_METHOD(getDriftList, evalDrift, TEST)
  FORWARD_METHOD(getDriftList, evalDriftBySample)
  FORWARD_METHOD(getDriftList, evalDriftBySampleInPlace)
  FORWARD_METHOD(getDriftList, hasDrift, false)

  FORWARD_METHOD(getDriftList, getMean, TEST)
  FORWARD_METHOD(getDriftList, getMeans)
  FORWARD_METHOD(getDriftList, evalDriftVarCoef,TEST)
  FORWARD_METHOD(getDriftList, evalDriftVarCoefs)

  FORWARD_METHOD_NON_CONST(_getDriftListModify, setFlagLinked)
  FORWARD_METHOD_NON_CONST(_getDriftListModify, setBetaHat)
  FORWARD_METHOD_NON_CONST(_getDriftListModify, setFiltered)
  FORWARD_METHOD_NON_CONST(_getDriftListModify, delDrift)
  FORWARD_METHOD_NON_CONST(_getDriftListModify, delAllDrifts)
  FORWARD_METHOD_NON_CONST(_getDriftListModify, copyCovContext)
  FORWARD_METHOD_NON_CONST(_getDriftListModify, setMeans)
  FORWARD_METHOD_NON_CONST(_getDriftListModify, setMean)
  
  // Forwarding the methods from _ctxt
  FORWARD_METHOD(getContext, getNVar, -1)
  FORWARD_METHOD(getContext, getNDim, -1)
  FORWARD_METHOD(getContext, getSpace)

  
  FORWARD_METHOD(getContext, getCovar0)
  FORWARD_METHOD(getContext, getField, TEST)

  FORWARD_METHOD_NON_CONST(_getContextModify, setField)
  FORWARD_METHOD_NON_CONST(_getContextModify, setCovar0s)
  FORWARD_METHOD_NON_CONST(_getContextModify, setCovar0)
  
  void setField(double field);
  bool isValid() const;

  void   setDriftList(const DriftList* driftlist);
  void   setDriftIRF(int order = 0, int nfex = 0);
  void   addDrift(const ADrift* drift);  // TODO: check that the same driftM has not been already defined
  void   setDrifts(const VectorString& driftSymbols);
  
  double computeLogLikelihood(const Db* db, bool verbose = false);  

private :
  virtual bool _isValid() const;

protected:               // TODO : pass into private to finish clean
  ACov* _cova;           /* Generic Covariance structure */
  DriftList* _driftList; /* Series of Drift functions */
  CovContext _ctxt;      /* Context */
};
