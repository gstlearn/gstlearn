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

#include "geoslib_define.h"
#include "gstlearn_export.hpp"

#include "Enum/ECov.hpp"
#include "Covariances/ACov.hpp"
#include "Covariances/CovContext.hpp"
#include "Matrix/MatrixSquareSymmetric.hpp"
#include "Drifts/DriftList.hpp"
#include "Basic/ICloneable.hpp"

class Model;
class Db;
class CovInternal;
class CovCalcMode;
class CovAnisoList;
class Vario;
class ADrift;
class AnamContinuous;
class AnamHermite;


typedef std::vector<ECov> VectorECov;

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
  ACov*             getCovModify()             { return  _cova;     }
  CovContext*       getContextModify()         { return &_ctxt;     }
  DriftList*        getDriftListModify()       { return  _driftList;}
  
  // Forwarding the methods from _cova
  FORWARD_METHOD(getCov, evalCovMatOptim)
  FORWARD_METHOD(getCov, evalCovMatSym)
  FORWARD_METHOD(getCov, evalCovMatSymOptim)
  FORWARD_METHOD(getCov, evalCovMatOptimByRanks)
  FORWARD_METHOD(getCov, evalCovMatSymOptimByRanks)
  FORWARD_METHOD(getCov, evalCovMatSymByRanks)
  FORWARD_METHOD(getCov, eval0Mat)
  FORWARD_METHOD(getCov, evalCovMat)
  FORWARD_METHOD(getCov, evalCovMatByRanks)
  FORWARD_METHOD(getCov, evalCovMatSparse)
  
  // Forwarding the methods from _driftList
  FORWARD_METHOD(getDriftList, evalDriftMat)
  FORWARD_METHOD(getDriftList, evalDriftMatByRanks)
  FORWARD_METHOD(getDriftList, evalDriftMatByTarget)

  // Forwarding the methods from _ctxt
  FORWARD_METHOD(getContext, getNVar)
  FORWARD_METHOD(getContext, getNDim)
  FORWARD_METHOD(getContext, getMeans)

  void setField(double field);
  bool isValid() const;

  // Case of _driftList
  int getNDrift() const;
  int getNDriftEquation() const;
  int getNExtDrift() const;
  int getDriftMaxIRFOrder(void) const;
  void delAllDrifts();

    // Case of _cova
  const CovAnisoList* getCovAnisoList() const;
  CovAnisoList* getCovAnisoListModify() const;
  int getCovaMinIRFOrder() const;
  int getNCov(bool skipNugget = false) const;
  void setActiveFactor(int iclass);
  int  getActiveFactor() const;

protected:               // TODO : pass into private to finish clean
  ACov* _cova;           /* Generic Covariance structure */
  DriftList* _driftList; /* Series of Drift functions */
  CovContext _ctxt;      /* Context */
};
