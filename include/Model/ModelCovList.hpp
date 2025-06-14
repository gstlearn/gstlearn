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

#include "Covariances/CovList.hpp"
#include "Model/ModelGeneric.hpp"
#include "geoslib_define.h"
#include "gstlearn_export.hpp"

/**
 * \brief
 * Class containing the ModelCovList Information describing the formal Spatial (or Temporal) Characteristics
 * of the (set of) random variable(s) under study.
 *
 * The ModelCovList is essentially a container with two main contents:
 * - the **covariance** part: see CovList.hpp for more information
 */
class GSTLEARN_EXPORT ModelCovList : public ModelGeneric
{
public:
  ModelCovList(const CovContext& ctxt = CovContext());
  ModelCovList(const ModelCovList &m);
  ModelCovList& operator= (const ModelCovList &m);
  virtual ~ModelCovList();

  const CovList* getCovList() const { return (const CovList*)getCov(); }
  CovList* getCovListModify() const { return  (CovList*)getCov(); }

  FORWARD_METHOD_NON_CONST(getCovListModify, delCov)
  FORWARD_METHOD_NON_CONST(getCovListModify, delAllCov)
  FORWARD_METHOD_NON_CONST(getCovListModify, setCovFiltered);
  FORWARD_METHOD_NON_CONST(getCovListModify, makeSillNoStatDb)
  FORWARD_METHOD_NON_CONST(getCovListModify, makeSillStationary)
  FORWARD_METHOD_NON_CONST(getCovListModify, makeSillsStationary)
  FORWARD_METHOD_NON_CONST(getCovListModify, makeSillNoStatFunctional)
  FORWARD_METHOD_NON_CONST(getCovListModify, setSill)
  FORWARD_METHOD_NON_CONST(getCovListModify, setSills)
  FORWARD_METHOD_NON_CONST(getCovListModify, normalize)

  FORWARD_METHOD(getCovList, getNCov)
  FORWARD_METHOD(getCovList, getSills)
  FORWARD_METHOD(getCovList, getSill, TEST) 
  FORWARD_METHOD(getCovList, getTotalSill)
  FORWARD_METHOD(getCovList, getTotalSills)
  FORWARD_METHOD(getCovList, isAllActiveCovList)
  
  void setCovList(CovList* covs);
  virtual void addCov(const CovBase* cov);

public:
  mutable AModelFitSills* _modelFitSills; /* Model fitting procedure for Sills */
};
