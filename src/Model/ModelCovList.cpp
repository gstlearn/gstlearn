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
#include "Model/ModelCovList.hpp"
#include "Covariances/CovBase.hpp"

ModelCovList::ModelCovList(const CovContext& ctxt)
  : ModelGeneric(ctxt)
  , _modelFitSills(nullptr)
{
  _cova = nullptr;
}

ModelCovList::ModelCovList(const ModelCovList &m)
  : ModelGeneric(m)
{
  _modelFitSills = (m._modelFitSills != nullptr) ? (AModelFitSills*)m._modelFitSills->clone() : nullptr;

}
ModelCovList& ModelCovList:: operator= (const ModelCovList &m)
{
  if (this != &m)
  {
    ModelGeneric::operator=(m);
    _modelFitSills = (m._modelFitSills != nullptr) ? (AModelFitSills*)m._modelFitSills->clone() : nullptr;

  }
  return *this;
}

void ModelCovList::setCovList(CovList* covs)
{
  setCov(covs);
}

ModelCovList::~ModelCovList() 
{
  delete _modelFitSills;
  _modelFitSills = nullptr;
}

void ModelCovList::addCov(const CovBase* cov)
{
  if (cov == nullptr)
  {
    messerr("Error: Covariance is nullptr");
    return;
  }

  if (!cov->getContext().isEqual(_ctxt))
  {
    messerr("Error: Covariance should share the same Context as 'Model'");
    messerr("Operation is cancelled");
    return;
  }
  if (getCovList() == nullptr)
  {
    messerr("Error: Covariance List is nullptr");
    return;
  }
  getCovListModify()->addCov(cov);
}
