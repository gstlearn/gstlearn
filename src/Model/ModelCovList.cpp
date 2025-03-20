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
{
  _cova = nullptr;
}

ModelCovList::ModelCovList(const ModelCovList &m)
  : ModelGeneric(m)
{

}
ModelCovList& ModelCovList:: operator= (const ModelCovList &m)
{
  if (this != &m)
  {
    ModelGeneric::operator=(m);
  }
  return *this;
}

void ModelCovList::setCovList(CovList* covs)
{
  setCov(covs);
}

ModelCovList::~ModelCovList() 
{

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
