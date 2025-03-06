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
  _cova = _covList = nullptr;
}

void ModelCovList::setCovList(CovList* covs)
{
  _covList = covs;
  _cova    = _covList;
}

ModelCovList::~ModelCovList() {}

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
  if (_covList == nullptr)
  {
    messerr("Error: Covariance List is nullptr");
    return;
  }
  _covList->addCov(cov);
}
