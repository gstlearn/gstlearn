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
#include "Model/ModelFitSillsVario.hpp"
#include "Model/ModelFitSillsVMap.hpp"
#include "Covariances/CovBase.hpp"

namespace gstlrn
{
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

void ModelCovList::setCovList(const CovList* covs)
{
  setCov(covs);
}

ModelCovList::~ModelCovList() 
{

}

void ModelCovList::addCov(const CovBase& cov)
{
  if (!cov.getContext().isEqual(_ctxt))
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

void ModelCovList::fitSills(Vario* vario,
                            const DbGrid* dbmap,
                            Constraints* constraints,
                            const ModelOptimParam& mop,
                            bool verbose,
                            bool trace)
{
  if (vario != nullptr)
  {
    setFitSills(ModelFitSillsVario::createForOptim(vario, this, constraints, mop));
  }
  else if (dbmap != nullptr)
  {
    setFitSills(ModelFitSillsVMap::createForOptim(dbmap, this, constraints, mop));
  }

  AModelFitSills* amf = getFitSills();
  if (amf == nullptr) return;

  amf->setVerbose(verbose);
  amf->setTrace(trace);

  _cova->updateCov();

  // Cancel the structure possibly used for Goulard (to be improved)
  deleteFitSills();
}
}