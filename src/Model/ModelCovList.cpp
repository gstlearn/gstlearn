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
