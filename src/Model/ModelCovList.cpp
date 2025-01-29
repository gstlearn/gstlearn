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


ModelCovList::ModelCovList(const CovContext &ctxt)
: ModelGeneric(ctxt)
{
    _cova = _covList = nullptr;
}


void ModelCovList::delAllCovas()
{
  if (_covList == nullptr) return;
  _covList->delAllCov();
}

const MatrixSquareSymmetric& ModelCovList::getSillValues(int icov) const
{
  if (_cova == nullptr) return _dummy;
  return _covList->getSill(icov);
}

double ModelCovList::getSill(int icov, int ivar, int jvar) const
{
  if (_covList == nullptr) return TEST;
  return _covList->getSill(icov, ivar, jvar);
}

double ModelCovList::getTotalSill(int ivar, int jvar) const
{
  return _covList->getTotalSill(ivar, jvar);
}

void ModelCovList::setCovList(CovList* covs)
{
    _covList = covs;
    _cova = _covList;
}

VectorInt ModelCovList::getActiveCovList() const
{
  if (_covList == nullptr) return VectorInt();
  return _covList->getActiveCovList();
}
VectorInt ModelCovList::getAllActiveCovList() const
{
  if (_covList == nullptr) return VectorInt();
  return _covList->getAllActiveCovList();
}

bool ModelCovList::isAllActiveCovList() const
{
  if (_covList == nullptr) return false;
  return _covList->isAllActiveCovList();
}

MatrixSquareSymmetric ModelCovList::getTotalSills() const
{
  return _covList->getTotalSills();
}

ModelCovList::~ModelCovList()
{
}
