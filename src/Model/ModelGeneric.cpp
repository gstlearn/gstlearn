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
#include "Model/ModelGeneric.hpp"
#include "Covariances/CovAnisoList.hpp"
#include "Covariances/CovLMCAnamorphosis.hpp"

ModelGeneric::ModelGeneric(const CovContext &ctxt)
    : _cova(nullptr),
      _driftList(nullptr),
      _ctxt(ctxt)
{
}

ModelGeneric::~ModelGeneric()
{
}


void ModelGeneric::setField(double field)
{
  _ctxt.setField(field);
  CovAnisoList* covalist = dynamic_cast<CovAnisoList*>(_cova);
  if (covalist != nullptr) covalist->copyCovContext(_ctxt);
  if (_driftList != nullptr) _driftList->copyCovContext(_ctxt);
}

// Pipes method to _ACov
const CovAnisoList* ModelGeneric::getCovAnisoList() const
{
  if (_cova == nullptr) return nullptr;
  const CovAnisoList* covalist = dynamic_cast<const CovAnisoList*>(_cova);
  return covalist;
}
CovAnisoList* ModelGeneric::getCovAnisoListModify() const
{
  if (_cova == nullptr) return nullptr;
  CovAnisoList* covalist = dynamic_cast<CovAnisoList*>(_cova);
  return covalist;
}
int ModelGeneric::getCovaMinIRFOrder() const
{
  const CovAnisoList* covalist = getCovAnisoList();
  if (covalist == nullptr) return ITEST;
  return covalist->getCovaMinIRFOrder();
}
int ModelGeneric::getNCov(bool skipNugget) const
{
  const CovAnisoList* covalist = getCovAnisoList();
  if (covalist == nullptr) return ITEST;
  return covalist->getNCov(skipNugget);
}
void ModelGeneric::setActiveFactor(int iclass)
{
  if (_cova == nullptr) return;
  CovLMCAnamorphosis* covalist = dynamic_cast<CovLMCAnamorphosis*>(_cova);
  if (covalist == nullptr) return;
  covalist->setActiveFactor(iclass);
}
int ModelGeneric::getActiveFactor() const
{
  if (_cova == nullptr) return 0;
  CovLMCAnamorphosis* covalist = dynamic_cast<CovLMCAnamorphosis*>(_cova);
  if (covalist == nullptr) return 0;
  return covalist->getActiveFactor();
}

bool ModelGeneric::isValid() const
{
  return _isValid();
}

bool ModelGeneric::_isValid() const
{
  return true;
}