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

MatrixRectangular ModelGeneric::evalDriftMat(const Db* db,
                                                int ivar0,
                                                const VectorInt& nbgh,
                                                const ECalcMember& member) const
{
  if (_driftList == nullptr) return MatrixRectangular();
  return _driftList->evalDriftMat(db, ivar0, nbgh, member);
}

MatrixRectangular ModelGeneric::evalDriftMatByRanks(const Db* db,
                                                    const VectorVectorInt& sampleRanks,
                                                    int ivar0,
                                                    const ECalcMember& member) const
{
  if (_driftList == nullptr) return MatrixRectangular();
  return _driftList->evalDriftMatByRanks(db, sampleRanks, ivar0, member);
}

MatrixRectangular
ModelGeneric::evalDriftMatByTarget(const Db* db,
                                    int ivar0,
                                    int iech2,
                                    const ECalcMember& member) const
{
  if (_driftList == nullptr) return MatrixRectangular();
  return _driftList->evalDriftMatByTarget(db, ivar0, iech2, member);
}






void ModelGeneric::setField(double field)
{
  _ctxt.setField(field);
  CovAnisoList* covalist = dynamic_cast<CovAnisoList*>(_cova);
  if (covalist != nullptr) covalist->copyCovContext(_ctxt);
  if (_driftList != nullptr) _driftList->copyCovContext(_ctxt);
}

bool ModelGeneric::isValid() const
{
  // Covariances: there should be some defined
  if (_cova == nullptr)
  {
    messerr("Model is not valid: no covariance has been defined");
    return false;
  }

  // Drifts: there should be valid
  if (_driftList != nullptr)
  {
    if (!_driftList->isValid()) return false;
  }

  // Check the consistency between the Covariance and the Drift parts
  int irf_drift = getDriftMaxIRFOrder();
  int irf_cova  = getCovaMinIRFOrder();
  if (irf_cova > irf_drift)
  {
    messerr("Model if invalid due to IRF degree inconsistency");
    messerr("- Covariance implies a order >= %d", irf_cova);
    messerr("- Drift implies a order %d", irf_drift);
    messerr("(Order -1 stands for strict stationarity)");
    return false;
  }
  return true;
}

// Pipes methods to _ctxt
int ModelGeneric::getNVar() const
{
  return _ctxt.getNVar();
}
int ModelGeneric::getNDim() const
{
  return _ctxt.getNDim();
}
const VectorDouble& ModelGeneric::getMeans() const
{
  return _ctxt.getMeans();
}

// Pipes method to _driftList
int ModelGeneric::getNDrift() const
{
  if (_driftList == nullptr) return 0;
  return _driftList->getNDrift();
}
int ModelGeneric::getNDriftEquation() const
{
  if (_driftList == nullptr) return 0;
  return _driftList->getNDriftEquation();
}
int ModelGeneric::getNExtDrift() const
{
  if (_driftList == nullptr) return 0;
  return _driftList->getNExtDrift();
}
int ModelGeneric::getDriftMaxIRFOrder(void) const
{
  if (_driftList == nullptr) return -1;
  return _driftList->getDriftMaxIRFOrder();
}
void ModelGeneric::delAllDrifts()
{
  if (_driftList == nullptr) return;
  _driftList->delAllDrifts();
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
