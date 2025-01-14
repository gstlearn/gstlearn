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

MatrixRectangular ModelGeneric::evalDriftMatrix(const Db* db,
                                                int ivar0,
                                                const VectorInt& nbgh,
                                                const ECalcMember& member) const
{
  if (_driftList == nullptr) return MatrixRectangular();
  return _driftList->evalDriftMatrix(db, ivar0, nbgh, member);
}

MatrixSquareGeneral ModelGeneric::eval0Mat(const CovCalcMode* mode) const
{
  if (_cova == nullptr) return MatrixSquareGeneral();
  return _cova->eval0Mat(mode);
}

MatrixRectangular ModelGeneric::evalCovMatrix(Db* db1,
                                              Db* db2,
                                              int ivar0,
                                              int jvar0,
                                              const VectorInt& nbgh1,
                                              const VectorInt& nbgh2,
                                              const CovCalcMode* mode)
{
  if (_cova == nullptr) return MatrixRectangular();
  return _cova->evalCovMatrix(db1, db2, ivar0, jvar0, nbgh1, nbgh2, mode);
}

MatrixSquareSymmetric ModelGeneric::evalCovMatrixSymmetric(
  const Db* db1, int ivar0, const VectorInt& nbgh1, const CovCalcMode* mode)
{
  if (_cova == nullptr) return MatrixSquareSymmetric();
  return _cova->evalCovMatrixSymmetric(db1, ivar0, nbgh1, mode);
}

MatrixRectangular ModelGeneric::evalCovMatrixOptim(Db* db1,
                                                   Db* db2,
                                                   int ivar0,
                                                   int jvar0,
                                                   const VectorInt& nbgh1,
                                                   const VectorInt& nbgh2,
                                                   const CovCalcMode* mode)
{
  if (_cova == nullptr) return MatrixRectangular();
  return _cova->evalCovMatrixOptim(db1, db2, ivar0, jvar0, nbgh1, nbgh2, mode);
}

MatrixSquareSymmetric ModelGeneric::evalCovMatrixSymmetricOptim(
  const Db* db1, int ivar0, const VectorInt& nbgh1, const CovCalcMode* mode)
{
  if (_cova == nullptr) return MatrixSquareSymmetric();
  return _cova->evalCovMatrixSymmetricOptim(db1, ivar0, nbgh1, mode);
}

MatrixSparse* ModelGeneric::evalCovMatrixSparse(Db* db1,
                                                Db* db2,
                                                int ivar0,
                                                int jvar0,
                                                const VectorInt& nbgh1,
                                                const VectorInt& nbgh2,
                                                const CovCalcMode* mode,
                                                double eps)
{
  if (_cova == nullptr) return nullptr;
  return _cova->evalCovMatrixSparse(db1, db2, ivar0, jvar0, nbgh1, nbgh2, mode, eps);
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
int ModelGeneric::getVariableNumber() const
{
  return _ctxt.getNVar();
}
int ModelGeneric::getDimensionNumber() const
{
  return _ctxt.getNDim();
}
const VectorDouble& ModelGeneric::getMeans() const
{
  return _ctxt.getMean();
}

// Pipes method to _driftList
int ModelGeneric::getDriftNumber() const
{
  if (_driftList == nullptr) return 0;
  return _driftList->getDriftNumber();
}
int ModelGeneric::getDriftEquationNumber() const
{
  if (_driftList == nullptr) return 0;
  return _driftList->getDriftEquationNumber();
}
int ModelGeneric::getExternalDriftNumber() const
{
  if (_driftList == nullptr) return 0;
  return _driftList->getExternalDriftNumber();
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
double ModelGeneric::getTotalSill(int ivar, int jvar) const
{
  const CovAnisoList* covalist = getCovAnisoList();
  if (covalist == nullptr) return ITEST;
  return covalist->getTotalSill(ivar, jvar);
}
int ModelGeneric::getCovaNumber(bool skipNugget) const
{
  const CovAnisoList* covalist = getCovAnisoList();
  if (covalist == nullptr) return ITEST;
  return covalist->getCovaNumber(skipNugget);
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
