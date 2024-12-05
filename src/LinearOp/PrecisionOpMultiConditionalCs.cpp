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
#include "LinearOp/PrecisionOpMultiConditionalCs.hpp"
#include "LinearOp/PrecisionOpMatrix.hpp"

#include "Basic/VectorHelper.hpp"
#include "Matrix/MatrixSparse.hpp"
#include "Matrix/MatrixFactory.hpp"

#include <math.h>
#include <vector>

PrecisionOpMultiConditionalCs::PrecisionOpMultiConditionalCs()
    : _Q(nullptr)
    , _chol(nullptr)
{

}

PrecisionOpMultiConditionalCs::~PrecisionOpMultiConditionalCs()
{
  _clear();
}

void PrecisionOpMultiConditionalCs::_clear()
{
  delete _chol;
  _chol = nullptr;
  delete _Q;
  _Q = nullptr;
}

int PrecisionOpMultiConditionalCs::push_back(PrecisionOp* pmatElem, IProj* projDataElem)
{
  _clear();
  PrecisionOpMatrix* pmatElemCs = dynamic_cast<PrecisionOpMatrix*>(pmatElem);
  if (pmatElemCs == nullptr)
  {
    messerr("The first argument of 'push_back' should be a pointer to PrecisionOpMatrix");
    return 1;
  }
  return PrecisionOpMultiConditional::push_back(pmatElem, projDataElem);
}

double PrecisionOpMultiConditionalCs::computeLogDetOp(int nbsimu) const
{
  DECLARE_UNUSED(nbsimu);

  if (_chol == nullptr)
    _chol = new CholeskySparse(_Q);
  return _chol->computeLogDeterminant();
}

MatrixSparse* PrecisionOpMultiConditionalCs::_buildQmult() const
{
  MatrixSparse* Qmult = nullptr;
  int number = sizes();
  if (number <= 0)
  {
    messerr("This method requires at least one registered covariance");
    return Qmult;
  }

  // Particular case of a single registered covariance
  if (number == 1)
  {
    const PrecisionOpMatrix* pmatElem = dynamic_cast<const PrecisionOpMatrix*>(getMultiPrecisionOp(0));
    if (pmatElem != nullptr) Qmult = pmatElem->getQ()->clone();
  }
  else
  {
    const PrecisionOpMatrix* pmat1 = dynamic_cast<const PrecisionOpMatrix*>(getMultiPrecisionOp(0));
    const MatrixSparse* Qref = pmat1->getQ();

    for (int is = 1; is < number; is++)
    {
      const PrecisionOpMatrix* pmataux = dynamic_cast<const PrecisionOpMatrix*>(getMultiPrecisionOp(is));
      delete Qmult;
      Qmult = dynamic_cast<MatrixSparse*>(MatrixFactory::createGlue(Qref, pmataux->getQ(), true, true));
      Qref = Qmult;
    }
  }
  return Qmult;
}

ProjMatrix* PrecisionOpMultiConditionalCs::_buildAmult() const
{
  ProjMatrix* Pmult = nullptr;
  int number = sizes();
  if (number <= 0)
  {
    messerr("This method requires at least one registered projection matrix");
    return Pmult;
  }

  // Particular case of a single registered covariance
  if (number == 1)
  {
    const ProjMatrix* projElem = getProjMatrix(0);
    if (projElem != nullptr) Pmult = new ProjMatrix(*projElem);
  }
  else
  {
    MatrixSparse* mstemp = nullptr;
    const MatrixSparse* msref = dynamic_cast<const MatrixSparse*>(getProjMatrix(0));
    for (int is = 1; is < number; is++)
    {
      const MatrixSparse* msaux = dynamic_cast<const MatrixSparse*>(getProjMatrix(is));
      delete mstemp;
      mstemp = dynamic_cast<MatrixSparse*>(MatrixFactory::createGlue(msref, msaux, false, true));
      msref = mstemp;
    }
    Pmult = new ProjMatrix(mstemp);
    delete mstemp;
  }
  return Pmult;
}

int PrecisionOpMultiConditionalCs::_buildQpAtA()
{
  if (_Q != nullptr) return 0;

  // Build the multiple projection matrix 'Amult'
  ProjMatrix* Amult = _buildAmult();
  if (Amult == nullptr) return 1;

  // Build the multiple precision matrix 'Qmult'
  MatrixSparse* Qmult = _buildQmult();
  if (Qmult == nullptr) return 1;

  // Create the conditional multiple precision matrix 'Q'
  VectorDouble invsigma = VectorHelper::inverse(getAllVarianceData());
  MatrixSparse* AtAsVar = prodNormMat(Amult, invsigma, true);
  _Q = MatrixSparse::addMatMat(Qmult, AtAsVar, 1., 1.);

  // Free core allocated
  delete Amult;
  delete Qmult;
  delete AtAsVar;

  return 0;
}

void PrecisionOpMultiConditionalCs::evalInverse(const std::vector<std::vector<double>> &vecin,
                                                std::vector<std::vector<double>> &vecout) const
{
  if (_chol == nullptr)
    _chol = new CholeskySparse(_Q);
  std::vector<double> locVecin = VH::flatten(vecin);
  std::vector<double> locVecout(locVecin.size());
  _chol->solve(locVecin, locVecout);
  VH::unflattenInPlace(locVecout, vecout);
}

void PrecisionOpMultiConditionalCs::makeReady()
{
  // Perform Cholesky decomposition (if not already performed)
  _buildQpAtA();
}
