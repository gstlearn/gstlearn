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
#include "LinearOp/PrecisionOpCs.hpp"

#include "Basic/Law.hpp"
#include "Basic/VectorHelper.hpp"
#include "Matrix/MatrixSquareSymmetric.hpp"
#include "Matrix/MatrixSparse.hpp"
#include "Polynomials/Chebychev.hpp"
#include "LinearOp/Cholesky.hpp"

#include <functional>

#include <math.h>

PrecisionOpMultiConditionalCs::PrecisionOpMultiConditionalCs()
    : _Q(nullptr)
{

}

PrecisionOpMultiConditionalCs::~PrecisionOpMultiConditionalCs()
{
}

int PrecisionOpMultiConditionalCs::push_back(PrecisionOp* pmatElem, IProjMatrix* projDataElem)
{
  PrecisionOpCs* pmatElemCs = dynamic_cast<PrecisionOpCs*>(pmatElem);
  if (pmatElemCs == nullptr)
  {
    messerr("The first argument of 'push_back' should be a pointer to PrecisionOpCs");
    return 1;
  }
  return PrecisionOpMultiConditional::push_back(pmatElem, projDataElem);
}

double PrecisionOpMultiConditionalCs::computeLogDetOp(int nbsimu, int seed) const
{
  DECLARE_UNUSED(nbsimu);
  DECLARE_UNUSED(seed);

  return _Q->getCholeskyLogDeterminant();
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
    const PrecisionOpCs* pmatElem = dynamic_cast<const PrecisionOpCs*>(getMultiPrecisionOp(0));
    if (pmatElem != nullptr) Qmult = pmatElem->getQ()->clone();
  }
  else
  {
    const PrecisionOpCs* pmat1 = dynamic_cast<const PrecisionOpCs*>(getMultiPrecisionOp(0));
    const MatrixSparse* Qref = pmat1->getQ();

    for (int is = 1; is < number; is++)
    {
      const PrecisionOpCs* pmataux = dynamic_cast<const PrecisionOpCs*>(getMultiPrecisionOp(is));
      if (Qmult != nullptr) delete Qmult;
      Qmult = MatrixSparse::glue(Qref, pmataux->getQ(), true, true);
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
      if (mstemp != nullptr) delete mstemp;
      mstemp = MatrixSparse::glue(msref, msaux, false, true);
      msref = mstemp;
    }
    Pmult = new ProjMatrix(mstemp);
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
  MatrixSparse* AtAsVar = prodNormMat(*Amult, invsigma, true);
  _Q = MatrixSparse::addMatMat(Qmult, AtAsVar, 1., 1.);

  // Prepare the Cholesky decomposition
  _Q->computeCholesky();

  return 0;
}

void PrecisionOpMultiConditionalCs::evalInverse(const VectorVectorDouble &vecin,
                                                VectorVectorDouble &vecout) const
{
  VectorDouble locVecin = VectorHelper::flatten(vecin);
  VectorDouble locVecout(locVecin.size());
  _Q->solveCholesky(locVecin, locVecout);
  VectorHelper::unflattenInPlace(locVecout, vecout);
}

void PrecisionOpMultiConditionalCs::makeReady()
{
  // Perform Cholesky decomposition (if not already performed)
  _buildQpAtA();
}
