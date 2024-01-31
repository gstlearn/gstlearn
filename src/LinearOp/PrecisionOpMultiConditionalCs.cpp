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

#include "csparse_f.h"

#include <functional>

#include <math.h>

PrecisionOpMultiConditionalCs::PrecisionOpMultiConditionalCs()
    : _qChol()
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

  return _qChol.computeLogDet();
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
    pmat1->getSize();

    for (int is = 1; is < number; is++)
    {
      const PrecisionOpCs* pmataux = dynamic_cast<const PrecisionOpCs*>(getMultiPrecisionOp(is));
      Qmult = matCS_glue(Qref, pmataux->getQ(), true, true);
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
    const ProjMatrix* projElem = dynamic_cast<const ProjMatrix*>(getProjMatrix(0));
    if (projElem != nullptr) Pmult = new ProjMatrix(*projElem);
  }

  else
  {
    MatrixSparse* Amult = nullptr;
    const ProjMatrix* Pref = dynamic_cast<const ProjMatrix*>(getProjMatrix(0));
    int npoint = Pref->getPointNumber();
    int napices = Pref->getApexNumber();
    const MatrixSparse* Aref = Pref->getAproj();

    for (int is = 1; is < number; is++)
    {
      const ProjMatrix* Paux = dynamic_cast<const ProjMatrix*>(getProjMatrix(is));
      const MatrixSparse* Aaux = Paux->getAproj();
      Amult = matCS_glue(Aref, Aaux, false, true);
      napices += Paux->getApexNumber();
      Aref = Amult;
    }
    Pmult = new ProjMatrix(npoint, napices, Amult);
  }
  return Pmult;
}

int PrecisionOpMultiConditionalCs::_buildQpAtA()
{
  if (_qChol.isCholeskyDecomposed()) return 0;

  // Build the multiple projection matrix 'Amult'
  ProjMatrix* Amult = _buildAmult();
  if (Amult == nullptr) return 1;

  // Build the multiple precision matrix 'Qmult'
  MatrixSparse* Qmult = _buildQmult();
  if (Qmult == nullptr) return 1;

  // Create the conditional multiple precision matrix 'Q'
  VectorDouble invsigma = VectorHelper::inverse(getAllVarianceData());
  MatrixSparse* AtAsVar = matCS_prod_norm_diagonal(1, Amult->getAproj(), invsigma);
  MatrixSparse* Q = matCS_add(Qmult, AtAsVar, 1., 1.);

  // Prepare the Cholesky decomposition
  _qChol.reset(Q, true);
  _qChol.mustShowStats(getLogStats().isMustPrint());

  return 0;
}

void PrecisionOpMultiConditionalCs::evalInverse(const VectorVectorDouble &vecin,
                                                VectorVectorDouble &vecout) const
{
  VectorDouble locVecin = VectorHelper::flatten(vecin);
  VectorDouble locVecout(locVecin.size());
  _qChol.evalInverse(locVecin, locVecout);
  VectorHelper::unflattenInPlace(locVecout, vecout);
}

void PrecisionOpMultiConditionalCs::makeReady()
{
  // Perform Cholesky decomposition (if not already performed)
  if (_buildQpAtA()) return;
}
