/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c3) MINES Paris / ARMINES                                        */
/* Authors: gstlearn Team                                                     */
/* Website: https://gstlearn.org                                              */
/* License: BSD 3-clause                                                      */
/*                                                                            */
/******************************************************************************/

#include "Covariances/CovProportional.hpp"
#include "Basic/AStringable.hpp"
#include "Covariances/ACov.hpp"
#include "Covariances/CovBase.hpp"
#include "Matrix/MatrixSquareSymmetric.hpp"

CovProportional::CovProportional(ACov* cor,
                const MatrixSquareSymmetric &sill)
: CovBase(cor,sill)
{
  _ctxt.setNVar(sill.getNCols());
  _workMat.resize(_ctxt.getNVar(), _ctxt.getNVar());
  _workMat.setIdentity();
  if (cor != nullptr)
    if (cor->getNVar() != 1)
    {
        messerr("Correlation function should have only 1 variable");
        messerr("You should use CovBase instead of CovProportional");
        messerr("Undefined behaviour");
        return;
    }
}

CovProportional::~CovProportional()
{

}

void CovProportional::setCor(ACov* cor)
{
  if (cor->getNVar() != 1)
  {
    messerr("Correlation function should have only 1 variable");
    return;
  }
  CovBase::setCor(cor);
}

/**
 * Calculate the Matrix of covariance between two space points
 * @param p1 Reference of the first space point
 * @param p2 Reference of the second space point
 * @param mat   Covariance matrix (Dimension: nvar * nvar)
 * @param mode  Calculation Options
 *
 * @remarks: Matrix 'mat' should be dimensioned and initialized beforehand
 */
void CovProportional::_addEvalCovMatBiPointInPlace(MatrixSquareGeneral &mat,
                                          const SpacePoint &p1,
                                          const SpacePoint &p2,
                                          const CovCalcMode *mode) const
{
  double cor = getCor()->eval(p1, p2, 0, 0, mode);

  if (mode == nullptr || ! mode->getUnitary())
    mat.addMatInPlace(_sill, 1., cor);
  else
  {
    mat.addMatInPlace(_workMat, 1., cor);
  }
}