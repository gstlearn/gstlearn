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
#include "Covariances/ACor.hpp"
#include "Covariances/ACov.hpp"
#include "Covariances/CovBase.hpp"
#include "Matrix/MatrixSquareSymmetric.hpp"

CovProportional::CovProportional(ACor* cor,
                const MatrixSquareSymmetric &sill)
: CovBase(cor,sill)
{
  _workMat.resize(_ctxt.getNVar(), _ctxt.getNVar());
  _workMat.setIdentity();
}

CovProportional::~CovProportional()
{

}

void CovProportional::setCor(ACor* cor)
{
  if (cor->getNVariables() != 1)
  {
    messerr("Correlation function should have only 1 variable");
    return;
  }
  setCor(cor);
}
