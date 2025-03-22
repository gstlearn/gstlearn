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

CovProportional::CovProportional(const CovProportional &r)
: CovBase(r)
{
  _workMat = r._workMat;
}

CovProportional& CovProportional::operator=(const CovProportional &r)
{
  if (this != &r)
  {
    CovBase::operator=(r);
    _workMat = r._workMat;
  }
  return *this;
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


double CovProportional::_eval(const SpacePoint& p1, 
  const SpacePoint& p2,
  int ivar, 
  int jvar, 
  const CovCalcMode* mode) const
{
return _sillCur.getValue(ivar,jvar) * getCor()->evalCov(p1, p2,0, 0, mode);
}