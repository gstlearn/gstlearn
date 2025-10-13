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
#include "LinearOp/CGParam.hpp"

namespace gstlrn
{
CGParam::CGParam(Id nitermax, double eps)
  : AStringable()
  , _nIterMax(nitermax)
  , _eps(eps)
  , _x0()
  , _precondStatus(0)
  , _precond(nullptr)
{
}

CGParam::CGParam(const CGParam& m)
  : AStringable(m)
  , _nIterMax(m._nIterMax)
  , _eps(m._eps)
  , _x0(m._x0)
  , _precondStatus(m._precondStatus)
  , _precond(m._precond)
{
}

CGParam& CGParam::operator=(const CGParam& m)
{
  if (this != &m)
  {
    AStringable::operator=(m);
    _nIterMax      = m._nIterMax;
    _eps           = m._eps;
    _x0            = m._x0;
    _precondStatus = m._precondStatus;
    _precond       = m._precond;
  }
  return *this;
}

CGParam::~CGParam()
{
}

/*****************************************************************************/
/*!
**  Define the Pre-Conditioner facility
**
** \param[in]  precond  Pointer to a ALinearOp operator
** \param[in]  status   Status of this Pre-conditioner
** \li                  0 : not defined and therefore not used
** \li                 -1 : Pre-conditioner is the Q_{-1}
** \li                  1 : Pre-conditioner is the Q
**
** \remarks When 'precond' argument is not provided, 'status' is forced to 0
**
*****************************************************************************/
void CGParam::setPrecond(const ALinearOp* precond, Id status)
{
  _precond       = precond;
  _precondStatus = status;
  if (precond == nullptr) _precondStatus = 0;
}

String CGParam::toString(const AStringFormat* strfmt) const
{
  DECLARE_UNUSED(strfmt);
  std::stringstream sstr;

  sstr << "Maximum number of Conjugate Gradient iterations = " << _nIterMax << std::endl;
  sstr << "Numerical tolerance = " << _eps << std::endl;
  sstr << "Initial value = " << _x0 << std::endl;
  sstr << "Using a Pre-conditioner = " << _precondStatus << std::endl;
  return sstr.str();
}
} // namespace gstlrn