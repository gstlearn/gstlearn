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
#include "geoslib_old_f.h"

#include "LinearOp/CGParam.hpp"

CGParam::CGParam(int nitermax, double eps)
  : _nIterMax(nitermax)
  , _eps(eps)
  , _x0()
  , _precondStatus(0)
  , _precond(nullptr)
{
}

CGParam::CGParam(const CGParam &m)
    : _nIterMax(m._nIterMax),
      _eps(m._eps),
      _x0(m._x0),
      _precondStatus(m._precondStatus),
      _precond(m._precond)
{
}

CGParam& CGParam::operator=(const CGParam &m)
{
  if (this != &m)
  {
    _nIterMax = m._nIterMax;
    _eps = m._eps;
    _x0 = m._x0;
    _precondStatus = m._precondStatus;
    _precond = m._precond;
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
void CGParam::setPrecond(const ALinearOp* precond, int status)
{ 
  _precond = precond; 
  _precondStatus = status;
  if (precond == NULL) _precondStatus = 0;
}
