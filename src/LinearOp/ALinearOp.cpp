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

#include "LinearOp/ALinearOp.hpp"
#include "LinearOp/Identity.hpp"
#include "Basic/AException.hpp"
#include "Basic/OptDbg.hpp"

#include <iostream>

ALinearOp::ALinearOp(int nitermax, double eps)
  : _nIterMax(nitermax)
  , _eps(eps)
  , _x0()
  , _precondStatus(0)
  , _precond(nullptr)
{
}

ALinearOp::ALinearOp(const ALinearOp &m)
    : _nIterMax(m._nIterMax),
      _eps(m._eps),
      _x0(m._x0),
      _precondStatus(m._precondStatus),
      _precond(m._precond)
{
}

ALinearOp& ALinearOp::operator=(const ALinearOp &m)
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

ALinearOp::~ALinearOp() 
{
}

/*****************************************************************************/
/*!
**  Evaluate the product: 'outv' = Q * 'inv'
**
** \param[in]  inv     Array of input values
**
** \param[out] outv    Array of output values
**
*****************************************************************************/
void ALinearOp::evalDirect(const VectorDouble& inv, VectorDouble& outv) const
{
  try
  {
    _evalDirect(inv,outv);
  }
  catch(const char * str)
  {
    // TODO : Check if std::exception can be used
    messerr("%s", str);
  }
}

/*****************************************************************************/
/*!
**  Evaluate the product: 'outv' = Q^{-1} * 'inv'
**
** \param[in]  inv     Array of input values
**
** \param[out] outv    Array of output values
**
*****************************************************************************/
void ALinearOp::evalInverse(const VectorDouble &inv, VectorDouble &outv) const
{
  int n = getSize();
  if (n <= 0) my_throw("ALinearOp size not defined. Call setSize before");

  VectorDouble z    = VectorDouble(n);
  VectorDouble r    = VectorDouble(n);
  VectorDouble temp = VectorDouble(n);
  VectorDouble p    = VectorDouble(n);

  if (! _x0.empty())
    for (int i=0; i<n; i++) outv[i] = _x0[i];
  else
    for (int i=0; i<n; i++) outv[i] = 0.;

  evalDirect(outv,temp);
  for(int i=0; i<n; i++) r[i] = inv[i] - temp[i];
  
  if (_precondStatus == 0)
    for (int i=0; i<n; i++) z[i] = r[i];
  else if (_precondStatus == -1)
    _precond->evalInverse(r,z);
  else
    _precond->evalDirect(r,z);
  double critold = _prod(r, z);

  for (int i=0; i<n; i++) p[i] = z[i];

  int niter = 0;
  while(niter < _nIterMax && _prod(r,r) > _eps)
  {
    niter++;
    evalDirect(p,temp);
    double alpha = _prod(r,z) / _prod(temp,p);

    for(int i=0; i<n; i++)
    {
      outv[i] += alpha * p[i];
      r[i]   -= alpha * temp[i];
    }
    if (_precondStatus == 0)
      for (int i=0; i<n; i++) z[i] = r[i];
    else if (_precondStatus == -1)
      _precond->evalInverse(r,z);
    else 
      _precond->evalDirect(r,z);

    double critnew = _prod(r, z);
    double beta    = critnew / critold;

    for(int i=0; i<n; i++)
      p[i] = z[i] + beta * p[i];
    critold = critnew;
  }

  if (OptDbg::query(EDbg::CONVERGE))
  {
    message("-- Conjugate Gradient (precond=%d) : %d iterations (max=%d) (eps=%lg)\n",
            _precondStatus,niter,_nIterMax,_eps);
  }
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
void ALinearOp::setPrecond(const ALinearOp* precond, int status)
{ 
  _precond = precond; 
  _precondStatus = status;
  if (precond == NULL) _precondStatus = 0;
}

/*****************************************************************************/
/*!
**  Returns the scalar product between 'x' and 'y'
**
** \param[in]  x      First array
** \param[in]  y      Second array
**
*****************************************************************************/
double ALinearOp::_prod(const VectorDouble &x, const VectorDouble &y) const
{
  double prod = 0.;
  int n = getSize();
  
  for (int i=0; i<n; i++)
    prod += x[i] * y[i];
  return prod;
}
