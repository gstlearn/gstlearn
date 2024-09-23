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
#include "LinearOp/ALinearOpOld.hpp"
#include "Basic/AException.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/OptDbg.hpp"
#include "Basic/Timer.hpp"

ALinearOpOld::ALinearOpOld(const CGParam& params)
    : _params(params),
      _logStats()
{
}

ALinearOpOld::ALinearOpOld(const ALinearOpOld &m)
    : _params(m._params),
      _logStats(m._logStats)
{
}

ALinearOpOld& ALinearOpOld::operator=(const ALinearOpOld &m)
{
  if (this != &m)
  {
    _params = m._params;
    _logStats = m._logStats;
  }
  return *this;
}

ALinearOpOld::~ALinearOpOld() 
{
  _logStats.statsShow();
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
void ALinearOpOld::evalDirect(const VectorDouble& inv, VectorDouble& outv) const
{
  try
  {
    _evalDirect(inv,outv);
  }
  catch(const std::string& str)
  {
    // TODO : Check if std::exception can be used
    messerr("%s", str.c_str());
  }
}

/*****************************************************************************/
/*!
**  Evaluate the product: 'outv' = Q^{-1} * 'inv' by Conjugate Gradient
**
** \param[in]  inv     Array of input values
**
** \param[out] outv    Array of output values
**
*****************************************************************************/
void ALinearOpOld::evalInverse(const VectorDouble &inv, VectorDouble &outv) const
{
  int n = getSize();
  if (n <= 0) my_throw("ALinearOpOld size not defined. Call setSize before");

  Timer time;

  VectorDouble z    = VectorDouble(n);
  VectorDouble r    = VectorDouble(n);
  VectorDouble temp = VectorDouble(n);
  VectorDouble p    = VectorDouble(n);

  if (! _params.getX0().empty())
    for (int i=0; i<n; i++) outv[i] = _params.getX0(i);
  else
    for (int i=0; i<n; i++) outv[i] = 0.;

  evalDirect(outv,temp);
  for(int i=0; i<n; i++) r[i] = inv[i] - temp[i];
  
  if (_params.getPrecondStatus() == 0)
    for (int i=0; i<n; i++) z[i] = r[i];
  else if (_params.getPrecondStatus() == -1)
    _params.getPrecond()->evalInverse(r,z);
  else
    _params.getPrecond()->evalDirect(r,z);
  double critold = _prod(r, z);

  for (int i=0; i<n; i++) p[i] = z[i];

  int niter = 0;
  while(niter < _params.getNIterMax() && _prod(r,r) > _params.getEps())
  {
    niter++;
    evalDirect(p,temp);
    double alpha = _prod(r,z) / _prod(temp,p);

    for(int i=0; i<n; i++)
    {
      outv[i] += alpha * p[i];
      r[i]   -= alpha * temp[i];
    }
    if (_params.getPrecondStatus() == 0)
      for (int i=0; i<n; i++) z[i] = r[i];
    else if (_params.getPrecondStatus() == -1)
      _params.getPrecond()->evalInverse(r,z);
    else 
      _params.getPrecond()->evalDirect(r,z);

    double critnew = _prod(r, z);
    double beta    = critnew / critold;

    for(int i=0; i<n; i++)
      p[i] = z[i] + beta * p[i];
    critold = critnew;
  }

  if (OptDbg::query(EDbg::CONVERGE))
  {
    message("-- Conjugate Gradient (precond=%d) : %d iterations (max=%d) (eps=%lg)\n",
            _params.getPrecondStatus(),niter,_params.getNIterMax(),_params.getEps());
  }

  getLogStats().incrementStatsInverseCG(niter, time.getIntervalSeconds());
}

/*****************************************************************************/
/*!
**  Returns the scalar product between 'x' and 'y'
**
** \param[in]  x      First array
** \param[in]  y      Second array
**
*****************************************************************************/
double ALinearOpOld::_prod(const VectorDouble &x, const VectorDouble &y) const
{
  double prod = 0.;
  for (int i=0, n = getSize(); i<n; i++)
    prod += x[i] * y[i];
  return prod;
}
