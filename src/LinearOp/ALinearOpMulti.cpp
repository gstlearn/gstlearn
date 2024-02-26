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

#include "LinearOp/ALinearOpMulti.hpp"
#include "LinearOp/Identity.hpp"
#include "Basic/AException.hpp"
#include "Basic/VectorHelper.hpp"
#include "Basic/Timer.hpp"
#include "Basic/OptDbg.hpp"

#include <iostream>

ALinearOpMulti::ALinearOpMulti(int nitermax, double eps)
    : _nIterMax(nitermax),
      _eps(eps),
      _precondStatus(false),
      _userInitialValue(false),
      _precond(nullptr),
      _initialized(false),
      _r(VectorVectorDouble()),
      _temp(),
      _p(),
      _z(),
      _logStats()
{
}

ALinearOpMulti::ALinearOpMulti(const ALinearOpMulti &m)
    : _nIterMax(m._nIterMax),
      _eps(m._eps),
      _precondStatus(m._precondStatus),
      _userInitialValue(m._userInitialValue),
      _precond(m._precond),
      _initialized(m._initialized),
      _r(m._r),
      _temp(m._temp),
      _p(m._p),
      _z(m._z),
      _logStats(m._logStats)
{
}

ALinearOpMulti& ALinearOpMulti::operator=(const ALinearOpMulti &m)
{
  if (this != &m)
  {
    _nIterMax = m._nIterMax;
    _eps = m._eps;
    _precondStatus = m._precondStatus;
    _userInitialValue = m._userInitialValue;
    _precond = m._precond;
    _initialized = m._initialized;
    _r = m._r;
    _temp = m._temp;
    _p = m._p;
    _z = m._z;
    _logStats = m ._logStats;
  }
  return *this;
}

ALinearOpMulti::~ALinearOpMulti()
{
  _logStats.statsShow();
}

/**
 * This method intends to resize the different working arrays.
 * It is considered as const in order to avoid breaking constness of calling function
 */
void ALinearOpMulti::prepare() const
{
  if (_initialized) return;

  _initialized = true;
  int ns = sizes();
  _r.resize(ns);
  _p.resize(ns);
  _temp.resize(ns);
  _z.resize(ns);

  for (int is = 0; is < ns; is++)
  {
    int n = size(is);
    _r[is].resize(n);
    _p[is].resize(n);
    _temp[is].resize(n);
    _z[is].resize(n);
  }
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
void ALinearOpMulti::evalDirect(const VectorVectorDouble &inv,
                                VectorVectorDouble &outv) const
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
 **  Evaluate the product: 'outv' = Q^{-1} * 'inv' by conjugate gradient
 **
 ** \param[in]  vecin   Array of input values
 **
 ** \param[out] vecout  Array of output values. Will be used as initial value if
 **                    _userInitialValue is true.
 **
 *****************************************************************************/
void ALinearOpMulti::evalInverse(const VectorVectorDouble &vecin,
                                 VectorVectorDouble &vecout) const
{
  prepare();
  int n = sizes();
  if (n <= 0) my_throw("ALinearOpMulti size not defined. Call setSize before");

  double rsnew;
  double rsold;
  double nb;
  double crit, alpha;

  Timer time;
  nb = VH::innerProduct(vecin, vecin);

  if (_userInitialValue)
  {
    evalDirect(vecout, _temp); //temp = Ax0 (x0 est stocké dans outv)
    VH::subtractInPlace(_temp, vecin, _r);    //r=b-Ax0
  }
  else
  {
    VH::fill(vecout, 0.);
    VH::fill(_temp, 0.); // temp = Ax0=0
    VH::copy(vecin, _r);   // r = b
  }

  if (OptDbg::query(EDbg::CONVERGE))
    message("initial crit %lg \n", VH::innerProduct(_r, _r));

  if (_precondStatus)
  {
    _precond->evalDirect(_r, _temp); //z=Mr
    VH::copy(_temp, _p); //p=z
    rsold = VH::innerProduct(_r, _temp); //<r, z>
    crit = VH::innerProduct(_r, _r);  //<r,r>
  }
  else
  {
    VH::copy(_r, _p); //p=r (=z)
    crit = rsold = VH::innerProduct(_r, _r);
  }

  crit /= nb;

  int niter = 0;

  while (niter < _nIterMax && crit > _eps)
  {
    niter++;
    evalDirect(_p, _temp);                                // temp = Ap
    alpha = rsold / VH::innerProduct(_temp, _p);          // r'r/p'Ap
    VH::linearCombVVD(1., vecout, alpha, _p, vecout);     // x = x + alpha * p
    VH::linearCombVVD(1., _r, -alpha, _temp, _r);         // r = r - alpha * Ap

    if (_precondStatus)
    {
      _precond->evalDirect(_r, _temp);                     // z = Mr
      rsnew = VH::innerProduct(_r, _temp);                 // r'z
      VH::linearCombVVD(1., _temp, rsnew / rsold, _p, _p); // p = z+beta p
    }
    else
    {
      rsnew = VH::innerProduct(_r, _r);
      crit = rsnew / nb;
      VH::linearCombVVD(1., _r, rsnew / rsold, _p, _p);    // p = r+beta p
    }

    if (OptDbg::query(EDbg::CONVERGE))
      message("%d iterations (max=%d)  crit %lg \n", niter, _nIterMax, crit);
    rsold = rsnew;
  }

  if (OptDbg::query(EDbg::CONVERGE))
  {
    message("-- Conjugate Gradient (precond=%d) : %d iterations (max=%d) (eps=%lg)\n",
            _precondStatus, niter, _nIterMax, _eps);
  }

  getLogStats().incrementStatsInverseCG(niter, time.getIntervalSeconds());
}

void ALinearOpMulti::initLk(const VectorVectorDouble &inv,
                            VectorVectorDouble &outv) const
{
  prepare();
  int n = sizes();
  if (n <= 0) my_throw("ALinearOpMulti size not defined. Call setSize before");

  VH::fill(outv,0.);
  VH::fill(_temp,0.);    // temp = Ax0=0
  VH::copy(inv,_p);     // p = r (=z)
  evalDirect(_p,_temp); // temp = Ap
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
void ALinearOpMulti::setPrecond(const ALinearOpMulti* precond, int status)
{ 
  _precond = precond; 
  _precondStatus = status;
  if (_precond == nullptr) _precondStatus = false;
}

void ALinearOpMulti::_updated()const
{
  _initialized = false;
}
