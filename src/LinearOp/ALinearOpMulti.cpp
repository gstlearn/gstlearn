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
#include "LinearOp/ALinearOpMulti.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/AException.hpp"
#include "Basic/Timer.hpp"
#include "Basic/OptDbg.hpp"
#include "Basic/Utilities.hpp"
#include "Basic/VectorHelper.hpp"
#include "Matrix/VectorEigen.hpp"
#include <Eigen/src/Core/Matrix.h>
#include <vector>

ALinearOpMulti::ALinearOpMulti(int nitermax, double eps)
    : _nIterMax(nitermax),
      _nIterRestart(0),
      _eps(eps),
      _precondStatus(false),
      _userInitialValue(false),
      _precond(nullptr),
      _initialized(false),
      _r(std::vector<std::vector<double>>()),
      _temp(),
      _p(),
      _z(),
      _nb(getNA<double>()),
      _logStats()
{
}

ALinearOpMulti::ALinearOpMulti(const ALinearOpMulti &m)
    : _nIterMax(m._nIterMax),
      _nIterRestart(0),
      _eps(m._eps),
      _precondStatus(m._precondStatus),
      _userInitialValue(m._userInitialValue),
      _precond(m._precond),
      _initialized(m._initialized),
      _r(m._r),
      _temp(m._temp),
      _p(m._p),
      _z(m._z),
      _nb(m._nb),
      _logStats(m._logStats)
{
}

ALinearOpMulti& ALinearOpMulti::operator=(const ALinearOpMulti &m)
{
  if (this != &m)
  {
    _nIterMax = m._nIterMax;
    _nIterRestart = m._nIterRestart;
    _eps = m._eps;
    _precondStatus = m._precondStatus;
    _userInitialValue = m._userInitialValue;
    _precond = m._precond;
    _initialized = m._initialized;
    _r = m._r;
    _temp = m._temp;
    _p = m._p;
    _z = m._z;
    _nb = m._nb;
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
void ALinearOpMulti::evalDirect(const std::vector<std::vector<double>> &inv,
                                std::vector<std::vector<double>> &outv) const
{
  try
  {
    _evalDirect(inv, outv);
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
void ALinearOpMulti::evalInverse(const std::vector<std::vector<double>> &vecin,
                                 std::vector<std::vector<double>> &vecout) const
{
  prepare();
  int n = sizes();
  if (n <= 0) my_throw("ALinearOpMulti size not defined. Call setSize before");

  double rsnew;
  double rsold = 0.;
  double nb;
  double crit = 0., alpha;

  Timer time;
  nb = 0.;
  for (const auto &e:vecin)
  {
    nb += VH::norm(e);
  }

  if (_userInitialValue)
  {
    evalDirect(vecout, _temp); //temp = Ax0 (x0 est stocké dans outv)
    VectorHelper::substractInPlace(_temp, vecin, _r);    //r=b-Ax0
    nb = VectorHelper::innerProduct(_r, _r);

    // If _nb is not set, then initialize the internal state from scratch.
    // If _nb is set, reuse the internal state of the solver (_p) to add
    // iterations. _nb is only needed for crit i.e. the stopping criterion.
    if (!isNA(_nb))
    {
      crit = rsold = nb;
      nb           = _nb;
    }
  }
  else
  {
    for (auto &e:vecout)
      std::fill(e.begin(),e.end(),0.);
    for (auto &e :_temp)
      std::fill(e.begin(),e.end(),0.); // temp = Ax0=0
  
    VectorHelper::copy(vecin, _r);   // r = b
  }

  if (OptDbg::query(EDbg::CONVERGE))
    message("initial crit %lg \n", VectorHelper::innerProduct(_r, _r));

  if (_precondStatus)
  {
    _precond->evalDirect(_r, _temp); //z=Mr
    VectorHelper::copy(_temp, _p); //p=z
    rsold = VectorHelper::innerProduct(_r, _temp); //<r, z>
    crit = VectorHelper::innerProduct(_r, _r);  //<r,r>
  }
  else if (!_userInitialValue || isNA(_nb)) // _p, rsold and crit are already set (see above)
  {
    VectorHelper::copy(_r, _p); //p=r (=z)
    crit = rsold = VectorHelper::innerProduct(_r, _r);
  }

  crit /= nb;

  int niter = 0;

  while (niter < _nIterMax && crit > _eps)
  {
    niter++;
    evalDirect(_p, _temp);                                // temp = Ap
    alpha = rsold / VectorHelper::innerProduct(_temp, _p);          // r'r/p'Ap
    VectorHelper::linearCombinationVVDInPlace(1., vecout, alpha, _p, vecout);     // x = x + alpha * p

    if (_nIterRestart > 0 && (niter + 1) % _nIterRestart == 0)
    {
      evalDirect(vecout, _temp);               // temp = Ax
      VectorHelper::substractInPlace(_temp, vecin, _r);   // r = b - Ax
      if (OptDbg::query(EDbg::CONVERGE))
        message("Recomputing exact residuals after %d iterations (max=%d)\n", niter, _nIterMax);
    }
    else
      VectorHelper::linearCombinationVVDInPlace(1., _r, -alpha, _temp, _r);         // r = r - alpha * Ap

    if (_precondStatus)
    {
      _precond->evalDirect(_r, _temp);                     // z = Mr
      rsnew = VectorHelper::innerProduct(_r, _temp);                 // r'z
      VectorHelper::linearCombinationVVDInPlace(1., _temp, rsnew / rsold, _p, _p); // p = z+beta p
    }
    else
    {
      rsnew = VectorHelper::innerProduct(_r, _r);
      VectorHelper::linearCombinationVVDInPlace(1., _r, rsnew / rsold, _p, _p);    // p = r+beta p
    }
    crit = rsnew / nb;

    if (OptDbg::query(EDbg::CONVERGE))
      message("%d iterations (max=%d)  crit %lg \n", niter, _nIterMax, crit);
    rsold = rsnew;
  }

  // Store _nb for further iterations (this also uses _p).
  _nb = nb;

  if (OptDbg::query(EDbg::CONVERGE))
  {
    message("-- Conjugate Gradient (precond=%d) : %d iterations (max=%d) (eps=%lg)\n",
            _precondStatus, niter, _nIterMax, _eps);
  }

  getLogStats().incrementStatsInverseCG(niter, time.getIntervalSeconds());
}

void ALinearOpMulti::initLk(const std::vector<std::vector<double>> &inv,
                            std::vector<std::vector<double>> &outv) const
{
  prepare();
  int n = sizes();
  if (n <= 0) my_throw("ALinearOpMulti size not defined. Call setSize before");

  for (auto &e: outv)
    std::fill(e.begin(),e.end(),0.);
  for (auto &e : _temp)
    std::fill(e.begin(),e.end(),0.);    // temp = Ax0=0
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
