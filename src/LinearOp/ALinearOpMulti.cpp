/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
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
    : _temp(VectorVectorDouble()),
      _p(VectorVectorDouble()),
      _nIterMax(nitermax),
      _eps(eps),
      _precondStatus(false),
      _userInitialValue(false),
      _precond(nullptr),
      _initialized(false),
      _r(VectorVectorDouble()),
      _timeCG(0),
      _niterCG(0),
      _numberCG(0)
{
}

ALinearOpMulti::ALinearOpMulti(const ALinearOpMulti &m)
    : _temp(m._temp),
      _p(m._p),
      _nIterMax(m._nIterMax),
      _eps(m._eps),
      _precondStatus(m._precondStatus),
      _userInitialValue(m._userInitialValue),
      _precond(m._precond),
      _initialized(m._initialized),
      _r(m._r),
      _timeCG(m._timeCG),
      _niterCG(m._niterCG),
      _numberCG(m._numberCG)
{
}

ALinearOpMulti& ALinearOpMulti::operator=(const ALinearOpMulti &m)
{
  if (this != &m)
  {
    _temp = m._temp;
    _p = m._p;
    _nIterMax = m._nIterMax;
    _eps = m._eps;
    _precondStatus = m._precondStatus;
    _userInitialValue = m._userInitialValue;
    _precond = m._precond;
    _initialized = m._initialized;
    _r = m._r;
    _timeCG = m._timeCG;
    _niterCG = m._niterCG;
    _numberCG = m._numberCG;
  }
  return *this;

}

ALinearOpMulti::~ALinearOpMulti()
{
}

/**
 * This method intends to size the different working arrays.
 */

void ALinearOpMulti::_init() const
{
  if(_initialized)
  {
    return;
  }
  _initialized = true;
  int n;
  int ns = sizes();
  _r.resize(ns);
  _p.resize(ns);
  _temp.resize(ns);
  _z.resize(ns);

  for(int is=0;is<ns;is++)
  {
    n = size(is);
    _r[is].resize(n);
    _p[is].resize(n);
    _temp[is].resize(n);
    _z[is].resize(n);
  }
}

void ALinearOpMulti::_initPublic() const{
  _init();
};

/*****************************************************************************/
/*!
 **  Evaluate the product: 'outv' = Q * 'inv'
 **
 ** \param[in]  inv     Array of input values
 **
 ** \param[out] outv    Array of output values
 **
 *****************************************************************************/
void ALinearOpMulti::evalDirect(const VectorVectorDouble& inv,
                                VectorVectorDouble& outv) const
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
 **  Evaluate the product: 'outv' = Q^{-1} * 'inv' by conjugate gradient
 **
 ** \param[in]  inv     Array of input values
 **
 ** \param[out] outv    Array of output values. Will be used as initial value if
 **                    _userInitialValue is true.
 **
 *****************************************************************************/
void ALinearOpMulti::evalInverse(const VectorVectorDouble &inv,
                                 VectorVectorDouble &outv) const
{
  _init();
  int n = sizes();
  if (n <= 0) my_throw("ALinearOpMulti size not defined. Call setSize before");

  double rsnew;
  double rsold;
  double nb;
  double crit,alpha;

  Timer time;

  nb = innerProduct(inv,inv);


  if(_userInitialValue)
  {

    evalDirect(outv,_temp); //temp = Ax0 (x0 est stocké dans outv)
    diff(_temp,inv,_r);    //r=b-Ax0
  }
  else
  {
    fillVal(outv,0.);
    fillVal(_temp,0.); // temp = Ax0=0
    _copyVals(inv,_r);   // r = b
  }

  if (OptDbg::query(EDbg::CONVERGE))
    message("initial crit %lg \n",innerProduct(_r,_r));

  if(_precondStatus)
  {
    _precond->evalDirect(_r,_temp); //z=Mr
    _copyVals(_temp,_p); //p=z
    rsold=innerProduct(_r,_temp); //<r, z>
    crit=innerProduct(_r,_r);  //<r,r>

  }
  else
  {
    _copyVals(_r,_p); //p=r (=z)
    crit=rsold=innerProduct(_r,_r);
  }

  crit/=nb;

  int niter = 0;

  while(niter < _nIterMax && crit > _eps)
  {
    niter++;
    evalDirect(_p,_temp); //temp = Ap
    alpha = rsold / innerProduct(_temp,_p); // r'r/p'Ap
    _linearComb(1.,outv,alpha,_p,outv);//x=x+alpha p
    _linearComb(1.,_r,-alpha,_temp,_r); //r=r-alpha*Ap

    if(_precondStatus)
    {
      _precond->evalDirect(_r,_temp); //z=Mr
      rsnew=innerProduct(_r,_temp); //r'z
      _linearComb(1.,_temp,rsnew/rsold,_p,_p); //p=z+beta p
    }
    else
    {
      rsnew=innerProduct(_r,_r);
      crit=rsnew/nb;
      _linearComb(1.,_r,rsnew/rsold,_p,_p);//p=r+beta p

    }

    if (OptDbg::query(EDbg::CONVERGE))
      message("%d iterations (max=%d)  crit %lg \n",niter,_nIterMax,crit);
    rsold = rsnew;

  }

  if (OptDbg::query(EDbg::CONVERGE))
  {
    message("-- Conjugate Gradient (precond=%d) : %d iterations (max=%d) (eps=%lg)\n",
            _precondStatus,niter,_nIterMax,_eps);
  }

  _timeCG   += time.getIntervalSeconds();
  _niterCG  += niter;
  _numberCG ++;

}

void ALinearOpMulti::initLk(const VectorVectorDouble& inv,
                                 VectorVectorDouble& outv) const
{
  _init();
  int n = sizes();
  if (n <= 0) my_throw("ALinearOpMulti size not defined. Call setSize before");

  fillVal(outv,0.);
  fillVal(_temp,0.); // temp = Ax0=0
  _copyVals(inv,_p); //p=r (=z)

  evalDirect(_p,_temp); //temp = Ap

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

/*****************************************************************************/
/*!
 **  Reset the statistics
 **
 ** \remarks These statistics concern:
 ** \li The number of calls to Conjugate Gradient algorithm
 ** \li The number of iterations in Conjugate Gradient algorithm
 ** \li The time spent in Conjugate Gradient algorithm
 **
 *****************************************************************************/
void ALinearOpMulti::resetStatCG(void) const
{
  _timeCG   = 0;
  _niterCG  = 0;
  _numberCG = 0;
}

/*****************************************************************************/
/*!
 **  Trigger the printout of the statistics
 **
 ** \remarks These statistics concern:
 ** \li The number of calls to Conjugate Gradient algorithm
 ** \li The number of iterations in Conjugate Gradient algorithm
 ** \li The time spent in Conjugate Gradient algorithm
 **
 *****************************************************************************/
void ALinearOpMulti::printStatCG(void) const
{
  message("Conjugate Gradient:\n");
  message("- Number of calls      = %d \n",_numberCG);
  message("- Number of iterations = %d \n",_niterCG);
  message("- Time spent           = %lf (s)\n",_timeCG);
}

/*****************************************************************************/
/*!
 **  Returns the inner product between 'x' and 'y'
 **
 ** \param[in]  x      First array
 ** \param[in]  y      Second array
 **
 *****************************************************************************/
double ALinearOpMulti::innerProduct(const VectorDouble& x,
                             const VectorDouble& y) const
{
  double prod = 0.;
  int n = static_cast<int> (x.size());
  for(int i=0; i<n;i++)
  {
      prod += x[i] * y[i];
  }
  return prod;
}



/*****************************************************************************/
/*!
 **  Returns the inner product between 'x' and 'y'
 **
 ** \param[in]  x      First array
 ** \param[in]  y      Second array
 **
 *****************************************************************************/
double ALinearOpMulti::innerProduct(const VectorVectorDouble& x,
                             const VectorVectorDouble& y) const
{
  double s = 0.;

  for (int i = 0; i<(int)x.size();i++)
    s+= VH::innerProduct(x[i],y[i]);

  return s;
}


void ALinearOpMulti::diff(const VectorVectorDouble& in1,
                           const VectorVectorDouble& in2,
                           VectorVectorDouble& outv) const
{
  for(int is = 0;is<sizes();is++)
  {
    for(int i = 0;i<size(is);i++)
    {
      outv[is][i] = in2[is][i] - in1[is][i];
    }
  }
}


void ALinearOpMulti::sum(const VectorVectorDouble& in1,
         const VectorVectorDouble& in2,
         VectorVectorDouble& outv) const
{
  for(int is = 0;is<sizes();is++)
  {
    for(int i = 0;i<size(is);i++)
    {
      outv[is][i] = in2[is][i] + in1[is][i];
    }
  }
}
/**
 * outv = val1 * in1 + val2  in2
 */

void ALinearOpMulti::_linearComb(double val1,
                                 const VectorVectorDouble& in1,
                                 double val2,
                                 const VectorVectorDouble& in2,
                                 VectorVectorDouble& outv) const
{
  for (int is = 0; is < sizes(); is++)
  {
    for (int i = 0; i < size(is); i++)
    {
      outv[is][i] = val1 * in1[is][i] + val2 * in2[is][i];
    }
  }
}

void ALinearOpMulti::prodScalar(double val1,
                                const VectorVectorDouble& in1,
                                VectorVectorDouble& outv) const
{
  for (int is = 0; is < sizes(); is++)
  {
    for (int i = 0; i < size(is); i++)
    {
      outv[is][i] = val1 * in1[is][i];
    }
  }
}

void ALinearOpMulti::addProdScalar(double val1,
                                   const VectorVectorDouble& in1,
                                   VectorVectorDouble& outv) const
{
  for (int is = 0; is < sizes(); is++)
  {
    for (int i = 0; i < size(is); i++)
    {
      outv[is][i] += val1 * in1[is][i];
    }
  }
}

void ALinearOpMulti::fillVal(VectorVectorDouble& vect, double val) const
{
  VH::fill(vect, val);
}

double ALinearOpMulti::max(const VectorVectorDouble& vect) const
{
  return VH::maximum(vect);
}

void ALinearOpMulti::_updated()const
{
  _initialized=false;
}

void ALinearOpMulti::_copyVals(const VectorVectorDouble& inv,
                               VectorVectorDouble& outv) const
{
  for (int is = 0; is < sizes(); is++)
  {
    for (int i = 0; i < size(is); i++)
    {
      outv[is][i] = inv[is][i];
    }
  }
}
