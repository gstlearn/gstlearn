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

#include "Basic/AException.hpp"
#include "Basic/AFunction.hpp"
#include "Polynomials/Chebychev.hpp"
#include "LinearOp/ALinearOpMulti.hpp"
#include "Matrix/LinkMatrixSparse.hpp"

#include <math.h>
#include <functional>

// External library /// TODO : Dependency to csparse to be removed
#include "csparse_d.h"
#include "csparse_f.h"

Chebychev::Chebychev()
  : _ncMax(10001)
  , _nDisc(100)
  , _a(0.)
  , _b(1.)
  , _verbose(false)
{
  // TODO Auto-generated constructor stub

}

Chebychev::~Chebychev()
{
  // TODO Auto-generated destructor stub
}

Chebychev* Chebychev::createFromCoeffs(const VectorDouble coeffs)
{
  Chebychev* cheb = new Chebychev();
  cheb->setCoeffs(coeffs);
  return cheb;
}

void Chebychev::init(int ncMax,int nDisc,double a,double b,bool verbose)
{
  _ncMax=ncMax;
  _nDisc=nDisc;
  _a = a;
  _b = b;
  _verbose = verbose;
}

int Chebychev::fit2(AFunction* f,
                   double a ,
                   double b ,
                   double tol )
{
  std::function<double(double)> func = [f](double val){return f->eval(val);};
  return fit(func,a,b,tol);

}

int Chebychev::fit(std::function<double(double)> f, double a, double b, double tol)
{

  _coeffs.resize(_ncMax,0.);

   /* Evaluate the coefficients of the AChebychev approximation */

   _fillCoeffs(f,a,b);

   /* Loop on some discretized samples of the interval */

   int number = 0;
   double incr = (b - a) / (_nDisc + 1);
   double value = a;
   while(value <= b)
   {
     int numloc = _countCoeffs(f,value,a,b);
     number = MAX(number,numloc);
     value  += incr;
   }


   /* Optional printout  */

   if (_verbose)
   {
     message("AChebychev Polynomial Approximation:\n");
     message("- Performed using %d terms\n",number);
     message("- between %lf and %lf (Nb. discretization steps=%d)\n",a,b,_nDisc);
     message("- with a tolerance of %lg\n",tol);
   }

   /* Core Reallocation */

   _coeffs.resize(number);

   // Optional printout

   if (_verbose)
     for (int i=0; i<(int) _coeffs.size(); i++)
       message("Chebychev coefficient[%d] = %lf\n",i+1,_coeffs[i]);

  return 0;
}

double Chebychev::eval(double x) const
{
  double y, xc, Tx, Tm1, Tm2;

   /* Calculate the approximate value until tolerance is reached */

   xc = 2 * (x - _a) / (_b - _a) - 1.;
   y = _coeffs[0] + _coeffs[1] * xc;
   Tm2 = 1.;
   Tm1 = xc;
   for (int i = 2; i < _ncMax; i++)
   {
     Tx = 2. * xc * Tm1 - Tm2;
     y += _coeffs[i] * Tx;
     Tm2 = Tm1;
     Tm1 = Tx;
   }

   return y;
}


int Chebychev::_countCoeffs(std::function<double(double)> f,double x,double a,double b,double tol)const
{
  double y, y0, y_2, xc, Tx, Tm1, Tm2;

  /* Get the true value  */

  y0 = f(x);
  double y0_2 = y0 * y0;
  double normy0 = y0_2 + EPSILON2;
  tol *= normy0;
  /* Calculate the approximate value until tolerance is reached */
  xc = 2 * (x - a) / (b - a) - 1.;
  y = _coeffs[0] + _coeffs[1] * xc;
  y_2 = y * y;
  if (ABS(y_2 - y0_2) < tol) return (2);
  Tm2 = 1.;
  Tm1 = xc;

  for (int i = 2; i < _ncMax; i++)
  {
    Tx = 2. * Tm1 * xc - Tm2;
    y += _coeffs[i] * Tx;
    if (ABS(y * y - y0 * y0)  < tol)
      return (i + 1);

    Tm2 = Tm1;
    Tm1 = Tx;
  }

  return (_ncMax);
}

void Chebychev::_fillCoeffs(std::function<double(double)> f,double a, double b)
{
  VectorDouble coeffs, x1, y1, x2, y2;
  int n;

  /* Initializations */

  double minsubdiv = pow(2., 20.);
  if (minsubdiv >= (_ncMax + 1) / 2)
    n = static_cast<int> (minsubdiv);
  else
    n = static_cast<int> (ceil((double) (_ncMax + 1) / 2));

  /* Core allocation */

  x1.resize(n);
  y1.resize(n);
  x2.resize(n);
  y2.resize(n);

  /* Filling the arrays */

  for (int i = 0; i < n; i++)
  {
    double theta = 2. * GV_PI * ((double) i) / ((double) n);
    double ct    = cos(theta / 2.);
    double val1  = f(((b + a) + (b - a) * ct) / 2.);
    double val2  = f(((b + a) - (b - a) * ct) / 2.);
    x1[i] = 0.5 * (val1 + val2);
    y1[i] = 0.;
    x2[i] = 0.5 * (val1 - val2) * cos(-theta / 2.);
    y2[i] = 0.5 * (val1 - val2) * sin(-theta / 2.);
  }

  /* Perform the FFT transform */

  if (fftn(1, &n, x1.data(), y1.data(),  1, 1.))
    my_throw("Problem when calculating 1-D Fast Fourrier Transform");
  if (fftn(1, &n, x2.data(), y2.data(), -1, 1.))
    my_throw("Problem when calculating 1-D Fast Fourrier Transform");

  /* Store the coefficients */

  double value = 2. / (double) n;
  for (int i = 0; i < n; i++)
  {
    if (2 * i >= _ncMax) break;
    _coeffs[2 * i    ] = value * x1[i];
    if (2 * i + 1 >= _ncMax) break;
    _coeffs[2 * i + 1] = value * x2[i];
  }
  _coeffs[0] /= 2.;
}


void Chebychev::evalOp(const ALinearOpMulti* Op,const VectorVectorDouble& inv, VectorVectorDouble& outv) const
{
  double v1 = 2. / (_b - _a);
  double v2 = -(_b + _a) / (_b - _a);
  // Initialization


  VectorVectorDouble* tm2 = &Op->_z;
  VectorVectorDouble* tm1 = &Op->_temp;
  VectorVectorDouble* t0 = &Op->_p;
  VectorVectorDouble* swap;

  Op->_copyVals(inv,*tm2);
  // tm1 = v1 Op tm2 + v2 tm2
    Op->evalDirect(*tm2,*tm1);
    Op->_linearComb(v1,*tm1,v2,*tm2,*tm1);

    Op->_linearComb(_coeffs[0],*tm2,_coeffs[1],*tm1,outv);


   /* Loop on the AChebychev polynomials */
    // Op * = 2
    v1 *= 2.;
    v2 *= 2.;

    for (int ib=2; ib<(int) _coeffs.size(); ib++)
    {
      // t0 = (v1 Op + v2 I) tm1
      Op->evalDirect(*tm1,*t0);
      Op->_linearComb(v1, *t0, v2, *tm1, *t0);

      // t0 = 2 * t0 - tm2

      Op->diff(*tm2,*t0,*t0);

      // outv += coeff * y
      Op->addProdScalar(_coeffs[ib],*t0,outv);

      // swap
      swap = tm2;
      tm2 = tm1;
      tm1 = t0;
      t0 = swap;
    }
}
#ifndef SWIG
void Chebychev::evalOp(cs* S,const VectorDouble& x,VectorDouble& y) const
{
  VectorDouble tm1, tm2, px, tx;
  int nvertex;
  cs *T1;

/* Initializations */

  if (!_isReady())
    my_throw("You must use 'initCoeffs' before 'operate'");
  nvertex = cs_getncol(S);
  double v1 = 2. / (_b - _a);
  double v2 = -(_b + _a) / (_b - _a);
  tm1.resize(nvertex);
  tm2.resize(nvertex);
  px.resize(nvertex);
  tx.resize(nvertex);

/* Create the T1 sparse matrix */

  T1 = cs_eye(nvertex, 1.);
  if (T1 == nullptr) my_throw("Problem in cs_eye");
  T1 = cs_add_and_release(T1, S, v2, v1, 1);
  if (T1 == nullptr) my_throw("Problem in cs_add");

/* Initialize the simulation */

  for (int i=0; i<nvertex; i++)
  {
    tm1[i] = 0.;
    y[i]   = x[i];
  }
  if (! cs_gaxpy(T1, y.data(), tm1.data())) my_throw("Problem in cs_gaxpy");
  for (int i=0; i<nvertex; i++)
  {
    px[i]  = _coeffs[0] * y[i] + _coeffs[1] * tm1[i];
    tm2[i] = y[i];
  }

  /* Loop on the AChebychev polynomials */

  for (int ib=2; ib<(int) _coeffs.size(); ib++)
  {
    cs_mulvec(T1, nvertex, tm1.data(), tx.data());
    for (int i=0; i<nvertex; i++)
    {
      tx[i]  = 2. * tx[i] - tm2[i];
      px[i] += _coeffs[ib] * tx[i];
      tm2[i] = tm1[i];
      tm1[i] = tx[i];
    }
  }

/* Return the results */

  for (int i=0; i<nvertex; i++)
    y[i] = px[i];
}
#endif
