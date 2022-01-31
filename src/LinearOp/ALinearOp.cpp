/******************************************************************************/
/* COPYRIGHT ARMINES, ALL RIGHTS RESERVED                                     */
/*                                                                            */
/* THE CONTENT OF THIS WORK CONTAINS CONFIDENTIAL AND PROPRIETARY             */
/* INFORMATION OF ARMINES. ANY DUPLICATION, MODIFICATION,                     */
/* DISTRIBUTION, OR DISCLOSURE IN ANY FORM, IN WHOLE, OR IN PART, IS STRICTLY */
/* PROHIBITED WITHOUT THE PRIOR EXPRESS WRITTEN PERMISSION OF ARMINES         */
/*                                                                            */
/* TAG_SOURCE_CG                                                              */
/******************************************************************************/
#include "geoslib_f.h"
#include "geoslib_old_f.h"

#include "LinearOp/ALinearOp.hpp"
#include "LinearOp/Identity.hpp"
#include "Basic/AException.hpp"
#include "Basic/OptDbg.hpp"

#include <iostream>

ALinearOp::ALinearOp()
  : _nIterMax(100)
  , _eps(1.e-06)
  , _x0()
  , _precondStatus(0)
  , _precond(nullptr)
{
}

ALinearOp::~ALinearOp() 
{
}

/*****************************************************************************/
/*!
**  Evaluate the product: 'out' = Q * 'in'
**
** \param[in]  in     Array of input values
**
** \param[out] out    Array of output values
**
*****************************************************************************/
void ALinearOp::evalDirect(const VectorDouble& in, VectorDouble& out) const
{
  try
  {
    _evalDirect(in,out);
  }
  catch(const char * str)
  {
    // TODO : Check if std::exception can be used
    std::cout << str << std::endl;
  }
}

/*****************************************************************************/
/*!
**  Evaluate the product: 'out' = Q^{-1} * 'in'
**
** \param[in]  in     Array of input values
**
** \param[out] out    Array of output values
**
*****************************************************************************/
void ALinearOp::evalInverse(const VectorDouble& in,
                            VectorDouble& out) const
{
  int n = getSize();
  if (n <= 0) my_throw("ALinearOp size not defined. Call setSize before");

	VectorDouble z    = VectorDouble(n);
	VectorDouble r    = VectorDouble(n);
	VectorDouble temp = VectorDouble(n);
	VectorDouble p    = VectorDouble(n);

	if (! _x0.empty())
    for (int i=0; i<n; i++) out[i] = _x0[i];
  else
    for (int i=0; i<n; i++) out[i] = 0.;

  evalDirect(out,temp);
	for(int i=0; i<n; i++) r[i] = in[i] - temp[i];
  
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
			out[i] += alpha * p[i];
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
double ALinearOp::_prod(const VectorDouble& x,
                        const VectorDouble& y) const
{
  double prod = 0.;
  int n = getSize();
  
  for (int i=0; i<n; i++)
    prod += x[i] * y[i];
  return prod;
}
