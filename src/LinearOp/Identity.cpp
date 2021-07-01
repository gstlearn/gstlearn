/******************************************************************************/
/* COPYRIGHT ARMINES, ALL RIGHTS RESERVED                                     */
/*                                                                            */
/* THE CONTENT OF THIS WORK CONTAINS CONFIDENTIAL AND PROPRIETARY             */
/* INFORMATION OF ARMINES. ANY DUPLICATION, MODIFICATION,                     */
/* DISTRIBUTION, OR DISCLOSURE IN ANY FORM, IN WHOLE, OR IN PART, IS STRICTLY */
/* PROHIBITED WITHOUT THE PRIOR EXPRESS WRITTEN PERMISSION OF ARMINES         */
/*                                                                            */
/* Created on: 9 avr. 2019 by N. Desassis                                     */
/*                                                                            */
/* TAG_SOURCE_CG                                                              */
/******************************************************************************/
#include <LinearOp/Identity.hpp>
#include <iostream>

Identity::Identity(int n) 
  : ALinearOp()
  , _n(n)
{
}

Identity::~Identity() 
{
}

/*****************************************************************************/
/*!
**  Evaluate the product (by the Identity) : 'out' = I * 'in' = 'in'
**
** \param[in]  in     Array of input values
**
** \param[out] out    Array of output values
**
*****************************************************************************/
void Identity::_evalDirect(const VectorDouble& in, VectorDouble& out) const
{
	for(int i=0, n=_n; i<n; i++)
		out[i] = in[i];
}

void Identity::evalInverse(const VectorDouble& in,
                           VectorDouble& out) const
{
  evalDirect(in,out);
}
