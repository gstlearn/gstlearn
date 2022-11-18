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
**  Evaluate the product (by the Identity) : 'outv' = I * 'inv' = 'inv'
**
** \param[in]  inv     Array of input values
**
** \param[out] outv    Array of output values
**
*****************************************************************************/
void Identity::_evalDirect(const VectorDouble& inv, VectorDouble& outv) const
{
  for(int i=0, n=_n; i<n; i++)
    outv[i] = inv[i];
}

void Identity::evalInverse(const VectorDouble& inv,
                           VectorDouble& outv) const
{
  evalDirect(inv,outv);
}
