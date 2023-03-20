/******************************************************************************/
/* gstlearn C++ Library                                                       */
/*                                                                            */
/* Authors: <authors>                                                         */
/*                                                                            */
/* License: BSD 3 Clause                                                      */
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

void Identity::evalInverse(const VectorDouble &inv, VectorDouble &outv) const
{
  evalDirect(inv,outv);
}
