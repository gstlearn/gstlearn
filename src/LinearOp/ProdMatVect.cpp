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
#include <LinearOp/ProdMatVect.hpp>

ProdMatVect::ProdMatVect(int nx,
                         int ny,
                         double* A) 
  : ALinearOp()
  , _nx(nx)
  , _ny(ny)
  , _A(A) 
{
}

ProdMatVect::~ProdMatVect()
{
}

void ProdMatVect::_evalDirect(const VectorDouble& in,
                              VectorDouble& out) const
{
	double s = 0;
	for(int j=0; j<_ny; j++)
		out[j] = 0;

	for(int i=0, nx=_nx; i<nx; i++)
	{
		s = 0;
		for(int j=0; j<_ny; j++)
			s += _A[i+_ny*j] * in[j];
		out[i] = s;
	}
}


