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
#include "Basic/AException.hpp"
#include "LinearOp/PrecisionOpCs.hpp"
#include "Polynomials/APolynomial.hpp"
#include "Basic/Vector.hpp"
#include "Model/Model.hpp"
#include "geoslib_e.h"
#include "csparse_d.h"
#include "LinearOp/ShiftOpCs.hpp"
#include "Polynomials/ClassicalPolynomial.hpp"

PrecisionOpCs::PrecisionOpCs(ShiftOpCs* shiftop,
                             const CovAniso* cova,
                             ENUM_POPTS power,
                             bool verbose)
    : PrecisionOp(shiftop, cova, power, verbose)
{
}

PrecisionOpCs::~PrecisionOpCs()
{
  // TODO Auto-generated destructor stub
}

VectorDouble PrecisionOpCs::getCoeffs()
{
  VectorDouble coeffs = getPoly(POPT_ONE)->getCoeffs();
  return coeffs;
}


void PrecisionOpCs::gradYQX(const VectorDouble & X, const VectorDouble &Y,VectorDouble& result)
{
  if (_work2.empty()) _work2.resize(getSize());
  if (_work3.empty()) _work3.resize(getSize());
  eval(X,_work3);
  double temp,val;
  int iadress;

  for(int igparam = 0;igparam<getShiftOp()->getNModelGradParam();igparam++)
  {
    for(int iapex=0;iapex<getSize();iapex++)
    {
      iadress = getShiftOp()->getSGradAddress(iapex,igparam);
      if(igparam < getShiftOp()->getDim()) // range parameters
      {
        val = getShiftOp()->getLambda(iapex);
        temp = getShiftOp()->getLambdaGrad(igparam,iapex);
        result[iadress]= 2 * Y[iapex] * temp * _work3[iapex] / val;

      }
      else
      {
        result[iadress] = 0.;
      }
      evalDeriv(X,_work2,iapex,igparam);
      for(int i = 0;i<getSize();i++)
      {
        result[iadress] += _work2[i]*Y[i];
      }
    }
  }
}

void PrecisionOpCs::evalDeriv(const VectorDouble& in, VectorDouble& out,int iapex,int igparam)
{
  const VectorDouble* inPtr = &in;
  if (_work.empty()) _work.resize(getSize());

  if(getPower() == POPT_MINUSONE)
     my_throw("'evalDeriv' is not yet implemented for 'POPT_MINUSONE'");
  if(getPower() == POPT_MINUSHALF)
     my_throw("'evalDeriv' is not yet implemented for 'POPT_MINUSHALF'");
  if(getPower() == POPT_LOG)
     my_throw("'evalDeriv' is not yet implemented for 'POPT_LOG'");

  // Pre-processing

  if (getPower() == POPT_ONE)
  {
    getShiftOp()->prodTildeC(in, _work, POPT_HALF);
    inPtr = &_work;
  }
  else if (getPower() == POPT_MINUSONE)
  {
    getShiftOp()->prodTildeC(in, _work, POPT_MINUSHALF);
    inPtr = &_work;

  }

  // Polynomial evaluation


    ((ClassicalPolynomial*)getPoly(getPower()))->evalDerivOp(getShiftOp(),
                                                             *inPtr,
                                                             out,
                                                             iapex,
                                                             igparam);

    // Post-processing

    if (getPower() == POPT_ONE)
    {
       getShiftOp()->prodTildeC(out, out, POPT_HALF);
       getShiftOp()->prodLambdaOnSqrtTildeC(out, out, 2.);
    }
    else if (getPower() == POPT_MINUSONE)
    {
      getShiftOp()->prodTildeC(out, out, POPT_MINUSHALF);
      getShiftOp()->prodLambdaOnSqrtTildeC(out, out, -2.);
    }
    else if (getPower() == POPT_MINUSHALF)
    {
      getShiftOp()->prodLambda(out, out, POPT_MINUSONE);
    }

}

void PrecisionOpCs::evalDerivPoly(const VectorDouble& in, VectorDouble& out,int iapex,int igparam)
{

  if(getPower() == POPT_ONE)
     my_throw("'evalDeriv' is not yet implemented for 'POPT_ONE'");
  if(getPower() == POPT_MINUSONE)
     my_throw("'evalDeriv' is not yet implemented for 'POPT_MINUSONE'");
  if(getPower() == POPT_MINUSHALF)
     my_throw("'evalDeriv' is not yet implemented for 'POPT_MINUSHALF'");
  if(getPower() == POPT_LOG)
     my_throw("'evalDeriv' is not yet implemented for 'POPT_LOG'");

}

cs *PrecisionOpCs::getQ()
{
  VectorDouble blin = getPoly(POPT_ONE)->getCoeffs();
  cs* Q = spde_build_Q(getShiftOp()->getS(), getShiftOp()->getLambda(),
                       static_cast<int> (blin.size()), blin.data());
  return Q;
}
