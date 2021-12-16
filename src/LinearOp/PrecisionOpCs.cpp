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
#include "geoslib_f_private.h"
#include "Basic/AException.hpp"
#include "LinearOp/PrecisionOpCs.hpp"
#include "Polynomials/APolynomial.hpp"
#include "Basic/Vector.hpp"
#include "Model/Model.hpp"
#include "csparse_d.h"
#include "LinearOp/ShiftOpCs.hpp"
#include "Polynomials/ClassicalPolynomial.hpp"

PrecisionOpCs::PrecisionOpCs(ShiftOpCs* shiftop,
                             const CovAniso* cova,
                             const EPowerPT& power,
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
  VectorDouble coeffs = getPoly(EPowerPT::ONE)->getCoeffs();
  return coeffs;
}


void PrecisionOpCs::gradYQX(const VectorDouble & X, const VectorDouble &Y,VectorDouble& result)
{
  if (_work2.empty()) _work2.resize(getSize());
  if (_work3.empty()) _work3.resize(getSize());
  if (_work4.empty()) _work4.resize(getSize());

  eval(X,_work3);
  eval(Y,_work4);
  double temp,val;
  int iadress;

  for(int igparam = 0;igparam<getShiftOp()->getNModelGradParam();igparam++)
  {
    for(int iapex=0;iapex<getSize();iapex++)
    {
      iadress = getShiftOp()->getSGradAddress(iapex,igparam);
      if(igparam < getShiftOp()->getLambdaGradSize()) // range parameters
      {
        val = getShiftOp()->getLambda(iapex);
        temp = getShiftOp()->getLambdaGrad(igparam,iapex);
        result[iadress]= (X[iapex] * _work4[iapex] + Y[iapex] * _work3[iapex]) * temp / val;

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


void PrecisionOpCs::gradYQXOptim(const VectorDouble & X, const VectorDouble &Y,VectorDouble& result)
{
  if (_work2.empty()) _work2.resize(getSize());
  if (_work3.empty()) _work3.resize(getSize());
  if (_work4.empty()) _work4.resize(getSize());

  setTraining(false);
  eval(Y,_work3);
  setTraining(true);
  eval(X,_work4);

  double temp,val;
  int iadress;

  for(int igparam = 0;igparam<getShiftOp()->getNModelGradParam();igparam++)
  {
    for(int iapex=0;iapex<getSize();iapex++)
    {
      iadress = getShiftOp()->getSGradAddress(iapex,igparam);
      result[iadress] = 0.;
      if(igparam < getShiftOp()->getLambdaGradSize())
      {
        val = getShiftOp()->getLambda(iapex);
        temp = getShiftOp()->getLambdaGrad(igparam,iapex);
        result[iadress]=  (Y[iapex] * _work4[iapex] + X[iapex] * _work3[iapex]) * temp / val;

      }

      evalDerivOptim(_work2,iapex,igparam);
      for(int i = 0;i<getSize();i++)
      {
        result[iadress] += _work2[i]*Y[i];
      }
    }
  }
}

void PrecisionOpCs::evalDeriv(const VectorDouble& in, VectorDouble& out,int iapex,int igparam)
{
  if (_work.empty()) _work.resize(getSize());

  if(getPower() == EPowerPT::MINUSONE)
     my_throw("'evalDeriv' is not yet implemented for 'EPowerPT::MINUSONE'");
  if(getPower() == EPowerPT::MINUSHALF)
     my_throw("'evalDeriv' is not yet implemented for 'EPowerPT::MINUSHALF'");
  if(getPower() == EPowerPT::LOG)
     my_throw("'evalDeriv' is not yet implemented for 'EPowerPT::LOG'");

  // Pre-processing

    getShiftOp()->prodLambda(in, _work ,EPowerPT::ONE);

  // Polynomial evaluation


    ((ClassicalPolynomial*)getPoly(getPower()))->evalDerivOp(getShiftOp(),
                                                             _work,
                                                             out,
                                                             iapex,
                                                             igparam);

    // Post-processing

       getShiftOp()->prodLambda(out, out ,EPowerPT::ONE);
}

void PrecisionOpCs::evalDerivOptim(VectorDouble& out,
                                   int iapex,
                                   int igparam)
{
  if (_work.empty()) _work3.resize(getSize());
  if (_work5.empty()) _work4.resize(getSize());

  if(getPower() == EPowerPT::MINUSONE)
     my_throw("'evalDeriv' is not yet implemented for 'POPT_MINUSONE'");
  if(getPower() == EPowerPT::MINUSHALF)
     my_throw("'evalDeriv' is not yet implemented for 'POPT_MINUSHALF'");
  if(getPower() == EPowerPT::LOG)
     my_throw("'evalDeriv' is not yet implemented for 'POPT_LOG'");


  ((ClassicalPolynomial*) getPoly(getPower()))->evalDerivOpOptim(
      getShiftOp(), _work,_work5,out,_workPoly, iapex, igparam);

    // Post-processing
       getShiftOp()->prodLambda(out, out ,EPowerPT::ONE);

}


//void PrecisionOpCs::evalDerivPoly(const VectorDouble& in, VectorDouble& out,int iapex,int igparam)
//{
//
//  if(getPower() == EPowerPT::ONE)
//     my_throw("'evalDeriv' is not yet implemented for 'EPowerPT::ONE'");
//  if(getPower() == EPowerPT::MINUSONE)
//     my_throw("'evalDeriv' is not yet implemented for 'EPowerPT::MINUSONE'");
//  if(getPower() == EPowerPT::MINUSHALF)
//     my_throw("'evalDeriv' is not yet implemented for 'EPowerPT::MINUSHALF'");
//  if(getPower() == EPowerPT::LOG)
//     my_throw("'evalDeriv' is not yet implemented for 'EPowerPT::LOG'");
//
//}

cs *PrecisionOpCs::getQ()
{
  VectorDouble blin = getPoly(EPowerPT::ONE)->getCoeffs();
  cs* Q = _spde_build_Q(getShiftOp()->getS(), getShiftOp()->getLambdas(),
                        static_cast<int> (blin.size()), blin.data());
  return Q;
}
