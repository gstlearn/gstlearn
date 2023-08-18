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
#include "geoslib_f_private.h"
#include "Basic/AException.hpp"
#include "LinearOp/PrecisionOpCs.hpp"
#include "Polynomials/APolynomial.hpp"
#include "Model/Model.hpp"
#include "Mesh/AMesh.hpp"
#include "LinearOp/ShiftOpCs.hpp"
#include "Model/Model.hpp"
#include "Polynomials/ClassicalPolynomial.hpp"
#include "Matrix/LinkMatrixSparse.hpp"

// External library /// TODO : Dependency to csparse to be removed
#include "csparse_d.h"
#include "csparse_f.h"

PrecisionOpCs::PrecisionOpCs(ShiftOpCs* shiftop,
                             const CovAniso* cova,
                             bool verbose)
    : PrecisionOp(shiftop, cova, verbose),
      _Q(nullptr),
      _qChol()
{
  _buildQ();
}

PrecisionOpCs::PrecisionOpCs(const AMesh* mesh,
                             Model* model,
                             int icov,
                             bool verbose)
    : PrecisionOp(mesh, model, icov, verbose),
      _Q(nullptr),
      _qChol()
{
  _buildQ();
}

PrecisionOpCs::~PrecisionOpCs()
{
  _Q = cs_spfree(_Q);
}

Triplet PrecisionOpCs::getQToTriplet(bool flag_from_1) const
{
  return csToTriplet(getQ(), flag_from_1);
}

void PrecisionOpCs::gradYQX(const VectorDouble & X, const VectorDouble &Y, VectorDouble& result, const EPowerPT& power)
{
  if (_work2.empty()) _work2.resize(getSize());
  if (_work3.empty()) _work3.resize(getSize());
  if (_work4.empty()) _work4.resize(getSize());

  evalPower(X,_work3, power);
  evalPower(Y,_work4, power);
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
      evalDeriv(X,_work2,iapex,igparam, power);
      for(int i = 0;i<getSize();i++)
      {
        result[iadress] += _work2[i]*Y[i];
      }
    }
  }
}


void PrecisionOpCs::gradYQXOptim(const VectorDouble & X, const VectorDouble &Y,VectorDouble& result, const EPowerPT& power)
{
  if (_work2.empty()) _work2.resize(getSize());
  if (_work3.empty()) _work3.resize(getSize());
  if (_work4.empty()) _work4.resize(getSize());

  setTraining(false);
  evalPower(Y,_work3, power);
  setTraining(true);
  evalPower(X,_work4, power);

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

      evalDerivOptim(_work2,iapex,igparam,power);
      for(int i = 0;i<getSize();i++)
      {
        result[iadress] += _work2[i]*Y[i];
      }
    }
  }
}

void PrecisionOpCs::eval(const VectorDouble &inv, VectorDouble &outv)
{
  _qChol.evalDirect(inv, outv);
}

void PrecisionOpCs::simulateOneInPlace(VectorDouble& whitenoise, VectorDouble& result)
{
  _qChol.simulate(whitenoise, result);
}

void PrecisionOpCs::evalInvVect(VectorDouble& in, VectorDouble& result)
{
  _qChol.evalInverse(in, result);
}

double PrecisionOpCs::computeLogDet(int nsimus, int seed)
{
  DECLARE_UNUSED(nsimus);
  DECLARE_UNUSED(seed);

  return _qChol.computeLogDet();
}

void PrecisionOpCs::evalDeriv(const VectorDouble& inv, VectorDouble& outv,int iapex,int igparam, const EPowerPT& power)
{
  if (_work.empty()) _work.resize(getSize());

  if(power == EPowerPT::MINUSONE)
     my_throw("'evalDeriv' is not yet implemented for 'EPowerPT::MINUSONE'");
  if(power == EPowerPT::MINUSHALF)
     my_throw("'evalDeriv' is not yet implemented for 'EPowerPT::MINUSHALF'");
  if(power == EPowerPT::LOG)
     my_throw("'evalDeriv' is not yet implemented for 'EPowerPT::LOG'");

  // Pre-processing

    getShiftOp()->prodLambda(inv, _work ,EPowerPT::ONE);

  // Polynomial evaluation


    ((ClassicalPolynomial*)getPoly(power))->evalDerivOp(getShiftOp(),
                                                             _work,
                                                             outv,
                                                             iapex,
                                                             igparam);

    // Post-processing

       getShiftOp()->prodLambda(outv, outv ,EPowerPT::ONE);
}

void PrecisionOpCs::evalDerivOptim(VectorDouble& outv,
                                   int iapex,
                                   int igparam,
                                   const EPowerPT& power)
{
  if (_work.empty()) _work3.resize(getSize());
  if (_work5.empty()) _work4.resize(getSize());

  if(power == EPowerPT::MINUSONE)
     my_throw("'evalDeriv' is not yet implemented for 'POPT_MINUSONE'");
  if(power == EPowerPT::MINUSHALF)
     my_throw("'evalDeriv' is not yet implemented for 'POPT_MINUSHALF'");
  if(power == EPowerPT::LOG)
     my_throw("'evalDeriv' is not yet implemented for 'POPT_LOG'");


  ((ClassicalPolynomial*) getPoly(power))->evalDerivOpOptim(
      getShiftOp(), _work,_work5,outv,_workPoly, iapex, igparam);

    // Post-processing
       getShiftOp()->prodLambda(outv, outv ,EPowerPT::ONE);
}


//void PrecisionOpCs::evalDerivPoly(const VectorDouble& inv, VectorDouble& outv,int iapex,int igparam)
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

void PrecisionOpCs::_buildQ()
{
  if (_Q != nullptr) _Q = cs_spfree(_Q);
  if (! isCovaDefined()) return;

  // Calculate the Vector of coefficients (blin)
  VectorDouble blin = getPoly(EPowerPT::ONE)->getCoeffs();

  // Calculate the Precision matrix Q
  _Q = _spde_build_Q(getShiftOp()->getS(), getShiftOp()->getLambdas(),
                     static_cast<int>(blin.size()), blin.data());

  // Store the Precision matrix for Cholesky decomposition
  _qChol.reset(_Q, false);
}
