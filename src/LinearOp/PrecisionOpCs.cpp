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
#include "Basic/AStringable.hpp"
#include "LinearOp/AShiftOp.hpp"
#include "geoslib_define.h"

#include "Covariances/CovAniso.hpp"
#include "LinearOp/PrecisionOpCs.hpp"
#include "LinearOp/ShiftOpMatrix.hpp"
#include "Polynomials/APolynomial.hpp"
#include "Polynomials/ClassicalPolynomial.hpp"
#include "Mesh/AMesh.hpp"
#include "Basic/AException.hpp"

PrecisionOpCs::PrecisionOpCs(ShiftOpMatrix* shiftop,
                             const CovAniso* cova,
                             bool verbose)
  : PrecisionOp(shiftop, cova, verbose)
  , _Q(nullptr)
  , _chol(nullptr)
{
  _buildQ();
}

PrecisionOpCs::PrecisionOpCs(const AMesh* mesh, CovAniso* cova, bool verbose)
  : PrecisionOp(mesh, cova, verbose)
  , _Q(nullptr)
  , _chol(nullptr)
{
  _buildQ();
}

PrecisionOpCs::~PrecisionOpCs()
{
  delete _Q;
  delete _chol;
}

void PrecisionOpCs::gradYQX(const constvect X,
                            const constvect Y,
                            vect result,
                            const EPowerPT& power)
{
  if (_work2.size() == 0) _work2.resize(getSize());
  if (_work3.size() == 0) _work3.resize(getSize());
  if (_work4.size() == 0) _work4.resize(getSize());

  vect w2s(_work2);
  vect w3s(_work3);
  vect w4s(_work4);
  evalPower(X, w3s, power);
  evalPower(Y, w4s, power);
  double temp, val;
  int iadress;
  const AShiftOp* sopt = getShiftOp();
  const auto *soptmat = dynamic_cast<const ShiftOpMatrix*>(sopt);
  if (soptmat == nullptr)
  {
    messerr("Method only available for ShiftOpMatrix\n");
    return;
  }

  for (int igparam = 0; igparam < soptmat->getNCovAnisoGradParam();
       igparam++)
  {
    for (int iapex = 0; iapex < getSize(); iapex++)
    {
      iadress = soptmat->getSGradAddress(iapex, igparam);
      if (igparam < soptmat->getLambdaGradSize()) // range parameters
      {
        val  = getShiftOp()->getLambda(iapex);
        temp = soptmat->getLambdaGrad(igparam, iapex);
        result[iadress] =
          (X[iapex] * _work4[iapex] + Y[iapex] * _work3[iapex]) * temp / val;
      }
      else
      {
        result[iadress] = 0.;
      }
      evalDeriv(X, w2s, iapex, igparam, power);
      for (int i = 0; i < getSize(); i++)
      {
        result[iadress] += _work2[i] * Y[i];
      }
    }
  }
}

void PrecisionOpCs::gradYQXOptim(const constvect X,
                                 const constvect Y,
                                 vect result,
                                 const EPowerPT& power)
{
  if (_work2.size() == 0) _work2.resize(getSize());
  if (_work3.size() == 0) _work3.resize(getSize());
  if (_work4.size() == 0) _work4.resize(getSize());

  vect w2s(_work2);
  vect w3s(_work3);
  vect w4s(_work4);
  setTraining(false);
  evalPower(Y, w3s, power);
  setTraining(true);
  evalPower(X, w4s, power);

  double temp, val;
  int iadress;

  const AShiftOp* sopt = getShiftOp();
  const auto *soptmat = dynamic_cast<const ShiftOpMatrix*>(sopt);
  if (soptmat == nullptr)
  {
    messerr("Method only available for ShiftOpMatrix\n");
    return;
  }
  for (int igparam = 0; igparam < soptmat->getNCovAnisoGradParam();
       igparam++)
  {
    for (int iapex = 0; iapex < getSize(); iapex++)
    {
      iadress         = soptmat->getSGradAddress(iapex, igparam);
      result[iadress] = 0.;
      if (igparam < soptmat->getLambdaGradSize())
      {
        val  = getShiftOp()->getLambda(iapex);
        temp = soptmat->getLambdaGrad(igparam, iapex);
        result[iadress] =
          (Y[iapex] * _work4[iapex] + X[iapex] * _work3[iapex]) * temp / val;
      }

      evalDerivOptim(w2s, iapex, igparam, power);
      for (int i = 0; i < getSize(); i++)
      {
        result[iadress] += _work2[i] * Y[i];
      }
    }
  }
}

int PrecisionOpCs::_addToDest(const constvect inv, vect outv) const
{
  return _Q->addToDest(inv, outv);
}

int PrecisionOpCs::_addSimulateToDest(const constvect whitenoise,
                                      vect outv) const
{
  if (_chol == nullptr) _chol = new CholeskySparse(_Q);
  _chol->addSimulateToDest(whitenoise, outv);
  return 0;
}

void PrecisionOpCs::evalInverse(const constvect vecin,
                                std::vector<double>& vecout)
{
  if (_chol == nullptr) _chol = new CholeskySparse(_Q);
  _chol->solve(vecin, vecout);
}

double PrecisionOpCs::getLogDeterminant(int nbsimu)
{
  DECLARE_UNUSED(nbsimu);
  if (_chol == nullptr) _chol = new CholeskySparse(_Q);
  return _chol->computeLogDeterminant();
}

void PrecisionOpCs::evalDeriv(
  const constvect inv, vect outv, int iapex, int igparam, const EPowerPT& power)
{
  DECLARE_UNUSED(iapex,igparam)
  if (_work.size()==0) _work.resize(getSize());

  if (power == EPowerPT::MINUSONE)
  my_throw("'evalDeriv' is not yet implemented for 'EPowerPT::MINUSONE'");
  if (power == EPowerPT::MINUSHALF)
  my_throw("'evalDeriv' is not yet implemented for 'EPowerPT::MINUSHALF'");
  if (power == EPowerPT::LOG)
  my_throw("'evalDeriv' is not yet implemented for 'EPowerPT::LOG'");

  // Pre-processing
  vect ws(_work);
  getShiftOp()->prodLambda(inv, ws, EPowerPT::ONE);

  // Polynomial evaluation

  // ((ClassicalPolynomial*) getPoly(power))->evalDerivOp(getShiftOp(), ws,
  //                                                      outv, iapex, igparam);

  // Post-processing

  getShiftOp()->prodLambda(outv, outv, EPowerPT::ONE);
}

void PrecisionOpCs::evalDerivOptim(vect outv,
                                   int iapex,
                                   int igparam,
                                   const EPowerPT& power)
{
  DECLARE_UNUSED(iapex,igparam)
  if (_work.size()  == 0) _work3.resize(getSize());
  if (_work5.size() == 0) _work5.resize(getSize());

  vect ws(_work);
  vect w5s(_work5);
  if (power == EPowerPT::MINUSONE)
  my_throw("'evalDeriv' is not yet implemented for 'POPT_MINUSONE'");
  if (power == EPowerPT::MINUSHALF)
  my_throw("'evalDeriv' is not yet implemented for 'POPT_MINUSHALF'");
  if (power == EPowerPT::LOG)
  my_throw("'evalDeriv' is not yet implemented for 'POPT_LOG'");

  // ((ClassicalPolynomial*) getPoly(power))->evalDerivOpOptim(getShiftOp(), ws,
  //                                                           w5s, outv,
  //                                                           _workPoly, iapex,
  //                                                           igparam);

  // Post-processing
  getShiftOp()->prodLambda(outv, outv, EPowerPT::ONE);
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
  delete _Q;
  if (!isCovaDefined()) return;

  // Calculate the Vector of coefficients (blin)
  //VectorDouble blin = getPoly(EPowerPT::ONE)->getCoeffs();

  // Calculate the Precision matrix Q

  _Q = _build_Q();
}

/****************************************************************************/
/*!
 **  Construct the final sparse matrix Q from the Model
 **
 ** \return Error return code
 **
 ** \param[in] S        Shift operator
 ** \param[in] Lambda   Lambda vector
 ** \param[in] nblin    Number of blin coeffbuicients
 ** \param[in] blin     Array of coefficients for Linear combinaison
 **
 *****************************************************************************/
MatrixSparse* PrecisionOpCs::_build_Q()
{
  // Preliminary checks
  auto *S = ((ShiftOpMatrix*)getShiftOp())->getS();
  auto Lambda = getShiftOp()->getLambdas();
  VectorDouble blin = getPoly(EPowerPT::ONE)->getCoeffs();
  int nblin = static_cast<int>(blin.size());
  int nvertex = S->getNCols();
  if (nvertex <= 0)
  {
    messerr("You must define a valid Meshing beforehand");
    return nullptr;
  }
  if (nblin <= 0)
  {
    messerr("You must have a set of already available 'blin' coefficients");
    messerr("These coefficients come from the decomposition in series for Q");
    messerr("This decomposition is available only if 'alpha' is an integer");
    messerr("where: alpha = param + ndim/2");
    return nullptr;
  }

  /* First step */

  MatrixSparse* Q  = MatrixSparse::diagConstant(nvertex, blin[0]);
  MatrixSparse* Bi = S->clone();

  /* Loop on the different terms */

  for (int iterm = 1; iterm < nblin; iterm++)
  {
    Q->addMatInPlace(*Bi, 1., blin[iterm]);
    if (iterm < nblin - 1) Bi->prodMatInPlace(S);
  }
  delete Bi;

  /* Final scaling */

  Q->prodNormDiagVecInPlace(Lambda, 1);
  return Q;
}

VectorDouble PrecisionOpCs::extractDiag() const
{
  return _Q->extractDiag();
}
