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
#include "Basic/Utilities.hpp"
#include "Basic/VectorHelper.hpp"
#include "Basic/AException.hpp"
#include "Basic/Law.hpp"
#include "LinearOp/PrecisionOp.hpp"
#include "LinearOp/ShiftOpCs.hpp"
#include "Polynomials/APolynomial.hpp"
#include "Polynomials/ClassicalPolynomial.hpp"
#include "Polynomials/Chebychev.hpp"
#include "Covariances/CovAniso.hpp"
#include "Model/Model.hpp"
#include "Mesh/AMesh.hpp"

#include <math.h>

PrecisionOp::PrecisionOp()
  : _shiftOp(nullptr)
  , _cova(nullptr)
  , _polynomials()
  , _verbose(false)
  , _training(false)
  , _destroyShiftOp(false)
  , _userPoly(false)
  , _work()
  , _work2()
  , _work3()
{
}

PrecisionOp::PrecisionOp(ShiftOpCs* shiftop,
                         const CovAniso* cova,
                         bool verbose)
  : _shiftOp(shiftop)
  , _cova(cova)
  , _polynomials()
  , _verbose(verbose)
  , _training(false)
  , _destroyShiftOp(false)
  , _userPoly(false)
  , _work()
  , _work2()
  , _work3()
{
  if (_shiftOp != nullptr)
  {
    _work.resize(_shiftOp->getSize());
    _work2.resize(_shiftOp->getSize());
    _work3.resize(_shiftOp->getSize());
  }
}

PrecisionOp::PrecisionOp(const AMesh* mesh,
                         Model* model,
                         int icov,
                         bool verbose)
  : _shiftOp(nullptr)
  , _cova(model->getCova(icov))
  , _polynomials()
  , _verbose(verbose)
  , _training(false)
  , _destroyShiftOp(true)
  , _userPoly(false)
  , _work()
  , _work2()
  , _work3()
{
  _shiftOp = new ShiftOpCs(mesh,model,nullptr,0,icov,verbose);

  _work.resize(_shiftOp->getSize());
  _work2.resize(_shiftOp->getSize());
  _work3.resize(_shiftOp->getSize());
}

PrecisionOp::PrecisionOp(const PrecisionOp &pmat)
  : _shiftOp(pmat._shiftOp)
  , _cova(pmat._cova)
  , _polynomials(pmat._polynomials)
  , _verbose(pmat._verbose)
  , _training(false)
  , _destroyShiftOp(pmat._destroyShiftOp)
  , _userPoly(false)
  , _work(pmat._work)
  , _work2(pmat._work2)
  , _work3(pmat._work3)
{
}

PrecisionOp& PrecisionOp::operator= (const PrecisionOp &pmat)
{
  if (this != &pmat)
  {
    _shiftOp = pmat._shiftOp;
    _cova = pmat._cova;
    _polynomials = pmat._polynomials;
    _verbose = pmat._verbose;
    _training = pmat._training;
    _destroyShiftOp = pmat._destroyShiftOp;
    _userPoly = pmat._userPoly;
    _work = pmat._work;
    _work2 = pmat._work2;
    _work3 = pmat._work3;
  }
  return *this;
}

void PrecisionOp::_purge()
{
  for (auto& e: _polynomials)
  {
    if(e.first!=EPowerPT::ONE || !_userPoly)
      delete e.second;
  }
  _polynomials.clear();

}

PrecisionOp::~PrecisionOp()
{
  _purge();
  if (_destroyShiftOp)
  {
    delete _shiftOp;
  }
}

PrecisionOp* PrecisionOp::createFromShiftOp(ShiftOpCs *shiftop,
                                            const CovAniso *cova,
                                            bool verbose)
{
  return new PrecisionOp(shiftop, cova, verbose);
}

PrecisionOp* PrecisionOp::create(const AMesh *mesh,
                           Model *model,
                           int icov,
                           bool verbose)
{
  return new PrecisionOp(mesh, model, icov, verbose);
}

int PrecisionOp::_preparePoly(const EPowerPT& power,bool force)
{
  // Polynomial already exists. Nothing to be done
  if (_polynomials.count(power) && !force) return 0;

  // Prepare Polynomial for EPowerPT::ONE
  if (_preparePrecisionPoly() && !force) return 1;

  // Prepare polynomials for other powers than 1
  if (power != EPowerPT::ONE)
  {
    if (_prepareChebychev(power)) return 1;
  }
  return 0;
}

VectorDouble PrecisionOp::getPolyCoeffs(EPowerPT power)
{
  return _polynomials[power]->getCoeffs();
}

void PrecisionOp::setPolynomialFromPoly(APolynomial* polynomial)
  {
    _purge();
    _userPoly = true;
    _polynomials[EPowerPT::ONE] = polynomial;
    _preparePoly(EPowerPT::MINUSONE,true);
    _preparePoly(EPowerPT::MINUSHALF,true);
    _preparePoly(EPowerPT::LOG,true);
  }

int PrecisionOp::_prepareChebychev(const EPowerPT& power)
{
  if (_cova == nullptr && _polynomials.count(EPowerPT::ONE)==0) return 1;
  if (_shiftOp == nullptr) return 1;

  double b = _shiftOp->getMaxEigenValue();
  APolynomial* chebMatern = new Chebychev();
  ((Chebychev*)chebMatern)->setA(0);
  ((Chebychev*)chebMatern)->setB(b);

  std::function<double(double)> f;

  if(power == EPowerPT::LOG)
  {
    f = [this](double val){return log(_polynomials[EPowerPT::ONE]->eval(val));};
  }
  else if (power == EPowerPT::MINUSONE)
  {
    f = [this](double val){return pow(_polynomials[EPowerPT::ONE]->eval(val),-1.);};
  }
  else if (power == EPowerPT::MINUSHALF)
  {
    f = [this](double val){return pow(_polynomials[EPowerPT::ONE]->eval(val),-0.5);};
  }

  chebMatern->fit(f,0,b);
  _polynomials[power] = chebMatern;
  return 0;
}

/**
 * Compute the Logarithm of the Determinant
 * @param nsimus
 * @param seed
 * @return The computed value or TEST if problem
 */
double PrecisionOp::computeLogDet(int nsimus,int seed)
{
  law_set_random_seed(seed);

  VectorDouble gauss;
  VectorDouble result;
  gauss.resize(getSize());
  result.resize(getSize());

  double val1 = 0.;
  for (int isimu = 0; isimu < nsimus; isimu++)
  {
    for (auto &e : gauss)
    {
      e = law_gaussian();
    }
    if (_evalPoly(EPowerPT::LOG, gauss, result)) return TEST;

    for (int i = 0; i < getSize(); i++)
    {
      val1 += gauss[i] * result[i];
    }
  }

  val1 /= nsimus;

  double val2 = 0.;

  for (auto &e : _shiftOp->getLambdas())
  {
    val2 += log(e);
  }

  val1 += 2. * val2;

  return val1;
}

int PrecisionOp::reset(const ShiftOpCs* shiftop,
                       const CovAniso* cova,
                       bool verbose)
{
  // Initializations

  int error = 0;

  try
  {
    // Store the pointer to the ShiftOp
    if (shiftop == (ShiftOpCs *) NULL)
      my_throw("The argument 'shiftop'must be provided");

    // Store the members

    _cova    = cova;
    _verbose = verbose;
    _shiftOp = new ShiftOpCs(*shiftop);

    _purge();
  }

  catch(const char * str)
  {
    error = 1;
    messerr("%s", str);
  }
  return error;
}

/**
 * Evaluate with power = ONE
 * @param inv Input array
 * @param outv Output array
 */
void PrecisionOp::eval(const VectorDouble &inv, VectorDouble &outv)
{
  evalPower(inv, outv, EPowerPT::ONE);
}

void PrecisionOp::evalPower(const VectorDouble& inv, VectorDouble& outv, const EPowerPT& power)
{
  const VectorDouble* inPtr = &inv;
  if (_work.empty()) _work.resize(getSize());

  // Pre-processing

  if (power == EPowerPT::ONE || power == EPowerPT::MINUSONE)
  {
    _shiftOp->prodLambda(inv, _work,power);
    inPtr = &_work;
  }

  // Polynomial evaluation

  if (_evalPoly(power,*inPtr,outv))
    my_throw("Computation in 'eval' interrupted due to problem in '_evalPoly'");

  // Post-processing

  if (power == EPowerPT::ONE || power == EPowerPT::MINUSONE)
  {
    _shiftOp->prodLambda(outv, outv, power);
  }
  else if (power == EPowerPT::MINUSHALF)
  {
    _shiftOp->prodLambda(outv, outv, EPowerPT::MINUSONE);
  }
}

int PrecisionOp::_evalPoly(const EPowerPT& power,
                           const VectorDouble& inv,
                           VectorDouble& outv)
{
  if (_preparePoly(power)) return 1;
  if(getTraining())
  {
    int degree = _polynomials[power]->getDegree();

    if(_workPoly.empty())
    {
      _workPoly = VectorVectorDouble(degree);
      for(auto &e: _workPoly)
      {
        e = VectorDouble(inv.size());
      }
    }

    if (_work5.empty()) _work5.resize(getSize());
    ((ClassicalPolynomial*)_polynomials[power])->evalOpTraining(_shiftOp->getS(),inv,_workPoly,_work5);

    for(int i=0;i<(int)inv.size();i++)
    {
      outv[i] = _workPoly[0][i];
    }
  }
  else
  {
    _polynomials[power]->evalOp(_shiftOp->getS(),inv,outv);
  }
  return 0;
}

void PrecisionOp::evalInvVect(VectorDouble& in, VectorDouble& result)
{
  _shiftOp->prodLambda(in,result,EPowerPT::MINUSONE);
  _evalPoly(EPowerPT::MINUSONE, result, in);
  _shiftOp->prodLambda(in, result, EPowerPT::MINUSONE);
}

VectorDouble PrecisionOp::evalCov(int imesh)
{

  int n = getSize();
  VectorDouble ei(n);
  VectorDouble result(n);
  VH::fill(ei,0.,n);
  ei[imesh] = 1.;
  _shiftOp->prodLambda(ei,result,EPowerPT::MINUSONE);
  _evalPoly(EPowerPT::MINUSONE,result,ei);
  _shiftOp->prodLambda(ei, result, EPowerPT::MINUSONE);

  return result;
}

void PrecisionOp::simulateOneInPlace(VectorDouble& whitenoise, VectorDouble& result)
{
  _evalPoly(EPowerPT::MINUSHALF,whitenoise,result);
  _shiftOp->prodLambda(result, result, EPowerPT::MINUSONE);
}

VectorVectorDouble PrecisionOp::simulate(int nbsimus)
{
  int n = getSize();
  VectorVectorDouble vect(nbsimus);
  VectorDouble whitenoise(n);

  for(auto &e : vect)
  {
    e.resize(n);
    whitenoise = VH::simulateGaussian(n);
    _evalPoly(EPowerPT::MINUSHALF,whitenoise,e);
    _shiftOp->prodLambda(e, e, EPowerPT::MINUSONE);
  }
  return vect;
}

VectorDouble PrecisionOp::simulateOne()
{
  int n = getSize();
  VectorDouble vect(n);
  VectorDouble whitenoise = VH::simulateGaussian(n);
  _evalPoly(EPowerPT::MINUSHALF,whitenoise,vect);
  _shiftOp->prodLambda(vect, vect, EPowerPT::MINUSONE);
  return vect;
}

int PrecisionOp::_preparePrecisionPoly()
{
  if (_cova == nullptr) return 1;
  if (!_cova->hasMarkovCoeffs()) return 1;

  _polynomials[EPowerPT::ONE] = new ClassicalPolynomial(_cova->getMarkovCoeffs());

  if (_polynomials.count(EPowerPT::ONE)) return 0;

  return 0;
}

/* Evaluation of the polynomial of the precision over the interval [0, lambda_max(ShiftOp)] */
/* discretized over ndiscr points */
/* and return the min and the max */

std::pair<double,double> PrecisionOp::getRangeEigenVal(int ndiscr)
{
  std::pair<double,double> rangeVals;

  double sill = _cova->getSill(0,0);
  double sMax = _shiftOp->getMaxEigenValue();
  double x = 0;
  double delta = sMax/(ndiscr-1);

  double val =  _polynomials[EPowerPT::ONE]->eval(x);
  rangeVals.first  = val;
  rangeVals.second = val;

  for(int i = 1; i < ndiscr; i++)
  {
    x += delta;
    val =  _polynomials[EPowerPT::ONE]->eval(x) ;
    rangeVals.first  = MIN(val,rangeVals.first);
    rangeVals.second = MAX(val,rangeVals.first);
  }

  rangeVals.first  = rangeVals.first/sill;
  rangeVals.second = rangeVals.second/sill;

  return rangeVals;
}

APolynomial* PrecisionOp::getPoly(const EPowerPT& power)
{
  if (_preparePoly(power))
    my_throw("Problem in function getPoly");
  return _polynomials[power];
}

void PrecisionOp::setNIterMax(int nitermax)
{
  if (_shiftOp != nullptr)
    _shiftOp->setNIterMax(nitermax);
}

void PrecisionOp::setEps(double eps)
{
  if (_shiftOp != nullptr)
    _shiftOp->setEps(eps);
}

VectorDouble PrecisionOp::getCoeffs()
{
  VectorDouble coeffs = getPoly(EPowerPT::ONE)->getCoeffs();
  return coeffs;
}

