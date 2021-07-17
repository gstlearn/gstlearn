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
#include "Basic/Utilities.hpp"
#include "LinearOp/PrecisionOp.hpp"
#include "Polynomials/ClassicalPolynomial.hpp"
#include "Polynomials/Chebychev.hpp"
#include "Basic/Vector.hpp"
#include "Basic/AException.hpp"
#include "Covariances/CovAniso.hpp"
#include "LinearOp/ShiftOpCs.hpp"
#include "Basic/Law.hpp"
#include "geoslib_e.h"

PrecisionOp::PrecisionOp(ShiftOpCs* shiftop,
                         const CovAniso* cova,
                         ENUM_POPTS power,
                         bool verbose)
  : _shiftOp(shiftop)
  , _cova(cova)
  , _power(power)
  , _polynomials()
  , _verbose(verbose)
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

PrecisionOp::PrecisionOp(const PrecisionOp &pmat)
  : _shiftOp(pmat._shiftOp)
  , _cova(pmat._cova)
  , _power(pmat._power)
  , _polynomials(pmat._polynomials)
  , _verbose(pmat._verbose)
  , _work(pmat._work)
  , _work2(pmat._work2)
  , _work3(pmat._work3)
{
}

PrecisionOp& PrecisionOp::operator= (const PrecisionOp &pmat)
{
  _shiftOp       = pmat._shiftOp;
  _cova          = pmat._cova;
  _power         = pmat._power;
  _polynomials   = pmat._polynomials;
  _verbose       = pmat._verbose;
  _work          = pmat._work;
  _work2         = pmat._work2;
  _work3         = pmat._work3;
  return *this;
}

void PrecisionOp::_purge()
{
  for (auto& e: _polynomials)
  {
    delete e.second;
  }
  _polynomials.clear();
}

PrecisionOp::~PrecisionOp()
{
  _purge();
}

int PrecisionOp::_preparePoly(ENUM_POPTS power)
{
  // Polynomial already exists. Nothing to be done
  if (_polynomials.count(power)) return 0;

  // Prepare Polynomial for POPT_ONE
  if (_preparePrecisionPoly()) return 1;

  // Prepare polynomials for other powers than 1
  if (power != POPT_ONE)
  {
    if (_prepareChebychev(power)) return 1;
  }
  return 0;
}

int PrecisionOp::_prepareChebychev(ENUM_POPTS power)
{
  if (_cova == nullptr) return 1;
  if (_shiftOp == nullptr) return 1;

  double b = _shiftOp->getMaxEigenValue();
  APolynomial* chebMatern = new Chebychev();
  ((Chebychev*)chebMatern)->setA(0);
  ((Chebychev*)chebMatern)->setB(b);

  std::function<double(double)> f;

  if(power == POPT_LOG)
  {
    f = [this](double val){return log(_polynomials[POPT_ONE]->eval(val));};
  }
  else if (power == POPT_MINUSONE)
  {
    f = [this](double val){return pow(_polynomials[POPT_ONE]->eval(val),-1.);};
  }
  else if (power == POPT_MINUSHALF)
  {
    f = [this](double val){return pow(_polynomials[POPT_ONE]->eval(val),-0.5);};
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

  double val1 = 0.;

  VectorDouble gauss;
  VectorDouble result;
  gauss.resize(getSize());
  result.resize(getSize());

  for (int i = 0; i < nsimus; i++)
  {
    for (auto &e : gauss)
    {
      e = law_gaussian();
    }
    if (_evalPoly(POPT_LOG, gauss, result)) return TEST;

    for (int i = 0; i < getSize(); i++)
    {
      val1 += gauss[i] * result[i];
    }
  }

  val1 /= nsimus;

  double val2 = 0.;

  for (auto &e : _shiftOp->getLambda())
  {
    val2 += log(e);
  }

  val1 += 2. * val2;

  return val1;
}

int PrecisionOp::init(const ShiftOpCs* shiftop,
                      const CovAniso*  cova,
                      ENUM_POPTS       power,
                      bool             verbose)
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
    _power   = power;
    _verbose = verbose;
    _shiftOp = new ShiftOpCs(*shiftop);

    _purge();
  }

  catch(const char * str)
  {
    error = 1;
    std::cout << str << std::endl;
  }
  return error;
}

void PrecisionOp::eval(const VectorDouble& in, VectorDouble& out)
{
  const VectorDouble* inPtr = &in;

  // Pre-processing

  if (_power == POPT_ONE)
  {
    _shiftOp->prodTildeC(in, _work, POPT_HALF);
    inPtr = &_work;
  }
  else if (_power == POPT_MINUSONE)
  {
    _shiftOp->prodTildeC(in, _work, POPT_MINUSHALF);
    inPtr = &_work;
  }

  // Polynomial evaluation

  if (_evalPoly(_power,*inPtr,out))
    my_throw("Computation in 'eval' interrupted due to problem in '_evalPoly'");

  // Post-processing

  if (_power == POPT_ONE)
  {
    _shiftOp->prodTildeC(out, out, POPT_HALF);
    _shiftOp->prodLambdaOnSqrtTildeC(out, out, 2.);
  }
  else if (_power == POPT_MINUSONE)
  {
    _shiftOp->prodTildeC(out, out, POPT_MINUSHALF);
    _shiftOp->prodLambdaOnSqrtTildeC(out, out, -2.);
  }
  else if (_power == POPT_MINUSHALF)
  {
    _shiftOp->prodLambda(out, out, POPT_MINUSONE);
  }
}

int PrecisionOp::_evalPoly(ENUM_POPTS power,const VectorDouble& in, VectorDouble& out)
{
  if (_preparePoly(power)) return 1;
  _polynomials[power]->evalOp(_shiftOp->getS(),in,out);
  return 0;
}

int PrecisionOp::_preparePrecisionPoly()
{
  if (_polynomials.count(POPT_ONE)) return 0;
  VectorDouble blin;

  if (_cova == nullptr) return 1;
  int ndim = _cova->getNDim();
  double param = _cova->getParam();
  double ndims2 = ((double) ndim) / 2.;
  double alpha = param + ndims2;
  if (! isInteger(alpha,EPSILON6)) return 1;

  double correc = spde_compute_correc(ndim, param);
  int p = getClosestInteger(alpha);
  int ndimp = p + 1;
  blin.resize(ndimp);
  for (int i = 0; i < ndimp; i++)
  {
    blin[i] = ut_cnp(p, i) * correc;
    if (_verbose) message("Coefficient blin[%d] = %lf\n",i+1,blin[i]);
  }
  _polynomials[POPT_ONE] = new ClassicalPolynomial(blin);
  return 0;
}

APolynomial* PrecisionOp::getPoly(ENUM_POPTS power)
{
  if (_preparePoly(power))
    my_throw("Problem in function getPoly");
  return _polynomials[power];
}
