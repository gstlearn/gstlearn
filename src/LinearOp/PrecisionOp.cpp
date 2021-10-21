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
#include "geoslib_old_f.h"

PrecisionOp::PrecisionOp(ShiftOpCs* shiftop,
                         const CovAniso* cova,
                         const EPowerPT& power,
                         bool verbose)
  : _shiftOp(shiftop)
  , _cova(cova)
  , _power(power)
  , _polynomials()
  , _verbose(verbose)
  , _training(false)
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
  , _training(false)
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
  _training      = pmat._training;
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

int PrecisionOp::_preparePoly(const EPowerPT& power)
{
  // Polynomial already exists. Nothing to be done
  if (_polynomials.count(power)) return 0;

  // Prepare Polynomial for EPowerPT::ONE
  if (_preparePrecisionPoly()) return 1;

  // Prepare polynomials for other powers than 1
  if (power != EPowerPT::ONE)
  {
    if (_prepareChebychev(power)) return 1;
  }
  return 0;
}

int PrecisionOp::_prepareChebychev(const EPowerPT& power)
{
  if (_cova == nullptr) return 1;
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
    if (_evalPoly(EPowerPT::LOG, gauss, result)) return TEST;

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
                      const EPowerPT&  power,
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
  if (_work.empty()) _work.resize(getSize());

  // Pre-processing

  if (_power == EPowerPT::ONE)
  {
    _shiftOp->prodTildeC(in, _work, EPowerPT::HALF);
    inPtr = &_work;
  }
  else if (_power == EPowerPT::MINUSONE)
  {
    _shiftOp->prodTildeC(in, _work, EPowerPT::MINUSHALF);
    inPtr = &_work;
  }

  // Polynomial evaluation

  if (_evalPoly(_power,*inPtr,out))
    my_throw("Computation in 'eval' interrupted due to problem in '_evalPoly'");

  // Post-processing

  if (_power == EPowerPT::ONE)
  {
    _shiftOp->prodTildeC(out, out, EPowerPT::HALF);
    _shiftOp->prodLambdaOnSqrtTildeC(out, out, 2.);
  }
  else if (_power == EPowerPT::MINUSONE)
  {
    _shiftOp->prodTildeC(out, out, EPowerPT::MINUSHALF);
    _shiftOp->prodLambdaOnSqrtTildeC(out, out, -2.);
  }
  else if (_power == EPowerPT::MINUSHALF)
  {
    _shiftOp->prodLambda(out, out, EPowerPT::MINUSONE);
  }
}

int PrecisionOp::_evalPoly(const EPowerPT& power,
                           const VectorDouble& in,
                           VectorDouble& out)
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
        e = VectorDouble(in.size());
      }
    }


    ((ClassicalPolynomial*)_polynomials[power])->evalOpTraining(_shiftOp->getS(),in,_workPoly,_work4);
    for(int i=0;i<(int)in.size();i++)
    {
       out[i] = _workPoly[0][i];
    }
   }
    else
    {
      _polynomials[power]->evalOp(_shiftOp->getS(),in,out);
    }
  return 0;
}

int PrecisionOp::_preparePrecisionPoly()
{
  if (_polynomials.count(EPowerPT::ONE)) return 0;
  VectorDouble blin;

  if (_cova == nullptr) return 1;
  int ndim = _cova->getNDim();
  double param = _cova->getParam();
  double ndims2 = ((double) ndim) / 2.;
  double alpha = param + ndims2;
  if (! isInteger(alpha,EPSILON3)) return 1;

  double correc = spde_compute_correc(ndim, param);
  int p = getClosestInteger(alpha);
  int ndimp = p + 1;
  blin.resize(ndimp);
  for (int i = 0; i < ndimp; i++)
  {
    blin[i] = ut_cnp(p, i) * correc;
    if (_verbose) message("Coefficient blin[%d] = %lf\n",i+1,blin[i]);
  }
  _polynomials[EPowerPT::ONE] = new ClassicalPolynomial(blin);
  return 0;
}

APolynomial* PrecisionOp::getPoly(const EPowerPT& power)
{
  if (_preparePoly(power))
    my_throw("Problem in function getPoly");
  return _polynomials[power];
}
