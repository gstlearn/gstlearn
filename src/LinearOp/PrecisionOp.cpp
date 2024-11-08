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
#include "Basic/VectorNumT.hpp"
#include "Basic/AException.hpp"
#include "LinearOp/PrecisionOp.hpp"
#include "LinearOp/ShiftOpCs.hpp"
#include "Polynomials/APolynomial.hpp"
#include "Polynomials/ClassicalPolynomial.hpp"
#include "Polynomials/Chebychev.hpp"
#include "Covariances/CovAniso.hpp"
#include "Mesh/AMesh.hpp"
#include <algorithm>

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
    _work.resize (_shiftOp->getSize());
    _work2.resize(_shiftOp->getSize());
    _work3.resize(_shiftOp->getSize());
  }
}

PrecisionOp::PrecisionOp(const AMesh* mesh,
                         CovAniso* cova,
                         bool verbose)
  : _shiftOp(nullptr)
  , _cova(cova)
  , _polynomials()
  , _verbose(verbose)
  , _training(false)
  , _destroyShiftOp(true)
  , _userPoly(false)
  , _work()
  , _work2()
  , _work3()
{ 
  
  _shiftOp = new ShiftOpCs(mesh,cova,nullptr,verbose);
  if (_cova->getNVariables() == 1)
  {
    _shiftOp->normalizeLambdaBySills(mesh);
  }
  _work.resize(_shiftOp->getSize());
  _work2.resize(_shiftOp->getSize());
  _work3.resize(_shiftOp->getSize());
}

PrecisionOp::PrecisionOp(const PrecisionOp& pmat)
  : _shiftOp(nullptr)
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
  if (_destroyShiftOp)
    _shiftOp = new ShiftOpCs(*pmat._shiftOp);
  else
    _shiftOp = pmat._shiftOp;
}

PrecisionOp& PrecisionOp::operator= (const PrecisionOp &pmat)
{
  if (this != &pmat)
  {
    _cova = pmat._cova;
    _polynomials = pmat._polynomials;
    _verbose = pmat._verbose;
    _training = pmat._training;
    _destroyShiftOp = pmat._destroyShiftOp;
    _userPoly = pmat._userPoly;
    _work = pmat._work;
    _work2 = pmat._work2;
    _work3 = pmat._work3;

    if (_destroyShiftOp)
      _shiftOp = new ShiftOpCs(*pmat._shiftOp);
    else
      _shiftOp = pmat._shiftOp;
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

PrecisionOp* PrecisionOp::create(const AMesh* mesh,
                                 CovAniso* cova,
                                 bool verbose)
{
  return new PrecisionOp(mesh, cova, verbose);
}

int PrecisionOp::_addToDest(const constvect inv, vect outv) const
{
    _addEvalPower(inv, outv, EPowerPT::ONE);
    return 0;
}

int PrecisionOp::_addSimulateToDest(const constvect whitenoise, vect outv) const
{
    _addEvalPower(whitenoise, outv, EPowerPT::MINUSHALF);
    return 0;
}

int PrecisionOp::_preparePoly(const EPowerPT& power,bool force) const
{
  // Polynomial already exists. Nothing to be done
  if (_polynomials.contains(power)  && !force) return 0;

  // Prepare Polynomial for EPowerPT::ONE
  if (_preparePrecisionPoly() != 0 && !force) return 1;

  // Prepare polynomials for other powers than 1
  if (power != EPowerPT::ONE)
  {
    if (_prepareChebychev(power) != 0) return 1;
  }
  return 0;
}

VectorDouble PrecisionOp::getPolyCoeffs(const EPowerPT& power)
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

int PrecisionOp::_prepareChebychev(const EPowerPT& power) const
{
  if (_cova == nullptr && _polynomials.contains(EPowerPT::ONE)) return 1;
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
 * @param nbsimu Number of simulations
 * @return The computed value or TEST if problem
 */
double PrecisionOp::getLogDeterminant(int nbsimu)
{
  VectorDouble gauss;
  VectorDouble result;
  gauss.resize(getSize());
  result.resize(getSize());

  double val1 = 0.;
  for (int isimu = 0; isimu < nbsimu; isimu++)
  {
    VH::simulateGaussianInPlace(gauss);
    vect results(result);
    if (_evalPoly(EPowerPT::LOG, gauss, results) != 0) return TEST;

    for (int i = 0; i < getSize(); i++)
    {
      val1 += gauss[i] * result[i];
    }
  }
  val1 /= nbsimu;

  double val2 = 0.;
  for (const auto &e : _shiftOp->getLambdas())
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
  catch(const std::string& str)
  {
    // TODO : Check if std::exception can be used
    error = 1;
    messerr("%s", str.c_str());
  }

  return error;
}

/**
 * Evaluate with power = ONE
 * @param inm Input array
 * @param outm Output array
 * @param power Power of the operation
 */
/* void PrecisionOp::evalDirect(const VectorDouble &vecin, VectorDouble &vecout)
{
  evalPower(vecin, vecout, EPowerPT::ONE);
}
 */
void PrecisionOp::evalPower(const constvect inm,
                            vect outm,
                            const EPowerPT& power)
{
  std::fill(outm.begin(),outm.end(),0.);
  for (int i = 0; i < (int)outm.size(); i++)
  _addEvalPower(inm, outm, power);
}

void PrecisionOp::_addEvalPower(const constvect inv,
                                vect outv,
                                const EPowerPT& power) const
{
  const constvect* inPtr = &inv;
  if (_work.size() == 0) _work.resize(getSize());
  vect worksp(_work);
  
  // Pre-processing

  if (power == EPowerPT::ONE || power == EPowerPT::MINUSONE)
  {
    _shiftOp->prodLambda(inv, worksp, power);
    inPtr = (constvect*)&worksp;
  }

  // Polynomial evaluation

  if (_evalPoly(power, *inPtr, outv) != 0)
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
                           const constvect inv,
                           vect outv) const
{
  constvect invs(inv);
  if (_preparePoly(power) != 0) return 1;
  if (getTraining())
  {
    int degree = _polynomials[power]->getDegree();

    if(_workPoly.empty())
    {
      _workPoly = std::vector<std::vector<double>>(degree);
      for(auto &e: _workPoly)
      {
        e.resize(inv.size());
      }
    }

    if (_work5.size() == 0) _work5.resize(getSize());
    
    ((ClassicalPolynomial*)_polynomials[power])->evalOpTraining(_shiftOp->getS(),
                                                                invs,_workPoly,
                                                                _work5);

    for (int i = 0; i < (int)inv.size(); i++)
    {
      outv[i] = _workPoly[0][i];
    }
  }
  else
  {
    vect outvs(outv);
    _polynomials[power]->evalOp(_shiftOp->getS(), invs, outvs);
  }
  return 0;
}

void PrecisionOp::evalInverse(const constvect vecin,
                              std::vector<double>& vecout)
{
  if (_work.size() != vecin.size()) _work.resize(vecin.size());
  vect vecouts(vecout);
  _shiftOp->prodLambda(vecin,vecouts,EPowerPT::MINUSONE);
  vect works(_work);
  _evalPoly(EPowerPT::MINUSONE, vecout, works);
  _shiftOp->prodLambda(works, vecouts, EPowerPT::MINUSONE);
}

VectorDouble PrecisionOp::evalCov(int imesh)
{

  int n = getSize();
  VectorDouble result(n);
  std::vector<double> ei(n);
  vect eis(ei);
  vect results(result);

  std::fill(ei.begin(),ei.end(),0.);
  ei[imesh] = 1.;
  _shiftOp->prodLambda(eis,results,EPowerPT::MINUSONE);
  _evalPoly(EPowerPT::MINUSONE,result,eis);
  _shiftOp->prodLambda(eis, results, EPowerPT::MINUSONE);
  return result;
}

VectorVectorDouble PrecisionOp::simulate(int nbsimu)
{
  int n = getSize();
  VectorVectorDouble vectv(nbsimu);
  VectorDouble whitenoise(n);

  for(auto &e : vectv)
  {
    e.resize(n);
    VH::simulateGaussianInPlace(whitenoise);
    vect es(e);
    _evalPoly(EPowerPT::MINUSHALF,whitenoise,es);
    _shiftOp->prodLambda(e, e, EPowerPT::MINUSONE);
  }
  return vectv;
}

VectorDouble PrecisionOp::simulateOne()
{
  int n = getSize();
  VectorDouble vectv(n);
  VectorDouble whitenoise(n);
  VH::simulateGaussianInPlace(whitenoise);
  vect vects(vectv);
  _evalPoly(EPowerPT::MINUSHALF,whitenoise,vects);
  _shiftOp->prodLambda(vectv, vectv, EPowerPT::MINUSONE);
  return vectv;
}

int PrecisionOp::_preparePrecisionPoly() const
{
  if (_cova == nullptr) return 1;
  if (!_cova->hasMarkovCoeffs()) return 1;

  _polynomials[EPowerPT::ONE] = new ClassicalPolynomial(_cova->getMarkovCoeffs());

  return 0;
}

/* Evaluation of the polynomial of the precision over the interval [0, lambda_max(ShiftOp)] */
/* discretized over ndiscr points */
/* and return the min and the max */

std::pair<double,double> PrecisionOp::getRangeEigenVal(int ndiscr)
{
  std::pair<double,double> rangeVals;
  double sill = _cova->getSill(0,0); //TODO handle non constant sill
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

  rangeVals.first  = rangeVals.first  / sill;
  rangeVals.second = rangeVals.second / sill;

  return rangeVals;
}

APolynomial* PrecisionOp::getPoly(const EPowerPT& power)
{
  if (_preparePoly(power) != 0)
    my_throw("Problem in function getPoly");
  return _polynomials[power];
}

VectorDouble PrecisionOp::getCoeffs()
{
  VectorDouble coeffs = getPoly(EPowerPT::ONE)->getCoeffs();
  return coeffs;
}

VectorDouble PrecisionOp::extractDiag() const
{
  int size = getSize();
  VectorDouble vec(size, 0.);
  const EPowerPT& power = EPowerPT::ONE;

  VectorDouble lambdas = _shiftOp->getLambdas();

  if (_preparePoly(power) != 0) return vec;

  for (int i = 0; i < size; i++)
  {
    double lambda = lambdas[i];
    vec[i] = _polynomials[power]->evalOpByRank(_shiftOp->getS(), i) * lambda * lambda;
  }
  return vec;
}
