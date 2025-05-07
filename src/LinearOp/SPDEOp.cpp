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
#include "LinearOp/SPDEOp.hpp"
#include "Basic/VectorNumT.hpp"
#include "LinearOp/ProjMulti.hpp"
#include "LinearOp/PrecisionOpMulti.hpp"
#include "Matrix/MatrixDense.hpp"
#include "geoslib_define.h"

ASPDEOp::ASPDEOp(const PrecisionOpMulti* const popKriging,
                 const ProjMulti* const projKriging,
                 const ASimulable* const invNoise,
                 const PrecisionOpMulti* const popSimu,
                 const ProjMulti* const projSimu,
                 bool todelete)
  : _QKriging(popKriging)
  , _projKriging(projKriging)
  , _invNoise(invNoise)
  , _QSimu(popSimu == nullptr ? popKriging : popSimu)
  , _projSimu(projSimu == nullptr ? projKriging : projSimu)
  , _solver(nullptr)
  , _noiseToDelete(todelete)
  , _ndat(0)
{
  if (_projKriging == nullptr) return;
  if (_invNoise == nullptr) return;
  if (_QKriging == nullptr) return;
  _ndat = _projKriging->getNPoint();
  _prepare(true, true);
}

ASPDEOp::~ASPDEOp()
{
  if (_noiseToDelete) delete _invNoise;
  delete _solver;
}

int ASPDEOp::getSize() const
{
  return _QKriging->getSize();
}

void ASPDEOp::_prepare(bool w1, bool w2) const
{
  if (w1) _workdat1.resize(_getNDat());
  if (w2) _workdat2.resize(_getNDat());
}

int ASPDEOp::_addToDest(const constvect inv, vect outv) const
{
  return _addToDestImpl(inv, outv);
}

int ASPDEOp::getSizeSimu() const
{
  return _QSimu->getSize();
}

void ASPDEOp::simCond(const constvect data, vect outv) const
{
  // Resize if necessary
  _workdat3.resize(_getNDat());
  _workdat4.resize(_getNDat());
  _workmesh.resize(getSizeSimu());
  _workNoiseMesh.resize(getSizeSimu());
  _workNoiseData.resize(_getNDat());

  //Non conditional simulation on mesh
  VH::simulateGaussianInPlace(_workNoiseMesh);
  _QSimu->evalSimulate(_workNoiseMesh, outv); 
  
  //Simulation at data locations (projection + noise)
  
  _projSimu->mesh2point(outv, _workdat3); //Projection on data locations
  VH::simulateGaussianInPlace(_workNoiseData);
  _invNoise->addSimulateToDest(_workNoiseData, _workdat3); //Add noise
  
  //compute residual _workdat4 = data - outv
  VH::subtractInPlace(_workdat3, data, _workdat4);

  //Co-Kriging of the residual on the mesh
  _solver->setTolerance(1e-5);
  kriging(_workdat4,_workmesh); 
  
  //Add the kriging to the non conditional simulation
  VH::addInPlace(_workmesh,outv); 
}

VectorDouble ASPDEOp::kriging(const VectorDouble& dat) const
{
  constvect datm(dat.data(), dat.size());
  VectorDouble outv(_QKriging->getSize());
  vect outvs(outv);
  int err = kriging(datm, outvs);
  if (err) return VectorDouble();
  return outv;
}

VectorDouble ASPDEOp::simCond(const VectorDouble& dat) const
{
  constvect datm(dat.data(), dat.size());
  VectorDouble outv(_QSimu->getSize());
  vect outvs(outv);
  simCond(datm, outvs);
  return outv;
}

VectorDouble ASPDEOp::krigingWithGuess(const VectorDouble& dat,
                                       const VectorDouble& guess) const
{
  constvect datm(dat.data(), dat.size());
  constvect guessm(guess.data(), guess.size());

  VectorDouble outv(_QKriging->getSize());
  vect outvs(outv);
  int err = krigingWithGuess(datm, guessm, outvs);
  if (err) return VectorDouble();
  return outv;
}

int ASPDEOp::kriging(const constvect inv, vect out) const
{
  _buildRhs(inv);
  return _solve(_rhs, out);
}

int ASPDEOp::krigingWithGuess(const constvect inv,
                              const constvect guess,
                              vect out) const
{
  _buildRhs(inv);
  return _solveWithGuess(_rhs, guess, out);
}

int ASPDEOp::_solve(const constvect in, vect out) const
{
  _solver->solve(in, out);
  return 0;
}

int ASPDEOp::_solveWithGuess(const constvect in,
                             const constvect guess,
                             vect out) const
{
  _solver->solveWithGuess(in, guess, out);
  return 0;
}

int ASPDEOp::_buildRhs(const constvect inv) const
{
  _rhs.resize(_QKriging->getSize());
  vect w1(_workdat1);
  _invNoise->evalDirect(inv, w1);
  _projKriging->point2mesh(_workdat1, _rhs);
  return 0;
}

/*****************************************************************************/
/*!
**  Evaluate the product (by the SPDEOp) :
**  'outv' = (_Q + _Proj' * _invNoise * Proj) * 'inv'
**
** \param[in]  inv     Array of input values
**
** \param[out] outv    Array of output values
**
*****************************************************************************/
int ASPDEOp::_addToDestImpl(const constvect inv, vect outv) const
{
  _prepare();
  vect w1s(_workdat1);
  vect w2s(_workdat2);
  _projKriging->mesh2point(inv, w1s);
  _invNoise->evalDirect(w1s, w2s);
  _projKriging->addPoint2mesh(w2s, outv);
  return _QKriging->addToDest(inv, outv);
}

void ASPDEOp::evalInvCov(const constvect inv, vect result) const
{
  // InvNoise - InvNoise * Proj' * (Q + Proj * InvNoise * Proj')^-1 * Proj * InvNoise
  
  _rhs.resize(getSize());
  _workmesh.resize(getSize());
  _workdat2.resize(_getNDat());
  vect rhss(_rhs);
  vect wms(_workmesh);
  vect w2s(_workdat2);

  _invNoise->evalDirect(inv,result);
  _projKriging->point2mesh(result,rhss);
  _solve(rhss,wms);
  _projKriging->mesh2point(wms,w2s);
  //VectorHelper::multiplyConstant(w2s,-1);
  _invNoise->addToDest(w2s,result);
}

VectorDouble ASPDEOp::computeDriftCoeffs(const VectorDouble& Z,
                                         const MatrixDense& drifts) const
{
  int xsize = (int)(drifts.getNCols());
  VectorDouble XtInvSigmaZ(xsize);
  MatrixSymmetric XtInvSigmaX(xsize);
  VectorDouble result(xsize);
  _workdat1.resize(_getNDat());
  vect w1s(_workdat1);
  for(int i = 0; i< xsize; i++)
  {
    auto xm = drifts.getColumnPtr(i);
    evalInvCov(xm,w1s);

    constvect ym(Z.data(),Z.size());
    constvect wd1(_workdat1.data(),_workdat1.size());
    XtInvSigmaZ[i] = VH::innerProduct(ym,wd1);
    for(int j = i; j < xsize;j++)
    {
      constvect xmj = drifts.getViewOnColumn(j);
      double prod = VH::innerProduct(xmj,w1s);
      XtInvSigmaX.setValue(i,j,prod);
    }
  }

  XtInvSigmaX.solve(XtInvSigmaZ,result);

  return result;
}
