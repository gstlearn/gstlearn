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
#include "Matrix/MatrixRectangular.hpp"
#include "geoslib_define.h"
#include <vector>

SPDEOp::SPDEOp(const PrecisionOpMulti* const popkriging,
               const ProjMulti* const proj,
               const ASimulable* const invNoise,
               const PrecisionOpMulti* const popsimu,
                const ProjMulti* const projSimu,
               bool todelete)
  : _QKriging(popkriging)
  , _projKriging(proj)
  , _invNoise(invNoise)
  , _QSimu(popsimu == nullptr?_QKriging:popsimu)
  , _projSimu(projSimu == nullptr?_projKriging:projSimu)
  , _noiseToDelete(todelete)
  , _ndat(0)
  , _solver(this)
{
   if (_projKriging == nullptr) return;
   if (_invNoise == nullptr) return;
   if (_QKriging == nullptr) return;
   _ndat = _projKriging->getNPoint();
  _prepare(true, true);
}

SPDEOp::~SPDEOp()
{
  if (_noiseToDelete) delete _invNoise;
}

int SPDEOp::getSize() const
{
  return _QKriging->getSize();
}

void SPDEOp::_prepare(bool w1, bool w2) const
{
  if (w1) _workdat1.resize(_getNDat());
  if (w2) _workdat2.resize(_getNDat());
}

int SPDEOp::_addToDest(const constvect inv, vect outv) const
{
  return _addToDestImpl(inv, outv);
}

int SPDEOp::getSizeSimu() const
{
  return _QSimu->getSize();
}

void SPDEOp::simCond(const constvect data, vect outv) const
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
  _solver.setTolerance(1e-5);
  kriging(_workdat4,_workmesh); 
  
  //Add the kriging to the non conditional simulation
  VH::addInPlace(_workmesh,outv); 
}

VectorDouble SPDEOp::kriging(const VectorDouble& dat) const
{
  constvect datm(dat.data(), dat.size());
  VectorDouble outv(_QKriging->getSize());
  vect outvs(outv);
  int err = kriging(datm, outvs);
  if (err) return VectorDouble();
  return outv;
}

VectorDouble SPDEOp::simCond(const VectorDouble& dat) const
{
  constvect datm(dat.data(), dat.size());
  VectorDouble outv(_QSimu->getSize());
  vect outvs(outv);
  simCond(datm, outvs);
  return outv;
}

VectorDouble SPDEOp::krigingWithGuess(const VectorDouble& dat,
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

int SPDEOp::kriging(const constvect inv, vect out) const
{
  _buildRhs(inv);
  return _solve(_rhs, out);
}

int SPDEOp::krigingWithGuess(const constvect inv,
                             const constvect guess,
                             vect out) const
{
  _buildRhs(inv);
  return _solveWithGuess(_rhs, guess, out);
}

int SPDEOp::_solve(const constvect in, vect out) const
{
  _solver.solve(in, out);
  return 0;
}

int SPDEOp::_solveWithGuess(const constvect in,
                            const constvect guess,
                            vect out) const
{
  _solver.solveWithGuess(in, guess, out);
  return 0;
}

int SPDEOp::_buildRhs(const constvect inv) const
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
int SPDEOp::_addToDestImpl(const constvect inv, vect outv) const
{
  _prepare();
  vect w1s(_workdat1);
  vect w2s(_workdat2);
  _projKriging->mesh2point(inv, w1s);
  _invNoise->evalDirect(w1s, w2s);
  _projKriging->addPoint2mesh(w2s, outv);
  return _QKriging->addToDest(inv, outv);
}

void SPDEOp::evalInvCov(const constvect inv, vect result) const
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


VectorDouble SPDEOp::computeDriftCoeffs(const VectorDouble& Z,
                                        const MatrixRectangular& drifts) const
{
  int xsize = (int)(drifts.size());
  VectorDouble XtInvSigmaZ(xsize);
  MatrixSquareSymmetric XtInvSigmaX(xsize);
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
