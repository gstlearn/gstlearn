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

SPDEOp::SPDEOp(const PrecisionOpMulti* const pop,
               const ProjMulti* const proj,
               const ASimulable* const invNoise,
               bool todelete)
  : _Q(pop)
  , _Proj(proj)
  , _invNoise(invNoise)
  , _noiseToDelete(todelete)
  , _ndat(0)
  , _solver(this)
{
   if (_Proj == nullptr) return;
   if (_invNoise == nullptr) return;
   if (_Q == nullptr) return;
   _ndat = _Proj->getPointNumber();
  _prepare(true, true);
}

SPDEOp::~SPDEOp()
{
  if (_noiseToDelete) delete _invNoise;
}

int SPDEOp::getSize() const
{
  return _Q->getSize();
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

int SPDEOp::_addSimulateToDest(const constvect whitenoise, vect outv) const
{
  DECLARE_UNUSED(whitenoise);
  DECLARE_UNUSED(outv);
  return 1;
}

VectorDouble SPDEOp::kriging(const VectorDouble& dat) const
{
  constvect datm(dat.data(), dat.size());
  VectorDouble outv(_Q->getSize());
  vect outvs(outv);
  int err = kriging(datm, outvs);
  if (err) return VectorDouble();
  return outv;
}

VectorDouble SPDEOp::krigingWithGuess(const VectorDouble& dat,
                                      const VectorDouble& guess) const
{
  constvect datm(dat.data(), dat.size());
  constvect guessm(guess.data(), guess.size());

  VectorDouble outv(_Q->getSize());
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
  _rhs.resize(_Q->getSize());
  vect w1(_workdat1);
  _invNoise->evalDirect(inv, w1);
  _Proj->point2mesh(_workdat1, _rhs);
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
  _Proj->mesh2point(inv, w1s);
  _invNoise->evalDirect(w1s, w2s);
  _Proj->addPoint2mesh(w2s, outv);
  return _Q->addToDest(inv, outv);
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
  _Proj->point2mesh(result,rhss);
  _solve(rhss,wms);
  _Proj->mesh2point(wms,w2s);
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
      // Eigen::Map<const Eigen::VectorXd> xmj(drifts[j].data(),drifts[j].size());
      // XtInvSigmaX.setValue(i,j,  xmj.adjoint() * wd1);
    }
  }

  XtInvSigmaX.solve(XtInvSigmaZ,result);

  return result;
}
