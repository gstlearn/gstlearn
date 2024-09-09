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
#include <Eigen/src/Core/Matrix.h>
#include "Matrix/VectorEigen.hpp"
#include "geoslib_define.h"

SPDEOp::SPDEOp(const PrecisionOpMulti* const pop,
               const ProjMulti* const proj,
               const ASimulable* const invNoise,
               bool todelete)
  : _Q(pop)
  , _Proj(proj)
  , _invNoise(invNoise)
  , _noiseToDelete(todelete)
  , _solver(this)
{
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
  if (w1) _workdat1.resize(_Proj->getPointNumber());
  if (w2) _workdat2.resize(_Proj->getPointNumber());
}

int SPDEOp::_addToDest(const Eigen::VectorXd& inv, Eigen::VectorXd& outv) const
{
  return _addToDestImpl(inv, outv);
}

int SPDEOp::_addSimulateToDest(const Eigen::VectorXd& whitenoise,
                               Eigen::VectorXd& outv) const
{
  DECLARE_UNUSED(whitenoise);
  DECLARE_UNUSED(outv);
  return 1;
}

VectorDouble SPDEOp::kriging(const VectorDouble& dat) const
{
  Eigen::Map<const Eigen::VectorXd> datm(dat.data(), dat.size());
  Eigen::VectorXd outv(_Q->getSize());
  int err = kriging(datm, outv);
  if (err) return VectorDouble();
  return VectorEigen::copyIntoVD(outv);
}

VectorDouble SPDEOp::krigingWithGuess(const VectorDouble& dat,
                                      const VectorDouble& guess) const
{
  Eigen::Map<const Eigen::VectorXd> datm(dat.data(), dat.size());
  Eigen::Map<const Eigen::VectorXd> guessm(guess.data(), guess.size());

  Eigen::VectorXd outv(_Q->getSize());
  int err = krigingWithGuess(datm, guessm, outv);
  if (err) return VectorDouble();
  return VectorEigen::copyIntoVD(outv);
}

int SPDEOp::kriging(const Eigen::VectorXd& inv, Eigen::VectorXd& out) const
{
  out.resize(_Q->getSize());
  _buildRhs(inv);
  return _solve(_rhs, out);
}

int SPDEOp::krigingWithGuess(const Eigen::VectorXd& inv,
                             const Eigen::VectorXd& guess,
                             Eigen::VectorXd& out) const
{
  out.resize(_Q->getSize());
  _buildRhs(inv);
  return _solveWithGuess(_rhs, guess, out);
}

int SPDEOp::_solve(const Eigen::VectorXd& in, Eigen::VectorXd& out) const
{
  _solver.solve(in, out);
  return 0;
}

int SPDEOp::_solveWithGuess(const Eigen::VectorXd& in,
                            const Eigen::VectorXd& guess,
                            Eigen::VectorXd& out) const
{
  _solver.solveWithGuess(in, guess, out);
  return 0;
}

int SPDEOp::_buildRhs(const Eigen::VectorXd& inv) const
{
  _rhs.resize(_Q->getSize());
  _invNoise->evalDirect(inv, _workdat1);
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
int SPDEOp::_addToDestImpl(const Eigen::VectorXd& inv,
                           Eigen::VectorXd& outv) const
{
  _prepare();
  _Proj->mesh2point(inv, _workdat1);
  _invNoise->evalDirect(_workdat1, _workdat2);
  _Proj->addPoint2mesh(_workdat2, outv);
  return _Q->addToDest(inv, outv);
}
