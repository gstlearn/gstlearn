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
#include "Basic/Law.hpp"
#include "Basic/VectorNumT.hpp"
#include "LinearOp/ASimulable.hpp"
#include "LinearOp/PrecisionOpMulti.hpp"
#include "LinearOp/ProjMulti.hpp"
#include "Matrix/MatrixDense.hpp"
#include "Polynomials/Chebychev.hpp"
#include "geoslib_define.h"

namespace gstlrn
{
ASPDEOp::ASPDEOp(const PrecisionOpMulti* const popKriging,
                 const ProjMulti* const projInKriging,
                 const ASimulable* invNoise,
                 const PrecisionOpMulti* const popSimu,
                 const ProjMulti* const projInSimu)
  : _QKriging(popKriging)
  , _projInKriging(projInKriging)
  , _invNoise(invNoise)
  , _QSimu(popSimu == nullptr ? popKriging : popSimu)
  , _projInSimu(projInSimu == nullptr ? projInKriging : projInSimu)
  , _solver(nullptr)
  , _verbose(false)
  , _ndat(0)
{
  if (_projInKriging == nullptr) return;
  if (_invNoise == nullptr) return;
  if (_QKriging == nullptr) return;
  _ndat = _projInKriging->getNPoint();
  _prepare(true, true);
}

ASPDEOp::~ASPDEOp()
{
  delete _solver;
}

Id ASPDEOp::getSize() const
{
  return _QKriging->getSize();
}

Id ASPDEOp::getSizeSimu() const
{
  return _QSimu->getSize();
}

void ASPDEOp::_prepare(bool w1, bool w2) const
{
  if (w1) _workdat1.resize(_getNDat());
  if (w2) _workdat2.resize(_getNDat());
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
// Id ASPDEOp::_addToDest(const constvect inv, vect outv) const
// {
//   _prepare();

//   vect w1s(_workdat1);
//   vect w2s(_workdat2);
//   _projKriging->mesh2point(inv, w1s);
//   _invNoise->evalDirect(w1s, w2s);
//   _projKriging->addPoint2mesh(w2s, outv);
//   return _QKriging->addToDest(inv, outv);
// }

Id ASPDEOp::_addToDest(const constvect inv, vect outv) const
{
  _prepare();

  Id status = _QKriging->addToDest(inv, outv); // TODO: find why outv is set to zero in multistructure case
  if (status) return status;
  vect w1s(_workdat1);
  vect w2s(_workdat2);
  _projInKriging->mesh2point(inv, w1s);
  _invNoise->evalDirect(w1s, w2s);
  _projInKriging->addPoint2mesh(w2s, outv);

  return status;
}

void ASPDEOp::_simCond(const constvect data, vect outvK, vect outvS) const
{
  // Resize if necessary
  _workdat3.resize(_getNDat());
  _workdat4.resize(_getNDat());
  _workmesh.resize(getSizeSimu());
  _workNoiseMesh.resize(getSizeSimu());
  _workNoiseData.resize(_getNDat());

  // Non conditional simulation on Simulation mesh
  VH::simulateGaussianInPlace(_workNoiseMesh);
  _QSimu->evalSimulate(_workNoiseMesh, outvS);

  // Simulation at data locations (projection + noise)
  _projInSimu->mesh2point(outvS, _workdat3); // Projection on data locations
  VH::simulateGaussianInPlace(_workNoiseData);
  _invNoise->addSimulateToDest(_workNoiseData, _workdat3); // Add noise

  // compute residual _workdat4 = data - outv
  VH::subtractInPlace(_workdat3, data, _workdat4);

  // Co-Kriging of the residual on the Kriging mesh
  _solver->setTolerance(1e-5);
  _kriging(_workdat4, outvK);
}

void ASPDEOp::_simNonCond(vect outv) const
{
  // Resize if necessary
  _workmesh.resize(getSizeSimu());
  _workNoiseMesh.resize(getSizeSimu());

  // Non conditional simulation on mesh
  VH::simulateGaussianInPlace(_workNoiseMesh);
  _QSimu->evalSimulate(_workNoiseMesh, outv);
}

VectorDouble ASPDEOp::kriging(const VectorDouble& dat, const ProjMulti* proj) const
{
  constvect datm(dat.data(), dat.size());
  VectorDouble outMeshK(getSize());
  vect outvs(outMeshK);
  Id err = _kriging(datm, outvs);
  if (err) return VectorDouble();

  // Project the result on the output mesh (optional)
  if (proj == nullptr) return outMeshK;
  VectorDouble result(proj->getNPoint());
  proj->mesh2point(outvs, result);
  return result;
}

Id ASPDEOp::centerDataByDriftMat(VectorDouble& Z,
                                 const MatrixDense& driftMat,
                                 const VectorDouble& driftCoeffs)
{
  auto nrows = driftMat.getNRows();
  auto ncols = driftMat.getNCols();
  if (nrows != static_cast<Id>(Z.size()))
  {
    messerr("Error in number of Rows of drift matrix (%d) and size of data vector (%d)",
            nrows, Z.size());
    return 1;
  }
  if (ncols != static_cast<Id>(driftCoeffs.size()))
  {
    messerr("Error in number of Columns of drift matrix (%d) and size of drift coefficients (%d)",
            ncols, driftCoeffs.size());
    return 1;
  }

  // Center the data set
  for (Id i = 0; i < nrows; i++)
  {
    double sum = 0.;
    for (Id j = 0; j < ncols; j++)
    {
      sum += driftCoeffs[j] * driftMat.getValue(i, j);
    }
    Z[i] -= sum;
  }
  return 0;
}

Id ASPDEOp::centerDataByMeanVec(VectorDouble& Z,
                                const VectorDouble& meanVec)
{
  if (static_cast<Id>(Z.size()) != static_cast<Id>(meanVec.size()))
  {
    messerr("Error in size of data vector (%d) and size of mean vector (%d)",
            Z.size(), meanVec.size());
    return 1;
  }

  // Center the data set
  for (Id i = 0, nrows = static_cast<Id>(Z.size()); i < nrows; i++)
    Z[i] -= meanVec[i];
  return 0;
}

VectorDouble ASPDEOp::simNonCond(const ProjMulti* proj) const
{
  VectorDouble outMeshS(getSizeSimu());
  vect outMeshSv(outMeshS);
  _simNonCond(outMeshSv);

  // Project the result on the output mesh (optional)
  if (proj == nullptr) return outMeshS;
  VectorDouble result(proj->getNPoint());
  proj->mesh2point(outMeshSv, result);
  return result;
}

VectorDouble ASPDEOp::simCond(const VectorDouble& dat,
                              const ProjMulti* projK,
                              const ProjMulti* projS) const
{
  constvect datv(dat.data(), dat.size());
  VectorDouble outMeshK(getSize());
  vect outMeshKv(outMeshK);
  VectorDouble outMeshS(getSizeSimu());
  vect outMeshSv(outMeshS);

  // Perform the conditional simulation
  _simCond(datv, outMeshKv, outMeshSv);

  // Project the result on the output mesh (optional)
  if (projK == nullptr || projS == nullptr)
  {
    VH::addInPlace(outMeshSv, outMeshKv);
    return outMeshK;
  }
  VectorDouble result(projS->getNPoint());
  projK->mesh2point(outMeshKv, result);
  projS->addMesh2point(outMeshSv, result);
  return result;
}

/**
 * @brief Computing Standard deviation of the estimation error using MonteCarlo
 * on conditional simulations
 *
 * @param dat Vector of Data
 * @param nMC  Number of Monte-Carlo simulations
 * @param seed Random seed for the Monte-Carlo simulations
 * @param projK Projection Matrix used for Kriging
 * @param projS Projection matrix used for Simulations
 * @return VectorDouble
 */
VectorDouble ASPDEOp::stdev(const VectorDouble& dat,
                            Id nMC,
                            Id seed,
                            const ProjMulti* projK,
                            const ProjMulti* projS) const
{
  auto memo = law_get_random_seed();
  law_set_random_seed(seed);

  // Standard Deviation using Monte-Carlo simulations
  Id nout = projS->getNPoint();
  VectorDouble temp_mean(nout, 0.);
  VectorDouble temp_mean2(nout, 0.);

  for (Id iMC = 0; iMC < nMC; iMC++)
  {
    VectorDouble temp = simCond(dat, projK, projS);
    VH::addInPlace(temp_mean, temp);
    VH::addSquareInPlace(temp_mean2, temp);
  }
  VH::mean1AndMean2ToStdev(temp_mean, temp_mean2, temp_mean, nMC);

  law_set_random_seed(memo);

  return temp_mean;
}

VectorDouble ASPDEOp::krigingWithGuess(const VectorDouble& dat,
                                       const VectorDouble& guess) const
{
  constvect datv(dat.data(), dat.size());
  constvect guessv(guess.data(), guess.size());

  VectorDouble outv(getSize());
  vect outvs(outv);
  Id err = krigingWithGuess(datv, guessv, outvs);
  if (err) return VectorDouble();
  return outv;
}

Id ASPDEOp::_kriging(const constvect inv, vect out) const
{
  _buildRhs(inv);
  Id status = _solve(_rhs, out);
  return status;
}

Id ASPDEOp::krigingWithGuess(const constvect inv,
                             const constvect guess,
                             vect out) const
{
  _buildRhs(inv);
  return _solveWithGuess(_rhs, guess, out);
}

Id ASPDEOp::_solve(const constvect in, vect out) const
{
  _solver->solve(in, out);
  return 0;
}

Id ASPDEOp::_solveWithGuess(const constvect in,
                            const constvect guess,
                            vect out) const
{
  _solver->solveWithGuess(in, guess, out);
  return 0;
}

Id ASPDEOp::_buildRhs(const constvect inv) const
{
  _rhs.resize(getSize());
  vect w1(_workdat1);
  _invNoise->evalDirect(inv, w1);
  _projInKriging->point2mesh(_workdat1, _rhs);
  return 0;
}

void ASPDEOp::evalInvCov(const constvect inv, vect result) const
{
  // InvNoise - InvNoise * Proj' * (Q + Proj * InvNoise * Proj')^-1 * Proj * InvNoise

  _rhs.resize(getSize());
  _workmesh.resize(getSize());
  _workdat2.resize(_getNDat());
  _workdat3.resize(_getNDat());

  _invNoise->evalDirect(inv, result);
  _projInKriging->point2mesh(result, _rhs);
  _solve(_rhs, _workmesh);
  _projInKriging->mesh2point(_workmesh, _workdat2);
  _invNoise->evalDirect(_workdat2, _workdat3);
  VectorHelper::subtractInPlace(_workdat3, result, result);
}

VectorDouble ASPDEOp::computeDriftCoeffs(const VectorDouble& Z,
                                         const MatrixDense& driftMat,
                                         bool verbose) const
{
  Id xsize = (driftMat.getNCols());
  VectorDouble XtInvSigmaZ(xsize);
  MatrixSymmetric XtInvSigmaX(xsize);
  VectorDouble result(xsize);

  _workdat1.resize(_getNDat());
  vect w1s(_workdat1);
  for (Id i = 0; i < xsize; i++)
  {
    auto xm = driftMat.getColumnPtr(i);
    evalInvCov(xm, w1s);

    constvect ym(Z.data(), Z.size());
    constvect wd1(_workdat1.data(), _workdat1.size());
    XtInvSigmaZ[i] = VH::innerProduct(ym, wd1);

    for (Id j = i; j < xsize; j++)
    {
      constvect xmj = driftMat.getViewOnColumn(j);
      double prod   = VH::innerProduct(xmj, w1s);
      XtInvSigmaX.setValue(i, j, prod);
    }
  }

  XtInvSigmaX.solve(XtInvSigmaZ, result);

  // Optional printout
  if (verbose)
    VH::dump("Drift coefficients", result);

  return result;
}

std::pair<double, double> ASPDEOp::_computeRangeEigenVal() const
{
  std::pair<double, double> result = _QKriging->rangeEigenValQ();
  // result.second += getMaxEigenValProj();
  return result;
}

void ASPDEOp::_preparePoly(Chebychev& logPoly) const
{
  std::pair<double, double> ranges = _computeRangeEigenVal();
  double a                         = ranges.first;
  double b                         = ranges.second;
  logPoly.setA(a);
  logPoly.setB(b);
  logPoly.setNcMax(1500);
  logPoly.fit([](double val)
              { return log(val); }, a, b, 2 * EPSILON4 / (a + b));
}

double ASPDEOp::computeLogDetOp(Id nbsimu) const
{
  Chebychev logPoly;
  _preparePoly(logPoly);
  _workNoiseMesh.resize(getSize());
  _workmesh.resize(getSize());
  double val = 0.;
  for (Id i = 0; i < nbsimu; i++)
  {
    VH::simulateGaussianInPlace(_workNoiseMesh);
    std::fill(_workmesh.begin(), _workmesh.end(), 0.);
    logPoly.addEvalOp(this, _workNoiseMesh, _workmesh);
    val += VH::innerProduct(_workNoiseMesh, _workmesh);
  }
  return val / nbsimu;
}

double ASPDEOp::computeQuadratic(const std::vector<double>& x) const
{
  _workdat4.resize(_getNDat());
  vect w1s(_workdat4);
  constvect xm(x);
  evalInvCov(xm, w1s);
  return VH::innerProduct(w1s, xm);
}

double ASPDEOp::computeLogDetQ(Id nMC) const
{
  return _QKriging->computeLogDet(nMC);
}

double ASPDEOp::computeLogDetInvNoise() const
{
  return _invNoise->computeLogDet();
}

// We use the fact that log|Sigma| = log |Q + A^t diag^(-1) (sigma) A|- log|Q| + log|Noise|
double ASPDEOp::computeTotalLogDet(Id nMC, Id seed) const
{
  auto memo = law_get_random_seed();

  law_set_random_seed(seed);
  double a1 = computeLogDetOp(nMC);
  double a2 = computeLogDetQ(nMC);
  double a3 = computeLogDetInvNoise();
  law_set_random_seed(memo);
  if (_verbose)
  {
    message("LogDet of Q + ADA': %lf\n", a1);
    message("LogDet of Q: %lf\n", a2);
    message("LogDet of InvNoise: %lf\n", a3);
  }
  double result = TEST;
  if (!FFFF(a1) && !FFFF(a2) && !FFFF(a3)) result = a1 - a2 - a3; // -a3 since a3 is log|invNoise|
  return result;
}
} // namespace gstlrn