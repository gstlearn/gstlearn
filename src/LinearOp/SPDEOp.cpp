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
#include "Basic/Law.hpp"
#include "LinearOp/ProjMulti.hpp"
#include "LinearOp/PrecisionOpMulti.hpp"
#include "Matrix/MatrixDense.hpp"
#include "geoslib_define.h"

ASPDEOp::ASPDEOp(const PrecisionOpMulti* const popKriging,
                 const ProjMulti* const projInKriging,
                 const ASimulable* const invNoise,
                 const PrecisionOpMulti* const popSimu,
                 const ProjMulti* const projInSimu,
                 const ProjMulti* const projOutKriging,
                 const ProjMulti* const projOutSimu,
                 bool todelete)
  : _QKriging(popKriging)
  , _projInKriging(projInKriging)
  , _invNoise(invNoise)
  , _QSimu(popSimu == nullptr ? popKriging : popSimu)
  , _projInSimu(projInSimu == nullptr ? projInKriging : projInSimu)
  , _projOutKriging(projOutKriging)
  , _projOutSimu(projOutSimu == nullptr ? projOutKriging : projOutSimu)
  , _solver(nullptr)
  , _noiseToDelete(todelete)
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
// int ASPDEOp::_addToDest(const constvect inv, vect outv) const
// {
//   _prepare();
  
//   vect w1s(_workdat1);
//   vect w2s(_workdat2);
//   _projKriging->mesh2point(inv, w1s);
//   _invNoise->evalDirect(w1s, w2s);
//   _projKriging->addPoint2mesh(w2s, outv);
//   return _QKriging->addToDest(inv, outv);
// }

int ASPDEOp::_addToDest(const constvect inv, vect outv) const
{
  _prepare();

  int status = _QKriging->addToDest(inv, outv); // TODO: find why outv is set to zero in multistructure case
  if (status) return status;
  vect w1s(_workdat1);
  vect w2s(_workdat2);
  _projInKriging->mesh2point(inv, w1s);
  _invNoise->evalDirect(w1s, w2s);
  _projInKriging->addPoint2mesh(w2s, outv);

  return status;
}

int ASPDEOp::getSizeSimu() const
{
  return _QSimu->getSize();
}

void ASPDEOp::_simCond(const constvect data, vect outvK, vect outvS) const
{
  // Resize if necessary
  _workdat3.resize(_getNDat());
  _workdat4.resize(_getNDat());
  _workmesh.resize(getSizeSimu());
  _workNoiseMesh.resize(getSizeSimu());
  _workNoiseData.resize(_getNDat());

  // Non conditional simulation on mesh
  VH::simulateGaussianInPlace(_workNoiseMesh);
  _QSimu->evalSimulate(_workNoiseMesh, outvS); 
  
  // Simulation at data locations (projection + noise)
  _projInSimu->mesh2point(outvS, _workdat3); //Projection on data locations
  VH::simulateGaussianInPlace(_workNoiseData);
  _invNoise->addSimulateToDest(_workNoiseData, _workdat3); //Add noise
  
  // compute residual _workdat4 = data - outv
  VH::subtractInPlace(_workdat3, data, _workdat4);

  //Co-Kriging of the residual on the mesh
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

VectorDouble ASPDEOp::kriging(const VectorDouble& dat) const
{
  constvect datm(dat.data(), dat.size());
  VectorDouble outMeshK(_QKriging->getSize());
  vect outvs(outMeshK);
  int err = _kriging(datm, outvs);
  if (err) return VectorDouble();

  if (_projOutKriging == nullptr)
    return outMeshK;

  // Project the result on the output mesh
  VectorDouble result(_projOutKriging->getNPoint());
  _projOutKriging->mesh2point(outvs, result);
  return result;
}

int ASPDEOp::centerDataByDriftMat(VectorDouble& Z,
                                  const MatrixDense& driftMat,
                                  const VectorDouble& driftCoeffs)
{
  int nrows = driftMat.getNRows();
  int ncols = driftMat.getNCols();
  if (nrows != (int) Z.size())
  {
    messerr("Error in number of Rows of drift matrix (%d) and size of data vector (%d)",
            nrows, Z.size());
    return 1;
  }
  if (ncols != (int) driftCoeffs.size())
  {
    messerr("Error in number of Columns of drift matrix (%d) and size of drift coefficients (%d)",
            ncols, driftCoeffs.size());
    return 1;
  }

  // Center the data set
  for (int i = 0; i < nrows; i++)
  {
    double sum = 0.;
    for (int j = 0; j < ncols; j++)
    {
      sum += driftCoeffs[j] * driftMat.getValue(i, j);
    }
    Z[i] -= sum;
  }
  return 0;
}

int ASPDEOp::centerDataByMeanVec(VectorDouble& Z,
                                 const VectorDouble& meanVec)
{
  if ((int)Z.size() != (int)meanVec.size())
  {
    messerr("Error in size of data vector (%d) and size of mean vector (%d)",
            Z.size(), meanVec.size());
    return 1;
  }

  // Center the data set
  for (int i = 0, nrows = (int)Z.size(); i < nrows; i++)
    Z[i] -= meanVec[i];
  return 0;
}

VectorDouble ASPDEOp::simCond(const VectorDouble& dat) const
{
  constvect datm(dat.data(), dat.size());
  VectorDouble outMeshK(_QKriging->getSize());
  vect outvK(outMeshK);
  VectorDouble outMeshS(_QSimu->getSize());
  vect outvS(outMeshS);
  _simCond(datm, outvK, outvS);

  if (_projOutKriging == nullptr && _projOutSimu == nullptr)
  {
    VH::addInPlace(outvS, outvK);
    return outMeshK;
  }
  VectorDouble result(_projOutSimu->getNPoint());
  _projOutKriging->mesh2point(outvK, result);
  _projOutSimu->addMesh2point(outvS, result);
  return result;
}

VectorDouble ASPDEOp::simNonCond() const
{
  VectorDouble outMeshS(_QSimu->getSize());
  vect outvs(outMeshS);
  _simNonCond(outvs);

  if (_projOutSimu == nullptr)
    return outMeshS;
  VectorDouble result(_projOutSimu->getNPoint());
  _projOutSimu->mesh2point(outvs, result);
  return result;
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

int ASPDEOp::_kriging(const constvect inv, vect out) const
{
  _buildRhs(inv);
  int status = _solve(_rhs, out);
  return status;
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
  _projInKriging->point2mesh(_workdat1, _rhs);
  return 0;
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
  _projInKriging->point2mesh(result,rhss);
  _solve(rhss,wms);
  _projInKriging->mesh2point(wms,w2s);
  _invNoise->evalDirect(w2s, wms);
  VectorHelper::subtractInPlace(wms, result, result);
}

VectorDouble ASPDEOp::computeDriftCoeffs(const VectorDouble& Z,
                                         const MatrixDense& driftMat,
                                         bool verbose) const
{
  int xsize = (int)(driftMat.getNCols());
  VectorDouble XtInvSigmaZ(xsize);
  MatrixSymmetric XtInvSigmaX(xsize);
  VectorDouble result(xsize);

  _workdat1.resize(_getNDat());
  vect w1s(_workdat1);
  for(int i = 0; i< xsize; i++)
  {
    auto xm = driftMat.getColumnPtr(i);
    evalInvCov(xm, w1s);

    constvect ym(Z.data(),Z.size());
    constvect wd1(_workdat1.data(),_workdat1.size());
    XtInvSigmaZ[i] = VH::innerProduct(ym,wd1);

    for(int j = i; j < xsize;j++)
    {
      constvect xmj = driftMat.getViewOnColumn(j);
      double prod   = VH::innerProduct(xmj, w1s);
      XtInvSigmaX.setValue(i,j,prod);
    }
  }

  XtInvSigmaX.solve(XtInvSigmaZ,result);

  // Optional printout
  if (verbose)
    VH::dump("Drift coefficients", result);
  
  return result;
}

double ASPDEOp::computeLogDetOp(int nbsimu) const
{
  DECLARE_UNUSED(nbsimu);
  messerr("Not implemented yet in Matrix-free version");
  return TEST;
}

double ASPDEOp::computeQuadratic(const std::vector<double>& x) const
{
  _workdat1.resize(_getNDat());
  vect w1s(_workdat1);
  constvect xm(x);
  evalInvCov(xm, w1s);
  return VH::innerProduct(w1s, xm);
}

double ASPDEOp::computeLogDetQ(int nMC) const
{
  return _QKriging->computeLogDetQ(nMC);
}

double ASPDEOp::computeLogDetNoise() const
{
  return _invNoise->computeLogDet();
}

// We use the fact that log|Sigma| = log |Q + A^t diag^(-1) (sigma) A|- log|Q| + Sum(log sigma_i^2)
double ASPDEOp::computeTotalLogDet(int nMC, int seed) const
{
  int memo = law_get_random_seed();

  law_set_random_seed(seed);
  double a1 = computeLogDetOp(nMC);
  double a2 = computeLogDetQ(nMC);
  double a3 = computeLogDetNoise();
  law_set_random_seed(memo);

  double result = TEST;
  if (! FFFF(a1) && ! FFFF(a2) && ! FFFF(a3)) result = a1 - a2 + a3;
  return result;
}
