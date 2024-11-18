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
#include "Covariances/ACov.hpp"
#include "Covariances/CovCalcMode.hpp"
#include "Matrix/MatrixSquareGeneral.hpp"
#include "Matrix/MatrixRectangular.hpp"
#include "Matrix/MatrixSparse.hpp"
#include "Matrix/MatrixSquareSymmetric.hpp"
#include "Matrix/NF_Triplet.hpp"
#include "Db/Db.hpp"
#include "Db/DbGrid.hpp"
#include "Basic/AException.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/VectorNumT.hpp"
#include "Basic/VectorHelper.hpp"
#include "Basic/Law.hpp"
#include "Space/ASpace.hpp"
#include "Space/SpacePoint.hpp"
#include "geoslib_define.h"

#include <vector>
#include <math.h>

ACov::ACov(const ASpace *space)
    : ASpaceObject(space),
      _optimEnabled(true),
      _isOptimPreProcessed(false),
      _p1As(),
      _p2A(space)
{
}

ACov::ACov(const ACov &r)
    : ASpaceObject(r),
      _optimEnabled(r._optimEnabled),
      _isOptimPreProcessed(false),
      _p1As(),
      _p2A(r.getSpace())
{
}

ACov& ACov::operator=(const ACov &r)
{
  if (this != &r)
  {
    ASpaceObject::operator=(r);
    _optimEnabled = r._optimEnabled;
    _isOptimPreProcessed = r._isOptimPreProcessed;
    _p2A = SpacePoint(r._space);
  }
  return *this;
}

ACov::~ACov()
{
}

void ACov::optimizationPostProcess() const
{
  _optimizationPostProcess();
  _isOptimPreProcessed = false;
}

void ACov::optimizationPreProcess(const Db* db) const
{
  if (_isOptimPreProcessed) return;
  db->getSamplesAsSP(_p1As,_space);
  _optimizationPreProcess(_p1As);
  
}

void ACov::optimizationSetTarget(const SpacePoint &pt) const
{
  _optimizationSetTarget(pt);
}

void ACov::_optimizationSetTarget(const SpacePoint &pt) const
{
  _p2A = pt;
}


void ACov::optimizationPreProcess(const std::vector<SpacePoint>& p) const
{
  _optimizationPreProcess(p);
  _isOptimPreProcessed = isOptimEnabled();
}

//Preprocess by default (only copy the points)
void ACov::_optimizationPreProcess(const std::vector<SpacePoint>& p) const
{
  if (!_p1As.empty() && _isOptimPreProcessed) return;
  _p1As.clear();
  for (const auto &e : p)
  {
   _p1As.push_back(e);
  }
}

void ACov::_optimizationPostProcess() const
{
    _p1As.clear();
}

MatrixSquareGeneral ACov::eval0Mat(const CovCalcMode* mode) const
{
  int nvar = getNVariables();
  MatrixSquareGeneral mat(nvar);
  eval0CovMatBiPointInPlace(mat,mode);
  return mat;
}

/**
 * Calculate the Matrix of covariance for zero distance
 * @param mat   Covariance matrix (Dimension: nvar * nvar)
 * @param mode  Calculation Options
 *
 * @remarks: Matrix 'mat' should be dimensioned and initialized beforehand
 */
void ACov::eval0CovMatBiPointInPlace(MatrixSquareGeneral& mat, const CovCalcMode* mode) const
{
  mat.fill(0.);
  addEval0CovMatBiPointInPlace(mat, mode);
}

/**
 * Add the contribution of the Matrix of covariance for zero distance
 * @param mat   Covariance matrix (Dimension: nvar * nvar)
 * @param mode  Calculation Options
 *
 * @remarks: Matrix 'mat' should be dimensioned and initialized beforehand
 */
void ACov::addEval0CovMatBiPointInPlace(MatrixSquareGeneral& mat, const CovCalcMode* mode) const
{
  int nvar = getNVariables();

  for (int ivar=0; ivar<nvar; ivar++)
    for (int jvar=0; jvar<nvar; jvar++)
        mat.addValue(ivar, jvar, eval0(ivar, jvar, mode)); // pure virtual method
}

VectorDouble ACov::eval(const std::vector<SpacePoint>& vec_p1,
                        const std::vector<SpacePoint>& vec_p2,
                        int ivar,
                        int jvar,
                        const CovCalcMode* mode) const
{
  VectorDouble vec;
  if (vec_p1.size() != vec_p2.size())
    my_throw ("Error: 'p1' and 'p2' should have same dimension");
  for (int i=0, n=static_cast<int> (vec_p1.size()); i < n; i++)
    vec.push_back(eval(vec_p1[i], vec_p2[i], ivar, jvar, mode)); // pure virtual method
  return vec;
}

double ACov::eval0(int ivar,
                   int jvar,
                   const CovCalcMode* mode) const
{
  SpacePoint p1(getSpace()->getOrigin(),-1);
  return eval(p1,p1,ivar,jvar,mode); // pure virtual method
}
MatrixSquareGeneral ACov::evalMat(const SpacePoint& p1,
                                  const SpacePoint& p2,
                                  const CovCalcMode* mode) const
{
  int nvar = getNVariables();
  MatrixSquareGeneral mat(nvar);
  evalCovMatBiPointInPlace(mat,p1, p2,mode);
  return mat;
}
/**
 * Calculate the Matrix of covariance between two space points
 * @param p1 Reference of the first space point
 * @param p2 Reference of the second space point
 * @param mat   Covariance matrix (Dimension: nvar * nvar)
 * @param mode  Calculation Options
 *
 * @remarks: Matrix 'mat' should be dimensioned and initialized beforehand
 */
void ACov::evalCovMatBiPointInPlace(MatrixSquareGeneral &mat,
                                    const SpacePoint &p1,
                                    const SpacePoint &p2,
                                    const CovCalcMode* mode) const
{
  mat.fill(0.);
  addEvalCovMatBiPointInPlace(mat,p1,p2,mode);
}

/**
 * Covariance from a given point (center) in a given direction (dir *step)
 * @param ivar   Rank of the first variable
 * @param jvar   Rank of the second variable
 * @param step   Step value
 * @param dir    Direction definition
 * @param mode   CovCalcMode structure
 * @return
 */
double ACov::evalIvarIpas(double step,
                          const VectorDouble& dir,
                          int ivar,
                          int jvar,
                          const CovCalcMode* mode) const
{
  // Define the point in the ACov space (center will be checked)
  const ASpace* space = getSpace();
  SpacePoint p1(space);
  SpacePoint p2(space);

  if (dir.empty())
  {
    VectorDouble vec(getNDim(),0.);
    vec[0] = 1.;
    VH::multiplyConstant(vec, step);
    p2.move(vec);
  }
  else
  {
    VectorDouble vec(dir);
    VH::multiplyConstant(vec, step);
    p2.move(vec);
  }

  return eval(p1, p2, ivar, jvar, mode); // pure virtual method
}

double ACov::evalIvarIpasIncr(const VectorDouble& dincr,
                              int ivar,
                              int jvar,
                              const CovCalcMode* mode) const
{
  // Define the point in the ACov space (center will be checked)
  SpacePoint p1(VectorDouble(_space->getNDim()),-1,_space);
  SpacePoint p2(VectorDouble(_space->getNDim()),-1,_space);
  p2.move(dincr);
  return eval(p1, p2, ivar, jvar, mode); // pure virtual method
}


/**
 * Covariance vector from a given point (center) in a given direction (dir * steps)
 * for a pair of variables and a set of steps
 * @param ivar      Rank of the first variable
 * @param jvar      Rank of the second variable
 * @param vec_step  Vector of step values
 * @param dir       Direction definition
 * @param mode      CovCalcMode structure
 * @return
 */
VectorDouble ACov::evalIvarNpas(const VectorDouble& vec_step,
                                const VectorDouble& dir,
                                int ivar,
                                int jvar,
                                const CovCalcMode* mode) const
{
  VectorDouble vec;
  for (int i=0, n=static_cast<int> (vec_step.size()); i < n; i++)
    vec.push_back(evalIvarIpas(vec_step[i], dir, ivar, jvar, mode));
  return vec;
}

/**
 * Covariance Matrix from a given point (center) in a given direction (dir * step)
 * for a set of variables and a given step
 * @param step   Step value
 * @param dir    Direction definition
 * @param mode   CovCalcMode structure
 * @return
 */
MatrixSquareGeneral ACov::evalNvarIpas(double step,
                                       const VectorDouble& dir,
                                       const CovCalcMode* mode) const
{
  int nvar = getNVariables();
  MatrixSquareGeneral mat(nvar);
  for (int ivar=0; ivar<nvar; ivar++)
    for (int jvar=0; jvar<nvar; jvar++)
      mat.setValue(ivar, jvar, evalIvarIpas(step, dir, ivar, jvar, mode));
  return mat;
}

MatrixSquareGeneral ACov::evalNvarIpasIncr(const VectorDouble& dincr,
                                           const CovCalcMode* mode) const
{
  int nvar = getNVariables();
  MatrixSquareGeneral mat(nvar);
  for (int ivar=0; ivar<nvar; ivar++)
    for (int jvar=0; jvar<nvar; jvar++)
      mat.setValue(ivar, jvar, evalIvarIpasIncr(dincr, ivar, jvar, mode));
  return mat;
}

/**
 * Covariance for a given unit global distance (without anisotropy)
 * for a pair of variables and a single step
 * @param ivar Rank of the first variable
 * @param jvar Rank of the second variable
 * @param step Step value
 * @param mode CovCalcMode structure
 * @return
 */
double ACov::evalIsoIvarIpas(double step,
                             int ivar,
                             int jvar,
                             const CovCalcMode* mode) const
{
  /// TODO : Not true whatever the space
  VectorDouble center = getOrigin();
  VectorDouble dir = getUnitaryVector();
  return evalIvarIpas(step, dir, ivar, jvar, mode);
}

/**
 * Covariance for a given unit global distance (without anisotropy)
 * for a pair of variables and a set of steps
 * @param ivar
 * @param jvar
 * @param vec_step
 * @param mode
 * @return
 */
VectorDouble ACov::evalIsoIvarNpas(const VectorDouble& vec_step,
                                   int ivar,
                                   int jvar,
                                   const CovCalcMode* mode) const
{
  VectorDouble vec;
  VectorDouble dir = getUnitaryVector();
  for (const auto& h : vec_step)
    vec.push_back(evalIvarIpas(h, dir, ivar, jvar, mode));
  return vec;
}

/**
 * Covariance for a given unit global distance (without anisotropy)
 * for a set of variables and a single step
 * @param step Step value
 * @param mode CovCalcMode structure
 * @return
 */
MatrixSquareGeneral ACov::evalIsoNvarIpas(double step,
                                          const CovCalcMode* mode) const
{
  int nvar = getNVariables();
  VectorDouble dir = getUnitaryVector();
  MatrixSquareGeneral mat(nvar);
  for (int ivar=0; ivar<nvar; ivar++)
    for (int jvar=0; jvar<nvar; jvar++)
      mat.setValue(ivar, jvar, evalIvarIpas(step, dir, ivar, jvar, mode));
  return mat;
}

/**
 * Calculate the (weighted) average Covariance between samples of two Dbs,
 * for a pair of variables
 * @param db1  Pointer to the first Db
 * @param db2  Pointer to the second Db
 * @param ivar Rank of the first variables
 * @param jvar Rank of the second variable
 * @param eps  Epsilon used for randomization in calculation of CVV (optional)
 * @param seed Seed for the randomization
 * @param mode CovCalcMode structure
 * @return
 */
double ACov::evalAverageDbToDb(const Db* db1,
                               const Db* db2,
                               int ivar,
                               int jvar,
                               double eps,
                               int seed,
                               const CovCalcMode* mode) const
{
  int memo = law_get_random_seed();
  if (eps > 0. && seed > 0)
    law_set_random_seed(seed);

  /* Loop on the first sample */

  double norme = 0.;
  double total = 0.;
  for (int iech1 = 0; iech1 < db1->getSampleNumber(); iech1++)
  {
    if (!db1->isActive(iech1)) continue;
    double w1 = db1->getWeight(iech1);
    if (isZero(w1)) continue;
    SpacePoint p1(db1->getSampleCoordinates(iech1));

    /* Loop on the second sample */

    for (int iech2 = 0; iech2 < db2->getSampleNumber(); iech2++)
    {
      if (!db2->isActive(iech2)) continue;
      double w2 = db2->getWeight(iech2);
      if (isZero(w2)) continue;
      VectorDouble coord2 = db2->getSampleCoordinates(iech2);
      if (eps > 0)
      {
        for (int idim = 0, ndim = getNDim(); idim < ndim; idim++)
          coord2[idim] += eps * law_uniform(-0.5,  0.5);
      }
      SpacePoint p2(coord2);

      /* Loop on the dimension of the space */

      total += w1 * w2 * eval(p1, p2, ivar, jvar, mode);
      norme += w1 * w2;
    }
  }

  // Scaling
  if (!isZero(norme)) total /= norme;

  if (eps > 0. && seed > 0)
    law_set_random_seed(memo);

  return total;
}

double ACov::evalAverageIncrToIncr(const VectorVectorDouble &d1,
                                   const VectorVectorDouble &d2,
                                   int ivar,
                                   int jvar,
                                   const CovCalcMode *mode) const
{
  int nech1 = (int) d1.size();
  int nech2 = (int) d2.size();

  /* Loop on the first sample */

  double total = 0.;
  for (int iech1 = 0; iech1 < nech1; iech1++)
  {
    SpacePoint p1(d1[iech1],-1,getSpace());

    /* Loop on the second sample */

    for (int iech2 = 0; iech2 < nech2; iech2++)
    {
      SpacePoint p2(d2[iech2],-1,getSpace());
      total += eval(p1, p2, ivar, jvar, mode);
    }
  }

  // Scaling
  total /= (double) (nech1 * nech2);

  return total;
}

/**
 * Calculate the (weighted) average Covariance between a point and a Db
 * for a pair of variables
 * @param p1   Point location
 * @param db2  Pointer to the second Db
 * @param ivar Rank of the first variables
 * @param jvar Rank of the second variable
 * @param mode CovCalcMode structure
 * @return
 */
double ACov::evalAveragePointToDb(const SpacePoint& p1,
                                  const Db* db2,
                                  int ivar,
                                  int jvar,
                                  const CovCalcMode* mode) const
{
  /* Loop on the first sample */

  double norme = 0.;
  double total = 0.;

  /* Loop on the second sample */

  for (int iech2 = 0; iech2 < db2->getSampleNumber(); iech2++)
  {
    if (!db2->isActive(iech2)) continue;
    double w2 = db2->getWeight(iech2);
    if (isZero(w2)) continue;
    SpacePoint p2(db2->getSampleCoordinates(iech2),iech2,getSpace());

    /* Loop on the dimension of the space */

    total += w2 * eval(p1, p2, ivar, jvar, mode);
    norme += w2;
  }

  // Scaling
  if (isZero(norme)) total /= norme;

  return total;
}

VectorDouble ACov::evalPointToDbAsSP(const std::vector<SpacePoint> &p1s,
                                     const SpacePoint &p2,
                                     int ivar,
                                     int jvar,
                                     const CovCalcMode *mode) const
{
  int nech1 = (int) p1s.size();
  VectorDouble values(nech1);

  /* Loop on the second sample */

  for (int iech1 = 0; iech1 < nech1; iech1++)
  {
    p1s[iech1].setIech(iech1);
    values[iech1] = eval(p1s[iech1], p2, ivar, jvar, mode);
  }
  return values;
}

/**
 * Calculate the Covariance vector between a Point and all the samples
 * of a Db, for a pair of variables
 * @param p1   Point location
 * @param db2  Pointer to the second Db
 * @param ivar Rank of the first variables
 * @param jvar Rank of the second variable
 * @param useSel When TRUE, the returned vector is reduced to active samples
 *               Otherwise, returns TEST for masked samples
 * @param nbgh2 Vector of indices of active samples in db2 (optional)
 * @param mode CovCalcMode structure
 * @return
 */
VectorDouble ACov::evalPointToDb(const SpacePoint& p1,
                                 const Db* db2,
                                 int ivar,
                                 int jvar,
                                 bool useSel,
                                 const VectorInt& nbgh2,
                                 const CovCalcMode* mode) const
{
  VectorDouble values;
  int nech2;
  if (nbgh2.empty())
    nech2 = db2->getSampleNumber();
  else
    nech2 = (int) nbgh2.size();

  /* Loop on the second sample */

  for (int kech2 = 0; kech2 < nech2; kech2++)
  {
    int iech2 = (nbgh2.empty()) ? kech2 : nbgh2[kech2];
    if (! nbgh2.empty())
    {
      SpacePoint p2(db2->getSampleCoordinates(iech2),iech2,getSpace());
      values.push_back(eval(p1, p2, ivar, jvar, mode));
    }
    else
    {
      if (db2->isActive(iech2))
      {
        SpacePoint p2(db2->getSampleCoordinates(iech2),iech2, getSpace());
        values.push_back(eval(p1, p2, ivar, jvar, mode));

      }
      else
      {
        if (!useSel) values.push_back(TEST);
      }
    }
  }

  return values;
}

/**
 * Average covariance over a block
 * @param ext    Vector of Block extensions
 * @param ndisc  Vector of Block discretization
 * @param angles Vector of rotation angles
 * @param ivar   Rank of the first variable
 * @param jvar   Rank of the second variable
 * @param mode   CovCalcMode structure
 * @return
 */
double ACov::evalCvv(const VectorDouble& ext,
                     const VectorInt& ndisc,
                     const VectorDouble& angles,
                     int ivar,
                     int jvar,
                     const CovCalcMode* mode) const
{
  int ndim = getNDim();
  if (ndim != (int) ext.size())
  {
    messerr("Block Extension (%d) should have same dimension as the Model %d)",
            (int) ext.size(),ndim);
    return TEST;
  }
  if (ndim != (int) ndisc.size())
  {
    messerr("Discretization (%d) should have same dimension as the Model (%d)",
            (int) ndisc.size(), ndim);
    return TEST;
  }

  DbGrid* dbgrid = _discretizeBlock(ext, ndisc, angles);
  if (dbgrid == nullptr) return TEST;
  Db* db = _discretizeBlockRandom(dbgrid);
  if (db == nullptr) return TEST;

  double total = evalAverageDbToDb(dbgrid,  db, ivar, jvar, 0., 0, mode);
  delete dbgrid;
  return total;
}

/**
 * Average covariance between a block and the same block shifted
 * @param ext    Vector of Block extensions
 * @param ndisc  Vector of Block discretization
 * @param angles Vector of rotation angles
 * @param shift  Shift between the two blocks
 * @param ivar   Rank of the first variable
 * @param jvar   Rank of the second variable
 * @param mode   CovCalcMode structure
 * @return
 */
double ACov::evalCvvShift(const VectorDouble& ext,
                          const VectorInt& ndisc,
                          const VectorDouble& shift,
                          const VectorDouble& angles,
                          int ivar,
                          int jvar,
                          const CovCalcMode* mode) const
{
  int ndim = getNDim();
  if (ndim != (int) ext.size())
  {
    messerr("Block Extension (%d) should have same dimension as the Model %d)",
            (int) ext.size(),ndim);
    return TEST;
  }
  if (ndim != (int) ndisc.size())
  {
    messerr("Discretization (%d) should have same dimension as the Model (%d)",
            (int) ndisc.size(), ndim);
    return TEST;
  }
  if (ndim != (int) shift.size())
  {
    messerr("Shift (%d) should have the same dimension as the Model (%d)",
            (int) shift.size(), ndim);
    return TEST;
  }

  DbGrid* dbgrid1 = _discretizeBlock(ext, ndisc, angles);
  if (dbgrid1 == nullptr) return TEST;
  DbGrid* dbgrid2 = _discretizeBlock(ext, ndisc, angles, shift);
  if (dbgrid2 == nullptr) return TEST;

  double total = evalAverageDbToDb(dbgrid1,  dbgrid2, ivar, jvar, 0., 0, mode);
  delete dbgrid1;
  delete dbgrid2;
  return total;
}

MatrixSquareGeneral ACov::evalCvvM(const VectorDouble& ext,
                                   const VectorInt& ndisc,
                                   const VectorDouble& angles,
                                   const CovCalcMode* mode) const
{
  int nvar = getNVariables();
  MatrixSquareGeneral mat(nvar);
  for (int ivar=0; ivar<nvar; ivar++)
    for (int jvar=0; jvar<nvar; jvar++)
      mat.setValue(ivar, jvar, evalCvv(ext, ndisc, angles, ivar, jvar, mode));
  return mat;
}

/**
 * Average covariance over a block
 * @param p1     Point location
 * @param ext    Vector of Block extensions
 * @param ndisc  Vector of Block discretization
 * @param angles Vector of rotation angles
 * @param x0     Vector for origin of block
 * @param ivar   Rank of the first variable
 * @param jvar   Rank of the second variable
 * @param mode   CovCalcMode structure
 * @return
 */
double ACov::evalCxv(const SpacePoint& p1,
                     const VectorDouble& ext,
                     const VectorInt& ndisc,
                     const VectorDouble& angles,
                     const VectorDouble& x0,
                     int ivar,
                     int jvar,
                     const CovCalcMode* mode) const
{
  int ndim = getNDim();
  if (ndim != (int) ext.size())
  {
    messerr("Block Extension (%d) should have same dimension as the Model %d)",
            (int) ext.size(),ndim);
    return TEST;
  }
  if (ndim != (int) ndisc.size())
  {
    messerr("Discretization (%d) should have same dimension as the Model (%d)",
            (int) ndisc.size(), ndim);
    return TEST;
  }

  double total = TEST;
  DbGrid* dbgrid = _discretizeBlock(ext, ndisc, angles, x0);
  if (dbgrid != nullptr)
    total = evalAveragePointToDb(p1, dbgrid, ivar, jvar, mode);
  delete dbgrid;

  return total;
}

double ACov::evalCxv(const Db* db,
                     const VectorDouble& ext,
                     const VectorInt& ndisc,
                     const VectorDouble& angles,
                     const VectorDouble& x0,
                     int ivar,
                     int jvar,
                     const CovCalcMode* mode) const
{
  int ndim = getNDim();
  if (db == nullptr)
  {
    messerr("Argument 'db' should be defined");
    return TEST;
  }
  if (ndim != db->getNDim())
  {
    messerr("Db (%d) should have the seame dimension as the Model(%d)",
            db->getNDim(), ndim);
    return TEST;
  }
  if (ndim != (int) ext.size())
  {
    messerr("Block Extension (%d) should have same dimension as the Model %d)",
            (int) ext.size(),ndim);
    return TEST;
  }
  if (ndim != (int) ndisc.size())
  {
    messerr("Discretization (%d) should have same dimension as the Model (%d)",
            (int) ndisc.size(), ndim);
    return TEST;
  }

  double total = TEST;
  DbGrid* dbgrid = _discretizeBlock(ext, ndisc, angles, x0);
  if (dbgrid != nullptr)
    total = evalAverageDbToDb(db, dbgrid, ivar, jvar, 0., 0, mode);
  delete dbgrid;

  return total;
}

MatrixSquareGeneral ACov::evalCxvM(const SpacePoint& p1,
                                   const VectorDouble& ext,
                                   const VectorInt& ndisc,
                                   const VectorDouble& angles,
                                   const VectorDouble& x0,
                                   const CovCalcMode* mode) const
{
  int nvar = getNVariables();
  MatrixSquareGeneral mat(nvar);
  for (int ivar=0; ivar<nvar; ivar++)
    for (int jvar=0; jvar<nvar; jvar++)
      mat.setValue(ivar, jvar, evalCxv(p1, ext, ndisc, angles, x0, ivar, jvar, mode));
  return mat;
}

/**
 * Creates the discretization grid
 * @param ext    Vecto of Block extensions
 * @param ndisc  Vector of Discretizations
 * @param angles Vector of rotation angles
 * @param x0     Vector of Discretization origin
 * @return
 *
 * @remark If block origin is not defined, it is set so that the
 * @remark center of the block is one the point (0,0)
 */
DbGrid* ACov::_discretizeBlock(const VectorDouble& ext,
                               const VectorInt& ndisc,
                               const VectorDouble& angles,
                               const VectorDouble& x0) const
{
  int ndim = getNDim();
  VectorDouble x0loc = x0;
  if (x0loc.empty() || ndim != (int) x0loc.size())
    x0loc.resize(ndim, 0.);
  for (int idim = 0; idim < ndim; idim++)
    x0loc[idim] -= ext[idim] / 2.;
  VectorDouble dx(ndim, 0.);
  for (int idim = 0; idim < ndim; idim++)
    dx[idim] = ext[idim] / ndisc[idim];
  DbGrid* dbgrid = DbGrid::create(ndisc, dx, x0loc, angles);
  return dbgrid;
}

Db* ACov::_discretizeBlockRandom(const DbGrid* dbgrid, int seed) const
{
  int ndim = getNDim();
  int nech = dbgrid->getSampleNumber();
  Db* db = Db::createFromSamples(nech);
  VectorString names = generateMultipleNames("x",ndim);
  law_set_random_seed(seed);

  for (int idim = 0; idim < ndim; idim++)
  {
    double taille = dbgrid->getDX(idim);
    VectorDouble vec = dbgrid->getCoordinates(idim, false);
    for (int i = 0; i < (int) vec.size(); i++)
      vec[i] += taille * law_uniform(-0.5, 0.5);
    db->addColumns(vec, names[idim], ELoc::X, idim);
  }
  return db;
}

VectorInt ACov::_getActiveVariables(int ivar0) const
{
  int nvar = getNVariables();

  VectorInt ivars;
  if (ivar0 >= 0)
  {
    if (!checkArg("Argument 'ivar0'", ivar0, nvar)) return VectorInt();
    ivars.push_back(ivar0);
  }
  else
  {
    ivars = VH::sequence(nvar);
  }
  return ivars;
}

/****************************************************************************/
/*!
 **  Establish the covariance matrix between two Dbs
 **  Takes into account selection and heterotopy
 **
 ** \return Dense matrix containing the covariance matrix
 **
 ** \param[in]  db1   First Db
 ** \param[in]  db2   Second Db (= db1 if absent)
 ** \param[in]  ivar0 Rank of the first variable (-1 for all variables)
 ** \param[in]  jvar0 Rank of the second variable (-1 for all variables)
 ** \param[in]  nbgh1 Vector of indices of active samples in db1 (optional)
 ** \param[in]  nbgh2 Vector of indices of active samples in db2 (optional)
 ** \param[in]  mode  CovCalcMode structure
 **
 ** \remarks If a Db does not contain any Z-variable defined, the covariance
 ** \remarks cannot treat possible heterotopy and therefore uses all samples
 **
 ** \remarks The returned matrix if dimension to nrows * ncols where
 ** \remarks each term is the product of the number of active samples
 ** \remarks by the number of samples where the variable is defined
 **
 ** \note 'dbin' and 'dbout' cannot be made 'const' as they can be updated
 ** \note due to the presence of 'nostat'
 **
 *****************************************************************************/
MatrixRectangular ACov::evalCovMatrix(const Db* db1,
                                      const Db* db2,
                                      int ivar0,
                                      int jvar0,
                                      const VectorInt& nbgh1,
                                      const VectorInt& nbgh2,
                                      const CovCalcMode* mode)
{
  MatrixRectangular mat;

  // Preliminary checks
  if (db2 == nullptr) db2 = db1;
  if (db1 == nullptr || db2 == nullptr) return MatrixRectangular();
  VectorInt ivars = _getActiveVariables(ivar0);
  if (ivars.empty()) return mat;
  VectorInt jvars = _getActiveVariables(jvar0);
  if (jvars.empty()) return mat;

  // Play the non-stationarity (if needed)
  manage(db1,db2);
  
  // Create the sets of Vector of valid sample indices per variable (not masked and defined)
  VectorVectorInt index1 = db1->getMultipleRanksActive(ivars, nbgh1);
  VectorVectorInt index2 = db2->getMultipleRanksActive(jvars, nbgh2);

  // Creating the matrix
  int neq1 = VH::count(index1);
  int neq2 = VH::count(index2);
  if (neq1 <= 0 || neq2 <= 0)
  {
    messerr("The returned matrix does not have any valid sample for any valid variable");
    return mat;
  }
  mat.resize(neq1, neq2);

  // Define the two space points
  SpacePoint p1(getSpace());
  SpacePoint p2(getSpace());

  // Loop on the first variable
  int irow = 0;
  for (int ivar = 0, nvar1 = (int) ivars.size(); ivar < nvar1; ivar++)
  {
    int ivar1 = ivars[ivar];

    // Loop on the first sample
    int nech1s = (int) index1[ivar].size();
    for (int jech1 = 0; jech1 < nech1s; jech1++)
    {
      int iech1 = index1[ivar][jech1];
      p1.setIech(iech1);
      db1->getSampleAsSPInPlace(p1);

      // Loop on the second variable
      int icol = 0;
      for (int jvar = 0, nvar2 = (int) jvars.size(); jvar < nvar2; jvar++)
      {
        int jvar2 = jvars[jvar];

        // Loop on the second sample
        int nech2s = (int) index2[jvar].size();
        for (int jech2 = 0; jech2 < nech2s; jech2++)
        {
          int iech2 = index2[jvar][jech2];
          p2.setIech(iech2);
          db2->getSampleAsSPInPlace(p2);

          // Modify the covariance (if non stationary)
          updateCovByPoints(1, iech1, 2, iech2);

          /* Loop on the dimension of the space */
          double value = eval(p1, p2, ivar1, jvar2, mode);
          mat.setValue(irow, icol, value);
          icol++;
        }
      }
      irow++;
    }
  }

  return mat;
}

void ACov::_updateCovMatrixSymmetricVerr(const Db *db1,
                                         AMatrix *mat,
                                         const VectorInt &ivars,
                                         const VectorVectorInt &index1)
{
  // Check if the correction can take place at all
  if (! db1->hasLocVariable(ELoc::V)) return;

  // Initializations
  int icolVerr = -1;
  double verr = 0.;

  // Loop on the first variable
  int irow = 0;
  for (int rvar1 = 0, nvar1 = (int) ivars.size(); rvar1 < nvar1; rvar1++)
  {
    int ivar1 = ivars[rvar1];
    icolVerr = db1->getColIdxByLocator(ELoc::V, ivar1);

    // Loop on the first sample
    int nech1s = (int) index1[rvar1].size();
    for (int rech1 = 0; rech1 < nech1s; rech1++)
    {
      int iech1 = index1[rvar1][rech1];

      // Update the Diagonal due to the presence of Variance of Measurement Error
      if (icolVerr >= 0)
      {
        verr = db1->getValueByColIdx(iech1, icolVerr);
        mat->updValue(irow, irow, EOperator::ADD, verr);
      }
      irow++;
    }
  }
}

void ACov::addEvalCovMatBiPointInPlace(MatrixSquareGeneral &mat,
                        const SpacePoint& pwork1, 
                        const SpacePoint& pwork2, 
                        const CovCalcMode *mode) const
{
  _addEvalCovMatBiPointInPlace(mat, pwork1, pwork2,mode);
}

void ACov::_addEvalCovMatBiPointInPlace(MatrixSquareGeneral &mat,
                        const SpacePoint& pwork1, 
                        const SpacePoint& pwork2, 
                        const CovCalcMode *mode) const
{
  for (int ivar = 0, nvar = getNVariables(); ivar < nvar; ivar++)
    for (int jvar = 0; jvar < nvar; jvar++)
      mat.addValue(ivar, jvar, eval(pwork1,pwork2,ivar,jvar,mode)); 
}

void ACov::evalCovKriging(MatrixSquareGeneral &mat,
                          SpacePoint &pwork1,
                          SpacePoint& pout, 
                          const CovCalcMode *mode) const
{
  mat.fill(0.);
  loadAndAddEvalCovMatBiPointInPlace(mat,pwork1,pout,mode);
}

void ACov::loadAndAddEvalCovMatBiPointInPlace(MatrixSquareGeneral &mat,const SpacePoint& p1,const SpacePoint&p2,
                                              const CovCalcMode *mode) const
{
  _loadAndAddEvalCovMatBiPointInPlace(mat,p1,p2,mode);
}

double ACov::loadAndEval(const SpacePoint& p1,
                          const SpacePoint&p2,
                          int ivar,
                          int jvar,
                          const CovCalcMode *mode) const
{
  return _loadAndEval(p1,p2,ivar,jvar,mode);
}

double ACov::_loadAndEval(const SpacePoint& p1,
                          const SpacePoint&p2,
                          int ivar,
                          int jvar,
                          const CovCalcMode *mode) const
{
  load(p1,true);
  load(p2,false);
  return eval(*_pw1,*_pw2,ivar,jvar,mode);

}
void ACov::_loadAndAddEvalCovMatBiPointInPlace(MatrixSquareGeneral &mat,
                                              const SpacePoint& p1,
                                              const SpacePoint&p2,
                                              const CovCalcMode *mode) const
{
  load(p1,true);
  load(p2,false);
  _addEvalCovMatBiPointInPlace(mat,*_pw1,*_pw2,mode);
}

void ACov::load(const SpacePoint& p,bool case1) const
{
  const SpacePoint** pp = case1? &_pw1:&_pw2;
  *pp = p.isTarget()? &_p2A:&_p1As[p.getIech()];
}


/****************************************************************************/
/*!
 **  Establish the covariance matrix within a Db
 **  Takes into account selection and heterotopy
 **  This method takes advantage of calculating covariance between a Db and itself
 **
 ** \return Dense matrix containing the covariance matrix
 **
 ** \param[in]  db1   First Db
 ** \param[in]  ivar0 Rank of the first variable (-1 for all variables)
 ** \param[in]  nbgh1 Vector of indices of active samples in db1 (optional)
 ** \param[in]  mode  CovCalcMode structure
 **
 ** \remarks If a Db does not contain any Z-variable defined, the covariance
 ** \remarks cannot treat possible heterotopy and therefore uses all samples
 **
 ** \remarks The returned matrix if dimension to nrows * ncols where
 ** \remarks each term is the product of the number of active samples
 ** \remarks by the number of samples where the variable is defined
 **
 *****************************************************************************/
MatrixSquareSymmetric ACov::evalCovMatrixSymmetric(const Db *db1,
                                                   int ivar0,
                                                   const VectorInt &nbgh1,
                                                   const CovCalcMode *mode)
{
  MatrixSquareSymmetric mat;

  // Preliminary checks
  if (db1 == nullptr) return MatrixRectangular();
  VectorInt ivars = _getActiveVariables(ivar0);
  if (ivars.empty()) return mat;

  // Play the non-stationarity (if needed)
  manage(db1,nullptr);
  

  // Create the sets of Vector of valid sample indices per variable (not masked and defined)
  VectorVectorInt index1 = db1->getMultipleRanksActive(ivars, nbgh1, true, true);

  // Creating the matrix
  int neq1 = VH::count(index1);
  if (neq1 <= 0)
  {
    messerr("The returned matrix does not have any valid sample for any valid variable");
    return mat;
  }
  mat.resize(neq1, neq1);

  // Define the two space points
  SpacePoint p1(getSpace());
  SpacePoint p2(getSpace());

  // Loop on the first variable
  int irow = 0;
  for (int rvar1 = 0, nvar1 = (int) ivars.size(); rvar1 < nvar1; rvar1++)
  {
    int ivar1 = ivars[rvar1];

    // Loop on the first sample
    int nech1s = (int) index1[rvar1].size();
    for (int rech1 = 0; rech1 < nech1s; rech1++)
    {
      int iech1 = index1[rvar1][rech1];
      p1.setIech(iech1);
      db1->getSampleAsSPInPlace(p1);

      // Loop on the second variable
      int icol = 0;
      for (int rvar2 = 0, nvar2 = (int) ivars.size(); rvar2 < nvar2; rvar2++)
      {
        int ivar2 = ivars[rvar2];

        // Loop on the second sample
        int nech2s = (int) index1[rvar2].size();
        for (int rech2 = 0; rech2 < nech2s; rech2++)
        {
          // Perform calculation only in upper triangle of the Symmetric Matrix
          if (icol >= irow)
          {
            int iech2 = index1[rvar2][rech2];
            p2.setIech(iech2);
            db1->getSampleAsSPInPlace(p2);

            // Modify the covariance (if non stationary)
            updateCovByPoints(1, iech1, 1, iech2);

            /* Loop on the dimension of the space */
            double value = eval(p1, p2, ivar1, ivar2, mode);
            mat.setValue(irow, icol, value);
          }
          icol++;
        }
      }
      irow++;
    }
  }

  // Update the matrix due to presence of Variance of Measurement Error
  _updateCovMatrixSymmetricVerr(db1, &mat, ivars, index1);

  return mat;
}

/****************************************************************************/
/*!
 **  Establish the covariance matrix between two Dbs where samples are selected by ranks
 **  The output is stored in a Sparse Matrix
 **
 ** \return Sparse matrix containing the covariance matrix
 **
 ** \param[in]  db1     First Db
 ** \param[in]  db2     Second Db
 ** \param[in]  ivar0   Rank of the first variable (-1: all variables)
 ** \param[in]  jvar0   Rank of the second variable (-1: all variables)
 ** \param[in]  nbgh1   Array giving ranks of selected samples (optional)
 ** \param[in]  nbgh2   Array giving ranks of selected samples (optional)
 ** \param[in]  mode    CovCalcMode structure
 ** \param[in]  eps     Tolerance for discarding a covariance value
 **
 ** \remarks The covariance matrix (returned) must be freed by calling routine
 ** \remarks The covariance matrix is established for the first variable
 ** \remarks and returned as a covariance
 ** \remarks As the ranks are used, no test is performed on any selection
 ** \remarks but only ranks positive or null are considered
 **
 *****************************************************************************/
MatrixSparse* ACov::evalCovMatrixSparse(const Db *db1,
                                        const Db *db2,
                                        int ivar0,
                                        int jvar0,
                                        const VectorInt &nbgh1,
                                        const VectorInt &nbgh2,
                                        const CovCalcMode *mode,
                                        double eps)
{
  MatrixSparse* mat = nullptr;
  if (db2 == nullptr) db2 = db1;
  if (db1 == nullptr || db2 == nullptr) return mat;
  bool flagSameDb = (db1 == db2);
  VectorInt ivars = _getActiveVariables(ivar0);
  if (ivars.empty()) return mat;
  VectorInt jvars = _getActiveVariables(jvar0);
  if (jvars.empty()) return mat;

  // Play the non-stationarity (if needed)

  manage(db1, db2);
  

  // Create the sets of Vector of valid sample indices per variable (not masked and defined)
  VectorVectorInt index1 = db1->getMultipleRanksActive(ivars, nbgh1, true, flagSameDb);
  VectorVectorInt index2 = db2->getMultipleRanksActive(jvars, nbgh2, true, flagSameDb);

  // Evaluate the matrix of sills
  int nvar1 = (int) ivars.size();
  int nvar2 = (int) jvars.size();
  MatrixRectangular mat0(nvar1, nvar2);
  for (int ivar = 0; ivar < nvar1; ivar++)
  {
    int ivar1 = ivars[ivar];
    for (int jvar = 0; jvar < nvar2; jvar++)
    {
      int jvar2 = jvars[jvar];
      double value = eval0(ivar1, jvar2, mode);
      mat0.setValue(ivar1, jvar2, value);
    }
  }

  // Constitute the triplet
  NF_Triplet NF_T;

  // Define the two space points
  SpacePoint p1(getSpace());
  SpacePoint p2(getSpace());

  // Loop on the first variable
  int irow = 0;
  for (int ivar = 0; ivar < nvar1; ivar++)
  {
    int ivar1 = ivars[ivar];

    // Loop on the first sample
    int nech1s = (int) index1[ivar].size();
    for (int jech1 = 0; jech1 < nech1s; jech1++)
    {
      int iech1 = index1[ivar][jech1];
      p1.setIech(iech1);
      db1->getSampleAsSPInPlace(p1);

      // Loop on the second variable
      int icol = 0;
      for (int jvar = 0; jvar < nvar2; jvar++)
      {
        int jvar2 = jvars[jvar];

        // Loop on the second sample
        int nech2s = (int) index2[jvar].size();
        for (int jech2 = 0; jech2 < nech2s; jech2++)
        {
          int iech2 = index2[jvar][jech2];
          p2.setIech(iech2);
          db2->getSampleAsSPInPlace(p2);

          // Modify the covariance (if non stationary)
          updateCovByPoints(1, iech1, 2, iech2);

          /* Loop on the dimension of the space */
          double value = eval(p1, p2, ivar1, jvar2, mode);

          if (ABS(value) >= eps * mat0.getValue(ivar1, jvar2))
            NF_T.add(irow, icol, value);
          icol++;
        }
      }
      irow++;
    }
  }

  // Convert from triplet to sparse matrix

  mat = MatrixSparse::createFromTriplet(NF_T);

  // Update the matrix due to presence of Variance of Measurement Error
  if (flagSameDb)
    _updateCovMatrixSymmetricVerr(db1, mat, ivars, index1);
  return mat;
}


/**
 * Variance of Extension of a set of points and the block
 * @param db      Reference Data Base
 * @param ext     Vector giving the extensions of the target block
 * @param ndisc   Vector giving the discretization
 * @param angles  Vector for the rotation angles of the block (optional)
 * @param x0      Optional origin of the Block
 * @param ivar    Rank of the first variable
 * @param jvar    Rank of the second variable
 * @return
 */
double ACov::extensionVariance(const Db* db,
                               const VectorDouble& ext,
                               const VectorInt&    ndisc,
                               const VectorDouble& angles,
                               const VectorDouble& x0,
                               int ivar,
                               int jvar) const
{
  double sigmaE = TEST;
  DbGrid* dbgrid = _discretizeBlock(ext, ndisc, angles, x0);
  if (dbgrid != nullptr)
  {
    double GxV = evalAverageDbToDb(db, dbgrid, ivar, jvar);
    double Gxx = evalAverageDbToDb(db, db, ivar, jvar);
    double GVV = evalAverageDbToDb(dbgrid, dbgrid, ivar, jvar);
    sigmaE = -2. * GxV + Gxx + GVV;
  }
  delete dbgrid;

  return sigmaE;
}

/**
 * Calculate the Sampling Density Variance
 * @param db      Set of data points
 * @param ext     Block extension
 * @param ndisc   Discretization
 * @param angles  Optional rotation angles for the Block
 * @param x0      Optional origin of the block
 * @param ivar    Rank of the first variable
 * @param jvar    Rank of the second variable
 * @return
 */
double ACov::samplingDensityVariance(const Db *db,
                                     const VectorDouble &ext,
                                     const VectorInt &ndisc,
                                     const VectorDouble &angles,
                                     const VectorDouble &x0,
                                     int ivar,
                                     int jvar) const
{
  double sigmaE = extensionVariance(db, ext, ndisc, angles, x0, ivar, jvar);
  double maille = _getVolume(ext);
  return sigmaE * maille;
}

/**
 * Calculate the Specific Volume
 * @param db     Set of data points
 * @param mean   Value of the Mean
 * @param ext    Target Block extension
 * @param ndisc  Vector of discretization
 * @param angles Optional rotation angle for block
 * @param x0     Optional origin of the Block
 * @param ivar   Rank of the first variable
 * @param jvar   Rank of the second variable
 * @return
 */
double ACov::specificVolume(const Db *db,
                            double mean,
                            const VectorDouble &ext,
                            const VectorInt &ndisc,
                            const VectorDouble &angles,
                            const VectorDouble &x0,
                            int ivar,
                            int jvar) const
{
  if (FFFF(mean) || mean <= 0.)
  {
    messerr("Argument 'mean'  must be defined and positive");
    return TEST;
  }
  return samplingDensityVariance(db, ext, ndisc, angles, x0, ivar, jvar)
      / (mean * mean);
}

/**
 * Calculate the Coefficient of Variation
 * @param db     Set of data points
 * @param volume Specific production volume
 * @param mean   Value of the Mean
 * @param ext    Target Block extension
 * @param ndisc  Vector of discretization
 * @param angles Optional rotation angle for block
 * @param x0     Optional origin of the Block
 * @param ivar   Rank of the first variable
 * @param jvar   Rank of the second variable
 * @return
 */
double ACov::coefficientOfVariation(const Db *db,
                                    double volume,
                                    double mean,
                                    const VectorDouble &ext,
                                    const VectorInt &ndisc,
                                    const VectorDouble &angles,
                                    const VectorDouble &x0,
                                    int ivar,
                                    int jvar) const
{
  if (FFFF(mean) || mean <= 0.)
  {
    messerr("Argument 'mean'  must be defined and positive");
    return TEST;
  }
  if (FFFF(volume) || volume <= 0.)
  {
    messerr("Argument 'volume'  must be defined and positive");
    return TEST;
  }
  double V0 = specificVolume(db, mean, ext, ndisc, angles, x0, ivar, jvar);
  return sqrt(V0 / volume);
}

/**
 * Derive the Specific volume for a given CoV
 * @param db     Set of data points
 * @param cov    Target Coefficient of Variation
 * @param mean   Value of the Mean
 * @param ext    Target Block extension
 * @param ndisc  Vector of discretization
 * @param angles Optional rotation angle for block
 * @param x0     Optional origin of the Block
 * @param ivar   Rank of the first variable
 * @param jvar   Rank of the second variable
 * @return
 */
double ACov::specificVolumeFromCoV(Db *db,
                                   double cov,
                                   double mean,
                                   const VectorDouble &ext,
                                   const VectorInt &ndisc,
                                   const VectorDouble &angles,
                                   const VectorDouble &x0,
                                   int ivar,
                                   int jvar) const
{
  double V0 = specificVolume(db, mean, ext, ndisc, angles, x0, ivar, jvar);
  return V0 / (cov * cov);
}

double ACov::_getVolume(const VectorDouble& ext) const
{
  double maille = 1.;
  int ndim = getNDim();
  for (int idim = 0; idim < ndim; idim++) maille *= ext[idim];
  return maille;
}

