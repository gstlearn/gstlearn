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
#include "Covariances/CovContext.hpp"
#include "Enum/ECalcMember.hpp"
#include "Enum/EKrigOpt.hpp"
#include "Matrix/MatrixSquareGeneral.hpp"
#include "Matrix/MatrixRectangular.hpp"
#include "Matrix/MatrixSparse.hpp"
#include "Matrix/MatrixSquareSymmetric.hpp"
#include "Matrix/NF_Triplet.hpp"
#include "Db/Db.hpp"
#include "Geometry/GeometryHelper.hpp"
#include "Db/DbGrid.hpp"
#include "Basic/AException.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/VectorNumT.hpp"
#include "Basic/VectorHelper.hpp"
#include "Basic/Law.hpp"
#include "Space/ASpace.hpp"
#include "Space/SpacePoint.hpp"
#include "geoslib_define.h"
#include "Covariances/NoStatArray.hpp"
#include "Covariances/NoStatFunctional.hpp"
#include "Variogram/Vario.hpp"
#include "Estimation/KrigOpt.hpp"

#include <vector>
#include <math.h>

ACov::ACov(const CovContext& ctxt)
  : ASpaceObject(ctxt.getSpace())
  , _ctxt(ctxt)
  , _optimEnabled(false)
  , _optimPreProcessedData(false)
  , _p1As()
  , _p2As()
  , _p2A(ctxt.getSpace())
  , _pAux(ctxt.getSpace())
  , _tabNoStat(nullptr)
{
  createNoStatTab();
}

ACov::ACov(const ACov& r)
  : ASpaceObject(r)
  , _ctxt(r._ctxt)
  , _optimEnabled(r._optimEnabled)
  , _optimPreProcessedData(r._optimPreProcessedData)
  , _p1As(r._p1As)
  , _p2As(r._p2As)
  , _p2A(r._p2A)
  , _pAux(r._pAux)
  , _pw1(r._pw1)
  , _pw2(r._pw2)
  , _tabNoStat(r._tabNoStat == nullptr ? nullptr : r._tabNoStat->clone())
{
}

ACov& ACov::operator=(const ACov& r)
{
  if (this != &r)
  {
    ASpaceObject::operator=(r);
    _ctxt                  = r._ctxt;
    _optimEnabled          = r._optimEnabled;
    _optimPreProcessedData = r._optimPreProcessedData;
    _p1As                  = r._p1As;
    _p2As                  = r._p2As;
    _p2A                   = r._p2A;
    _pAux                  = r._pAux;
    _pw1                   = r._pw1;
    _pw2                   = r._pw2;
    _tabNoStat             = r._tabNoStat->clone();
  }
  return *this;
}

ACov::~ACov()
{
  delete _tabNoStat;
}

double ACov::evalCov(const SpacePoint& p1,
                     const SpacePoint& p2,
                     int ivar,
                     int jvar,
                     const CovCalcMode* mode) const
{
  return _eval(p1, p2, ivar, jvar, mode);
}
void ACov::optimizationPostProcess() const
{
  _p1As.clear();
  _p2As.clear();
  _optimizationPostProcess();
  _optimPreProcessedData = false;
}

void ACov::optimizationPreProcessForData(const Db* db1) const
{
  // Do not proceed to optimization when this is forbidden
  // if (! _optimEnabled) return;

  // Do not proceed if already done
  if (_optimPreProcessedData) return;

  std::vector<SpacePoint> ps;

  // Add projected samples from db1 (optional)
  if (db1 != nullptr)
  {
    db1->getSamplesAsSP(ps, getSpace());
    _optimizationPreProcess(1, ps);
  }

  _optimPreProcessedData = true;
}

void ACov::_optimizationPreProcessForTarget(const Db* db2, const VectorInt& nbgh2) const
{
  std::vector<SpacePoint> ps;

  // Add projected samples from db2 (optional)
  if (nbgh2.empty())
    db2->getSamplesAsSP(ps, getSpace(), true);
  else
    db2->getSamplesFromNbghAsSP(ps, nbgh2);
  _optimizationPreProcess(2, ps);
}

void ACov::optimizationPreProcess(int mode, const std::vector<SpacePoint>& ps) const
{
  _optimizationPreProcess(mode, ps);
}

SpacePoint& ACov::optimizationLoadInPlace(int iech, int mode, int rank) const
{
  return _optimizationLoadInPlace(iech, mode, rank);
}

void ACov::optimizationSetTarget(SpacePoint& pt) const
{
  _optimizationSetTarget(pt);
}

void ACov::_optimizationSetTarget(SpacePoint& pt) const
{
  DECLARE_UNUSED(pt);
}

/**
 * @brief Preprocess the input Data.
 * By default, this method only copies the SpacePoints
 * In the Anisoptric version, the samples are projected along with the Covariance
 *
 * @param mode 1 for p1As, 2 for p2As
 * @param ps Set of SpacePoints to be copied
 */
void ACov::_optimizationPreProcess(int mode, const std::vector<SpacePoint>& ps) const
{
  if (mode == 1)
  {
    _p1As.clear();
    for (const auto& e: ps)
      _p1As.push_back(e);
  }
  else
  {
    _p2As.clear();
    for (const auto& e: ps)
      _p2As.push_back(e);
  }
}

void ACov::createNoStatTab()
{
  delete _tabNoStat;
  _tabNoStat = _createNoStatTab();
}

void ACov::attachNoStatDb(const Db* db)
{
  _tabNoStat->setDbNoStatRef(db);
}

VectorDouble ACov::informCoords(const VectorVectorDouble& coords,
                                const EConsElem& econs,
                                int iv1,
                                int iv2) const
{
  VectorDouble result(coords[0].size(), getValue(econs, iv1, iv2));
  _tabNoStat->informCoords(coords, econs, iv1, iv2, result);
  return result;
}

bool ACov::checkAndManageNoStatDb(const Db*& db, const String& namecol)
{
  if (_tabNoStat->getDbNoStatRef() == nullptr && db == nullptr)
  {
    messerr("You have to define a Db (with attachNoStatDb or by specifying a Db here)");
    return false;
  }
  setNoStatDbIfNecessary(db);

  if (db->getUID(namecol) < 0)
  {
    messerr("You have to specified a name of a column of the reference Db");
    return false;
  }
  return true;
}

TabNoStat* ACov::_createNoStatTab()
{
  return new TabNoStat();
}

bool ACov::_checkDims(int idim, int jdim) const
{
  int ndim = getNDim();
  if ((idim > ndim) || (jdim > ndim))
  {
    messerr("Your model is only in dimension %d.", ndim);
    return false;
  }
  return true;
}

void ACov::_optimizationPostProcess() const
{
}

MatrixSquareSymmetric ACov::eval0Mat(const CovCalcMode* mode) const
{
  int nvar = getNVar();
  MatrixSquareSymmetric mat(nvar);
  mat.fill(0.);
  eval0CovMatBiPointInPlace(mat, mode);
  return mat;
}

/**
 * Calculate the Matrix of covariance for zero distance
 * @param mat   Covariance matrix (Dimension: nvar * nvar)
 * @param mode  Calculation Options
 *
 * @remarks: Matrix 'mat' should be dimensioned and initialized beforehand
 */
void ACov::eval0CovMatBiPointInPlace(MatrixSquareSymmetric& mat,
                                     const CovCalcMode* mode) const
{
  int nvar = getNVar();
  for (int ivar = 0; ivar < nvar; ivar++)
    for (int jvar = 0; jvar <= ivar; jvar++)
    {
      double value = eval0(ivar, jvar, mode);
      mat.addValue(ivar, jvar, value);
    }
}

MatrixSquareSymmetric ACov::evalCovMat0(const Db* db, int iech, const KrigOpt& krigopt) const
{
  MatrixSquareSymmetric mat;

  int error = evalCovMat0InPlace(mat, db, iech, krigopt);
  return (error == 0) ? mat : MatrixSquareSymmetric();
}

int ACov::evalCovMat0InPlace(MatrixSquareSymmetric& mat,
                             const Db* db,
                             int iech,
                             const KrigOpt& krigopt) const
{
  const EKrigOpt& calcul = krigopt.getCalcul();
  if (calcul == EKrigOpt::DGM)
  {
    messerr("This method is not designed for DGM Krigopt option");
    return 1;
  }

  int nvar = getNVar();
  mat.resize(nvar, nvar);
  mat.fill(0.);

  if (calcul == EKrigOpt::DRIFT) return 1;
  bool flagBlock         = calcul == EKrigOpt::BLOCK;
  const CovCalcMode mode = CovCalcMode(ECalcMember::VAR);
  bool isNoStatLocal     = isNoStat();

  SpacePoint p0(getSpace());
  db->getSampleAsSPInPlace(p0, iech);

  // Modify the covariance (if non stationary)
  if (isNoStatLocal)
    updateCovByPoints(2, iech, 2, iech);

  if (!flagBlock)
    eval0CovMatBiPointInPlace(mat, &mode);
  else
  {
    if (krigopt.isFlagCell()) krigopt.blockDiscretize(0, true);

    VectorVectorDouble d1 = krigopt.getDisc1VVD();
    VectorVectorDouble d2 = krigopt.getDisc2VVD();

    for (int ivar = 0; ivar < nvar; ivar++)
      for (int jvar = 0; jvar <= ivar; jvar++)
        mat.setValue(ivar, jvar, evalAverageIncrToIncr(d1, d2, ivar, jvar, &mode));
  };

  // In case of combined R.H.S., modify the output matrix
  if (krigopt.isMatLC()) mat = mat.compress0MatLC(*krigopt.getMatLC());

  return 0;
}

VectorDouble ACov::eval(const std::vector<SpacePoint>& vec_p1,
                        const std::vector<SpacePoint>& vec_p2,
                        int ivar,
                        int jvar,
                        const CovCalcMode* mode) const
{
  VectorDouble vec;
  if (vec_p1.size() != vec_p2.size())
    my_throw("Error: 'p1' and 'p2' should have same dimension");
  for (int i = 0, n = static_cast<int>(vec_p1.size()); i < n; i++)
    vec.push_back(evalCov(vec_p1[i], vec_p2[i], ivar, jvar, mode)); // pure virtual method
  return vec;
}

double ACov::eval0(int ivar,
                   int jvar,
                   const CovCalcMode* mode) const
{
  SpacePoint p1(getSpace()->getOrigin(), -1);
  return evalCov(p1, p1, ivar, jvar, mode); // pure virtual method
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
  ASpaceSharedPtr space = getSpace();
  SpacePoint p1(space);
  SpacePoint p2(space);

  if (dir.empty())
  {
    VectorDouble vec(getNDim(), 0.);
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

  return evalCov(p1, p2, ivar, jvar, mode); // pure virtual method
}

double ACov::evalIvarIpasIncr(const VectorDouble& dincr,
                              int ivar,
                              int jvar,
                              const CovCalcMode* mode) const
{
  // Define the point in the ACov space (center will be checked)
  SpacePoint p1(VectorDouble(_space->getNDim()), -1, getSpace());
  SpacePoint p2(VectorDouble(_space->getNDim()), -1, getSpace());
  p2.move(dincr);
  return evalCov(p1, p2, ivar, jvar, mode); // pure virtual method
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
VectorDouble ACov::evalIvarNlag(const VectorDouble& vec_step,
                                const VectorDouble& dir,
                                int ivar,
                                int jvar,
                                const CovCalcMode* mode) const
{
  VectorDouble vec;
  for (int i = 0, n = static_cast<int>(vec_step.size()); i < n; i++)
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
  int nvar = getNVar();
  MatrixSquareGeneral mat(nvar);
  for (int ivar = 0; ivar < nvar; ivar++)
    for (int jvar = 0; jvar < nvar; jvar++)
      mat.setValue(ivar, jvar, evalIvarIpas(step, dir, ivar, jvar, mode));
  return mat;
}

MatrixSquareGeneral ACov::evalNvarIpasIncr(const VectorDouble& dincr,
                                           const CovCalcMode* mode) const
{
  int nvar = getNVar();
  MatrixSquareGeneral mat(nvar);
  for (int ivar = 0; ivar < nvar; ivar++)
    for (int jvar = 0; jvar < nvar; jvar++)
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
  VectorDouble dir    = getUnitaryVector();
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
VectorDouble ACov::evalIsoIvarNlag(const VectorDouble& vec_step,
                                   int ivar,
                                   int jvar,
                                   const CovCalcMode* mode) const
{
  VectorDouble vec;
  VectorDouble dir = getUnitaryVector();
  for (const auto& h: vec_step)
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
  int nvar         = getNVar();
  VectorDouble dir = getUnitaryVector();
  MatrixSquareGeneral mat(nvar);
  for (int ivar = 0; ivar < nvar; ivar++)
    for (int jvar = 0; jvar < nvar; jvar++)
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
  for (int iech1 = 0; iech1 < db1->getNSample(); iech1++)
  {
    if (!db1->isActive(iech1)) continue;
    double w1 = db1->getWeight(iech1);
    if (isZero(w1)) continue;
    SpacePoint p1(db1->getSampleCoordinates(iech1));

    /* Loop on the second sample */

    for (int iech2 = 0; iech2 < db2->getNSample(); iech2++)
    {
      if (!db2->isActive(iech2)) continue;
      double w2 = db2->getWeight(iech2);
      if (isZero(w2)) continue;
      VectorDouble coord2 = db2->getSampleCoordinates(iech2);
      if (eps > 0)
      {
        for (int idim = 0, ndim = getNDim(); idim < ndim; idim++)
          coord2[idim] += eps * law_uniform(-0.5, 0.5);
      }
      SpacePoint p2(coord2);

      /* Loop on the dimension of the space */

      total += w1 * w2 * evalCov(p1, p2, ivar, jvar, mode);
      norme += w1 * w2;
    }
  }

  // Scaling
  if (!isZero(norme)) total /= norme;

  if (eps > 0. && seed > 0)
    law_set_random_seed(memo);

  return total;
}

double ACov::evalAverageIncrToIncr(const VectorVectorDouble& d1,
                                   const VectorVectorDouble& d2,
                                   int ivar,
                                   int jvar,
                                   const CovCalcMode* mode) const
{
  int nincr1 = (int)d1.size();
  int nincr2 = (int)d2.size();

  /* Loop on the first sample */

  double total = 0.;
  for (int incr1 = 0; incr1 < nincr1; incr1++)
  {
    SpacePoint p1(d1[incr1], -1, getSpace());

    /* Loop on the second sample */

    for (int incr2 = 0; incr2 < nincr2; incr2++)
    {
      SpacePoint p2(d2[incr2], -1, getSpace());
      total += evalCov(p1, p2, ivar, jvar, mode);
    }
  }

  // Scaling
  total /= (double)(nincr1 * nincr2);

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

  for (int iech2 = 0; iech2 < db2->getNSample(); iech2++)
  {
    if (!db2->isActive(iech2)) continue;
    double w2 = db2->getWeight(iech2);
    if (isZero(w2)) continue;
    SpacePoint p2(db2->getSampleCoordinates(iech2), iech2, getSpace());

    /* Loop on the dimension of the space */

    total += w2 * evalCov(p1, p2, ivar, jvar, mode);
    norme += w2;
  }

  // Scaling
  if (isZero(norme)) total /= norme;

  return total;
}

void ACov::evalPointToDbAsSP(VectorDouble& values,
                             const std::vector<SpacePoint>& p1s,
                             const SpacePoint& p2,
                             int ivar,
                             int jvar,
                             const CovCalcMode* mode) const
{
  int nech1 = (int)p1s.size();
  if (nech1 != (int)values.size()) values.resize(nech1);

  /* Loop on the second sample */

  for (int iech1 = 0; iech1 < nech1; iech1++)
  {
    const SpacePoint& p1 = p1s[iech1];
    p1.setIech(iech1);
    values[iech1] = evalCov(p1, p2, ivar, jvar, mode);
  }
}

/**
 * Calculate the Covariance vector between a Point and all the samples
 * of a Db, for a pair of variables
 * @param values Array of returned values (possible resized)
 * @param p1   Point location
 * @param db2  Pointer to the second Db
 * @param ivar Rank of the first variables
 * @param jvar Rank of the second variable
 * @param useSel When TRUE, the returned vector is reduced to active samples
 *               Otherwise, returns TEST for masked samples
 * @param nbgh2 Vector of indices of active samples in db2 (optional)
 * @param mode CovCalcMode structure
 */
void ACov::evalPointToDb(VectorDouble& values,
                         const SpacePoint& p1,
                         const Db* db2,
                         int ivar,
                         int jvar,
                         bool useSel,
                         const VectorInt& nbgh2,
                         const CovCalcMode* mode) const
{
  SpacePoint p2(getSpace());

  VectorInt index2 = db2->getRanksActive(nbgh2, jvar, useSel);
  int nech2        = (int)index2.size();
  if (nech2 != (int)values.size()) values.resize(nech2);

  int irow = 0;
  for (const auto i: index2.getVector())
  {
    int iabs2 = (nbgh2.empty()) ? i : nbgh2[i];
    db2->getSampleAsSPInPlace(p2, iabs2);
    values[irow++] = evalCov(p1, p2, ivar, jvar, mode);
  }
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
  if (ndim != (int)ext.size())
  {
    messerr("Block Extension (%d) should have same dimension as the Model %d)",
            (int)ext.size(), ndim);
    return TEST;
  }
  if (ndim != (int)ndisc.size())
  {
    messerr("Discretization (%d) should have same dimension as the Model (%d)",
            (int)ndisc.size(), ndim);
    return TEST;
  }

  DbGrid* dbgrid = _discretizeBlock(ext, ndisc, angles);
  if (dbgrid == nullptr) return TEST;
  Db* db = _discretizeBlockRandom(dbgrid);
  if (db == nullptr) return TEST;

  double total = evalAverageDbToDb(dbgrid, db, ivar, jvar, 0., 0, mode);
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
  if (ndim != (int)ext.size())
  {
    messerr("Block Extension (%d) should have same dimension as the Model %d)",
            (int)ext.size(), ndim);
    return TEST;
  }
  if (ndim != (int)ndisc.size())
  {
    messerr("Discretization (%d) should have same dimension as the Model (%d)",
            (int)ndisc.size(), ndim);
    return TEST;
  }
  if (ndim != (int)shift.size())
  {
    messerr("Shift (%d) should have the same dimension as the Model (%d)",
            (int)shift.size(), ndim);
    return TEST;
  }

  DbGrid* dbgrid1 = _discretizeBlock(ext, ndisc, angles);
  if (dbgrid1 == nullptr) return TEST;
  DbGrid* dbgrid2 = _discretizeBlock(ext, ndisc, angles, shift);
  if (dbgrid2 == nullptr) return TEST;

  double total = evalAverageDbToDb(dbgrid1, dbgrid2, ivar, jvar, 0., 0, mode);
  delete dbgrid1;
  delete dbgrid2;
  return total;
}

MatrixSquareGeneral ACov::evalCvvM(const VectorDouble& ext,
                                   const VectorInt& ndisc,
                                   const VectorDouble& angles,
                                   const CovCalcMode* mode) const
{
  int nvar = getNVar();
  MatrixSquareGeneral mat(nvar);
  for (int ivar = 0; ivar < nvar; ivar++)
    for (int jvar = 0; jvar < nvar; jvar++)
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
  if (ndim != (int)ext.size())
  {
    messerr("Block Extension (%d) should have same dimension as the Model %d)",
            (int)ext.size(), ndim);
    return TEST;
  }
  if (ndim != (int)ndisc.size())
  {
    messerr("Discretization (%d) should have same dimension as the Model (%d)",
            (int)ndisc.size(), ndim);
    return TEST;
  }

  double total   = TEST;
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
  if (ndim != (int)ext.size())
  {
    messerr("Block Extension (%d) should have same dimension as the Model %d)",
            (int)ext.size(), ndim);
    return TEST;
  }
  if (ndim != (int)ndisc.size())
  {
    messerr("Discretization (%d) should have same dimension as the Model (%d)",
            (int)ndisc.size(), ndim);
    return TEST;
  }

  double total   = TEST;
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
  int nvar = getNVar();
  MatrixSquareGeneral mat(nvar);
  for (int ivar = 0; ivar < nvar; ivar++)
    for (int jvar = 0; jvar < nvar; jvar++)
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
  int ndim           = getNDim();
  VectorDouble x0loc = x0;
  if (x0loc.empty() || ndim != (int)x0loc.size())
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
  int ndim           = getNDim();
  int nech           = dbgrid->getNSample();
  Db* db             = Db::createFromSamples(nech);
  VectorString names = generateMultipleNames("x", ndim);
  law_set_random_seed(seed);

  for (int idim = 0; idim < ndim; idim++)
  {
    double taille    = dbgrid->getDX(idim);
    VectorDouble vec = dbgrid->getOneCoordinate(idim, false);
    for (int i = 0; i < (int)vec.size(); i++)
      vec[i] += taille * law_uniform(-0.5, 0.5);
    db->addColumns(vec, names[idim], ELoc::X, idim);
  }
  return db;
}

VectorInt ACov::_getActiveVariables(int ivar0) const
{
  int nvar = getNVar();

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
 ** \param[in]  cleanOptim When True, clean optimization internal when ended
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
MatrixRectangular ACov::evalCovMat(const Db* db1,
                                   const Db* db2,
                                   int ivar0,
                                   int jvar0,
                                   const VectorInt& nbgh1,
                                   const VectorInt& nbgh2,
                                   const CovCalcMode* mode,
                                   bool cleanOptim) const
{
  MatrixRectangular mat;

  int error = evalCovMatInPlace(mat, db1, db2, ivar0, jvar0, nbgh1, nbgh2, mode, cleanOptim);
  return (error) == 0 ? mat : MatrixRectangular();
}

int ACov::evalCovMatInPlace(MatrixRectangular& mat,
                            const Db* db1,
                            const Db* db2,
                            int ivar0,
                            int jvar0,
                            const VectorInt& nbgh1,
                            const VectorInt& nbgh2,
                            const CovCalcMode* mode,
                            bool cleanOptim) const
{
  // Preliminary checks
  if (db2 == nullptr) db2 = db1;
  if (db1 == nullptr || db2 == nullptr) return 1;
  VectorInt ivars = _getActiveVariables(ivar0);
  if (ivars.empty()) return 1;
  VectorInt jvars = _getActiveVariables(jvar0);
  if (jvars.empty()) return 1;

  VectorVectorInt index1 = db1->getSampleRanks(ivars, nbgh1);
  VectorVectorInt index2 = db2->getSampleRanks(jvars, nbgh2);

  return evalCovMatInPlaceFromIdx(mat, db1, db2, index1, index2, nbgh2, mode, cleanOptim);
}

int ACov::evalCovMatInPlaceFromIdx(MatrixRectangular& mat,
                                   const Db* db1,
                                   const Db* db2,
                                   const VectorVectorInt& index1,
                                   const VectorVectorInt& index2,
                                   const VectorInt& nbgh2,
                                   const CovCalcMode* mode,
                                   bool cleanOptim) const
{
  // Prepare Non-stationarity (if needed)
  manage(db1, db2);

  // Prepare Optimization for covariance calculation (if not forbidden or already done)
  optimizationPreProcessForData(db1);
  _optimizationPreProcessForTarget(db2, nbgh2);

  // Creating the matrix
  int nvar1 = (int)index1.size();
  int nvar2 = (int)index2.size();
  int neq1  = VH::count(index1);
  int neq2  = VH::count(index2);
  if (neq1 <= 0 || neq2 <= 0)
  {
    messerr("The returned matrix has no valid sample and no valid variable");
    return 1;
  }
  mat.resize(neq1, neq2);

  // Define the two space points
  bool isNoStatLocal = isNoStat();

  // Loop on Data
  int icol = 0;
  for (int ivar2 = 0; ivar2 < nvar2; ivar2++)
  {
    const VectorInt& index2i = index2[ivar2];
    const int* ptr2          = index2i.data();
    for (int irel2 = 0, n2 = (int)index2i.size(); irel2 < n2; irel2++)
    {
      int iabs2      = *ptr2++;
      SpacePoint& p2 = optimizationLoadInPlace(irel2, 2, 2);

      int irow = 0;
      for (int ivar1 = 0; ivar1 < nvar1; ivar1++)
      {
        const VectorInt& index1i = index1[ivar1];
        for (const auto iabs1: index1i.getVector())
        {
          SpacePoint& p1 = optimizationLoadInPlace(iabs1, 1, 1);

          // Modify the covariance (if non stationary)
          if (isNoStatLocal)
            updateCovByPoints(1, iabs1, 2, iabs2);

          // Calculate the covariance between two points
          double value = evalCov(p1, p2, ivar1, ivar2, mode);
          mat.setValue(irow, icol, value);

          irow++;
        }
      }
      icol++;
    }
  }

  if (cleanOptim) optimizationPostProcess();

  return 0;
}

/**
 * @brief Returns the references in 'pt' and set the local pointer '_pw2'
 *
 * @param iech Relative index of the target sample (within 'pXAs')
 * @param mode 1 for _p1As, 2 for _p2As and 3 for _p2A
 * @param rank 1 for the first point and 2 for the second
 */
SpacePoint& ACov::_optimizationLoadInPlace(int iech,
                                           int mode,
                                           int rank) const
{
  if (mode == 1)
  {
    if (rank == 1)
      _pw1 = &_p1As[iech];
    else
      _pw2 = &_p1As[iech];
    return _p1As[iech];
  }

  if (mode == 2)
  {
    if (rank == 1)
      _pw1 = &_p2As[iech];
    else
      _pw2 = &_p2As[iech];
    return _p2As[iech];
  }

  if (rank == 1)
    _pw1 = &_p2A;
  else
    _pw2 = &_p2A;
  return _p2A;
}

/****************************************************************************/
/*!
 **  Establish covariance matrix between one Db and one sample of a Target Db
 **
 ** \return Dense matrix containing the covariance matrix
 **
 ** \param[in]  mat Matrix (possibly resized)
 ** \param[in]  db1   First Db
 ** \param[in]  db2   Second Db
 ** \param[in]  index1 Vector of vector indices of active samples in db1
 ** \param[in]  iech2 Sample rank within db2
 ** \param[in]  krigopt KrigOpt structure
 ** \param[in]  cleanOptim When True, clean optimization internal when ended
 **
 ** \remarks If a Db does not contain any Z-variable defined, the covariance
 ** \remarks cannot treat possible heterotopy and therefore uses all samples
 **
 ** \remarks The returned matrix if dimension to nrows * 1 where
 ** \remarks each 'nrows' is the number of active samples
 ** \remarks by the number of samples where the variable is defined
 **
 ** \note 'dbin' and 'dbout' cannot be made 'const' as they can be updated
 ** \note due to the presence of 'nostat'
 **
 *****************************************************************************/
int ACov::evalCovMatRHSInPlaceFromIdx(MatrixRectangular& mat,
                                      const Db* db1,
                                      const Db* db2,
                                      const VectorVectorInt& index1,
                                      int iech2,
                                      const KrigOpt& krigopt,
                                      bool cleanOptim) const
{
  // Preliminary checks
  if (db1 == nullptr || db2 == nullptr) return 1;
  const EKrigOpt& calcul = krigopt.getCalcul();
  if (calcul == EKrigOpt::DGM)
  {
    messerr("This method is not designed for DGM Krigopt option");
    return 1;
  }
  VectorInt ivars = VH::sequence(getNVar());
  if (ivars.empty()) return 1;

  // Create the sets of Vectors of valid sample indices per variable
  // (not masked and defined)
  VectorInt nbgh2        = {iech2};
  VectorVectorInt index2 = db2->getSampleRanks(ivars, nbgh2, true, false, false);

  // Creating the matrix
  int neq1 = VH::count(index1);
  int neq2 = VH::count(index2);
  if (neq1 <= 0 || neq2 <= 0)
  {
    messerr("The returned matrix has no valid sample and no valid variable");
    return 1;
  }

  // Dimension the returned matrix
  mat.resize(neq1, neq2);
  mat.fill(0.);

  // Prepare covariance optimization and non-stationarity (except for Mean)
  if (calcul != EKrigOpt::DRIFT)
  {
    // Play the non-stationarity (if needed)
    manage(db1, db2);

    // Prepare the Optimization for covariance calculation
    optimizationPreProcessForData(db1);
    _optimizationPreProcessForTarget(db2, nbgh2);
  }

  // Dispatch according to the type of estimation
  if (calcul == EKrigOpt::POINT)
    _evalCovMatRHSInPlacePoint(mat, index1, index2, krigopt);
  else if (calcul == EKrigOpt::BLOCK)
    _evalCovMatRHSInPlaceBlock(mat, db2, index1, index2, krigopt);
  else if (calcul == EKrigOpt::DRIFT)
  {
    // No calculation needed for Large scale drift estimation
  }
  else
  {
    messerr("Unknown Calculation type");
    return 1;
  }

  // In case of combined R.H.S., modify the output matrix
  if (krigopt.isMatLC()) mat = mat.compressMatLC(*krigopt.getMatLC());

  if (cleanOptim) optimizationPostProcess();
  return 0;
}

/****************************************************************************/
/*!
 **  Establish covariance matrix between one Db and one sample of a Target Db
 **
 ** \return Dense matrix containing the covariance matrix
 **
 ** \param[in]  mat Matrix (possibly resized)
 ** \param[in]  db1   First Db
 ** \param[in]  db2   Second Db
 ** \param[in]  index1 Vector of vector indices of active samples in db1
 ** \param[in]  iech2 Sample rank within db2
 ** \param[in]  krigopt KrigOpt structure
 ** \param[in]  cleanOptim When True, clean optimization internal when ended
 **
 ** \remarks If a Db does not contain any Z-variable defined, the covariance
 ** \remarks cannot treat possible heterotopy and therefore uses all samples
 **
 ** \remarks The returned matrix if dimension to nrows * 1 where
 ** \remarks each 'nrows' is the number of active samples
 ** \remarks by the number of samples where the variable is defined
 **
 ** \note 'dbin' and 'dbout' cannot be made 'const' as they can be updated
 ** \note due to the presence of 'nostat'
 **
 *****************************************************************************/
int ACov::evalCovVecRHSInPlace(vect vect,
                               const Db* db2,
                               const VectorInt& index1,
                               int iech2,
                               const KrigOpt& krigopt,
                               SpacePoint& pin,
                               SpacePoint& pout,
                               VectorDouble& tabwork,
                               double lambda) const
{
  db2->getSampleAsSPInPlace(pin, iech2);
  for (int i = 0; i < (int)vect.size(); i++)
    vect[i] = 0.;
  return addEvalCovVecRHSInPlace(vect, index1, iech2, krigopt, pin, pout, tabwork, lambda);
}

int ACov::addEvalCovVecRHSInPlace(vect vect,
                                  const VectorInt& index1,
                                  int iech2,
                                  const KrigOpt& krigopt,
                                  SpacePoint& pin,
                                  SpacePoint& pout,
                                  VectorDouble& tabwork,
                                  double lambda) const
{ 
  DECLARE_UNUSED(pout,tabwork);
  optimizationSetTarget(pin);
  bool flagNoStat = isNoStat();
  const CovCalcMode& mode = krigopt.getMode();
  const int* inds = index1.data();
  for (int i = 0; i < (int)vect.size();i++)
  {
    if (flagNoStat)
        updateCovByPoints(1, *inds, 2, iech2);
    SpacePoint& p1 = optimizationLoadInPlace(*inds++, 1, 1);

    vect[i] += lambda * evalCov(p1, pin, 0, 0, &mode);
  }
  return 0;
}

int ACov::_evalCovMatRHSInPlaceBlock(MatrixRectangular& mat,
                                     const Db* db2,
                                     const VectorVectorInt& index1,
                                     const VectorVectorInt& index2,
                                     const KrigOpt& krigopt) const
{
  int ndisc = krigopt.getNDisc();
  SpacePoint p2(getSpace());
  SpacePoint p2aux(getSpace());
  const CovCalcMode& mode = krigopt.getMode();
  bool isNoStatLocal      = isNoStat();

  // Loop on Data

  int icol = 0;
  for (int ivar2 = 0, nvar2 = (int)index2.size(); ivar2 < nvar2; ivar2++)
  {
    const VectorInt& index2i = index2[ivar2];
    const int* ptr2          = index2i.data();
    for (int irel2 = 0, n2 = (int)index2i.size(); irel2 < n2; irel2++)
    {
      int iabs2 = *ptr2++;

      // Identify the center of the block
      db2->getSampleAsSPInPlace(p2, iabs2);

      // Discretize the block if adapted to the cell dimensions
      if (krigopt.isFlagCell()) krigopt.blockDiscretize(p2.getIech());

      // Loop on the discretization points
      for (int idisc = 0; idisc < ndisc; idisc++)
      {
        // Move the target to the discretization point
        p2aux = p2;
        p2aux.move(krigopt.getDisc1VD(idisc));
        optimizationSetTarget(p2aux);

        _loopOnData(mat, p2aux, ivar2, iabs2, icol,
                    true, isNoStatLocal && idisc == 0,
                    index1, mode);
      }

      // Scale the matrix
      _scaleOnData(mat, icol, ndisc);
      icol++;
    }
  }

  return 0;
}

int ACov::_evalCovMatRHSInPlacePoint(MatrixRectangular& mat,
                                     const VectorVectorInt& index1,
                                     const VectorVectorInt& index2,
                                     const KrigOpt& krigopt) const
{
  // Local shortcuts for parameters
  const CovCalcMode& mode = krigopt.getMode();
  bool isNoStatLocal      = isNoStat();

  // Loop on Target
  int icol = 0;
  for (int ivar2 = 0, nvar2 = (int)index2.size(); ivar2 < nvar2; ivar2++)
  {
    const VectorInt& index2i = index2[ivar2];
    const int* ptr2          = index2i.data();
    for (int irel2 = 0, n2 = (int)index2i.size(); irel2 < n2; irel2++)
    {
      int iabs2 = *ptr2++;

      SpacePoint& p2 = optimizationLoadInPlace(irel2, 2, 2);

      _loopOnData(mat, p2, ivar2, iabs2, icol,
                  false, isNoStatLocal, index1, mode);
      icol++;
    }
  }
  return 0;
}

void ACov::_loopOnData(MatrixRectangular& mat,
                       const SpacePoint& p2,
                       int ivar2,
                       int iabs2,
                       int icol,
                       bool flagUpdate,
                       bool flagNoStat,
                       const VectorVectorInt& index1,
                       const CovCalcMode& mode) const
{
  double value;

  int irow = 0;
  for (int ivar1 = 0, nvar1 = (int)index1.size(); ivar1 < nvar1; ivar1++)
  {
    const VectorInt& index1i = index1[ivar1];
    for (const auto iabs1: index1i.getVector())
    {
      SpacePoint& p1 = optimizationLoadInPlace(iabs1, 1, 1);

      // Modify the covariance (if non stationary)
      if (flagNoStat)
        updateCovByPoints(1, iabs1, 2, iabs2);

      // Loop on the discretization points
      value = evalCov(p1, p2, ivar1, ivar2, &mode);
      if (flagUpdate)
        mat.updValue(irow, icol, EOperator::ADD, value);
      else
        mat.setValue(irow, icol, value);
      irow++;
    }
  }
}

void ACov::_scaleOnData(MatrixRectangular& mat, int icol, int ndisc)
{
  int nrows = mat.getNRows();
  double value;

  for (int irow = 0; irow < nrows; irow++)
  {
    value = mat.getValue(irow, icol);
    mat.setValue(irow, icol, value / ndisc);
  }
}

void ACov::_updateCovMatrixSymmetricForVerr(const Db* db1,
                                            AMatrix* mat,
                                            const VectorVectorInt& index1)
{
  // Check if the correction can take place at all
  if (!db1->hasLocVariable(ELoc::V)) return;

  // Loop on Data
  int irow = 0;
  for (int ivar1 = 0, nvar1 = (int)index1.size(); ivar1 < nvar1; ivar1++)
  {
    int icolVerr             = db1->getColIdxByLocator(ELoc::V, ivar1);
    const VectorInt& index1i = index1[ivar1];
    for (const auto iabs1: index1i.getVector())
    {
      double verr = 0.;

      // Update the Diagonal due to the presence of Variance of Measurement Error
      if (icolVerr >= 0)
        verr = db1->getValueByColIdx(iabs1, icolVerr);

      // Update the Covariance matrix
      if (verr > 0) mat->updValue(irow, irow, EOperator::ADD, verr);
      irow++;
    }
  }
}

void ACov::load(const SpacePoint& p, bool case1) const
{
  _load(p, case1);
}

void ACov::_load(const SpacePoint& p, bool option) const
{
  DECLARE_UNUSED(p);
  DECLARE_UNUSED(option);
}

/****************************************************************************/
/*!
 **  Establish the covariance matrix within a Db
 **  Takes into account selection and heterotopy
 **  This method takes advantage of calculating covariance between
 **  a Db and itself
 **
 ** \return Dense matrix containing the covariance matrix
 **
 ** \param[in]  db1   First Db
 ** \param[in]  ivar0 Rank of the first variable (-1 for all variables)
 ** \param[in]  nbgh1 Vector of indices of active samples in db1 (optional)
 ** \param[in]  mode  CovCalcMode structure
 ** \param[in]  cleanOptim When True, clean optimization internal arrays at end
 **
 ** \remarks If a Db does not contain any Z-variable defined, the covariance
 ** \remarks cannot treat possible heterotopy and therefore uses all samples
 **
 ** \remarks The returned matrix if dimension to nrows * ncols where
 ** \remarks each term is the product of the number of active samples
 ** \remarks by the number of samples where the variable is defined
 **
 *****************************************************************************/
MatrixSquareSymmetric ACov::evalCovMatSym(const Db* db1,
                                          const VectorInt& nbgh1,
                                          int ivar0,
                                          const CovCalcMode* mode,
                                          bool cleanOptim) const
{
  MatrixSquareSymmetric mat;

  int error = evalCovMatSymInPlace(mat, db1, nbgh1, ivar0, mode, cleanOptim);
  return (error == 0) ? mat : MatrixSquareSymmetric();
}

int ACov::evalCovMatSymInPlace(MatrixSquareSymmetric& mat,
                               const Db* db1,
                               const VectorInt& nbgh1,
                               int ivar0,
                               const CovCalcMode* mode,
                               bool cleanOptim) const
{
  // Preliminary checks
  if (db1 == nullptr) return 1;
  VectorInt ivars = _getActiveVariables(ivar0);
  if (ivars.empty()) return 1;

  // Create the sets of Vector of valid sample indices per variable (not masked and defined)
  VectorVectorInt index1 = db1->getSampleRanks(ivars, nbgh1, true, true, true);

  return evalCovMatSymInPlaceFromIdx(mat, db1, index1, mode, cleanOptim);
}

int ACov::evalCovMatSymInPlaceFromIdx(MatrixSquareSymmetric& mat,
                                      const Db* db1,
                                      const VectorVectorInt& index1,
                                      const CovCalcMode* mode,
                                      bool cleanOptim) const
{
  // Creating the matrix
  int neq1 = VH::count(index1);
  if (neq1 <= 0)
  {
    messerr("The returned matrix has no valid sample and no valid variable");
    return 1;
  }
  mat.resize(neq1, neq1);

  // Prepare Non-stationarity (if needed)
  manage(db1, nullptr);

  // Prepare Optimization for covariance calculation (if not forbidden or already done)
  optimizationPreProcessForData(db1);

  // Define the two space points
  bool isNoStatLocal = isNoStat();

  // Loop on Data
  double value;
  int icol = 0;
  for (int ivar2 = 0, nvar2 = (int)index1.size(); ivar2 < nvar2; ivar2++)
  {
    const VectorInt& index2i = index1[ivar2];
    for (const auto iabs2: index2i.getVector())
    {
      SpacePoint& p2 = optimizationLoadInPlace(iabs2, 1, 2);

      int irow = 0;
      for (int ivar1 = 0, nvar1 = (int)index1.size(); ivar1 < nvar1; ivar1++)
      {
        const VectorInt& index1i = index1[ivar1];
        for (const auto iabs1: index1i.getVector())
        {
          if (icol >= irow)
          {
            SpacePoint& p1 = optimizationLoadInPlace(iabs1, 1, 1);

            // Modify the covariance (if non stationary)
            if (isNoStatLocal)
              updateCovByPoints(1, iabs1, 1, iabs2);

            // Calculate the covariance between two points
            if (iabs1 == iabs2)
              value = eval0(ivar1, ivar2, mode);
            else
              value = evalCov(p1, p2, ivar1, ivar2, mode);
            mat.setValue(irow, icol, value);
          }
          irow++;
        }
      }
      icol++;
    }
  }

  // Update the matrix due to presence of Variance of Measurement Error
  _updateCovMatrixSymmetricForVerr(db1, &mat, index1);

  if (cleanOptim) optimizationPostProcess();
  return 0;
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
 ** \param[in]  cleanOptim When True, clean optimization internal when ended
 ** \param[in]  eps     Tolerance for discarding a covariance value
 **
 ** \remarks The covariance matrix (returned) must be freed by calling routine
 ** \remarks The covariance matrix is established for the first variable
 ** \remarks and returned as a covariance
 ** \remarks As the ranks are used, no test is performed on any selection
 ** \remarks but only ranks positive or null are considered
 **
 *****************************************************************************/
MatrixSparse* ACov::evalCovMatSparse(const Db* db1,
                                     const Db* db2,
                                     int ivar0,
                                     int jvar0,
                                     const VectorInt& nbgh1,
                                     const VectorInt& nbgh2,
                                     const CovCalcMode* mode,
                                     bool cleanOptim,
                                     double eps) const
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
  bool isNoStatLocal = isNoStat();
  manage(db1, db2);

  // Prepare the Optimization for covariance calculation
  optimizationPreProcessForData(db1);
  _optimizationPreProcessForTarget(db2, nbgh2);

  // Create the sets of Vector of valid sample indices per variable (not masked and defined)
  VectorVectorInt index1 = db1->getSampleRanks(ivars, nbgh1, true, true, flagSameDb);
  VectorVectorInt index2 = db2->getSampleRanks(jvars, nbgh2, true, true, flagSameDb);

  // Evaluate the matrix of sills
  int nvar1 = (int)ivars.size();
  int nvar2 = (int)jvars.size();
  MatrixRectangular mat0(nvar1, nvar2);
  for (int ivar = 0; ivar < nvar1; ivar++)
  {
    int ivar1 = ivars[ivar];
    for (int jvar = 0; jvar < nvar2; jvar++)
    {
      int jvar2    = jvars[jvar];
      double value = eval0(ivar1, jvar2, mode);
      mat0.setValue(ivar1, jvar2, value);
    }
  }

  // Constitute the triplet
  NF_Triplet NF_T;

  // Loop on Data
  int icol = 0;
  for (int ivar2 = 0, nvar2 = (int)index2.size(); ivar2 < nvar2; ivar2++)
  {
    const VectorInt& index2i = index2[ivar2];
    const int* ptr2          = index2i.data();
    for (int irel2 = 0, n2 = (int)index2i.size(); irel2 < n2; irel2++)
    {
      int iabs2      = *ptr2++;
      SpacePoint& p2 = optimizationLoadInPlace(irel2, 2, 2);

      int irow = 0;
      for (int ivar1 = 0, nvar1 = (int)index1.size(); ivar1 < nvar1; ivar1++)
      {
        const VectorInt& index1i = index1[ivar1];
        for (const auto iabs1: index1i.getVector())
        {
          SpacePoint& p1 = optimizationLoadInPlace(iabs1, 1, 1);

          // Modify the covariance (if non stationary)
          if (isNoStatLocal)
            updateCovByPoints(1, iabs1, 2, iabs2);

          /* Loop on the dimension of the space */
          double value = evalCov(p1, p2, ivar1, ivar2, mode);

          if (ABS(value) >= eps * mat0.getValue(ivar1, ivar2))
            NF_T.add(irow, icol, value);

          irow++;
        }
      }
      icol++;
    }
  }

  // Convert from triplet to sparse matrix

  mat = MatrixSparse::createFromTriplet(NF_T);

  // Update the matrix due to presence of Variance of Measurement Error
  if (flagSameDb)
    _updateCovMatrixSymmetricForVerr(db1, mat, index1);

  if (cleanOptim) optimizationPostProcess();
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
                               const VectorInt& ndisc,
                               const VectorDouble& angles,
                               const VectorDouble& x0,
                               int ivar,
                               int jvar) const
{
  double sigmaE  = TEST;
  DbGrid* dbgrid = _discretizeBlock(ext, ndisc, angles, x0);
  if (dbgrid != nullptr)
  {
    double GxV = evalAverageDbToDb(db, dbgrid, ivar, jvar);
    double Gxx = evalAverageDbToDb(db, db, ivar, jvar);
    double GVV = evalAverageDbToDb(dbgrid, dbgrid, ivar, jvar);
    sigmaE     = -2. * GxV + Gxx + GVV;
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
double ACov::samplingDensityVariance(const Db* db,
                                     const VectorDouble& ext,
                                     const VectorInt& ndisc,
                                     const VectorDouble& angles,
                                     const VectorDouble& x0,
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
double ACov::specificVolume(const Db* db,
                            double mean,
                            const VectorDouble& ext,
                            const VectorInt& ndisc,
                            const VectorDouble& angles,
                            const VectorDouble& x0,
                            int ivar,
                            int jvar) const
{
  if (FFFF(mean) || mean <= 0.)
  {
    messerr("Argument 'mean'  must be defined and positive");
    return TEST;
  }
  return samplingDensityVariance(db, ext, ndisc, angles, x0, ivar, jvar) / (mean * mean);
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
double ACov::coefficientOfVariation(const Db* db,
                                    double volume,
                                    double mean,
                                    const VectorDouble& ext,
                                    const VectorInt& ndisc,
                                    const VectorDouble& angles,
                                    const VectorDouble& x0,
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
double ACov::specificVolumeFromCoV(Db* db,
                                   double cov,
                                   double mean,
                                   const VectorDouble& ext,
                                   const VectorInt& ndisc,
                                   const VectorDouble& angles,
                                   const VectorDouble& x0,
                                   int ivar,
                                   int jvar) const
{
  double V0 = specificVolume(db, mean, ext, ndisc, angles, x0, ivar, jvar);
  return V0 / (cov * cov);
}

double ACov::_getVolume(const VectorDouble& ext) const
{
  double maille = 1.;
  int ndim      = getNDim();
  for (int idim = 0; idim < ndim; idim++) maille *= ext[idim];
  return maille;
}

/////////////  Functions to attach no stat information on various supports ////////
void ACov::informMeshByMesh(const AMesh* amesh) const
{
  _tabNoStat->informMeshByMesh(amesh);
}
void ACov::informMeshByApex(const AMesh* amesh) const
{
  _tabNoStat->informMeshByApex(amesh);
}
void ACov::informDbIn(const Db* dbin) const
{
  _tabNoStat->informDbIn(dbin);
}
void ACov::informDbOut(const Db* dbout) const
{
  _tabNoStat->informDbOut(dbout);
}

void ACov::setNoStatDbIfNecessary(const Db*& db)
{
  if (_tabNoStat->getDbNoStatRef() == nullptr)
    attachNoStatDb(db);
  if (db == nullptr)
    db = _tabNoStat->getDbNoStatRef();
}

void ACov::makeStationary()
{
  _tabNoStat->clear();
}
int ACov::makeElemNoStat(const EConsElem& econs, int iv1, int iv2, const AFunctional* func, const Db* db, const String& namecol)
{
  std::shared_ptr<ANoStat> ns;
  if (func == nullptr)
  {
    if (!checkAndManageNoStatDb(db, namecol)) return 1;
    ns = std::shared_ptr<ANoStat>(new NoStatArray(db, namecol));
  }
  else
  {
    ns = std::unique_ptr<ANoStat>(new NoStatFunctional(func));
  }
  return _tabNoStat->addElem(ns, econs, iv1, iv2);
}

/*****************************************************************************/
/*!
 **  Returns the covariance for an increment
 **  This is the generic internal function
 **  It can be called for stationary or non-stationary case
 **
 ** \param[in]  covint       Internal structure for non-stationarityAddress for the next term after the drift
 **                          or NULL (for stationary case)
 ** \param[in]  mode         CovCalcMode structure
 ** \param[in]  weight       Multiplicative weight
 ** \param[in]  d1           Distance vector
 **
 *****************************************************************************/
double ACov::evaluateOneGeneric(const CovInternal* covint,
                                const VectorDouble& d1,
                                double weight,
                                const CovCalcMode* mode) const
{
  // Load the non-stationary parameters if needed

  if (covint != nullptr)
    updateCovByPoints(covint->getIcas1(), covint->getIech1(),
                      covint->getIcas2(), covint->getIech2());

  // Return the (weighted) Model value

  return (weight * evalIvarIpas(1, d1, 0, 0, mode));
}

/*****************************************************************************/
/*!
 **  Returns the standard deviation at a given increment for a given model
 **  between two samples of two Dbs
 **
 ** \param[in]  db1         First Db
 ** \param[in]  iech1       Rank in the first Db
 ** \param[in]  db2         Second Db
 ** \param[in]  iech2       Rank in the second Db
 ** \param[in]  verbose     Verbose flag
 ** \param[in]  factor      Multiplicative factor for standard deviation
 ** \param[in]  mode        CovCalcMode structure
 **
 *****************************************************************************/
double ACov::calculateStDev(Db* db1,
                            int iech1,
                            Db* db2,
                            int iech2,
                            bool verbose,
                            double factor,
                            const CovCalcMode* mode) const
{

  /* Covariance at origin */

  int ndim = db1->getNDim();
  VectorDouble dd(ndim, 0.);
  double c00 = evaluateOneGeneric(nullptr, dd, 1., mode);

  /* Covariance at increment */

  if (db1->getDistanceVecInPlace(iech1, iech2, dd, db2) != 0) return TEST;
  double cov   = evaluateOneGeneric(nullptr, dd, 1., mode);
  double stdev = factor * sqrt(c00 - cov);

  if (verbose)
  {
    message("Db1(%d) - Db2(%d)", iech1 + 1, iech2 + 1);
    message(" - Incr=");
    for (int idim = 0; idim < ndim; idim++)
      message(" %lf", dd[idim]);
    message(" - c(0)=%lf cov=%lf stdev=%lf\n", c00, cov, stdev);
  }
  return stdev;
}

/****************************************************************************/
/*!
 **  Evaluate the model on a Db
 **
 ** \param[in]  db         Db structure
 ** \param[in]  ivar       Rank of the first variable
 ** \param[in]  jvar       Rank of the second variable
 ** \param[in]  mode       CovCalcMode structure
 **
 *****************************************************************************/
VectorDouble ACov::evaluateFromDb(Db* db,
                                  int ivar,
                                  int jvar,
                                  const CovCalcMode* mode) const
{
  if ((int)getNDim() != db->getNDim())
  {
    messerr("Dimension of the Db (%d) does not match dimension of the Model (%d)",
            db->getNDim(), getNDim());
    return VectorDouble();
  }
  int ndim = getNDim();
  int nvar = getNVar();
  int nech = db->getNSample();

  /* Core allocation */

  VectorDouble d1(ndim, 0.);
  MatrixSquareGeneral covtab(nvar);
  VectorDouble gg(nech, TEST);

  /* Loop on the lags */

  for (int iech = 0; iech < nech; iech++)
  {
    if (!db->isActive(iech)) continue;
    db->getCoordinatesInPlace(d1, iech);
    evaluateMatInPlace(nullptr, d1, covtab, true, 1., mode);
    gg[iech] = covtab.getValue(ivar, jvar);
  }
  return gg;
}

/*****************************************************************************/
/*!
 **  Returns the covariances for an increment
 **  This is the generic internal function
 **  It can be called for stationary or non-stationary case
 **
 ** \param[in]  covint       Internal structure for non-stationarityAddress for the next term after the drift
 **                          or NULL (for stationary case)
 ** \param[in]  mode         CovCalcMode structure
 ** \param[in]  flag_init    Initialize the array beforehand
 ** \param[in]  weight       Multiplicative weight
 ** \param[in]  d1           Distance vector
 ** \param[out] covtab       Covariance array
 **
 *****************************************************************************/
void ACov::evaluateMatInPlace(const CovInternal* covint,
                              const VectorDouble& d1,
                              MatrixSquareGeneral& covtab,
                              bool flag_init,
                              double weight,
                              const CovCalcMode* mode) const
{
  // Load the non-stationary parameters if needed
  if (isNoStat() && covint != nullptr)
    updateCovByPoints(covint->getIcas1(), covint->getIech1(),
                      covint->getIcas2(), covint->getIech2());

  // Evaluate the Model
  MatrixSquareGeneral mat = evalNvarIpas(1., d1, mode);

  int nvar = getNVar();
  for (int ivar = 0; ivar < nvar; ivar++)
    for (int jvar = 0; jvar < nvar; jvar++)
    {
      double value = weight * mat.getValue(ivar, jvar);
      if (flag_init)
        covtab.setValue(ivar, jvar, value);
      else
        covtab.updValue(ivar, jvar, EOperator::ADD, value);
    }
}

/****************************************************************************/
/*!
 **  Calculate the variogram map from a Model
 **  (presented as Variogram, not Covariance)
 **
 ** \return  Error return code
 **
 ** \param[in]  dbgrid      Grid structure
 ** \param[in]  namconv     Naming convention
 **
 *****************************************************************************/
int ACov::buildVmapOnDbGrid(DbGrid* dbgrid, const NamingConvention& namconv) const
{
  if (dbgrid == nullptr) return 1;

  /* Initializations */

  int ndim = dbgrid->getNDim();
  int nvar = dbgrid->getNLoc(ELoc::Z);
  int nv2  = nvar * (nvar + 1) / 2;

  /* Create the variables in the Variogram Map file */

  int iptr = dbgrid->addColumnsByConstant(nv2, 0.);
  if (iptr < 0) return 1;

  /* Loop on the grid nodes */

  CovCalcMode mode(ECalcMember::LHS);
  mode.setAsVario(true);
  VectorInt center = dbgrid->getCenterIndices();
  VectorDouble dincr(ndim);
  VectorInt indices(ndim);
  MatrixSquareGeneral mat;
  for (int iech = 0; iech < dbgrid->getNSample(); iech++)
  {
    if (!dbgrid->isActive(iech)) continue;
    dbgrid->rankToIndice(iech, indices);

    for (int idim = 0; idim < ndim; idim++)
      dincr[idim] = (indices[idim] - center[idim]) * dbgrid->getDX(idim);

    // Evaluate the variogram map
    mat = evalNvarIpasIncr(dincr, &mode);

    int ecr = 0;
    for (int ivar = 0; ivar < nvar; ivar++)
      for (int jvar = 0; jvar <= ivar; jvar++, ecr++)
        dbgrid->setArray(iech, iptr + ecr, mat.getValue(ivar, jvar));
  }

  /* Set the error return code */

  namconv.setNamesAndLocators(dbgrid, iptr, "Model", nv2);
  return 0;
}

/****************************************************************************/
/*!
 **  Calculate the value of the model for a set of distances
 **
 ** \return  The model value
 **
 ** \param[in]  ivar       Rank of the first variable
 ** \param[in]  jvar       Rank of the second variable
 ** \param[in]  mode       CovCalcMode structure
 ** \param[in]  codir      Array giving the direction coefficients (optional)
 ** \param[in]  hh         Vector of increments
 **
 *****************************************************************************/
double ACov::evaluateOneIncr(double hh,
                             const VectorDouble& codir,
                             int ivar,
                             int jvar,
                             const CovCalcMode* mode) const
{
  int ndim = getNDim();
  int nvar = getNVar();

  /* Core allocation */

  VectorDouble d1(ndim);
  MatrixSquareGeneral covtab(nvar);

  /* Normalize the direction vector codir */

  /* Get the normalized direction vector */

  VectorDouble codir_loc = codir;
  if (codir_loc.empty())
  {
    (void)GH::rotationGetDirectionDefault(ndim, codir_loc);
  }
  else
  {
    VH::normalizeCodir(ndim, codir_loc);
  }

  for (int idim = 0; idim < ndim; idim++)
    d1[idim] = hh * codir_loc[idim];
  evaluateMatInPlace(nullptr, d1, covtab, true, 1., mode);
  return covtab.getValue(ivar, jvar);
}

/****************************************************************************/
/*!
 **  Calculate the value of the model for a set of distances
 **
 ** \return  Array containing the model values
 **
 ** \param[in]  ivar       Rank of the first variable
 ** \param[in]  jvar       Rank of the second variable
 ** \param[in]  codir      Array giving the direction coefficients (optional)
 ** \param[in]  h          Vector of increments
 ** \param[in]  mode       CovCalcMode structure
 ** \param[in]  covint     Non-stationary parameters
 **
 *****************************************************************************/
VectorDouble ACov::sample(const VectorDouble& h,
                          const VectorDouble& codir,
                          int ivar,
                          int jvar,
                          const CovCalcMode* mode,
                          const CovInternal* covint) const
{
  int nh   = (int)h.size();
  int ndim = getNDim();
  int nvar = getNVar();

  /* Core allocation */

  VectorDouble d1(ndim);
  MatrixSquareGeneral covtab(nvar);

  /* Get the normalized direction vector */

  VectorDouble codir_loc = codir;
  if (codir_loc.empty())
  {
    (void)GH::rotationGetDirectionDefault(ndim, codir_loc);
  }
  else
  {
    VH::normalizeCodir(ndim, codir_loc);
  }

  /* Loop on the lags */

  VectorDouble g(nh);
  for (int ih = 0; ih < nh; ih++)
  {
    double hh = h[ih];
    for (int idim = 0; idim < ndim; idim++)
      d1[idim] = hh * codir_loc[idim];
    evaluateMatInPlace(covint, d1, covtab, true, 1., mode);
    g[ih] = covtab.getValue(ivar, jvar);
  }
  return g;
}

/**
 * Returns the value of the normalized covariance (by the variance/covariance value)
 * for a given pair of variables
 * @param hh    Vector of distances
 * @param ivar  Rank of the first variable
 * @param jvar  Rank of the second variable
 * @param codir Direction coefficients
 * @param mode  CovCalcMode structure
 * @return
 */
VectorDouble ACov::sampleUnitary(const VectorDouble& hh,
                                 int ivar,
                                 int jvar,
                                 VectorDouble codir,
                                 const CovCalcMode* mode) const
{
  if (ivar < 0 || ivar >= getNVar()) return VectorDouble();
  if (jvar < 0 || jvar >= getNVar()) return VectorDouble();
  if (ivar == jvar) return VectorDouble();
  int ndim = getNDim();
  if (codir.empty())
  {
    (void)GH::rotationGetDirectionDefault(ndim, codir);
  }
  int nh = (int)hh.size();

  double c00      = eval0(ivar, ivar, mode);
  double c11      = eval0(jvar, jvar, mode);
  c00             = sqrt(c00 * c11);
  VectorDouble gg = sample(hh, codir, ivar, jvar, mode);

  for (int i = 0; i < nh; i++)
    gg[i] /= c00;

  return gg;
}

VectorDouble ACov::envelop(const VectorDouble& hh,
                           int ivar,
                           int jvar,
                           int isign,
                           VectorDouble codir,
                           const CovCalcMode* mode) const
{
  if (ivar < 0 || ivar >= getNVar()) return VectorDouble();
  if (jvar < 0 || jvar >= getNVar()) return VectorDouble();
  if (ivar == jvar) return VectorDouble();
  if (isign != -1 && isign != 1) return VectorDouble();
  int ndim = getNDim();
  if (codir.empty())
  {
    (void)GH::rotationGetDirectionDefault(ndim, codir);
  }
  int nh = (int)hh.size();
  VectorDouble gg(nh);
  VectorDouble g1 = sample(hh, codir, ivar, ivar, mode);
  VectorDouble g2 = sample(hh, codir, jvar, jvar, mode);

  for (int i = 0; i < nh; i++)
    gg[i] = isign * sqrt(abs(g1[i] * g2[i]));

  return gg;
}

/**
 * Evaluate the Goodness-of_fit of the Model on the Experimental Variogram
 * It is expressed as the average departure between Model and Variogram
 * scaled to the variance.
 * As this variance may be poorly calculated (< gmax / 5), it may be replaced
 * by the largest value (gmax) divided by 2 (highly non_stationary cases).
 * @param vario Experimental variogram
 * @param verbose Verbose flag

 * @return Value for the Goodness-of_fit (as percentage of the total sill)
 */
double ACov::gofToVario(const Vario* vario, bool verbose) const
{
  int nvar = getNVar();
  int ndir = vario->getNDir();

  double total = 0.;

  // Loop on the pair of variables

  CovCalcMode mode(ECalcMember::LHS);
  mode.setAsVario(true);
  for (int ivar = 0; ivar < nvar; ivar++)
    for (int jvar = 0; jvar < nvar; jvar++)
    {
      double varij  = vario->getVar(ivar, jvar);
      double varmax = vario->getGmax(ivar, jvar);
      // Modify the normalization as variance seems not consistent
      if (ABS(varij) < varmax / 5)
      {
        if (verbose)
          messerr("Variance seems erroneous. It is replaced by Gmax / 2.");
        varij = varmax / 2.;
      }

      // Loop on the variogram directions

      double totdir = 0.;
      for (int idir = 0; idir < ndir; idir++)
      {

        // Read information from Experimental Variogram

        VectorDouble codir = vario->getCodirs(idir);
        VectorDouble sw    = vario->getSwVec(idir, ivar, jvar);
        VectorDouble hh    = vario->getHhVec(idir, ivar, jvar);
        VectorDouble gexp  = vario->getGgVec(idir, ivar, jvar);

        // Evaluate the Model

        int nlag          = (int)gexp.size();
        VectorDouble gmod = sample(hh, codir, ivar, jvar, &mode);

        // Evaluate the score

        double totpas = 0;
        double scale  = 0.;
        for (int ilag = 0; ilag < nlag; ilag++)
        {
          if (sw[ilag] <= 0 || hh[ilag] <= 0.) continue;
          double ecart = sw[ilag] * ABS(gexp[ilag] - gmod[ilag]) / hh[ilag];
          totpas += ecart;
          scale += sw[ilag] / hh[ilag];
        }
        totpas = totpas / scale;
        totdir += totpas;
      }
      totdir /= (double)ndir;
      totdir /= varij;
      total += ABS(totdir);
    }
  total = 100. * total / (double)(nvar * nvar);
  return total;
}

/**
 * Printout of statement concerning the Quality of the GOF
 * @param gof        Value of the Gof
 * @param byValue    true: display GOF value; false: print its quality level
 * @param thresholds Vector giving the Quality thresholds
 */
void ACov::gofDisplay(double gof, bool byValue, const VectorDouble& thresholds)
{
  message("Goodness-of-fit (as a percentage of the variance)");
  if (byValue)
  {
    message(" = %5.2lf\n", gof);
    return;
  }
  int nclass = (int)thresholds.size();
  for (int iclass = 0; iclass < nclass; iclass++)
  {
    if (gof < thresholds[iclass])
    {
      message(" corresponds to level #%d (1 for very good)\n", iclass + 1);
      return;
    }
  }
}

void ACov::setContext(const CovContext& ctxt)
{
  _ctxt = ctxt;
  _setContext(ctxt);
}
