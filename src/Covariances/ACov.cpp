/******************************************************************************/
/* COPYRIGHT ARMINES, ALL RIGHTS RESERVED                                     */
/*                                                                            */
/* THE CONTENT OF THIS WORK CONTAINS CONFIDENTIAL AND PROPRIETARY             */
/* INFORMATION OF ARMINES. ANY DUPLICATION, MODIFICATION,                     */
/* DISTRIBUTION, OR DISCLOSURE IN ANY FORM, IN WHOLE, OR IN PART, IS STRICTLY */
/* PROHIBITED WITHOUT THE PRIOR EXPRESS WRITTEN PERMISSION OF ARMINES         */
/*                                                                            */
/* TAG_SOURCE_CG                                                              */
/******************************************************************************/
#include "Covariances/ACov.hpp"
#include "Matrix/MatrixSquareGeneral.hpp"
#include "Db/Db.hpp"
#include "Db/DbGrid.hpp"
#include "Basic/AException.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/Vector.hpp"
#include "Basic/Law.hpp"

ACov::ACov(const ASpace* space)
: ASpaceObject(space)
{
}

ACov::ACov(const ACov &r)
: ASpaceObject(r)
{
}

ACov& ACov::operator=(const ACov &r)
{
  if (this != &r)
  {
    ASpaceObject::operator=(r);
  }
  return *this;
}

ACov::~ACov()
{
}

MatrixSquareGeneral ACov::eval0Nvar(const CovCalcMode& mode) const
{
  int nvar = getNVariables();
  MatrixSquareGeneral mat(nvar);
  for (int ivar=0; ivar<nvar; ivar++)
    for (int jvar=0; jvar<nvar; jvar++)
      mat.setValue(ivar, jvar, eval0(ivar, jvar, mode)); // pure virtual method
  return mat;
}

VectorDouble ACov::eval(int ivar,
                        int jvar,
                        const std::vector<SpacePoint>& vec_p1,
                        const std::vector<SpacePoint>& vec_p2,
                        const CovCalcMode& mode) const
{
  VectorDouble vec;
  if (vec_p1.size() != vec_p2.size())
    my_throw ("Error: 'p1' and 'p2' should have same dimension");
  for (int i=0, n=static_cast<int> (vec_p1.size()); i < n; i++)
    vec.push_back(eval(ivar, jvar, vec_p1[i], vec_p2[i], mode)); // pure virtual method
  return vec;
}

MatrixSquareGeneral ACov::eval(const SpacePoint& p1,
                               const SpacePoint& p2,
                               const CovCalcMode& mode) const
{
  int nvar = getNVariables();
  MatrixSquareGeneral mat(nvar);
  for (int ivar=0; ivar<nvar; ivar++)
    for (int jvar=0; jvar<nvar; jvar++)
      mat.setValue(ivar, jvar, eval(ivar, jvar, p1, p2, mode)); // pure virtual method
  return mat;
}

/**
 * Covariance from a given point (center) in a given direction (dir *step)
 * @param ivar   Rank of the first variable
 * @param jvar   Rank of the second variable
 * @param step   Step value
 * @param dir    Direction definition
 * @param center Coordinates of the center
 * @param mode   CovCalcMode structure
 * @return
 */
double ACov::evalIvarIpas(int ivar,
                          int jvar,
                          double step,
                          const VectorDouble& dir,
                          const VectorDouble& center,
                          const CovCalcMode& mode) const
{
  // Define the point in the ACov space (center will be checked)
  SpacePoint p1(center,getSpace());
  SpacePoint p2(center,getSpace());
  VectorDouble dirloc(dir);
  if (dirloc.empty())
  {
    dirloc.resize(getNDim(),0.);
    dirloc[0] = 1.;
  }
  VectorDouble vec(dirloc);
  ut_vector_multiply_inplace(vec, step);
  p2.move(vec);
  return eval(ivar, jvar, p1, p2, mode); // pure virtual method
}

/**
 * Covariance vector from a given point (center) in a given direction (dir * steps)
 * for a pair of variables and a set of steps
 * @param ivar      Rank of the first variable
 * @param jvar      Rank of the second variable
 * @param vec_step  Vector of step values
 * @param dir       Direction definition
 * @param center    Coordinates of the Center
 * @param mode      CovCalcMode structure
 * @return
 */
VectorDouble ACov::evalIvarNpas(int ivar,
                                int jvar,
                                const VectorDouble& vec_step,
                                const VectorDouble& dir,
                                const VectorDouble& center,
                                const CovCalcMode& mode) const
{
  VectorDouble vec;
  for (int i=0, n=static_cast<int> (vec_step.size()); i < n; i++)
    vec.push_back(evalIvarIpas(ivar, jvar, vec_step[i], dir, center, mode));
  return vec;
}

/**
 * Covariance Matrix from a given point (center) in a given direction (dir * step)
 * for a set of variables and a given step
 * @param step   Step value
 * @param dir    Direction definition
 * @param center Coordinates of the center
 * @param mode   CovCalcMode structure
 * @return
 */
MatrixSquareGeneral ACov::evalNvarIpas(double step,
                                       const VectorDouble& dir,
                                       const VectorDouble& center,
                                       const CovCalcMode& mode) const
{
  int nvar = getNVariables();
  MatrixSquareGeneral mat(nvar);
  for (int ivar=0; ivar<nvar; ivar++)
    for (int jvar=0; jvar<nvar; jvar++)
      mat.setValue(ivar, jvar, evalIvarIpas(ivar, jvar, step, dir, center, mode));
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
double ACov::evalIsoIvarIpas(int ivar,
                             int jvar,
                             double step,
                             const CovCalcMode& mode) const
{
  /// TODO : Not true whatever the space
  VectorDouble center = getOrigin();
  VectorDouble dir = getUnitaryVector();
  return evalIvarIpas(ivar, jvar, step, dir, center, mode);
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
VectorDouble ACov::evalIsoIvarNpas(int ivar,
                                   int jvar,
                                   const VectorDouble& vec_step,
                                   const CovCalcMode& mode) const
{
  VectorDouble vec;
  /// TODO : Not true whatever the space
  VectorDouble center = getOrigin();
  VectorDouble dir = getUnitaryVector();
  for (const auto& h : vec_step)
    vec.push_back(evalIvarIpas(ivar, jvar, h, dir, center, mode));
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
                                          const CovCalcMode& mode) const
{
  int nvar = getNVariables();
  /// TODO : Not true whatever the space
  VectorDouble center = getOrigin();
  VectorDouble dir = getUnitaryVector();
  MatrixSquareGeneral mat(nvar);
  for (int ivar=0; ivar<nvar; ivar++)
    for (int jvar=0; jvar<nvar; jvar++)
      mat.setValue(ivar, jvar, evalIvarIpas(ivar, jvar, step, dir, center, mode));
  return mat;
}

/**
 * Calculate the (weighted) average Covariance between samples of two Dbs,
 * for a pair of variables
 * @param db1  Pointer to the first Db
 * @param db2  Pointer to the second Db
 * @param ivar Rank of the first variables
 * @param jvar Rank of the second variable
 * @param mode CovCalcMode structure
 * @return
 */
double ACov::evalAverageDbToDb(const Db* db1,
                               const Db* db2,
                               int ivar,
                               int jvar,
                               const CovCalcMode& mode)
{
  /* Loop on the first sample */

  double norme = 0.;
  double total = 0.;
  for (int iech1 = 0; iech1 < db1->getSampleNumber(); iech1++)
  {
    if (!db1->isActive(iech1)) continue;
    double w1 = db1->getWeight(iech1);
    if (w1 == 0.) continue;
    SpacePoint p1(db1->getSampleCoordinates(iech1),getSpace());

    /* Loop on the second sample */

    for (int iech2 = 0; iech2 < db2->getSampleNumber(); iech2++)
    {
      if (!db2->isActive(iech2)) continue;
      double w2 = db2->getWeight(iech2);
      if (w2 == 0.) continue;
      SpacePoint p2(db2->getSampleCoordinates(iech2),getSpace());

      /* Loop on the dimension of the space */

      total += w1 * w2 * eval(ivar, jvar, p1, p2, mode);
      norme += w1 * w2;
    }
  }

  // Scaling
  if (norme != 0.) total /= norme;

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
                                  const CovCalcMode& mode)
{
  /* Loop on the first sample */

  double norme = 0.;
  double total = 0.;

  /* Loop on the second sample */

  for (int iech2 = 0; iech2 < db2->getSampleNumber(); iech2++)
  {
    if (!db2->isActive(iech2)) continue;
    double w2 = db2->getWeight(iech2);
    if (w2 == 0.) continue;
    SpacePoint p2(db2->getSampleCoordinates(iech2),getSpace());

    /* Loop on the dimension of the space */

    total += w2 * eval(ivar, jvar, p1, p2, mode);
    norme += w2;
  }

  // Scaling
  if (norme != 0.) total /= norme;

  return total;
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
 * @param mode CovCalcMode structure
 * @return
 */
VectorDouble ACov::evalPointToDb(const SpacePoint& p1,
                                 const Db* db2,
                                 int ivar,
                                 int jvar,
                                 bool useSel,
                                 const CovCalcMode& mode)
{
  VectorDouble values;

  /* Loop on the second sample */

  for (int iech2 = 0; iech2 < db2->getSampleNumber(); iech2++)
  {
    if (db2->isActive(iech2))
    {
      SpacePoint p2(db2->getSampleCoordinates(iech2),getSpace());
      values.push_back(eval(ivar, jvar, p1, p2, mode));
    }
    else
    {
      if (! useSel) values.push_back(TEST);
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
                     const CovCalcMode& mode)
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

  double total = evalAverageDbToDb(dbgrid,  db, ivar, jvar, mode);
  delete dbgrid;
  return total;
}

/**
 * Average covariance over a block
 * @param p1     Point location
 * @param ext    Vector of Block extensions
 * @param ndisc  Vector of Block discretization
 * @param angles Vector of rotation angles
 * @param ivar   Rank of the first variable
 * @param jvar   Rank of the second variable
 * @param mode   CovCalcMode structure
 * @return
 */
double ACov::evalCxv(const SpacePoint& p1,
                     const VectorDouble& ext,
                     const VectorInt& ndisc,
                     const VectorDouble& angles,
                     int ivar,
                     int jvar,
                     const CovCalcMode& mode)
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
  DbGrid* dbgrid = _discretizeBlock(ext, ndisc, angles);
  if (dbgrid != nullptr)
    total = evalAveragePointToDb(p1, dbgrid, ivar, jvar, mode);
  delete dbgrid;

  return total;
}

/**
 * Creates the discretization grid
 * @param ext    Vecto of Block extensions
 * @param ndisc  Vector of Discretizations
 * @param angles Vector of rotation angles
 * @return
 */
DbGrid* ACov::_discretizeBlock(const VectorDouble& ext,
                               const VectorInt& ndisc,
                               const VectorDouble& angles)
{
  int ndim = getNDim();
  VectorDouble x0(ndim, 0.);
  VectorDouble dx(ndim, 0.);
  for (int idim = 0; idim < ndim; idim++)
    dx[idim] = ext[idim] / ndisc[idim];
  DbGrid* dbgrid = DbGrid::create(ndisc, dx, x0, angles);
  return dbgrid;
}

Db* ACov::_discretizeBlockRandom(const DbGrid* dbgrid)
{
  int ndim = getNDim();
  Db* db = Db::create();
  VectorString names = generateMultipleNames("x",ndim);

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
