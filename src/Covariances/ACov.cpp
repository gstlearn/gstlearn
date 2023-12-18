/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://github.com/gstlearn                                       */
/* License: BSD 3 clauses                                                     */
/*                                                                            */
/******************************************************************************/
#include "Covariances/ACov.hpp"
#include "Matrix/MatrixSquareGeneral.hpp"
#include "Matrix/MatrixRectangular.hpp"
#include "Db/Db.hpp"
#include "Db/DbGrid.hpp"
#include "Basic/AException.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/VectorNumT.hpp"
#include "Basic/VectorHelper.hpp"
#include "Basic/Law.hpp"

#include <vector>
#include <math.h>

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

VectorDouble ACov::eval(const std::vector<SpacePoint>& vec_p1,
                        const std::vector<SpacePoint>& vec_p2,
                        int ivar,
                        int jvar,
                        const CovCalcMode& mode) const
{
  VectorDouble vec;
  if (vec_p1.size() != vec_p2.size())
    my_throw ("Error: 'p1' and 'p2' should have same dimension");
  for (int i=0, n=static_cast<int> (vec_p1.size()); i < n; i++)
    vec.push_back(eval(vec_p1[i], vec_p2[i], ivar, jvar, mode)); // pure virtual method
  return vec;
}

MatrixSquareGeneral ACov::evalMat(const SpacePoint& p1,
                                  const SpacePoint& p2,
                                  const CovCalcMode& mode) const
{
  int nvar = getNVariables();
  MatrixSquareGeneral mat(nvar);
  for (int ivar=0; ivar<nvar; ivar++)
    for (int jvar=0; jvar<nvar; jvar++)
      mat.setValue(ivar, jvar, eval(p1, p2, ivar, jvar, mode)); // pure virtual method
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
double ACov::evalIvarIpas(double step,
                          const VectorDouble& dir,
                          int ivar,
                          int jvar,
                          const VectorDouble& center,
                          const CovCalcMode& mode) const
{
  // Define the point in the ACov space (center will be checked)
  const ASpace* space = getSpace();
  SpacePoint p1(center,space);
  SpacePoint p2(center,space);
  VectorDouble dirloc(dir);
  if (dirloc.empty())
  {
    dirloc.resize(getNDim(),0.);
    dirloc[0] = 1.;
  }
  VectorDouble vec(dirloc);
  VH::multiplyConstant(vec, step);
  p2.move(vec);
  return eval(p1, p2, ivar, jvar, mode); // pure virtual method
}

double ACov::evalIvarIpasIncr(const VectorDouble& dincr,
                              int ivar,
                              int jvar,
                              const CovCalcMode& mode) const
{
  // Define the point in the ACov space (center will be checked)
  SpacePoint p1(VectorDouble(),getSpace());
  SpacePoint p2(VectorDouble(),getSpace());
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
 * @param center    Coordinates of the Center
 * @param mode      CovCalcMode structure
 * @return
 */
VectorDouble ACov::evalIvarNpas(const VectorDouble& vec_step,
                                const VectorDouble& dir,
                                int ivar,
                                int jvar,
                                const VectorDouble& center,
                                const CovCalcMode& mode) const
{
  VectorDouble vec;
  for (int i=0, n=static_cast<int> (vec_step.size()); i < n; i++)
    vec.push_back(evalIvarIpas(vec_step[i], dir, ivar, jvar, center, mode));
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
      mat.setValue(ivar, jvar, evalIvarIpas(step, dir, ivar, jvar, center, mode));
  return mat;
}

MatrixSquareGeneral ACov::evalNvarIpasIncr(const VectorDouble& dincr,
                                           const CovCalcMode& mode) const
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
                             const CovCalcMode& mode) const
{
  /// TODO : Not true whatever the space
  VectorDouble center = getOrigin();
  VectorDouble dir = getUnitaryVector();
  return evalIvarIpas(step, dir, ivar, jvar, center, mode);
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
                                   const CovCalcMode& mode) const
{
  VectorDouble vec;
  /// TODO : Not true whatever the space
  VectorDouble center = getOrigin();
  VectorDouble dir = getUnitaryVector();
  for (const auto& h : vec_step)
    vec.push_back(evalIvarIpas(h, dir, ivar, jvar, center, mode));
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
      mat.setValue(ivar, jvar, evalIvarIpas(step, dir, ivar, jvar, center, mode));
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
                               const CovCalcMode& mode) const
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

      total += w1 * w2 * eval(p1, p2, ivar, jvar, mode);
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
                                  const CovCalcMode& mode) const
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

    total += w2 * eval(p1, p2, ivar, jvar, mode);
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
                                 const CovCalcMode& mode) const
{
  VectorDouble values;

  /* Loop on the second sample */

  for (int iech2 = 0; iech2 < db2->getSampleNumber(); iech2++)
  {
    if (db2->isActive(iech2))
    {
      SpacePoint p2(db2->getSampleCoordinates(iech2),getSpace());
      values.push_back(eval(p1, p2, ivar, jvar, mode));
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
                     const CovCalcMode& mode) const
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
                          const CovCalcMode& mode) const
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

  double total = evalAverageDbToDb(dbgrid1,  dbgrid2, ivar, jvar, mode);
  delete dbgrid1;
  delete dbgrid2;
  return total;
}

MatrixSquareGeneral ACov::evalCvvM(const VectorDouble& ext,
                                   const VectorInt& ndisc,
                                   const VectorDouble& angles,
                                   const CovCalcMode& mode) const
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
                     const CovCalcMode& mode) const
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
                     const CovCalcMode& mode) const
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
    total = evalAverageDbToDb(db, dbgrid, ivar, jvar, mode);
  delete dbgrid;

  return total;
}

MatrixSquareGeneral ACov::evalCxvM(const SpacePoint& p1,
                                   const VectorDouble& ext,
                                   const VectorInt& ndisc,
                                   const VectorDouble& angles,
                                   const VectorDouble& x0,
                                   const CovCalcMode& mode) const
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

/****************************************************************************/
/*!
 **  Establish the covariance matrix between two Dbs
 **
 ** \param[in]  db1   First Db
 ** \param[in]  db2   Second Db (=db1 if absent)
 ** \param[in]  ivar  Rank of the first variable (-1: all variables)
 ** \param[in]  jvar  Rank of the second variable (-1: all variables)
 ** \param[in]  mode  CovCalcMode structure
 **
 ** \remarks The returned matrix if dimension to nrows * ncols where
 ** \remarks nrows is the number of active samples in db1
 ** \remarks ncols is the number of active samples in db2
 **
 *****************************************************************************/
MatrixRectangular ACov::evalCovMatrix(const Db* db1,
                                      const Db* db2,
                                      int ivar,
                                      int jvar,
                                      const CovCalcMode& mode) const
{
  if (db2 == nullptr) db2 = db1;
  int nechtot1 = db1->getSampleNumber(false);
  int nechtot2 = db2->getSampleNumber(false);
  int nech1 = db1->getSampleNumber(true);
  int nech2 = db2->getSampleNumber(true);
  MatrixRectangular mat(nech1, nech2);

  /* Loop on the first sample */

  int jech1 = 0;
  for (int iech1 = 0; iech1 < nechtot1; iech1++)
  {
    if (!db1->isActive(iech1)) continue;
    SpacePoint p1(db1->getSampleCoordinates(iech1),getSpace());

    /* Loop on the second sample */

    int jech2 = 0;
    for (int iech2 = 0; iech2 < nechtot2; iech2++)
    {
      if (!db2->isActive(iech2)) continue;
      SpacePoint p2(db2->getSampleCoordinates(iech2),getSpace());

      /* Loop on the dimension of the space */

      double value = eval(p1, p2, ivar, jvar, mode);
      mat.setValue(jech1, jech2, value);
      jech2++;
    }
    jech1++;
  }
  return mat;
}

void ACov::preProcessFromDb(const Db* db) const
{
	int nechtot1 = db->getSampleNumber(false);
	int nech1 = db->getSampleNumber(true);
	std::vector<SpacePoint> pvect(nech1);
	VectorDouble temp(nech1);

	int jech1 = 0;

	for (int iech1 = 0; iech1 <nechtot1;iech1++)
	{
	  if (!db->isActive(iech1)) continue;
	  pvect[jech1] = db->getSampleCoordinates(iech1);
	  jech1++;
	}

	preProcess(pvect);

}

MatrixEigen ACov::evalCovMatrixEigen(const Db* db1,
                                     const Db* db2,
									 bool preprocess,
									 int ivar,
									 int jvar,
									 const CovCalcMode& mode) const
{

	int nech1 = db1->getSampleNumber(true);
	VectorDouble temp(nech1);

	MatrixEigen mat;
	if (preprocess)
		preProcessFromDb(db1);

	SpacePoint p1(db1->getSampleCoordinates(0),getSpace());
	SpacePoint ptemp(db1->getSampleCoordinates(0),getSpace());




	if (db2 == nullptr)
	{
		mat = MatrixEigen(nech1,nech1);
		evalOptimEigen(mat,temp,ivar, jvar, mode);
	}
	else
	{
		int nechtot2 = db2->getSampleNumber(false);
		int nech2 = db2->getSampleNumber(true);
		mat = MatrixEigen(nech1,nech2);
		int jech2 = 0;
		for (int iech2 =  0; iech2 < nechtot2; iech2++)
		{
			if (!db2->isActive(iech2)) continue;
			db2->getSampleCoordinates(iech2, p1.getCoordM());
			evalOptimEigen(p1,ptemp,mat,jech2,jech2,temp,ivar, jvar, mode);
			jech2++;
		}
	}
	 if (preprocess)
		 cleanPreProcessInfo();
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
                               const VectorDouble &ext,
                               const VectorInt &ndisc,
                               const VectorDouble &angles,
                               const VectorDouble& x0,
                               int ivar,
                               int jvar) const
{
  CovCalcMode mode = CovCalcMode();
  mode.setAsVario(true);

  double sigmaE = TEST;
  DbGrid* dbgrid = _discretizeBlock(ext, ndisc, angles, x0);
  if (dbgrid != nullptr)
  {
    double GxV = evalAverageDbToDb(db, dbgrid, ivar, jvar, mode);
    double Gxx = evalAverageDbToDb(db, db, ivar, jvar, mode);
    double GVV = evalAverageDbToDb(dbgrid, dbgrid, ivar, jvar, mode);
    sigmaE = 2. * GxV - Gxx - GVV;
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
