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
#include "Covariances/ACor.hpp"
#include "Covariances/CovCalcMode.hpp"
#include "Matrix/MatrixSquareGeneral.hpp"
#include "Matrix/MatrixRectangular.hpp"
#include "Matrix/MatrixSparse.hpp"
#include "Matrix/MatrixSquareSymmetric.hpp"
#include "Matrix/NF_Triplet.hpp"
#include "Db/Db.hpp"
#include "Basic/AException.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/VectorNumT.hpp"
#include "Basic/VectorHelper.hpp"
#include "Space/ASpace.hpp"
#include "Space/SpacePoint.hpp"
#include "geoslib_define.h"

#include <vector>
#include <math.h>

ACor::ACor(const ASpace *space,int nvar)
    : ASpaceObject(space),
      _nvar(nvar),
      _optimEnabled(true),
      _isOptimPreProcessed(false),
      _p1As(),
      _p2A(space)
{
}

ACor::ACor(const ACor &r)
    : ASpaceObject(r),
      _nvar(r._nvar),
      _optimEnabled(r._optimEnabled),
      _isOptimPreProcessed(false),
      _p1As(),
      _p2A(r.getSpace())
{
}

ACor& ACor::operator=(const ACor &r)
{
  if (this != &r)
  {
    ASpaceObject::operator=(r);
    _nvar = r._nvar;
    _optimEnabled = r._optimEnabled;
    _isOptimPreProcessed = r._isOptimPreProcessed;
    _p2A = SpacePoint(r._space);
  }
  return *this;
}

ACor::~ACor()
{
}

void ACor::optimizationPostProcess() const
{
  _optimizationPostProcess();
  _isOptimPreProcessed = false;
}

void ACor::optimizationPreProcess(const Db* db) const
{
  if (_isOptimPreProcessed) return;
  db->getSamplesAsSP(_p1As,_space);
  _optimizationPreProcess(_p1As);
  
}

void ACor::optimizationSetTarget(const SpacePoint &pt) const
{
  _optimizationSetTarget(pt);
}

void ACor::_optimizationSetTarget(const SpacePoint &pt) const
{
  _p2A = pt;
}


void ACor::optimizationPreProcess(const std::vector<SpacePoint>& p) const
{
  _optimizationPreProcess(p);
  _isOptimPreProcessed = isOptimEnabled();
}

//Preprocess by default (only copy the points)
void ACor::_optimizationPreProcess(const std::vector<SpacePoint>& p) const
{
  if (!_p1As.empty() && _isOptimPreProcessed) return;
  _p1As.clear();
  for (const auto &e : p)
  {
   _p1As.push_back(e);
  }
}

void ACor::_optimizationPostProcess() const
{
    _p1As.clear();
}

MatrixSquareGeneral ACor::eval0Mat(const CovCalcMode* mode) const
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
void ACor::eval0CovMatBiPointInPlace(MatrixSquareGeneral& mat, const CovCalcMode* mode) const
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
void ACor::addEval0CovMatBiPointInPlace(MatrixSquareGeneral& mat, const CovCalcMode* mode) const
{
  int nvar = getNVariables();

  for (int ivar=0; ivar<nvar; ivar++)
    for (int jvar=0; jvar<nvar; jvar++)
        mat.addValue(ivar, jvar, eval0(ivar, jvar, mode)); // pure virtual method
}

VectorDouble ACor::eval(const std::vector<SpacePoint>& vec_p1,
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

double ACor::eval0(int ivar,
                   int jvar,
                   const CovCalcMode* mode) const
{
  DECLARE_UNUSED(ivar,jvar,mode)
  return 1.;
}
MatrixSquareGeneral ACor::evalMat(const SpacePoint& p1,
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
void ACor::evalCovMatBiPointInPlace(MatrixSquareGeneral &mat,
                                    const SpacePoint &p1,
                                    const SpacePoint &p2,
                                    const CovCalcMode* mode) const
{
  mat.fill(0.);
  addEvalCovMatBiPointInPlace(mat,p1,p2,mode);
}

VectorInt ACor::_getActiveVariables(int ivar0) const
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
MatrixRectangular ACor::evalCovMatrix(const Db* db1,
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

void ACor::_updateCovMatrixSymmetricVerr(const Db *db1,
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

void ACor::addEvalCovMatBiPointInPlace(MatrixSquareGeneral &mat,
                        const SpacePoint& pwork1, 
                        const SpacePoint& pwork2, 
                        const CovCalcMode *mode) const
{
  _addEvalCovMatBiPointInPlace(mat, pwork1, pwork2,mode);
}

void ACor::_addEvalCovMatBiPointInPlace(MatrixSquareGeneral &mat,
                        const SpacePoint& pwork1, 
                        const SpacePoint& pwork2, 
                        const CovCalcMode *mode) const
{
  for (int ivar = 0, nvar = getNVariables(); ivar < nvar; ivar++)
    for (int jvar = 0; jvar < nvar; jvar++)
      mat.addValue(ivar, jvar, eval(pwork1,pwork2,ivar,jvar,mode)); 
}

void ACor::evalCovKriging(MatrixSquareGeneral &mat,
                          SpacePoint &pwork1,
                          SpacePoint& pout, 
                          const CovCalcMode *mode) const
{
  mat.fill(0.);
  loadAndAddEvalCovMatBiPointInPlace(mat,pwork1,pout,mode);
}

void ACor::loadAndAddEvalCovMatBiPointInPlace(MatrixSquareGeneral &mat,const SpacePoint& p1,const SpacePoint&p2,
                                              const CovCalcMode *mode) const
{
  _loadAndAddEvalCovMatBiPointInPlace(mat,p1,p2,mode);
}

double ACor::loadAndEval(const SpacePoint& p1,
                          const SpacePoint&p2,
                          int ivar,
                          int jvar,
                          const CovCalcMode *mode) const
{
  return _loadAndEval(p1,p2,ivar,jvar,mode);
}

double ACor::_loadAndEval(const SpacePoint& p1,
                          const SpacePoint&p2,
                          int ivar,
                          int jvar,
                          const CovCalcMode *mode) const
{
  load(p1,true);
  load(p2,false);
  return eval(*_pw1,*_pw2,ivar,jvar,mode);

}
void ACor::_loadAndAddEvalCovMatBiPointInPlace(MatrixSquareGeneral &mat,
                                              const SpacePoint& p1,
                                              const SpacePoint&p2,
                                              const CovCalcMode *mode) const
{
  load(p1,true);
  load(p2,false);
  _addEvalCovMatBiPointInPlace(mat,*_pw1,*_pw2,mode);
}

void ACor::load(const SpacePoint& p,bool case1) const
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
MatrixSquareSymmetric ACor::evalCovMatrixSymmetric(const Db *db1,
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
MatrixSparse* ACor::evalCovMatrixSparse(const Db *db1,
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

