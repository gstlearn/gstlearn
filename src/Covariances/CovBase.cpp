/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c3) MINES Paris / ARMINES                                       */
/* Authors: gstlearn Team                                                     */
/* Website: https://gstlearn.org                                              */
/* License: BSD 3-clause                                                      */
/*                                                                            */
/******************************************************************************/

#include "Covariances/CovBase.hpp"
#include "Basic/ParamInfo.hpp"
#include "Basic/VectorNumT.hpp"
#include "Covariances/ACov.hpp"
#include "Covariances/CovContext.hpp"
#include "Matrix/MatrixRectangular.hpp"
#include "Matrix/MatrixSquareSymmetric.hpp"
#include "Db/Db.hpp"
#include "Covariances/NoStatArray.hpp"
#include "Covariances/NoStatFunctional.hpp"
#include "geoslib_define.h"
#include <cstddef>
#include <functional>

ParamInfo CovBase::createParamInfoForCholSill(int ivar, int jvar)
{
  std::function<void(double)> setCholSill = [this, ivar, jvar](double value)
  {
    this->setCholSill(ivar, jvar, value);
  };
  ParamInfo pinf(String("Cholesky sill"),
                 TEST,
                 {-INF, INF},
                 String("Term of the Cholesky decomposition of the sill matrix"));
  return pinf;
}
CovBase::CovBase(ACov* cor,
                 const MatrixSquareSymmetric& sill)
  : ACov(cor == nullptr ? CovContext() : cor->getContext())
  , _cholSillsInfo(MatrixT<ParamInfo>(sill.getNRows(), sill.getNCols(), createParamInfoForCholSill()))
  , _cholSills(MatrixRectangular(sill.getNRows(), sill.getNCols()))
  , _sillCur(sill)
  , _cor(cor)
{
  _ctxt.setNVar(sill.getNSize());
  for (size_t i = 0, n = getNVar(); i < n; i++)
  {

    for (size_t j = 0; j <= n; j++)
    {
    }
    for (size_t j = i + 1; j < n; j++)
    {
      _cholSillsInfo(i, j).setFixed(true);
    }
  }
  if (cor != nullptr)
  {
    _ctxt = cor->getContextCopy();
  }
  _ctxt.setNVar(sill.getNSize());
  _workMat.resize(_ctxt.getNVar(), _ctxt.getNVar());
  _workMat.setIdentity();
}

CovBase::CovBase(const CovBase& r)
  : ACov(r)
{
  _cholSillsInfo = r._cholSillsInfo;
  _cholSills     = r._cholSills;
  _tabNoStat     = r._tabNoStat;
  _sillCur       = r._sillCur;
  _workMat       = r._workMat;
  _cor           = (ACov*)r._cor->clone();
}

CovBase& CovBase::operator=(const CovBase& r)
{
  if (this != &r)
  {
    _cholSillsInfo = r._cholSillsInfo;
    _cholSills     = r._cholSills;
    _tabNoStat     = r._tabNoStat;
    _sillCur       = r._sillCur;
    _workMat       = r._workMat;
    _cor           = (ACov*)r._cor->clone();
  }
  return *this;
}

CovBase::~CovBase()
{
}

void CovBase::loadInfoValues()
{
  for (size_t ivar = 0, n = getNVar(); ivar < n; ivar++)
  {
    for (size_t jvar = 0; jvar < n; jvar++)
    {
      // _cholSills.setValue(ivar, jvar, _cholSillsInfo(ivar, jvar).getValue());
    }
  }
  _sillCur.prodMatMatInPlace(&_cholSills, &_cholSills, false, true);
  _cor->loadInfoValues();
}
void CovBase::setCor(ACov* cor)
{
  _cor     = cor;
  int nvar = getNVar();
  if (cor != nullptr)
  {
    _ctxt = cor->getContextCopy();
    _ctxt.setNVar(nvar);
  }
}
void CovBase::_setContext(const CovContext& ctxt)
{
  _cor->setContext(ctxt);
  _updateFromContext();
}

void CovBase::setSill(double sill) const
{
  int nvar = getNVar();
  if (nvar > 0 && nvar != 1)
  {
    messerr("Number of provided sill doesn't match number of variables");
    return;
  }
  _sillCur.resetFromValue(1, 1, sill);
}

void CovBase::setSill(const MatrixSquareSymmetric& sill) const
{
  int nvar = getNVar();
  if (nvar > 0 && nvar != sill.getNCols())
  {
    messerr("Number of provided sills doesn't match number of variables");
    return;
  }
  _sillCur = sill;
}

void CovBase::setSill(const VectorDouble& sill) const
{
  int size = static_cast<int>(sill.size());
  int nvar = getNVar();
  if (size != nvar * nvar)
  {
    messerr("Number of provided sills doesn't match number of variables");
    return;
  }
  _sillCur.setValues(sill);
}

void CovBase::setSill(int ivar, int jvar, double sill) const
{
  if (!_isVariableValid(ivar)) return;
  if (!_isVariableValid(jvar)) return;
  /// TODO : Test if sill matrix is positive definite (if not, generate a warning)
  if (!_sillCur.isValid(ivar, jvar)) return;
  _sillCur.setValue(ivar, jvar, sill);
}

void CovBase::setCholSill(int ivar, int jvar, double val) const
{
  if (!_isVariableValid(ivar)) return;
  if (!_isVariableValid(jvar)) return;
  if (ivar < jvar)
  {
    messerr("The Cholesky decomposition of the sill matrix is lower triangular");
    return;
  }
  _cholSills.setValue(ivar, jvar, val);
}

bool CovBase::_isVariableValid(int ivar) const
{
  return checkArg("Rank of the Variable", ivar, getNVar());
}

void CovBase::_initFromContext()
{
  _cor->initFromContext();
  _sillCur.reset(_ctxt.getNVar(), _ctxt.getNVar());
  setOptimEnabled(true);

}
void CovBase::initSill(double value)
{
  _sillCur.fill(value);
}

bool CovBase::isConsistent(const ASpace* space) const
{
  return _cor->isConsistent(space);
}

int CovBase::addEvalCovVecRHSInPlace(vect vect,
                                     const VectorInt& index1,
                                     SpacePoint& pin,
                                     SpacePoint& pout,
                                     const int iech2) const
{
  return getSill(0, 0) * _cor->addEvalCovVecRHSInPlace(vect, index1, pin, pout, iech2);
}

double CovBase::_eval(const SpacePoint& p1,
                      const SpacePoint& p2,
                      int ivar,
                      int jvar,
                      const CovCalcMode* mode) const
{
  return getSill(ivar, jvar) * _cor->evalCov(p1, p2, ivar, jvar, mode);
}

double CovBase::getSill(int ivar, int jvar) const
{
  return _sillCur.getValue(ivar, jvar);
}

/*****************************************************************************/
/*!
 **  Update the Model in the case of Non-stationary parameters
 **  This requires the knowledge of the two end-points
 **
 ** \param[in]  covint       Internal structure for non-stationarity
 **                          or NULL (for stationary case)
 **
 *****************************************************************************/
void CovBase::nostatUpdate(CovInternal* covint) const
{
  if (covint == NULL) return;
  updateCovByPoints(covint->getIcas1(), covint->getIech1(),
                    covint->getIcas2(), covint->getIech2());
}

void CovBase::_copyCovContext(const CovContext& ctxt)
{
  _ctxt.copyCovContext(ctxt);
  _cor->copyCovContext(ctxt);
}

/**
 * Transform a set of Space Points using the anisotropy tensor
 * The set of resulting Space Points are stored as private member of this.
 * Note that ALL samples are processed, independently from the presence of a selection
 * or checking for heterotopy.
 * @param mode 1 for p1As; 2for p2As
 * @param ps vector of SpacePoints
 */
void CovBase::_optimizationPreProcess(int mode, const std::vector<SpacePoint>& ps) const
{
  _cor->optimizationPreProcess(mode, ps);
}

SpacePoint& CovBase::_optimizationLoadInPlace(int iech, int mode, int rank) const
{
  return _cor->optimizationLoadInPlace(iech, mode, rank);
}

/**
 * Checks that the Optimization has already been initiated, by:
 * - checking that the storage (for Sample Points projected in the Covariance
 * rotation system) is already allocated
 * - checking that the dimension of this storage is correct (only if 'db' is provided):
 * in particular, this check is not necessary when freeing this storage.
 */
bool CovBase::isOptimizationInitialized(const Db* db) const
{
  if (_p1As.empty()) return false;
  if (db == nullptr) return true;
  int n = (int)_p1As.size();
  return n == db->getNSample();
}

// Set of functions to make parameters no stationary (or to make them back stationary).
// There is to types of non stationarities : NoStatDb in which the parameters are read in a
// DbGrid or NoStatFunctional for which you have to provide a function of the coordinates.
// Each parameter can have its own type of No stationarity and its own DbGrid in case
// of NoStatDb.
// For specifying the NoStat DbGrid, you can first attach it by using attachNoStatDb.
// If not, you have to specify the DbGrid when you make the first parameter non stationary.

void CovBase::attachNoStatDb(const Db* db)
{
  _tabNoStat.setDbNoStatRef(db);
  _cor->attachNoStatDb(db);
}

bool CovBase::_checkAndManageNoStatDb(const Db*& db, const String& namecol)
{
  if (_tabNoStat.getDbNoStatRef() == nullptr && db == nullptr)
  {
    messerr("You have to define a Db (with attachNoStatDb or by specifying a Db here)");
    return false;
  }
  _setNoStatDbIfNecessary(db);

  if (db->getUID(namecol) < 0)
  {
    messerr("You have to specified a name of a column of the reference Db");
    return false;
  }
  return true;
}

void CovBase::_setNoStatDbIfNecessary(const Db*& db)
{
  if (_tabNoStat.getDbNoStatRef() == nullptr)
    attachNoStatDb(db);
  if (db == nullptr)
    db = _tabNoStat.getDbNoStatRef();
}

void CovBase::_makeElemNoStat(const EConsElem& econs, int iv1, int iv2, const AFunctional* func, const Db* db, const String& namecol)
{
  if (func == nullptr)
  {
    if (!_checkAndManageNoStatDb(db, namecol)) return;
  }

  if (econs != EConsElem::SILL)
  {
    _cor->makeElemNoStat(econs, iv1, iv2, func, db, namecol);
    return;
  }

  std::shared_ptr<ANoStat> ns;
  if (func == nullptr)
  {
    ns = std::shared_ptr<ANoStat>(new NoStatArray(db, namecol));
  }
  else
  {
    ns = std::unique_ptr<ANoStat>(new NoStatFunctional(func));
  }

  _tabNoStat.addElem(ns, econs, iv1, iv2);
}

///////////////////// Sill ////////////////////////

void CovBase::makeSillNoStatDb(const String& namecol, int ivar, int jvar, const Db* db)
{
  if (!_checkSill(ivar, jvar)) return;
  _makeElemNoStat(EConsElem::SILL, ivar, jvar, nullptr, db, namecol);
  _cor->checkAndManageNoStatDb(db, namecol);
}

void CovBase::makeSillNoStatFunctional(const AFunctional* func, int ivar, int jvar)
{
  if (!_checkSill(ivar, jvar)) return;
  _makeElemNoStat(EConsElem::SILL, ivar, jvar, func);
}

void CovBase::makeSillsStationary(bool silent)
{
  if (_tabNoStat.getNSills() == 0 && !silent)
  {
    messerr("All the sills are already stationary!");
    return;
  }
  _tabNoStat.clear();
}
void CovBase::makeSillStationary(int ivar, int jvar)
{
  if (!_checkSill(ivar, jvar)) return;
  if (_tabNoStat.removeElem(EConsElem::SILL, ivar, jvar) == 0)
  {
    messerr("This parameter was already stationary!");
  }
}

/////////////////////////// Check functions ////////////////////:

bool CovBase::_checkSill(int ivar, int jvar) const
{
  int nvar = getNVar();
  if ((ivar > nvar) || (jvar > nvar))
  {
    messerr("Your model has only %d variables.", nvar);
    return false;
  }
  return true;
}

bool CovBase::_checkDims(int idim, int jdim) const
{
  int ndim = getNDim();
  if ((idim > ndim) || (jdim > ndim))
  {
    messerr("Your model is only in dimension %d.", ndim);
    return false;
  }
  return true;
}

/////////////  Functions to attach no stat information on various supports ////////
void CovBase::informMeshByMesh(const AMesh* amesh) const
{
  _tabNoStat.informMeshByMesh(amesh);
  _cor->informMeshByMesh(amesh);
}
void CovBase::informMeshByApex(const AMesh* amesh) const
{
  _tabNoStat.informMeshByApex(amesh);
  _cor->informMeshByApex(amesh);
}
void CovBase::informDbIn(const Db* dbin) const
{
  _tabNoStat.informDbIn(dbin);
  _cor->informDbIn(dbin);
}
void CovBase::informDbOut(const Db* dbout) const
{
  _tabNoStat.informDbOut(dbout);
  _cor->informDbOut(dbout);
}

double CovBase::getValue(const EConsElem& econs, int iv1, int iv2) const
{
  double val = _cor->getValue(econs, iv1, iv2);
  if (val == TEST)
  {
    if (econs == EConsElem::SILL)
      return getSill(iv1, iv2);
  }
  return val;
}

VectorDouble CovBase::informCoords(const VectorVectorDouble& coords,
                                   const EConsElem& econs,
                                   int iv1,
                                   int iv2) const
{
  if (econs == EConsElem::SILL)
  {
    VectorDouble result(coords[0].size(), getValue(econs, iv1, iv2));
    _tabNoStat.informCoords(coords, econs, iv1, iv2, result);
    return result;
  }

  return _cor->informCoords(coords, econs, iv1, iv2);
}

void CovBase::informMeshByMeshForSills(const AMesh* amesh) const
{
  _tabNoStat.informMeshByMesh(amesh, EConsElem::SILL);
}

void CovBase::informMeshByApexForSills(const AMesh* amesh) const
{
  _tabNoStat.informMeshByApex(amesh, EConsElem::SILL);
}

void CovBase::informDbInForSills(const Db* dbin) const
{
  _tabNoStat.informDbIn(dbin, EConsElem::SILL);
}

void CovBase::informDbOutForSills(const Db* dbout) const
{
  _tabNoStat.informDbOut(dbout, EConsElem::SILL);
}

/**
 * Update the Model according to the Non-stationary parameters
 * @param icas1 Type of first Db: 1 for Input; 2 for Output
 * @param iech1 Rank of the target within Db1 (or -1)
 * @param icas2 Type of first Db: 1 for Input; 2 for Output
 * @param iech2 Rank of the target within Dbout (or -2)
 */
void CovBase::updateCovByPoints(int icas1, int iech1, int icas2, int iech2) const
{
  // If no non-stationary parameter is defined, simply skip
  if (!isNoStat()) return;
  double val1, val2;

  const auto paramsnostat = _tabNoStat.getTable();
  // Loop on the elements that can be updated one-by-one

  for (const auto& e: paramsnostat)
  {
    EConsElem type = e.first.getType();
    e.second->getValuesOnDb(icas1, iech1, &val1, icas2, iech2, &val2);

    if (type == EConsElem::SILL)
    {
      int iv1 = e.first.getIV1();
      int iv2 = e.first.getIV2();
      setSill(iv1, iv2, sqrt(val1 * val2));
    }
  }
  _cor->updateCovByPoints(icas1, iech1, icas2, iech2);
}

void CovBase::updateCovByMesh(int imesh, bool aniso) const
{
  // If no non-stationary parameter is defined, simply skip
  if (!isNoStat()) return;

  // Loop on the elements that can be updated one-by-one
  if (!aniso)
  {
    const auto paramsnostat = _tabNoStat.getTable();
    for (const auto& e: paramsnostat)
    {
      EConsElem type = e.first.getType();
      if (type == EConsElem::SILL)
      {
        double sill = e.second->getValueOnMeshByApex(imesh);
        int iv1     = e.first.getIV1();
        int iv2     = e.first.getIV2();
        setSill(iv1, iv2, sill);
      }
    }
  }
  _cor->updateCovByMesh(imesh, aniso);
}

void CovBase::makeStationary()
{
  _cor->makeStationary();
  makeSillsStationary(true);
}

void CovBase::_manage(const Db* db1, const Db* db2) const
{
  if (db1 != nullptr)
    informDbIn(db1);
  if (db2 != nullptr)
    informDbOut(db2);
  _cor->manage(db1, db2);
}

void CovBase::_optimizationPostProcess() const
{
  _cor->optimizationPostProcess();
}

void CovBase::_updateFromContext()
{
  _cor->updateFromContext();
}

void CovBase::_load(const SpacePoint& p, bool case1) const
{
  _cor->load(p, case1);
}

void CovBase::_optimizationSetTarget(SpacePoint& pt) const
{
  _cor->optimizationSetTarget(pt);
}
