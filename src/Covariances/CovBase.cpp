/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c3) MINES Paris / ARMINES                                        */
/* Authors: gstlearn Team                                                     */
/* Website: https://gstlearn.org                                              */
/* License: BSD 3-clause                                                      */
/*                                                                            */
/******************************************************************************/

#include "Covariances/CovBase.hpp"
#include "Covariances/ACov.hpp"
#include "Covariances/CovContext.hpp"
#include "Matrix/MatrixSquareSymmetric.hpp"
#include "Db/Db.hpp"
#include "Covariances/NoStatArray.hpp"
#include "Covariances/NoStatFunctional.hpp"

CovBase::CovBase(ACov* cor,
                const MatrixSquareSymmetric &sill)
: ACov(cor == nullptr? CovContext() : cor->getContext())
, _sill(sill)
, _cor(cor)
{
  if (cor != nullptr)
  {
    _ctxt = cor->getContextCopy();
  }
  _ctxt.setNVar(sill.getNSize());
  _workMat.resize(_ctxt.getNVar(), _ctxt.getNVar());
  _workMat.setIdentity();
}

CovBase::~CovBase()
{

}

void CovBase::setCor(ACov* cor)
{
  _cor = cor;
  int nvar = getNVariables();
  if (cor != nullptr)
  {
    _ctxt = cor->getContextCopy();
    _ctxt.setNVar(nvar);
  }
}
void CovBase::setContext(const CovContext &ctxt)
{
  _ctxt = ctxt;
  _updateFromContext();
}
void CovBase::setSill(double sill) const
{
  int nvar = getNVariables();
  if (nvar > 0 && nvar!= 1)
  {
    messerr("Number of provided sill doesn't match number of variables");
    return;
  }
  _sill.resetFromValue(1, 1, sill);
}

void CovBase::setSill(const MatrixSquareSymmetric &sill) const
{
  int nvar = getNVariables();
  if (nvar > 0 && nvar != sill.getNCols())
  {
    messerr("Number of provided sills doesn't match number of variables");
    return;
  }
  _sill = sill;
}

void CovBase::setSill(const VectorDouble &sill) const 
{
  int size = static_cast<int>(sill.size());
  int nvar = getNVariables();
  if (size != nvar * nvar)
  {
    messerr("Number of provided sills doesn't match number of variables");
    return;
  }
  _sill.setValues(sill);
}

void CovBase::setSill(int ivar, int jvar, double sill) const
{
  if (!_isVariableValid(ivar)) return;
  if (!_isVariableValid(jvar)) return;
  /// TODO : Test if sill matrix is positive definite (if not, generate a warning)
  if (!_sill.isValid(ivar, jvar)) return;
  _sill.setValue(ivar, jvar, sill);
}

bool CovBase::_isVariableValid(int ivar) const
{
  return checkArg("Rank of the Variable", ivar, getNVariables());
}

void CovBase::_initFromContext()
{
  _cor->initFromContext();
  _sill.reset(_ctxt.getNVar(), _ctxt.getNVar());
}
void CovBase::initSill(double value)
{
  _sill.fill(value);
}

bool CovBase::isConsistent(const ASpace* space) const
{
  return _cor->isConsistent(space);
}



double CovBase::getSill(int ivar, int jvar) const
{
  return _sill.getValue(ivar, jvar);
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
void CovBase::nostatUpdate(CovInternal *covint) const
{
  if (covint == NULL) return;
  updateCovByPoints(covint->getIcas1(), covint->getIech1(),
                    covint->getIcas2(), covint->getIech2());
}



void CovBase::_copyCovContext(const CovContext &ctxt)
{
  _ctxt.copyCovContext(ctxt);
  _cor->copyCovContext(ctxt);
}

/**
 * Define the Second Space Point as coinciding with the Input Space Point 'iech'.
 * Note that, as the Input Space Points are already transformed in the basis
 * of the current structure, it is just an assignment.
 *
 * @param iech Rank of the sample among the recorded Space Points
 */
void CovBase::optimizationSetTargetByIndex(int iech) const
{
  if (_isOptimPreProcessed)
  {
    _p2A = _p1As[iech];
    _p2A.setTarget(true);
  }
}

/**
 * Transform a set of Space Points using the anisotropy tensor
 * The set of resulting Space Points are stored as private member of this.
 * Note that ALL samples are processed, independently from the presence of a selection
 * or checking for heterotopy.
 * @param p vector of SpacePoints
 */
void CovBase::_optimizationPreProcess(const std::vector<SpacePoint>& p) const
{
  if (!isOptimEnabled())
  {
     ACov::_optimizationPreProcess(p);
     return;
  }
  _cor->optimizationPreProcess(p,_p1As);
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
  int n = (int) _p1As.size();
  return n == db->getSampleNumber();
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

bool CovBase::_checkAndManageNoStatDb(const Db*&  db, const String& namecol)
{
 if (_tabNoStat.getDbNoStatRef() == nullptr && db == nullptr)
 {
  messerr("You have to define a Db (with attachNoStatDb or by specifying a Db here)");  
  return false;
 }
  _setNoStatDbIfNecessary(db);

 if (db->getUID(namecol)< 0)
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

void CovBase::_makeElemNoStat(const EConsElem &econs, int iv1, int iv2,const AFunctional* func, const Db* db, const String& namecol)
{
  if (func == nullptr)
  {
    if(!_checkAndManageNoStatDb(db,namecol)) return;
  }

  if (econs != EConsElem::SILL)
  {
    _cor->makeElemNoStat(econs,iv1,iv2,func,db,namecol);
    return;
  }
  

  std::shared_ptr<ANoStat> ns;
  if (func == nullptr)
  {
    ns = std::shared_ptr<ANoStat>(new NoStatArray(db,namecol));
  }
  else 
  {
    ns = std::unique_ptr<ANoStat>(new NoStatFunctional(func));
  }
  
  _tabNoStat.addElem(ns, econs,iv1,iv2);


}

///////////////////// Sill ////////////////////////

void CovBase::makeSillNoStatDb(const String &namecol, int ivar, int jvar,const Db* db)
{
  if (!_checkSill(ivar,jvar)) return;
  _makeElemNoStat(EConsElem::SILL, ivar, jvar,nullptr,db, namecol);
  _cor->checkAndManageNoStatDb(db,namecol);
}

void CovBase::makeSillNoStatFunctional(const AFunctional  *func, int ivar, int jvar)
{
  if (!_checkSill(ivar,jvar)) return;
  _makeElemNoStat(EConsElem::SILL, ivar, jvar,func);

}
  
void CovBase::makeSillStationary(int ivar, int jvar)
{
  if (!_checkSill(ivar,jvar)) return;
  if(_tabNoStat.removeElem(EConsElem::SILL, ivar,jvar) == 0)
  {
    messerr("This parameter was already stationary!");
  }
}

/////////////////////////// Check functions ////////////////////:

bool CovBase::_checkSill(int ivar, int jvar) const
{
  int nvar = getNVariables();
  if ((ivar > nvar) || (jvar > nvar))
  {
    messerr("Your model has only %d variables.",nvar);
    return false;
  }
  return true;
}

bool CovBase::_checkDims(int idim, int jdim) const
{
  int ndim = getNDim();
  if ((idim > ndim) || (jdim > ndim))
  {
    messerr("Your model is only in dimension %d.",ndim);
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

double CovBase::getValue(const EConsElem &econs,int iv1,int iv2) const
{
  double val = _cor->getValue(econs,iv1,iv2);
  if (val == TEST)
  {
    if (econs == EConsElem::SILL)
    return getSill(iv1,iv2);
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
    VectorDouble result(coords[0].size(),getValue(econs,iv1,iv2));
    _tabNoStat.informCoords(coords,econs,iv1,iv2,result);
    return result;
  }
 
  return _cor->informCoords(coords,econs,iv1,iv2);
 
}


void CovBase::informMeshByMeshForSills(const AMesh* amesh) const
{
   _tabNoStat.informMeshByMesh(amesh,EConsElem::SILL);
}

void CovBase::informMeshByApexForSills(const AMesh* amesh) const
{
   _tabNoStat.informMeshByApex(amesh,EConsElem::SILL);
}

void CovBase::informDbInForSills(const Db* dbin) const
{
   _tabNoStat.informDbIn(dbin,EConsElem::SILL);
}

void CovBase::informDbOutForSills(const Db* dbout) const
{
  _tabNoStat.informDbOut(dbout,EConsElem::SILL);
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
  if (! isNoStat()) return;
  double val1, val2;

  const auto paramsnostat = _tabNoStat.getTable();
  // Loop on the elements that can be updated one-by-one

  for (const auto &e : paramsnostat)
  {
    EConsElem type = e.first.getType();
    e.second->getValuesOnDb( icas1, iech1, &val1, icas2, iech2, &val2);

    if (type == EConsElem::SILL)
    {
      int iv1 = e.first.getIV1();
      int iv2 = e.first.getIV2();
      setSill(iv1, iv2, sqrt(val1 * val2));
    }  
  }
  _cor->updateCovByPoints(icas1, iech1 , icas2, iech2); 
}


void CovBase::updateCovByMesh(int imesh,bool aniso) const
{
  // If no non-stationary parameter is defined, simply skip
  if (! isNoStat()) return;

  // Loop on the elements that can be updated one-by-one
  if (!aniso)
  {
    const auto paramsnostat = _tabNoStat.getTable();
    for (const auto &e : paramsnostat)
    {
      EConsElem type = e.first.getType();
      if (type == EConsElem::SILL)
      {
        double sill = e.second->getValueOnMeshByApex(imesh);
        int iv1 = e.first.getIV1();
        int iv2 = e.first.getIV2();
        setSill(iv1, iv2, sill);
      }
    }
  }
 _cor->updateCovByMesh(imesh,aniso);
}

void CovBase::makeStationary()
{
  _cor->makeStationary();
}

void CovBase::_manage(const Db* db1,const Db* db2) const
{
  if (db1!=nullptr)
    informDbIn(db1);
  if (db2!=nullptr)
    informDbOut(db2);
  _cor->manage(db1,db2);
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
void CovBase::_addEvalCovMatBiPointInPlace(MatrixSquareGeneral &mat,
                                          const SpacePoint &p1,
                                          const SpacePoint &p2,
                                          const CovCalcMode *mode) const
{
  int nvar = getNVariables();
  for (int ivar = 0; ivar < nvar; ivar++)
    for (int jvar = 0; jvar < nvar; jvar++)
    {
      double cor = _cor->eval(p1,p2,ivar,jvar,mode);
      mat.addValue(ivar, jvar, _sill.getValue(ivar, jvar) * cor);
    }
}

void CovBase::_optimizationPostProcess() const 
{
  _cor->optimizationPostProcess();
}

void CovBase::_updateFromContext()
{
  _cor->updateFromContext();
}