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
#include "Covariances/CovContext.hpp"
#include "Covariances/NoStatArray.hpp"
#include "Covariances/NoStatFunctional.hpp"
#include "Covariances/TabNoStat.hpp"
#include "Db/Db.hpp"

#include "Space/ASpace.hpp"
#include "geoslib_define.h"

#include <math.h>

ACor::ACor(const CovContext &ctxt)
    : ASpaceObject(ctxt.getASpace()),
      _nvar(ctxt.getNVar()),
      _tabNoStat(nullptr),
      _ctxt(ctxt)
{
  createNoStatTab();
}


ACor::ACor(const ACor &r)
    : ASpaceObject(r),
      _nvar(r._nvar),
      _tabNoStat(r._tabNoStat == nullptr? nullptr:new TabNoStat(*r._tabNoStat))
{
}

ACor& ACor::operator=(const ACor &r)
{
  if (this != &r)
  {
    ASpaceObject::operator=(r);
    _nvar = r._nvar;
    _tabNoStat = r._tabNoStat->clone();
  }
  return *this;
}

ACor::~ACor()
{
  delete _tabNoStat;
}

void ACor::createNoStatTab()
{
  delete _tabNoStat;
  _tabNoStat = _createNoStatTab();
}

TabNoStat* ACor::_createNoStatTab()
{
  return new TabNoStat();
}
double ACor::eval0(int ivar,
                   int jvar,
                   const CovCalcMode* mode) const
{
  DECLARE_UNUSED(ivar,jvar,mode)
  return 1.;
}

void ACor::copyCovContext(const CovContext &ctxt)
{
  _ctxt.copyCovContext(ctxt);
  _copyCovContext(ctxt);
}

/////////////  Functions to attach no stat information on various supports ////////
void ACor::informMeshByMesh(const AMesh* amesh) const
{
  _tabNoStat->informMeshByMesh(amesh);
}
void ACor::informMeshByApex(const AMesh* amesh) const
{
  _tabNoStat->informMeshByApex(amesh);
}
void ACor::informDbIn(const Db* dbin) const
{
  _tabNoStat->informDbIn(dbin);
}
void ACor::informDbOut(const Db* dbout) const
{
  _tabNoStat->informDbOut(dbout);
}


void ACor::setNoStatDbIfNecessary(const Db*& db)
{
  if (_tabNoStat->getDbNoStatRef() == nullptr)
    attachNoStatDb(db);
  if (db == nullptr)
    db = _tabNoStat->getDbNoStatRef();
}

void ACor::makeStationary()
{
  _tabNoStat->clear();
}
int ACor::makeElemNoStat(const EConsElem &econs, int iv1, int iv2,const AFunctional* func, const Db* db, const String& namecol)
{
  std::shared_ptr<ANoStat> ns;
  if (func == nullptr)
  {
    if(!checkAndManageNoStatDb(db,namecol)) return 1;
    ns = std::shared_ptr<ANoStat>(new NoStatArray(db,namecol));
  }
  else 
  {
    ns = std::unique_ptr<ANoStat>(new NoStatFunctional(func));
  }
   return _tabNoStat->addElem(ns, econs,iv1,iv2);
  
}

bool ACor::_checkDims(int idim, int jdim) const
{
  int ndim = getNDim();
  if ((idim > ndim) || (jdim > ndim))
  {
    messerr("Your model is only in dimension %d.",ndim);
    return false;
  }
  return true;
}


// Set of functions to make parameters no stationary (or to make them back stationary).
// There is to types of non stationarities : NoStatDb in which the parameters are read in a
// DbGrid or NoStatFunctional for which you have to provide a function of the coordinates.
// Each parameter can have its own type of No stationarity and its own DbGrid in case
// of NoStatDb. 
// For specifying the NoStat DbGrid, you can first attach it by using attachNoStatDb.
// If not, you have to specify the DbGrid when you make the first parameter non stationary.

void ACor::attachNoStatDb(const Db* db)
{
  _tabNoStat->setDbNoStatRef(db);
}


bool ACor::checkAndManageNoStatDb(const Db*&  db, const String& namecol)
{
 if (_tabNoStat->getDbNoStatRef() == nullptr && db == nullptr)
 {
  messerr("You have to define a Db (with attachNoStatDb or by specifying a Db here)");  
  return false;
 }
  setNoStatDbIfNecessary(db);

 if (db->getUID(namecol)< 0)
 {
    messerr("You have to specified a name of a column of the reference Db");
    return false;
 }
 return true;
}


VectorDouble ACor::informCoords(const VectorVectorDouble& coords, 
                                    const EConsElem& econs,
                                    int iv1,
                                    int iv2) const
{
  VectorDouble result(coords[0].size(),getValue(econs,iv1,iv2));
  _tabNoStat->informCoords(coords,econs,iv1,iv2,result);
  return result;
}