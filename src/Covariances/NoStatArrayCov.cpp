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
#include "Covariances/NoStatArrayCov.hpp"

#include "Calculators/CalcMigrate.hpp"
#include "Basic/AException.hpp"
#include "Matrix/MatrixRectangular.hpp"
#include "Basic/VectorHelper.hpp"
#include "Basic/Utilities.hpp"
#include "Basic/String.hpp"
#include "Covariances/CovAniso.hpp"
#include "Covariances/ANoStatCov.hpp"
#include "Db/Db.hpp"
#include "Db/DbGrid.hpp"
#include "Mesh/MeshETurbo.hpp"
#include "Mesh/AMesh.hpp"

NoStatArrayCov::NoStatArrayCov()
: ANoStatCov(),
  _dbnostat(nullptr),
  _tab()
{
}

NoStatArrayCov::NoStatArrayCov(const VectorString& codes, const Db* dbnostat)
: ANoStatCov(codes),
  _dbnostat(dbnostat),
  _tab()
{
  if (! _checkValid())
  {
    messerr("Error in the Definition of Non-Stationarity Parameters");
    messerr("The 'NoStat' element is not valid");
  }
}

NoStatArrayCov::NoStatArrayCov(const NoStatArrayCov &m)
: ANoStatCov(m),
  _dbnostat(m._dbnostat),
  _tab(m._tab)
{

}

NoStatArrayCov& NoStatArrayCov::operator= (const NoStatArrayCov &m)
{
  if (this != &m)
  {
    ANoStatCov::operator=(m);
    _dbnostat = m._dbnostat;
    _tab = m._tab;
  }
  return *this;
}

NoStatArrayCov::~NoStatArrayCov()
{
}

bool NoStatArrayCov::_checkValid() const
{
  // Get the number of non-stationary parameters from codes
  int nparcd = getNoStatElemNumber();
  if (nparcd <= 0) return false;

  // Get the number of non-stationary parameters from Dbnostat
  if (_dbnostat == nullptr) return false;
  int npardb = _dbnostat->getLocatorNumber(ELoc::NOSTAT);
  if (npardb < nparcd)
  {
    messerr("The Non-Stationary codes require %d parameters", nparcd);
    messerr("The Dbnostat only contains %d ELoc::NOSTAT attributes", npardb);
    return false;
  }
  return true;
}

int NoStatArrayCov::attachToMesh(const AMesh* mesh, bool center, bool verbose) const
{
  if (_dbnostat == nullptr)
  {
    messerr("dbNoStat must be defined beforehand");
    return 1;
  }

  // Create the array of coordinates

  ANoStatCov::attachToMesh(mesh,center,verbose);
  VectorVectorDouble coords;
  
  if (center)
  {
    coords = mesh->getAllCenterCoordinates();
  }
  else 
  {
    coords = mesh->getAllCoordinates();
  }

  int npoint = coords[0].size();
  VectorDouble loctab(npoint,0);

  // Create the internal array

  int npar = getNoStatElemNumber();
  _tab.reset(npoint, npar);

  /* Evaluate the non-stationary parameters */

  for (int ipar=0; ipar<npar; ipar++)
  {
    // Evaluate the non-stationary attribute at the target points
    if (_informField(ipar, coords, loctab, verbose)) return 1;

    // Store the local vector within the Matrix
    _tab.setColumn(ipar, loctab);
  }

  return 0;
}

void NoStatArrayCov::detachFromMesh() const
{
  ANoStatCov::detachFromMesh();
  _tab.reset(0,0);
}

/**
 * Attaching the current Non-Stationary parameters to the Db
 * This function creates new fields set to the locator ELoc::NOSTAT.
 * They will be deleted using the detachFromDb() method.
 * @param db      Db where the new ELoc::NOSTAT fields must be added
 * @param icas    1 for Dbin and 2 for Dbout
 * @param verbose Verbose flag
 * @return
 */
int NoStatArrayCov::attachToDb(Db* db, int icas, bool verbose) const
{
  if (db == nullptr) return 0;

  // Preliminary checks
  if (_dbnostat == nullptr)
  {
    messerr("Dbnostat must be defined beforehand");
    return 1;
  }

  // Store the reference to the Db
  ANoStatCov::attachToDb(db,icas,verbose);

  // If the Db to be attached coincides with _dbnostat, do nothing
  if (db == _dbnostat) return 0;

  // Create the array of coordinates

  int nech = db->getSampleNumber(true);
  VectorDouble tab(nech,0);
  VectorVectorDouble coords = db->getAllCoordinates(true);

  /* Identify the non-stationary parameters within data base(s) */

  int npar = getNoStatElemNumber();
  VectorString names = generateMultipleNames("Nostat",npar);
  for (int ipar=0; ipar<npar; ipar++)
  {
    // Evaluate the non-stationary attribute at the target points
    if (_informField(ipar, coords, tab, verbose)) return 1;

    // Store the local vector within the Db as a new field
    db->addColumns(tab,names[ipar],ELoc::UNKNOWN,0,true);
  }

  // Set locators to the newly created variables
  db->setLocators(names,ELoc::NOSTAT);

  return 0;
}

void NoStatArrayCov::detachFromDb(Db* db, int icas) const
{
  if (db == nullptr) return;
  if (db == _dbnostat) return;
  ANoStatCov::detachFromDb(db,icas);
  db->deleteColumnsByLocator(ELoc::NOSTAT);
}

/**
 * Check if the non-stationary values defined or not
 * @param icas Type of information (0: meshing; 1: Dbin; 2: Dbout)
 * @return
 */
bool NoStatArrayCov::isEmpty(int icas) const
{

  // Dispatch

  if (icas == 0)
  {
    if (_tab.empty()) return true;
  }
  if (icas == 1)
  {
    if (_getDbin() == nullptr) return true;
  }
  if (icas == 2)
  {
    if (_getDbout() == nullptr) return true;
  }
  return false;
}

/**
 * Returns the value of a non-stationary parameter at a target sample
 * @param type  Type of non-stationary element (EConsElem)
 * @param icas  Additional identifier (0 for Meshing; 1 for Dbin; 2 for Dbout)
 * @param rank  Rank of the target (in Meshing (0); in Dbin (1) or in Dbout (2)
 * @param iv1   Rank of the first variable (optional)
 * @param iv2   Rank of the second variable (optional)
 * @return
 */
double NoStatArrayCov::getValue(const EConsElem &type,
                             int icas,
                             int rank,
                             int iv1,
                             int iv2) const
{
  if (! _isValid(icas, rank)) return TEST;
  int ipar = getRank(type, iv1, iv2);
  return getValueByParam(ipar, icas, rank);
}

/**
 * Return the value of the non-stationary parameter (ipar) at target (rank)
 * @param ipar  Rank of the non-stationary parameter
 * @param icas  Source definition:
 *              0 : from Meshing (rank: absolute rank to be converted into relative)
 *              1 : from Dbin
 *              2 : from Dbout
 * @param rank  Rank of the target
 * @return
 */
double NoStatArrayCov::getValueByParam(int ipar, int icas, int rank) const
{
  if (!_isValid(icas, rank)) return TEST;
  if (icas == 0)
  {

    // From Meshing

    return _tab(rank, ipar);
  }
  if (icas == 1)
  {

    // From Dbin

    return _getDbin()->getFromLocator(ELoc::NOSTAT, rank, ipar);
  }
  if (icas == 2)
  {

    // From Dbout

    return _getDbout()->getFromLocator(ELoc::NOSTAT, rank, ipar);
  }
  my_throw("Invalid argument 'icas'");
  return 0.;
}

String NoStatArrayCov::_displayStats(int ipar, int icas) const
{
  std::stringstream sstr;

  VectorDouble vec;
  if (icas == 0)
  {
    if (_tab.empty()) return sstr.str();;
    for (int iech = 0; iech < _tab.getNRows(); iech++)
      vec.push_back(_tab.getValue(iech, ipar));
  }
  else if (icas == 1)
  {
    if (_getDbin() == nullptr) return sstr.str();
    vec = _getDbin()->getColumnByLocator(ELoc::NOSTAT,ipar,true);
  }
  else
  {
    if (_getDbout() == nullptr) return sstr.str();
    vec = _getDbout()->getColumnByLocator(ELoc::NOSTAT,ipar,true);
  }

  // Produce the statistics
  if (vec.size() > 0)
    sstr << toVector("Non-stationary Parameter",vec);

  return sstr.str();
}

String NoStatArrayCov::_displayStats(int icas) const
{
  std::stringstream sstr;

  for (int ipar = 0; ipar < getNoStatElemNumber(); ipar++)
    sstr << _displayStats(ipar, icas);

  return sstr.str();
}

String NoStatArrayCov::toString(const AStringFormat* strfmt) const
{
  std::stringstream sstr;
  sstr << ANoStatCov::toString(strfmt);

  AStringFormat sf;
  if (strfmt != nullptr) sf = *strfmt;
  if (sf.getLevel() > 1)
    sstr << _displayStats(0);
  return sstr.str();
}

/**
 * Inform the value of non-stationary parameter 'ipar' at samples
 * defined by their coordinates
 * @param ipar Rank of the non-stationary parameter
 * @param coords Coordinates of the samples (dimension: ndim, npar)
 * @param tab Return array of non-stationary values at samples
 * @param verbose Verbose flag
 * @return
 */
int NoStatArrayCov::_informField(int ipar,
                              const VectorVectorDouble& coords,
                              VectorDouble& tab,
                              bool verbose) const
{
  // Identify the attribute in the Db

  int iatt = _dbnostat->getUIDByLocator(ELoc::NOSTAT, ipar);
  if (iatt < 0)
  {
    messerr("The Non-stationary attribute (%d) is not defined in _dbnostat",
            ipar);
    return 1;
  }

  // Migrate the information from Db onto the Vertex locations

  if (_dbnostat->isGrid())
  {
    const DbGrid* dbgrid = dynamic_cast<const DbGrid*>(_dbnostat);
    if (migrateGridToCoor(dbgrid, iatt, coords, tab)) return 1;
  }
  else
  {
    if (expandPointToCoor(_dbnostat, iatt, coords, tab)) return 1;
  }

  int ndef = VH::countUndefined(tab);
  if (ndef > 0)
  {

    // Calculate local statistics

    double mean = VH::mean(tab);
    if (FFFF(mean))
    {
      messageFlush(getItems(ipar).toString());
      messerr("This Non-Stationary parameter is not valid");
      return 1;
    }

    if (verbose)
    {
      message("For Non-Stationary Parameter (%d), there are %d undefined values\n",
              ipar + 1, ndef);
      message("They have been replaced by its average value (%lf)\n", mean);
    }

    // Modify the TEST values to the mean value

    VH::fillUndef(tab, mean);
  }

  // Printout some statistics (optional)

  if (verbose)
  {
    char str[LONG_SIZE];
    (void) gslSPrintf(str,
                      "Statistics for Non-Stationary Parameter #%d on Mesh",
                      ipar + 1);
    VH::displayStats(str,tab);
  }

  return 0;
}