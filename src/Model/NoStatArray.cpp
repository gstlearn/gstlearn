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
#include "Model/NoStatArray.hpp"
#include "Basic/AException.hpp"
#include "Matrix/MatrixRectangular.hpp"
#include "Basic/VectorHelper.hpp"
#include "Basic/Tensor.hpp"
#include "Basic/Utilities.hpp"
#include "Basic/String.hpp"
#include "Basic/File.hpp"
#include "Covariances/CovAniso.hpp"
#include "Model/ANoStat.hpp"
#include "Db/Db.hpp"
#include "Db/DbGrid.hpp"
#include "Mesh/MeshETurbo.hpp"
#include "geoslib_old_f.h"

#include <math.h>

NoStatArray::NoStatArray()
: ANoStat(),
  _dbnostat(nullptr),
  _tab()
{
}

NoStatArray::NoStatArray(const VectorString& codes, const Db* dbnostat)
: ANoStat(codes),
  _dbnostat(dbnostat),
  _tab()
{
  if (! _checkValid())
  {
    messerr("Error in the Definition of Non-Stationarity Parameters");
    messerr("The 'NoStat' element is not valid");
  }
}

NoStatArray::NoStatArray(const NoStatArray &m)
: ANoStat(m),
  _dbnostat(m._dbnostat),
  _tab(m._tab)
{

}

NoStatArray& NoStatArray::operator= (const NoStatArray &m)
{
  if (this != &m)
  {
    ANoStat::operator=(m);
    _dbnostat = m._dbnostat;
    _tab = m._tab;
  }
  return *this;
}

NoStatArray::~NoStatArray()
{
}

bool NoStatArray::_checkValid() const
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

int NoStatArray::attachToMesh(const AMesh* mesh, bool verbose) const
{
  if (_dbnostat == nullptr)
  {
    messerr("dbNoStat must be defined beforehand");
    return 1;
  }

  // Create the array of coordinates

  ANoStat::attachToMesh(mesh,verbose);
  int nmeshes = mesh->getNMeshes();
  VectorDouble loctab(nmeshes,0);
  VectorVectorDouble coords = mesh->getAllCenterCoordinates();

  // Create the internal array

  int npar = getNoStatElemNumber();
  _tab.reset(nmeshes, npar);

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

void NoStatArray::detachFromMesh() const
{
  ANoStat::detachFromMesh();
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
int NoStatArray::attachToDb(Db* db, int icas, bool verbose) const
{
  if (db == nullptr) return 0;

  // Preliminary checks
  if (_dbnostat == nullptr)
  {
    messerr("Dbnostat must be defined beforehand");
    return 1;
  }

  // Store the reference to the Db
  ANoStat::attachToDb(db,icas,verbose);

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

void NoStatArray::detachFromDb(Db* db, int icas) const
{
  if (db == nullptr) return;
  if (db == _dbnostat) return;
  ANoStat::detachFromDb(db,icas);
  db->deleteColumnsByLocator(ELoc::NOSTAT);
}

/**
 * Check if the non-stationary values defined or not
 * @param icas Type of information (0: meshing; 1: Dbin; 2: Dbout)
 * @return
 */
bool NoStatArray::isEmpty(int icas) const
{

  // Dispatch

  if (icas == 0)
  {
    if (_tab.empty()) return true;
  }
  if (icas == 1)
  {
    if (_dbin == nullptr) return true;
  }
  if (icas == 2)
  {
    if (_dbout == nullptr) return true;
  }
  return false;
}

/**
 * Returns the value of a non-stationary parameter at a target sample
 * @param type  Type of non-stationary element (EConsElem)
 * @param icas  Additional identifier (0 for Meshing; 1 for Dbin; 2 for Dbout)
 * @param rank  Rank of the target (in Meshing (0); in Dbin (1) or in Dbout (2)
 * @param icov  Rank of the Covariance
 * @param iv1   Rank of the first variable (optional)
 * @param iv2   Rank of the second variable (optional)
 * @param igrf  Rank of the GRF
 * @return
 */
double NoStatArray::getValue(const EConsElem &type,
                             int icas,
                             int rank,
                             int icov,
                             int iv1,
                             int iv2,
                             int igrf) const
{
  if (! _isValid(icas, rank)) return TEST;
  int ipar = getRank(type, icov, iv1, iv2, igrf);
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
double NoStatArray::getValueByParam(int ipar, int icas, int rank) const
{
  if (! _isValid(icas, rank)) return TEST;
  if (icas == 0)
  {

    // From Meshing

    return _tab(rank, ipar);
  }
  else if (icas == 1)
  {

    // From Dbin

    return _dbin->getFromLocator(ELoc::NOSTAT, rank, ipar);
  }
  else if (icas == 2)
  {

    // From Dbout

    return _dbout->getFromLocator(ELoc::NOSTAT, rank, ipar);
  }
  else
  {
    my_throw("Invalid argument 'icas'");
  }
  return 0.;
}

String NoStatArray::_displayStats(int ipar, int icas) const
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
    if (_dbin == nullptr) return sstr.str();
    vec = _dbin->getColumnByLocator(ELoc::NOSTAT,ipar,true);
  }
  else
  {
    if (_dbout == nullptr) return sstr.str();
    vec = _dbout->getColumnByLocator(ELoc::NOSTAT,ipar,true);
  }

  // Produce the statistics
  if (vec.size() > 0)
    sstr << toVector("Non-stationary Parameter",vec);

  return sstr.str();
}

String NoStatArray::_displayStats(int icas) const
{
  std::stringstream sstr;

  for (int ipar = 0; ipar < getNoStatElemNumber(); ipar++)
    sstr << _displayStats(ipar, icas);

  return sstr.str();
}

String NoStatArray::toString(const AStringFormat* strfmt) const
{
  std::stringstream sstr;
  sstr << ANoStat::toString(strfmt);

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
int NoStatArray::_informField(int ipar,
                              const VectorVectorDouble& coords,
                              VectorDouble& tab,
                              bool verbose) const
{
  // Identify the attribute in the Db

  int iatt = db_attribute_identify(_dbnostat, ELoc::NOSTAT, ipar);
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
    if (migrate_grid_to_coor(dbgrid, iatt, coords, tab)) return 1;
  }
  else
  {
    if (expand_point_to_coor(_dbnostat, iatt, coords, tab)) return 1;
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
