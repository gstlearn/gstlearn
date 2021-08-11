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
#include "Model/NoStatArray.hpp"

#include "geoslib_e.h"

#include "Basic/AException.hpp"
#include "MatrixC/MatrixCRectangular.hpp"
#include "Basic/Vector.hpp"
#include "Basic/Tensor.hpp"
#include "Basic/Utilities.hpp"
#include "Basic/String.hpp"
#include "Covariances/CovAniso.hpp"
#include "Model/ANoStat.hpp"

NoStatArray::NoStatArray()
    : ANoStat(),
      _dbnostat(nullptr),
      _dbin(nullptr),
      _dbout(nullptr),
      _tab()
{
}

NoStatArray::NoStatArray(const VectorString& codes, const Db* dbnostat)
    : ANoStat(codes),
      _dbnostat(dbnostat),
      _dbin(nullptr),
      _dbout(nullptr),
      _tab()
{
  if (! _checkValid())
    my_throw("Error in the Definition of Non-Stationarity Parameters");
}

NoStatArray::NoStatArray(const NoStatArray &m)
    : ANoStat(m),
      _dbnostat(m._dbnostat),
      _dbin(m._dbin),
      _dbout(m._dbout),
      _tab(m._tab)
{

}

NoStatArray& NoStatArray::operator= (const NoStatArray &m)
{
  if (this != &m)
  {
    ANoStat::operator=(m);
    _dbnostat = m._dbnostat;
    _dbin = m._dbin;
    _dbout = m._dbout;
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
  int npardb = _dbnostat->getLocatorNumber(LOC_NOSTAT);
  if (npardb < nparcd)
  {
    messerr("The Non-Stationary codes require %d parameters", nparcd);
    messerr("The Dbnostat only contains %d LOC_NOSTAT attributes", npardb);
    return false;
  }
  return true;
}

int NoStatArray::attachToMesh(const AMesh* mesh, bool verbose) const
{
  double* coorloc[3];

  // Preliminary checks
  if (_dbnostat == nullptr)
  {
    messerr("Dbnostat must be defined beforehand");
    return 1;
  }

  // Create the array of coordinates

  int ndim    = _dbnostat->getNDim();
  int nvertex = mesh->getNApices();
  VectorDouble tab(nvertex,0);
  double* coor = (double *) mem_alloc(sizeof(double) * ndim * nvertex, 1);
  int ecr = 0;
  for (int idim=0; idim<ndim; idim++)
    for (int ip=0; ip<nvertex; ip++,ecr++)
      coor[ecr] = mesh->getApexCoor(ip,idim);
  for (int idim=0; idim<3; idim++)
    coorloc[idim] = (idim < ndim) ? &coor[idim * nvertex] : NULL;

  // Create the internal array

  int npar = getNoStatElemNumber();
  _tab.reset(nvertex, npar);

  /* Evaluate the non-stationary parameters */

  for (int ipar=0; ipar<npar; ipar++)
  {
    // Evaluate the non-stationary attribute at the target points
    if (_informField(ipar, nvertex, coorloc, tab, verbose)) return 1;

    // Store the local vector within the Matrix
    _tab.setColumn(ipar, tab);
  }

  coor = (double *) mem_free((char *) coor);
  return 0;
}

void NoStatArray::detachFromMesh(bool verbose) const
{
  _tab.reset(0,0);
}

/**
 * Attaching the current Non-Stationary parameters to the Db
 * This function creates new fields set to the locator LOC_NOSTAT.
 * They will be deleted using the detachFromDb() method.
 * @param db      Db where the new LOC_NOSTAT fields must be added
 * @param icas    1 for Dbin and 2 for Dbout
 * @param verbose Verbose flag
 * @return
 */
int NoStatArray::attachToDb(Db* db, int icas, bool verbose) const
{
  double* coorloc[3];
  if (db == nullptr) return 0;

  // Preliminary checks
  if (_dbnostat == nullptr)
  {
    messerr("Dbnostat must be defined beforehand");
    return 1;
  }

  // Store the reference to the Db
  if (icas == 1)
    _dbin = db;
  else
    _dbout = db;

  // If the Db to be attached coincides with _dbnostat, do nothing
  if (db == _dbnostat) return 0;

  // Create the array of coordinates

  int ndim = _dbnostat->getNDim();
  int nech = db->getActiveSampleNumber();

  VectorDouble tab(nech,0);
  double* coor = (double *) mem_alloc(sizeof(double) * ndim * nech, 1);
  int ecr = 0;
  for (int idim=0; idim<ndim; idim++)
    for (int iech=0; iech<db->getSampleNumber(); iech++)
    {
      if (! db->isActive(iech)) continue;
      coor[ecr++] = db->getCoordinate(iech,idim);
    }
  for (int idim=0; idim<3; idim++)
    coorloc[idim] = (idim < ndim) ? &coor[idim * nech] : NULL;

  /* Identify the non-stationary parameters within data base(s) */

  int npar = getNoStatElemNumber();
  VectorString names = generateMultipleNames("Nostat",npar);
  for (int ipar=0; ipar<npar; ipar++)
  {
    // Evaluate the non-stationary attribute at the target points
    if (_informField(ipar, nech, coorloc, tab, verbose)) return 1;

    // Store the local vector within the Db as a new field
    db->addFields(tab,names[ipar],LOC_UNKNOWN,true);
  }

  // Set locators to the newly created variables
  db->setLocator(names,LOC_NOSTAT);

  coor = (double *) mem_free((char *) coor);
  return 0;
}

void NoStatArray::detachFromDb(Db* db, int icas, bool verbose) const
{
  if (db == nullptr) return;
  if (db == _dbnostat) return;
  db->deleteFieldByLocator(LOC_NOSTAT);
  if (icas == 1)
    _dbin = nullptr;
  else
    _dbout = nullptr;
}

/**
 * Returns the value of a non-stationary parameter at a target sample
 * @param igrf  Rank of the GRF
 * @param icov  Rank of the Covariance
 * @param type  Type of non-stationary element
 * @param iv1   Rank of the first variable (optional)
 * @param iv2   Rank of the second variable (optional)
 * @param icas  Additional identifier (0 for Meshing; 1 for Dbin; 2 for Dbout)
 * @param rank  Rank of the target (in Meshing (0); in Dbin (1) or in Dbout (2)
 * @return
 */
double NoStatArray::getValue(int igrf,
                             int icov,
                             ENUM_CONS type,
                             int iv1,
                             int iv2,
                             int icas,
                             int rank) const
{
  int ipar = getRank(igrf, icov, type, iv1, iv2);
  return getValue(ipar, icas, rank);
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
    if (_tab.isEmpty()) return true;
  }
  if (icas == 1)
  {
    if (_dbin == (Db *) NULL) return true;
  }
  if (icas == 2)
  {
    if (_dbout == (Db *) NULL) return true;
  }
  return false;
}

/**
 * Return the value of the non-stationary parameter (ipar) at target (rank)
 * @param ipar  Rank of the non-stationary parameter
 * @param icas  Additional identifier
 * @param rank  Rank of the target
 * @return
 */
double NoStatArray::getValue(int ipar, int icas, int rank) const
{
  if (ipar < 0)
    my_throw("Invalid rank when searching for Non-stationary parameter");
  if (isEmpty(icas))
    my_throw("The Non-Stationary storage must be defined beforehand");

  // Dispatch

  if (icas == 0)
  {

    // From Meshing

    if (rank < 0 || rank > _tab.getNRows()) return TEST;
    return _tab(rank, ipar);
  }
  else if (icas == 1)
  {

    // From Dbin

    if (_dbin == nullptr) return TEST;
    if (rank < 0 || rank > _dbin->getSampleNumber()) return TEST;
    return _dbin->getFromLocator(LOC_NOSTAT, rank, ipar);
  }
  else if (icas == 2)
  {

    // From Dbout

    if (_dbout == nullptr) return TEST;
    if (rank < 0 || rank > _dbout->getSampleNumber()) return TEST;
    return _dbout->getFromLocator(LOC_NOSTAT, rank, ipar);
  }
  else
  {
    my_throw("Invalid argument 'icas'");
  }
  return 0.;
}

double NoStatArray::_interpolate(int ipar,
                                 int icas1,
                                 int iech1,
                                 int icas2,
                                 int iech2) const
{
  double val1 = getValue(ipar, icas1, iech1);
  double val2 = getValue(ipar, icas2, iech2);

  if (! FFFF(val1) && ! FFFF(val2))
    return sqrt(val1 * val2);

  else if (! FFFF(val1))
    return val2;
  else
    return val1;
}

/**
 * Get the information from the storage in Dbin and/or Dbout
 * @param ipar  Rank of the non-stationary parameter
 * @param icas1 Type of first Db: 1 for Input; 2 for Output
 * @param iech1 Rank of the first sample (in Dbin)
 * @param icas2 Type of first Db: 1 for Input; 2 for Output
 * @param iech2 Rank of the second sample (in Dbout)
 * @param val1  Returned value at first sample
 * @param val2  Returned value at the second sample
 */
void NoStatArray::_getInfoFromDb(int ipar,
                                 int icas1,
                                 int iech1,
                                 int icas2,
                                 int iech2,
                                 double *val1,
                                 double *val2) const
{
  *val1 = getValue(ipar, icas1, iech1);
  *val2 = getValue(ipar, icas2, iech2);

  if (FFFF(*val1) && FFFF(*val2)) return;

  if (! FFFF(*val1))
    *val2 = *val1;
  if (! FFFF(*val2))
    *val1 = *val2;
}

/**
 * Update the Model according to the Non-stationary parameters
 * @param model Model to be patched
 * @param icas1 Type of first Db: 1 for Input; 2 for Output
 * @param iech1 Rank of the target within Db1 (or -1)
 * @param icas2 Type of first Db: 1 for Input; 2 for Output
 * @param iech2 Rank of the target within Dbout (or -2)
 */
void NoStatArray::updateModel(Model* model,
                              int icas1,
                              int iech1,
                              int icas2,
                              int iech2) const
{
  double val1, val2;

  // If no non-stationary parameter is defined, simply skip
  const ANoStat* nostat = model->getNoStat();
  if (nostat == nullptr) return;

  // Loop on the elements that can be updated one-by-one

  for (int ipar = 0; ipar < nostat->getNoStatElemNumber(); ipar++)
  {
    int icov = nostat->getICov(ipar);
    int type = nostat->getType(ipar);

    if (type == CONS_SILL)
    {
      _getInfoFromDb(ipar, icas1, iech1, icas2, iech2, &val1, &val2);
      int iv1  = nostat->getIV1(ipar);
      int iv2  = nostat->getIV2(ipar);
      model->setSill(icov, iv1, iv2, sqrt(val1 * val2));
    }
    else if (type == CONS_PARAM)
    {
      _getInfoFromDb(ipar, icas1, iech1, icas2, iech2, &val1, &val2);
      model->getCova(icov)->setParam(0.5 * (val1 + val2));
    }
  }

  // Loop on the other parameters (Anisotropy) that must be processed globally

  for (int icov = 0; icov < model->getCovaNumber(); icov++)
  {
    if (! isDefinedforAnisotropy(-1, icov)) continue;
    CovAniso* cova = model->getCova(icov);

    VectorDouble angle0(cova->getAnisoAngles());
    VectorDouble angle1(angle0);
    VectorDouble angle2(angle0);

    VectorDouble scale0(cova->getScales());
    VectorDouble scale1(scale0);
    VectorDouble scale2(scale0);

    VectorDouble range0(cova->getRanges());
    VectorDouble range1(range0);
    VectorDouble range2(range0);

    // Define the angles (for all space dimensions)
    bool flagRot = false;
    if (isDefined(-1, icov, CONS_ANGLE, -1, -1))
    {
      flagRot = true;
      for (int idim = 0; idim < model->getDimensionNumber(); idim++)
      {
        if (isDefined(-1, icov, CONS_ANGLE, idim, 0))
        {
          int ipar = getRank(-1, icov, CONS_ANGLE, idim, -1);
          if (ipar < 0) continue;
          _getInfoFromDb(ipar, icas1, iech1, icas2, iech2,
                         &angle1[idim], &angle2[idim]);
        }
      }
    }

    // Define the Theoretical ranges (for all space dimensions)

    bool flagScale = false;
    if (isDefined(-1, icov, CONS_SCALE, -1, -1))
    {
      flagScale = true;
      for (int idim = 0; idim < model->getDimensionNumber(); idim++)
      {
        if (isDefined(-1, icov, CONS_SCALE, idim, -1))
        {
          int ipar = getRank(-1, icov, CONS_SCALE, idim, -1);
          if (ipar < 0) continue;
          _getInfoFromDb(ipar, icas1, iech1, icas2, iech2,
                         &scale1[idim], &scale2[idim]);
        }
      }
    }

    // Define the Practical ranges (for all space dimensions)

    bool flagRange = false;
    if (isDefined(-1, icov, CONS_RANGE, -1, -1))
    {
      flagRange = true;
      for (int idim = 0; idim < model->getDimensionNumber(); idim++)
      {
        if (isDefined(-1, icov, CONS_RANGE, idim, -1))
        {
          int ipar = getRank(-1, icov, CONS_RANGE, idim, -1);
          if (ipar < 0) continue;
          _getInfoFromDb(ipar, icas1, iech1, icas2, iech2,
                         &range1[idim], &range2[idim]);
        }
      }
    }

    // Create and exploit the Tensor for First location
    if (flagRot || flagRange || flagScale)
    {
      if (flagRot)   cova->setAnisoAngles(angle1);
      if (flagRange) cova->setRanges(range1);
      if (flagScale) cova->setScales(scale1);
      const MatrixCSGeneral direct1 = cova->getAniso().getTensorDirect();

      // Create and exploit the Tensor for the second location
      if (flagRot)   cova->setAnisoAngles(angle2);
      if (flagRange) cova->setRanges(range2);
      if (flagScale) cova->setScales(scale2);

      // Build the new Tensor (as average of tensors at end-points)
      Tensor tensor = cova->getAniso();
      MatrixCSGeneral direct = tensor.getTensorDirect();
      direct.linearCombination(0.5, 0.5, direct1);
      tensor.setTensorDirect(direct);
      cova->setAniso(tensor);
    }
  }
}

/**
 * Update the Model according to the Non-stationary parameters
 * @param model Model to be patched
 * @param vertex Rank of the target vertex
 */
void NoStatArray::updateModel(Model* model, int vertex) const
{
  // If no non-stationary parameter is defined, simply skip
  const ANoStat* nostat = model->getNoStat();
  if (nostat == nullptr) return;

  // Loop on the elements that can be updated one-by-one

  for (int ipar = 0; ipar < model->getNoStatElemNumber(); ipar++)
  {
    int icov = nostat->getICov(ipar);
    ENUM_CONS type = nostat->getType(ipar);

    if (type == CONS_SILL)
    {
      double sill = getValue(ipar, 0, vertex);
      int iv1  = nostat->getIV1(ipar);
      int iv2  = nostat->getIV2(ipar);
      model->setSill(icov, iv1, iv2, sill);
    }
  }

  // Loop on the other parameters (Anisotropy) that must be processed globally

  for (int icov = 0; icov < model->getCovaNumber(); icov++)
  {
    if (! isDefinedforAnisotropy(-1, icov)) continue;
    CovAniso* cova = model->getCova(icov);

    VectorDouble angle(cova->getAnisoAngles());
    VectorDouble scale(cova->getScales());
    VectorDouble range(cova->getRanges());

    // Define the angles (for all space dimensions)
    bool flagRot = false;
    if (isDefined(-1, icov, CONS_ANGLE, -1, -1))
    {
      flagRot = true;
      for (int idim = 0; idim < model->getDimensionNumber(); idim++)
      {
        if (isDefined(-1, icov, CONS_ANGLE, idim, 0))
        {
          int ipar = getRank(-1, icov, CONS_ANGLE, idim, -1);
          if (ipar < 0) continue;
          angle[idim] = getValue(ipar, 0, vertex);
        }
      }
    }

    // Define the Theoretical ranges (for all space dimensions)

    bool flagScale = false;
    if (isDefined(-1, icov, CONS_SCALE, -1, -1))
    {
      flagScale = true;
      for (int idim = 0; idim < model->getDimensionNumber(); idim++)
      {
        if (isDefined(-1, icov, CONS_SCALE, idim, -1))
        {
          int ipar = getRank(-1, icov, CONS_SCALE, idim, -1);
          if (ipar < 0) continue;
          scale[idim] = getValue(ipar, 0, vertex);
        }
      }
    }

    // Define the Practical ranges (for all space dimensions)

    bool flagRange = false;
    if (isDefined(-1, icov, CONS_RANGE, -1, -1))
    {
      flagRange = true;
      for (int idim = 0; idim < model->getDimensionNumber(); idim++)
      {
        if (isDefined(-1, icov, CONS_RANGE, idim, -1))
        {
          int ipar = getRank(-1, icov, CONS_RANGE, idim, -1);
          if (ipar < 0) continue;
          range[idim] = getValue(ipar, 0, vertex);
        }
      }
    }

    // Exploit the Anisotropy
    if (flagRot)   cova->setAnisoAngles(angle);
    if (flagRange) cova->setRanges(range);
    if (flagScale) cova->setScales(scale);
  }
}

String NoStatArray::displayStats(int ipar, int icas) const
{
  std::stringstream sstr;

  VectorDouble vec;
  if (icas == 0)
  {
    if (_tab.isEmpty()) return sstr.str();;
    for (int iech = 0; iech < _tab.getNRows(); iech++)
      vec.push_back(_tab.getValue(iech, ipar));
  }
  else if (icas == 1)
  {
    if (_dbin == (Db *) NULL) return sstr.str();
    vec = _dbin->getFromLocatorVector(LOC_NOSTAT,ipar,true);
  }
  else
  {
    if (_dbout == (Db *) NULL) return sstr.str();
    vec = _dbout->getFromLocatorVector(LOC_NOSTAT,ipar,true);
  }

  // Produce the statistics
  if (vec.size() > 0)
    sstr << toVector("Non-stationary Parameter",vec);

  return sstr.str();
}

String NoStatArray::displayStats(int icas) const
{
  std::stringstream sstr;

  for (int ipar = 0; ipar < getNoStatElemNumber(); ipar++)
    sstr << displayStats(ipar, icas);

  return sstr.str();
}

String NoStatArray::toString(int level) const
{
  std::stringstream sstr;
  sstr << ANoStat::toString(level);
  if (level > 0)
    sstr << displayStats(0);
  return sstr.str();
}

int NoStatArray::_informField(int ipar,
                              int nech,
                              double* coor[3],
                              VectorDouble& tab,
                              bool verbose) const
{
  // Identify the attribute in the Db

  int iatt = db_attribute_identify(_dbnostat, LOC_NOSTAT, ipar);
  if (iatt < 0)
  {
    messerr("The Non-stationary attribute (%d) is not defined in Dbnostat",
            ipar);
    return 1;
  }

  // Migrate the information from Db onto the Vertex locations

  if (is_grid(_dbnostat))
  {
    if (migrate_grid_to_coor(_dbnostat, iatt, nech, coor[0], coor[1],
                             coor[2], tab.data())) return 1;
  }
  else
  {
    if (expand_point_to_coor(_dbnostat, iatt, nech, coor[0], coor[1],
                             coor[2], tab.data())) return 1;
  }

  int ndef = ut_vector_count_undefined(tab);
  if (ndef > 0)
  {

    // Calculate local statistics

    double mean = ut_vector_mean(tab);
    if (FFFF(mean))
    {
      messageFlush(getItems(ipar).toString());
      messerr("This Non-Stationary parameter is not valid");
      return 1;
    }

    message("For Non-Stationary Parameter (%d), there are some undefined values (%d)\n",
        ipar + 1, ndef);
    message("They have been replaced by its average value (%lf)\n", mean);

    // Modify the TEST values to the mean value

    for (int ip = 0; ip < nech; ip++)
    {
      if (FFFF(tab[ip])) tab[ip] = mean;
    }
  }

  // Printout some statistics (optional)

  if (verbose)
    ut_vector_display_stats(stringCompose("Statistics for Non-Stationary Parameter #%d on Mesh",
                                          ipar + 1),tab);

  return 0;
}
