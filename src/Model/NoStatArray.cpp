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
      _dbin(nullptr),
      _attIn(),
      _dbout(nullptr),
      _attOut(),
      _tab()
{
}

NoStatArray::NoStatArray(const VectorString& codes)
    : ANoStat(codes),
      _dbin(nullptr),
      _attIn(),
      _dbout(nullptr),
      _attOut(),
      _tab()
{
}

NoStatArray::NoStatArray(const NoStatArray &m)
    : ANoStat(m),
      _dbin(m._dbin),
      _attIn(m._attIn),
      _dbout(m._dbout),
      _attOut(m._attOut),
      _tab(m._tab)
{

}

NoStatArray& NoStatArray::operator= (const NoStatArray &m)
{
  if (this != &m)
  {
    ANoStat::operator=(m);
    _dbin = m._dbin;
    _attIn = m._attIn;
    _dbout = m._dbout;
    _attOut = m._attOut;
    _tab = m._tab;
  }
  return *this;
}

NoStatArray::~NoStatArray()
{

}

int NoStatArray::attachMesh(Db* db, const AMesh* mesh, bool verbose)
{
  double* coorloc[3];

  int error = 1;

  // Preliminary checks

  if (db == (Db *) NULL)
  {
    messerr("Db must be defined");
    return 1;
  }
  int ndim = db->getNDim();
  int flag_grid = is_grid(db);
  int npar = getNoStatElemNumber();

  // Local array

  int nvertex = mesh->getNApices();
  VectorDouble tab(nvertex,0);
  double* coor = (double *) mem_alloc(sizeof(double) * ndim * nvertex, 1);

  // Create the array of coordinates

  int ecr = 0;
  for (int idim=0; idim<ndim; idim++)
    for (int ip=0; ip<nvertex; ip++,ecr++)
      coor[ecr] = mesh->getApexCoor(ip,idim);

  for (int idim=0; idim<3; idim++)
    coorloc[idim] = (idim < ndim) ? &coor[idim * nvertex] : NULL;

  // Create the internal array

  _tab.reset(nvertex, npar);

  /* Evaluate the non-stationary parameters */

  for (int ipar=0; ipar<npar; ipar++)
  {
    // Identify the attribute in the Db

    int iatt = db_attribute_identify(db,LOC_NOSTAT,ipar);
    if (iatt < 0)
    {
      messerr("The Non-stationary attribute (%d) is not defined in Db", ipar);
      goto label_end;
    }

    // Migrate the information from Db onto the Vertex locations

    if (flag_grid)
    {
      if (migrate_grid_to_coor(db,iatt,nvertex,
                               coorloc[0],coorloc[1],coorloc[2],
                               tab.data())) goto label_end;
    }
    else
    {
      if (expand_point_to_coor(db,iatt,nvertex,
                               coorloc[0],coorloc[1],coorloc[2],
                               tab.data())) goto label_end;
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

      for (int ip = 0; ip < nvertex; ip++)
      {
        if (FFFF(tab[ip])) tab[ip] = mean;
      }
    }

    // Store the local vector within the Matrix

    _tab.setColumn(ipar, tab);

    // Printout some statistics (optional)

    if (verbose)
      ut_vector_display_stats(stringCompose("Statistics for Non-Stationary Parameter #%d on Mesh",ipar+1),tab);
  }

  // Set the error return code

  error = 0;

  label_end:
  coor = (double *) mem_free((char *) coor);
  return error;
}

/**
 * Attaching the current Non-Stationary parameters to the Db
 * @param db      Db containing the LOC_NOSTAT fields
 * @param icas    Type of Db: 1 for 'dbin' and 2 for 'dbout
 * @param verbose Verbose flag
 * @return
 */
int NoStatArray::attachDb(Db* db, int icas, bool verbose)
{
  int npar = getNoStatElemNumber();

  // Preliminary checks

  if (db == (Db *) NULL)
  {
    messerr("Db must be defined");
    return 1;
  }
  if (icas != 1 && icas != 2)
  {
    messerr("Argument 'icas' should be 1 or 2");
    return 1;
  }

  // Save the pointers to the Data Bases

  if (icas == 1)
  {
    _dbin = db;
    _attIn.clear();
  }
  else
  {
    _dbout = db;
    _attOut.clear();
  }

  /* Identify the non-stationary parameters within data base(s) */

  for (int ipar=0; ipar<npar; ipar++)
  {
    // Identify the attribute in the input Db

    int attIn = db_attribute_identify(db, LOC_NOSTAT, ipar);
    if (attIn < 0)
    {
      messerr("The Non-stationary attribute (%d) is not defined in Db", ipar+1);
      return 1;
    }
    _attIn.push_back(attIn);

    // Printout some statistics (optional)

    if (verbose)
    {
      VectorInt atts;
      atts.push_back(attIn);
      VectorString opers;
      opers.push_back("mean");
      db_stats_print(db,atts,opers,0,0,
                     stringCompose("Statistics for Non-Stationary Parameter #%d on Mesh",
                                   ipar+1));
    }
  }

  return 0;
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
 * @param ipar Rank of the Non-Stationary parameter
 * @param icas Type of information (0: meshing; 1: Dbin; 2: Dbout; -1: Any)
 * @return
 */
bool NoStatArray::isEmpty(int ipar, int icas) const
{

  // Dispatch

  if (icas < 0 || icas == 0)
  {
    if (_tab.isEmpty())
      return true;
    if (ipar >= 0 && ipar >= _tab.getNCols())
      return true;
  }
  if (icas < 0 || icas == 1)
  {
    if (_dbin == (Db *) NULL)
      return true;
    if (_attIn.empty())
      return true;
    if (ipar >= 0 && ipar >= (int) _attIn.size())
      return true;
  }
  if (icas < 0 || icas == 2)
  {
    if (_dbout == (Db *) NULL)
      return true;
    if (_attOut.empty())
      return true;
    if (ipar >= 0 && ipar >= (int) _attOut.size())
      return true;
  }
  return false;
}

/**
 * Return the value of the non-stationary parameter (ipar) at target (rank)
 * @param ipar  Rank of the non-stationary parameter
 * @param icas  Additional identifier (useless here)
 * @param rank  Rank of the target
 * @return
 */
double NoStatArray::getValue(int ipar, int icas, int rank) const
{
  if (ipar < 0)
    my_throw("Invalid rank when searching for Non-stationary parameter");
  if (isEmpty(ipar, icas))
    my_throw("The Non-Stationary storage must be defined beforehand");

  // Dispatch

  if (icas == 0)
  {

    // From Meshing

    if (rank < 0 || rank > _tab.getNRows())
      my_throw("Invalid Vertex index when searching for Non-stationary parameter");
    return _tab(rank, ipar);
  }
  else if (icas == 1)
  {

    // From Dbin

    if (rank < 0 || rank > get_NECH(_dbin))
      my_throw(
          "Error: Invalid Rank in Dbin when searching for NonStat parameter");
    return get_ARRAY(_dbin, rank, _attIn[ipar]);
  }
  else if (icas == 2)
  {

    // From Dbout

    if (rank < 0 || rank > get_NECH(_dbout))
      my_throw(
          "Error: Invalid Rank in Dbout when searching for NonStat parameter");
    return get_ARRAY(_dbout, rank, _attOut[ipar]);
  }
  else
  {
    my_throw("Invalid argument 'icas'");
  }
  return 0.;
}

double NoStatArray::_interpolate(int ipar, int iech1, int iech2)
{
  double val1 = TEST;
  double val2 = TEST;
  if (_dbin != (Db *) NULL)
    val1 = getValue(ipar, 1, iech1);
  if (_dbout != (Db *) NULL)
    val2 = getValue(ipar, 2, iech2);

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
 * @param iech1 Rank of the first sample (in Dbin)
 * @param iech2 Rank of the second sample (in Dbout)
 * @param val1  Returned value at first sample
 * @param val2  Returned value at the second sample
 */
void NoStatArray::_getInfoFromDb(int ipar,
                                 int iech1,
                                 int iech2,
                                 double *val1,
                                 double *val2)
{
  *val1 = *val2 = TEST;
  *val2 = TEST;
  if (_dbin != (Db *) NULL && ! _attIn.empty())
    *val1 = getValue(ipar, 1, iech1);
  if (_dbout != (Db *) NULL && ! _attOut.empty())
    *val2 = getValue(ipar, 2, iech2);

  if (FFFF(*val1) && FFFF(*val2))
    my_throw("Non-stationary information may not be undefined");

  if (! FFFF(*val1))
    *val2 = *val1;
  if (! FFFF(*val2))
    *val1 = *val2;
}

/**
 * Update the Model according to the Non-stationary parameters
 * @param model Model to be patched
 * @param iech1 Rank of the target within Db1 (or -1)
 * @param iech2 Rank of the target within Dbout (or -2)
 */
void NoStatArray::updateModel(Model* model, int iech1, int iech2)
{
  double val1, val2;

  // Loop on the elements that can be updated one-by-one

  for (int ipar = 0; ipar < model->getNoStatElemNumber(); ipar++)
  {
    int icov = model->getNoStat().getICov(ipar);
    int type = model->getNoStat().getType(ipar);

    if (type == CONS_SILL)
    {
      _getInfoFromDb(ipar, iech1, iech2, &val1, &val2);
      int iv1  = model->getNoStat().getIV1(ipar);
      int iv2  = model->getNoStat().getIV2(ipar);
      model->setSill(icov, iv1, iv2, sqrt(val1 * val2));
    }
    else if (type == CONS_PARAM)
    {
      _getInfoFromDb(ipar, iech1, iech2, &val1, &val2);
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
          _getInfoFromDb(ipar, iech1, iech2, &angle1[idim], &angle2[idim]);
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
          _getInfoFromDb(ipar, iech1, iech2, &scale1[idim], &scale2[idim]);
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
          _getInfoFromDb(ipar, iech1, iech2, &range1[idim], &range2[idim]);
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
void NoStatArray::updateModel(Model* model, int vertex)
{

  // Loop on the elements that can be updated one-by-one

  for (int ipar = 0; ipar < model->getNoStatElemNumber(); ipar++)
  {
    int icov = model->getNoStat().getICov(ipar);
    ENUM_CONS type = model->getNoStat().getType(ipar);

    if (type == CONS_SILL)
    {
      double sill = getValue(ipar, 0, vertex);
      int iv1  = model->getNoStat().getIV1(ipar);
      int iv2  = model->getNoStat().getIV2(ipar);
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
    int iatt = _attIn[ipar];
    for (int iech = 0; iech < get_NECH(_dbin); iech++)
    {
      if (! get_ACTIVE(_dbin, iech)) continue;
      vec.push_back(get_ARRAY(_dbin,iech,iatt));
    }
  }
  else
  {
    if (_dbout == (Db *) NULL) return sstr.str();
    int iatt = _attOut[ipar];
    for (int iech = 0; iech < get_NECH(_dbout); iech++)
    {
      if (! get_ACTIVE(_dbout, iech)) continue;
      vec.push_back(get_ARRAY(_dbout,iech,iatt));
    }
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

