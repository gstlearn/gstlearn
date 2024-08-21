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
#include "Enum/ELoc.hpp"

#include "Basic/NamingConvention.hpp"
#include "Basic/OptDbg.hpp"
#include "Calculators/CalcMigrate.hpp"
#include "Calculators/ACalcDbToDb.hpp"
#include "Geometry/GeometryHelper.hpp"
#include "Db/DbGrid.hpp"
#include "Db/Db.hpp"
#include "Morpho/Morpho.hpp"
#include "Tree/Ball.hpp"

#include "geoslib_old_f.h"

#include <math.h>

#define RES(nval,idim)      (res[(idim) + (ndim+1) * (nval)])

// The next functions must be static as they ar called in GLOBAL functions
// defined in this class

/*****************************************************************************
 *!
 ** Calculate the shift of a grid node (rank i)
 **
 ** \param[in]  rank   Rank of the shift
 ** \param[in]  dbgrid Pointer to the Db
 ** \param[in]  ndim   Space dimension
 ** \param[in]  indg1  Input index array
 ** \param[in]  prop   Array of distance (proportions)
 **
 ** \param[out] indg2  Output index array
 ** \param[out] weight Output weight
 **
 *****************************************************************************/
static void st_shift(int rank,
                     Db *dbgrid,
                     const VectorInt &indg1,
                     const VectorDouble &prop,
                     VectorInt &indg2,
                     double *weight)
{
  int idim, ndiv, ival, ndim;
  double wgt;

  wgt = 1.;
  ndim = dbgrid->getNDim();
  ndiv = (int) pow(2., ndim - 1);
  for (idim = ndim - 1; idim >= 0; idim--)
  {
    ival = rank / ndiv;
    rank = rank - ndiv * ival;
    ndiv /= 2;
    indg2[idim] = indg1[idim] + ival;
    wgt *= (ival > 0) ? prop[idim] : (1. - prop[idim]);
  }
  *weight = wgt;
}

/*****************************************************************************
 *!
 ** Evaluate the value and the weight for multilinear interpolation
 **
 ** \return  Error return code (1 if the sample is outside the grid or if
 ** \return  the sample value is not defined)
 **
 ** \param[in]  db_grid   descriptor of the grid parameters
 ** \param[in]  indg      Working array
 ** \param[in]  iatt      Rank of the grid attribute
 **
 ** \param[out]  value    Output value
 **
 *****************************************************************************/
static int st_multilinear_evaluate(DbGrid *db_grid,
                                   const VectorInt &indg,
                                   int iatt,
                                   double *value)
{
  int jech = db_grid->indiceToRank(indg);
  if (jech < 0) return (1);
  if (!db_grid->isActive(jech)) return (1);
  *value = db_grid->getArray(jech, iatt);
  if (FFFF(*value)) return (1);
  return (0);
}

/*****************************************************************************/
/*!
 ** Perform the multi-linear interpolation from a regular grid Db
 **
 ** \return  Interpolated value
 **
 ** \param[in]  dbgrid    descriptor of the grid parameters
 ** \param[in]  iatt      rank of the target variable in dbgrid
 ** \param[in]  distType  Type of distance for calculating maximum distance
 **                       1 for L1 and 2 for L2 distance
 ** \param[in]  dmax      Array of maximum distances (optional)
 ** \param[in]  coor      Coordinates of the target point
 **
 *****************************************************************************/
static double st_multilinear_interpolation(DbGrid *dbgrid,
                                           int iatt,
                                           int distType,
                                           const VectorDouble& dmax,
                                           const VectorDouble& coor)
{
  int ndim = dbgrid->getNDim();
  int number = (int) pow(2., ndim);
  VectorInt iwork2(ndim);
  VectorInt indg(ndim);
  VectorDouble aux(ndim);
  VectorDouble prop(ndim);

  /* Identify the closest grid node */

  if (point_to_grid(dbgrid, coor.data(), 0, indg.data()) != 0) return TEST;
  dbgrid->indicesToCoordinateInPlace(indg, aux);

  /* Calculate distance to lower corner as proportion of the grid mesh */

  double rtot = 0.;
  bool dmaxEmpty = dmax.empty();
  for (int idim = 0; idim < ndim; idim++)
  {
    double mesh = dbgrid->getDX(idim);
    double delta = coor[idim] - aux[idim];
    if (delta < 0)
    {
      indg[idim]--;
      delta += mesh;
    }
    if (! dmaxEmpty && distType == 1)
    {
      if (dmax[idim] <= 0) return TEST;
      double ratio = delta / dmax[idim];
      rtot += ratio * ratio;
    }
    if (! dmaxEmpty && delta > dmax[idim]) return TEST;
    prop[idim] = delta / mesh;
  }
  if (! dmaxEmpty && distType == 1 && rtot > 1.) return TEST;

  /* Calculate the estimation */

  double estim = 0.;
  double wgt_tot = 0.;
  for (int i = 0; i < number; i++)
  {
    /* Get the shift */

    double weight;
    st_shift(i, dbgrid, indg, prop, iwork2, &weight);
    if (ABS(weight) < EPSILON6) continue;

    /* Get the sample value */

    double value;
    if (!st_multilinear_evaluate(dbgrid, iwork2, iatt, &value))
    {
      estim += weight * value;
      wgt_tot += weight;
    }
    else
    {
      estim = TEST;
      break;
    }
  }
  if (!FFFF(estim)) estim /= wgt_tot;

  return estim;
}

static double st_distance_modify(DbGrid* dbgrid,
                                 int ig,
                                 Db* dbpoint,
                                 int ip,
                                 VectorDouble& dvect,
                                 int flag_aniso,
                                 int iatt_time,
                                 int iatt_angle,
                                 int iatt_scaleu,
                                 int iatt_scalev,
                                 int iatt_scalew)
{
  double dd = TEST;

  // Case of anisotropic distances
  if (flag_aniso)
  {
    double cosa, sina, x, y, dloc;

    int ndim = dbgrid->getNDim();
    double angle = 0.;
    double scaleu = 1.;
    double scalev = 1.;
    double scalew = 1.;

    // Rotation angle in degrees (optional)
    if (iatt_angle >= 0 && ndim >= 2)
    {
      angle = dbgrid->getArray(ig, iatt_angle);
      GH::rotationGetSinCos(angle, &cosa, &sina);
      x = dvect[0];
      y = dvect[1];
      dvect[0] = x * cosa + y * sina;
      dvect[1] = y * cosa - x * sina;
    }

    // Scaled distance (optional)
    dd = 0;
    if (ndim >= 1)
    {
      if (iatt_scaleu >= 0) scaleu = dbgrid->getArray(ig, iatt_scaleu);
      if (!FFFF(scaleu) && scaleu > 0)
      {
        dloc = dvect[0] / scaleu;
        dd += dloc * dloc;
      }
    }
    if (ndim >= 2)
    {
      if (iatt_scalev >= 0) scalev = dbgrid->getArray(ig, iatt_scalev);
      if (!FFFF(scalev) && scalev > 0.)
      {
        dloc = dvect[1] / scalev;
        dd += dloc * dloc;
      }
    }
    if (ndim >= 3)
    {
      if (iatt_scalew >= 0) scalew = dbgrid->getArray(ig, iatt_scalew);
      if (!FFFF(scalew) && scalew > 0.)
      {
        dloc = dvect[2] / scalew;
        dd += dloc * dloc;
      }
    }

    // Complementary increments
    for (int idim = 4; idim < ndim; idim++)
    {
      dloc = dvect[idim - 1];
      dd += dloc * dloc;
    }
    dd = sqrt(dd);
  }

  // Add the Time Shift penalty (optional)
  if (iatt_time >= 0)
  {
    double time = dbpoint->getArray(ip, iatt_time);
    if (time > 0)
      dd = (FFFF(dd)) ? time : dd + time;
  }
  return dd;
}

/*****************************************************************************/
/*!
 **  Update minimum distance and rank of the corresponding sample
 **
 ** \param[in]  dbgrid      Descriptor of the grid parameters
 ** \param[in]  ig          Rank of the sample in the Grid file
 ** \param[in]  dbpoint     Descriptor of the point parameters
 ** \param[in]  ip          Rank of the sample in the Point file
 ** \param[in]  flag_aniso  1 if anisotropic distance must be calculated
 ** \param[in]  iatt_time   Optional variable for Time shift
 ** \param[in]  iatt_angle  Optional variable for anisotropy angle (around Z)
 ** \param[in]  iatt_scaleu Optional variable for anisotropy scale factor (U)
 ** \param[in]  iatt_scalev Optional variable for anisotropy scale factor (V)
 ** \param[in]  iatt_scalew Optional variable for anisotropy scale factor (W)
 **
 ** \param[in,out]  ipmin   Rank of the Point sample
 ** \param[in,out]  ddmin   Minimum distance
 ** \param[out]     dvect   Working vector
 **
 ** \remarks The Time Shift is an optional variable which increases the
 ** \remarks distance (the time-to-distance conversion is assumed to be 1)
 ** \remarks Only positive Time Shifts are considered
 **
 *****************************************************************************/
void st_get_closest_sample(DbGrid *dbgrid,
                           int ig,
                           Db *dbpoint,
                           int ip,
                           int flag_aniso,
                           int iatt_time,
                           int iatt_angle,
                           int iatt_scaleu,
                           int iatt_scalev,
                           int iatt_scalew,
                           int *ipmin,
                           double *ddmin,
                           VectorDouble& dvect)
{
  // Calculate the euclidean distance
  double dd = distance_inter(dbgrid, dbpoint, ig, ip, nullptr);

  // Distance modifier
  if (flag_aniso || iatt_time >= 0)
  {
    double dd2 = st_distance_modify(dbgrid, ig, dbpoint, ip, dvect, flag_aniso,
                                    iatt_time, iatt_angle,
                                    iatt_scaleu, iatt_scalev, iatt_scalew);
    if (! FFFF(dd2)) dd = dd2;
  }

  // Evaluate the closest distance
  if (dd < (*ddmin))
  {
    (*ddmin) = dd;
    (*ipmin) = ip;
  }
}

/**
 * Get the rank of the next sample just above the target
 * To speed up the process, this operation is performed
 * starting from the rank assigned to the previous sample
 * (this assumes that samples are treated by increasing coordinate)
 * @param ip0_init Rank of the previous sample
 * @param rank     Array of sample ordering
 * @param xtab     Array of sample coordinates
 * @param xtarget  Target coordinate
 * @return Rank of the sample just above the target (or equal)
 *
 * @note: The use of 'ip0_init' which could be different from 0 has been abandoned temporarily
 */
int st_next_sample(int ip0_init,
                   const VectorInt    &rank,
                   const VectorDouble &xtab,
                   double xtarget)
{
  int np = (int) rank.size();
  for (int ip = 0; ip < np; ip++)
  {
    int jp = ip0_init + ip;
    if (jp >= np) jp -= np;
    if (xtarget <= xtab[rank[jp]]) return jp;
  }
  return ip0_init;
}

/*****************************************************************************/
/*!
 **  Locate a set of points on a grid
 **
 ** \return  Number of samples located on the grid
 **
 ** \param[in]  db_point  descriptor of the point parameters
 ** \param[in]  db_grid   descriptor of the grid parameters
 **
 ** \param[out]  coor     Working array
 ** \param[out]  tab      Output array (Dimension: Number of point samples)
 **
 ** \remark  The array tab contains the index of the closest grid node
 ** \remark  even if the sample does not lie within the grid
 **
 *****************************************************************************/
static int st_locate_point_on_grid(const Db *db_point,
                                   const DbGrid *db_grid,
                                   VectorDouble &coor,
                                   VectorDouble &tab)
{
  int ndim = db_grid->getNDim();

  int number = 0;
  for (int iech = 0; iech < db_point->getSampleNumber(); iech++)
  {
    tab[iech] = TEST;
    if (!db_point->isActive(iech)) continue;
    for (int idim = 0; idim < ndim; idim++)
      coor[idim] = db_point->getCoordinate(iech, idim);
    int iad = db_grid->getGrid().coordinateToRank(coor);
    if (iad >= 0)
    {
      tab[iech] = iad;
      number++;
    }
  }
  return (number);
}

/*****************************************************************************/
/*!
 **  Locate a set of points on a grid
 **
 ** \return  Number of samples located on the grid
 **
 ** \param[in]  np        Number of samples
 ** \param[in]  coords    Array of coordinates (dimension: ndim, np)
 ** \param[in]  db_grid   descriptor of the grid parameters
 **
 ** \param[out]  tab      Output array (Dimension: Number of discretized points)
 **
 ** \remark  The array tab contains the index of the closest grid node
 ** \remark  even if the sample does not lie within the grid
 ** \remark  This function is limited to 3D space maximum. The consistency
 ** \remark  of space dimension must have been performed beforehand
 **
 *****************************************************************************/
static int st_locate_coor_on_grid(int np,
                                  const VectorVectorDouble& coords,
                                  const DbGrid *db_grid,
                                  VectorDouble& tab)
{
  int iech, iad, number;
  int ndim = (int) coords.size();
  VectorDouble local(ndim,0.);

  /* Loop on the point samples */

  for (iech = number = 0; iech < np; iech++)
  {
    tab[iech] = TEST;
    for (int idim = 0; idim < ndim; idim++)
      local[idim] = coords[idim][iech];
    iad = db_grid->coordinateToRank(local);
    if (iad >= 0)
    {
      tab[iech] = iad;
      number++;
    }
  }
  return (number);
}

/*****************************************************************************/
/*!
 **  Check if vector is out of range by comparing each component to the
 **  maximum value defined per direction
 **
 ** \return  Error return code
 **
 ** \param[in]  ndim      Space dimension
 ** \param[in]  dvect     Vector
 ** \param[in]  distType  Type of distance for calculating maximum distance
 **                       1 for L1 and 2 for L2 distance
 ** \param[in]  dmax      Array of maximum distances (optional)
 **
 ** \remarks If 'dmax' is not defined, the function always returns 0
 ** \remarks The function returns 1 as soon as one vector component is
 ** \remarks larger than the maximum value for this direction
 **
 *****************************************************************************/
static int st_larger_than_dmax(int ndim,
                               const VectorDouble &dvect,
                               int distType,
                               const VectorDouble &dmax)
{
  double ratio, rtot;

  if (dmax.empty()) return 0;

  /* Dispatch according to the type of distance */

  if (distType == 1)
  {

    // L1 distance

    for (int idim = 0; idim < ndim; idim++)
    {
      if (ABS(dvect[idim]) > dmax[idim])
        return 1;
    }
  }
  else
  {

    // L2 Distance

    rtot = 0.;
    for (int idim = 0; idim < ndim; idim++)
    {
      if (dmax[idim] <= 0) return (1);
      ratio = dvect[idim] / dmax[idim];
      rtot += ratio * ratio;
    }
    if (rtot > 1)
      return 1;
  }

  return 0;
}

/*****************************************************************************/
/*!
 **  Expand the joins at each cell in its vicinity, the radius is given
 **  per pixel in the array 'tab1'
 **
 ** \param[in]  flag_size   when 0, the norder pixels are painted with 1
 **                         when 1, the border pixels report the border thckness
 ** \param[in]  dbgrid      Descriptor of the grid parameters
 ** \param[in]  tab1        Array containing expansion radius
 **
 ** \param[in]  indg0       Array used for encoding/decoding
 ** \param[in]  indg        Array used for encoding/decoding
 ** \param[in]  tab2        Returned array
 **
 *****************************************************************************/
static void st_expand(int flag_size,
                      DbGrid *dbgrid,
                      VectorDouble &tab1,
                      VectorInt& indg0,
                      VectorInt& indg,
                      VectorDouble &tab2)
{
  int nech = dbgrid->getSampleNumber();
  int ndim = dbgrid->getNDim();

  /* Loop on the grid nodes */

  for (int iech = 0; iech < nech; iech++)
  {
    if (tab1[iech] <= 0.) continue;

    if (tab1[iech] <= 1.)
      tab2[iech] = 1.;
    else
    {
      int radius = (int)tab1[iech];
      dbgrid->rankToIndice(iech, indg0);

      for (int idim = 0; idim < ndim; idim++)
        for (int ifois = -1; ifois <= 1; ifois += 2)
          for (int irad = 0; irad < radius; irad++)
          {
            for (int jdim = 0; jdim < ndim; jdim++)
              indg[jdim] = indg0[jdim];
            indg[idim] = indg0[idim] + irad * ifois;
            int jech   = dbgrid->indiceToRank(indg);
            if (jech >= 0) tab2[jech] = (flag_size) ? radius : 1.;
          }
    }
  }
}

/*****************************************************************************/
/*!
 **  Migrates a variable from grid structure into a variable in point structure
 **
 ** \return  Error return code
 **
 ** \param[in]  db_grid   descriptor of the grid parameters
 ** \param[in]  db_point  descriptor of the point parameters
 ** \param[in]  iatt      rank of the grid attribute
 ** \param[in]  distType  Type of distance for calculating maximum distance
 **                       1 for L1 and 2 for L2 distance
 ** \param[in]  dmax      Array of maximum distances (optional)
 **
 ** \param[out]  tab      Output array (Dimension: number of samples in db_point)
 **
 *****************************************************************************/
int CalcMigrate::_migrateGridToPoint(DbGrid* db_grid,
                                     Db* db_point,
                                     int iatt,
                                     int distType,
                                     const VectorDouble& dmax,
                                     VectorDouble& tab)
{
  if (!db_grid->hasLargerDimension(db_point)) return 1;
  int ndim_min = MIN(db_grid->getNDim(), db_point->getNDim());
  int ndim_max = MAX(db_grid->getNDim(), db_point->getNDim());
  VectorDouble dvect(ndim_max);
  VectorDouble coor(ndim_max);

  /* Define the default values for 'coor'*/

  db_point->getCoordinatesPerSampleInPlace(0, coor);

  /* Locate the samples on the grid */

  (void) st_locate_point_on_grid(db_point, db_grid, coor, tab);

  /* Loop on the point samples */

  for (int iech = 0; iech < db_point->getSampleNumber(); iech++)
  {
    if (FFFF(tab[iech])) continue;
    int rank = (int) tab[iech];
    if (! dmax.empty())
    {
      (void) distance_inter(db_grid, db_point, rank, iech, dvect.data());
      if (st_larger_than_dmax(ndim_min, dvect, distType, dmax)) continue;
    }
    tab[iech] = db_grid->getArray(rank, iatt);
  }
  return 0;
}

CalcMigrate::CalcMigrate()
    : ACalcDbToDb(false),
      _iattOut(-1),
      _iuids(),
      _distType(1),
      _dmax(),
      _flagFill(false),
      _flagInter(false),
      _flagLocate(false),
      _locatorType(ELoc::Z)
{
}

CalcMigrate::~CalcMigrate()
{
}

bool CalcMigrate::_check()
{
  if (! ACalcDbToDb::_check()) return false;

  if (! hasDbin()) return false;
  if (! hasDbout()) return false;

  if (_iuids.empty())
  {
    messerr("At least one variable should be defined");
    return false;
  }
  if (_distType != 1 && _distType != 2)
  {
    messerr("Argument 'dist_type'(%d)  should be 1 (for L1 distance) or 2 (for L2 distance",_distType);
    return false;
  }
  return true;
}

bool CalcMigrate::_preprocess()
{
  int nvar = _getNVar();
  _iattOut = _addVariableDb(2, 1, ELoc::UNKNOWN, 0, nvar, 0.);
  return (_iattOut >= 0);
}

bool CalcMigrate::_postprocess()
{
  /* Free the temporary variables */
  _cleanVariableDb(2);

  int nvar = _getNVar();
  _renameVariable(2, getDbin()->getNamesByUID(_iuids), ELoc::UNKNOWN, nvar,
                    _iattOut, String(), 1);

  if (_flagLocate)
    getDbout()->setLocatorsByUID(nvar, _iattOut, _locatorType);

  return true;
}

void CalcMigrate::_rollback()
{
  _cleanVariableDb(1);
}

int CalcMigrate::_getNVar() const
{
  return (int) _iuids.size();
}

/****************************************************************************/
/*!
 **  Standard Kriging
 **
 ** \return  Error return code
 **
 *****************************************************************************/
bool CalcMigrate::_run()
{
  int nvar = _getNVar();

  // Perform the migrations

  for (int i = 0; i < nvar; i++)
  {
    int iatt1 = _iuids[i];
    int iatt2 = _iattOut + i;
    if (_migrate(getDbin(), getDbout(), iatt1, iatt2, _distType, _dmax, _flagFill,
                 _flagInter, _flagBall)) return false;
  }

  return true;
}

/*****************************************************************************/
/*!
 **  Migrates a variable from one Db to another one
 **
 ** \return  Error return code
 **
 ** \param[in]  dbin       Descriptor of the input Db
 ** \param[in]  dbout      Descriptor of the output Db
 ** \param[in]  name       Name of the attribute to be migrated
 ** \param[in]  dist_type  Type of distance for calculating maximum distance
 **                        1 for L1 and 2 for L2 distance
 ** \param[in]  dmax       Array of maximum distances (optional)
 ** \param[in]  flag_fill  Filling option
 ** \param[in]  flag_inter Interpolation
 ** \param[in]  flag_ball  Use BallTree sorting algorithm when available
 ** \param[in]  namconv    Naming convention
 **
 *****************************************************************************/
int migrate(Db *dbin,
            Db *dbout,
            const String &name,
            int dist_type,
            const VectorDouble& dmax,
            bool flag_fill,
            bool flag_inter,
            bool flag_ball,
            const NamingConvention &namconv)
{
  CalcMigrate migrate;
  migrate.setDbin(dbin);
  migrate.setDbout(dbout);
  migrate.setNamingConvention(namconv);

  VectorInt iuids(1);
  iuids[0] = dbin->getUID(name);
  migrate.setIuids(iuids);
  migrate.setDistType(dist_type);
  migrate.setDmax(dmax);
  migrate.setFlagFill(flag_fill);
  migrate.setFlagInter(flag_inter);
  migrate.setFlagBall(flag_ball);

  // Run the calculator
  int error = (migrate.run()) ? 0 : 1;

  return error;
}

/*****************************************************************************/
/*!
 **  Migrates a set of variables from one Db to another one
 **
 ** \return  Error return code
 **
 ** \param[in]  dbin       Descriptor of the input Db
 ** \param[in]  dbout      Descriptor of the output Db
 ** \param[in]  names      Name of the attribute to be migrated
 ** \param[in]  dist_type  Type of distance for calculating maximum distance
 **                        1 for L1 and 2 for L2 distance
 ** \param[in]  dmax       Array of maximum distances (optional)
 ** \param[in]  flag_fill  Filling option
 ** \param[in]  flag_inter Interpolation
 ** \param[in]  flag_ball  Use BallTree sorting algorithm when available
 ** \param[in]  namconv    Naming convention
 **
 *****************************************************************************/
int migrateMulti(Db *dbin,
                 Db *dbout,
                 const VectorString &names,
                 int dist_type,
                 const VectorDouble& dmax,
                 bool flag_fill,
                 bool flag_inter,
                 bool flag_ball,
                 const NamingConvention &namconv)
{
  CalcMigrate migrate;
  migrate.setDbin(dbin);
  migrate.setDbout(dbout);
  migrate.setNamingConvention(namconv);

  VectorInt iuids = dbin->getUIDs(names);
  migrate.setIuids(iuids);
  migrate.setDistType(dist_type);
  migrate.setDmax(dmax);
  migrate.setFlagFill(flag_fill);
  migrate.setFlagInter(flag_inter);
  migrate.setFlagBall(flag_ball);

  // Run the calculator
  int error = (migrate.run()) ? 0 : 1;
  return error;
}

/*****************************************************************************/
/*!
 **  Migrates a variable from one Db to another one
 **
 ** \return  Error return code
 **
 ** \param[in]  dbin       Descriptor of the input Db
 ** \param[in]  dbout      Descriptor of the output Db
 ** \param[in]  atts       Array of attributes to be migrated
 ** \param[in]  dist_type  Type of distance for calculating maximum distance
 **                        1 for L1 and 2 for L2 distance
 ** \param[in]  dmax       Array of maximum distances (optional)
 ** \param[in]  flag_fill  Filling option
 ** \param[in]  flag_inter Interpolation
 ** \param[in]  flag_ball  Use BallTree sorting algorithm when available
 ** \param[in]  namconv    Naming Convention
 **
 *****************************************************************************/
int migrateByAttribute(Db *dbin,
                       Db *dbout,
                       const VectorInt& atts,
                       int dist_type,
                       const VectorDouble& dmax,
                       bool flag_fill,
                       bool flag_inter,
                       bool flag_ball,
                       const NamingConvention &namconv)
{
  CalcMigrate migrate;
  migrate.setDbin(dbin);
  migrate.setDbout(dbout);
  migrate.setNamingConvention(namconv);

  VectorInt iuids = atts;
  if (iuids.empty()) iuids = dbin->getAllUIDs();

  migrate.setIuids(iuids);
  migrate.setDistType(dist_type);
  migrate.setDmax(dmax);
  migrate.setFlagFill(flag_fill);
  migrate.setFlagInter(flag_inter);
  migrate.setFlagBall(flag_ball);

  // Run the calculator
  int error = (migrate.run()) ? 0 : 1;
  return error;
}

/*****************************************************************************/
/*!
 **  Migrates all z-locator variables from one Db to another one
 **
 ** \return  Error return code
 **
 ** \param[in]  dbin        Descriptor of the input Db
 ** \param[in]  dbout       Descriptor of the output Db
 ** \param[in]  locatorType Locator Type
 ** \param[in]  dist_type   Type of distance for calculating maximum distance
 **                         1 for L1 and 2 for L2 distance
 ** \param[in]  dmax        Array of maximum distances (optional)
 ** \param[in]  flag_fill   Filling option
 ** \param[in]  flag_inter  Interpolation
 ** \param[in]  flag_ball  Use BallTree sorting algorithm when available
 ** \param[in]  namconv     Naming convention
 **
 ** \remark The output variable receive the same locator as the input variables
 **
 *****************************************************************************/
int migrateByLocator(Db *dbin,
                     Db *dbout,
                     const ELoc& locatorType,
                     int dist_type,
                     const VectorDouble& dmax,
                     bool flag_fill,
                     bool flag_inter,
                     bool flag_ball,
                     const NamingConvention &namconv)
{
  CalcMigrate migrate;
  migrate.setDbin(dbin);
  migrate.setDbout(dbout);
  migrate.setNamingConvention(namconv);

  VectorString names = dbin->getNamesByLocator(locatorType);
  VectorInt iuids = dbin->getUIDs(names);
  migrate.setIuids(iuids);
  migrate.setDistType(dist_type);
  migrate.setDmax(dmax);
  migrate.setFlagFill(flag_fill);
  migrate.setFlagInter(flag_inter);
  migrate.setFlagBall(flag_ball);
  migrate.setFlagLocate(true);
  migrate.setLocatorType(locatorType);

  // Run the calculator
  int error = (migrate.run()) ? 0 : 1;
  return error;
}

/*****************************************************************************/
/*!
 **  Derive the external information(s) from the Output db (if Grid)
 **  to the Input Db
 **
 ** \return  Error return code
 **
 ** \param[in]  mode        1 for allocation; -1 for deallocation
 ** \param[in]  locatorType Type of the pointer (ELoc)
 ** \param[in]  dbin        Descriptor of the input Db
 ** \param[in]  dbout       Descriptor of the output Db
 **
 ** \param[out] flag_created True if variables have been created
 **
 ** \remark This function only functions when the Output Db is a grid
 ** \remark However, in case of a Point output Db, this function should not
 ** \remark be used: the external drift functions should already be present
 ** \remark in the output Db.
 ** \remark If this is not the case, an error is issued.
 **
 *****************************************************************************/
int manageExternalInformation(int mode,
                              const ELoc &locatorType,
                              Db *dbin,
                              Db *dbout,
                              bool *flag_created)
{
  //  VectorDouble tab;

  if (dbin == nullptr) return 0;
  int ninfo = get_LOCATOR_NITEM(dbout, locatorType);
  if (ninfo <= 0) return 0;

  /* Case when the Output Db is not a grid */

  if (!dbout->isGrid())
  {
    if (get_LOCATOR_NITEM(dbin, locatorType) == ninfo) return 0;
    messerr("The Output Db is not a Grid file");
    messerr("The Input Db does not contain the %d External Drifts");
    return 1;
  }
  DbGrid* dbgrid = dynamic_cast<DbGrid*>(dbout);

  /* Dispatch */

  if (mode > 0)
  {
    /* If the drift vector is present in the input file, skip the rest */

    if (dbin->getLocatorNumber(locatorType) >= ninfo)
    {
      *flag_created = false;
      return 0;
    }

    /* Creating variables */

    *flag_created = true;
    for (int info = 0; info < ninfo; info++)
    {
      String name = dbgrid->getNameByLocator(locatorType, info);
      if (migrate(dbgrid, dbin, name, 0, VectorDouble(), false, false, false))
        continue;
    }
    return 0;
  }
  if (*flag_created)
  {
    // If no variable has been created (when 'mode' == 0), then do nothing

    return 0;
  }
  for (int info = 0; info < ninfo; info++)
  {
    int jatt = dbin->getUIDByLocator(locatorType, info);
    dbin->deleteColumnByUID(jatt);
  }
  return 0;
}

/*****************************************************************************/
/*!
 ** Interpolate a variable from a grid Db on discretization points
 **
 ** \param[in]  db_grid   Descriptor of the grid parameters
 ** \param[in]  iatt      Rank of the attribute in db_grid
 ** \param[in]  np        Number of discretized points
 ** \param[in]  xp        Array of first coordinates
 ** \param[in]  yp        Array of second coordinates
 ** \param[in]  zp        Array of third coordinates
 **
 ** \param[out]  tab      Output array
 **
 ** \remark The arguments 'xp', 'yp' and 'zp' must be defined in accordance
 ** \remark with the space dimension in argument 'db_grid'
 **
 ** \remark A point which does not lie between two valuated grid nodes
 ** \remark (in all space dimensions) is always set to FFFF
 **
 *****************************************************************************/
int interpolateVariableToPoint(DbGrid *db_grid,
                               int iatt,
                               int np,
                               const double *xp,
                               const double *yp,
                               const double *zp,
                               double *tab)
{
  int error, ndim;
  VectorDouble coor(3);

  /* Initializations */

  error = 1;
  for (int idim = 0; idim < 3; idim++)
    coor[idim] = 0.;
  ndim = db_grid->getNDim();
  if (ndim > 3)
  {
    messerr("This procedure is limited to 3-D grid");
    goto label_end;
  }
  if ((ndim >= 1 && xp == nullptr) || (ndim >= 2 && yp == nullptr)
      || (ndim >= 3 && zp == nullptr))
  {
    messerr("The Grid space dimension (%d) must be in accordance with", ndim);
    messerr("the definition of arguments 'xp', 'yp' and 'zp'");
    goto label_end;
  }

  /* Loop on the point samples */

  for (int ip = 0; ip < np; ip++)
  {
    if (ndim >= 1) coor[0] = xp[ip];
    if (ndim >= 2) coor[1] = yp[ip];
    if (ndim >= 3) coor[2] = zp[ip];
    tab[ip] = st_multilinear_interpolation(db_grid, iatt, 0, VectorDouble(), coor);
  }

  /* Set the error return code */

  error = 0;

  /* Core deallocation */

  label_end: return (error);
}

/*****************************************************************************/
/*!
 **  Sampling vertices within a Grid between two points
 **
 ** \return  Array of sampled vertices
 **
 ** \param[in]  dbgrid      reference Grid
 ** \param[in]  x1          Array giving the coordinates of the first point
 ** \param[in]  x2          Array giving the coordinates of the second point
 ** \param[in]  ndisc       Number of discretized points in the segment
 ** \param[in]  ncut        Number of cutoffs
 ** \param[in]  cuts        Array of cutoffs
 **
 ** \param[out] nval_ret    Number of samples in the output array
 **
 ** \remarks This function considers the segment [x1,x2] and subdivises it
 ** \remarks into 'ndisc' intervals. The endpoints of each interval correspond
 ** \remarks to two points in the space
 ** \remarks At each endpoint, the target variable is interpolated from the grid
 ** \remarks If the target variable values cross a cutoff, the coordinates of
 ** \remarks the intersection are calculated.
 ** \remarks The program returns the list of all these intersection coordinates
 **
 ** TODO FUTURE_REFACTOR
 *****************************************************************************/
double* dbgridLineSampling(DbGrid *dbgrid,
                           const double *x1,
                           const double *x2,
                           int ndisc,
                           int ncut,
                           const double *cuts,
                           int *nval_ret)
{
  double *res, delta, vi1, vi2, cut, v1, v2;
  int ndim, iatt, nval;
  VectorDouble xi1;
  VectorDouble xi2;

  /* Initializations */

  *nval_ret = 0;
  res = nullptr;
  ndim = dbgrid->getNDim();
  iatt = dbgrid->getColIdxByLocator(ELoc::Z, 0);

  /* Preliminary checks */

  if (ndisc <= 1)
  {
    messerr("The number of discretization points must be larger than 1");
    goto label_end;
  }
  if (iatt < 0)
  {
    messerr("You need a target variable on the grid");
    goto label_end;
  }

  /* Core allocation */

  xi1.resize(ndim, 0);
  xi2.resize(ndim, 0);

  /* Loop on the discretized points */

  nval = 0;
  for (int idisc = 0; idisc < ndisc; idisc++)
  {

    /* Calculate the discretization segment */

    for (int idim = 0; idim < ndim; idim++)
    {
      delta = (x2[idim] - x1[idim]) / (double) ndisc;
      xi1[idim] = x1[idim] + delta * (double) (idisc);
      xi2[idim] = x1[idim] + delta * (double) (idisc + 1);
    }

    /* Calculate the target variable value at segment endpoints */

    vi1 = st_multilinear_interpolation(dbgrid, iatt, 0, VectorDouble(), xi1);
    vi2 = st_multilinear_interpolation(dbgrid, iatt, 0, VectorDouble(), xi2);
    v1 = MIN(vi1, vi2);
    v2 = MAX(vi1, vi2);

    /* Loop on the cutoffs */

    for (int icut = 0; icut < ncut; icut++)
    {
      cut = cuts[icut];
      if (cut < v1 || cut > v2) continue;
      res = (double*) mem_realloc((char* ) res,
                                  (ndim + 1) * (nval + 1) * sizeof(double), 0);
      if (res == nullptr) goto label_end;

      for (int idim = 0; idim < ndim; idim++)
      {
        delta = (v2 > v1) ? (xi2[idim] - xi1[idim]) / (v2 - v1) : 0.;
        RES(nval,idim) = xi1[idim] + delta * (cut - v1);
      }
      RES(nval,ndim) = icut + 1;
      nval++;
    }
  }
  *nval_ret = nval;

  label_end:
  return (res);
}

/*****************************************************************************/
/*!
 **  Expands a variable from the point structure
 **  into a variable in the grid structure
 **
 ** \return  Error return code
 **
 ** \param[in]  db_point    Descriptor of the point parameters
 ** \param[in]  db_grid     Descriptor of the grid parameters
 ** \param[in]  iatt        Rank of the point attribute
 ** \param[in]  iatt_time   Optional variable for Time shift
 ** \param[in]  iatt_angle  Optional variable for anisotropy angle (around Z)
 ** \param[in]  iatt_scaleu Optional variable for anisotropy scale factor (U)
 ** \param[in]  iatt_scalev Optional variable for anisotropy scale factor (V)
 ** \param[in]  iatt_scalew Optional variable for anisotropy scale factor (W)
 ** \param[in]  flag_index  1 if the Index must be assigned to grid node
 **                         0 the 'iatt' attribute is assigned instead
 ** \param[in]  distType    Type of distance for calculating maximum distance
 **                         1 for L1 and 2 for L2 distance
 ** \param[in]  dmax        Array of maximum distances (optional)
 **
 ** \param[out]  tab        Output array
 **
 ** \remarks When a Time Shift is present, this corresponds to Johnson-Mehl
 ** \remarks The Time Shift is an optional variable which increases the
 ** \remarks distance (the time-to-distance conversion is assumed to be 1)
 ** \remarks Only positive Time Shifts are considered
 **
 *****************************************************************************/
int expandPointToGrid(Db *db_point,
                      DbGrid *db_grid,
                      int iatt,
                      int iatt_time,
                      int iatt_angle,
                      int iatt_scaleu,
                      int iatt_scalev,
                      int iatt_scalew,
                      int flag_index,
                      int distType,
                      const VectorDouble &dmax,
                      VectorDouble &tab)
{
  double dd1d, ddmin;
  int jpmin;

  /* Preliminary checks */

  if (!db_point->hasLargerDimension(db_grid)) return 1;
  int ndim_min = MIN(db_point->getNDim(), db_grid->getNDim());
  int ndim_max = MAX(db_point->getNDim(), db_grid->getNDim());
  bool flag_aniso = 0;
  if (ndim_min >= 1 && iatt_scaleu >= 0) flag_aniso = 1;
  if (ndim_min >= 2 && iatt_scalev >= 0) flag_aniso = 1;
  if (ndim_min >= 3 && iatt_scalew >= 0) flag_aniso = 1;
  int idim_ref = ndim_min - 1;
  double dmax_ref = 1.e30;
  if (! dmax.empty()) dmax_ref = dmax[idim_ref];

  // Core allocation

  int ng = db_grid->getSampleNumber();
  int np = db_point->getSampleNumber(true);
  VectorDouble dvect(ndim_max);
  VectorDouble xtab(np);

  /* Sort the point samples according to their coordinate ranked 'idim_ref' */

  for (int ip = np = 0; ip < db_point->getSampleNumber(); ip++)
  {
    if (!db_point->isActive(ip)) continue;
    if (FFFF(db_point->getArray(ip, iatt))) continue;
    xtab[np] = db_point->getCoordinate(ip, idim_ref);
    np++;
  }
  VectorInt rank = VH::orderRanks(xtab, true);

  /* Calculate the maximum time (if defined) */

  double time_max = 0.;
  if (iatt_time >= 0)
  {
    for (int ip = 0; ip < np; ip++)
    {
      double time = db_point->getArray(ip, iatt_time);
      if (time > time_max) time_max = time;
    }
  }

  /* Loop on the grid nodes */

  int ip0 = 0;
  for (int ig = 0; ig < ng; ig++)
  {
    if (!db_grid->isActive(ig)) continue;
    double xtarget = db_grid->getCoordinate(ig, idim_ref);

    /* Locate the grid node within the ordered list (1D coordinate) */

    ip0 = st_next_sample(0, rank, xtab, xtarget);

    /* Calculate minimum distance between the two closest ordered samples */

    jpmin = -1;
    ddmin = dmax_ref;

    /* Look for closer points for samples located below rank[ip0] */

    for (int ip = ip0; ip >= 0; ip--)
    {
      int jp = rank[ip];
      dd1d = ABS(xtab[jp] - xtarget);
      if (dd1d >= ddmin)
        break;
      st_get_closest_sample(db_grid, ig, db_point, jp, flag_aniso, iatt_time,
                            iatt_angle, iatt_scaleu, iatt_scalev, iatt_scalew,
                            &jpmin, &ddmin, dvect);
    }

    /* Look for closer points for samples located above rank[ip0+1] */

    for (int ip = ip0 + 1; ip < np; ip++)
    {
      int jp = rank[ip];
      dd1d = ABS(xtab[jp] - xtarget);
      if (dd1d >= ddmin)
        break;
      st_get_closest_sample(db_grid, ig, db_point, jp, flag_aniso, iatt_time,
                            iatt_angle, iatt_scaleu, iatt_scalev, iatt_scalew,
                            &jpmin, &ddmin, dvect);
    }

    /* Truncation by 'dmax' if provided */

    if (! dmax.empty() && jpmin >= 0)
    {
      (void) distance_inter(db_grid, db_point, ig, jpmin, dvect.data());
      if (st_larger_than_dmax(ndim_min, dvect, distType, dmax)) continue;
    }

    /* Set the value */

    if (jpmin < 0)
      tab[ig] = TEST;
    else
    {
      if (flag_index)
        tab[ig] = (double) jpmin;
      else
        tab[ig] = db_point->getArray(jpmin, iatt);
    }
  }
  return 0;
}

/*****************************************************************************/
/*!
 **  Expands a variable from one point Db
 **  into a variable at points defined by coordinate vectors (maximum 3D)
 **
 ** \return  Error return code
 **
 ** \param[in]  db1   descriptor of the input parameters
 ** \param[in]  iatt  rank of the input attribute
 ** \param[in]  coords    Array of coordinates
 **
 ** \param[out]  tab      Output array (Dimension: number of discretized points)
 **
 *****************************************************************************/
int expandPointToCoor(const Db *db1,
                      int iatt,
                      const VectorVectorDouble &coords,
                      VectorDouble &tab)
{
  double *tab1 = nullptr;
  double *tab2 = nullptr;

  /* Preliminary checks */

  int ndim = db1->getNDim();
  int ndimp = (int) coords.size();
  int np = (int) coords[0].size();
  if (ndim != ndimp)
  {
    messerr("The Space Dimension of the First Db (%d)", ndim);
    messerr("must be equal to the Space Dimension of the coordinate arrays",
            ndimp);
    return 1;
  }

  ut_distance_allocated(ndim, &tab1, &tab2);

  /* Loop on the output structure */

  for (int iech2 = 0; iech2 < np; iech2++)
  {
    /* Store the coordinates */

    if (ndim >= 1) tab2[0] = coords[0][iech2];
    if (ndim >= 2) tab2[1] = coords[1][iech2];
    if (ndim >= 3) tab2[2] = coords[2][iech2];

    /* Loop on the input structure */

    double distmin = 1.e30;
    int iechmin = -1;
    for (int iech1 = 0; iech1 < db1->getSampleNumber(); iech1++)
    {
      if (!db1->isActive(iech1)) continue;
      for (int idim = 0; idim < ndim; idim++)
        tab1[idim] = db1->getCoordinate(iech1, idim);

      double dist = ut_distance(ndim, tab1, tab2);
      if (dist < distmin)
      {
        distmin = dist;
        iechmin = iech1;
      }
    }
    if (iechmin >= 0) tab[iech2] = db1->getArray(iechmin, iatt);
  }
  return 0;
}

/*****************************************************************************/
/*!
 **  Plunge a set of isolated points within a discretization grid
 **  in order to compute the voronoi of the points and derive:
 **  - the statistics on the volume and perimeter of the cells
 **  - the edge between cells
 **
 ** \return  Error return code
 **
 ** \param[in]  dbpoint     Descriptor of the point parameters
 ** \param[in]  dbgrid      Descriptor of the grid parameters
 ** \param[in]  option      Connectivity option (0 for cross and 1 for block)
 ** \param[in]  flag_size   When 1, the border pixels report the border thickness
 **                         When 0, the border pixels are painted in 1
 ** \param[in]  iatt_time   Attribute of 'dbpoint'for Time shift (optional)
 ** \param[in]  iatt_size   Attribute of 'dbpoint' giving size (optional)
 ** \param[in]  iatt_angle  Optional variable for anisotropy angle (around Z)
 ** \param[in]  iatt_scaleu Optional variable for anisotropy scale factor (U)
 ** \param[in]  iatt_scalev Optional variable for anisotropy scale factor (V)
 ** \param[in]  iatt_scalew Optional variable for anisotropy scale factor (W)
 **
 ** \remarks The value of 'flag_index' can be turned on for assigning
 ** \remarks the sample index to the grid cell (instead of the 'iatt' value)
 ** \remarks using: set_keypair("PTB_flag_index")
 **
 *****************************************************************************/
int pointToBlock(Db *dbpoint,
                 DbGrid *dbgrid,
                 int option,
                 int flag_size,
                 int iatt_time,
                 int iatt_size,
                 int iatt_angle,
                 int iatt_scaleu,
                 int iatt_scalev,
                 int iatt_scalew)
{
  int iatt_edge, iatt_rank, iatt_surf, iatt_vol, iatt_code;
  int val_iech, val_jech, jech, ndim, nvois, lec, error, flag_index;
  VectorInt indret;
  VectorDouble tab1, tab2;

  /* Initializations */

  error = 1;
  iatt_rank = -1;
  if (!dbgrid->hasSameDimension(dbpoint)) return 1;
  ndim = dbgrid->getNDim();
  flag_index = (int) get_keypone("PTB_Flag_Index", 0.);

  /* Core allocation */

  VectorInt indg0(ndim);
  VectorInt indg(ndim);;
  tab1.resize(dbgrid->getSampleNumber());
  tab2.resize(dbgrid->getSampleNumber(), -1.);

  /* Variable allocation */

  iatt_edge = dbgrid->addColumnsByConstant(1, 0.);
  if (iatt_edge < 0) goto label_end;
  iatt_rank = dbpoint->addColumnsByConstant(1, 0.);
  if (iatt_rank < 0) goto label_end;
  iatt_surf = dbpoint->addColumnsByConstant(1, 0.);
  if (iatt_surf < 0) goto label_end;
  iatt_vol = dbpoint->addColumnsByConstant(1, 0.);
  if (iatt_vol < 0) goto label_end;
  iatt_code = dbpoint->addColumnsByConstant(1, 1.);
  if (iatt_code < 0) goto label_end;

  /* Create the sample rank attribute and expand it over the grid */

  for (int iech = 0; iech < dbpoint->getSampleNumber(); iech++)
    dbpoint->setArray(iech, iatt_rank, (double) iech);
  if (expandPointToGrid(dbpoint, dbgrid, iatt_rank, iatt_time, iatt_angle,
                        iatt_scaleu, iatt_scalev, iatt_scalew, flag_index, 0,
                        VectorDouble(), tab1)) goto label_end;

  /* When sample index needed, work is over, except copying 'tab1' into 'tab2'*/

  if (flag_index)
  {
    for (int i = 0; i < dbgrid->getSampleNumber(); i++)
      tab2[i] = tab1[i];
    goto label_suite;
  }

  /* Define the neighboring elements */

  indret = gridcell_neigh(ndim, option, 1, 1, 0);
  nvois = (int) indret.size() / ndim;

  /* Loop on the grid nodes */

  for (int iech = 0; iech < dbgrid->getSampleNumber(); iech++)
  {
    if (!dbgrid->isActive(iech)) continue;

    /* Identify the sample within the grid */

    val_iech = (int)tab1[iech];
    dbgrid->rankToIndice(iech, indg0);

    /* Increment the volume by one */

    dbpoint->updArray(val_iech, iatt_vol, EOperator::ADD, 1.);

    /* Loop on the neighborhood nodes */

    lec = 0;
    for (int ivois = 0; ivois < nvois; ivois++)
    {
      for (int idim = 0; idim < ndim; idim++)
        indg[idim] = indg0[idim] + indret[lec++];
      jech = dbgrid->indiceToRank(indg);
      if (jech < 0)
      {

        /* The center point is located on the edge of the field */

        dbpoint->setArray(val_iech, iatt_code, 2.);
      }
      else
      {
        if (!dbgrid->isActive(jech)) continue;
        val_jech = (int) tab1[jech];
        if (val_iech == val_jech) continue;

        /* Increment the perimeter by one */

        dbpoint->updArray(val_iech, iatt_surf, EOperator::ADD, 1.);

        /* Set the edge */

        tab2[iech] = (iatt_size >= 0) ? dbpoint->getArray(val_iech, iatt_size) : 1.;
      }
    }
  }

  /* If size is defined, expand the perimeter */

  if (iatt_size >= 0)
  {
    for (int i = 0; i < dbgrid->getSampleNumber(); i++)
      tab1[i] = tab2[i];
    st_expand(flag_size, dbgrid, tab1, indg0, indg, tab2);
  }

  /* Transform values into 0 and 1 */

  for (int i = 0; i < dbgrid->getSampleNumber(); i++)
  {
    if (flag_size)
      tab2[i] = (tab2[i] < 0) ? 0 : tab2[i];
    else
      tab2[i] = (tab2[i] < 0) ? 0 : 1;
  }

  /* Save the array 'tab' in the Grid Db file */

  label_suite: dbgrid->setColumnByUID(tab2, iatt_edge);

  /* Set the error return code */

  error = 0;

  /* Core deallocation */

  label_end:
    if (iatt_rank >= 0) dbpoint->deleteColumnByUID(iatt_rank);
  return (error);
}

/*****************************************************************************/
/*!
 **  Migrates a variable from the grid structure
 **  into a variable at points defined by coordinate vectors
 **
 ** \return  Error return code
 **
 ** \param[in]  db_grid   descriptor of the grid parameters
 ** \param[in]  iatt      rank of the grid attribute
 ** \param[in]  coords    Array of coordinates (dimension: ndim, np)
 **
 ** \param[out]  tab      Output array (Dimension: number of discretized points)
 **
 *****************************************************************************/
int migrateGridToCoor(const DbGrid *db_grid,
                      int iatt,
                      const VectorVectorDouble &coords,
                      VectorDouble &tab)
{
  int ndim = (int) coords.size();
  int np = (int) coords[0].size();
  if (db_grid->getNDim() != ndim)
  {
    messerr("The Space Dimension of the First Db (%d)", db_grid->getNDim());
    messerr("must be equal to the Space Dimension of the coordinate arrays",
            ndim);
    return 1;
  }

  /* Locate the samples on the grid */

  (void) st_locate_coor_on_grid(np, coords, db_grid, tab);

  /* Loop on the point samples */

  for (int iech = 0; iech < np; iech++)
  {
    if (FFFF(tab[iech])) continue;
    tab[iech] = db_grid->getArray((int) tab[iech], iatt);
  }
  return 0;
}

/*****************************************************************************/
/*!
 **  Expands a variable from one point Db
 **  into a variable at points defined by coordinate vectors (maximum 3D)
 **
 ** \return  Error return code
 **
 ** \param[in]  db1   descriptor of the input parameters
 ** \param[in]  iatt  rank of the input attribute
 ** \param[in]  coords    Array of coordinates
 **
 ** \param[out]  tab      Output array (Dimension: number of discretized points)
 **
 *****************************************************************************/
int expand_point_to_coor(const Db *db1,
                         int iatt,
                         const VectorVectorDouble& coords,
                         VectorDouble& tab)
{
  double *tab1 = nullptr;
  double *tab2 = nullptr;

  /* Preliminary checks */

  int ndim = db1->getNDim();
  int ndimp = (int) coords.size();
  int np = (int) coords[0].size();
  if (ndim != ndimp)
  {
    messerr("The Space Dimension of the First Db (%d)", ndim);
    messerr("must be equal to the Space Dimension of the coordinate arrays",
            ndimp);
    return 1;
  }

  ut_distance_allocated(ndim, &tab1, &tab2);

  /* Loop on the output structure */

  for (int iech2 = 0; iech2 < np; iech2++)
  {
    /* Store the coordinates */

    if (ndim >= 1) tab2[0] = coords[0][iech2];
    if (ndim >= 2) tab2[1] = coords[1][iech2];
    if (ndim >= 3) tab2[2] = coords[2][iech2];

    /* Loop on the input structure */

    double distmin = 1.e30;
    int iechmin = -1;
    for (int iech1 = 0; iech1 < db1->getSampleNumber(); iech1++)
    {
      if (!db1->isActive(iech1)) continue;
      for (int idim = 0; idim < ndim; idim++)
        tab1[idim] = db1->getCoordinate(iech1, idim);

      double dist = ut_distance(ndim, tab1, tab2);
      if (dist < distmin)
      {
        distmin = dist;
        iechmin = iech1;
      }
    }
    if (iechmin >= 0) tab[iech2] = db1->getArray(iechmin, iatt);
  }
  return 0;
}

/*****************************************************************************/
/*!
 **  Migrates a variable from one Db to another one
 **
 ** \return  Error return code
 **
 ** \param[in]  db1        descriptor of the input Db
 ** \param[in]  db2        descriptor of the output Db
 ** \param[in]  iatt1      Attribute in Db1 to be migrated
 ** \param[in]  iatt2      Attribute in Db2 where the result must be stored
 ** \param[in]  distType   Type of distance for calculating maximum distance
 **                        1 for L1 and 2 for L2 distance
 ** \param[in]  dmax       Array of maximum distances (optional)
 ** \param[in]  flag_fill  Filling option
 ** \param[in]  flag_inter Interpolation
 ** \param[in]  flag_ball  Use the BallTree sort for speeding up calculations
 **
 *****************************************************************************/
int CalcMigrate::_migrate(Db *db1,
                          Db *db2,
                          int iatt1,
                          int iatt2,
                          int distType,
                          const VectorDouble &dmax,
                          bool flag_fill,
                          bool flag_inter,
                          bool flag_ball)
{
  int size = db2->getSampleNumber();
  VectorDouble tab(size, TEST);

  if (db2->isGrid())

  {
    DbGrid* db2grid = dynamic_cast<DbGrid*>(db2);

    // To Grid
    if (db1->isGrid())
    {
      DbGrid* db1grid = dynamic_cast<DbGrid*>(db1);

      // Grid to Grid
      if (flag_fill)
      {
        // Grid to Grid (flag_fill = TRUE)
        if (_expandGridToGrid(db1grid, db2grid, iatt1, distType, dmax, tab))
          return 1;
      }
      else
      {
        // Grid to Grid (flag_fill = FALSE)
        if (_migrateGridToGrid(db1grid, db2grid, iatt1, distType, dmax, tab))
          return 1;
      }
    }
    else
    {
      // Point to Grid
      if (flag_fill)
      {
        // Point to Grid (flag_fill = TRUE)
        if (flag_ball)
        {
          // Note that we do not benefit from the fact that db2 is a Grid
          if (_expandPointToPointBall(db1, db2grid, iatt1, distType, dmax, tab))
            return 1;
        }
        else
        {
          if (expandPointToGrid(db1, db2grid, iatt1, -1, -1, -1, -1, -1, 0,
                                distType, dmax, tab))
            return 1;
        }
      }
      else
      {
        // Point to Grid (flag_fill = FALSE)
        // flag_ball option is not considered as it does not save time
        if (_migratePointToGrid(db1, db2grid, iatt1, distType, dmax, tab))
          return 1;
      }
    }
  }
  else if (db1->isGrid())
  {
    DbGrid* db1grid = dynamic_cast<DbGrid*>(db1);

    // Grid to Point
    if (flag_inter)
    {
      // Grid to Point (flag_inter = TRUE)
      if (_interpolateGridToPoint(db1grid, db2, iatt1, distType, dmax, tab))
        return 1;
    }
    else
    {
      // Grid to Point (flag_inter = FALSE)
      if (_migrateGridToPoint(db1grid, db2, iatt1, distType, dmax, tab))
        return 1;
    }
  }
  else
  {
    // Point to Point
    if (flag_ball)
    {
      if (_expandPointToPointBall(db1, db2, iatt1, distType, dmax, tab))
        return 1;
    }
    else
    {
      if (_expandPointToPoint(db1, db2, iatt1, distType, dmax, tab))
        return 1;
    }
  }

  // Store the resulting array in the output Db

  db2->setColumnByUID(tab, iatt2);
  return 0;
}

/*****************************************************************************/
/*!
 **  Migrate a variable from the point structure
 **  into a variable in the grid structure
 **
 ** \return  Error return code
 **
 ** \param[in]  db_point  Descriptor of the point parameters
 ** \param[in]  db_grid   Descriptor of the grid parameters
 ** \param[in]  iatt      Rank of the point attribute
 ** \param[in]  distType  Type of distance for calculating maximum distance
 **                       1 for L1 and 2 for L2 distance
 ** \param[in]  dmax      Array of maximum distances (optional)
 **
 ** \param[out] tab       Output array
 **
 *****************************************************************************/
int CalcMigrate::_migratePointToGrid(Db *db_point,
                                        DbGrid *db_grid,
                                        int iatt,
                                        int distType,
                                        const VectorDouble &dmax,
                                        VectorDouble &tab)
{
  if (!db_point->hasLargerDimension(db_grid)) return 1;
  int ndim_min = MIN(db_point->getNDim(), db_grid->getNDim());
  int ndim_max = MAX(db_point->getNDim(), db_grid->getNDim());

  /* Core allocation */

  VectorDouble local(db_point->getSampleNumber());
  VectorDouble dvect(ndim_max);
  VectorDouble coor(ndim_max);
  db_grid->getCoordinatesPerSampleInPlace(0, coor);

  /* Locate the samples on the grid */

  (void) st_locate_point_on_grid(db_point, db_grid, coor, local);

  /* Assign the index of the closest sample to each grid node */

  int nb_assign = 0;
  for (int iech = 0; iech < db_point->getSampleNumber(); iech++)
  {
    if (FFFF(local[iech])) continue;
    if (FFFF(db_point->getArray(iech, iatt))) continue;
    int inode = (int) local[iech];
    nb_assign++;
    if (FFFF(tab[inode]))
    {
      /* If the grid node is empty assign the sample to it */

      tab[inode] = iech;
    }
    else
    {
      /* If the grid is not empty, find the closest sample */

      int jech = (int) tab[inode];
      double dist1 = distance_inter(db_grid, db_point, inode, iech, dvect.data());
      if (st_larger_than_dmax(ndim_min, dvect, distType, dmax)) continue;
      double dist2 = distance_inter(db_grid, db_point, inode, jech, dvect.data());
      if (st_larger_than_dmax(ndim_min, dvect, distType, dmax)) continue;
      tab[inode] = (dist1 < dist2) ? iech : jech;
    }
  }
  if (OptDbg::query(EDbg::DB))
    message("Number of nodes directly assigned = %d/%d\n", nb_assign,
            db_grid->getSampleNumber());

  /* Convert into data values */

  for (int jnode = 0; jnode < db_grid->getSampleNumber(); jnode++)
  {
    if (FFFF(tab[jnode])) continue;
    tab[jnode] = db_point->getArray((int) tab[jnode], iatt);
  }
  return 0;
}

/*****************************************************************************/
/*!
 **  Expand a variable from the input structure
 **  into a variable in the output structure
 **
 ** \return  Error return code
 **
 ** \param[in]  db1       Descriptor of the input parameters
 ** \param[in]  db2       Descriptor of the output parameters
 ** \param[in]  iatt      Rank of the input attribute
 ** \param[in]  distType  Type of distance for calculating maximum distance
 **                       1 for L1 and 2 for L2 distance
 ** \param[in]  dmax      Array of maximum distances (optional)
 **
 ** \param[out] tab       Output array
 **
 ** \remarks Method is designed when the two 'Db' share the same space dimension
 **
 *****************************************************************************/
int CalcMigrate::_expandPointToPointBall(Db *db1,
                                         Db *db2,
                                         int iatt,
                                         int distType,
                                         const VectorDouble &dmax,
                                         VectorDouble &tab)
{
  if (! db1->hasSameDimension(db2)) return 1;
  int ndim = db1->getNDim();
  VectorDouble coor(ndim);
  VectorDouble dvect(ndim);

  // Establish the ball tree (on the grid)

  int leaf_size = 30;
  Ball ball(db1, leaf_size, 1);

  // Loop on the sample points

  for (int inode = 0; inode < db2->getSampleNumber(); inode++)
  {
    if (! db2->isActive(inode)) continue;
    db2->getCoordinatesPerSampleInPlace(inode, coor);
    int iech = ball.queryClosest(coor);

    if (! dmax.empty())
    {
      (void) distance_inter(db2, db1, inode, iech, dvect.data());
      if (st_larger_than_dmax(ndim, dvect, distType, dmax)) continue;
    }

    tab[inode] = db1->getArray(iech, iatt);
  }
  return 0;
}

/*****************************************************************************/
/*!
 **  Migrates a variable from the grid structure
 **  into a variable in another grid structure
 **
 ** \return  Error return code
 **
 ** \param[in]  db_gridin  descriptor of the grid parameters
 ** \param[in]  db_gridout descriptor of the point parameters
 ** \param[in]  iatt       rank of the grid attribute
 ** \param[in]  distType   Type of distance for calculating maximum distance
 **                        1 for L1 and 2 for L2 distance
 ** \param[in]  dmax       Array of maximum distances (optional)
 **
 ** \param[out]  tab      Output array
 **
 *****************************************************************************/
int CalcMigrate::_migrateGridToGrid(DbGrid *db_gridin,
                                    DbGrid *db_gridout,
                                    int iatt,
                                    int distType,
                                    const VectorDouble &dmax,
                                    VectorDouble &tab)
{
  if (!db_gridin->hasLargerDimension(db_gridout)) return 1;
  int ndim_min = MIN(db_gridin->getNDim(), db_gridout->getNDim());
  int ndim_max = MAX(db_gridin->getNDim(), db_gridout->getNDim());
  if (! db_gridin->isGrid()) return 1;
  if (! db_gridout->isGrid()) return 1;

  /* Core allocation */

  VectorDouble coor(ndim_max);
  VectorDouble dvect(ndim_max);
  VectorDouble dist(db_gridout->getSampleNumber(), 1.e30);

  // Initialize 'coor' as the first target sample
  db_gridout->rankToCoordinatesInPlace(0, coor);

  /* Loop on the input grid nodes */

  for (int iech = 0; iech < db_gridin->getSampleNumber(); iech++)
  {
    double value = db_gridin->getArray(iech, iatt);
    if (FFFF(value)) continue;

    /* Get the coordinates of the node from the input grid node */

    db_gridin->rankToCoordinatesInPlace(iech, coor);

    /* Locate in the output grid */

    int jech = db_gridout->coordinateToRank(coor);
    if (jech < 0) continue;
    double dist_loc = distance_inter(db_gridin, db_gridout, iech, jech, dvect.data());
    if (st_larger_than_dmax(ndim_min, dvect, distType, dmax)) continue;
    if (dist_loc > dist[jech]) continue;
    tab[jech] = value;
    dist[jech] = dist_loc;
  }
  return 0;
}

/*****************************************************************************/
/*!
 **  Expands a variable from one point Db
 **  into a variable into another point Db
 **
 ** \return  Error return code
 **
 ** \param[in]  db1       descriptor of the input parameters
 ** \param[in]  db2       descriptor of the output parameters
 ** \param[in]  iatt      rank of the input attribute
 ** \param[in]  distType  Type of distance for calculating maximum distance
 **                       1 for L1 and 2 for L2 distance
 ** \param[in]  dmax      Array of maximum distances (optional)
 **
 ** \param[out] tab       Output array
 **
 ** \remark: The method has been enlarged to give a valid answer even if
 ** \remark: Space Dimension of 'db2' >= Space Dimension of 'db1'
 ** \remark: In that case, 'dmax' must be defined for the smallest Space Dimension
 **
 *****************************************************************************/
int CalcMigrate::_expandPointToPoint(Db *db1,
                                    Db *db2,
                                    int iatt,
                                    int distType,
                                    const VectorDouble &dmax,
                                    VectorDouble &tab)
{
  if (!db1->hasLargerDimension(db2)) return 1;
  int ndim_min = MIN(db1->getNDim(), db2->getNDim());
  int ndim_max = MAX(db1->getNDim(), db2->getNDim());

  /* Core allocation (using the smallest possible space dimension: db2) */

  VectorDouble dvect(ndim_max);

  /* Loop on the output structure */

  for (int iech2 = 0; iech2 < db2->getSampleNumber(); iech2++)
  {
    if (!db2->isActive(iech2)) continue;

    /* Loop on the input structure */

    double distmin = 1.e30;
    int iechmin = -1;
    for (int iech1 = 0; iech1 < db1->getSampleNumber(); iech1++)
    {
      if (!db1->isActive(iech1)) continue;
      double dist = distance_inter(db1, db2, iech1, iech2, dvect.data());
      if (st_larger_than_dmax(ndim_min, dvect, distType, dmax)) continue;
      if (dist < distmin)
      {
        distmin = dist;
        iechmin = iech1;
      }
    }
    if (iechmin >= 0) tab[iech2] = db1->getArray(iechmin, iatt);
  }
  return 0;
}

/*****************************************************************************/
/*!
 **  Expands a variable from the grid structure
 **  into a variable in another grid structure
 **  All nodes of the Output Grid will be filled
 **
 ** \return  Error return code
 **
 ** \param[in]  db_gridin  descriptor of the grid parameters
 ** \param[in]  db_gridout descriptor of the point parameters
 ** \param[in]  iatt       rank of the grid attribute
 ** \param[in]  distType   Type of distance for calculating maximum distance
 **                        1 for L1 and 2 for L2 distance
 ** \param[in]  dmax       Array of maximum distance (optional)
 **
 ** \param[out]  tab      Output array
 **
 *****************************************************************************/
int CalcMigrate::_expandGridToGrid(DbGrid *db_gridin,
                                      DbGrid *db_gridout,
                                      int iatt,
                                      int distType,
                                      const VectorDouble &dmax,
                                      VectorDouble &tab)
{
  if (!db_gridin->hasLargerDimension(db_gridout)) return 1;
  if (!db_gridin->isGrid())
  {
    messerr("The 'db_gridin' file should be a Grid Db");
    return 1;
  }
  if (!db_gridout->isGrid())
  {
    messerr("The 'db_gridout' file should be a Grid Db");
    return 1;
  }
  int ndim_min = MIN(db_gridin->getNDim(), db_gridout->getNDim());
  int ndim_max = MAX(db_gridin->getNDim(), db_gridout->getNDim());

  /* Core allocation */

  VectorDouble coor(ndim_max);
  VectorDouble dvect(ndim_max);
  VectorDouble dist(db_gridout->getSampleNumber());
  for (int jech = 0; jech < db_gridout->getSampleNumber(); jech++)
    dist[jech] = 1.e30;

  /* Loop on the output grid nodes */

  for (int iech = 0; iech < db_gridout->getSampleNumber(); iech++)
  {
    if (!db_gridout->isActive(iech)) continue;

    db_gridout->rankToCoordinatesInPlace(iech, coor);
    int jech = db_gridin->coordinateToRank(coor);
    if (jech < 0) continue;

    double dist_loc = distance_inter(db_gridin, db_gridout, jech, iech,
                                     dvect.data());
    if (st_larger_than_dmax(ndim_min, dvect, distType, dmax)) continue;
    if (dist_loc > dist[iech]) continue;
    tab[iech] = db_gridin->getArray(jech, iatt);
    dist[iech] = dist_loc;
  }
  return 0;
}

/*****************************************************************************/
/*!
 ** Interpolate the value on a Db by interpolating the data
 ** from a regular grid Db
 **
 ** \return  Error return code
 **
 ** \param[in]  db_grid   descriptor of the grid parameters
 ** \param[in]  db_point  descriptor of the point parameters
 ** \param[in]  iatt      rank of the grid attribute
 ** \param[in]  distType  Type of distance for calculating maximum distance
 **                       1 for L1 and 2 for L2 distance
 ** \param[in]  dmax      Array of maximum distances (optional)
 **
 ** \param[out]  tab      Output array
 **
 ** \remark A point which does not lie between two valuated grid nodes
 ** \remark (in all space dimensions) is always set to FFFF
 **
 *****************************************************************************/
int CalcMigrate::_interpolateGridToPoint(DbGrid *db_grid,
                                         Db *db_point,
                                         int iatt,
                                         int distType,
                                         const VectorDouble &dmax,
                                         VectorDouble &tab)
{
  if (!db_grid->hasLargerDimension(db_point)) return 1;

  /* Core allocation */

  VectorDouble coor(db_point->getNDim());

  /* Loop on the point samples */

  for (int iech = 0; iech < db_point->getSampleNumber(); iech++)
  {
    if (!db_point->isActive(iech)) continue;
    db_point->getCoordinatesPerSampleInPlace(iech, coor);
    tab[iech] = st_multilinear_interpolation(db_grid, iatt, distType, dmax, coor);
  }
  return 0;
}
