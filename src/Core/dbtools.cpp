/* COPYRIGHT ARMINES, ALL RIGHTS RESERVED                                     */
/*                                                                            */
/* THE CONTENT OF THIS WORK CONTAINS CONFIDENTIAL AND PROPRIETARY             */
/* INFORMATION OF ARMINES. ANY DUPLICATION, MODIFICATION,                     */
/* DISTRIBUTION, OR DISCLOSURE IN ANY FORM, IN WHOLE, OR IN PART, IS STRICTLY */
/* PROHIBITED WITHOUT THE PRIOR EXPRESS WRITT%EN PERMISSION OF ARMINES        */
/*                                                                            */
/* TAG_SOURCE_CG                                                              */
/******************************************************************************/
#include "geoslib_e.h"
#include "Mesh/MeshETurbo.hpp"
#include "LinearOp/ShiftOpCs.hpp"
#include "LinearOp/PrecisionOp.hpp"
#include "LinearOp/ProjMatrix.hpp"
#include "LinearOp/OptimCostColored.hpp"
#include "Stats/Classical.hpp"
#include "Morpho/Morpho.hpp"
#include "Basic/NamingConvention.hpp"
#include "Basic/Utilities.hpp"
#include "Basic/String.hpp"
#include "Covariances/CovAniso.hpp"
#include "Model/ANoStat.hpp"
#include "Model/NoStatArray.hpp"
#include "Db/Db.hpp"
#include "Basic/Law.hpp"

/*! \cond */
#define TRACE(i,iseg)       (trace[(i) * nseg + (iseg)])
#define TAB(iech,idim)      (tab[ndim * (iech) + (idim)])
#define LINE(nbline,i)      (line[npline * (nbline) + (i)])
#define PROP1(iz,iprop)     (prop1[(iz) * nprop + (iprop)])
#define PROP2(iz,iprop)     (prop2[(iz) * nprop + (iprop)])
#define WTAB(iz,icode,ivar) (wtab[(ivar) + nvar * ((iz) + nz * (icode))])
#define WCOR(iz,icode,idim) (wcor[(idim) + ndim * ((iz) + nz * (icode))])
#define WCNT(iz,icode)      (wcnt[                 (iz) + nz * (icode)])
#define RES(nval,idim)      (res[(idim) + (ndim+1) * (nval)])
#define R(i,j)              (R[(i) * n + (j)])

typedef struct
{
  char key[5];
  char title[STRING_LENGTH];
  int flag_rank;
  int flag_bounds;
} Edit_Item;
static int N_EDIT = 10;
static Edit_Item EDIT[] = { { "P", "Define the Properties", 0, 0 },
                            { ".", "Same as previous Command", 0, 0 },
                            { "S", "Shift the Sample Rank", 1, 0 },
                            { "V", "Shift the Variable Rank", 1, 0 },
                            { "AS", "Set the Sample Rank", 1, 0 },
                            { "AV", "Set the Variable Rank", 1, 0 },
                            { "M", "Modify the current Value", 0, 0 },
                            { "FD", "Find Next value in Interval", 0, 1 },
                            { "FU", "Find Next value in Interval", 0, 1 },
                            { "D", "Current Display", 0, 0 } };

static Db *DB_GRID_FILL;
/*! \endcond */

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
                                   const Db *db_grid,
                                   VectorDouble& coor,
                                   VectorDouble& tab)
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
 ** \param[in]  np        Number of discretized points
 ** \param[in]  xp        Array of first coordinates
 ** \param[in]  yp        Array of second coordinates
 ** \param[in]  zp        Array of third coordinates
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
                                  double *xp,
                                  double *yp,
                                  double *zp,
                                  const Db *db_grid,
                                  double *tab)
{
  int iech, iad, number;
  VectorDouble coor(3);

  /* Initializations */

  for (int idim = 0; idim < 3; idim++) coor[idim] = 0.;

  /* Loop on the point samples */

  for (iech = number = 0; iech < np; iech++)
  {
    tab[iech] = TEST;
    if (xp != (double *) NULL) coor[0] = xp[iech];
    if (yp != (double *) NULL) coor[1] = yp[iech];
    if (zp != (double *) NULL) coor[2] = zp[iech];
    iad = db_grid->coordinateToRank(coor);
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
 ** \param[in]  ldmax     Type of distance for calculating maximum distance
 **                       1 for L1 and 2 for L2 distance
 ** \param[in]  dmax      Array of maximum distances (optional)
 **
 ** \remarks If 'dmax' is not defined, the function always returns 0
 ** \remarks The function returns 1 as soon as one vector component is
 ** \remarks larger than the maximum value for this direction
 **
 *****************************************************************************/
static int st_larger_than_dmax(int ndim,
                               const VectorDouble& dvect,
                               int ldmax,
                               const VectorDouble& dmax)
{
  double ratio, rtot;

  if (dmax.empty()) return (0);

  /* Dispatch according to the type of distance */

  if (ldmax == 1)
  {

    // L1 distance

    for (int idim = 0; idim < ndim; idim++)
    {
      if (ABS(dvect[idim]) > dmax[idim]) return (1);
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
    if (rtot > 1) return (1);
  }

  return (0);
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
 ** \param[in]  ldmax     Type of distance for calculating maximum distance
 **                       1 for L1 and 2 for L2 distance
 ** \param[in]  dmax      Array of maximum distances (optional)
 **
 ** \param[out]  tab      Output array (Dimension: number of samples in db_point)
 **
 *****************************************************************************/
static int st_migrate_grid_to_point(Db *db_grid,
                                    Db *db_point,
                                    int iatt,
                                    int ldmax,
                                    const VectorDouble& dmax,
                                    VectorDouble& tab)
{
  if (! db_grid->hasLargerDimension(db_point)) return 1;
  int ndim_min = MIN(db_grid->getNDim(),db_point->getNDim());
  int ndim_max = MIN(db_grid->getNDim(),db_point->getNDim());
  VectorDouble dvect(ndim_max);
  VectorDouble coor(ndim_max);

  /* Define the default values for 'coor'*/

  db_point->getCoordinate(0,coor);

  /* Locate the samples on the grid */

  (void) st_locate_point_on_grid(db_point, db_grid, coor, tab);

  /* Loop on the point samples */

  for (int iech = 0; iech < db_point->getSampleNumber(); iech++)
  {
    if (FFFF(tab[iech])) continue;
    int rank = (int) tab[iech];
    if (!dmax.empty())
    {
      (void) distance_inter(db_grid, db_point, rank, iech, dvect.data());
      if (st_larger_than_dmax(ndim_min, dvect, ldmax, dmax)) continue;
    }
    tab[iech] = db_grid->getArray(rank, iatt);
  }
  return 0;
}

/*****************************************************************************/
/*!
 **  Find the Space dimension for a set of coordinate arrays (limited to 3-D)
 **
 ** \return  Space Dimension
 **
 ** \param[in]  xp        Array of first coordinates
 ** \param[in]  yp        Array of second coordinates
 ** \param[in]  zp        Array of third coordinates
 **
 *****************************************************************************/
static int st_get_ndim(double *xp, double *yp, double *zp)
{
  int ndim = 0;
  if (xp != (double *) NULL) ndim = 1;
  if (yp != (double *) NULL) ndim = 2;
  if (zp != (double *) NULL) ndim = 3;
  return ndim;
}

/*****************************************************************************/
/*!
 **  Migrates a variable from the grid structure
 **  into a variable at points defined by coordinate vectors (maximum 3D)
 **
 ** \return  Error return code
 **
 ** \param[in]  db_grid   descriptor of the grid parameters
 ** \param[in]  iatt      rank of the grid attribute
 ** \param[in]  np        Number of discretized points
 ** \param[in]  xp        Array of first coordinates
 ** \param[in]  yp        Array of second coordinates
 ** \param[in]  zp        Array of third coordinates
 **
 ** \param[out]  tab      Output array (Dimension: number of discretized points)
 **
 *****************************************************************************/
GEOSLIB_API int migrate_grid_to_coor(const Db *db_grid,
                                     int iatt,
                                     int np,
                                     double *xp,
                                     double *yp,
                                     double *zp,
                                     double *tab)
{
  int ndim = st_get_ndim(xp, yp, zp);
  if (db_grid->getNDim() != ndim)
  {
    messerr("The Space Dimension of the First Db (%d)", db_grid->getNDim());
    messerr("must be equal to the Space Dimension of the coordinate arrays",
            ndim);
    return 1;
  }

  /* Locate the samples on the grid */

  (void) st_locate_coor_on_grid(np, xp, yp, zp, db_grid, tab);

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
 **  Migrate a variable from the point structure
 **  into a variable in the grid structure
 **
 ** \return  Error return code
 **
 ** \param[in]  db_point  Descriptor of the point parameters
 ** \param[in]  db_grid   Descriptor of the grid parameters
 ** \param[in]  iatt      Rank of the point attribute
 ** \param[in]  ldmax     Type of distance for calculating maximum distance
 **                       1 for L1 and 2 for L2 distance
 ** \param[in]  dmax      Array of maximum distances (optional)
 **
 ** \param[out] tab       Output array
 **
 *****************************************************************************/
static int st_migrate_point_to_grid(Db *db_point,
                                    Db *db_grid,
                                    int iatt,
                                    int ldmax,
                                    const VectorDouble& dmax,
                                    VectorDouble& tab)
{
  double dist1, dist2;
  int iech, jech, inode, error, nb_assign, ndim_min, ndim_max;
  VectorDouble dvect, coor, local;

  /* Initializations */

  error = 1;

  /* Preliminary checks */

  if (! db_point->hasLargerDimension(db_grid)) goto label_end;
  ndim_min = MIN(db_point->getNDim(),db_grid->getNDim());
  ndim_max = MIN(db_point->getNDim(),db_grid->getNDim());

  /* Core allocation */

  local.resize(db_point->getSampleNumber());
  dvect.resize(ndim_max);
  coor.resize(ndim_max);
  db_grid->getCoordinate(0,coor);

  /* Locate the samples on the grid */

  (void) st_locate_point_on_grid(db_point, db_grid, coor, local);

  /* Assign the index of the closest sample to each grid node */

  for (iech = nb_assign = 0; iech < db_point->getSampleNumber(); iech++)
  {
    if (FFFF(local[iech])) continue;
    if (FFFF(db_point->getArray(iech, iatt))) continue;
    inode = (int) local[iech];
    nb_assign++;
    if (FFFF(tab[inode]))
    {
      /* If the grid node is empty assign the sample to it */

      tab[inode] = iech;
    }
    else
    {
      /* If the grid is not empty, find the closest sample */

      jech = (int) tab[inode];
      dist1 = distance_inter(db_grid, db_point, inode, iech, dvect.data());
      if (st_larger_than_dmax(ndim_min, dvect, ldmax, dmax)) continue;
      dist2 = distance_inter(db_grid, db_point, inode, jech, dvect.data());
      if (st_larger_than_dmax(ndim_min, dvect, ldmax, dmax)) continue;
      tab[inode] = (dist1 < dist2) ? iech : jech;
    }
  }
  if (debug_query("db"))
    message("Number of nodes directly assigned = %d/%d\n", nb_assign,
            db_grid->getSampleNumber());

  /* Convert into data values */

  for (inode = 0; inode < db_grid->getSampleNumber(); inode++)
  {
    if (FFFF(tab[inode])) continue;
    tab[inode] = db_point->getArray((int) tab[inode], iatt);
  }

  /* Set the error return code */

  error = 0;

  /* Core deallocation */

  label_end:
  return (error);
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
 ** \param[in]  ldmax      Type of distance for calculating maximum distance
 **                        1 for L1 and 2 for L2 distance
 ** \param[in]  dmax       Array of maximum distances (optional)
 **
 ** \param[out]  tab      Output array
 **
 *****************************************************************************/
static int st_migrate_grid_to_grid(Db *db_gridin,
                                   Db *db_gridout,
                                   int iatt,
                                   int ldmax,
                                   const VectorDouble& dmax,
                                   VectorDouble& tab)
{
  int error, iech, jech, ndim_min, ndim_max;
  double *dist, dist_loc, value;
  VectorDouble dvect, coor;

  /* Initializations */

  error = 1;
  dist = (double *) NULL;

  /* Preliminary checks */

  if (! db_gridin->hasLargerDimension(db_gridout)) goto label_end;
  ndim_min = MIN(db_gridin->getNDim(),db_gridout->getNDim());
  ndim_max = MAX(db_gridin->getNDim(),db_gridout->getNDim());
  if (!is_grid(db_gridin)) goto label_end;
  if (!is_grid(db_gridout)) goto label_end;

  /* Core allocation */

  coor.resize(ndim_max);
  dvect.resize(ndim_max);
  dist = (double *) mem_alloc(sizeof(double) * db_gridout->getSampleNumber(), 0);
  if (dist == (double *) NULL) goto label_end;
  for (jech = 0; jech < db_gridout->getSampleNumber(); jech++) dist[jech] = 1.e30;

  // Initialize 'coor' as the first target sample
  db_gridout->rankToCoordinate(0, coor);

  /* Loop on the input grid nodes */

  for (iech = 0; iech < db_gridin->getSampleNumber(); iech++)
  {
    value = db_gridin->getArray(iech, iatt);
    if (FFFF(value)) continue;

    /* Get the coordinates of the node from the input grid node */

    db_gridin->rankToCoordinate(iech, coor);

    /* Locate in the output grid */

    jech = db_gridout->coordinateToRank(coor);
    if (jech < 0) continue;
    dist_loc = distance_inter(db_gridin, db_gridout, iech, jech, dvect.data());
    if (st_larger_than_dmax(ndim_min, dvect, ldmax, dmax)) continue;
    if (dist_loc > dist[jech]) continue;
    tab[jech] = value;
    dist[jech] = dist_loc;
  }

  /* Set the error return code */

  error = 0;

  label_end:
  dist = (double *) mem_free((char * ) dist);
  return (error);
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
 ** \param[in]  ldmax     Type of distance for calculating maximum distance
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
static int st_expand_point_to_point(Db *db1,
                                    Db *db2,
                                    int iatt,
                                    int ldmax,
                                    const VectorDouble& dmax,
                                    VectorDouble& tab)
{
  if (! db1->hasLargerDimension(db2)) return 1;
  int ndim_min = MIN(db1->getNDim(),db2->getNDim());
  int ndim_max = MAX(db1->getNDim(),db2->getNDim());

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
      if (! db1->isActive(iech1)) continue;
      double dist = distance_inter(db1, db2, iech1, iech2, dvect.data());
      if (st_larger_than_dmax(ndim_min, dvect, ldmax, dmax)) continue;
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
 **  Expands a variable from one point Db
 **  into a variable at points defined by coordinate vectors (maximum 3D)
 **
 ** \return  Error return code
 **
 ** \param[in]  db1   descriptor of the input parameters
 ** \param[in]  iatt  rank of the input attribute
 **
 ** \param[in]  np        Number of discretized points
 ** \param[in]  xp        Array of first coordinates
 ** \param[in]  yp        Array of second coordinates
 ** \param[in]  zp        Array of third coordinates
 **
 ** \param[out]  tab      Output array (Dimension: number of discretized points)
 **
 *****************************************************************************/
GEOSLIB_API int expand_point_to_coor(const Db *db1,
                                     int iatt,
                                     int np,
                                     double *xp,
                                     double *yp,
                                     double *zp,
                                     double *tab)
{
  double* tab1 = (double *) NULL;
  double* tab2 = (double *) NULL;

  /* Preliminary checks */

  int ndim = db1->getNDim();
  if (ndim > 3)
  {
    messerr("Function 'expand_point_to_coor' is limited to space <= 3");
    return 1;
  }
  int ndimp = st_get_ndim(xp, yp, zp);
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

    if (ndim >= 1) tab2[0] = xp[iech2];
    if (ndim >= 2) tab2[1] = yp[iech2];
    if (ndim >= 3) tab2[2] = zp[iech2];

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
 **  Expands a variable from the grid structure
 **  into a variable in another grid structure
 **  All nodes of the Output Grid will be filled
 **
 ** \return  Error return code
 **
 ** \param[in]  db_gridin  descriptor of the grid parameters
 ** \param[in]  db_gridout descriptor of the point parameters
 ** \param[in]  iatt       rank of the grid attribute
 ** \param[in]  ldmax      Type of distance for calculating maximum distance
 **                        1 for L1 and 2 for L2 distance
 ** \param[in]  dmax       Array of maximum distance (optional)
 **
 ** \param[out]  tab      Output array
 **
 *****************************************************************************/
static int st_expand_grid_to_grid(Db *db_gridin,
                                  Db *db_gridout,
                                  int iatt,
                                  int ldmax,
                                  const VectorDouble& dmax,
                                  VectorDouble& tab)
{
  if (! db_gridin->hasLargerDimension(db_gridout)) return 1;
  if (! is_grid(db_gridin,true)) return 1;
  if (! is_grid(db_gridout,true)) return 1;
  int ndim_min = MIN(db_gridin->getNDim(),db_gridout->getNDim());
  int ndim_max = MAX(db_gridin->getNDim(),db_gridout->getNDim());

  /* Core allocation */

  VectorDouble coor(ndim_max);
  VectorDouble dvect(ndim_max);
  VectorDouble dist(db_gridout->getSampleNumber());
  for (int jech = 0; jech < db_gridout->getSampleNumber(); jech++) dist[jech] = 1.e30;

  /* Loop on the output grid nodes */

  for (int iech = 0; iech < db_gridout->getSampleNumber(); iech++)
  {
    if (!db_gridout->isActive(iech)) continue;

    db_gridout->rankToCoordinate(iech, coor);
    int jech = db_gridin->coordinateToRank(coor);
    if (jech < 0) continue;

    double dist_loc = distance_inter(db_gridin, db_gridout, jech, iech, dvect.data());
    if (st_larger_than_dmax(ndim_min, dvect, ldmax, dmax)) continue;
    if (dist_loc > dist[iech]) continue;
    tab[iech] = db_gridin->getArray(jech, iatt);
    dist[iech] = dist_loc;
  }
  return 0;
}

/****************************************************************************/
/*!
 **  Look for duplicates
 **
 ** \return  Error return code
 **
 ** \param[in]  db1        First Db
 ** \param[in]  db2        Second Db
 ** \param[in]  flag_same  1 if the two Db files are the same
 ** \param[in]  flag_print 1 for verbose output
 ** \param[in]  opt_code   code selection option (if code is defined)
 ** \li                     0 : no use of the code selection
 ** \li                     1 : codes must be close enough
 ** \li                     2 : codes must be different
 ** \param[in]  tolcode    Code tolerance
 ** \param[in]  dist       Array of the minimum distance
 **
 ** \param[out]  sel       Array containing the selection
 **
 *****************************************************************************/
GEOSLIB_API int db_duplicate(Db *db1,
                             Db *db2,
                             int flag_same,
                             int flag_print,
                             int opt_code,
                             double tolcode,
                             double *dist,
                             double *sel)
{
  int idim, iech1, iech2, flag_diff, flag_code;
  double v1, v2;

  /* Initializations */

  flag_code = db1->hasCode() && db2->hasCode();

  /* Set the selection */

  for (iech2 = 0; iech2 < db2->getSampleNumber(); iech2++)
    sel[iech2] = 1;

  /* Loop on the samples of the second Db */

  for (iech2 = 0; iech2 < db2->getSampleNumber(); iech2++)
  {
    if (!db2->isActive(iech2)) continue;
    for (iech1 = 0; iech1 < db1->getSampleNumber(); iech1++)
    {
      if (!db1->isActive(iech1)) continue;
      if (flag_same)
      {
        if (iech1 == iech2) continue;
        if (!sel[iech1]) continue;
      }

      /* Check if the two points are collocated */

      for (idim = flag_diff = 0; idim < db1->getNDim() && flag_diff == 0; idim++)
      {
        v1 = db1->getCoordinate(iech1, idim);
        v2 = db2->getCoordinate(iech2, idim);
        if (flag_code)
        {
          if (code_comparable(db1, db2, iech1, iech2, opt_code, (int) tolcode))
            continue;
        }
        if (ABS(v1 - v2) > dist[idim]) flag_diff = 1;
      }
      if (flag_diff) continue;

      sel[iech2] = 0;

      /* Optional printout */

      if (flag_print)
      {
        message("Sample %d too close to sample %d\n", iech1 + 1, iech2 + 1);
        db_sample_print(db1, iech1, 1, 0, 0);
        db_sample_print(db2, iech2, 1, 0, 0);
        message("\n");
      }
    }
  }

  return 0;
}

/****************************************************************************/
/*!
 **  Calculate the (discretized) surface of influence
 **
 ** \return  Error returned code
 **
 ** \param[in]  db_point Db containing the data points
 ** \param[in]  db_grid  Db containing the discretization grid
 ** \param[in]  icol     Rank of the variable whose value at a data
 **                      is attached to the grid nodes of the influence polygon
 ** \param[in]  dlim     Maximum distance (TEST if not defined)
 **
 ** \param[out]  dtab    Array containing the surface of influence
 **                      (Dimension = Number of samples in db_point)
 ** \param[out]  gtab    Array containing the surface of influence of the
 **                      polygon to which it belongs (or TEST)
 **                      (Dimension = Number of samples in db_grid)
 **
 *****************************************************************************/
GEOSLIB_API int surface(Db *db_point,
                        Db *db_grid,
                        int icol,
                        double dlim,
                        double *dtab,
                        double *gtab)
{
  double v1, v2, dist, d2min, d2max, delta, maille;
  int idim, iech, igrid, ndim;

  /* Initializations */

  if (! db_grid->hasSameDimension(db_point)) return (1);
  ndim = db_point->getNDim();

  /* Preliminary calculations */

  d2max = (FFFF(dlim)) ? 1.e30 :
                         dlim * dlim;
  maille = db_grid->getCellSize();
  for (iech = 0; iech < db_point->getSampleNumber(); iech++)
    dtab[iech] = 0.;

  /* Loop on the target points */

  for (igrid = 0; igrid < db_grid->getSampleNumber(); igrid++)
  {
    gtab[igrid] = -1;
    if (!db_grid->isActive(igrid)) continue;

    /* Loop on the data points */

    d2min = d2max;
    for (iech = 0; iech < db_point->getSampleNumber(); iech++)
    {
      if (!db_point->isActive(iech)) continue;

      /* Calculate the distance between node and data */

      dist = 0.;
      for (idim = 0; idim < ndim && !FFFF(dist); idim++)
      {
        v1 = db_grid->getCoordinate(igrid, idim);
        v2 = db_point->getCoordinate(iech, idim);
        if (FFFF(v1) || FFFF(v2)) dist = TEST;
        delta = v1 - v2;
        dist += delta * delta;
      }
      if (FFFF(dist)) continue;
      if (dist > d2max) continue;

      /* Keep the closest sample */

      if (dist > d2min) continue;
      gtab[igrid] = iech;
      d2min = dist;
    }
  }

  /* Calculate the influence of each datum */

  for (igrid = 0; igrid < db_grid->getSampleNumber(); igrid++)
  {
    iech = (int) gtab[igrid];
    if (iech >= 0) dtab[iech]++;
  }
  for (iech = 0; iech < db_point->getSampleNumber(); iech++)
    dtab[iech] *= maille;

  /* Evaluate each grid node with the size of the influence polygon */
  /* to which it belongs                                            */

  for (igrid = 0; igrid < db_grid->getSampleNumber(); igrid++)
  {
    iech = (int) gtab[igrid];
    if (iech >= 0)
      gtab[igrid] = dtab[iech];
    else
      gtab[igrid] = TEST;
  }

  return (0);
}

/****************************************************************************/
/*!
 **  Print the Db contents
 **
 ** \param[in]  db    Db descriptor
 ** \param[in]  nrdv  Radius of the variable display
 ** \param[in]  nrds  Radius of the sample display
 ** \param[in]  ivar  Rank of the Target Variable
 ** \param[in]  iech  Rank of the Target Sample
 **
 *****************************************************************************/
static void st_edit_display(Db *db, int nrdv, int nrds, int ivar, int iech)
{
  int item, nvar, nech, ivar_deb, ivar_fin, iech_deb, iech_fin, jvar, jech;
  ENUM_LOCS locatorType;
  char string[5];

  /* Initializations */

  (void) strcpy(string, "NA");
  nech = db->getSampleNumber();
  nvar = db->getFieldNumber();

  ivar_deb = ivar - nrdv;
  ivar_fin = ivar + nrdv;
  if (ivar_deb < 0)
  {
    ivar_deb = 0;
    ivar_fin = MIN(2 * nrdv, nvar - 1);
  }
  if (ivar_fin >= nvar)
  {
    ivar_fin = nvar - 1;
    ivar_deb = MAX(0, ivar_fin - 2 * nrdv);
  }

  iech_deb = iech - nrds;
  iech_fin = iech + nrds;
  if (iech_deb < 0)
  {
    iech_deb = 0;
    iech_fin = MIN(2 * nrds, nech - 1);
  }
  if (iech_fin >= nech)
  {
    iech_fin = nech - 1;
    iech_deb = MAX(0, iech_fin - 2 * nrds);
  }

  /* Print the Header (Variable name) */

  tab_prints(NULL, 1, GD_J_RIGHT, " ");
  for (jvar = ivar_deb; jvar <= ivar_fin; jvar++)
  {
    if (db->getLocatorByColumn(jvar, &locatorType, &item))
    {
      String strloc = getLocatorName(locatorType, item);
      (void) strcpy(string, strloc.c_str());
    }
    else
      (void) strcpy(string, "NA");
    if (jvar == ivar) (void) strcat(string, "*");
    tab_prints(NULL, 1, GD_J_RIGHT, string);
  }
  message("\n");

  /* Print the Header (Variable rank) */

  tab_prints(NULL, 1, GD_J_RIGHT, " ");
  for (jvar = ivar_deb; jvar <= ivar_fin; jvar++)
    tab_print_rc(NULL, 1, GD_J_RIGHT, 2, jvar + 1);
  message("\n");

  /* Loop on the samples */

  for (jech = iech_deb; jech <= iech_fin; jech++)
  {
    tab_print_rc(NULL, 1, GD_J_RIGHT, 3, jech + 1);
    if (iech == jech)
      message("*");
    else
      message(" ");
    for (jvar = ivar_deb; jvar <= ivar_fin; jvar++)
      tab_printg(NULL, 1, GD_J_RIGHT, db->getArray(jech, jvar));
    message("\n");
  }
  return;
}

/****************************************************************************/
/*!
 **  Look for the next sample with a value within the interval
 **
 ** \return  Rank of the next sample
 **
 ** \param[in]  db      Db descriptor
 ** \param[in]  iech    Rank of the current sample
 ** \param[in]  ivar    Rank of the current variable
 ** \param[in]  orient  Orientation (1: downwards; -1: forwards)
 ** \param[in]  vmin    Minimum value
 ** \param[in]  vmax    Maximum value
 **
 *****************************************************************************/
static int st_edit_find(Db *db,
                        int iech,
                        int ivar,
                        int orient,
                        double vmin,
                        double vmax)
{
  double value;
  int i;

  /* Dispatch */

  if (orient > 0)
  {
    for (i = iech + 1; i < db->getSampleNumber(); i++)
    {
      value = db->getArray(i, ivar);
      if (FFFF(value)) continue;
      if (!FFFF(vmin) && value < vmin) continue;
      if (!FFFF(vmax) && value > vmax) continue;
      return (i);
    }
    messerr("--> String not found before the end-of-file");
    return (iech);
  }
  else
  {
    for (i = iech - 1; i >= 0; i--)
    {
      value = db->getArray(i, ivar);
      if (FFFF(value)) continue;
      if (!FFFF(vmin) && value < vmin) continue;
      if (!FFFF(vmax) && value > vmax) continue;
      return (i);
    }
    messerr("--> String not found before the top-of-file");
    return (iech);
  }
}

/****************************************************************************/
/*!
 **  Ask for the next keyword in the Editor
 **
 ** \return  Return code:
 ** \return  0 : A valid keyword has been found
 ** \return  1 : The 'stop' has been met
 ** \return -1 : The 'quit' has been met
 **
 ** \param[out]  item   Selected item
 ** \param[out]  rank   Value for the Shift
 ** \param[out]  vmin   Value for the lower bound
 ** \param[out]  vmax   Value for the upper bound
 **
 *****************************************************************************/
static int st_edit_ask(int *item, int *rank, double *vmin, double *vmax)
{
  int found, flag_skip, mem_long;
  char string[STRING_LENGTH], *decode;
  static int mem_item = 1;
  static int mem_rank = 1;
  static double mem_vmin = 0.;
  static double mem_vmax = 1.;

  label_loop: _lire_string("Enter Command (or 'stop' or 'quit' or '?')", 0,
                           NULL, string);

  /* Look for the string */

  found = -1;
  for (int i = 0; i < N_EDIT; i++)
    if (!strncasecmp(string, EDIT[i].key, strlen(EDIT[i].key))) found = i;

  /* Check for the special keyword */

  if (!strcasecmp(string, "STOP")) return (1);
  if (!strcasecmp(string, "QUIT")) return (-1);

  /* A valid keyword has not been found */

  if (found < 0)
  {
    mestitle(1, "List of the Valid Editor Keywords:");
    for (int i = 0; i < N_EDIT; i++)
      message("%2s : %s\n", EDIT[i].key, EDIT[i].title);
    goto label_loop;
  }

  /* A valid keyword has been encountered: Check for the rest */

  decode = &string[strlen(EDIT[found].key)];

  /* The 'same' command has been encountered */

  flag_skip = 0;
  if (found == 1)
  {
    found = mem_item;
    *rank = mem_rank;
    *vmin = mem_vmin;
    *vmax = mem_vmax;
    flag_skip = 1;
  }

  /* Ask for complementary information */

  if (!flag_skip)
  {
    /* A rank must be specified */

    if (EDIT[found].flag_rank)
    {
      string_strip_blanks(decode, 1);
      mem_long = static_cast<int> (strlen(decode));
      if (mem_long > 0)
      {
        *rank = strtol(decode, &decode, 0);
        if (mem_long == (int) strlen(decode))
        {
          messerr("Cannot convert '%s' into a valid Rank", decode);
          goto label_loop;
        }
      }
      else
        *rank = _lire_int("Value for the Shift", 1, mem_rank, ITEST, ITEST);
    }

    /* Bounds must be specified */

    if (EDIT[found].flag_bounds)
    {
      string_strip_blanks(decode, 1);
      mem_long = static_cast<int> (strlen(decode));
      if (mem_long > 0)
      {
        *vmin = strtod(decode, &decode);
        if (mem_long == (int) strlen(decode))
        {
          messerr("Cannot convert '%s' into a valid Minimum Bound", decode);
          goto label_loop;
        }
      }
      else
        *vmin = _lire_double("Minimum value", 1, mem_vmin, TEST, TEST);

      string_strip_blanks(decode, 1);
      mem_long = static_cast<int> (strlen(decode));
      if (mem_long > 0)
      {
        *vmax = strtod(decode, &decode);
        if (mem_long == (int) strlen(decode))
        {
          messerr("Cannot convert '%s' into a valid Maximum Bound", decode);
          goto label_loop;
        }
        if (*vmax < *vmin)
        {
          messerr("Upper bound (%lf) may not be smaller than Lower bound (%lf)",
                  (*vmax), (*vmin));
        }
      }
      else
        *vmax = _lire_double("Maximum value", 1, mem_vmax, *vmin, TEST);
    }
  }

  /* Return argument */

  *item = found;

  /* Store the answers for default values in next operation */

  mem_item = *item;
  mem_rank = *rank;
  mem_vmin = *vmin;
  mem_vmax = *vmax;

  return (0);
}

/****************************************************************************/
/*!
 **  Edit the Data Base Db
 **
 ** \return  Error return code
 **
 ** \param[in]  db   Db descriptor
 **
 ** \param[out] flag_valid: 1 for 'stop' and 0 for 'quit'
 **
 *****************************************************************************/
GEOSLIB_API int db_edit(Db *db, int *flag_valid)

{
  int nech, nvar, ivar, iech, incr, type, ok, nrds, nrdv, flag_inter;
  double vmin, vmax, value;

  /* Initializations */

  nech = db->getSampleNumber();
  nvar = db->getFieldNumber();
  ivar = iech = 0;
  ok = nrds = nrdv = incr = 1;
  vmin = vmax = TEST;
  if (nech < 1 || nvar < 1) return (1);

  while (ok)
  {
    st_edit_display(db, nrdv, nrds, ivar, iech);
    flag_inter = st_edit_ask(&type, &incr, &vmin, &vmax);
    if (flag_inter > 0)
    {
      *flag_valid = 1;
      ok = 0;
      break;
    }
    else if (flag_inter < 0)
    {
      *flag_valid = 0;
      ok = 0;
      break;
    }

    /* Dispatch */

    switch (type)
    {
      case 0: /* Set the Parameters */
        nrdv = _lire_int("Display Radius along Variable", 1, nrdv, 0, ITEST);
        nrds = _lire_int("Display Radius along Sample", 1, nrds, 0, ITEST);
        break;

      case 2: /* Relative Sample Rank */
        iech = MAX(0, MIN(iech + incr, nech - 1));
        break;

      case 3: /* Relative Variable Rank */
        ivar = MAX(0, MIN(ivar + incr, nvar - 1));
        break;

      case 4: /* Absolute Sample Rank */
        iech = MAX(0, MIN(incr - 1, nech - 1));
        break;

      case 5: /* Absolute Variable Rank */
        ivar = MAX(0, MIN(incr - 1, nvar - 1));
        break;

      case 6: /* Modify the Value */
        value = _lire_double("New value", 1, db->getArray(iech, ivar), TEST,TEST);
        db->setArray(iech, ivar, value);
        break;

      case 7: /* Next sample within an interval */
        iech = st_edit_find(db, iech, ivar, 1, vmin, vmax);
        break;

      case 8: /* Previous sample within an interval */
        iech = st_edit_find(db, iech, ivar, -1, vmin, vmax);
        break;

      case 9: /* Display the current selection */
        break;
    }
  }
  return (0);
}

/****************************************************************************/
/*!
 **  Normalize a set of variables
 **
 ** \return  Error returned code
 **
 ** \param[in]  db     Db structure
 ** \param[in]  oper   Name of the operator
 ** \li                "mean"  : Normalize the mean to the input value
 ** \li                "stdv"  : Normalize the st. dev. to the input value
 ** \li                "scal"  : Normalize the mean and st. dev.
 ** \li                "prop"  : Normalize the variables to proportions
 ** \param[in]  ncol   Number of variables
 ** \param[in]  cols   Ranks of the variables
 ** \param[in]  center Theoretical Mean value
 ** \param[in]  stdv   Theoretical Standard Deviation value
 **
 *****************************************************************************/
GEOSLIB_API int db_normalize(Db *db,
                             const char *oper,
                             int ncol,
                             int *cols,
                             double center,
                             double stdv)
{
  int iech, nech, icol, jcol, ndef, iptr;
  double *num, *mm, *vv, proptot, value;

  /* Initializations */

  nech = db->getSampleNumber();
  num = mm = vv = (double *) NULL;

  /* Check that all variables are defined */

  for (icol = 0; icol < ncol; icol++)
  {
    jcol = cols[icol];
    if (!db->isColumnIndexValid(jcol))
    {
      messerr("Column %d is not defined", cols[icol]);
      return (1);
    }
  }

  /* Core allocation */

  num = (double *) mem_alloc(sizeof(double) * ncol, 0);
  if (num == (double *) NULL) goto label_end;
  mm = (double *) mem_alloc(sizeof(double) * ncol, 0);
  if (mm == (double *) NULL) goto label_end;
  vv = (double *) mem_alloc(sizeof(double) * ncol, 0);
  if (vv == (double *) NULL) goto label_end;

  /* Initializations */

  for (icol = 0; icol < ncol; icol++)
    num[icol] = mm[icol] = vv[icol] = 0.;

  /* Printout */

  if (!strcmp(oper, "mean"))
  {
    if (FFFF(center)) center = 0.;
  }
  else if (!strcmp(oper, "stdv"))
  {
    if (FFFF(stdv)) stdv = 1.;
  }
  else if (!strcmp(oper, "scal"))
  {
    if (FFFF(center)) center = 0.;
    if (FFFF(stdv)) stdv = 1.;
  }
  else if (!strcmp(oper, "prop"))
  {
    center = 0.;                  // Useless line
  }
  else
  {
    messerr("Invalid operator name (%s)", oper);
    messerr("The list of operators available is:");
    messerr("mean  : Normalize the mean");
    messerr("stdv  : Normalize the st. dev.");
    messerr("scal  : Normalize the mean and st dev");
    messerr("prop  : Normalize the proportions");
    return (1);
  }

  /* Creation of the new attributes */

  iptr = db->addFields(ncol, TEST);
  if (iptr < 0) return (1);

  /* Loop on the samples */

  for (iech = 0; iech < nech; iech++)
  {
    if (!db->isActive(iech)) continue;

    /* Loop on the variables */

    ndef = 0;
    proptot = 0.;
    for (icol = 0; icol < ncol; icol++)
    {
      value = db->getArray(iech, cols[icol]);
      if (FFFF(value)) continue;
      if (!strcmp(oper, "prop")) value = MIN(1., MAX(0.,value));

      /* Update statistics */

      ndef += 1;
      num[icol] += 1;
      proptot += value;
      mm[icol] += value;
      vv[icol] += value * value;
    }

    /* Set the output array (in the case of proportions) */

    if (!strcmp(oper, "prop"))
    {
      for (icol = 0; icol < ncol; icol++)
      {
        value = db->getArray(iech, cols[icol]);
        value = MIN(1., MAX(0.,value));
        if (ndef == ncol && proptot > 0.)
          db->setArray(iech, iptr + icol, value / proptot);
      }
    }
  }

  if (strcmp(oper, "prop"))
  {

    /* Global Normation */

    for (icol = 0; icol < ncol; icol++)
    {
      if (num[icol] <= 0)
      {
        mm[icol] = TEST;
        vv[icol] = TEST;
      }
      else
      {
        mm[icol] /= num[icol];
        vv[icol] = vv[icol] / num[icol] - mm[icol] * mm[icol];
      }
    }

    /* Set the output array */

    for (iech = 0; iech < nech; iech++)
    {
      for (icol = 0; icol < ncol; icol++)
      {
        jcol = cols[icol];
        value = db->getArray(iech, jcol);
        if (!strcmp(oper, "mean"))
        {
          if (!FFFF(mm[icol]))
            db->setArray(iech, iptr + icol, value + center - mm[icol]);
        }
        else if (!strcmp(oper, "stdv"))
        {
          if (!FFFF(vv[icol]) && vv[icol] > 0.)
            db->setArray(iech, iptr + icol, value * stdv / sqrt(vv[icol]));
        }
        else if (!strcmp(oper, "scal"))
        {
          if (!FFFF(vv[icol]) && vv[icol] > 0. && !FFFF(mm[icol]))
            db->setArray(iech, iptr + icol,
                         center + (value - mm[icol]) * stdv / sqrt(vv[icol]));
        }
      }
    }
  }

  label_end: num = (double *) mem_free((char * ) num);
  mm = (double *) mem_free((char * ) mm);
  vv = (double *) mem_free((char * ) vv);

  return (0);
}

/****************************************************************************/
/*!
 **  Check if the cell is already filled
 **
 ** \return  1 if the cell (filled with facies) is already filled
 **
 ** \param[in]  ipos  Absolute grid index of the input grid node
 **
 *****************************************************************************/
static int st_grid_fill_already_filled(int ipos)

{
  int value;

  if (!DB_GRID_FILL->getSelection(ipos)) return (0);
  value = FFFF(DB_GRID_FILL->getVariable(ipos, 0)) ? 0 : 1;
  return (value);
}

/****************************************************************************/
/*!
 **  Check if the cell can be filled with fluid
 **
 ** \return  1 if the cell (filled with facies) can be filled with Fluid
 **
 ** \param[in]  ipos  Absolute grid index of the input grid node
 **
 *****************************************************************************/
static int st_grid_fill_to_be_filled(int ipos)

{
  int value;

  if (!DB_GRID_FILL->getSelection(ipos)) return (0);
  value = FFFF(DB_GRID_FILL->getVariable(ipos, 0)) ? 1 : 0;
  return (value);
}

/****************************************************************************/
/*!
 **  Find the neighborhood of the current cell
 **
 ** \return  1 if the neighborhood criterion is not fulfilled; 0 otherwise
 **
 ** \param[in]  ipos    Absolute grid index of the input grid node
 ** \param[in]  ndim    Space dimension
 ** \param[in]  radius  Radius of the neighborhood
 **
 ** \param[out] nech_loc Number of samples in the neighborhood
 ** \param[out] tabind   Index of the neighboring sample
 ** \param[out] tabval   Value of the neighboring sample
 **
 *****************************************************************************/
static void st_grid_fill_neigh(int ipos,
                               int ndim,
                               int radius,
                               int *nech_loc,
                               int *tabind,
                               double *tabval)
{
  double value;
  int i, ix, iy, iz, jpos, nech, iwork1[3], iwork2[3], nrx, nry, nrz, nmx, nmy,
      nmz;

  /* Initializations */

  nech = 0;
  nrx = (ndim >= 1) ? radius : 0;
  nry = (ndim >= 2) ? radius : 0;
  nrz = (ndim >= 3) ? radius : 0;
  nmx = (ndim >= 1) ? DB_GRID_FILL->getNX(0) : 1;
  nmy = (ndim >= 2) ? DB_GRID_FILL->getNX(1) : 1;
  nmz = (ndim >= 3) ? DB_GRID_FILL->getNX(2) : 1;

  /* Locate the central cell */

  for (i = 0; i < 3; i++)
    iwork1[i] = iwork2[i] = 0;
  db_index_sample_to_grid(DB_GRID_FILL, ipos, iwork1);

  /* Loop on the neighborhood cells */

  for (ix = -nrx; ix <= nrx; ix++)
  {
    iwork2[0] = iwork1[0] + ix;
    if (iwork2[0] < 0 || iwork2[0] >= nmx) continue;
    for (iy = -nry; iy <= nry; iy++)
    {
      iwork2[1] = iwork1[1] + iy;
      if (iwork2[1] < 0 || iwork2[1] >= nmy) continue;
      for (iz = -nrz; iz <= nrz; iz++)
      {
        iwork2[2] = iwork1[2] + iz;
        if (iwork2[2] < 0 || iwork2[2] >= nmz) continue;
        jpos = db_index_grid_to_sample(DB_GRID_FILL, iwork2);
        value = DB_GRID_FILL->getVariable(jpos, 0);
        if (FFFF(value)) continue;
        tabind[nech] = jpos;
        tabval[nech] = value;
        nech++;
      }
    }
  }

  /* Returning arguments */
  *nech_loc = nech;
  return;
}

/****************************************************************************/
/*!
 **  Calculate the extrapolation
 **
 ** \return  Error return code
 **
 ** \param[in]  ipos    Absolute grid index of the input grid node
 ** \param[in]  mode    Type of interpolation
 ** \li                 0 : Moving average
 ** \li                 1 : Inverse squared distance
 ** \li                 2 : Interpolation by a linear plane
 ** \li                 3 : Distance to the initial grains
 ** \param[in]  nech    Number of samples in the neighborhood
 ** \param[in]  tabind  Index of the neighboring sample
 ** \param[in]  tabval  Value of the neighboring sample
 **
 *****************************************************************************/
static int st_grid_fill_calculate(int ipos,
                                  int mode,
                                  int nech,
                                  int *tabind,
                                  double *tabval)
{
  double result, dist, dist2, top, bot, a[6], b[3], sol[3], f[4], dmin;
  int iech, indg[3], j, k, l, pivot;
  static int neq = 3;

  /* Initializations */

  result = top = bot = 0.;

  /* Dispatch according to the extrapolation mode */

  switch (mode)
  {
    case 0:
      for (iech = 0; iech < nech; iech++)
        top += tabval[iech];
      result = top / (double) nech;
      break;

    case 1:
      for (iech = 0; iech < nech; iech++)
      {
        dist = distance_intra(DB_GRID_FILL, ipos, tabind[iech], NULL);
        dist2 = dist * dist;
        top += tabval[iech] / dist2;
        bot += 1. / dist2;
      }
      result = top / bot;
      break;

    case 2:
      if (nech < 3) return (1);
      f[0] = 1.;
      for (iech = 0; iech < nech; iech++)
      {
        db_index_sample_to_grid(DB_GRID_FILL, tabind[iech], indg);
        grid_to_point(DB_GRID_FILL, indg, NULL, &f[1]);
        for (j = l = 0; j < neq; j++)
        {
          for (k = 0; k <= j; k++, l++)
            a[l] += f[j] * f[k];
          b[j] += tabval[iech] * f[j];
        }
      }
      if (matrix_solve(0, a, b, sol, neq, 1, &pivot)) return (1);
      if (pivot != 0) return (1);
      db_index_sample_to_grid(DB_GRID_FILL, ipos, indg);
      grid_to_point(DB_GRID_FILL, indg, NULL, &f[1]);
      for (j = 0; j < neq; j++)
        result += sol[j] * f[j];
      break;

    case 3:
      dmin = 0.;
      for (iech = 0; iech < nech; iech++)
      {
        dist = tabval[iech];
        if (dist > dmin) dmin = dist;
      }
      result = dmin + 1.;
      break;
  }

  /* Assign the result */

  DB_GRID_FILL->setVariable(ipos, 0, result);
  return (0);
}

/****************************************************************************/
/*!
 **  Fill an incomplete grid
 **
 ** \return  Error returned code
 **
 ** \param[in]  dbgrid  Db grid structure
 ** \param[in]  mode    Type of interpolation
 ** \li                  0 : Moving average
 ** \li                  1 : Inverse squared distance
 ** \li                  2 : Interpolation by a linear plane
 ** \li                  3 : Distance to the initial grains
 ** \param[in]  seed    Seed used for the random number generation
 ** \param[in]  radius  Radius of the neighborhood
 ** \param[in]  verbose Verbose flag
 ** \param[in]  namconv Naming convention
 **
 *****************************************************************************/
GEOSLIB_API int db_grid_fill(Db *dbgrid,
                             int mode,
                             int seed,
                             int radius,
                             bool verbose,
                             NamingConvention namconv)
{
  Skin   *skin;
  double *tabval;
  int    *tabind, error, rank, ipos, ndim, count, nech;

  /* Initializations */

  error = 1;

  /* Preliminary checks */

  if (!is_grid(dbgrid))
  {
    messerr("This function is limited to Grid Db");
    return (1);
  }
  if (! dbgrid->isVariableNumberComparedTo(1)) return(1);
  ndim = dbgrid->getNDim();
  if (ndim > 3)
  {
    messerr("This function is limited to a maximum 3-D space");
    return (1);
  }
  if (mode < 0 || mode > 3)
  {
    messerr("The argument 'mode' should lie between 0 and 3");
    return (1);
  }
  if (mode == 2 && radius < 1)
  {
    messerr("The linear interpolation requires a neighborhood radius > 1");
    return (1);
  }

  // Create the new variable and duplicate the Z-locator variable

  int iatt_in  = dbgrid->getAttribute(LOC_Z, 0);
  int iatt_out = dbgrid->addFields(1);
  dbgrid->duplicateColumnByAttribute(iatt_in, iatt_out);
  dbgrid->setLocatorByAttribute(iatt_out,LOC_Z);

  /* Global variables */

  DB_GRID_FILL = dbgrid;
  skin   = (Skin   *) NULL;
  tabval = (double *) NULL;
  tabind = (int    *) NULL;
  count  = (int) pow(2. * radius + 1., (double) ndim) - 1;

  /* Core allocation */

  law_set_random_seed(seed);
  tabval = (double *) mem_alloc(sizeof(double) * count, 0);
  if (tabval == (double *) NULL) goto label_end;
  tabind = (int *) mem_alloc(sizeof(int) * count, 0);
  if (tabind == (int *) NULL) goto label_end;
  skin   = skin_define(dbgrid, st_grid_fill_already_filled,st_grid_fill_to_be_filled, NULL);
  if (skin == (Skin *) NULL) goto label_end;

  if (skin_init(skin, verbose))
  {
    error = 0;
    goto label_end;
  }

  /* Implicit loop on the cells to be filled */

  while (skin_remains(skin))
  {

    /* Find the next cell to be processed */

    skin_next(skin, &rank, &ipos);

    /* Find the neighborhood */

    st_grid_fill_neigh(ipos, ndim, radius, &nech, tabind, tabval);

    /* Calculate the extrapolated value */

    if (st_grid_fill_calculate(ipos, mode, nech, tabind, tabval)) continue;

    /* Deduce the initial influence of the central cell */

    if (skin_unstack(skin, rank, ipos)) goto label_end;
  }

  // Optional printout

  if (verbose) skin_print(skin);

  /* Set the error return code */

  error = 0;
  namconv.setNamesAndLocators(dbgrid, LOC_Z, -1, dbgrid, iatt_out);

label_end:
  skin = skin_undefine(skin);
  tabval = (double *) mem_free((char * ) tabval);
  tabind = (int *) mem_free((char * ) tabind);
  return (error);
}

/****************************************************************************/
/*!
 **  Check consistency between the different bounds vectors
 **  and returns the number of classes
 **
 ** \return  Error return code
 **
 ** \param[in]  mini    Array containing the minima per class
 ** \param[in]  maxi    Array containing the maxima per class
 ** \param[in]  incmini Array containing the inclusive flag for minima per class
 ** \param[in]  incmaxi Array containing the inclusive flag for maxima per class
 **
 ** \param[out] nclass_arg Number of classes
 **
 *****************************************************************************/
static int st_check_bound_consistency(const VectorDouble& mini,
                                      const VectorDouble& maxi,
                                      const VectorBool& incmini,
                                      const VectorBool& incmaxi,
                                      int *nclass_arg)
{
  int nclass = 0;
  if (!mini.empty())
  {
    if (nclass > 0 && nclass != (int) mini.size())
    {
      messerr("Wrong dimension of 'mini'(%d). It should be %d", mini.size(),
              nclass);
      return 1;
    }
    nclass = static_cast<int> (mini.size());
  }
  if (!maxi.empty())
  {
    if (nclass > 0 && nclass != (int) maxi.size())
    {
      messerr("Wrong dimension of 'maxi'(%d). It should be %d", maxi.size(),
              nclass);
      return 1;
    }
    nclass = static_cast<int> (maxi.size());
  }
  if (!incmini.empty())
  {
    if (nclass > 0 && nclass != (int) incmini.size())
    {
      messerr("Wrong dimension of 'incmini'(%d). It should be %d",
              incmini.size(), nclass);
      return 1;
    }
    nclass = static_cast<int> (incmini.size());
  }
  if (!incmaxi.empty())
  {
    if (nclass > 0 && nclass != (int) incmaxi.size())
    {
      messerr("Wrong dimension of 'incmaxi'(%d). It should be %d",
              incmaxi.size(), nclass);
      return 1;
    }
    nclass = static_cast<int> (incmaxi.size());
  }
  if (nclass <= 0)
  {
    messerr("You must define at least one valid limit");
    return 1;
  }
  *nclass_arg = nclass;
  return 0;
}

/****************************************************************************/
/*!
 **  Convert the contents of the ivar-th continuous variable
 **  into a categorical array
 **
 ** \return  Error return code
 **
 ** \param[in]  db      Db structure
 ** \param[in]  iatt    Rank of the attribute
 ** \param[in]  mini    Array containing the minima per class
 ** \param[in]  maxi    Array containing the maxima per class
 ** \param[in]  incmini Array containing the inclusive flag for minima per class
 ** \param[in]  incmaxi Array containing the inclusive flag for maxima per class
 ** \param[in]  namconv Naming convention
 **
 ** \remark When the array mini or maxi are not provided, then
 ** \remark mini[iclass] = int(iclass-1)
 ** \remark maxi[iclass[ = int(iclass)
 **
 *****************************************************************************/
GEOSLIB_API int db_category(Db *db,
                            int iatt,
                            const VectorDouble& mini,
                            const VectorDouble& maxi,
                            const VectorBool& incmini,
                            const VectorBool& incmaxi,
                            NamingConvention namconv)
{
  // Determination of the number of classes

  int nclass;
  if (st_check_bound_consistency(mini, maxi, incmini, incmaxi, &nclass))
    return 1;

  /* Create the variable */

  int iptr = db->addFields(1, TEST);
  if (iptr < 0) return (1);

  /* Loop on the samples */

  for (int iech = 0; iech < db->getActiveSampleNumber(); iech++)
  {
    if (!db->isActive(iech)) continue;
    double value = db->getArray(iech, iatt);
    if (FFFF(value)) continue;

    /* Loop on the limits classes */

    int ival = 0;
    for (int iclass = 0; iclass < nclass; iclass++)
    {
      double minival = (mini.empty()) ? iclass + 0.5 : mini[iclass];
      double maxival = (maxi.empty()) ? iclass + 1.5 : maxi[iclass];
      if (!FFFF(minival))
      {
        int flag = (incmini.empty()) ? 1 :
                                       (int) incmini[iclass];
        if ((flag == 0 && value <= minival) || (flag == 1 && value < minival))
          continue;
      }
      if (!FFFF(maxival))
      {
        int flag = (incmaxi.empty()) ? 0 :
                                       (int) incmaxi[iclass];
        if ((flag == 0 && value >= maxival) || (flag == 1 && value > maxival))
          continue;
      }
      ival = iclass + 1;
    }

    /* Set the returning value */

    db->setArray(iech, iptr, (double) ival);
  }

  namconv.setNamesAndLocators(db, iatt, db, iptr);

  return (0);
}

/****************************************************************************/
/*!
 **  Create indicator variables
 **
 ** \return  Error returned code
 **
 ** \param[in]  db      Db structure
 ** \param[in]  iatt    Rank of the target variable
 ** \param[in]  flag_indic Type of variable(s) to be stored:
 **                     1 the indicator variable
 **                     0 the discretized mean variable per class
 ** \param[in]  mini    Array containing the minima per class
 ** \param[in]  maxi    Array containing the maxima per class
 ** \param[in]  incmini Array containing the inclusive flag for minima per class
 ** \param[in]  incmaxi Array containing the inclusive flag for maxima per class
 ** \param[in]  namconv Naming convention
 **
 ** \remark When both arrays mini and maxi are not provided, then:
 ** \remark mini[iclass] = iclass-1
 ** \remark maxi[iclass[ = iclass
 **
 *****************************************************************************/
GEOSLIB_API int db_indicator(Db *db,
                             int iatt,
                             int flag_indic,
                             const VectorDouble& mini,
                             const VectorDouble& maxi,
                             const VectorBool& incmini,
                             const VectorBool& incmaxi,
                             NamingConvention namconv)
{
  int nactive = 0;
  int iptr = 0;

  // Determination of the number of classes

  int nclass;
  if (st_check_bound_consistency(mini, maxi, incmini, incmaxi, &nclass))
    return 1;

  /* Core allocation */

  VectorInt count(nclass);
  VectorInt flagmin(nclass);
  VectorInt flagmax(nclass);
  VectorDouble mean(nclass);
  VectorDouble minival(nclass);
  VectorDouble maxival(nclass);
  for (int iclass = 0; iclass < nclass; iclass++)
  {
    count[iclass] = 0;
    mean[iclass] = 0.;
    minival[iclass] = (mini.empty()) ? iclass + 0.5 :
                                       mini[iclass];
    maxival[iclass] = (maxi.empty()) ? iclass + 1.5 :
                                       maxi[iclass];
    flagmin[iclass] = (incmini.empty()) ? 1 :
                                          (int) incmini[iclass];
    flagmax[iclass] = (incmaxi.empty()) ? 0 :
                                          (int) incmaxi[iclass];
  }
  int nbelow = 0;
  int nabove = 0;
  int mbelow = 0;
  int mabove = 0;

  /* Find extrema of all classes to sort too small or too large samples */

  double zmini = 1.e30;
  double zmaxi = -1.e30;
  for (int iclass = 0; iclass < nclass; iclass++)
  {
    if (!FFFF(minival[iclass])) zmini = MIN(zmini, minival[iclass]);
    if (!FFFF(maxival[iclass])) zmaxi = MAX(zmaxi, maxival[iclass]);
  }
  if (FFFF(zmini)) zmini = -1.e30;
  if (FFFF(zmaxi)) zmaxi = 1.e30;

  /* Create the variables */

  if (flag_indic)
  {
    iptr = db->addFields(nclass, 0.);
    if (iptr < 0) return 1;
  }
  else
  {
    iptr = db->addFields(1, TEST);
    if (iptr < 0) return 1;
  }

  /* Loop on the samples */

  for (int iech = 0; iech < db->getActiveSampleNumber(); iech++)
  {
    if (!db->isActive(iech)) continue;
    double value = db->getArray(iech, iatt);
    if (FFFF(value)) continue;

    /* Loop on the limit classes */

    int found = -1;
    for (int iclass = 0; iclass < nclass; iclass++)
    {
      int belong = 1;
      if (!FFFF(minival[iclass]))
      {
        if ((flagmin[iclass] == 0 && value <= minival[iclass]) || (flagmin[iclass]
            == 1
                                                                   && value < minival[iclass]))
          belong = 0;
      }
      if (!FFFF(maxival[iclass]))
      {
        if ((flagmax[iclass] == 0 && value >= maxival[iclass]) || (flagmax[iclass]
            == 1
                                                                   && value > maxival[iclass]))
          belong = 0;
      }

      /* Store the indicator (if required) */

      if (flag_indic) db->setArray(iech, iptr + iclass, belong);

      if (belong)
      {
        mean[iclass] += value;
        count[iclass]++;
        found = iclass;
      }
    }

    /* Update the statistics for classes below and above */

    int iclass;
    if (found < 0)
    {
      if (value < zmini)
      {
        nbelow += 1;
        iclass = -1;
        mbelow += (int) value;
      }
      else
      {
        nabove += 1;
        iclass = nclass;
        mabove += (int) value;
      }
    }
    else
    {
      iclass = found;
    }

    /* Store the rank of the class the sample belongs to */

    if (!flag_indic) db->setArray(iech, iptr, (double) iclass);
    nactive++;
  }

  /* Calculate the statistics */

  for (int iclass = 0; iclass < nclass; iclass++)
  {
    if (count[iclass] <= 0)
      mean[iclass] = TEST;
    else
      mean[iclass] /= (double) count[iclass];
  }
  if (nbelow > 0) mbelow = (int) (mbelow / (double) nbelow);
  if (nabove > 0) mabove = (int) (mabove / (double) nabove);

  /* Calculate the discretized variable */

  if (!flag_indic)
  {
    for (int iech = 0; iech < db->getActiveSampleNumber(); iech++)
    {
      double value = db->getArray(iech, iptr);
      if (FFFF(value)) continue;
      int iclass = (int) value;
      if (iclass < 0)
        db->setArray(iech, iptr, mbelow);
      else if (iclass < nclass)
        db->setArray(iech, iptr, mean[iclass]);
      else
        db->setArray(iech, iptr, mabove);
    }
  }

  // Naming convention
  if (flag_indic == 1)
    namconv.setNamesAndLocators(db, iatt, db, iptr, "Class", nclass);
  else
    namconv.setNamesAndLocators(db, iatt, db, iptr, "Mean", 1);

  return (0);
}

/*****************************************************************************/
/*!
 **  Select samples from a file according to the 2-D convex hull
 **  computed over the active samples of a second file
 **
 ** \return  Error returned code
 **
 ** \param[in]  db1     descriptor of the Db serving for convex hull calculation
 ** \param[in]  db2     descriptor of the Db where the mask must be performed
 ** \param[in]  verbose Verbose flag
 ** \param[in]  namconv Naming convention
 **
 ** \remark The Naming Convention locator Type is overwritten to LOC_SEL
 **
 *****************************************************************************/
GEOSLIB_API int db_selhull(Db *db1, Db *db2, bool verbose, NamingConvention namconv)
{
  Polygons* polygons = (Polygons *) NULL;

  // Create the variable in the output Db

  int isel = db2->addFields(1, 1.);

  /* Create the polygon as the convex hull of first Db */

  polygons = polygon_hull(db1);
  if (polygons == (Polygons *) NULL) return 1;

  /* Loop on the samples of the second Db */

  int ntotal = db2->getSampleNumber();
  int nactive = 0;
  int nout = 0;
  int nin = 0;
  for (int iech = 0; iech < ntotal; iech++)
  {
    if (!db2->isActive(iech)) continue;
    nactive++;
    if (!polygon_inside(db2->getCoordinate(iech, 0),
                        db2->getCoordinate(iech, 1),
                        db2->getCoordinate(iech, 2), 0, polygons))
    {
      db2->setArray(iech, isel, 0.);
      nout++;
    }
    else
    {
      nin++;
    }
  }

  // Verbose optional output
  if (verbose)
  {
    mestitle(1,"Convex Hull calculation");
    message("- Number of target samples = %d\n",ntotal);
    message("- Number of active samples = %d\n",nactive);
    message("- Number of masked samples = %d\n",nout);
    message("- Number of valid samples  = %d\n",nin);
  }

  // Set the Naming Convention
  namconv.setNamesAndLocators(db2, isel);
  return 0;
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
static int st_multilinear(Db *db_grid,
                          const VectorInt& indg,
                          int iatt,
                          double *value)
{
  int jech = db_index_grid_to_sample(db_grid, indg.data());
  if (jech < 0) return (1);
  if (!db_grid->isActive(jech)) return (1);
  *value = db_grid->getArray(jech, iatt);
  if (FFFF(*value)) return (1);
  return (0);
}

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
                     const VectorInt& indg1,
                     const VectorDouble& prop,
                     VectorInt& indg2,
                     double *weight)
{
  int idim, ndiv, ival, ndim;
  double wgt;

  wgt  = 1.;
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

/*****************************************************************************/
/*!
 ** Perform the multilinear interpolation from a regular grid Db
 **
 ** \return  Interpolated value
 **
 ** \param[in]  dbgrid    descriptor of the grid parameters
 ** \param[in]  iatt      rank of the target variable in dbgrid
 ** \param[in]  ldmax     Type of distance for calculating maximum distance
 **                       1 for L1 and 2 for L2 distance
 ** \param[in]  dmax      Array of maximum distances (optional)
 ** \param[in]  coor      Coordinates of the target point
 **
 *****************************************************************************/
static double st_multilinear_interpolation(Db *dbgrid,
                                           int iatt,
                                           int ldmax,
                                           const VectorDouble& dmax,
                                           double *coor)
{
  int ndim = dbgrid->getNDim();
  int number = (int) pow(2., ndim);
  VectorInt iwork2(ndim);
  VectorInt indg(ndim);
  VectorDouble aux(ndim);
  VectorDouble prop(ndim);

  /* Identify the closest grid node */

  if (point_to_grid(dbgrid, coor, 0, indg.data()) != 0) return TEST;
  grid_to_point(dbgrid, indg.data(), NULL, aux.data());

  /* Calculate distance to lower corner as proportion of the grid mesh */

  double rtot = 0.;
  for (int idim = 0; idim < ndim; idim++)
  {
    double mesh = dbgrid->getDX(idim);
    double delta = coor[idim] - aux[idim];
    if (delta < 0)
    {
      indg[idim]--;
      delta += mesh;
    }
    if (! dmax.empty() && ldmax == 1)
    {
      if (dmax[idim] <= 0) return TEST;
      double ratio = delta / dmax[idim];
      rtot += ratio * ratio;
    }
    if (! dmax.empty() && delta > dmax[idim]) return TEST;
    prop[idim] = delta / mesh;
  }
  if (! dmax.empty() && ldmax == 1 && rtot > 1.) return TEST;

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
    if (!st_multilinear(dbgrid, iwork2, iatt, &value))
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
 ** \param[in]  ldmax     Type of distance for calculating maximum distance
 **                       1 for L1 and 2 for L2 distance
 ** \param[in]  dmax      Array of maximum distances (optional)
 **
 ** \param[out]  tab      Output array
 **
 ** \remark A point which does not lie between two valuated grid nodes
 ** \remark (in all space dimensions) is always set to FFFF
 **
 *****************************************************************************/
static int st_interpolate_grid_to_point(Db *db_grid,
                                        Db *db_point,
                                        int iatt,
                                        int ldmax,
                                        const VectorDouble& dmax,
                                        VectorDouble& tab)
{
  int iech, error;
  double *coor;

  /* Initializations */

  error = 1;
  coor = (double *) NULL;

  /* Preliminary checks */

  if (! db_grid->hasLargerDimension(db_point)) goto label_end;

  /* Core allocation */

  coor = db_sample_alloc(db_point, LOC_X);
  if (coor == (double *) NULL) goto label_end;

  /* Loop on the point samples */

  for (iech = 0; iech < db_point->getSampleNumber(); iech++)
  {
    if (!db_point->isActive(iech)) continue;
    db_sample_load(db_point, LOC_X, iech, coor);
    tab[iech] = st_multilinear_interpolation(db_grid, iatt, ldmax, dmax, coor);
  }

  /* Set the error return code */

  error = 0;

  /* Core deallocation */

  label_end: coor = db_sample_free(coor);
  return (error);
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
GEOSLIB_API int interpolate_variable_to_point(Db *db_grid,
                                              int iatt,
                                              int np,
                                              double *xp,
                                              double *yp,
                                              double *zp,
                                              double *tab)
{
  int error, ndim;
  double coor[3];

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
  if ((ndim >= 1 && xp == (double *) NULL) || (ndim >= 2
      && yp == (double *) NULL)
      || (ndim >= 3 && zp == (double *) NULL))
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
 **  Generates the discretized points along the trace
 **
 ** \param[in]  nseg   Number of vertices along the trace
 ** \param[in]  trace  Array defining the trace
 **                    (Dimension: 2 * nseg)
 ** \param[in]  disc   Discretization distance
 **
 ** \param[out] np_arg   Number of discretized points
 ** \param[out] xp_arg   Array of first coordinates
 ** \param[out] yp_arg   Array of second coordinates
 ** \param[out] dd_arg   Array of distances between discretized points
 ** \param[out] del_arg  Array of distances between vertices
 ** \param[out] dist_arg Total distance of the trace
 **
 *****************************************************************************/
GEOSLIB_API void ut_trace_discretize(int nseg,
                                     double *trace,
                                     double disc,
                                     int *np_arg,
                                     double **xp_arg,
                                     double **yp_arg,
                                     double **dd_arg,
                                     double **del_arg,
                                     double *dist_arg)
{
  double *xp, *yp, *dd, *del, deltax, deltay, x0, y0, x1, y1, dist;
  int iseg, np, ecr, nloc, ip;

  /* Initializations */

  xp = yp = dd = (double *) NULL;
  (*np_arg) = np = 0;
  (*dist_arg) = x0 = y0 = x1 = y1 = 0.;
  del = (double *) mem_alloc(sizeof(double) * nseg, 1);
  del[0] = 0.;

  /* Loop on the trace segments */

  for (iseg = ecr = 0; iseg < nseg - 1; iseg++)
  {

    /* Consider a segment trace */

    x0 = TRACE(0, iseg);
    y0 = TRACE(1, iseg);
    x1 = TRACE(0, iseg + 1);
    y1 = TRACE(1, iseg + 1);
    deltax = x1 - x0;
    deltay = y1 - y0;
    dist = sqrt(deltax * deltax + deltay * deltay);
    (*dist_arg) += dist;
    del[iseg + 1] = (*dist_arg);

    /* Discretize the trace segment */

    nloc = (int) floor(dist / disc);
    if (ABS(nloc * disc - dist) < dist / 1000) nloc--;
    np += nloc;
    xp = (double *) mem_realloc((char * ) xp, sizeof(double) * np, 1);
    yp = (double *) mem_realloc((char * ) yp, sizeof(double) * np, 1);

    for (ip = 0; ip < nloc; ip++, ecr++)
    {
      xp[ecr] = x0 + deltax * ip / nloc;
      yp[ecr] = y0 + deltay * ip / nloc;
    }
  }

  /* Adding the last vertex */

  np++;
  xp = (double *) mem_realloc((char * ) xp, sizeof(double) * np, 1);
  yp = (double *) mem_realloc((char * ) yp, sizeof(double) * np, 1);
  xp[ecr] = x1;
  yp[ecr] = y1;
  ecr++;

  /* Elaborate the vector of distances */

  dd = (double *) mem_alloc(sizeof(double) * np, 1);
  dd[0] = 0.;
  for (ip = 0; ip < np - 1; ip++)
  {
    deltax = xp[ip + 1] - xp[ip];
    deltay = yp[ip + 1] - yp[ip];
    dd[ip + 1] = dd[ip] + sqrt(deltax * deltax + deltay * deltay);
  }

  /* Returning arguments */

  (*np_arg) = np;
  (*xp_arg) = xp;
  (*yp_arg) = yp;
  (*dd_arg) = dd;
  (*del_arg) = del;
  return;
}

/*****************************************************************************/
/*!
 **  Sample the point Db close to discretized points of the trace
 **
 ** \param[in]  db     Db to be sampled
 ** \param[in]  ptype  Type of locator (::ENUM_LOCS)
 ** \param[in]  np     Number of discretized points
 ** \param[in]  xp     Array of first coordinates
 ** \param[in]  yp     Array of second coordinates
 ** \param[in]  dd     Array of distances
 ** \param[in]  radius Neighborhood radius
 **
 ** \param[out] ns_arg     Number of sampled points
 ** \param[out] xs_arg     Array of first coordinates of sampled points
 ** \param[out] ys_arg     Array of second coordinates of sampled points
 ** \param[out] rks_arg    Array of sample indices (starting from 1)
 ** \param[out] lys_arg    Array of layer indices of sampled points
 ** \param[out] typ_arg    Array of data type
 **                        1 for hard data in Z or TIME
 **                        2 for lower bound
 **                        3 for upper bound
 **
 *****************************************************************************/
GEOSLIB_API void ut_trace_sample(Db *db,
                                 ENUM_LOCS ptype,
                                 int np,
                                 double *xp,
                                 double *yp,
                                 double *dd,
                                 double radius,
                                 int *ns_arg,
                                 double **xs_arg,
                                 double **ys_arg,
                                 int **rks_arg,
                                 int **lys_arg,
                                 int **typ_arg)
{
  int *lys, *typ, *rks, iech, ip, ns, ipmin, nvar;
  double *xs, *ys, cote, layer, bound[2];
  double radcarre, xx, yy, delx, dely, dist, ddmin;

  /* Initializations */

  radcarre = radius * radius;
  xs = ys = (double *) NULL;
  lys = typ = rks = (int *) NULL;
  ns = 0;
  nvar = db->getIntervalNumber();

  /* Loop on the samples */

  for (iech = 0; iech < db->getSampleNumber(); iech++)
  {
    if (!db->isActive(iech)) continue;

    /* Coordinates of the sample point */

    xx = db->getCoordinate(iech, 0);
    yy = db->getCoordinate(iech, 1);

    /* Loop on the discretized samples */

    ipmin = -1;
    ddmin = 1.e30;
    for (ip = 0; ip < np; ip++)
    {
      delx = xx - xp[ip];
      dely = yy - yp[ip];
      dist = delx * delx + dely * dely;
      if (dist > radcarre || dist > ddmin) continue;
      ddmin = dist;
      ipmin = ip;
    }
    if (ipmin < 0) continue;

    /* Keep sample defined by locator */

    cote = get_LOCATOR_ITEM(db, ptype, 0, iech);
    if (!FFFF(cote))
    {
      layer = get_LOCATOR_ITEM(db, LOC_LAYER, 0, iech);
      xs = (double *) mem_realloc((char * ) xs, sizeof(double) * (ns + 1), 1);
      ys = (double *) mem_realloc((char * ) ys, sizeof(double) * (ns + 1), 1);
      lys = (int *) mem_realloc((char * ) lys, sizeof(int) * (ns + 1), 1);
      typ = (int *) mem_realloc((char * ) typ, sizeof(int) * (ns + 1), 1);
      rks = (int *) mem_realloc((char * ) rks, sizeof(int) * (ns + 1), 1);
      xs[ns] = dd[ipmin];
      ys[ns] = cote;
      lys[ns] = (FFFF(layer)) ? 1 :
                                (int) layer + 1;
      typ[ns] = 1;
      rks[ns] = iech + 1;
      ns++;
    }

    /* Keep sample defined by locator UP or LOW */

    for (int ivar = 0; ivar < nvar; ivar++)
    {
      bound[0] = db->getLowerBound(iech, ivar);
      bound[1] = db->getUpperBound(iech, ivar);
      for (int ib = 0; ib < 2; ib++)
      {
        if (FFFF(bound[ib])) continue;
        xs = (double *) mem_realloc((char * )xs, sizeof(double) * (ns + 1), 1);
        ys = (double *) mem_realloc((char * )ys, sizeof(double) * (ns + 1), 1);
        lys = (int *) mem_realloc((char * )lys, sizeof(int) * (ns + 1), 1);
        typ = (int *) mem_realloc((char * )typ, sizeof(int) * (ns + 1), 1);
        rks = (int *) mem_realloc((char * )rks, sizeof(int) * (ns + 1), 1);
        xs[ns] = dd[ipmin];
        ys[ns] = bound[ib];
        lys[ns] = ivar + 1;
        typ[ns] = ib + 2;
        rks[ns] = iech + 1;
        ns++;
      }
    }
  }

  /* Returning arguments */

  *ns_arg = ns;
  *xs_arg = xs;
  *ys_arg = ys;
  *lys_arg = lys;
  *typ_arg = typ;
  *rks_arg = rks;
  return;
}

/*****************************************************************************/
/*!
 **  Derive the external information(s) from the Output db (if Grid)
 **  to the Input Db
 **
 ** \return  Error return code
 **
 ** \param[in]  mode        1 for allocation; -1 for deallocation
 ** \param[in]  locatorType      Type of the pointer (::ENUM_LOCS)
 ** \param[in]  dbin        Descriptor of the input Db
 ** \param[in]  dbout       Descriptor of the output Db
 ** \param[in]  istart      Address of the first allocated external information
 **
 ** \param[out] istart Address of the first external information variable
 **
 ** \remark This function only functions when the Output Db is a grid
 ** \remark However, in case of a Point output Db, this function should not
 ** \remark be used: the external drift functions should already be present
 ** \remark in the output Db.
 ** \remark If this is not the case, an error is issued.
 **
 *****************************************************************************/
GEOSLIB_API int manage_external_info(int mode,
                                     ENUM_LOCS locatorType,
                                     Db *dbin,
                                     Db *dbout,
                                     int *istart)
{
  int info, jstart, iatt, jatt, nechin, ninfo;
  VectorDouble tab;

  if (dbin == (Db *) NULL) return (0);
  nechin = dbin->getSampleNumber();
  ninfo = get_LOCATOR_NITEM(dbout, locatorType);
  if (ninfo <= 0) return (0);

  /* Case when the Output Db is not a grid */

  if (!is_grid(dbout))
  {
    if (get_LOCATOR_NITEM(dbin, locatorType) == ninfo) return (0);
    messerr("The Output Db is not a Grid file");
    messerr("The Input Db does not contain the %d External Drifts");
    return (1);
  }

  /* Dispatch */

  if (mode > 0)
  {

    /* Core allocation */

    *istart = -1;
    tab.resize(nechin, 0.);

    /* Allocation */

    for (info = 0; info < ninfo; info++)
    {

      /* If the drift vector is present in the input file, skip the rest */

      if (dbin->getLocatorNumber(locatorType) > 0)
      {
        jatt = db_attribute_identify(dbin, locatorType, info);
        if (jatt >= 0) continue;
      }

      /* Add the Drift vector in the Input file */

      jatt = dbin->addFields(1, 0.);
      if (jatt < 0) return (1);
      if (*istart < 0) *istart = jatt;

      /* Locate the variable in the Grid file */

      iatt = db_attribute_identify(dbout, locatorType, info);
      if (iatt < 0) return (1);
      dbin->setLocatorByAttribute(jatt, locatorType, info);

      /* Perform the migration */

      if (st_migrate_grid_to_point(dbout, dbin, iatt, 0, VectorDouble(), tab))
        continue;

      /* Save the migrated array */

      dbin->setFieldByAttribute(tab, jatt);
    }
  }
  else
  {
    jstart = *istart;
    if (jstart < 0) return (0);
    for (info = 0; info < ninfo; info++)
    {
      jatt = db_attribute_identify(dbin, locatorType, info);
      if (jatt >= *istart) dbin->deleteFieldByAttribute(jatt);
    }
  }
  return (0);
}

/*****************************************************************************/
/*!
 **  Derive the non-stationary information(s) from the Output db (if Grid)
 **  to the Input Db
 **
 ** \return  Error return code
 **
 ** \param[in]  mode        1 for allocation; -1 for deallocation
 ** \param[in]  model       Descriptor of the Model
 ** \param[in]  dbin        Descriptor of the input Db
 ** \param[in]  dbout       Descriptor of the output Db
 **
 *****************************************************************************/
GEOSLIB_API int manage_nostat_info(int mode, Model* model, Db *dbin, Db *dbout)
{
  VectorDouble tab;

  if (! model->isNoStat()) return 0;
  const ANoStat* nostat = model->getNoStat();

  /* Dispatch */

  if (mode > 0)
  {

    // Attach the Input Db
    if (nostat->attachToDb(dbin,1)) return 1;

    // Attach the Output Db
    if (nostat->attachToDb(dbout,2)) return 1;
  }
  else
  {

    // Detach the Input Db
    nostat->detachFromDb(dbin,1);

    // Detach the output Db
    nostat->detachFromDb(dbout,2);
  }
  return (0);
}

/*****************************************************************************/
/*!
 **  Migrates the samples of a Db to the center of blocks of a grid Db
 **
 ** \return  Error return code
 **
 ** \param[in]  db_point   descriptor of the point parameters
 ** \param[in]  db_grid    descriptor of the grid parameters
 ** \param[in]  eps_random Randomisation Epsilon
 **
 ** \remark The argument 'eps_random' allows perturbating the centered
 ** \remark coordinate so that it does not lie exactly on the node.
 ** \remark This possibility makes sense in order to identify migrated data
 ** \remark from data actually located on the grid center (before migration)
 ** \remark The perturbation is calculated as DX(i) * eps
 **
 *****************************************************************************/
GEOSLIB_API int db_center_point_to_grid(Db *db_point,
                                        Db *db_grid,
                                        double eps_random)
{
  int *indg, iech, error, idim, ndim;
  double *coor;

  /* Initializations */

  error = 1;
  indg = (int *) NULL;
  coor = (double *) NULL;

  /* Preliminary checks */

  if (db_point == (Db *) NULL) goto label_end;
  if (db_grid == (Db *) NULL) goto label_end;
  if (!is_grid(db_grid)) goto label_end;
  if (! db_point->hasSameDimension(db_grid)) goto label_end;
  ndim = db_point->getNDim();

  /* Core allocation */

  coor = db_sample_alloc(db_point, LOC_X);
  if (coor == (double *) NULL) goto label_end;
  indg = db_indg_alloc(db_grid);
  if (indg == (int *) NULL) goto label_end;

  /* Loop on the samples of the Point Db */

  for (iech = 0; iech < db_point->getSampleNumber(); iech++)
  {

    /* Read the coordinates of the point sample */

    for (idim = 0; idim < ndim; idim++)
      coor[idim] = db_point->getCoordinate(iech, idim);

    /* Get the indices of the grid node */

    (void) point_to_grid(db_grid, coor, -1, indg);

    /* Get the coordinates of the grid center */

    grid_to_point(db_grid, indg, NULL, coor);

    /* Randomize the migrated center */

    if (eps_random > 0) for (idim = 0; idim < ndim; idim++)
      coor[idim] += db_grid->getDX(idim) * law_uniform(0., eps_random);

    /* Correct the sample locations */

    for (idim = 0; idim < ndim; idim++)
      db_point->setCoordinate(iech, idim, coor[idim]);
  }

  /* Set the error return code */

  error = 0;

  label_end: coor = db_sample_free(coor);
  indg = db_indg_free(indg);
  return (error);
}

/*****************************************************************************/
/*!
 **  Sample a grid into a finer subgrid (all variables)
 **
 ** \return  Error return code
 **
 ** \param[in]  dbin  Descriptor of the grid parameters
 ** \param[in]  nmult Array of multiplicity coefficients
 **
 *****************************************************************************/
GEOSLIB_API Db *db_grid_sample(Db *dbin, const VectorInt& nmult)
{
  Db *dbout;
  VectorDouble coor;
  int ncol, icol, iech, iad, item, rank, ndim;
  ENUM_LOCS locatorType;

  /* Initializations */

  dbout = (Db *) NULL;
  ncol = dbin->getFieldNumber();
  ndim = dbin->getNDim();

  /* Core allocation */

  coor.resize(ndim);

  /* Create the subgrid */

  dbout = db_create_grid_multiple(dbin, nmult, 1);
  if (dbout == (Db *) NULL) goto label_end;
  rank = dbout->addFields(ncol, TEST);
  if (rank < 0) goto label_end;
  for (icol = 0; icol < ncol; icol++)
  {
    (void) dbin->getLocatorByColumn(icol, &locatorType, &item);
    dbout->setLocatorByAttribute(icol, locatorType, item);
  }

  /* Loop on the samples of the output grid */

  for (iech = 0; iech < dbout->getSampleNumber(); iech++)
  {
    if (!dbout->isActive(iech)) continue;
    db_sample_load(dbout, LOC_X, iech, coor.data());
    iad = dbin->coordinateToRank(coor);
    if (iad < 0) continue;

    /* Loop on the variables of the input grid */

    for (icol = 0; icol < ncol; icol++)
      dbout->setByColumn(iech, icol, dbin->getByColumn(iad, icol));
  }

  label_end:
  return (dbout);
}

/*****************************************************************************/
/*!
 **  Check distance in 1D from Data point to the target
 **
 ** \returns The Distance
 **
 ** \param[in]  dbpoint     Descriptor of the point parameters
 ** \param[in]  ip          Rank of the sample in the Point file
 ** \param[in]  idim_ref    Space dimension where the quick test is performed
 ** \param[in]  xtarget     Coordinate of the target along X
 **
 *****************************************************************************/
static double st_get_1d_distance(Db *dbpoint, int ip, int idim_ref, double xtarget)
{
  return ABS(dbpoint->getCoordinate(ip, idim_ref) - xtarget);
}

/*****************************************************************************/
/*!
 **  Update minimum distance and rank of the corresponding sample
 **
 ** \return - 1 If the distance has not been minimized
 ** \return - otherwise
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
 ** \param[out]     dvect   Vector for distance increments (Dimension: ndim)
 ** \param[out]     dvmin   Vector for minimum distance increment (Dim: ndim)
 **
 ** \remarks The Time Shift is an optional variable which increases the
 ** \remarks distance (the time-to-distance conversion is assumed to be 1)
 ** \remarks Only positive Time Shifts are considered
 **
 *****************************************************************************/
static int st_get_closest_sample(Db *dbgrid,
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
                                 VectorDouble& dvect,
                                 VectorDouble& dvmin)
{
  double dd, time, angle, scaleu, scalev, scalew, x, y, dloc, cosa, sina, dmem0, dmem1;
  int ndim;

  // Default values
  ndim = dbgrid->getNDim();
  angle = 0.;
  scaleu = 1.;
  scalev = 1.;
  scalew = 1.;

  // Calculate the euclidean distance
  dd = distance_inter(dbgrid, dbpoint, ig, ip, dvect.data());
  dmem0 = dvect[0];
  dmem1 = dvect[1];

  // Case of anisotropic distances
  if (flag_aniso)
  {

    // Rotation angle in degrees (optional)
    if (iatt_angle >= 0 && ndim >= 2)
    {
      angle = dbgrid->getArray(ig, iatt_angle);
      ut_rotation_sincos(angle, &cosa, &sina);
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

  // Set the initial values in the first two elements of 'dvect'

  dvect[0] = dmem0;
  dvect[1] = dmem1;

  // Add the Time Shift penalty (optional)
  if (iatt_time >= 0)
  {
    time = dbpoint->getArray(ip, iatt_time);
    if (time > 0) dd += time;
  }

  // Evaluate the closest distance
  if (dd >= (*ddmin)) return (1);
  (*ddmin) = dd;
  (*ipmin) = ip;
  for (int idim = 0; idim < ndim; idim++)
    dvmin[idim] = dvect[idim];
  return (0);
}

/*****************************************************************************/
/*!
 **  Migrate a variable from the point structure
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
 ** \param[in]  ldmax       Type of distance for calculating maximum distance
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
GEOSLIB_API int expand_point_to_grid(Db *db_point,
                                     Db *db_grid,
                                     int iatt,
                                     int iatt_time,
                                     int iatt_angle,
                                     int iatt_scaleu,
                                     int iatt_scalev,
                                     int iatt_scalew,
                                     int flag_index,
                                     int ldmax,
                                     const VectorDouble& dmax,
                                     VectorDouble& tab)
{
  /* Preliminary checks */

  if (! db_point->hasLargerDimension(db_grid)) return 1;
  int ndim_min = MIN(db_point->getNDim(),db_grid->getNDim());
  int ndim_max = MAX(db_point->getNDim(),db_grid->getNDim());
  bool flag_aniso = 0;
  if (ndim_min >= 1 && iatt_scaleu >= 0) flag_aniso = 1;
  if (ndim_min >= 2 && iatt_scalev >= 0) flag_aniso = 1;
  if (ndim_min >= 3 && iatt_scalew >= 0) flag_aniso = 1;
  int idim_ref = ndim_min - 1;

  // Core allocation

  int ng = db_grid->getSampleNumber();
  int np = db_point->getSampleNumber();
  VectorDouble xtab(np);
  VectorDouble dvect(ndim_max);
  VectorDouble dvmin(ndim_max);
  VectorInt    rank(np);

  /* Sort the point samples according to their coordinate ranked 'idim_ref' */

  for (int ip = np = 0; ip < db_point->getSampleNumber(); ip++)
  {
    if (!db_point->isActive(ip)) continue;
    if (FFFF(db_point->getArray(ip, iatt))) continue;
    rank[np] = ip;
    xtab[np] = db_point->getCoordinate(ip, idim_ref);
    np++;
  }
  ut_sort_double(1, np, rank.data(), xtab.data());

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

  int npin = 0;
  for (int ig = 0; ig < ng; ig++)
  {
    if (! db_grid->isActive(ig)) continue;
    double xtarget = db_grid->getCoordinate(ig, idim_ref);

    /* Locate the grid node within the ordered list (1D coordinate) */

    int ip0 = 0;
    for (int ip = 0; ip < np; ip++)
    {
      int jp = rank[ip];
      if (xtarget < xtab[jp]) break;
      ip0 = ip;
    }

    /* Calculate minimum distance between the two closest ordered samples */

    int ipmin = -1;
    double ddmin = 1.e30;
    if (ip0 >= 0)
      st_get_closest_sample(db_grid, ig, db_point, rank[ip0], flag_aniso,
                            iatt_time, iatt_angle, iatt_scaleu, iatt_scalev,
                            iatt_scalew, &ipmin, &ddmin, dvect, dvmin);
    if (ip0 + 1 < np - 1)
      st_get_closest_sample(db_grid, ig, db_point, rank[ip0 + 1], flag_aniso,
                            iatt_time, iatt_angle, iatt_scaleu, iatt_scalev,
                            iatt_scalew, &ipmin, &ddmin, dvect, dvmin);

    if (! dmax.empty() && dvmin[idim_ref] > dmax[idim_ref]) continue;

    /* Look for closer points for samples located below rank[ip0] */

    for (int ip = ip0 - 1; ip >= 0; ip--)
    {
      int jp = rank[ip];
      double dd1d = st_get_1d_distance(db_point, jp, idim_ref, xtarget);
      if (! flag_aniso && dd1d > ddmin + time_max) break;
      if (! dmax.empty() && dd1d > dmax[idim_ref]) break;
      if (st_get_closest_sample(db_grid, ig, db_point, jp, flag_aniso,
                                iatt_time, iatt_angle, iatt_scaleu, iatt_scalev,
                                iatt_scalew, &ipmin, &ddmin, dvect, dvmin))
        continue;
    }

    /* Look for closer points for samples located above rank[ip0+1] */

    for (int ip = ip0 + 1; ip < np; ip++)
    {
      int jp = rank[ip];
      double dd1d = st_get_1d_distance(db_point, jp, idim_ref, xtarget);
      if (! flag_aniso && dd1d > ddmin + time_max) break;
      if (! dmax.empty() && dd1d > dmax[idim_ref]) break;
      if (st_get_closest_sample(db_grid, ig, db_point, jp, flag_aniso,
                                iatt_time, iatt_angle, iatt_scaleu, iatt_scalev,
                                iatt_scalew, &ipmin, &ddmin, dvect, dvmin))
        continue;
    }

    /* Truncation by 'dmax' if provided */

    if (! dmax.empty())
    {
      (void) distance_inter(db_grid, db_point, ig, ipmin, dvmin.data());
      if (st_larger_than_dmax(ndim_min, dvmin, ldmax, dmax)) continue;
    }

    /* Set the value */

    npin++;
    if (flag_index)
      tab[ig] = (ipmin < 0) ? TEST : (double) ipmin;
    else
      tab[ig] = (ipmin < 0) ? TEST : db_point->getArray(ipmin, iatt);
  }
  return 0;
}

/*****************************************************************************/
/*!
 **  Write a sample in the output file
 **
 ** \param[in]  db        Db characteristics
 ** \param[in]  iech      Rank of the sample
 ** \param[in]  nvar      Number of attributes
 ** \param[in]  iatt      Array of the output attribute
 ** \param[in]  tab       Array of values
 **
 *****************************************************************************/
static void st_write_active_sample(Db *db,
                                   int iech,
                                   int nvar,
                                   int *iatt,
                                   double *tab)
{
  int ivar;

  for (ivar = 0; ivar < nvar; ivar++)
    db->setArray(iech, iatt[ivar], tab[ivar]);
}

/*****************************************************************************/
/*!
 **  Check if a sample must be kept or not
 **
 ** \return  1 if the sample must be kept
 **
 ** \param[in]  db        Db characteristics
 ** \param[in]  flag_zero 1 to forbid zero values
 ** \param[in]  iech      Rank of the sample
 ** \param[in]  nvar      Number of attributes
 ** \param[in]  iatt      Array of the input attribute
 ** \param[in]  eps       Minimum value assigned to a zero
 **
 ** \param[out] tab       Array of values
 **
 *****************************************************************************/
static int st_read_active_sample(Db *db,
                                 int flag_zero,
                                 int iech,
                                 int nvar,
                                 int *iatt,
                                 double eps,
                                 double *tab)
{
  int ivar, number;

  /* Initialize the output array */

  for (ivar = 0; ivar < nvar; ivar++)
    tab[ivar] = TEST;

  /* Check for the presence of a Test value */

  for (ivar = 0; ivar < nvar; ivar++)
  {
    tab[ivar] = db->getArray(iech, iatt[ivar]);
    if (FFFF(tab[ivar])) return (0);
  }

  /* Count the number of zeroes */

  if (!flag_zero) return (1);
  number = 0;
  for (ivar = 0; ivar < nvar; ivar++)
    if (tab[ivar] <= 0.) number++;

  /* If zeroes have been found, correct values */

  if (number > 0)
  {
    for (ivar = 0; ivar < nvar; ivar++)
    {
      if (tab[ivar] <= 0.)
        tab[ivar] = eps;
      else
        tab[ivar] = tab[ivar] - number * eps / (nvar - number);
    }
  }
  return (1);
}

/*****************************************************************************/
/*!
 **  Translate a set of compositional variables into auxiliary variables
 **
 ** \return  Error return code
 **
 ** \param[in]  db       Db characteristics
 ** \param[in]  verbose  1 for a Verbose option
 ** \param[in]  mode     1 for Forward and -1 for Backward transformation
 ** \param[in]  type     Type of conversion
 ** \li                  0 : Simple transformation in proportion
 ** \li                  1 : Additive logratio
 ** \li                  2 : Centered logratio
 ** \li                  3 : Isometric logratio
 ** \param[in]  number   Number of input attributes
 ** \param[in]  iatt_in  Array of the input attribute
 ** \param[in]  iatt_out Array of the output attribute
 **
 ** \param[out] numout   Number of variables in output
 **
 ** \remarks  The additive and the isometric logratio transformations
 ** \remarks  transform N+1 compositional variables into N elements (Forward)
 ** \remarks  and from N transformed elements into N+1 compositional variables
 ** \remarks  (Backwards).
 **
 ** \remarks  The zero-values are replaced by a conventional small value
 ** \remarks  This is defined by a variable that can be corrected using
 ** \remarks  the keypair facility with the keyword 'CompositionalEps'
 **
 ** \remarks  The arguments 'iatt_in' and 'iatt_out' can coincide.
 **
 ** \remarks  Outlier Detection for Compositional Data Using Robust Methods
 ** \remarks  Math Geosciences (2008) 40: 233-248
 **
 *****************************************************************************/
GEOSLIB_API int db_compositional_transform(Db *db,
                                           int verbose,
                                           int mode,
                                           int type,
                                           int number,
                                           int *iatt_in,
                                           int *iatt_out,
                                           int *numout)
{
  int error, nech, number1, iech, ivar, jvar;
  double *tabin, *tabout, sum, eps;

  /* Initializations */

  error = 1;
  nech = db->getSampleNumber();
  tabin = tabout = (double *) NULL;
  eps = get_keypone("CompositionalEps", EPSILON3);

  /* Core allocation (may be one more than needed, but general) */

  number1 = number + 1;
  tabin = (double *) mem_alloc(sizeof(double) * number1, 0);
  if (tabin == (double *) NULL) goto label_end;
  tabout = (double *) mem_alloc(sizeof(double) * number1, 0);
  if (tabout == (double *) NULL) goto label_end;

  /* Verbose output */

  if (verbose) mestitle(0, "Compositional transformation");

  /* Dispatch */

  switch (type)
  {
    case 0: /* No action */
      if (verbose)
      {
        if (mode > 0)
          message("- Scaling to Proportions (Forwards)\n");
        else
          message("- Scaling to Proportions (Backwards)\n");
      }
      (*numout) = number;
      for (iech = 0; iech < nech; iech++)
      {
        if (!db->isActive(iech)) continue;
        (void) st_read_active_sample(db, 1, iech, number, iatt_in, eps, tabin);
        sum = 0.;
        for (ivar = 0; ivar < (*numout); ivar++)
          sum += tabin[ivar];
        for (ivar = 0; ivar < (*numout); ivar++)
          tabout[ivar] = (sum > 0) ? tabin[ivar] / sum :
                                     TEST;
        st_write_active_sample(db, iech, *numout, iatt_out, tabout);
      }
      break;

    case 1: /* Additive Logratio */
      if (mode > 0)
      {
        /* Forward transformation */

        if (verbose) message("- using Additive Log-Ratio (Forwards)\n");
        (*numout) = number - 1;
        for (iech = 0; iech < nech; iech++)
        {
          if (!db->isActive(iech)) continue;
          if (st_read_active_sample(db, 1, iech, number, iatt_in, eps, tabin))
          {
            for (ivar = 0; ivar < (*numout); ivar++)
              tabout[ivar] = log(tabin[ivar] / tabin[(*numout)]);
          }
          st_write_active_sample(db, iech, *numout, iatt_out, tabout);
        }
      }
      else
      {
        /* Backward transformation */

        message("- using Additive Log-Ratio (Backwards)\n");
        (*numout) = number + 1;
        for (iech = 0; iech < nech; iech++)
        {
          if (!db->isActive(iech)) continue;
          if (st_read_active_sample(db, 0, iech, number, iatt_in, eps, tabin))
          {
            sum = 1.;
            for (ivar = 0; ivar < number; ivar++)
            {
              tabout[ivar] = exp(tabin[ivar]);
              sum += tabout[ivar];
            }
            for (ivar = 0; ivar < number; ivar++)
              tabout[ivar] /= sum;
            tabout[number] = 1 / sum;
          }
          st_write_active_sample(db, iech, *numout, iatt_out, tabout);
        }
      }
      break;

    case 2: /* Centered Logratio */
      if (mode > 0)
      {
        /* Forward transformation */

        if (verbose) message("- using Centered Log-Ratio (Forwards)\n");
        (*numout) = number;
        for (iech = 0; iech < nech; iech++)
        {
          if (!db->isActive(iech)) continue;
          if (st_read_active_sample(db, 1, iech, number, iatt_in, eps, tabin))
          {
            sum = 0.;
            for (ivar = 0; ivar < number; ivar++)
            {
              tabout[ivar] = log(tabin[ivar]);
              sum += tabout[ivar];
            }
            for (ivar = 0; ivar < number; ivar++)
              tabout[ivar] = tabout[ivar] - 0.5 * sum;
          }
          st_write_active_sample(db, iech, *numout, iatt_out, tabout);
        }
      }
      else
      {
        /* Backward transformation */

        if (verbose) message("- using Centered Log-Ratio (Backwards)\n");
        (*numout) = number;
        for (iech = 0; iech < nech; iech++)
        {
          if (!db->isActive(iech)) continue;
          if (st_read_active_sample(db, 0, iech, number, iatt_in, eps, tabin))
          {
            sum = 0.;
            for (ivar = 0; ivar < number; ivar++)
            {
              tabout[ivar] = exp(tabin[ivar]);
              sum += tabout[ivar];
            }
            for (ivar = 0; ivar < number; ivar++)
              tabout[ivar] /= sum;
          }
          st_write_active_sample(db, iech, *numout, iatt_out, tabout);
        }
      }
      break;

    case 3: /* Isometric Logratio */
      if (mode > 0)
      {
        /* Forward transformation */

        if (verbose) message("- using Isometric Log-Ratio (Forwards)\n");
        (*numout) = number - 1;
        for (iech = 0; iech < nech; iech++)
        {
          if (!db->isActive(iech)) continue;
          if (st_read_active_sample(db, 1, iech, number, iatt_in, eps, tabin))
          {
            for (ivar = 0; ivar < number; ivar++)
              tabin[ivar] = log(tabin[ivar]);
            sum = 0.;
            for (ivar = 0; ivar < (*numout); ivar++)
            {
              sum += tabin[ivar];
              tabout[ivar] = (sum / (ivar + 1.) - tabin[ivar + 1])
                  * sqrt((ivar + 1.) / (ivar + 2.));
            }
          }
          st_write_active_sample(db, iech, *numout, iatt_out, tabout);
        }
      }
      else
      {
        /* Backward transformation */

        if (verbose) message("- using Isometric Log-Ratio (Backwards)\n");
        (*numout) = number + 1;
        for (iech = 0; iech < nech; iech++)
        {
          if (!db->isActive(iech)) continue;
          if (st_read_active_sample(db, 0, iech, number, iatt_in, eps, tabin))
          {
            tabin[number] = 0.;
            for (ivar = 0; ivar < (*numout); ivar++)
            {
              sum = 0.;
              for (jvar = ivar; jvar < (*numout); jvar++)
                sum += tabin[jvar] / sqrt((jvar + 1.) * (jvar + 2.));
              tabout[ivar] = sum;
              if (ivar > 0)
                tabout[ivar] -= tabin[ivar - 1] * sqrt(ivar / (ivar + 1.));
            }

            sum = 0.;
            for (jvar = 0; jvar < (*numout); jvar++)
            {
              tabout[jvar] = exp(tabout[jvar]);
              sum += tabout[jvar];
            }
            for (ivar = 0; ivar < (*numout); ivar++)
              tabout[ivar] /= sum;
          }
          st_write_active_sample(db, iech, *numout, iatt_out, tabout);
        }
      }
      break;
  }

  /* Set the error return code */

  error = 0;

  label_end: tabin = (double *) mem_free((char * ) tabin);
  tabout = (double *) mem_free((char * ) tabout);
  return (error);
}

/*****************************************************************************/
/*!
 **  Unfold a 2-D Db with respect to a polyline
 **
 ** \return  Error return code
 **
 ** \param[in]  db      Db structure
 ** \param[in]  nvert   Number of segments in the polyline
 ** \param[in]  xl      Array of X-coordinates of the polyline
 ** \param[in]  yl      Array of Y-coordinates of the polyline
 **
 *****************************************************************************/
GEOSLIB_API int db_unfold_polyline(Db *db, int nvert, double *xl, double *yl)
{
  PL_Dist *pldist, *pldist0;
  double xx, yy, newx, newy;
  int error, iptr;

  /* Initializations */

  error = 1;
  pldist = pldist0 = (PL_Dist *) NULL;

  /* Preliminary checks */

  if (db->getNDim() != 2)
  {
    messerr("This function is restricted to 2-D Db");
    goto label_end;
  }
  if (nvert <= 1)
  {
    messerr("This function requires a Polyline with at least one segment");
    goto label_end;
  }

  /* Add the variables */

  iptr = db->addFields(2, 0.);
  if (iptr < 0) goto label_end;

  /* Define the internal structures */

  pldist0 = pldist_manage(1, NULL, 2, nvert);
  pldist = pldist_manage(1, NULL, 2, nvert);

  /* Project the starting point */

  distance_point_to_polyline(xl[0], yl[0], nvert, xl, yl, pldist0);

  /* Loop on the samples of the Db */

  for (int iech = 0; iech < db->getSampleNumber(); iech++)
  {
    if (!db->isActive(iech)) continue;
    xx = db->getCoordinate(iech, 0);
    yy = db->getCoordinate(iech, 1);
    distance_point_to_polyline(xx, yy, nvert, xl, yl, pldist);
    newx = pldist->dist;
    newy = distance_along_polyline(pldist0, pldist, xl, yl);
    db->setArray(iech, iptr, newx);
    db->setArray(iech, iptr + 1, newy);
  }

  /* Set the error return code */

  error = 0;

  label_end: pldist0 = pldist_manage(-1, pldist0, 2, nvert);
  pldist = pldist_manage(-1, pldist, 2, nvert);
  return (error);
}

/*****************************************************************************/
/*!
 **  Fold an input Db into an output Db with respect to a polyline
 **
 ** \return  Error return code
 **
 ** \param[in]  dbin    Input Db structure
 ** \param[in]  dbout   Output Db structure
 ** \param[in]  ncol    Number of target variables
 ** \param[in]  cols    Ranks of the target variables
 ** \param[in]  nvert   Number of segments in the polyline
 ** \param[in]  xl      Array of X-coordinates of the polyline
 ** \param[in]  yl      Array of Y-coordinates of the polyline

 *****************************************************************************/
GEOSLIB_API int db_fold_polyline(Db *dbin,
                                 Db *dbout,
                                 int ncol,
                                 int *cols,
                                 int nvert,
                                 double *xl,
                                 double *yl)
{
  PL_Dist *pldist, *pldist0;
  double xx, yy, value;
  int error, iptr, iad;
  VectorDouble coor(2);

  /* Initializations */

  error = 1;
  pldist = pldist0 = (PL_Dist *) NULL;

  /* Preliminary checks */

  if (dbin->getNDim() != 2 || !is_grid(dbin))
  {
    messerr("This function is restricted to 2-D Input Grid Db");
    goto label_end;
  }
  if (dbout->getNDim() != 2)
  {
    messerr("This function is restricted to 2-D Output Db");
    goto label_end;
  }
  if (nvert <= 1)
  {
    messerr("This function requires a Polyline with at least one segment");
    goto label_end;
  }

  /* Add the variables */

  iptr = dbout->addFields(ncol, TEST);
  if (iptr < 0) goto label_end;

  /* Define the internal structures */

  pldist0 = pldist_manage(1, NULL, 2, nvert);
  pldist = pldist_manage(1, NULL, 2, nvert);

  /* Project the starting point */

  distance_point_to_polyline(xl[0], yl[0], nvert, xl, yl, pldist0);

  /* Loop on the samples of the output Db */

  for (int iech = 0; iech < dbout->getSampleNumber(); iech++)
  {
    if (!dbout->isActive(iech)) continue;
    xx = dbout->getCoordinate(iech, 0);
    yy = dbout->getCoordinate(iech, 1);

    /* Project the target point according to the line */

    distance_point_to_polyline(xx, yy, nvert, xl, yl, pldist);
    coor[0] = pldist->dist;
    coor[1] = distance_along_polyline(pldist0, pldist, xl, yl);

    /* Locate the sample on the Input Grid */

    iad = dbin->coordinateToRank(coor);
    if (iad < 0) continue;

    /* Loop on the variables */

    for (int icol = 0; icol < ncol; icol++)
    {
      value = dbin->getArray(iad, cols[icol]);
      dbout->setArray(iech, iptr + icol, value);
    }
  }

  /* Set the error return code */

  error = 0;

  label_end: pldist0 = pldist_manage(-1, pldist0, 2, nvert);
  pldist = pldist_manage(-1, pldist, 2, nvert);
  return (error);
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
                      Db *dbgrid,
                      VectorDouble& tab1,
                      int *indg0,
                      int *indg,
                      VectorDouble& tab2)
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
      int radius = (int) tab1[iech];
      db_index_sample_to_grid(dbgrid, iech, indg0);

      for (int idim = 0; idim < ndim; idim++)
        for (int ifois = -1; ifois <= 1; ifois += 2)
          for (int irad = 0; irad < radius; irad++)
          {
            for (int jdim = 0; jdim < ndim; jdim++)
              indg[jdim] = indg0[jdim];
            indg[idim] = indg0[idim] + irad * ifois;
            int jech = db_index_grid_to_sample(dbgrid, indg);
            if (jech >= 0) tab2[jech] = (flag_size) ? radius : 1.;
          }
    }
  }
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
GEOSLIB_API int points_to_block(Db *dbpoint,
                                Db *dbgrid,
                                int option,
                                int flag_size,
                                int iatt_time,
                                int iatt_size,
                                int iatt_angle,
                                int iatt_scaleu,
                                int iatt_scalev,
                                int iatt_scalew)
{
  int *indg, *indg0, iatt_edge, iatt_rank, iatt_surf, iatt_vol, iatt_code;
  int val_iech, val_jech, jech, ndim, nvois, lec, error, flag_index;
  VectorInt indret;
  VectorDouble tab1, tab2;

  /* Initializations */

  error = 1;
  indg = indg0 = (int *) NULL;
  iatt_edge = iatt_rank = iatt_surf = iatt_vol = iatt_code = -1;
  if (! dbgrid->hasSameDimension(dbpoint)) goto label_end;
  ndim = dbgrid->getNDim();
  flag_index = get_keypone("PTB_Flag_Index", 0.);

  /* Core allocation */

  indg0 = db_indg_alloc(dbgrid);
  if (indg0 == (int *) NULL) goto label_end;
  indg = db_indg_alloc(dbgrid);
  if (indg == (int *) NULL) goto label_end;
  tab1.resize(dbgrid->getSampleNumber());
  tab2.resize(dbgrid->getSampleNumber(),-1.);

  /* Variable allocation */

  iatt_edge = dbgrid->addFields(1, 0.);
  if (iatt_edge < 0) goto label_end;
  iatt_rank = dbpoint->addFields(1, 0.);
  if (iatt_rank < 0) goto label_end;
  iatt_surf = dbpoint->addFields(1, 0.);
  if (iatt_surf < 0) goto label_end;
  iatt_vol = dbpoint->addFields(1, 0.);
  if (iatt_vol < 0) goto label_end;
  iatt_code = dbpoint->addFields(1, 1.);
  if (iatt_code < 0) goto label_end;

  /* Create the sample rank attribute and expand it over the grid */

  for (int iech = 0; iech < dbpoint->getSampleNumber(); iech++)
    dbpoint->setArray(iech, iatt_rank, (double) iech);
  if (expand_point_to_grid(dbpoint, dbgrid, iatt_rank, iatt_time, iatt_angle,
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

  indret = gridcell_neigh(ndim, option, 1, 1, 0, &nvois);

  /* Loop on the grid nodes */

  for (int iech = 0; iech < dbgrid->getSampleNumber(); iech++)
  {
    if (!dbgrid->isActive(iech)) continue;

    /* Identify the sample within the grid */

    val_iech = (int) tab1[iech];
    db_index_sample_to_grid(dbgrid, iech, indg0);

    /* Increment the volume by one */

    dbpoint->updArray(val_iech, iatt_vol, 0, 1.);

    /* Loop on the neighborhood nodes */

    lec = 0;
    for (int ivois = 0; ivois < nvois; ivois++)
    {
      for (int idim = 0; idim < ndim; idim++)
        indg[idim] = indg0[idim] + indret[lec++];
      jech = db_index_grid_to_sample(dbgrid, indg);
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

        dbpoint->updArray(val_iech, iatt_surf, 0, 1.);

        /* Set the edge */

        tab2[iech] = (iatt_size >= 0) ? dbpoint->getArray(val_iech, iatt_size) : 1.;
      }
    }
  }

  /* If size is defined, expand the perimeter */

  if (iatt_size >= 0)
  {
    for (int i = 0; i < dbgrid->getSampleNumber(); i++) tab1[i] = tab2[i];
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

  label_suite: dbgrid->setFieldByAttribute(tab2, iatt_edge);

  /* Set the error return code */

  error = 0;

  /* Core deallocation */

  label_end:
  if (iatt_rank >= 0) dbpoint->deleteFieldByAttribute(iatt_rank);
  indg = db_indg_free(indg);
  indg0 = db_indg_free(indg0);
  return (error);
}

/*****************************************************************************/
/*!
 **  Get the number of Poisson points to fill a Volume
 **
 ** \return  Number of points
 **
 ** \param[in]  verbose     Verbose option
 ** \param[in]  ndim        Space dimension
 ** \param[in]  density     Average Poisson intensity
 ** \param[in]  extend      Vector of field extends
 **
 ** \remarks  In the case where the number of points generated is zero,
 ** \remarks  a message is issued
 **
 *****************************************************************************/
static int st_get_number_poisson(int verbose,
                                 int ndim,
                                 double *extend,
                                 double density)
{
  double volume;
  int number;

  /* Calculate the volume */

  volume = 1.;
  for (int idim = 0; idim < ndim; idim++)
    volume *= extend[idim];

  /* Draw the count of samples */

  number = law_poisson(density * volume);
  if (number <= 0)
  {
    messerr("For density (%lf) and Volume (%lf) no sample is generated",
            density, volume);
    return (0);
  }

  /* Optional printout */

  if (verbose)
  {
    mestitle(1, "Poisson Point Process");
    message("Average density  = %lf\n", density);
    message("Volume           = %lf\n", volume);
    message("Number of points = %d\n", number);
  }

  return (number);
}

/*****************************************************************************/
/*!
 **  Create a set of samples according to a Poisson process
 **
 ** \return  Array of returned values (or NULL)
 **
 ** \param[in]  verbose     Verbose option
 ** \param[in]  ndim        Space dimension
 ** \param[in]  seed        Seed for the random number generator
 ** \param[in]  density     Average Poisson intensity
 ** \param[in]  origin      Vector of field origin
 ** \param[in]  extend      Vector of field extends
 **
 ** \param[out] count       Number of points generated
 **
 *****************************************************************************/
static VectorDouble st_point_init_poisson(int verbose,
                                          int ndim,
                                          int seed,
                                          double density,
                                          VectorDouble origin,
                                          VectorDouble extend,
                                          int *count)
{
  VectorDouble tab;
  int number;

  /* Initializations */

  law_set_random_seed(seed);

  /* Calculate the number of samples to be generated */

  number = st_get_number_poisson(verbose, ndim, extend.data(), density);
  if (number > 0)
  {
    tab.resize(ndim * number);
    for (int iech = 0; iech < number; iech++)
      for (int idim = 0; idim < ndim; idim++)
        TAB(iech,idim) = origin[idim] + law_uniform(0., extend[idim]);
  }
  *count = number;

  return (tab);
}

/*****************************************************************************/
/*!
 **  Create a set of samples according to a Poisson Regionalized process
 **
 ** \return  Array of returned values (or NULL)
 **
 ** \param[in]  verbose     Verbose option
 ** \param[in]  seed        Seed for the random number generator
 ** \param[in]  dbgrid      Descriptor of the Db grid parameters
 **
 ** \param[out] count       Number of points generated
 **
 *****************************************************************************/
static VectorDouble st_point_init_poisreg(int verbose,
                                          int seed,
                                          Db *dbgrid,
                                          int *count)
{
  int *indg, ndim, number, nbloc, ind;
  double *extend, densloc, densmax, density, test;
  VectorDouble tab,coor;

  /* Initializations */

  *count = number = 0;
  extend = (double *) NULL;
  indg = (int *) NULL;
  ndim = dbgrid->getNDim();
  if (!is_grid(dbgrid))
  {
    messerr("This function requires the Db organized as a grid");
    goto label_end;
  }
  if (dbgrid->isVariableNumberComparedTo(0))
  {
    messerr("This function requires the density to be defined in the Db");
    goto label_end;
  }
  law_set_random_seed(seed);

  /* Core allocation */

  indg = db_indg_alloc(dbgrid);
  if (indg == (int *) NULL) goto label_end;
  extend = (double *) mem_alloc(sizeof(double) * ndim, 0);
  if (extend == (double *) NULL) goto label_end;
  coor.resize(ndim);

  /* Calculate the volume of the grid */

  for (int idim = 0; idim < ndim; idim++)
    extend[idim] = dbgrid->getNX(idim) * dbgrid->getDX(idim);

  /* Look for the maximum intensity */

  densmax = density = 0.;
  for (int iech = 0; iech < dbgrid->getSampleNumber(); iech++)
  {
    if (!dbgrid->isActive(iech)) continue;
    densloc = dbgrid->getVariable(iech, 0);
    if (FFFF(densloc) || densloc < 0) continue;
    if (densmax < densloc) densmax = densloc;
    density += densloc;
  }
  density /= dbgrid->getSampleNumber();

  /* Draw the count of samples */

  number = st_get_number_poisson(verbose, ndim, extend, density);
  if (number <= 0) goto label_end;

  /* Core allocation */

  tab.resize(ndim * number);

  /* Point generation */

  nbloc = 0;
  while (nbloc < number)
  {
    for (int idim = 0; idim < ndim; idim++)
      coor[idim] = dbgrid->getX0(idim) + law_uniform(0., extend[idim]);

    /* Locate the new sample in the grid */

    ind = dbgrid->coordinateToRank(coor);
    if (ind < 0) continue;

    /* Read the local intensity */

    densloc = dbgrid->getVariable(ind, 0);
    if (!dbgrid->isActive(ind)) continue;
    if (FFFF(densloc) || densloc < 0) continue;

    /* Draw the acceptation-rejection criterion */

    test = law_uniform(0., densmax);
    if (test > densloc) continue;

    /* The sample is accepted */

    for (int idim = 0; idim < ndim; idim++)
      TAB(nbloc,idim) = coor[idim];
    nbloc++;
  }

  /* Return error code */

  *count = number;

  label_end: indg = db_indg_free(indg);
  extend = (double *) mem_free((char * ) extend);
  return (tab);
}

/*****************************************************************************/
/*!
 **  Create a set of samples according to a Poisson Thinning process
 **
 ** \return  Array of returned values (or NULL)
 **
 ** \param[in]  verbose     Verbose option
 ** \param[in]  ndim        Space dimension
 ** \param[in]  seed        Seed for the random number generator
 ** \param[in]  density     Average Poisson intensity
 ** \param[in]  range       Repulsion range
 ** \param[in]  beta        Bending coefficient
 ** \param[in]  origin      Vector of field origin
 ** \param[in]  extend      Vector of field extends
 **
 ** \param[out] count       Number of points generated
 **
 *****************************************************************************/
static VectorDouble st_point_init_poisthin(int verbose,
                                           int ndim,
                                           int seed,
                                           double density,
                                           double range,
                                           double beta,
                                           VectorDouble origin,
                                           VectorDouble extend,
                                           int *count)
{
  int *keep, number, lec, ecr;
  double ddmin, dd, delta, alea, proba;
  VectorDouble tab;

  /* Initializations */

  *count = number = 0;
  keep = (int *) NULL;
  law_set_random_seed(seed);

  /* Preliminary check */

  if (range <= 0)
  {
    messerr("Argument 'range' must be non negative");
    goto label_end;
  }

  /* Generate the Poisson point process */

  tab = st_point_init_poisson(verbose, ndim, seed, density, origin, extend,
                              &number);

  /* Core allocation */

  keep = (int *) mem_alloc(sizeof(int) * number, 0);
  if (keep == (int *) NULL) goto label_end;

  /* Operate the thining algorithm */

  for (int ip = 0; ip < number; ip++)
  {

    /* Find the closest data */

    ddmin = 1.e30;
    for (int jp = 0; jp < number; jp++)
    {
      if (ip == jp) continue;
      dd = 0.;
      for (int idim = 0; idim < ndim; idim++)
      {
        delta = TAB(ip,idim) - TAB(jp, idim);
        dd += delta * delta;
      }
      if (dd > ddmin) continue;
      ddmin = dd;
    }

    /* Check the rejection criterion */

    proba = exp(-pow(sqrt(ddmin) / range, beta));
    alea = law_uniform(0., 1.);
    keep[ip] = (alea > proba);
  }

  /* Discard the rejected samples */

  ecr = lec = 0;
  for (int ip = 0; ip < number; ip++)
  {
    if (keep[ip])
    {
      for (int idim = 0; idim < ndim; idim++)
        TAB(ecr,idim) = TAB(lec, idim);
      ecr++;
    }
    lec++;
  }
  if (number <= 0)
  {
    messerr("After thinning, the Point Process does not have any sample left");
    goto label_end;
  }

  /* Core reallocation */

  tab.resize(ndim * ecr);

  /* Optional complementary printout */

  if (verbose)
  {
    message("Repulsion range  = %lf\n", range);
    message("Final number     = %d\n", ecr);
  }

  /* Return error code */

  *count = ecr;

  label_end: keep = (int *) mem_free((char * ) keep);
  return (tab);
}

/****************************************************************************/
/*!
 **  Create indicator residual variables
 **
 ** \return  Error returned code
 **
 ** \param[in]  db      Db structure
 ** \param[in]  rank    Rank of the target variable
 ** \param[in]  ncut    Number of cutoffs
 ** \param[in]  zcut    Array containing the cutoffs
 **
 ** \remarks The array 'zcut' must be provided in increasing order
 **
 *****************************************************************************/
GEOSLIB_API int db_resind(Db *db, int rank, int ncut, double *zcut)
{
  double *tonnage, value, zval, ind_cut0, ind_cut1, ton_cut0, ton_cut1, ir;
  int ntot, nech, iptr;

  /* Initializations */

  nech = db->getSampleNumber();
  iptr = 0;
  tonnage = (double *) NULL;
  for (int icut = 1; icut < ncut; icut++)
  {
    if (zcut[icut] > zcut[icut - 1]) continue;
    messerr("The cutoffs must be provided in increasing order");
    goto label_end;
  }

  /* Core allocation */

  tonnage = (double *) mem_alloc(sizeof(double) * ncut, 0);
  if (tonnage == (double *) NULL) goto label_end;
  for (int icut = 0; icut < ncut; icut++)
    tonnage[icut] = 0;

  /* Calculate the tonnages */

  ntot = 0;
  for (int iech = 0; iech < nech; iech++)
  {
    if (!db->isActive(iech)) continue;
    value = db->getArray(iech, rank);
    if (FFFF(value)) continue;
    ntot++;

    for (int icut = 0; icut < ncut; icut++)
      if (value >= zcut[icut]) tonnage[icut]++;
  }
  for (int icut = 0; icut < ncut; icut++)
    tonnage[icut] /= (double) ntot;

  /* Create the variables */

  iptr = db->addFields(ncut, TEST);
  if (iptr < 0) goto label_end;

  /* Loop on the samples */

  for (int iech = 0; iech < nech; iech++)
  {
    if (!db->isActive(iech)) continue;
    value = db->getArray(iech, rank);
    if (FFFF(value)) continue;

    /* Loop on the cutoffs */

    for (int icut = 0; icut < ncut; icut++)
    {
      zval = zcut[icut];
      ind_cut0 = (value > zval);
      zval = (icut > 0) ? zcut[icut - 1] :
                          0.;
      ind_cut1 = (value > zval);
      ton_cut0 = tonnage[icut];
      ton_cut1 = (icut > 0) ? tonnage[icut - 1] :
                              1.;
      ir = ind_cut0 / ton_cut0 - ind_cut1 / ton_cut1;
      db->setArray(iech, iptr + icut, ir);
    }
  }

  label_end: tonnage = (double *) mem_free((char * ) tonnage);
  return (0);
}

/****************************************************************************/
/*!
 **  Normalize the gradient components
 **
 ** \param[in]  dbgrid  Db structure (grid organized)
 **
 *****************************************************************************/
static void st_gradient_normalize(Db *dbgrid)

{
  double norme, grad;
  int ndim;

  /* Initializations */

  ndim = dbgrid->getNDim();

  /* Loop on the samples */

  for (int iech = 0; iech < dbgrid->getSampleNumber(); iech++)
  {

    norme = 0.;
    for (int idim = 0; idim < ndim; idim++)
    {
      grad = dbgrid->getGradient(iech, idim);
      norme += grad * grad;
    }

    if (norme <= 0) continue;
    norme = sqrt(norme);

    for (int idim = 0; idim < ndim; idim++)
    {
      grad = dbgrid->getGradient(iech, idim);
      dbgrid->setGradient(iech, idim, grad / norme);
    }
  }
}

/****************************************************************************/
/*!
 **  Calculate the gradient over a grid
 **
 ** \return  Rank of the newly created variables (or -1)
 **
 ** \param[in]  dbgrid  Db structure (grid organized)
 **
 *****************************************************************************/
GEOSLIB_API int db_gradient_components(Db *dbgrid)

{
  int *indg, iptrz, iptr, nx, ny, nz, nmax, error, ndim, j1, j2, number;
  double dinc, v1, v2, delta;

  /* Preliminary check */

  error = number = 1;
  iptrz = iptr = -1;
  indg = (int *) NULL;
  ndim = dbgrid->getNDim();
  if (!is_grid(dbgrid))
  {
    messerr("The Db should be organized as a Grid");
    goto label_end;
  }
  if (! dbgrid->isVariableNumberComparedTo(1)) goto label_end;
  if (ndim > 3)
  {
    messerr("This function is limited to Space Dimension <= 3");
    goto label_end;
  }

  /* Initializations */

  nx = dbgrid->getNX(0);
  ny = dbgrid->getNX(1);
  nz = dbgrid->getNX(2);
  indg = db_indg_alloc(dbgrid);
  if (indg == (int *) NULL) goto label_end;

  /* Create the new variable */

  iptrz = dbgrid->getColumnByLocator(LOC_Z, 0);
  if (iptrz < 0) goto label_end;
  iptr = dbgrid->addFields(ndim,TEST,String(),LOC_G);

  /* Calculate the Gradient components */

  for (int ix = 0; ix < nx; ix++)
    for (int iy = 0; iy < ny; iy++)
      for (int iz = 0; iz < nz; iz++)
      {
        for (int idim = 0; idim < ndim; idim++)
        {
          nmax = dbgrid->getNX(idim);
          dinc = dbgrid->getDX(idim);

          v1 = v2 = 0.;
          if (idim == 0)
          {
            j1 = (ix + 1 > nmax - 1) ? ix :
                                       ix + 1;
            v1 = get_grid_value(dbgrid, iptrz, indg, j1, iy, iz);
            if (FFFF(v1)) continue;
            j2 = (ix - 1 < 0) ? ix :
                                ix - 1;
            v2 = get_grid_value(dbgrid, iptrz, indg, j2, iy, iz);
            if (FFFF(v2)) continue;
            number = j1 - j2;
          }

          if (idim == 1)
          {
            j1 = (iy + 1 > nmax - 1) ? iy :
                                       iy + 1;
            v1 = get_grid_value(dbgrid, iptrz, indg, ix, j1, iz);
            if (FFFF(v1)) continue;
            j2 = (iy - 1 < 0) ? iy :
                                iy - 1;
            v2 = get_grid_value(dbgrid, iptrz, indg, ix, j2, iz);
            if (FFFF(v2)) continue;
            number = j1 - j2;
          }

          if (idim == 2)
          {
            j1 = (iz + 1 > nmax - 1) ? iz :
                                       iz + 1;
            v1 = get_grid_value(dbgrid, iptrz, indg, ix, iy, j1);
            if (FFFF(v1)) continue;
            j2 = (iz - 1 < 0) ? iz :
                                iz - 1;
            v2 = get_grid_value(dbgrid, iptrz, indg, ix, iy, j2);
            if (FFFF(v2)) continue;
            number = j1 - j2;
          }

          delta = (v1 - v2) / (number * dinc);
          set_grid_value(dbgrid, iptr + idim, indg, ix, iy, iz, delta);
        }
      }

  /* Set the error return code */

  error = 0;

  label_end: if (error)
  {
    (void) db_attribute_del_mult(dbgrid, iptr, ndim);
    iptr = -1;
  }
  indg = db_indg_free(indg);
  return (iptr);
}

/****************************************************************************/
/*!
 **  Check if one (at least) of the gradient components is undefined
 **
 ** \return  1 If one component (at least) is undefined
 **
 ** \param[in]  dbgrid    Db structure (grid organized)
 ** \param[in]  iptr_grad Rank of the first gradient component
 ** \param[in]  iech      Sample rank
 **
 *****************************************************************************/
static int st_is_undefined(Db *dbgrid, int iptr_grad, int iech)
{
  int ndim;

  ndim = dbgrid->getNDim();
  for (int idim = 0; idim < ndim; idim++)
  {
    if (FFFF(dbgrid->getArray(iech, iptr_grad + idim))) return (1);
  }
  return (0);
}

/****************************************************************************/
/*!
 **  Check if the gradient is null
 **
 ** \return  1 If gradient is null
 **
 ** \param[in]  dbgrid    Db structure (grid organized)
 ** \param[in]  iptr_grad Rank of the first gradient component
 ** \param[in]  iech      Sample rank
 **
 *****************************************************************************/
static int st_is_zero(Db *dbgrid, int iptr_grad, int iech)
{
  double grad, delta;
  int ndim;

  grad = 0.;
  ndim = dbgrid->getNDim();
  for (int idim = 0; idim < ndim; idim++)
  {
    delta = dbgrid->getArray(iech, iptr_grad + idim);
    grad += delta * delta;
  }
  return (grad < 1.e-5);
}

/****************************************************************************/
/*!
 **  Get the next gradient-based data
 **
 ** \return  1 Error return code
 **
 ** \param[in]  dbgrid    Db structure (grid organized)
 ** \param[in]  iptr_grad Rank of the first gradient component
 ** \param[in]  coor      Array of coordinates
 **
 ** \param[out] indg      Array of indices
 ** \param[out] knd       Absolute index
 ** \param[out] surf      Local value of the surface
 **
 *****************************************************************************/
static int st_get_next(Db *dbgrid,
                       int iptr_grad,
                       VectorDouble& coor,
                       int *indg,
                       int *knd,
                       double *surf)
{
  int knd_loc;
  double surf_loc;

  knd_loc = dbgrid->coordinateToRank(coor);
  if (knd_loc < 0) return 1;
  if (!dbgrid->isActive(knd_loc)) return 1;
  surf_loc = dbgrid->getVariable(knd_loc, 0);
  if (FFFF(surf_loc) || st_is_undefined(dbgrid, iptr_grad, knd_loc)) return (1);
  if (st_is_zero(dbgrid, iptr_grad, knd_loc)) return (1);
  *knd = knd_loc;
  *surf = surf_loc;
  return 0;
}

/****************************************************************************/
/*!
 **  Calculate the streamlines
 **
 ** \return  Error return code
 **
 ** \param[in]  dbgrid    Db structure (grid organized)
 ** \param[in]  dbpoint   Db structure for control points
 ** \param[in]  niter     Maximum number of iterations
 ** \param[in]  step      Progress step
 ** \param[in]  flag_norm 1 if the gradients must be normalized
 ** \param[in]  use_grad  1 if the existing gradients must be used
 **                       0 the gradients must be calculated here
 ** \param[in]  save_grad 1 save the gradients generated in this function
 **                       0 delete gradients when calculated here
 **
 ** \param[out] nbline_loc Number of streamline steps
 ** \param[out] npline_loc Number of values per line vertex
 ** \param[out] line_loc   Array of streamline steps (Dimension: 5 * nbline)
 **
 ** \remarks The returned array 'line_loc' must be freed by the calling function
 ** \remarks Use get_keypone("Streamline_Skip",1) to define the skipping ratio
 **
 *****************************************************************************/
GEOSLIB_API int db_streamline(Db *dbgrid,
                              Db *dbpoint,
                              int niter,
                              double step,
                              int flag_norm,
                              int use_grad,
                              int save_grad,
                              int *nbline_loc,
                              int *npline_loc,
                              double **line_loc)
{
  int *indg, error, npline, idim, ecr;
  int iptr_time, iptr_accu, iptr_grad, nbline, knd, nquant, nbyech, ndim;
  double *coor0, *line, surf, date;
  static int quant = 1000;
  VectorDouble coor;

  /* Initializations */

  error = 1;
  nbline = nquant = 0;
  iptr_grad = iptr_accu = iptr_time = -1;
  indg = (int *) NULL;
  coor0 = line = (double *) NULL;
  if (dbpoint == (Db *) NULL) dbpoint = dbgrid;
  nbyech = (int) get_keypone("Streamline_Skip", 1.);

  /* Preliminary checks */

  ndim = dbgrid->getNDim();
  if (ndim < 2 || ndim > 3)
  {
    messerr("This function is limited to 2-D or 3-D case");
    goto label_end;
  }
  npline = ndim + 3;

  /* Core allocation on the Grid Db */

  indg = db_indg_alloc(dbgrid);
  if (indg == (int *) NULL) goto label_end;
  coor.resize(ndim);
  coor0 = db_sample_alloc(dbgrid, LOC_X);
  if (coor0 == (double *) NULL) goto label_end;
  iptr_time = dbgrid->addFields(1, TEST);
  if (iptr_time < 0) goto label_end;
  iptr_accu = dbgrid->addFields(1, 0.);
  if (iptr_accu < 0) goto label_end;

  /* Calculate the gradient */

  if (use_grad)
  {
    if (dbgrid->getGradientNumber() != ndim)
    {
      messerr("When using the option 'use.grad'");
      messerr("the number of gradients should be %d (%d)", ndim,
              dbgrid->getGradientNumber());
      goto label_end;
    }
    iptr_grad = dbgrid->getColumnByLocator(LOC_G, 0);
  }
  else
  {
    iptr_grad = db_gradient_components(dbgrid);
  }
  if (iptr_grad < 0) goto label_end;

  /* Normalize the gradient (optional) */

  if (flag_norm) st_gradient_normalize(dbgrid);

  /* Loop on the drop points */

  for (int iech = 0; iech < dbpoint->getSampleNumber(); iech++)
  {
    if (!dbpoint->isActive(iech)) continue;
    if (iech % nbyech != 0) continue;
    db_sample_load(dbpoint, LOC_X, iech, coor.data());
    if (st_get_next(dbgrid, iptr_grad, coor, indg, &knd, &surf)) break;

    /* Store the new point in the Streamline */

    if (nbline >= nquant * quant)
    {
      nquant++;
      line = (double *) mem_realloc((char * ) line,
                                    sizeof(double) * npline * nquant * quant,1);
    }
    for (idim = ecr = 0; idim < ndim; idim++)
      LINE(nbline,ecr++) = coor[idim];
    LINE(nbline,ecr++) = surf;
    LINE(nbline,ecr++) = knd + 1.;
    LINE(nbline,ecr++) = 0.;
    nbline++;

    for (int i = 0; i < niter; i++)
    {
      for (idim = ecr = 0; idim < ndim; idim++)
        coor[idim] -= step * dbgrid->getArray(knd, iptr_grad + idim);
      if (st_get_next(dbgrid, iptr_grad, coor, indg, &knd, &surf)) break;

      /* Store the new point in the Streamline */

      if (nbline >= nquant * quant)
      {
        nquant++;
        line = (double *) mem_realloc((char * ) line,
                                      sizeof(double) * npline * nquant * quant,
                                      1);
      }
      for (idim = ecr = 0; idim < ndim; idim++)
        LINE(nbline,ecr++) = coor[idim];
      LINE(nbline,ecr++) = surf;
      LINE(nbline,ecr++) = knd + 1.;
      LINE(nbline,ecr++) = i + 1.;
      nbline++;

      /* Update variables in the grid Db */

      date = MIN(dbgrid->getArray(knd, iptr_time), i + 1.);
      dbgrid->setArray(knd, iptr_time, date);
      dbgrid->updArray(knd, iptr_accu, 0, 1.);
    }

    /* Add the endpoint */

    if (nbline >= nquant * quant)
    {
      nquant++;
      line = (double *) mem_realloc((char * ) line,
                                    sizeof(double) * npline * nquant * quant,
                                    1);
    }
    for (idim = ecr = 0; idim < ndim; idim++)
      LINE(nbline,ecr++) = TEST;
    LINE(nbline,ecr++) = TEST;
    LINE(nbline,ecr++) = TEST;
    LINE(nbline,ecr++) = TEST;
    nbline++;
  }

  /* Final reallocation */

  line = (double *) mem_realloc((char * ) line,
                                sizeof(double) * npline * nbline, 1);

  /* Set the error return code */

  *nbline_loc = nbline;
  *npline_loc = npline;
  *line_loc = line;
  error = 0;

  label_end: indg = db_indg_free(indg);
  coor0 = db_sample_free(coor0);
  if (!use_grad && !save_grad && iptr_grad >= 0)
    db_attribute_del_mult(dbgrid, iptr_grad, ndim);
  return (error);
}

/*****************************************************************************/
/*!
 **  Calculate and store new variables in the Db which contain
 **  the non-stationary Model component
 **
 ** \return  Distance value
 **
 ** \param[in]  db          Db structure
 ** \param[in]  model       Model structure
 ** \param[in]  icov        Rank of the Basic structure
 ** \param[in]  namconv     Naming convention
 **
 ** \remarks This procedure automatically creates several fields:
 ** \remarks ndim fields for storing the ranges
 ** \remarks ndim fields for storing the angles
 ** \remarks 1 field for storing the sill
 **
 *****************************************************************************/
GEOSLIB_API int db_model_nostat(Db *db,
                                Model *model,
                                int icov,
                                NamingConvention namconv)
{
  if (icov < 0 || icov >= model->getCovaNumber()) return 1;
  if (! model->isNoStat()) return 0;

  // The Non-stationary must be defined in the tabulated way
  if (manage_nostat_info(1, model, db, nullptr)) return 1;

  /* Create the new variables */

  int ndim = model->getDimensionNumber();
  CovInternal covint(1,-1,1,-1,ndim,db,db);
  int iptr = db->addFields(2 * ndim + 1, 0.);
  if (iptr < 0) return 1;

  /* Loop on the samples */

  for (int iech = 0; iech < db->getSampleNumber(); iech++)
  {
    if (!db->isActive(iech)) continue;

    /* Load the non_stationary parameters */

    covint.setIech1(iech);
    covint.setIech2(iech);
    model_nostat_update(&covint, model);
    CovAniso* cova = model->getCova(icov);

    /* Store the variables */

    int jptr = iptr;
    for (int idim = 0; idim < ndim; idim++)
    {
      db->setArray(iech, jptr, cova->getRange(idim));
      jptr++;
    }
    for (int idim = 0; idim < ndim; idim++)
    {
      db->setArray(iech, jptr, cova->getAnisoAngles(idim));
      jptr++;
    }
    db->setArray(iech, jptr++, cova->getSill(0, 0));
  }

  // Naming convention

  int jptr = iptr;
  for (int idim = 0; idim < ndim; idim++)
    namconv.setNamesAndLocators(nullptr, LOC_UNKNOWN, -1, db, jptr++,
                                concatenateStrings("-","Range",intToString(idim+1)));
  for (int idim = 0; idim < ndim; idim++)
    namconv.setNamesAndLocators(nullptr, LOC_UNKNOWN, -1, db, jptr++,
                                concatenateStrings("-","Angle",intToString(idim+1)));
  namconv.setNamesAndLocators(nullptr, LOC_UNKNOWN, -1, db, jptr++, "Sill");
  namconv.setLocators(db, iptr, 1, 2*ndim+1);

  (void) manage_nostat_info(-1, model, db, nullptr);
  return 0;
}

/*****************************************************************************/
/*!
 **  Smooth out the VPC
 **
 ** \return  Distance value
 **
 ** \param[in]  db          3-D Db structure containing the VPCs
 ** \param[in]  width       Width of the Filter
 ** \param[in]  range       Range of the Gaussian Weighting Function
 **
 ** \remarks Work is performed IN PLACE
 **
 *****************************************************************************/
GEOSLIB_API int db_smooth_vpc(Db *db, int width, double range)
{
  int iz, nz, nprop, ecr, nkern, jz, error;
  double *prop1, *prop2, *kernel, total, propval, dz, quant, quant0;

  /* Initializations */

  error = 1;
  nprop = db->getProportionNumber();
  nz = db->getNX(2);
  dz = db->getDX(2);
  prop1 = prop2 = kernel = (double *) NULL;

  /* Core allocation */

  quant0 = law_invcdf_gaussian(0.975);
  if (FFFF(range))
    range = dz * width / quant0;
  else if (IFFFF(width))
    width = (int) (range * quant0 / dz);
  else
  {
    messerr("You must define either 'width' or 'range'");
    goto label_end;
  }
  nkern = 2 * width + 1;
  prop1 = (double *) mem_alloc(sizeof(double) * nz * nprop, 1);
  prop2 = (double *) mem_alloc(sizeof(double) * nz * nprop, 1);
  kernel = (double *) mem_alloc(sizeof(double) * nkern, 1);

  /* Establish the Kernel */

  for (int i = -width; i <= width; i++)
  {
    quant = (i * dz) / range;
    kernel[i + width] = law_df_gaussian(quant) / range;
  }

  /* Preliminary checks */

  if (!is_grid(db) || db->getNDim() != 3) goto label_end;

  /* Loop on the 2-D grid cells */

  for (int ix = 0; ix < db->getNX(0); ix++)
    for (int iy = 0; iy < db->getNX(1); iy++)
    {

      /* Load the proportions */

      if (db_prop_read(db, ix, iy, prop1)) goto label_end;

      /* Loop on the proportions */

      for (int iprop = 0; iprop < nprop; iprop++)
      {

        /* Loop on the levels of the VPC */

        for (iz = ecr = 0; iz < nz; iz++)
        {

          /* Loop on the kernel items */

          total = 0.;
          for (int i = -width; i <= width; i++)
          {
            jz = get_mirror_sample(nz, iz + i);
            propval = PROP1(jz, iprop);
            total += kernel[i + width] * propval;
          }
          PROP2(iz,iprop) = total;
        }
      }
      if (db_prop_write(db, ix, iy, prop2)) goto label_end;
    }

  /* Set the error return code */

  error = 0;

  label_end: prop1 = (double *) mem_free((char * ) prop1);
  prop2 = (double *) mem_free((char * ) prop2);
  kernel = (double *) mem_free((char * ) kernel);
  return (error);
}

/*****************************************************************************/
/*!
 **  Extract a new Db from an old Db using a sample seleciton vector
 **
 ** \return  Pointer to the newly created Db
 **
 ** \param[in]  db          Initial Db
 ** \param[in]  ranks       Vector of selected ranks
 **
 ** \remarks The 'ranks' array is dimensionned to the number of samples of 'db'
 ** \remarks and ranks[i]=1 if the sample 'i' must be kept and 0 otherwise.
 ** \remarks All the variables contained in the input Db are copied into the
 ** \remarks output Db
 ** \remarks This operations is limited to non-grid Db
 ** \remarks Possible selection in the input Db is taken into account
 **
 *****************************************************************************/
GEOSLIB_API Db *db_extract(Db *db, int *ranks)
{
  int *iatts, nech, nech_all, natt, error, ecr;
  VectorDouble tab;
  Db *dbnew;

  // Initializations

  error = 1;
  dbnew = (Db *) NULL;
  iatts = (int *) NULL;
  if (db == (Db *) NULL || ranks == (int *) NULL) return (dbnew);
  nech_all = db->getSampleNumber();
  natt = db->getFieldNumber();

  // Count the number of samples in 'ranks'
  nech = 0;
  for (int iech = 0; iech < nech_all; iech++)
  {
    if (!db->isActive(iech)) continue;
    if (ranks[iech] != 0) nech++;
  }
  if (nech <= 0)
  {
    messerr("The extraction failed as no sample is selected");
    return (dbnew);
  }

  // Fill the array of selected variables
  // TODO Enable selecting variables by LOC type

  iatts = (int *) mem_alloc(sizeof(int) * natt, 0);
  if (iatts == (int *) NULL) goto label_end;
  for (int iatt = 0; iatt < natt; iatt++)
    iatts[iatt] = 1;

  // Core allocation 

  tab.resize(nech * natt);

  // Create the extracted array

  ecr = 0;
  for (int iech = 0; iech < nech_all; iech++)
  {
    if (!db->isActive(iech)) continue;
    if (ranks[iech] == 0) continue;
    for (int iatt = 0; iatt < natt; iatt++)
      tab[ecr++] = db->getArray(iech, iatt);
  }

  // Create the new db

  dbnew = db_create_point(nech, natt, LOAD_BY_SAMPLE, 1, tab);
  if (dbnew == (Db *) NULL) goto label_end;

  // Delete the unnecessary variables 

  for (int iatt = natt - 1; iatt >= 0; iatt--)
  {
    if (iatts[iatt] != 0) continue;
    dbnew->deleteFieldByAttribute(iatt + 1);
  }

  // Set the error return code

  error = 0;

  label_end: iatts = (int *) mem_free((char * ) iatts);
  if (error) dbnew = db_delete(db);
  return (dbnew);
}

/*****************************************************************************/
/*!
 **  Regularize variables along vertical wells
 **
 ** \return  Pointer to the newly created Db
 **
 ** \param[in]  db          Initial Db
 ** \param[in]  dbgrid      Reference Grid
 ** \param[in]  flag_center When TRUE, the sample is centered in the layer
 **                         to which it belongs
 **
 ** \remarks This function requires the input well ('db') and the grid to be
 ** \remarks defined in space >= 3D
 ** \remarks It requires a CODE variable to be defined in the input 'db'
 ** \remarks This function regularizes all the variables marked with a Z-locator
 ** \remarks This function takes a sample into account only if isotopic
 **
 *****************************************************************************/
GEOSLIB_API Db *db_regularize(Db *db, Db *dbgrid, int flag_center)
{
  int ncode, nvar, nz, ndim, not_defined, ecr, nech, size, ntot, iz, icode;
  double *wcnt, *wtab, *wcor, *coor, code, ratio;
  Db *dbnew;
  VectorDouble wecr;
  VectorDouble codes;

  // Initializations

  dbnew = (Db *) NULL;
  wtab = wcor = coor = wcnt = (double *) NULL;
  if (db == (Db *) NULL || dbgrid == (Db *) NULL) return (dbnew);

  // Preliminary checks */

  if (!is_grid(dbgrid))
  {
    messerr("This function requires 'dbgrid' to correspond to a Grid");
    return (dbnew);
  }

  if (db->getNDim() < 3)
  {
    messerr("This function requires the 'db' to be defined in 3D or more");
    return (dbnew);
  }

  if (dbgrid->getNDim() < 3)
  {
    messerr("This function requires the 'dbgrid' to be defined in 3D or more");
    return (dbnew);
  }

  if (!db->hasCode())
  {
    messerr("This function requires the definition of a CODE variable in 'db'");
    return (dbnew);
  }

  if (db->isVariableNumberComparedTo(0))
  {
    messerr("You should define some Z-variables in input 'db'");
    return (dbnew);
  }

  // Core allocation 

  nz = dbgrid->getNX(2);
  nvar = db->getVariableNumber();
  ndim = db->getNDim();
  size = ndim + nvar + 1;

  codes = db->getCodeList();
  ncode = static_cast<int> (codes.size());
  coor = (double *) mem_alloc(sizeof(double) * ndim, 0);
  if (coor == (double *) NULL) goto label_end;

  wcnt = (double *) mem_alloc(sizeof(double) * ncode * nz, 0);
  if (wcnt == (double *) NULL) goto label_end;

  wcor = (double *) mem_alloc(sizeof(double) * ncode * nz * ndim, 0);
  if (wcor == (double *) NULL) goto label_end;

  wtab = (double *) mem_alloc(sizeof(double) * ncode * nz * nvar, 0);
  if (wtab == (double *) NULL) goto label_end;

  /* Initialize the different arrays */

  for (int i = 0; i < ncode * nz; i++)
    wcnt[i] = 0.;
  for (int i = 0; i < ncode * nz * ndim; i++)
    wcor[i] = 0.;
  for (int i = 0; i < ncode * nz * nvar; i++)
    wtab[i] = 0.;

  // Loop on the different samples

  ntot = db->getSampleNumber();

  //message("Before regularization: ncode = %d, nz = %d, ntot = %d\n", (int)ncode, (int)nz, (int)ntot);

  for (int iech = 0; iech < ntot; iech++)
  {
    if (!db->isActive(iech)) continue;
    mes_process("Regularize Wells", ntot, iech);
    code = db->getCode(iech);

    // Identify the rank of the code

    icode = -1;
    for (int i = 0; i < ncode && icode < 0; i++)
      if (code == codes[i]) icode = i;
    if (icode < 0) continue;

    // Load the coordinates

    for (int idim = 0; idim < ndim; idim++)
      coor[idim] = db->getCoordinate(iech, idim);

    int err = point_to_bench(dbgrid, coor, 0, &iz);
    if (err < 0) continue;
    if (iz < 0 || iz >= nz) continue;

    // Check if all variables are defined

    not_defined = 0;
    for (int ivar = 0; ivar < nvar && not_defined == 0; ivar++)
      if (FFFF(db->getVariable(iech, ivar))) not_defined = 1;
    if (not_defined) continue;

    // Cumulate this sample

    WCNT(iz,icode) += 1.;
    for (int idim = 0; idim < ndim; idim++)
      WCOR(iz,icode,idim) += db->getCoordinate(iech, idim);
    for (int ivar = 0; ivar < nvar; ivar++)
      WTAB(iz,icode,ivar) += db->getVariable(iech, ivar);
  }

  // Normalization

  nech = 0;
  for (icode = 0; icode < ncode; icode++)
    for (iz = 0; iz < nz; iz++)
    {
      ratio = WCNT(iz, icode);
      if (ratio <= 0) continue;
      for (int idim = 0; idim < ndim; idim++)
        WCOR(iz,icode,idim) /= ratio;
      if (flag_center)
      WCOR(iz,icode,2) = dbgrid->getX0(2) + (0.5 + iz) * dbgrid->getDX(2);
      for (int ivar = 0; ivar < nvar; ivar++)
        WTAB(iz,icode,ivar) /= ratio;
      nech++;
    }

  //message("After normalization: nech = %d, size = %d\n", (int)nech, (int)size);

  // Load in storing array

  wecr.resize(size * nech);

  ecr = 0;
  for (icode = 0; icode < ncode; icode++)
    for (iz = 0; iz < nz; iz++)
    {
      ratio = WCNT(iz, icode);
      if (ratio <= 0) continue;
      for (int idim = 0; idim < ndim; idim++)
        wecr[ecr++] = WCOR(iz, icode, idim);
      wecr[ecr++] = codes[icode];
      for (int ivar = 0; ivar < nvar; ivar++)
        wecr[ecr++] = WTAB(iz, icode, ivar);
    }

  // Create the new db

  dbnew = db_create_point(nech, size, LOAD_BY_SAMPLE, 0, wecr);
  if (dbnew == (Db *) NULL) goto label_end;

  ecr = 0;
  dbnew->setLocatorsByAttribute(ndim, ecr, LOC_X);
  ecr += ndim;
  dbnew->setLocatorByAttribute(ecr, LOC_C);
  ecr += 1;
  dbnew->setLocatorsByAttribute(nvar, ecr, LOC_Z);
  ecr += nvar;

  label_end:
  coor = (double *) mem_free((char * ) coor);
  wcnt = (double *) mem_free((char * ) wcnt);
  wcor = (double *) mem_free((char * ) wcor);
  wtab = (double *) mem_free((char * ) wtab);
  return (dbnew);
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
 *****************************************************************************/
GEOSLIB_API double *db_grid_sampling(Db *dbgrid,
                                     double *x1,
                                     double *x2,
                                     int ndisc,
                                     int ncut,
                                     double *cuts,
                                     int *nval_ret)
{
  double *xi1, *xi2, *res, delta, vi1, vi2, cut, v1, v2;
  int ndim, iatt, nval;

  /* Initializations */

  *nval_ret = 0;
  res = xi1 = xi2 = (double *) NULL;
  ndim = dbgrid->getNDim();
  iatt = dbgrid->getColumnByLocator(LOC_Z, 0);

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

  xi1 = (double *) mem_alloc(sizeof(double) * ndim, 0);
  if (xi1 == (double *) NULL) goto label_end;
  xi2 = (double *) mem_alloc(sizeof(double) * ndim, 0);
  if (xi2 == (double *) NULL) goto label_end;

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
      res = (double *) mem_realloc((char * ) res,
                                   (ndim + 1) * (nval + 1) * sizeof(double), 0);
      if (res == (double *) NULL) goto label_end;

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

  label_end: xi1 = (double *) mem_free((char * ) xi1);
  xi2 = (double *) mem_free((char * ) xi2);
  return (res);
}

/*****************************************************************************/
/*!
 **  Sampling a fine grid in a coarser set of cells
 **
 ** \return  Error returned code
 **
 ** \param[in]  dbgrid      reference Grid
 ** \param[in]  nvar        Number of variables
 ** \param[in]  vars        Array of variable ranks
 ** \param[in]  npacks      Vector of packing factors
 ** \param[in]  npcell      Number of samples per cell
 ** \param[in]  nmini       Minimum number of nodes before drawing
 **
 ** \param[out] nech_ret    Number of selected samples
 ** \param[out] coor_ret    Array of coordinates
 ** \param[out] data_ret    Array of variables
 **
 ** \remarks The returned arrays 'coor' and 'data' must be freed by
 ** \remarks the calling function
 **
 *****************************************************************************/
GEOSLIB_API int db_grid2point_sampling(Db *dbgrid,
                                       int nvar,
                                       int *vars,
                                       int *npacks,
                                       int npcell,
                                       int nmini,
                                       int *nech_ret,
                                       double **coor_ret,
                                       double **data_ret)
{
  int ndim, ntotal, nech, indg[3], nret, nfine, iech, ecrc, ecrd, error;
  int *ranks, *retain;
  double *coor, *data, *rndval;

  // Initializations

  *nech_ret = 0;

  error = 1;
  nech = 0;
  coor = data = rndval = (double *) NULL;
  ranks = retain = (int *) NULL;
  ndim = dbgrid->getNDim();
  nfine = dbgrid->getSampleNumber();
  nmini = MAX(nmini, npcell);
  if (ndim > 3)
  {
    messerr("This function is limited to 3D or less");
    goto label_end;
  }

  // Core allocation 

  ntotal = 1;
  for (int idim = 0; idim < ndim; idim++) ntotal *= npacks[idim];
  rndval = (double *) mem_alloc(sizeof(double) * ntotal, 0);
  if (rndval == (double *) NULL) goto label_end;
  ranks = (int *) mem_alloc(sizeof(int) * ntotal, 0);
  if (ranks == (int *) NULL) goto label_end;
  retain = (int *) mem_alloc(sizeof(int) * nfine, 0);
  if (retain == (int *) NULL) goto label_end;

  // Dispatch

  nret = 0;
  if (ndim == 1)
  {
    for (int ixcell = 0; ixcell < dbgrid->getNX(0); ixcell += npacks[0])
    {

      // Collect eligible samples

      nech = 0;
      for (int ix = 0; ix < npacks[0]; ix++)
      {
        indg[0] = ixcell + ix;
        if (indg[0] >= dbgrid->getNX(0)) break;
        iech = db_index_grid_to_sample(dbgrid, indg);
        if (dbgrid->isActive(iech)) ranks[nech++] = iech;
      }
      if (nech < nmini) continue;

      // Draw sample(s) at random 

      for (int i = 0; i < nech; i++)
        rndval[i] = law_uniform(0., 1.);
      ut_sort_double(0, nech, ranks, rndval);
      for (int i = 0; i < npcell; i++)
        retain[nret++] = ranks[i];
    }
  }
  else if (ndim == 2)
  {
    for (int ixcell = 0; ixcell < dbgrid->getNX(0); ixcell += npacks[0])
      for (int iycell = 0; iycell < dbgrid->getNX(1); iycell += npacks[1])
      {

        // Collect eligible samples

        nech = 0;
        for (int ix = 0; ix < npacks[0]; ix++)
          for (int iy = 0; iy < npacks[1]; iy++)
          {
            indg[0] = ixcell + ix;
            if (indg[0] >= dbgrid->getNX(0)) break;
            indg[1] = iycell + iy;
            if (indg[1] >= dbgrid->getNX(1)) break;
            iech = db_index_grid_to_sample(dbgrid, indg);
            if (dbgrid->isActive(iech)) ranks[nech++] = iech;
          }
        if (nech < nmini) continue;

        // Draw sample(s) at random 

        for (int i = 0; i < nech; i++)
          rndval[i] = law_uniform(0., 1.);
        ut_sort_double(0, nech, ranks, rndval);
        for (int i = 0; i < npcell; i++)
          retain[nret++] = ranks[i];
      }
  }
  else
  {
    for (int ixcell = 0; ixcell < dbgrid->getNX(0); ixcell += npacks[0])
      for (int iycell = 0; iycell < dbgrid->getNX(1); iycell += npacks[1])
        for (int izcell = 0; izcell < dbgrid->getNX(2); izcell += npacks[2])
        {

          // Collect eligible samples

          nech = 0;
          for (int ix = 0; ix < npacks[0]; ix++)
            for (int iy = 0; iy < npacks[1]; iy++)
              for (int iz = 0; iz < npacks[2]; iz++)
              {
                indg[0] = ixcell + ix;
                if (indg[0] >= dbgrid->getNX(0)) break;
                indg[1] = iycell + iy;
                if (indg[1] >= dbgrid->getNX(1)) break;
                indg[2] = izcell + iz;
                if (indg[2] >= dbgrid->getNX(2)) break;
                iech = db_index_grid_to_sample(dbgrid, indg);
                if (dbgrid->isActive(iech)) ranks[nech++] = iech;
              }
          if (nech < nmini) continue;

          // Draw sample(s) at random 

          for (int i = 0; i < nech; i++)
            rndval[i] = law_uniform(0., 1.);
          ut_sort_double(0, nech, ranks, rndval);
          for (int i = 0; i < npcell; i++)
            retain[nret++] = ranks[i];
        }
  }
  rndval = (double *) mem_free((char * ) rndval);
  ranks = (int *) mem_free((char * ) ranks);

  // Allocate the array for coordinates and data

  coor = (double *) mem_alloc(sizeof(double) * ndim * nret, 0);
  if (coor == (double *) NULL) goto label_end;
  data = (double *) mem_alloc(sizeof(double) * nvar * nret, 0);
  if (data == (double *) NULL) goto label_end;

  // Load the returned arrays

  ecrc = ecrd = 0;
  for (int i = 0; i < nret; i++)
  {
    iech = retain[i];
    for (int idim = 0; idim < ndim; idim++)
      coor[ecrc++] = dbgrid->getCoordinate(iech, idim);
    for (int ivar = 0; ivar < nvar; ivar++)
      data[ecrd++] = dbgrid->getArray(iech, vars[ivar]);
  }

  // Set the error return code

  *nech_ret = nret;
  *coor_ret = coor;
  *data_ret = data;
  error = 0;

  // Core deallocation

  label_end: retain = (int *) mem_free((char * ) retain);
  ranks = (int *) mem_free((char * ) ranks);
  rndval = (double *) mem_free((char * ) rndval);
  return (error);
}

/*****************************************************************************/
/*!
 **  Determine the distance to a polyline
 **
 ** \return  Error returned code
 **
 ** \param[in]  db      Db structure
 ** \param[in]  polygon Polygons structure
 ** \param[in]  dmax    Maximum distance
 ** \param[in]  scale   Scaling option
 **                     0 : no scaling
 **                    >0 : scaling between 0 and 1
 **                    <0 : scaling between 1 and 0
 ** \param[in]  polin   Option for checking against the polygon
 **                     0 : no check
 **                    >0 : if sample is outside polygon, return TEST
 **                    <0 : if sample is inside polygon, return TEST
 **
 ** \remarks When patching values with respect to the polygon, when abs(polin):
 ** \remarks 1 : put NA
 ** \remarks 2 : put minimum distance
 ** \remarks 3 : put maximum distance
 **
 *****************************************************************************/
GEOSLIB_API int db_polygon_distance(Db *db,
                                    Polygons *polygon,
                                    double dmax,
                                    int scale,
                                    int polin)
{
  PL_Dist *pldist;
  double distmin, distloc, distmax, value, valtest;
  int iptr, nech, inside;

  // Initializations

  iptr = -1;
  nech = db->getSampleNumber();
  distmin = distmax = 0.;

  // Create a new attribute

  iptr = db->addFields(1, TEST);
  if (iptr < 0) return (1);

  // Loop on the polysets 

  distmin = TEST;
  for (int iset = 0; iset < polygon->getPolySetNumber(); iset++)
  {
    const PolySet& polyset = polygon->getPolySet(iset);
    pldist = pldist_manage(1, NULL, 2, polyset.getNVertices());

    // Loop on the samples

    for (int iech = 0; iech < nech; iech++)
    {
      if (!db->isActive(iech)) continue;
      distance_point_to_polyline(db->getCoordinate(iech, 0), db->getCoordinate(iech, 1),
                                 polyset.getNVertices(), polyset.getX().data(),
                                 polyset.getY().data(), pldist);
      distloc = pldist->dist;
      if (FFFF(distloc)) continue;
      distmin = db->getArray(iech, iptr);
      if (FFFF(distmin))
        distmin = distloc;
      else
        distmin = MIN(distmin, distloc);
      if (!FFFF(dmax) && distmin > dmax) distmin = dmax;
      db->setArray(iech, iptr, distmin);
    }
    pldist = pldist_manage(-1, pldist, 2, polyset.getNVertices());
  }

  // Calculate the extrema

  if (scale != 0 || polin != 0)
  {
    distmin = 1.e30;
    distmax = 0.;
    for (int iech = 0; iech < nech; iech++)
    {
      if (!db->isActive(iech)) continue;
      distloc = db->getArray(iech, iptr);
      if (FFFF(distloc)) continue;
      if (polin != 0)
      {
        inside = polygon_inside(db->getCoordinate(iech, 0), db->getCoordinate(iech, 1),
                                TEST, 0, polygon);
        if (polin > 0)
        {
          if (!inside) continue;
        }
        else
        {
          if (inside) continue;
        }
      }
      if (distloc > distmax) distmax = distloc;
      if (distloc < distmin) distmin = distloc;
    }
  }

  // Scaling option

  if (scale != 0)
  {
    if (scale > 0)
    {
      for (int iech = 0; iech < nech; iech++)
      {
        if (!db->isActive(iech)) continue;
        distloc = db->getArray(iech, iptr);
        if (FFFF(distloc)) continue;
        value = (distloc - distmin) / (distmax - distmin);
        db->setArray(iech, iptr, value);
      }
    }
    else
    {
      for (int iech = 0; iech < nech; iech++)
      {
        if (!db->isActive(iech)) continue;
        distloc = db->getArray(iech, iptr);
        if (FFFF(distloc)) continue;
        value = (distloc - distmax) / (distmin - distmax);
        db->setArray(iech, iptr, value);
      }
    }
    distmin = 0.;
    distmax = 1.;
  }

  // Check if the sample belongs to the polygon or not

  if (polin != 0)
  {
    valtest = TEST;
    if (ABS(polin) == 2) valtest = distmin;
    if (ABS(polin) == 3) valtest = distmax;
    for (int iech = 0; iech < nech; iech++)
    {
      inside = polygon_inside(db->getCoordinate(iech, 0), db->getCoordinate(iech, 1),
                              TEST, 0, polygon);
      if (polin > 0 && !inside) db->setArray(iech, iptr, valtest);
      if (polin < 0 && inside)  db->setArray(iech, iptr, valtest);
    }
  }

  return (0);
}

/*****************************************************************************/
/*!
 **  Create a new Data Base with points generated at random
 **
 ** \return  Pointer for the new Db structure
 **
 ** \param[in]  mode        Type of Point generation
 **                         0: Poisson
 ** \param[in]  verbose     Verbose option
 ** \param[in]  ndim        Space dimension
 ** \param[in]  seed        Seed for the random number generator
 ** \param[in]  density     Average Poisson intensity
 ** \param[in]  range       Repulsion range (mode=2)
 ** \param[in]  beta        Bending coefficient (mode=2)
 ** \param[in]  dbgrid      Descriptor of the Db grid parameters (mode=1)
 ** \param[in]  origin      Vector of field origin
 ** \param[in]  extend      Vector of field extends
 **
 *****************************************************************************/
GEOSLIB_API Db *db_point_init(int mode,
                              int verbose,
                              int ndim,
                              int seed,
                              double density,
                              double range,
                              double beta,
                              Db *dbgrid,
                              const VectorDouble& origin,
                              const VectorDouble& extend)
{
  VectorDouble tab;
  int count;
  Db* db;
  String string;
  static int flag_add_rank = 1;

  // Dispatch

  switch (mode)
  {
    case 0:
      tab = st_point_init_poisson(verbose, ndim, seed, density, origin, extend,
                                  &count);
      break;

    case 1:
      tab = st_point_init_poisreg(verbose, seed, dbgrid, &count);
      break;

    case 2:
      tab = st_point_init_poisthin(verbose, ndim, seed, density, range, beta,
                                   origin, extend, &count);
  }

  /* Allocate the main structure */

  db = db_create_point(count, ndim, 0, flag_add_rank, tab);

  /* Set the locators */

  for (int idim = 0; idim < ndim; idim++)
  {
    string = getLocatorName(LOC_X, idim);
    db_name_set(db, idim + flag_add_rank, string);
    db->setLocatorByAttribute(idim + flag_add_rank, LOC_X, idim);
  }

  return (db);
}

/****************************************************************************/
/*!
 **  Find the interval (among vector X) to which the sample (x) belongs
 **
 ** \return  Rank of the interval containing the target
 **
 ** \param[in]  x       Target coordinate
 ** \param[in]  ndef    Number of defined samples
 ** \param[in]  X       Vector of coordinate of valued samples
 **
 ** \remarks If x < min(X) or x > max(X) the returned index is -1
 ** \remarks Array X is assumed to be increasingly ordered
 **
 *****************************************************************************/
static int st_find_interval(double x, int ndef, double *X)
{
  if (ndef < 2) return -1;
  if (x < X[0] || x > X[ndef - 1]) return -1;

  // The target belongs to an interval
  for (int k = 0; k < ndef - 1; k++)
    if (x >= X[k] && x < X[k + 1]) return k;

  // The target matches the upper bound of the last interval
  if (x == X[ndef - 1]) return ndef - 1;
  return -1;
}

/****************************************************************************/
/*!
 **  Fill an incomplete 1-D grid by linear interpolation from a set of
 **  valued samples
 **
 ** \return  Error returned code
 **
 ** \param[in]  dbgrid  Db grid structure
 ** \param[in]  ivar    Rank of the variable to be filled
 ** \param[in]  ndef    Number of defined samples
 ** \param[in]  X       Vector of coordinate of valued samples
 ** \param[in]  Y       Vector of values of valued samples
 **
 *****************************************************************************/
static void st_grid1D_interpolate_linear(Db* dbgrid,
                                         int ivar,
                                         int ndef,
                                         double* X,
                                         double* Y)
{
  int nech = dbgrid->getSampleNumber();

  // Loop on the grid nodes

  for (int iech = 0; iech < nech; iech++)
  {
    if (!dbgrid->isActive(iech)) continue;
    double x = dbgrid->getCoordinate(iech, 0);
    int k = st_find_interval(x, ndef, X);
    if (k < 0) continue;
    double y = Y[k] + (Y[k + 1] - Y[k]) * (x - X[k]) / (X[k + 1] - X[k]);
    dbgrid->setVariable(iech, ivar, y);
  }
}

/****************************************************************************/
/*!
 **  Fill an incomplete 1-D grid by spline interpolation from a set of
 **  valued samples
 **
 ** \return  Error returned code
 **
 ** \param[in]  dbgrid  Db grid structure
 ** \param[in]  ivar    Rank of the variable to be filled
 ** \param[in]  ndef    Number of defined samples
 ** \param[in]  X       Vector of coordinate of valued samples
 ** \param[in]  Y       Vector of values of valued samples
 **
 *****************************************************************************/
static int st_grid1D_interpolate_spline(Db *dbgrid,
                                        int ivar,
                                        int ndef,
                                        double* X,
                                        double* Y)
{
  VectorDouble h, F, R, M, C, Cp;
  int nech = dbgrid->getSampleNumber();

  // Preliminary calculations

  int n = ndef;
  int nm1 = n - 1;

  h.resize(nm1);
  for (int i = 0; i < nm1; i++)
    h[i] = X[i + 1] - X[i];

  F.resize(n);
  for (int i = 1; i < nm1; i++)
    F[i] = (Y[i + 1] - Y[i]) / h[i] - (Y[i] - Y[i - 1]) / h[i - 1];
  F[0] = 0;
  F[nm1] = 0;

  R.resize(n * n, 0);
  R(0,0) = 1;
  R(nm1,nm1) = 1.;
  for (int i = 1; i < nm1; i++)
  {
    R(i,i) = (h[i - 1] + h[i]) / 3;
    R(i,i+1) = h[i] / 6;
    R(i,i-1) = h[i - 1] / 6;
  }

  M.resize(n, 0);
  if (matrix_invert(R.data(), n, -1)) return 1;
  matrix_product(n, n, 1, R.data(), F.data(), M.data());

  C.resize(nm1, 0);
  Cp.resize(nm1, 0);
  for (int i = 0; i < nm1; i++)
  {
    C[i] = (Y[i + 1] - Y[i]) / h[i] - (M[i + 1] - M[i]) * h[i] / 6;
    Cp[i] = Y[i] - M[i] * h[i] * h[i] / 6;
  }

  // Loop on the grid nodes

  for (int iech = 0; iech < nech; iech++)
  {
    double y = TEST;
    if (dbgrid->isActive(iech))
    {
      double x = dbgrid->getCoordinate(iech, 0);
      int k = st_find_interval(x, ndef, X);
      if (k >= 0)
      {
        double d1 = X[k + 1] - x;
        double d2 = x - X[k];
        double h6 = 6 * h[k];
        y = (M[k] * pow(d1, 3) + M[k + 1] * pow(d2, 3)) / h6 + C[k] * d2
            + Cp[k];
      }
    }
    dbgrid->setVariable(iech, ivar, y);
  }
  return 0;
}

/****************************************************************************/
/*!
 **  Fill an incomplete 1-D grid
 **
 ** \return  Error returned code
 **
 ** \param[in]  dbgrid  Db grid structure
 ** \param[in]  mode    Type of interpolation
 ** \li                  0 : Linear interpolation
 ** \li                  1 : Cubic Spline
 ** \param[in]  seed    Seed used for the random number generation
 ** \param[in]  namconv Naming convention
 **
 *****************************************************************************/
GEOSLIB_API int db_grid1D_fill(Db *dbgrid,
                               int mode,
                               int seed,
                               NamingConvention namconv)
{
  /* Preliminary checks */

  if (!is_grid(dbgrid))
  {
    messerr("This function is limited to Grid Db");
    return (1);
  }
  int ndim = dbgrid->getNDim();
  if (ndim != 1)
  {
    messerr("This function is limited to 1-D space");
    return (1);
  }
  if (mode < 0 || mode > 1)
  {
    messerr("The argument 'mode' should lie between 0 and 1");
    return (1);
  }
  int nvar = dbgrid->getVariableNumber();
  if (nvar <= 0)
  {
    messerr("You must have at least one Z-locator defined");
    return 1;
  }

  // Add the variables (they must be defined as LOC_Z) for following functions

  int iatt_out = dbgrid->addFields(nvar);
  for (int ivar = 0; ivar < nvar; ivar++)
  {
    int iatt_in = dbgrid->getAttribute(LOC_Z, ivar);
    dbgrid->duplicateColumnByAttribute(iatt_in, iatt_out + ivar);
  }

  int nech = dbgrid->getSampleNumber();

  /* Core allocation */

  law_set_random_seed(seed);
  VectorDouble X(nech, 0);
  VectorDouble Y(nech, 0);

  /* Loop on the variables to be filled */

  for (int ivar = 0; ivar < nvar; ivar++)
  {
    // Copy the input variable

    // Look for the defined values

    int ndef = 0;
    for (int iech = 0; iech < nech; iech++)
    {
      if (!dbgrid->isActive(iech)) continue;
      double value = dbgrid->getVariable(iech, ivar);
      if (FFFF(value)) continue;
      X[ndef] = dbgrid->getCoordinate(iech, 0);
      Y[ndef] = value;
      ndef++;
    }

    // Perform the 1-D interpolation

    if (mode == 0)
      st_grid1D_interpolate_linear(dbgrid, ivar, ndef, X.data(), Y.data());
    else
    {
      if (st_grid1D_interpolate_spline(dbgrid, ivar, ndef, X.data(), Y.data()))
        return 1;
    }
  }

  /* Set the error return code */

  namconv.setNamesAndLocators(dbgrid, LOC_Z, -1, dbgrid, iatt_out);

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
 ** \param[in]  ldmax      Type of distance for calculating maximum distance
 **                        1 for L1 and 2 for L2 distance
 ** \param[in]  dmax       Array of maximum distances (optional)
 ** \param[in]  flag_fill  Filling option
 ** \param[in]  flag_inter Interpolation
 **
 *****************************************************************************/
static int st_migrate(Db* db1,
                      Db* db2,
                      int iatt1,
                      int iatt2,
                      int ldmax,
                      const VectorDouble& dmax,
                      int flag_fill,
                      int flag_inter)
{
  int size = db2->getSampleNumber();
  VectorDouble tab(size, TEST);

  if (db2->isGrid())
  {
    // To Grid
    if (db1->isGrid())
    {
      // Grid to Grid
      if (flag_fill)
      {
        // Grid to Grid (flag_fill = TRUE)
        if (st_expand_grid_to_grid(db1, db2, iatt1, ldmax, dmax, tab))
          return 1;
      }
      else
      {
        // Grid to Grid (flag_fill = FALSE)
        if (st_migrate_grid_to_grid(db1, db2, iatt1, ldmax, dmax, tab))
          return 1;
      }
    }
    else
    {
      // Point to Grid
      if (flag_fill)
      {
        // Point to Grid (flag_fill = TRUE)
        if (expand_point_to_grid(db1, db2, iatt1, -1, -1, -1, -1, -1, 0, ldmax,
                                 dmax, tab)) return 1;
      }
      else
      {
        // Point to Grid (flag_fill = FALSE)
        if (st_migrate_point_to_grid(db1, db2, iatt1, ldmax, dmax, tab))
          return 1;
      }
    }
  }
  else if (db1->isGrid())
  {
    // Grid to Point
    if (flag_inter)
    {
      // Grid to Point (flag_inter = TRUE)
      if (st_interpolate_grid_to_point(db1, db2, iatt1, ldmax, dmax,
                                       tab)) return 1;
    }
    else
    {
      // Grid to Point (flag_inter = FALSE)
      if (st_migrate_grid_to_point(db1, db2, iatt1, ldmax, dmax, tab))
        return 1;
    }
  }
  else
  {
    // Point to Point
    if (st_expand_point_to_point(db1, db2, iatt1, ldmax, dmax, tab))
      return 1;
  }

  // Store the resulting array in the output Db

  db2->setFieldByAttribute(tab, iatt2);
  return 0;
}

/*****************************************************************************/
/*!
 **  Migrates a variable from one Db to another one
 **
 ** \return  Error return code
 **
 ** \param[in]  db1        Descriptor of the input Db
 ** \param[in]  db2        Descriptor of the output Db
 ** \param[in]  atts_arg   Array of attributes to be migrated
 ** \param[in]  ldmax      Type of distance for calculating maximum distance
 **                        1 for L1 and 2 for L2 distance
 ** \param[in]  dmax       Array of maximum distances (optional)
 ** \param[in]  flag_fill  Filling option
 ** \param[in]  flag_inter Interpolation
 ** \param[in]  namconv    Naming Convention
 **
 *****************************************************************************/
GEOSLIB_API int migrateByAttribute(Db* db1,
                                   Db* db2,
                                   const VectorInt& atts_arg,
                                   int ldmax,
                                   const VectorDouble& dmax,
                                   int flag_fill,
                                   int flag_inter,
                                   NamingConvention namconv)
{
  // CDesignate the input variables

  VectorInt atts = atts_arg;
  int ncol = static_cast<int> (atts.size());
  if (atts.empty())
  {
    atts = db1->getAttributes();
    ncol = static_cast<int> (atts.size());
  }

  // Create the output variables

  int iatt0 = db2->addFields(ncol, TEST);

  // Loop on the different variables obtained by various migrations

  for (int i = 0; i < ncol; i++)
  {
    int iatt1 = atts[i];
    int iatt2 = iatt0 + i;
    if (st_migrate(db1, db2, iatt1, iatt2, ldmax, dmax, flag_fill, flag_inter))
      return 1;
  }

  // Set the output variable names and locators
  namconv.setNamesAndLocators(db1, atts, db2, iatt0);
  return 0;
}

/*****************************************************************************/
/*!
 **  Migrates a variable from one Db to another one
 **
 ** \return  Error return code
 **
 ** \param[in]  db1        Descriptor of the input Db
 ** \param[in]  db2        Descriptor of the output Db
 ** \param[in]  name       Name of the attribute to be migrated
 ** \param[in]  ldmax      Type of distance for calculating maximum distance
 **                        1 for L1 and 2 for L2 distance
 ** \param[in]  dmax       Array of maximum distances (optional)
 ** \param[in]  flag_fill  Filling option
 ** \param[in]  flag_inter Interpolation
 ** \param[in]  namconv    Naming convention
 **
 *****************************************************************************/
GEOSLIB_API int migrate(Db* db1,
                        Db* db2,
                        const String& name,
                        int ldmax,
                        const VectorDouble& dmax,
                        int flag_fill,
                        int flag_inter,
                        NamingConvention namconv)
{
  VectorInt iatts = db1->ids(name, true);
  if (iatts.empty()) return 1;

  // Create the output variables

  int iatt0 = db2->addFields(1, TEST);

  // Perform the migration

  if (st_migrate(db1, db2, iatts[0], iatt0, ldmax, dmax, flag_fill, flag_inter))
    return 1;

  // Set the output variable names and locators

  namconv.setNamesAndLocators(name, db2, iatt0);
  return 0;
}

/*****************************************************************************/
/*!
 **  Migrates all z-locator variables from one Db to another one
 **
 ** \return  Error return code
 **
 ** \param[in]  db1         Descriptor of the input Db
 ** \param[in]  db2         Descriptor of the output Db
 ** \param[in]  locatorType Locator Type
 ** \param[in]  ldmax       Type of distance for calculating maximum distance
 **                         1 for L1 and 2 for L2 distance
 ** \param[in]  dmax        Array of maximum distances (optional)
 ** \param[in]  flag_fill   Filling option
 ** \param[in]  flag_inter  Interpolation
 ** \param[in]  namconv     Naming convention
 **
 *****************************************************************************/
GEOSLIB_API int migrateByLocator(Db* db1,
                                 Db* db2,
                                 ENUM_LOCS locatorType,
                                 int ldmax,
                                 const VectorDouble& dmax,
                                 int flag_fill,
                                 int flag_inter,
                                 NamingConvention namconv)
{
  VectorString names = db1->getNames(locatorType);
  int natt = static_cast<int> (names.size());
  if (natt <= 0) return 0;

  // Create the output variables

  int iatt0 = db2->addFields(natt, TEST);

  // Loop on the different variables obtained by various migrations

  for (int i = 0; i < natt; i++)
  {
    int iatt1 = db1->getAttribute(names[i]);
    int iatt2 = iatt0 + i;
    if (st_migrate(db1, db2, iatt1, iatt2, ldmax, dmax, flag_fill, flag_inter))
      return 1;
  }

  // Set the output variable names and locators
  namconv.setNamesAndLocators(db1, locatorType, -1, db2, iatt0);
  return 0;
}

/****************************************************************************/
/*!
**  Standard Kriging
**
** \return  Error return code
**
** \param[in]  dbin        Input Db structure
** \param[in]  dbout       Output Db structure
** \param[in]  model       Model structure
** \param[in]  niter       Number of iterations
** \param[in]  verbose     Verbose flag
** \param[in]  namconv     Naming convention
**
** \remarks The procedure uses the FIRST covariance of the Model
** \remarks to describe the spatial structure
**
*****************************************************************************/
GEOSLIB_API int db_proportion_estimate(Db *dbin,
                                       Db *dbout,
                                       Model *model,
                                       int niter,
                                       bool verbose,
                                       NamingConvention namconv)
{
  VectorVectorInt splits;

  // Preliminary checks

  if (dbin == nullptr)
  {
    messerr("This method requires a 'dbin' argument");
    return 1;
  }
  if (dbout == nullptr)
  {
    messerr("This method requires a 'dbout' argument");
    return 1;
  }
  if (model == nullptr)
  {
    messerr("This method requires a 'model' argument");
    return 1;
  }
  if (dbin->getVariableNumber() != 1)
  {
    messerr("The argument 'dbin' should have a single variable");
    return 1;
  }

  // Define the environment

  MeshETurbo mesh = MeshETurbo(*dbout);
  ShiftOpCs S = ShiftOpCs(&mesh, model, dbout);
  PrecisionOp Qprop =  PrecisionOp(&S, model->getCova(0),  POPT_ONE);
  ProjMatrix Aproj =  ProjMatrix(dbin, &mesh);

  // Invoke the calculation

  VectorDouble propGlob = dbStatisticsFacies(dbin);
  int ncat = propGlob.size();
  OptimCostColored Oc =  OptimCostColored(ncat,&Qprop,&Aproj);
  Oc.setCGParams(200,1.e-10);

  VectorDouble facies = dbin->getFieldByLocator(LOC_Z);
  VectorVectorDouble props = Oc.minimize(facies,splits,propGlob,verbose,niter);

  // Loading the resulting results in the output 'dbout'

  int iptr0 = -1;
  for (int i = 0; i < ncat; i++)
  {
    int iptr = dbout->addFields(props[i]);
    if (i == 0) iptr0 = iptr;
    namconv.setNamesAndLocators(
        nullptr, LOC_UNKNOWN, -1, dbout, iptr,
        concatenateStrings("-", intToString(i + 1)));
  }
  namconv.setLocators(dbout, iptr0, 1, ncat);

  return 0;
}
