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
#include "geoslib_old_f.h"

#include "Polygon/Polygons.hpp"
#include "Basic/Utilities.hpp"
#include "Basic/Grid.hpp"
#include "Db/Db.hpp"
#include "Db/DbGrid.hpp"
#include "Core/Memory.hpp"

#include <math.h>

/****************************************************************************/
/*!
 **  Checks if Space Dimension of the first Db is at least equal
 **  to the Space dimension of the second Db
 **
 ** \return  Error return code (1 for correct and 0 otherwise).
 ** \return  If an error is found, a message is issued
 **
 ** \param[in]  db1    first Db descriptor (usually dbin)
 ** \param[in]  db2    second Db descriptor (usually dbout)
 **
 ** \remarks This function checks that the we can perform an operation
 ** \remarks (i.e. migration) of a 2D data set onto a 3D one.
 ** \remarks The opposite is not correct.
 **
 *****************************************************************************/
int compat_NDIM(Db *db1, Db *db2)
{
  if (db1->getNDim() <= db2->getNDim()) return (1);
  messerr("The Space Dimension of the First Db (%d)", db1->getNDim());
  messerr("must not be smaller than the Space Dimension of the Second Db",
          db2->getNDim());
  return (0);
}

/****************************************************************************/
/*!
 **  Print the Summary of the Grid structure
 **
 ** \param[in]  db  Pointer to the Db structure (organized as a grid)
 **
 *****************************************************************************/
void db_grid_print(Db *db)
{
  if (db->isGrid()) message(db->toString().c_str());
}

/****************************************************************************/
/*!
 **  Returns the number of items for a given locator in the Db
 **
 ** \return  Number of items
 **
 ** \param[in]  db     Db descriptor
 ** \param[in]  locatorType Rank of the pointer (ELoc)
 **
 *****************************************************************************/
int get_LOCATOR_NITEM(const Db *db, const ELoc& locatorType)
{
  if (db == nullptr) return (0);
  if (db->isGrid() && locatorType == ELoc::X)
    return (db->getNDim());
  return (db->getFromLocatorNumber(locatorType));
}

/****************************************************************************/
/*!
 **  Frees the array for storing a sample
 **
 ** \return  A pointer to the array to be freed
 **
 ** \param[in]  tab  Sample array to be freed
 **
 *****************************************************************************/
double* db_sample_free(double *tab)

{
  tab = (double*) mem_free((char* ) tab);
  return (tab);
}

/****************************************************************************/
/*!
 **  Allocates the array for storing a sample
 **
 ** \return  A pointer to the allocated array
 **
 ** \param[in]  db      Db descriptor
 ** \param[in]  locatorType  vector type (ELoc)
 **
 ** \remark  The allocated array must be freed using db_sample_free()
 ** \remark  A fatal error occurs if the core allocation fails.
 **
 *****************************************************************************/
double* db_sample_alloc(const Db *db, const ELoc& locatorType)
{
  double *tab;
  int size;

  /* Initializations */

  tab = nullptr;
  size = get_LOCATOR_NITEM(db, locatorType);

  /* In the case of a grid, there may be no actual data vector */
  if (locatorType == ELoc::X && db->isGrid()) size = db->getNDim();
  if (size > 0) tab = (double*) mem_alloc(sizeof(double) * size, 1);
  return (tab);
}

/****************************************************************************/
/*!
 **  Calculates the distance between two points
 **
 ** \return  The calculated distance or TEST if one coordinate not defined
 **
 ** \param[in]  db1    Db structure for the first sample
 ** \param[in]  db2    Db structure for the second sample
 ** \param[in]  iech1  rank of the first sample
 ** \param[in]  iech2  rank of the second sample
 **
 ** \param[out] dist_vect  If the output vector is provided.
 **                        Returns the distance as a vector
 **
 ** \remark: If both Data Bases do not share the same Space Dimension
 ** \remark: the test is performed on their minimum.
 ** \remark: The returned array 'vect' must be dimension to that value
 **
 *****************************************************************************/
double distance_inter(const Db *db1,
                      const Db *db2,
                      int iech1,
                      int iech2,
                      double *dist_vect)
{
  double v1, v2, *tab1, *tab2;
  int idim, ndim;

  ndim = MIN(db1->getNDim(), db2->getNDim());
  ut_distance_allocated(ndim, &tab1, &tab2);

  for (idim = 0; idim < ndim; idim++)
  {
    v1 = db1->getCoordinate(iech1, idim);
    v2 = db2->getCoordinate(iech2, idim);
    if (FFFF(v1) || FFFF(v2)) return (TEST);
    tab1[idim] = v1;
    tab2[idim] = v2;
    if (dist_vect != nullptr) dist_vect[idim] = v1 - v2;
  }
  return (ut_distance(ndim, tab1, tab2));
}

/****************************************************************************/
/*!
 **  Calculates the distance between two points in the same Db
 **
 ** \return  The calculated distance or TEST if one coordinate not defined
 **
 ** \param[in]  db     Db structure
 ** \param[in]  iech1  rank of the first sample
 ** \param[in]  iech2  rank of the second sample
 **
 ** \param[out] dist_vect  If the output vector is provided.
 **                        Returns the distance as a vector
 **
 *****************************************************************************/
double distance_intra(const Db *db, int iech1, int iech2, double *dist_vect)
{
  double v1, v2, *tab1, *tab2;
  int idim, ndim;

  ndim = db->getNDim();
  ut_distance_allocated(ndim, &tab1, &tab2);

  for (idim = 0; idim < ndim; idim++)
  {
    v1 = db->getCoordinate(iech1, idim);
    v2 = db->getCoordinate(iech2, idim);
    if (FFFF(v1) || FFFF(v2)) return (TEST);
    tab1[idim] = v1;
    tab2[idim] = v2;
    if (dist_vect != nullptr) dist_vect[idim] = v1 - v2;
  }
  return (ut_distance(ndim, tab1, tab2));
}

/****************************************************************************/
/*!
 **  Calculates the distance between two points in the same grid Db
 **
 ** \return  The calculated distance
 **
 ** \param[in]  db          Db structure
 ** \param[in]  flag_moins1 1 to use the cell center
 ** \param[in]  iech1       rank of the first sample
 ** \param[in]  iech2       rank of the second sample
 **
 ** \param[out] dist_vect   If the output vector is provided.
 **                         Returns the distance as a vector
 **
 *****************************************************************************/
double distance_grid(DbGrid *db,
                     int flag_moins1,
                     int iech1,
                     int iech2,
                     double *dist_vect)
{
  int ndim = db->getNDim();
  VectorInt iwork1(ndim);
  VectorInt iwork2(ndim);

  /* Locate the two samples on the grid */

  if (iech1 == iech2)
  {
    if (dist_vect != nullptr) for (int idim = 0; idim < db->getNDim(); idim++)
      dist_vect[idim] = 0.;
    return (0.);
  }

  /* Find the grid indices */

  db->rankToIndice(iech1, iwork1);
  db->rankToIndice(iech2, iwork2);

  /* Calculate the distance */

  double dist = 0.;
  for (int idim = 0; idim < db->getNDim(); idim++)
  {
    int number = ABS(iwork1[idim] - iwork2[idim]);
    if (flag_moins1 != 0 && number > 1) number--;
    double delta = number * db->getDX(idim);
    if (dist_vect != nullptr) dist_vect[idim] = delta;
    dist += delta * delta;
  }

  return (sqrt(dist));
}

/****************************************************************************/
/*!
 **  Print a sample
 **
 ** \param[in]  db        Db structure
 ** \param[in]  iech      Rank of the sample
 ** \param[in]  flag_ndim 1 if the coordinates must be printed
 ** \param[in]  flag_nvar 1 if the variables must be printed
 ** \param[in]  flag_nerr 1 if the error measurement variance must be printed
 **
 *****************************************************************************/
void db_sample_print(Db *db,
                     int iech,
                     int flag_ndim,
                     int flag_nvar,
                     int flag_nerr)
{
  message("Sample #%d (from %d)\n", iech + 1, db->getSampleNumber());
  if (flag_ndim != 0)
  {
    for (int idim = 0; idim < db->getNDim(); idim++)
    {
      double value = db->getCoordinate(iech, idim);
      if (FFFF(value))
        message("Coordinate #%d = NA\n", idim + 1);
      else
        message("Coordinate #%d = %lf\n", idim + 1,
                db->getCoordinate(iech, idim));
    }
  }
  if (flag_nvar != 0)
  {
    for (int ivar = 0; ivar < db->getLocNumber(ELoc::Z); ivar++)
    {
      double value = db->getZVariable(iech, ivar);
      if (FFFF(value))
        message("Variable   #%d = NA\n", ivar + 1);
      else
        message("Variable   #%d = %lf\n", ivar + 1, db->getZVariable(iech, ivar));
    }
  }
  if (flag_nerr != 0)
  {
    for (int ierr = 0; ierr < db->getLocNumber(ELoc::V); ierr++)
    {
      double value = db->getLocVariable(ELoc::V,iech, ierr);
      if (FFFF(value))
        message("Variance   #%d = NA\n", ierr + 1);
      else
        message("Variance   #%d = %lf\n", ierr + 1,
                db->getLocVariable(ELoc::V,iech, ierr));
    }
  }
  if (db->hasLocVariable(ELoc::C))
  {
    double value = db->getLocVariable(ELoc::C,iech,0);
    if (FFFF(value))
      message("Code          = NA\n");
    else
      message("Code          = %d\n", (int) value);
  }
}

/****************************************************************************/
/*!
 **  Returns the center of a data set
 **
 ** \param[in]  db    Db structure
 **
 ** \param[out] center  Array containing the coordinates of the center
 **                     (Dimension = get_NDIM(db))
 **
 *****************************************************************************/
int db_center(Db *db, double *center)
{
  VectorDouble wgt;
  VectorDouble sel;

  /* Initializations */

  if (db->hasLocVariable(ELoc::SEL))
    sel = db->getColumnByLocator(ELoc::SEL);
  if (db->hasLocVariable(ELoc::W))
    wgt = db->getColumnByLocator(ELoc::W);

  /* Loop on the space dimension */

  for (int idim = 0; idim < db->getNDim(); idim++)
  {
    VectorDouble tab = db->getCoordinates(idim);
    StatResults stats = ut_statistics(db->getSampleNumber(), tab.data(), sel.data(), wgt.data());
    center[idim] = stats.mean;
  }

  return (0);
}

/****************************************************************************/
/*!
 **  Returns the range for a variable
 **
 ** \return  Error return code
 **
 ** \param[in]  db    Db structure
 ** \param[in]  iatt  Rank of the target attribute
 **
 ** \param[out]  mini   Minimum
 ** \param[out]  maxi   Maximum
 ** \param[out]  delta  Extension
 **
 *****************************************************************************/
int db_attribute_range(const Db *db,
                       int iatt,
                       double *mini,
                       double *maxi,
                       double *delta)
{
  StatResults stats;

  /* Initializations */

  *mini = TEST;
  *maxi = TEST;
  *delta = TEST;

  /* Load the variable */

  VectorDouble tab = db->getColumnByUID(iatt);
  VectorDouble sel;
  if (db->hasLocVariable(ELoc::SEL))
    sel = db->getColumnByLocator(ELoc::SEL);

  /* Calculate the statistics */

  stats  = ut_statistics(db->getSampleNumber(), tab.data(), sel.data());
  *mini  = stats.mini;
  *maxi  = stats.maxi;
  *delta = stats.delta;
  return 0;
}

/****************************************************************************/
/*!
 **  Define the coordinates in a Grid structure
 **
 ** \return Error return code
 **
 ** \param[in]  db        The Db structure
 **
 ** \remark This function considers the grid characteristics and updates
 ** \remark the locators dedicated to coordinates
 ** \remark This makes sense when a new grid is generated or when the grid
 ** \remark characteristics have changed
 **
 *****************************************************************************/
int db_grid_define_coordinates(DbGrid *db)

{
  if (db == nullptr) return (0);
  int ndim = db->getNDim();
  int nech = db->getSampleNumber();
  VectorInt ntab(ndim, 0);
  VectorDouble coor(ndim);
  VectorDouble cbis(ndim);

  /* Fill the coordinates */

  for (int iech = 0; iech < nech; iech++)
  {

    /* Set the array of coordinates */

    for (int idim = 0; idim < ndim; idim++)
      coor[idim] = db->getDX(idim) * ntab[idim];

    /* Process the rotation (optional) */

    if (db->isGridRotated())
    {
      db->getGrid().getRotation().rotateDirect(coor, cbis);
      coor = cbis;
    }

    for (int idim = 0; idim < ndim; idim++)
    {
      coor[idim] += db->getX0(idim);
      db->setFromLocator(ELoc::X, iech, idim, coor[idim]);
    }

    /* Calculate the new coordinates */

    int nval = 1;
    for (int idim = 0; idim < ndim; idim++)
    {
      if ((iech + 1) % nval == 0)
      {
        ntab[idim]++;
        if (ntab[idim] == db->getNX(idim)) ntab[idim] = 0;
      }
      nval *= db->getNX(idim);
    }
  }
  return (0);
}

/****************************************************************************/
/*!
 **  Update the Data Base for Kriging with Gradient components
 **  This update is limited to the 2-D Monovariate case
 **
 ** \return  Error return code
 **
 ** \param[in]  db Input/output Db
 **
 ** \remark  In the case of Kriging with Gradient, the gradient components
 ** \remark  are transformed into additional variables
 **
 *****************************************************************************/
int db_gradient_update(Db *db)

{
  int ndim = db->getNDim();
  int ngrad = db->getLocNumber(ELoc::G);
  int nvar = db->getLocNumber(ELoc::Z);

  /* Preliminary checks */

  if (nvar != 1)
  {
    messerr("Kriging with Gradients if limited to the Monovariate case");
    return (1);
  }
  if (ndim != ngrad)
  {
    messerr(
        "The number of Gradient components (%d) must coincide with Space dimension (%d)",
        ngrad, ndim);
    return (1);
  }

  // Turn the locators from Gradients into Variables
  db->switchLocator(ELoc::G, ELoc::Z);

  return (0);
}

/*****************************************************************************/
/*!
 **  Select a 2-D grid from a Grid Db
 **
 ** \return  Error returned code
 **
 ** \param[in]  ndim   Space dimension of the Db
 ** \param[in]  nx     Number of meshes of the grid
 ** \param[in]  ref    Array giving the indices of the reference planes
 ** \li                Grid indices are counted starting from 1
 ** \li                -1 for the directions of interest
 ** \param[in]  tabin  Input array
 **
 ** \param[out] tabout Output array
 **
 *****************************************************************************/
int db_selref(int ndim, const int* nx, const int* ref, const double* tabin, double* tabout)
{
  int *rank, *ind1, idim, jdim, ntotal, nval, lec, ecr, iech, ival, error, neff_ndim;

  /* Initializations */

  error = 1;
  rank = ind1 = nullptr;

  /* Core allocation */

  rank = (int*) mem_alloc(sizeof(int) * ndim, 0);
  if (rank == nullptr) goto label_end;
  ind1 = (int*) mem_alloc(sizeof(int) * ndim, 0);
  if (ind1 == nullptr) goto label_end;

  /* Set the indices */

  for (idim = 0; idim < ndim; idim++)
  {
    rank[idim] = -1;
    for (jdim = 0; jdim < ndim; jdim++)
      if (ref[jdim] == -(idim + 1)) rank[idim] = jdim;
  }

  /* Calculate the effective dimension of the resulting grid */

  neff_ndim = ndim;
  for (idim = ndim - 1; idim >= 0; idim--)
    if (rank[idim] < 0) neff_ndim--;

  /* Count the total number of samples */

  ntotal = 1;
  for (idim = 0; idim < ndim; idim++)
    ntotal *= nx[idim];

  bool skip;
  for (lec = ecr = 0; lec < ntotal; lec++)
  {
    nval = ntotal;
    iech = lec;
    skip = false;
    for (idim = ndim - 1; idim >= 0 && ! skip; idim--)
    {
      nval /= nx[idim];
      ival = iech / nval;
      if (ref[idim] > 0 && ival != ref[idim] - 1) skip = true;
      ind1[idim] = ival;
      iech -= ival * nval;
    }
    if (skip) continue;

    /* Find the writing index */

    idim = neff_ndim - 1;
    jdim = rank[idim];
    ecr = ind1[jdim];
    for (idim = neff_ndim - 2; idim >= 0; idim--)
    {
      jdim = rank[idim];
      ecr = ecr * nx[jdim] + ind1[jdim];
    }
    tabout[ecr] = tabin[lec];
  }

  /* Set the error code */

  error = 0;

  /* Core deallocation */

  label_end:
  mem_free((char* ) rank);
  mem_free((char* ) ind1);
  return (error);
}

/****************************************************************************/
/*!
 **  Add and initiate several attributes corresponding to a given locator
 **
 ** \return  Error return code
 **
 ** \param[in]  db      Db structure
 ** \param[in]  locatorType  Rank of the Pointer (ELoc)
 ** \param[in]  number  Number of locators to be defined
 ** \param[in]  r_tem   Rank of the first item in the pointer
 ** \param[in]  valinit Value to be used for initialization
 **
 ** \param[out] iptr    Rank of the first new attribute
 **
 *****************************************************************************/
int db_locator_attribute_add(Db *db,
                             const ELoc& locatorType,
                             int number,
                             int r_tem,
                             double valinit,
                             int *iptr)
{
  (*iptr) = db->addColumnsByConstant(number, valinit);
  if ((*iptr) < 0) return (1);
  db->setLocatorsByUID(number, (*iptr), locatorType, r_tem);

  /* Set the default names to the newly created variables */

  for (int i = 0; i < number; i++)
  {
    String string = getLocatorName(locatorType, r_tem + i);
    db->setNameByUID((*iptr) + i, string);
  }
  return (0);
}

/****************************************************************************/
/*!
 **  Copy a set of variables from a grid Db to another grid Db
 **
 ** \return  Error return code
 **
 ** \param[in]  db1     Input Grid Db structure
 ** \param[in]  db2     Output Grid Db structure
 ** \param[in]  ind1    Vector of indices in the first Db
 ** \param[in]  ind2    Vector of indices in the second Db
 ** \param[in]  ncol    Number of variables to be copied
 ** \param[in]  cols    Array of input variable columns
 **
 *****************************************************************************/
int db_grid_copy(DbGrid *db1,
                 DbGrid *db2,
                 const int *ind1,
                 const int *ind2,
                 int ncol,
                 int *cols)
{
  /* Set the constant indices in the input Grid Db */

  int ndim1 = db1->getNDim();
  VectorInt iwork1(ndim1);
  for (int idim = 0; idim < ndim1; idim++)
    iwork1[idim] = ind1[idim] - 1;

  /* Set the variable index to 1 only for subsequent test */

  int ndim2 = db2->getNDim();
  for (int idim = 0; idim < ndim2; idim++)
    if (ind2[idim] != 0) iwork1[ind2[idim] - 1] = 1;
  for (int idim = 0; idim < ndim1; idim++)
  {
    if (iwork1[idim] < 0 || iwork1[idim] >= db1->getNX(idim))
    {
      messerr("The index %d of the input Grid Db is not assigned", idim);
      messerr("Copy operation is cancelled");
      return 1;
    }
  }

  /* Add the variables */

  int iptr = db2->addColumnsByConstant(ncol, TEST);

  /* Loop on the output grid Db */

  for (int iech = 0; iech < db2->getSampleNumber(); iech++)
  {

    /* Find the indices of the target grid node */

    db2->rankToIndice(iech, iwork1);
    for (int idim = 0; idim < db2->getNDim(); idim++)
    {
      if (ind2[idim] > 0)
        iwork1[ind2[idim] - 1] = iwork1[idim];
      else
        iwork1[-ind2[idim] - 1] = db2->getNX(idim) - iwork1[idim] - 1;
    }

    /* Restrain the indices to the extension of the first Grid Db */

    for (int idim = 0; idim < db1->getNDim(); idim++)
    {
      int indice = iwork1[idim];
      if (FFFF(indice)) messageAbort("This error should not happen");
      if (indice < 0) indice = 0;
      if (indice >= db1->getNX(idim)) indice = db1->getNX(idim) - 1;
      iwork1[idim] = indice;
    }

    /* Convert into absolute index */

    int jech = db1->indiceToRank(iwork1);

    /* Loop on the variables to be copied */

    for (int icol = 0; icol < ncol; icol++)
      db2->setArray(iech, iptr + icol, db1->getArray(jech, cols[icol]));
  }

  return 0;
}

/****************************************************************************/
/*!
 **  Copy one variable from a grid Db to another grid Db
 **
 ** \return  Error return code
 **
 ** \param[in]  db1     Input Grid Db structure
 ** \param[in]  iatt1   Rank of the attribute in Db1
 ** \param[in]  db2     Output Grid Db structure
 ** \param[in]  iatt2   Rank of the attribute in Db2
 ** \param[in]  mode    1 for dilation; -1 for compression
 ** \param[in]  nshift  Vector of dilation/compression defined in cell number
 **                     along each space direction
 **
 *****************************************************************************/
int db_grid_copy_dilate(DbGrid *db1,
                        int iatt1,
                        DbGrid *db2,
                        int iatt2,
                        int mode,
                        const int *nshift)
{
  double value;

  /* Initializations */

  int ndim  = db1->getNDim();
  VectorInt indg(ndim);

  /* Check that the grids are compatible */

  if (!db1->hasSameDimension(db2)) return 1;
  if (! db1->isGrid() || ! db2->isGrid())
  {
    messerr("The function 'db_grid_copy_dilate' requires two grid Dbs");
    return 1;
  }

  /* Loop on the samples of the second Db */

  for (int iech2 = 0; iech2 < db2->getSampleNumber(); iech2++)
  {
    db2->rankToIndice(iech2, indg);
    for (int idim = 0; idim < ndim; idim++)
      indg[idim] += mode * nshift[idim];
    int iech1 = db1->indiceToRank(indg);

    if (iech1 < 0)
      value = TEST;
    else
      value = db1->getArray(iech1, iatt1);

    db2->setArray(iech2, iatt2, value);
  }
  return 0;
}

/*****************************************************************************/
/*!
 **  Returns the rank of the closest isolated point
 **
 ** \return  Output index
 **
 ** \param[in]  db    descriptor of the Db
 ** \param[in]  coor  array of coordinates of the point
 **
 *****************************************************************************/
int point_to_point(Db *db, const double *coor)
{
  double dist, distmin, delta, x;
  int idim, iech, iechmin;

  /* Loop on the input structure */

  iechmin = 0;
  distmin = 1.e30;
  for (iech = 0; iech < db->getSampleNumber(); iech++)
  {
    if (!db->isActive(iech)) continue;

    dist = 0.;
    for (idim = 0; idim < db->getNDim(); idim++)
    {
      x = db->getCoordinate(iech, idim);
      if (FFFF(x)) continue;
      delta = x - coor[idim];
      dist += delta * delta;
    }

    if (dist < distmin)
    {
      distmin = dist;
      iechmin = iech;
    }
  }

  /* Returning arguments */

  return (iechmin);
}

/*****************************************************************************/
/*!
 **  Converts from point coordinates to nearest grid node indices
 **
 ** \return  Error return code
 ** \return   0 if the point is inside the grid
 ** \return   1 if the point is outside the grid
 ** \return  -1 if one coordinate is undefined
 **
 ** \param[in]  db            descriptor of the grid parameters
 ** \param[in]  coor          array of coordinates of the point
 ** \param[in]  flag_outside  value returned for the point outside the grid
 ** \li                       1 the index is set to the closest grid node
 ** \li                       0 the index is set to -1
 ** \li                      -1 do not correct the index
 **
 ** \param[out] indg          indices of the closest grid node
 **
 *****************************************************************************/
int point_to_grid(const DbGrid *db, const double *coor, int flag_outside, int *indg)
{
  int ndim = db->getNDim();
  VectorDouble work1(ndim);
  VectorDouble work2(ndim);

  /* Check if all coordinates are defined */

  for (int idim = 0; idim < ndim; idim++)
    if (FFFF(coor[idim])) return (-1);

  /* Process the grid rotation (if any) */

  for (int idim = 0; idim < ndim; idim++)
    work1[idim] = coor[idim] - db->getX0(idim);

  db->getGrid().getRotation().rotateInverse(work1, work2);

  /* Calculate the grid indices */

  int out = 0;
  for (int idim = 0; idim < ndim; idim++)
  {
    int ix = (int) (floor(work2[idim] / db->getDX(idim) + 0.5));
    if (ix < 0)
    {
      if (flag_outside > 0)
        ix = 0;
      else if (flag_outside == 0) ix = -1;
      out = 1;
    }
    else if (ix >= db->getNX(idim))
    {
      if (flag_outside > 0)
        ix = db->getNX(idim) - 1;
      else if (flag_outside == 0) ix = -1;
      out = 1;
    }
    indg[idim] = ix;
  }
  return (out);
}

/*****************************************************************************/
/*!
 **  Converts from point coordinates to index of the bench to which it belongs
 **
 ** \return  Error return code
 ** \return   0 if the point is inside the grid
 ** \return   1 if the point is outside the grid
 ** \return  -1 if one coordinate is undefined
 ** \return  -2 if the grid is defined in space less than 3D
 **
 ** \param[in]  db            descriptor of the grid parameters
 ** \param[in]  coor          array of coordinates of the point
 ** \param[in]  flag_outside  value returned for the point outside the grid
 ** \li                       1 the index is set to the closest grid node
 ** \li                       0 the index is set to -1
 ** \li                      -1 do not correct the index
 **
 ** \param[out] indb         index of the bench
 **
 ** \remarks The bench corresponds to the third dimension of the grid provided
 ** \remarks as reference
 **
 *****************************************************************************/
int point_to_bench(const DbGrid *db, double *coor, int flag_outside, int *indb)
{
  int ndim = db->getNDim();
  VectorDouble work1(ndim);
  VectorDouble work2(ndim);
  int idim0 = 2;
  (*indb) = -1;

  /* Check that the grid is defined in 3-D (or more) space */

  if (! db->isGrid() || ndim <= 2) return (-2);

  /* Check if all coordinates are defined */

  for (int idim = 0; idim < ndim; idim++)
    if (FFFF(coor[idim])) return (-1);

  /* Process the grid rotation (if any) */

  for (int idim = 0; idim < ndim; idim++)
    work1[idim] = coor[idim] - db->getX0(idim);

  db->getGrid().getRotation().rotateInverse(work1, work2);

  /* Calculate the bench index */

  int out = 0;

  double z = work2[idim0];
  double dz = db->getDX(idim0);
  int nz = db->getNX(idim0);

  int iz;
  if (dz <= 0.)
    iz = 0;
  else
    iz = (int) (floor(z / dz + 0.5));

  //message("iz = %d & nz = %d\n", iz, nz);

  if (iz < 0)
  {
    if (flag_outside > 0)
      iz = 0;
    else if (flag_outside == 0) iz = -1;
    out = 1;
  }
  else if (iz >= nz)
  {
    if (flag_outside > 0)
      iz = nz - 1;
    else if (flag_outside == 0) iz = -1;
    out = 1;
  }
  (*indb) = iz;

  return (out);
}

/*****************************************************************************/
/*!
 **  Find the index of the output grid file which is the closest
 **  to the sample of the input file
 **
 ** \return  Error return code
 ** \return  >= index of the grid node
 ** \return  -1 if the point is outside the grid
 **
 ** \param[in]  dbin          descriptor of the input file
 ** \param[in]  iech          Index of the data point
 ** \param[in]  flag_outside  value returned for the point outside the grid
 ** \li                       1 the index is set to the closest grid node
 ** \li                       0 the index is set to -1
 ** \li                      -1 do not correct the index
 ** \param[in]  dbout         descriptor of the output grid file
 **
 ** \param[out] coor          Working array (dimension: ndim)
 **
 *****************************************************************************/
int index_point_to_grid(const Db *dbin,
                        int iech,
                        int flag_outside,
                        const DbGrid *dbout,
                        double *coor)
{
  int ndim = dbin->getNDim();
  int nech = dbin->getSampleNumber();
  VectorInt iwork1(ndim);

  /* Get the coordinates of the input sample */

  if (iech < 0 || iech >= nech) return (-1);
  for (int idim = 0; idim < ndim; idim++)
    coor[idim] = dbin->getCoordinate(iech, idim);

  /* Get the indices of the grid node */

  if (point_to_grid(dbout, coor, flag_outside, iwork1.data()) < 0) return (-1);

  /* Convert the indices into the absolute grid node */

  int jech = dbout->indiceToRank(iwork1);

  return (jech);
}

/*****************************************************************************/
/*!
 **  Check if a sample of a Db lies within a Grid Db
 **
 ** \return  1 if the point lies within the grid; 0 otherwise
 **
 ** \param[in]  db      Db structure
 ** \param[in]  iech    Rank of the target sample
 ** \param[in]  dbgrid  Grid Db structure
 **
 *****************************************************************************/
int point_inside_grid(Db *db, int iech, const DbGrid *dbgrid)
{
  int ndim = db->getNDim();
  VectorDouble work1(ndim);
  VectorDouble work2(ndim);

  /* Process the grid rotation (if any) */

  for (int idim = 0; idim < ndim; idim++)
    work1[idim] = db->getCoordinate(iech, idim) - dbgrid->getX0(idim);

  dbgrid->getGrid().getRotation().rotateInverse(work1, work2);

  /* Calculate the grid indices */

  for (int idim = 0; idim < ndim; idim++)
  {
    int ix = (int) (floor(work2[idim] / dbgrid->getDX(idim) + 0.5));
    if (ix < 0 || ix >= dbgrid->getNX(idim)) return (0);
  }
  return (1);
}

/****************************************************************************/
/*!
 **  Monovariate statistics
 **
 ** \param[in]  db     Db structure
 ** \param[in]  iatt   Rank of the attribute
 **
 ** \param[out]  wtot Sum of the weights
 ** \param[out]  mean Mean of the variable
 ** \param[out]  var  variance of the variable
 ** \param[out]  mini Minimum value
 ** \param[out]  maxi Maximum value
 **
 *****************************************************************************/
void db_monostat(Db *db,
                 int iatt,
                 double *wtot,
                 double *mean,
                 double *var,
                 double *mini,
                 double *maxi)
{
  int iech;
  double weight, value;

  /* Initializations */

  (*mini) = 1.e30;
  (*maxi) = -1.e30;
  (*wtot) = (*mean) = (*var) = 0.;

  /* Loop on the data */

  for (iech = 0; iech < db->getSampleNumber(); iech++)
  {
    if (!db->isActive(iech)) continue;
    value = db->getArray(iech, iatt);
    if (FFFF(value)) continue;
    weight = db->getWeight(iech);
    if (weight <= 0.) continue;
    (*wtot) += weight;
    (*mean) += weight * value;
    (*var) += weight * value * value;
    if (value < (*mini)) (*mini) = value;
    if (value > (*maxi)) (*maxi) = value;
  }

  /* Normation */

  if (*wtot <= 0.)
  {
    *mean = *var = *wtot = TEST;
  }
  else
  {
    (*mean) /= (*wtot);
    (*var) = (*var) / (*wtot) - (*mean) * (*mean);
  }
}


/*****************************************************************************/
/*!
 **  Calculates the proportions of facies within a grid
 **
 ** \return  Error return code
 **
 ** \param[in]  db        Input Db structure
 ** \param[in]  dbgrid    Output Grid Db structure
 ** \param[in]  nfac1max  Maximum number of facies for the first variable
 ** \param[in]  nfac2max  Maximum number of facies for the second variable
 **
 ** \param[out] nclout    Total number of classes
 **
 ** \remark  This procedure is designed for one or two Z-variables
 ** \remark  When used for 2 variables, the proportions are for pairs of
 ** \remark  variables
 **
 *****************************************************************************/
int db_proportion(Db *db, DbGrid *dbgrid, int nfac1max, int nfac2max, int *nclout)
{
  int nval, mini, invalid;

  /* Initializations */

  int nvar = db->getLocNumber(ELoc::Z);
  int nech = db->getSampleNumber();
  if (nvar <= 0 || nvar > 2)
  {
    messerr("This procedure is designed for 1 or 2 variables");
    return 1;
  }
  if (!dbgrid->isGrid())
  {
    messerr(" This procedure is designed for a Grid Output Db");
    return 1;
  }

  VectorInt nmax(2);
  nmax[0] = nfac1max;
  nmax[1] = nfac2max;

  /* Count the number of facies */

  for (int ivar = 0; ivar < nvar; ivar++)
  {
    if (nmax[ivar] <= 0)
    {
      VectorDouble tab = db->getColumnByLocator(ELoc::Z, ivar);
      VectorDouble sel;
      if (db->hasLocVariable(ELoc::SEL))
        sel = db->getColumnByLocator(ELoc::SEL);
      ut_facies_statistics(nech, tab.data(), sel.data(), &nval, &mini, &nmax[ivar]);
    }
  }
  int nclass = VH::product(nmax);

  /* Core allocation */

  VectorDouble coor(db->getLocatorNumber(ELoc::X));

  /* Allocate the variables */

  int iptr = dbgrid->addColumnsByConstant(nclass, 0.);
  if (iptr < 0) return 1;
  dbgrid->setLocatorsByUID(nclass, iptr, ELoc::P);

  /* Loop on the samples of the input data Db */

  VectorInt ifac(2, 0);
  for (int iech = 0; iech < nech; iech++)
  {
    /* Load the facies information */

    for (int ivar = invalid = 0; ivar < nvar && invalid == 0; ivar++)
    {
      ifac[ivar] = (int) db->getZVariable(iech, ivar);
      if (ifac[ivar] > nmax[ivar]) invalid = 1;
    }
    if (invalid) continue;

    /* Locate the data sample within the grid */

    int jech = index_point_to_grid(db, iech, -1, dbgrid, coor.data());
    if (jech < 0) continue;

    /* Compute the class index */

    int iclass;
    if (nvar == 1)
      iclass = ifac[0] - 1;
    else
      iclass = (ifac[1] - 1) * nmax[0] + (ifac[0] - 1);

    /* Update the number of samples in the cell */

    dbgrid->setLocVariable(ELoc::P,jech, iclass,
                          dbgrid->getLocVariable(ELoc::P,jech, iclass) + 1);
  }

  /* Normalization phase */

  for (int jech = 0; jech < dbgrid->getSampleNumber(); jech++)
  {
    /* Cumulate the proportions */

    double total = 0.;
    for (int iclass = 0; iclass < nclass; iclass++)
      total += dbgrid->getLocVariable(ELoc::P,jech, iclass);
    if (total == 1.) continue;
    if (total <= 0.)
    {
      /* No sample in the current cell */

      for (int iclass = 0; iclass < nclass; iclass++)
        dbgrid->setLocVariable(ELoc::P,jech, iclass, TEST);
    }
    else
    {
      for (int iclass = 0; iclass < nclass; iclass++)
        dbgrid->setLocVariable(ELoc::P,jech, iclass,
                              dbgrid->getLocVariable(ELoc::P,jech, iclass) / total);
    }
  }

  *nclout = nclass;
  return 0;
}

/****************************************************************************/
/*!
 **  Merge a set of variables in a Db
 **
 ** \return  Error return code
 **
 ** \param[in]  db      Db structure
 ** \param[in]  ncol    Number of variables to be merged
 ** \param[in]  cols    Array of input variable columns
 **
 *****************************************************************************/
int db_merge(Db *db, int ncol, int *cols)
{
  int iptr, iech, icol;
  double value = TEST;

  /* Preliminary check */

  if (ncol < 1)
  {
    messerr("This procedure requires at least one variable to be merged");
    return (1);
  }

  /* Add the new variable */

  iptr = db->addColumnsByConstant(1, TEST);

  /* Loop on the samples */

  for (iech = 0; iech < db->getSampleNumber(); iech++)
  {

    /* Loop on the other variables until a defined value is found */

    for (icol = 0; icol < ncol; icol++)
    {
      value = db->getArray(iech, cols[icol]);
      if (!FFFF(value)) break;
    }

    /* Save the current value */

    db->setArray(iech, iptr, value);
  }

  return (0);
}

/****************************************************************************/
/*!
 **  Count the number of defined values for a Column
 **
 ** \return  Number of defined values
 **
 ** \param[in]  db      Db structure
 ** \param[in]  icol    Rank of the Column
 **
 *****************************************************************************/
int db_count_defined(Db *db, int icol)
{
  int iech, number;
  double value;

  number = 0;
  if (db == nullptr) return (number);
  for (iech = 0; iech < db->getSampleNumber(); iech++)
  {
    value = db->getArray(iech, icol);
    if (!FFFF(value)) number++;
  }
  return (number);
}

/****************************************************************************/
/*!
 **  Define the array of locators
 **
 ** \param[in]  strings     Array of locators
 ** \param[in]  current     Array of ranks of the last modified locator
 ** \param[in]  flag_locnew Reset all locators
 **
 ** \remark  The elements of the array current are numbered starting from 1
 **
 *****************************************************************************/
void db_locators_correct(VectorString &strings,
                         const VectorInt &current,
                         int flag_locnew)
{
  int cur_item, ref_item, found, nmatch, ncount, nmult;
  ELoc cur_type, ref_type;

  /* Dispatch */

  int number = static_cast<int>(strings.size());
  int ncur = static_cast<int>(current.size());
  if (number <= 0 || ncur <= 0) return;
  VectorInt rank(number);
  VectorInt ind(number);

  if (flag_locnew)
  {

    /* Undefine all the variables with the newly fixed locators */

    for (int i = 0; i < number; i++)
    {

      /* Check if the variable is one of the newly defined locator */

      found = -1;
      for (int j = 0; j < ncur && found < 0; j++)
        if (i == current[j] - 1) found = j;
      if (found < 0) continue;

      /* Get the locator characteristics of the current variable */

      if (locatorIdentify(strings[i], &cur_type, &cur_item, &nmult)) continue;

      /* Undefine the locator if: */
      /* - it coincides with a newly defined locator */
      /* - it is not part of the newly defined locators */

      for (int ip = 0; ip < number; ip++)
      {
        found = -1;
        for (int j = 0; j < ncur && found < 0; j++)
          if (ip == current[j] - 1) found = j;
        if (found >= 0) continue;

        if (locatorIdentify(strings[ip], &ref_type, &ref_item, &nmult))
          continue;
        if (cur_type == ref_type) strings[ip] = "NA";
      }
    }
  }
  else
  {

    /* Check locators similar to the ones just defined */

    for (int i = 0; i < number; i++)
    {

      /* Bypass the locator if newly defined */

      found = -1;
      for (int j = 0; j < ncur && found < 0; j++)
        if (i == current[j] - 1) found = j;
      if (found >= 0) continue;

      /* Get the characteristics of the current locator */

      if (locatorIdentify(strings[i], &cur_type, &cur_item, &nmult)) continue;

      /* Undefine the locator if it coincides with a newly defined locator */

      for (int j = 0; j < ncur; j++)
      {
        if (locatorIdentify(strings[current[j] - 1], &ref_type, &ref_item,
                            &nmult)) continue;
        if (cur_type == ref_type && cur_item == ref_item) strings[i] = "NA";
      }
    }
  }

  /* Loop on the reference locator */

  auto it = ELoc::getIterator();
  while (it.hasNext())
  {
    if (*it != ELoc::UNKNOWN)
    {
      /* Store the ranks of the locators matching the reference locator */
      nmatch = 0;
      for (int i = 0; i < number; i++)
      {
        if (locatorIdentify(strings[i], &cur_type, &cur_item, &nmult)) continue;
        if (cur_type != *it) continue;
        rank[nmatch++] = cur_item;
      }
      // Do not forget to increment the iterator!
      // 'continue' keyword should be forbidden!!
      if (nmatch <= 0)
      {
        it.toNext();
        continue;
      }

      /* Sort the indices */

      ncount = nmatch;
      for (int i = 0; i < ncount; i++) ind[i] = i;
      VH::arrangeInPlace(0, ind, rank, true, ncount);

      /* Store the ranks of the locators matching the reference locator */

      nmatch = 0;
      for (int i = 0; i < number; i++)
      {
        if (locatorIdentify(strings[i], &cur_type, &cur_item, &nmult)) continue;
        if (cur_type != *it) continue;
        found = -1;
        for (int k = 0; k < ncount && found < 0; k++)
          if (ind[k] == nmatch) found = k;
        strings[i] = getLocatorName(cur_type, found + 1);
        nmatch++;
      }
    }
    it.toNext();
  }
}

/****************************************************************************/
/*!
 **  Read the proportions per VPC
 **
 ** \return  Error return code
 **
 ** \param[in]  db       Db structure
 ** \param[in]  ix       Rank of the cell along X
 ** \param[in]  iy       Rank of the cell along Y
 **
 ** \param[out] props    Array of proportions
 **
 ** \remark  This procedure is meant for a 3-D grid file
 **
 *****************************************************************************/
int db_prop_read(DbGrid *db, int ix, int iy, double *props)
{
  int ecr, flag_no;
  double value, total;

  /* Initializations */

  int nprop = db->getLocNumber(ELoc::P);
  int ndim  = db->getNDim();
  int nz    = db->getNX(2);
  for (int i = 0; i < nz * nprop; i++) props[i] = 0.;

  /* Preliminary checks */

  if (ndim != 3) return (1);
  if (ix < 0 || ix >= db->getNX(0)) return (1);
  if (iy < 0 || iy >= db->getNX(1)) return (1);

  /* Blank out the array */

  VectorInt indices(ndim, 0);
  indices[0] = ix;
  indices[1] = iy;

  /* Load the proportions */

  for (int iz = ecr = 0; iz < nz; iz++)
  {
    indices[2] = iz;
    int iech = db->indiceToRank(indices);

    /* Check if the proportions are ALL defined */

    total = 0.;
    for (int iprop = flag_no = 0; iprop < nprop && flag_no == 0; iprop++)
    {
      value = db->getLocVariable(ELoc::P,iech, iprop);
      if (FFFF(value))
        flag_no = 1;
      else
        total += value;
    }

    for (int iprop = 0; iprop < nprop; iprop++, ecr++)
      props[ecr] = (flag_no && total > 0)
                   ? TEST : db->getLocVariable(ELoc::P, iech, iprop) / total;
  }
  return (0);
}

/****************************************************************************/
/*!
 **  Write the proportions (per VPC)
 **
 ** \return  Error return code
 **
 ** \param[in]  db       Db structure
 ** \param[in]  ix       Rank of the cell along X
 ** \param[in]  iy       Rank of the cell along Y
 ** \param[in]  props    Array of proportions
 **
 ** \remark  This procedure is meant for a 3-D grid file
 **
 *****************************************************************************/
int db_prop_write(DbGrid *db, int ix, int iy, double *props)
{
  int nprop = db->getLocNumber(ELoc::P);
  int nz    = db->getNX(2);
  int ndim = db->getNDim();

  /* Preliminary checks */

  if (ndim != 3) return (1);
  if (ix < 0 || ix >= db->getNX(0)) return (1);
  if (iy < 0 || iy >= db->getNX(1)) return (1);

  /* Blank out the array */

  VectorInt indices(ndim,0);
  indices[0] = ix;
  indices[1] = iy;

  /* Load the proportions */

  int ecr;
  for (int iz = ecr = 0; iz < nz; iz++)
  {
    indices[2] = iz;
    int iech = db->indiceToRank(indices);
    for (int iprop = 0; iprop < nprop; iprop++, ecr++)
      db->setLocVariable(ELoc::P,iech, iprop, props[ecr]);
  }
  return (0);
}

/****************************************************************************/
/*!
 **  Evaluate the array of distances between samples of two data sets
 **  Takes the selection into account (if any)
 **
 ** \return Array of distances
 **
 ** \param[in]  db1     Db first descriptor
 ** \param[in]  db2     Db second descriptor
 ** \param[in]  niso    Number of variables tested for isotopy
 ** \param[in]  mode    Type of array returned
 **                     0 : extreme distances
 **                     1 : the distance to the closest sample
 **                     2 : all point-to-point distances
 ** \param[in]  flag_same 1 if both Db coincide
 **
 ** \param[out] n1      First dimension of the returned array
 ** \param[out] n2      Second dimension of the returned array
 ** \param[out] dmin    Minimum distance
 ** \param[out] dmax    Maximum distance
 **
 ** \remarks The returned array must be freed by calling routine
 ** \remarks When the two Dbs coincide, the distance calculation excludes
 ** \remarks the comparison between one sample and itself
 **
 *****************************************************************************/
double* db_distances_general(Db *db1,
                             Db *db2,
                             int niso,
                             int mode,
                             int flag_same,
                             int *n1,
                             int *n2,
                             double *dmin,
                             double *dmax)
{
  int nech1, nech2, iech1, iech2, ecr, max_all, nvalid;
  double *dist, dlocmin, dloc, dist_min, dist_max;

  /* Preliminary calculations */

  *n1 = 0;
  *n2 = 0;
  nech1 = db1->getSampleNumber(true);
  nech2 = db2->getSampleNumber(true);
  dist = nullptr;
  max_all = nech1 * nech2;

  /* Preliminary checks */

  if (niso > db1->getLocNumber(ELoc::Z) || niso > db2->getLocNumber(ELoc::Z))
  {
    messerr("You ask for distances between samples with %d variables defined",
            niso);
    messerr("But the input 'Db' have %d and %d variables defined",
            db1->getLocNumber(ELoc::Z), db2->getLocNumber(ELoc::Z));
    return (dist);
  }

  /* Core allocation */

  if (mode > 0)
  {
    dist = (double*) mem_alloc(sizeof(double) * max_all, 0);
    if (dist == nullptr) return (dist);
    for (int i = 0; i < max_all; i++)
      dist[i] = 0.;
  }

  /* Loop on the second point */

  dist_min = 1.e30;
  dist_max = -1.e30;
  ecr = nvalid = 0;
  for (iech2 = 0; iech2 < nech2; iech2++)
  {
    if (! db2->isActive(iech2)) continue;
    if (! db2->isIsotopic(iech2, niso)) continue;
    nvalid++;
    dlocmin = 1.e30;

    /* Loop on the first point */

    for (iech1 = 0; iech1 < nech1; iech1++)
    {
      if (mode != 2 && flag_same && iech1 == iech2) continue;
      if (! db1->isActive(iech1)) continue;
      if (! db1->isIsotopic(iech1, niso)) continue;

      /* Calculate distance */

      dloc = distance_inter(db1, db2, iech1, iech2, NULL);
      if (dloc < dist_min) dist_min = dloc;
      if (dloc > dist_max) dist_max = dloc;

      if (mode == 1)
      {
        if (dloc < dlocmin) dlocmin = dloc;
      }
      else if (mode == 2)
      {
        dist[ecr++] = dloc;
      }
    }

    if (mode == 1) dist[ecr++] = dlocmin;
  }

  /* Reallocate the distance array */

  if (mode > 0 && ecr < max_all)
  {
    dist = (double*) mem_realloc((char* ) dist, sizeof(double) * ecr, 0);
    if (dist == nullptr) return (dist);
  }

  /* Returned arguments */

  *dmin = dist_min;
  *dmax = dist_max;
  if (mode == 1)
  {
    *n1 = ecr;
    *n2 = 1;
  }
  else if (mode == 2)
  {
    *n1 = nvalid;
    *n2 = nvalid;
  }
  return (dist);
}

/****************************************************************************/
/*!
 **  Check if a grid is a multiple of the other grid
 **
 ** \return  1 if the two grid are multiple; 0 otherwise
 **
 ** \param[in]  db1   Db1 small grid structure
 ** \param[in]  db2   Db1 coarse grid structure
 **
 *****************************************************************************/
int is_grid_multiple(DbGrid *db1, DbGrid *db2)
{
  int idim, ndim, error;
  double ratio, delta;

  /* Initializations */

  error = 1;

  /* Preliminary checks */

  if (!db1->hasSameDimension(db2)) return 1;

  /* Core allocation */

  ndim = db1->getNDim();
  VectorInt indg(ndim);
  VectorDouble perc(ndim);
  VectorDouble coor1(ndim);
  VectorDouble coor2(ndim);

  /* Check that the grid meshes are multiple */

  for (idim = 0; idim < ndim; idim++)
  {
    ratio = db2->getDX(idim) / db1->getDX(idim);
    if (!isInteger(ratio)) goto label_end;
  }

  /* Get the lower left corners of the both grid */

  for (idim = 0; idim < ndim; idim++)
  {
    indg[idim] = 0;
    perc[idim] = -0.5;
  }
  db1->indicesToCoordinateInPlace(indg, coor1, perc);
  db2->indicesToCoordinateInPlace(indg, coor2, perc);

  /* Check that these corners are close enough */

  for (idim = 0; idim < ndim; idim++)
  {
    delta = (coor1[idim] - coor2[idim]) / db1->getDX(idim);
    if (ABS(delta) > EPSILON3) goto label_end;
  }

  /* Set the error return code */

  error = 0;

  label_end:
  return (1 - error);
}

/****************************************************************************/
/*!
 **  Transform a set of gradients defined by (modulus,angle) into (gx,gy)
 **  Only defined in the 2-D case
 **
 ** \return Error return code
 **
 ** \param[in]  db        Initial Db
 ** \param[in]  ang_conv  Convention on the angles
 ** \li                   1: Trigonometry (From East counter-clockwise)
 ** \li                   2: From North clockwise
 ** \param[in]  iad_mod   Rank of the 'modulus' attribute
 ** \param[in]  iad_ang   Rank of the 'angle' attribute (defined in degrees)
 ** \param[in]  iad_gx    Rank of the 'gx' attribute (defined in degrees)
 ** \param[in]  iad_gy    Rank of the 'gy' attribute (defined in degrees)
 **
 ** \remarks  Attributes 'gx' and 'gy' may coincide with 'modulus' and 'angle'
 **
 *****************************************************************************/
int db_gradient_modang_to_component(Db *db,
                                    int ang_conv,
                                    int iad_mod,
                                    int iad_ang,
                                    int iad_gx,
                                    int iad_gy)
{
  double angdeg, angrad, modulus;

  /* Load the information */

  VectorDouble v1 = db->getColumnByColIdx(iad_mod);
  VectorDouble v2 = db->getColumnByColIdx(iad_ang);

  /* Gradient conversion */

  for (int iech = 0; iech < db->getSampleNumber(); iech++)
  {
    if (FFFF(v1[iech]) || FFFF(v2[iech])) continue;
    modulus = v1[iech];
    angdeg = v2[iech];
    if (ang_conv == 2) angdeg = 90. - angdeg;
    angrad = ut_deg2rad(angdeg);
    v1[iech] = modulus * cos(angrad);
    v2[iech] = modulus * sin(angrad);
  }

  /* Save the information */

  db->setColumnByColIdx(v1, iad_gx);
  db->setColumnByColIdx(v2, iad_gy);

  return 0;
}

/****************************************************************************/
/*!
 **  Transform a set of gradients defined by (gx,gy) into (module,angle)
 **  Only defined in the 2-D case
 **
 ** \return Error return code
 **
 ** \param[in]  db        Initial Db
 ** \param[in]  verbose   1 for the verbose option
 ** \param[in]  iad_gx    Rank of the 'gx' attribute (defined in degrees)
 ** \param[in]  iad_gy    Rank of the 'gy' attribute (defined in degrees)
 ** \param[in]  iad_mod   Rank of the 'modulus' attribute
 ** \param[in]  iad_ang   Rank of the 'angle' attribute (defined in degrees)
 ** \param[in]  scale     Scaling factor for the modulus
 ** \param[in]  ve        Moderation factor
 **
 ** \remarks  Attributes 'gx' and 'gy' may coincide with 'modulus' and 'angle'
 ** \remarks  Nothing is done (more than the conversion if scale=1 and ve=0)
 **
 *****************************************************************************/
int db_gradient_component_to_modang(Db *db,
                                    int verbose,
                                    int iad_gx,
                                    int iad_gy,
                                    int iad_mod,
                                    int iad_ang,
                                    double scale,
                                    double ve)
{
  double norme, angle, vmax, surr, alpha, mini, maxi;

  /* Load the information */

  VectorDouble v1 = db->getColumnByColIdx(iad_gx);
  VectorDouble v2 = db->getColumnByColIdx(iad_gy);

  /* Convert gradient components into modulus and azimuth */

  vmax = 0.;
  for (int iech = 0; iech < db->getSampleNumber(); iech++)
  {
    if (FFFF(v1[iech]) || FFFF(v2[iech])) continue;
    norme = sqrt(v1[iech] * v1[iech] + v2[iech] * v2[iech]);
    angle = atan2(v2[iech], v1[iech]);
    v1[iech] = norme;
    v2[iech] = ut_rad2deg(angle);
    if (norme > vmax) vmax = norme;
  }

  /* Modify the local anisotropy */

  mini = 1.e30;
  maxi = -1.e30;
  for (int iech = 0; iech < db->getSampleNumber(); iech++)
  {
    if (!db->isActive(iech)) continue;
    alpha = 1. / (1. + ve);
    surr = 1 + v1[iech] * (1. - alpha) / (vmax * alpha);
    v1[iech] = surr * scale;
    if (surr < mini) mini = surr;
    if (surr > maxi) maxi = surr;
  }

  /* Save the information */

  db->setColumnByColIdx(v1, iad_mod);
  db->setColumnByColIdx(v2, iad_ang);

  /* Print statistics (optional) */

  if (verbose)
  {
    mestitle(1, "Range correction");
    message("Value of the vector effect = %lf\n", ve);
    message("Range correction varies between %lf and %lf\n", mini, maxi);
  }

  return 0;
}

/****************************************************************************/
/*!
 **  Returns the relative rank of a sample from its absolute rank
 **  These are different due to the presence of a selection
 **  Returns -1 if not found
 **
 ** \return  Relative rank of a sample
 **
 ** \param[in]  db    Db structure
 ** \param[in]  iech0 Absolute sample rank
 **
 *****************************************************************************/
int db_get_rank_absolute_to_relative(Db *db, int iech0)
{
  int iech, jech;

  if (!db->hasLocVariable(ELoc::SEL)) return (iech0);

  for (iech = jech = 0; iech < db->getSampleNumber(); iech++)
  {
    if (!db->isActive(iech)) continue;
    if (iech == iech0) return (jech);
    jech++;
  }
  return (-1);
}

/****************************************************************************/
/*!
 **  Returns the absolute rank of a sample from its relative rank
 **  These are different due to the presence of a selection
 **  Returns -1 if not found
 **
 ** \return  Relative rank of a sample
 **
 ** \param[in]  db    Db structure
 ** \param[in]  iech0 Relative sample rank
 **
 *****************************************************************************/
int db_get_rank_relative_to_absolute(Db *db, int iech0)
{
  int iech, jech;

  if (!db->hasLocVariable(ELoc::SEL)) return (iech0);

  for (iech = jech = 0; iech < db->getSampleNumber(); iech++)
  {
    if (!db->isActive(iech)) continue;
    if (jech == iech0) return (iech);
    jech++;
  }
  return (-1);
}

/****************************************************************************/
/*!
 **  Returns a value from the a 3-D (maximum) grid
 **
 ** \return  Returned value
 **
 ** \param[in]  dbgrid   Db Grid structure
 ** \param[in]  iptr     Rank of the column
 ** \param[in]  indg     Working index array (Dimension: get_NDIM(dbgrid))
 ** \param[in]  ix       Rank of the node along first dimension
 ** \param[in]  iy       Rank of the node along second dimension
 ** \param[in]  iz       Rank of the node along third dimension
 **
 *****************************************************************************/
double get_grid_value(DbGrid *dbgrid, int iptr, VectorInt& indg, int ix, int iy, int iz)
{
  int ndim, iad;
  double value;

  ndim = dbgrid->getNDim();
  if (ndim >= 1) indg[0] = ix;
  if (ndim >= 2) indg[1] = iy;
  if (ndim >= 3) indg[2] = iz;

  iad = dbgrid->indiceToRank(indg);
  value = dbgrid->getArray(iad, iptr);
  return (value);
}

/****************************************************************************/
/*!
 **  Set the value in the a 3-D (maximum) grid
 **
 ** \param[in]  dbgrid   Db Grid structure
 ** \param[in]  iptr     Rank of the column
 ** \param[in]  indg     Working index array (Dimension: get_NDIM(dbgrid))
 ** \param[in]  ix       Rank of the node along first dimension
 ** \param[in]  iy       Rank of the node along second dimension
 ** \param[in]  iz       Rank of the node along third dimension
 ** \param[in]  value    Assigned value
 **
 *****************************************************************************/
void set_grid_value(DbGrid *dbgrid,
                    int iptr,
                    VectorInt& indg,
                    int ix,
                    int iy,
                    int iz,
                    double value)
{
  int ndim, iad;

  ndim = dbgrid->getNDim();
  if (ndim >= 1) indg[0] = ix;
  if (ndim >= 2) indg[1] = iy;
  if (ndim >= 3) indg[2] = iz;

  iad = dbgrid->indiceToRank(indg);
  dbgrid->setArray(iad, iptr, value);
}

/****************************************************************************/
/*!
 **  Extract the subgrid (from a grid) which contains the only cells where
 **  the target variable lies wuthin the target interval
 **  Selection in the input grid is taken into account
 **
 ** \param[in]  db_grid       Db structure
 ** \param[in]  iptr          Rank of the column of the target variable
 ** \param[in]  margin        Array of margins (or NULL)
 ** \param[in]  limmin        Array of minimum dimensions (or NULL)
 ** \param[in]  flag_sel      Create the selection
 ** \param[in]  flag_copy     1 if the selection must be copied in sub-grid
 ** \param[in]  verbose       Verbose flag
 ** \param[in]  vmin          Lower bound (included)
 ** \param[in]  vmax          Upper bound (excluded)
 **
 *****************************************************************************/
DbGrid* db_grid_reduce(DbGrid *db_grid,
                       int iptr,
                       const int *margin,
                       const int *limmin,
                       int flag_sel,
                       int flag_copy,
                       int verbose,
                       double vmin,
                       double vmax)
{
  DbGrid *ss_grid;
  int error, ndim, nech, flag_refuse, isel, icopy, iech;
  int mini, maxi, ecart, size;
  double value, retval;
  VectorInt nx;
  VectorDouble x0;

  // Initializations */

  error = 1;
  ss_grid = nullptr;

  // Core allocation

  nech = db_grid->getSampleNumber();
  ndim = db_grid->getNDim();
  VectorInt indcur(ndim);
  VectorInt indmin(ndim, nech);
  VectorInt indmax(ndim, -1);
  VectorDouble coor(ndim, 0.);

  // Loop on the input grid

  for (int i = 0; i < db_grid->getSampleNumber(); i++)
  {
    if (!db_grid->isActive(i)) continue;
    value = db_grid->getArray(i, iptr);
    if (value < vmin || value >= vmax) continue;
    db_grid->rankToIndice(i, indcur);

    for (int idim = 0; idim < ndim; idim++)
    {
      if (indcur[idim] < indmin[idim]) indmin[idim] = indcur[idim];
      if (indcur[idim] > indmax[idim]) indmax[idim] = indcur[idim];
    }
  }

  // Calculate extension of the next grid

  flag_refuse = 0;
  for (int idim = 0; idim < ndim && !flag_refuse; idim++)
  {
    if (indmin[idim] > indmax[idim]) flag_refuse = 1;
    size = db_grid->getNX(idim);
    mini = indmin[idim];
    if (margin != nullptr) mini = MAX(0, mini - margin[idim]);
    maxi = indmax[idim];
    if (margin != nullptr) maxi = MIN(size - 1, maxi + margin[idim]);
    if (limmin != nullptr)
    {
      ecart = limmin[idim];
      if (ecart > 0)
      {
        mini -= ecart / 2;
        maxi += ecart / 2;
      }
    }
    indmin[idim] = mini;
    indmax[idim] = maxi - mini + 1;
  }

  // Accept the new grid or refuse it

  if (flag_refuse)
  {
    error = 0;
    goto label_end;
  }

  // Optional printout

  if (verbose)
  {
    mestitle(1, "Grid Extraction");
    message("From:");
    for (int idim = 0; idim < ndim; idim++)
      message(" %d", indmin[idim]);
    message("\n");
    message("Size:");
    for (int idim = 0; idim < ndim; idim++)
      message(" %d", indmax[idim]);
    message("\n");
  }

  // Create the new sub-grid

  db_grid->indicesToCoordinateInPlace(indmin, coor);
  nx = indmax;
  x0 = coor;
  ss_grid = DbGrid::create(nx, db_grid->getDXs(), x0, db_grid->getAngles(),
                           ELoadBy::COLUMN, VectorDouble(),
                           VectorString(), VectorString(), 0);

  // Create the selection (optional)

  if (flag_sel)
  {
    isel = ss_grid->addColumnsByConstant(1, 0., String(), ELoc::SEL);
    for (int i = 0; i < ss_grid->getSampleNumber(); i++)
    {
      ss_grid->rankToIndice(i, indcur);
      for (int idim = 0; idim < ndim; idim++) indcur[idim] += indmin[idim];
      iech = db_grid->indiceToRank(indcur);
      value = db_grid->getArray(iech, iptr);
      retval = (value >= vmin && value < vmax);
      ss_grid->setArray(i, isel, retval);
    }
  }

  // Copy the selection from initial grid to subgrid

  if (flag_copy)
  {
    icopy = ss_grid->addColumnsByConstant(1, 0., String(), ELoc::SEL);
    for (int i = 0; i < ss_grid->getSampleNumber(); i++)
    {
      ss_grid->rankToIndice(i, indcur);
      for (int idim = 0; idim < ndim; idim++) indcur[idim] += indmin[idim];
      iech = db_grid->indiceToRank(indcur);

      retval = 1.;
      if (!db_grid->isActive(iech)) retval = 0.;
      if (retval)
      {
        value = db_grid->getArray(iech, iptr);
        if (value < vmin || value >= vmax) retval = 0.;
      }
      ss_grid->setArray(i, icopy, retval);
    }
  }

  // Set the error return code

  error = 0;

  label_end:
  if (error)
  {
    delete ss_grid;
    ss_grid = nullptr;
  }
  return (ss_grid);
}

/****************************************************************************/
/*!
 **  Patch a sub-grid within a main grid
 **
 ** \return Error return code
 **
 ** \param[in]  ss_grid       Db sub-grid structure
 ** \param[in]  db_grid       Db main grid structure
 ** \param[in]  iptr_ss       Rank of the attribute in the sub-grid
 ** \param[in]  iptr_db       Rank of the attribute in the main grid
 ** \param[in]  iptr_rank     Rank of the attribute storing object rank
 **                           If <0, no check is performed: always patch
 ** \param[in]  new_rank      Rank of the current object to patch in main grid
 ** \param[in]  oper          >0 for larger; <0 for smaller
 ** \param[in]  verbose       Verbose flag
 **
 ** \remarks When 'iptr_rank' is defined (>=0), defined pixels of the current
 ** \remarks sub-grid overwrites the corresponding pixel within the main grid
 ** \remarks only if its rank ('new_rank') is larger (oper>0) or smaller
 ** \remarks (oper<0)) than the rank of the same pixel in the main grid
 ** \remarks (attribute 'iptr_rank')
 ** \remarks When 'iptr_rank' is undefined, arguments 'new_rank', 'oper' are
 ** \remarks useless
 **
 *****************************************************************************/
int db_grid_patch(DbGrid *ss_grid,
                  DbGrid *db_grid,
                  int iptr_ss,
                  int iptr_db,
                  int iptr_rank,
                  int new_rank,
                  int oper,
                  int verbose)
{
  int error, ndim, jech, nused, noused, nout, nundef, nmask, ndef, nbnomask, flag_save;
  double *coor1, *coor2, value, rank;
  VectorInt indg;
  VectorInt indg0;

  /* Initializations */

  error = 1;
  ndim = ss_grid->getNDim();
  coor1 = coor2 = nullptr;

  /* Check that the two grids are compatible */

  if (!db_grid->hasSameDimension(ss_grid)) goto label_end;
  if (!db_grid->isSameGridMesh(*ss_grid)) goto label_end;
  if (!db_grid->isSameGridRotation(*ss_grid)) goto label_end;

  /* Core allocation */

  coor1 = (double*) mem_alloc(sizeof(double) * ndim, 0);
  if (coor1 == nullptr) goto label_end;
  coor2 = (double*) mem_alloc(sizeof(double) * ndim, 0);
  if (coor2 == nullptr) goto label_end;
  indg0.resize(ndim);
  indg.resize(ndim);

  /* Find the coordinates of the origin of the sub-grid within the main grid */

  for (int idim = 0; idim < ndim; idim++)
    coor1[idim] = ss_grid->getX0(idim);
  (void) point_to_grid(db_grid, coor1, -1, indg0.data());
  if (point_to_grid(db_grid, coor1, -1, indg0.data()) == -1)
  {
    messerr("Subgrid origin does not lie within the main grid");
    db_grid_print(db_grid);
    messerr("Subgrid origin:");
    for (int idim = 0; idim < ndim; idim++)
      messerr("- Dimension #%d: Coordinate=%lf", idim + 1, coor1[idim]);
    goto label_end;
  }

  /* Patch the values */

  nused = noused = nundef = nout = nmask = 0;
  for (int iech = 0; iech < ss_grid->getSampleNumber(); iech++)
  {

    // Sample masked off in the subgrid
    if (!ss_grid->isActive(iech))
    {
      nmask++;
      continue;
    }

    // Get the value
    value = ss_grid->getArray(iech, iptr_ss);

    // The value of the subgrid is undefined
    if (FFFF(value))
    {
      nundef++;
      continue;
    }

    // Find the location of corresponding pixel in the main grid
    ss_grid->rankToIndice(iech, indg);
    for (int idim = 0; idim < ndim; idim++) indg[idim] += indg0[idim];
    jech = db_grid->indiceToRank(indg);

    if (jech < 0)
    {
      // The pixel is outside the main grid
      nout++;
    }
    else
    {
      // The pixel is inside the main grid
      flag_save = 0;

      // The rank attribute is defined
      if (iptr_rank >= 0)
      {
        rank = db_grid->getArray(jech, iptr_rank);
        if (FFFF(rank))
          flag_save = 1;
        else
        {
          if (oper > 0)
          {
            if (new_rank > rank) flag_save = 1;
          }
          else
          {
            if (new_rank < rank) flag_save = 1;
          }
        }
      }
      else
      {
        flag_save = 1;
      }

      // Is the new value saved
      if (flag_save)
      {
        db_grid->setArray(jech, iptr_db, value);
        nused++;
      }
      else
      {
        noused++;
      }
    }
  }

  /* Optional printout */

  if (verbose)
  {
    ndef = nbnomask = 0;
    for (int iech = 0; iech < db_grid->getSampleNumber(); iech++)
    {
      if (!db_grid->isActive(iech)) continue;
      value = db_grid->getArray(iech, iptr_db);
      nbnomask++;
      if (!FFFF(value)) ndef++;
    }

    message("Patching:\n");
    message("(Naming convention: *_S for subgrid and _G for main grid)\n");
    for (int idim = 0; idim < ndim; idim++)
      message("- Dimension %d: NX_S =%4d - NX_G =%4d - Shift =%4d\n", idim + 1,
              ss_grid->getNX(idim), db_grid->getNX(idim), indg0[idim]);
    message("Subgrid                               = %d\n",
            ss_grid->getSampleNumber());
    message("- Number of masked off samples        = %d\n", nmask);
    message("- Number of undefined values          = %d\n", nundef);
    message("- Number of samples outside main grid = %d\n", nout);
    message("- Number of valid values (save)       = %d\n", nused);
    message("- Number of valid values (skip)       = %d\n", noused);
    message("Main Grid                             = %d\n",
            db_grid->getSampleNumber());
    message("- Number of non-masked values         = %d\n", nbnomask);
    message("- Number of valid values              = %d\n", ndef);
  }

  /* Set the error returned code */

  error = 0;

  label_end:
  mem_free((char* ) coor1);
  mem_free((char* ) coor2);
  return (error);
}

/****************************************************************************/
/*!
 **  Identify the attribute by its name
 **
 ** \return  Rank of the variable starting from 0 (or -1 if not found)
 **
 ** \param[in]  db       Db descriptor
 ** \param[in]  string   attribute name
 **
 *****************************************************************************/
int db_name_identify(Db *db, const String &string)
{
  for (int iatt = 0, natt = db->getUIDMaxNumber(); iatt < natt; iatt++)
  {
    int icol = db->getColIdxByUID(iatt);
    if (string != db->getNameByColIdx(icol)) return iatt;
  }
  return (-1);
}

/****************************************************************************/
/*!
 **  Update the bounding box limits due to a point
 **
 ** \param[in]  ndim     Space dimension
 ** \param[in]  rotmat   Rotation matrix (optional)
 ** \param[in]  coor     Array of sample coordinates
 **
 ** \param[out]  mini    Array containing the minimum
 ** \param[out]  maxi    Array containing the maximum
 **
 ** \remarks This function is limited to the 3-D case
 **
 *****************************************************************************/
static void st_rotate(int ndim,
                      double *rotmat,
                      VectorDouble& coor,
                      VectorDouble& mini,
                      VectorDouble& maxi)
{
  double d2[3];

  matrix_product_safe(1, ndim, ndim, coor.data(), rotmat, d2);
  for (int idim = 0; idim < ndim; idim++)
  {
    mini[idim] = getMin(mini[idim], d2[idim]);
    maxi[idim] = getMax(maxi[idim], d2[idim]);
  }
}

/****************************************************************************/
/*!
 **  Returns the extension of the field (after rotation) along each axis
 **
 ** \param[in]   db        Db structure
 ** \param[in]   rotmat    Rotation matrix (optional)
 ** \param[out]  mini      Array containing the minimum (Dimension = ndim)
 ** \param[out]  maxi      Array containing the maximum (Dimension =  ndim)
 **
 ** \remarks This function does nothing if:
 ** \remarks - no rotation matrix is defined
 ** \remarks - space dimension is 1 or larger than 3
 **
 *****************************************************************************/
void db_extension_rotated(Db *db,
                          double *rotmat,
                          VectorDouble& mini,
                          VectorDouble& maxi)
{
  int ndim = db->getNDim();
  VectorDouble coor(ndim);
  VectorDouble minrot(ndim);
  VectorDouble maxrot(ndim);

  // Calculate the extension (without rotation)
  db->getExtensionInPlace(mini, maxi);

  // Bypass the calculations
  if (rotmat == nullptr) return;
  if (ndim <= 1 || ndim > 3) return;

  // Perform the rotation
  for (int idim = 0; idim < ndim; idim++)
    minrot[idim] = maxrot[idim] = TEST;

  if (ndim == 2)
  {

    // 2-D case

    coor[0] = mini[0];        // 0,0
    coor[1] = mini[1];
    st_rotate(ndim, rotmat, coor, minrot, maxrot);
    coor[0] = mini[0];        // 0,1
    coor[1] = maxi[1];
    st_rotate(ndim, rotmat, coor, minrot, maxrot);
    coor[0] = maxi[0];        // 1,0
    coor[1] = mini[1];
    st_rotate(ndim, rotmat, coor, minrot, maxrot);
    coor[0] = maxi[0];        // 1,1
    coor[1] = maxi[1];
    st_rotate(ndim, rotmat, coor, minrot, maxrot);
  }
  else
  {

    // 3-D case

    coor[0] = mini[0];        // 0,0,0
    coor[1] = mini[1];
    coor[2] = mini[2];
    st_rotate(ndim, rotmat, coor, minrot, maxrot);
    coor[0] = mini[0];        // 0,0,1
    coor[1] = mini[1];
    coor[2] = maxi[2];
    st_rotate(ndim, rotmat, coor, minrot, maxrot);
    coor[0] = mini[0];        // 0,1,0
    coor[1] = maxi[1];
    coor[2] = mini[2];
    st_rotate(ndim, rotmat, coor, minrot, maxrot);
    coor[0] = mini[0];        // 0,1,1
    coor[1] = maxi[1];
    coor[2] = maxi[2];
    st_rotate(ndim, rotmat, coor, minrot, maxrot);
    coor[0] = maxi[0];        // 1,0,0
    coor[1] = mini[1];
    coor[2] = mini[2];
    st_rotate(ndim, rotmat, coor, minrot, maxrot);
    coor[0] = maxi[0];        // 1,0,1
    coor[1] = mini[1];
    coor[2] = maxi[2];
    st_rotate(ndim, rotmat, coor, minrot, maxrot);
    coor[0] = maxi[0];        // 1,1,0
    coor[1] = maxi[1];
    coor[2] = mini[2];
    st_rotate(ndim, rotmat, coor, minrot, maxrot);
    coor[0] = maxi[0];        // 1,1,1
    coor[1] = maxi[1];
    coor[2] = maxi[2];
    st_rotate(ndim, rotmat, coor, minrot, maxrot);
  }

  for (int idim = 0; idim < ndim; idim++)
  {
    mini[idim] = minrot[idim];
    maxi[idim] = maxrot[idim];
  }
}

/****************************************************************************/
/*!
 **  Returns a vector containing the coordinates of a Grid
 **  along one space dimension
 **
 ** \param[in]   dbgrid    Db structure
 ** \param[in]   idim      Rank of the space dimension
 **
 ****************************************************************************/
VectorDouble db_get_grid_axis(DbGrid *dbgrid, int idim)
{
  VectorDouble vect;

  if (dbgrid == nullptr) return (vect);
  if (! dbgrid->isGrid()) return (vect);
  if (idim < 0 || idim >= dbgrid->getNDim()) return (vect);

  int nvect = dbgrid->getNX(idim);
  double origin = dbgrid->getX0(idim);
  double pas = dbgrid->getDX(idim);
  vect.resize(nvect);

  for (int i = 0; i < nvect; i++)
    vect[i] = origin + i * pas;
  return vect;
}

/****************************************************************************/
/*!
 **  Returns a vector containing an attribute from the Db
 **
 ** \param[in]   db        Db structure
 ** \param[in]   iatt      Attribute rank
 ** \param[in]   verbose   Verbose flag
 **
 ****************************************************************************/
VectorDouble db_get_attribute(Db *db, int iatt, bool verbose)
{
  VectorDouble vect;

  if (db == nullptr)
  {
    if (verbose) messerr("Function 'db_get_attribute' requires a valid 'Db'");
    return (vect);
  }
  int nvect = db->getSampleNumber();
  vect.resize(nvect);

  for (int i = 0; i < nvect; i++)
    vect[i] = db->getArray(i, iatt);

  return vect;
}

/****************************************************************************/
/*!
 **  Identify the variables of a Db where names match a criterion
 **
 ** \param[in]   db       Db structure
 ** \param[in]   pattern  Matching pattern
 **
 ****************************************************************************/
VectorInt db_identify_variables_by_name(Db *db, const String &pattern)
{
  VectorString names = db->getName(pattern);
  VectorInt ranks = db->getUIDs(names);
  return ranks;
}

/****************************************************************************/
/*!
 **  Initialize the Grid iterator
 **
 ****************************************************************************/
void grid_iterator_init(Grid *grid, const VectorInt &order)
{
  grid->iteratorInit(order);
}

/****************************************************************************/
/*!
 **  Returns 1 when the last element of the iteration is reached
 **
 **  Increment the Grid iterator
 **
 ****************************************************************************/
VectorInt grid_iterator_next(Grid *grid)
{
  VectorInt indices = grid->iteratorNext();
  return (indices);
}

