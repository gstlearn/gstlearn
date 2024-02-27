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
#include "geoslib_f.h"
#include "geoslib_old_f.h"

#include "Enum/EJustify.hpp"
#include "Enum/ECalcVario.hpp"

#include "Space/SpaceRN.hpp"
#include "Variogram/Vario.hpp"
#include "Variogram/VarioParam.hpp"
#include "Anamorphosis/AAnam.hpp"
#include "Anamorphosis/AnamHermite.hpp"
#include "Polynomials/Hermite.hpp"
#include "Polygon/Polygons.hpp"
#include "Morpho/Morpho.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/Utilities.hpp"
#include "Basic/NamingConvention.hpp"
#include "Basic/File.hpp"
#include "Basic/String.hpp"
#include "Basic/OptDbg.hpp"
#include "Db/Db.hpp"
#include "Model/Model.hpp"
#include "Stats/PCA.hpp"
#include "Stats/PCAStringFormat.hpp"
#include "Matrix/MatrixSquareGeneral.hpp"
#include "Geometry/GeometryHelper.hpp"
#include "Geometry/BiTargetCheckGeometry.hpp"

#include <string.h>
#include <math.h>

/****************************************************************************/
/*!
 **  Replace zero values by TEST values
 **
 ** \param[in]  db    Discretization Grid descriptor
 ** \param[in]  iptr  Pointer of the attribute to be modified
 **
 *****************************************************************************/
static void st_final_discretization_grid(Db *db, int iptr)
{
  int iech, nech;
  double value;

  nech = db->getSampleNumber();
  for (iech = 0; iech < nech; iech++)
  {
    value = db->getArray(iech, iptr);
    if (value != 0.) continue;
    db->setArray(iech, iptr, TEST);
  }
}

/****************************************************************************/
/*!
 **  Add one pick to the discretization grid
 **
 ** \return  Index of the grid cell (or -1)
 **
 ** \param[in]  db    Discretization Grid descriptor
 ** \param[in]  x     Coordinate along the first axis
 ** \param[in]  y     Coordinate along the first axis
 **
 *****************************************************************************/
static int st_update_discretization_grid(DbGrid *db, double x, double y)
{
  int iech, ix, iy, indg[2];

  ix = (int) floor((x - db->getX0(0)) / db->getDX(0) + 0.5);
  iy = (int) floor((y - db->getX0(1)) / db->getDX(1) + 0.5);
  if (ix < 0 || ix >= db->getNX(0)) return (-1);
  if (iy < 0 || iy >= db->getNX(1)) return (-1);
  indg[0] = ix;
  indg[1] = iy;
  iech = db_index_grid_to_sample(db, indg);
  return (iech);
}

/****************************************************************************/
/*!
 **  Evaluate the correlation
 **     Correl(Z1(x) , Z2(x))
 **
 ** \return  Array of the indices of pairs of samples (or VectorVectorInt())
 **
 ** \param[in]  db1          Db descriptor (first variable)
 ** \param[in]  db2          Db descriptor (second variable for flag.same=T)
 ** \param[in]  name1        Name of the first variable
 ** \param[in]  name2        Name of the second variable
 ** \param[in]  flagFrom1    Start numbering of indices from 1 if True
 ** \param[in]  verbose      Verbose flag
 **
 ** \remarks The two input Db must match exactly (same number of samples with
 ** \remarks same set of coordinates and same optional selection)
 **
 ** \remarks The returned Vector of Vector of integer 'indices' contain
 ** \remarks the set of indices of the pairs of samples.
 ** \remarks Its contents is i1,j1,i2,j2,...
 ** \remarks The indices are numbered starting from 0
 **
 *****************************************************************************/
VectorVectorInt correlationPairs(Db *db1,
                                 Db *db2,
                                 const String& name1,
                                 const String& name2,
                                 bool flagFrom1,
                                 bool verbose)
{
  VectorVectorInt indices;

  /* Initializations */

  if (db1 == nullptr) return indices;
  if (db2 == nullptr) return indices;
  if (db1->getNDim() != db2->getNDim() || db1->getActiveSampleNumber() != db2->getActiveSampleNumber())
  {
    messerr("The two input 'db' are not compatible");
    return indices;
  }

  int nech = db1->getSampleNumber();
  int ndim = db1->getNDim();
  int shift = (flagFrom1) ? 1 : 0;
  SpaceRN space(ndim);
  SpaceTarget T1(&space);
  SpaceTarget T2(&space);

  /* Regular correlation */

  indices.resize(2);
  int nb = 0;
  for (int iech = 0; iech < nech; iech++)
  {
    if (!db1->isActive(iech)) continue;
    double val1 = db1->getValue(name1, iech);
    if (FFFF(val1)) continue;
    double val2 = db2->getValue(name2, iech);
    if (FFFF(val2)) continue;

    indices[0].push_back(iech + shift);
    indices[1].push_back(iech + shift);
    nb++;
  }

  /* Messages */

  if (nb <= 0)
  {
    messerr("No sample found where all variables are defined");
    return indices;
  }
  else
  {
    if (verbose)
    {
      message("Total number of samples = %d\n", nech);
      message("Number of samples defined = %d\n", (int) nb);
    }
  }
  return indices;
}

/****************************************************************************/
/*!
 **  Evaluate the shifted correlation calculated as follows:
 **     Correl(Z1(x) , Z2(x+h))
 **
 ** \return  Vector of indices (or VectorVectorInt())
 **
 ** \param[in]  db           Db descriptor
 ** \param[in]  name1        Name of the first variable
 ** \param[in]  name2        Name of the second variable
 ** \param[in]  varioparam   pointer to a VarioParam structure
 ** \param[in]  ipas         Rank of the lag of interest
 ** \param[in]  idir         Rank of the direction of interest (within VarioParam)
 ** \param[in]  verbose      Verbose flag
 **
 ** \remarks The returned Vector of Vector of integer 'indices' contain
 ** \remarks the set of indices of the pairs of samples.
 ** \remarks Its contents is i1,j1,i2,j2,...
 ** \remarks The indices are numbered starting from 1
 **
 *****************************************************************************/
VectorVectorInt hscatterPairs(Db *db,
                              const String& name1,
                              const String& name2,
                              VarioParam *varioparam,
                              int ipas,
                              int idir,
                              bool verbose)
{
  VectorVectorInt indices;
  double dist = 0.;

  // Preliminary checks

  if (db == nullptr) return indices;
  if (varioparam == nullptr) return indices;
  if (idir < 0 || idir >= varioparam->getDirectionNumber()) return indices;

  /* Initializations */

  const DirParam dirparam = varioparam->getDirParam(idir);
  int nech = db->getSampleNumber();
  int ndim = db->getNDim();
  SpaceRN space(ndim);
  SpaceTarget T1(&space);
  SpaceTarget T2(&space);
  indices.resize(2);

  // Creating a local Vario structure (to constitute the BiTargetCheck list)

  Vario *vario = Vario::create(*varioparam);
  vario->setDb(db);
  if (vario->prepare()) return 1;

  // Local variables to speed up calculations
  bool hasSel = db->hasLocVariable(ELoc::SEL);

  int nb = 0;
  for (int iech = 0; iech < nech - 1; iech++)
  {
    if (hasSel && !db->isActive(iech)) continue;
    double val1 = db->getValue(name1, iech);
    if (FFFF(val1)) continue;
    db->getSampleAsST(iech, T1);

    for (int jech = iech + 1; jech < nech; jech++)
    {
      if (hasSel && !db->isActive(jech)) continue;
      double val2 = db->getValue(name2, jech);
      if (FFFF(val2)) continue;
      db->getSampleAsST(jech, T2);

      // Reject the point as soon as one BiTargetChecker is not correct
      if (!vario->keepPair(0, T1, T2, &dist)) continue;

      /* Get the rank of the lag */

      int ipasloc = dirparam.getLagRank(dist);
      if (IFFFF(ipasloc)) continue;
      if (ipas != ipasloc) continue;

      /* Point update */

      indices[0].push_back(iech);
      indices[1].push_back(jech);
      nb++;
    }
  }

  /* Messages */

  if (nb <= 0)
  {
    messerr("No sample found where all variables are defined");
  }
  else
  {
    if (verbose)
    {
      message("Total number of samples = %d\n", nech);
      message("Number of pairs used for translated correlation = %d\n", (int) nb);
    }
  }
  return indices;
}

/****************************************************************************/
/*!
 **  Identify samples from scatter plot when included within a polygon
 **
 ** \return  Error return code
 **
 ** \param[in]  db1          Db descriptor (first variable)
 ** \param[in]  db2          Db descriptor (second variable for flag.same=T)
 ** \param[in]  icol1        Rank of the first column
 ** \param[in]  icol2        Rank of the second column
 ** \param[in]  polygon      Polygons structure
 **
 ** \remarks The two input Db must match exactly (same number of samples with
 ** \remarks same set of coordinates and same optional selection)
 **
 *****************************************************************************/
int correlation_ident(Db *db1, Db *db2, int icol1, int icol2, Polygons *polygon)
{
  if (db1 == nullptr) return (1);
  if (db2 == nullptr) return (1);
  int nech = db1->getSampleNumber();
  int number = 0;

  /* Correlation */

  for (int iech = 0; iech < nech; iech++)
  {
    if (!db1->isActive(iech)) continue;
    double val1 = db1->getArray(iech, icol1);
    if (FFFF(val1)) continue;
    double val2 = db2->getArray(iech, icol2);
    if (FFFF(val2)) continue;

    /* Check of the sample belongs to the polygon */

    VectorDouble coor(3, TEST);
    coor[0] = val1;
    coor[1] = val2;
    if (!polygon->inside(coor, false)) continue;

    /* Print the reference of the sample */

    if (number == 0) mestitle(0, "Samples selected from scatter plot");
    message("Sample #%d - Variable #1=%lf - Variable #2=%lf\n", iech + 1, val1,
            val2);
    number++;
  }

  return (0);
}

/****************************************************************************/
/*!
 **  Evaluate the variogram cloud
 **
 ** \param[in]  db      Db descriptor
 ** \param[in]  dbgrid  Discretization Grid descriptor
 ** \param[in]  iptr    Pointer for the variogram cloud (direction)
 ** \param[in]  varioparam VarioParam structure
 ** \param[in]  idir    Rank of the Direction
 **
 *****************************************************************************/
static void st_variogram_cloud(Db *db,
                               DbGrid *dbgrid,
                               int iptr,
                               const VarioParam *varioparam,
                               int idir)
{
  double dist, value, z1, z2;
  int nech, iech, jech, igrid, ideb;
  SpaceTarget T1(varioparam->getSpace());
  SpaceTarget T2(varioparam->getSpace());

  /* Preliminary calculations */

  nech = db->getSampleNumber();

  // Creating a local Vario structure (to constitute the BiTargetCheck list
  Vario* vario = Vario::create(*varioparam);
  vario->setDb(db);
  if (vario->prepare()) return;

  // Local variables to speed up calculations
  bool hasSel = db->hasLocVariable(ELoc::SEL);
  bool hasWeight = db->hasLocVariable(ELoc::W);
  double w1 = 1.;
  double w2 = 1.;

  /* Loop on the first point */

  for (iech = 0; iech < nech - 1; iech++)
  {
    if (hasSel && !db->isActive(iech)) continue;
    if (hasWeight)
    {
      w1 = db->getWeight(iech);
      if (FFFF(w1)) continue;
    }
    z1 = db->getLocVariable(ELoc::Z, iech, 0);
    if (FFFF(z1)) continue;
    db->getSampleAsST(iech, T1);

    ideb = (varioparam->isDateUsed(db)) ? 0 : iech + 1;
    for (jech = ideb; jech < nech; jech++)
    {
      if (hasSel && !db->isActive(jech)) continue;
      if (hasWeight)
      {
        w2 = db->getWeight(jech);
        if (FFFF(w2)) continue;
      }
      z2 = db->getLocVariable(ELoc::Z, jech, 0);
      if (FFFF(z2)) continue;
      db->getSampleAsST(jech, T2);

      // Reject the point as soon as one BiTargetChecker is not correct
      if (! vario->keepPair(idir, T1, T2, &dist)) continue;

      value = w1 * w2 * (z2 - z1) * (z2 - z1) / 2.;
      igrid = st_update_discretization_grid(dbgrid, dist, value);
      if (igrid < 0) continue;
      dbgrid->updArray(igrid, iptr, 0, 1.);
    }
  }
  return;
}

/****************************************************************************/
/*!
 **  Check the samples which are involved in the pairs which are located
 **  within the polygon
 **
 ** \param[in]  db      Db descriptor
 ** \param[in]  dbgrid  Discretization Grid descriptor
 ** \param[in]  vario   Vario structure
 ** \param[in]  polygon Polygons structure
 **
 *****************************************************************************/
void variogram_cloud_ident(Db *db, DbGrid *dbgrid, Vario *vario, Polygons *polygon)
{
  double dist, z1, z2, value;
  int iech, jech, igrid, idir, ideb;
  VectorDouble coor;
  SpaceTarget T1(vario->getSpace());
  SpaceTarget T2(vario->getSpace());

  /* Initializations */

  int* indg = nullptr;
  int* rank = nullptr;
  double* ids = nullptr;
  const VarioParam &varioparam = vario->getVarioParam();

  // Local variables to speed up calculations
  bool hasSel = db->hasLocVariable(ELoc::SEL);
  bool hasWeight = db->hasLocVariable(ELoc::W);
  double w1 = 1.;
  double w2 = 1.;

  /* Core allocation */

  int nech = db->getSampleNumber();
  indg = db_indg_alloc(dbgrid);
  if (indg == nullptr) goto label_end;
  rank = (int*) mem_alloc(sizeof(int) * nech, 0);
  if (rank == nullptr) goto label_end;
  ids = db_vector_alloc(db);
  if (ids == nullptr) goto label_end;
  for (iech = 0; iech < nech; iech++)
    ids[iech] = 0.;
  coor.resize(dbgrid->getNDim());

  /* Loop on the first point */

  for (iech = 0; iech < nech - 1; iech++)
  {
    if (hasSel && !db->isActive(iech)) continue;
    if (hasWeight)
    {
      w1 = db->getWeight(iech);
     if (FFFF(w1)) continue;
    }
    z1 = db->getLocVariable(ELoc::Z, iech, 0);
    if (FFFF(z1)) continue;
    db->getSampleAsST(iech, T1);

    ideb = (varioparam.isDateUsed(db)) ? 0 : iech + 1;
    for (jech = ideb; jech < nech; jech++)
    {
      if (hasSel && !db->isActive(jech)) continue;
      if (hasWeight)
      {
        w2 = db->getWeight(jech);
        if (FFFF(w2)) continue;
      }
      z2 = db->getLocVariable(ELoc::Z, jech, 0);
      if (FFFF(z2)) continue;
      db->getSampleAsST(jech, T2);

      /* Loop on the directions */

      for (idir = 0; idir < vario->getDirectionNumber(); idir++)
      {
        // Reject the point as soon as one BiTargetChecker is not correct
        if (! vario->keepPair(idir, T1, T2, &dist)) continue;

        value = w1 * w2 * (z2 - z1) * (z2 - z1) / 2.;
        igrid = st_update_discretization_grid(dbgrid, dist, value);
        if (igrid < 0) continue;

        /* Check if the grid cell belongs to the polygon */

        db_index_sample_to_grid(dbgrid, igrid, indg);
        grid_to_point(dbgrid, indg, NULL, coor.data());
        if (!polygon->inside(coor, false)) continue;

        /* Add the references */

        ids[iech] += 1.;
        ids[jech] += 1.;
      }
    }
  }

  /* Printout the scores: they are ranked by decreasing number */

  mestitle(0, "Samples in variogram cloud (by decreasing order of occurence)");
  for (iech = 0; iech < nech; iech++)
    rank[iech] = iech;
  ut_sort_double(0, nech, rank, ids);

  for (iech = 0; iech < nech; iech++)
  {
    jech = nech - iech - 1;
    if (ids[jech] <= 0.) break;
    message("Sample #%3d: %d occurence(s)\n", rank[jech] + 1, (int) ids[jech]);
  }

  label_end: indg = db_indg_free(indg);
  ids = db_vector_free(ids);
  rank = (int*) mem_free((char* ) rank);
  return;
}

/****************************************************************************/
/*!
 **  Evaluate the variogram cloud
 **
 ** \param[in]  db     Db descriptor
 ** \param[in]  varioparam VarioParam structure
 ** \param[in]  idir   Rank of the Direction
 **
 ** \param[out] vmax   Maximum variogram value
 **
 *****************************************************************************/
static void st_variogram_cloud_dim(Db *db,
                                   const VarioParam *varioparam,
                                   int idir,
                                   double *vmax)
{
  double dist = 0;
  SpaceTarget T1(varioparam->getSpace());
  SpaceTarget T2(varioparam->getSpace());

  /* Preliminary calculations */

  const DirParam &dirparam = varioparam->getDirParam(idir);
  int nech = db->getSampleNumber();

  // Creating a local Vario structure (to constitute the BiTargetCheck list
  Vario* vario = Vario::create(*varioparam);
  vario->setDb(db);
  if (vario->prepare()) return;

  // Local variables to speed up calculations
  bool hasSel = db->hasLocVariable(ELoc::SEL);
  bool hasWeight = db->hasLocVariable(ELoc::W);
  bool hasDate = varioparam->isDateUsed(db);

  /* Loop on the first point */

  for (int iech = 0; iech < nech - 1; iech++)
  {
    if (hasSel && !db->isActive(iech)) continue;
    if (hasWeight && FFFF(db->getWeight(iech))) continue;
    db->getSampleAsST(iech, T1);

    int ideb = (hasDate) ? 0 : iech + 1;
    for (int jech = ideb; jech < nech; jech++)
    {
      if (hasSel && !db->isActive(jech)) continue;
      if (hasWeight && FFFF(db->getWeight(jech))) continue;
      db->getSampleAsST(jech, T2);

      // Reject the point as soon as one BiTargetChecker is not correct
      if (! vario->keepPair(idir, T1, T2, &dist)) continue;

      if (floor(dist / dirparam.getDPas() + 0.5) >= dirparam.getLagNumber())
        continue;

      double w1 = db->getWeight(iech);
      double w2 = db->getWeight(jech);
      if (FFFF(w1) || FFFF(w2)) continue;
      double z1 = db->getLocVariable(ELoc::Z, iech, 0);
      double z2 = db->getLocVariable(ELoc::Z, jech, 0);
      if (FFFF(z1) || FFFF(z2)) continue;
      double value = w1 * w2 * (z2 - z1) * (z2 - z1) / 2.;
      if (value > (*vmax)) (*vmax) = value;
    }
  }
  return;
}

/****************************************************************************/
/*!
 **  Evaluate the experimental variogram cloud on irregular data
 **
 ** \return  Error return code
 **
 ** \param[in]  db           Db descriptor
 ** \param[in]  varioparam   VarioParam structure
 ** \param[in]  dbgrid       Output grid for storing the variogram cloud
 ** \param[in]  namconv      Naming convention
 **
 *****************************************************************************/
int variogram_cloud(Db *db,
                    const VarioParam *varioparam,
                    DbGrid *dbgrid,
                    const NamingConvention& namconv)
{
  if (db == nullptr) return (1);
  if (dbgrid == nullptr) return (1);
  if (varioparam == (VarioParam*) NULL) return (1);

  /* Preliminary checks */

  if (db->getNDim() != varioparam->getDimensionNumber())
  {
    messerr("Inconsistent parameters:");
    messerr("Data Base: NDIM=%d", db->getNDim());
    messerr("Variogram: NDIM=%d", varioparam->getDimensionNumber());
    return (1);
  }
  if (!db->isVariableNumberComparedTo(1)) return 1;
  if (dbgrid->getNDim() != 2)
  {
    messerr("The output Db for storing the variogram cloud must be 2-D");
    return (1);
  }

  /* Allocate new variables */

  int ndir = varioparam->getDirectionNumber();
  int iptr = dbgrid->addColumnsByConstant(ndir, 0.);
  if (iptr < 0) return (1);

  /* Loop on the directions to evaluate */

  for (int idir = 0; idir < ndir; idir++)
  {
    st_variogram_cloud(db, dbgrid, iptr + idir, varioparam, idir);

    /* Convert zero values into TEST */

    st_final_discretization_grid(dbgrid, iptr + idir);
  }

  // Naming of the newly created variables

  namconv.setNamesAndLocators(db, VectorString(), ELoc::Z, -1, dbgrid, iptr, String(), ndir, false);

  return (0);
}

/****************************************************************************/
/*!
 **  Evaluate the bounds for the experimental variogram cloud
 **
 ** \return  Error return code
 **
 ** \param[in]  db           Db descriptor
 ** \param[in]  varioparam   VarioParam structure
 **
 ** \param[out] vmax         Maximum variogram value
 **
 *****************************************************************************/
int variogram_cloud_dim(Db *db, const VarioParam *varioparam, double *vmax)
{
  if (db == nullptr) return (1);
  if (varioparam == (VarioParam*) NULL) return (1);

  /* Preliminary checks */

  if (db->getNDim() != varioparam->getDimensionNumber())
  {
    messerr("Inconsistent parameters:");
    messerr("Data Base: NDIM=%d", db->getNDim());
    messerr("Variogram: NDIM=%d", varioparam->getDimensionNumber());
    return (1);
  }
  if (!db->isVariableNumberComparedTo(1)) return 1;

  /* Loop on the directions to evaluate */

  *vmax = 0.;
  for (int idir = 0; idir < varioparam->getDirectionNumber(); idir++)
    st_variogram_cloud_dim(db, varioparam, idir, vmax);

  return (0);
}

/****************************************************************************/
/*!
 **  Ask the characteristics of the Vario structure
 **
 ** \return  Error returned code
 **
 ** \param[in]  vario  Vario structure
 **
 ** \param[out]  calcul_type Type of calculation (ECalcVario)
 ** \param[out]  ndim        Space dimension
 ** \param[out]  nvar        Number of variables
 ** \param[out]  ndir        Number of calculation directions
 ** \param[out]  ndate       Number of Date Intervals
 ** \param[out]  scale       Scaling factor for the transitive covariogram
 ** \param[out]  dates       Array of bounds for Date Intervals
 **
 ** \remark The array 'dates' must be freed by calling function.
 ** \remark The following code shows how to extract the calculation results
 ** \remark from a variogram
 **
 *****************************************************************************/
int vario_extract(Vario *vario,
                  ECalcVario *calcul_type,
                  int *ndim,
                  int *nvar,
                  int *ndir,
                  int *ndate,
                  double *scale,
                  double **dates)
{
  double *date_loc;

  /* Returning arguments */

  *calcul_type = vario->getCalcul();
  *ndim = vario->getDimensionNumber();
  *nvar = vario->getVariableNumber();
  *ndir = vario->getDirectionNumber();
  *ndate = vario->getDateNumber();
  *scale = vario->getScale();
  date_loc = (double*) mem_alloc(sizeof(double) * vario->getDateNumber() * 2, 1);
  int ecr = 0;
  for (int i = 0; i < vario->getDateNumber(); i++)
    for (int icas = 0; icas < 2; icas++)
      date_loc[ecr++] = vario->getDates(i, icas);
  *dates = date_loc;

  return (0);
}

/****************************************************************************/
/*!
 **  Ask for the rank of the 'vardir' structure, given direction and date
 **
 ** \return  Absolute rank (or -1 for error)
 **
 ** \param[in]  vario  Vario structure
 ** \param[in]  idir   Rank for the direction (starting from 0)
 ** \param[in]  idate  Rank for the Date (starting from 0)
 **
 ** \remark  An error occurs if 'idir' is negative or larger than 'ndir'
 ** \remark  or if 'idate' is negative or larger than 'ndate'
 **
 *****************************************************************************/
int vario_get_rank(Vario *vario, int idir, int idate)
{
  int rank = idir;
  int ndir = vario->getDirectionNumber();
  int ndate = vario->getDateNumber();
  if (idir < 0 || idir >= ndir) return (-1);
  if (ndate > 0)
  {
    if (idate < 0 || idate >= ndate) return (-1);
    rank = rank * idir + idate;
  }
  return (rank);
}

/****************************************************************************/
/*!
 **  Copy a direction from one Vario to another Vario
 **
 ** \param[in]  vario_in     Input Vario structure
 ** \param[in]  idir_in      Rank of the Input Direction
 ** \param[in]  vario_out    Output Vario structure
 ** \param[in]  idir_out     Rank of the Output Direction
 **
 *****************************************************************************/
void vardir_copy(VarioParam *vario_in,
                 int idir_in,
                 VarioParam *vario_out,
                 int idir_out)
{
  if (vario_in == (VarioParam*) NULL) return;
  if (idir_in < 0 || idir_in >= vario_in->getDirectionNumber()) return;
  if (vario_out == (VarioParam*) NULL) return;
  if (idir_out < 0 || idir_out >= vario_in->getDirectionNumber()) return;

  DirParam dir_in = vario_in->getDirParam(idir_in);
  DirParam dir_out = vario_out->getDirParam(idir_out);
  dir_out = dir_in;
}

/****************************************************************************/
/*!
 **  Linear interpolation
 **
 ** \return  Interpolated value
 **
 ** \param[in]  n      Number of discretization steps
 ** \param[in]  x      Discretized X (sorted increasingly)
 ** \param[in]  y      Discretized Y
 ** \param[in]  x0     Origin
 **
 *****************************************************************************/
static double st_linear_interpolate(int n, double *x, double *y, double x0)
{
  int i;

  if (x0 < x[0]) return (y[0]);
  if (x0 > x[n - 1]) return (y[n - 1]);
  for (i = 1; i < n; i++)
  {
    if (x0 < x[i - 1]) continue;
    if (x0 > x[i]) continue;
    return (y[i - 1] + (y[i] - y[i - 1]) * (x0 - x[i - 1]) / (x[i] - x[i - 1]));
  }
  return (TEST);
}

/****************************************************************************/
/*!
 **  Calculate the experimental variogram of the completed variable starting
 **  from the experimental variogram of the truncated variable
 **
 ** \param[in,out] vario  Vario structure
 ** \param[in]  nh     Number of Hermite polynomials
 ** \param[in]  ycut   Truncation (lowest) value
 **
 *****************************************************************************/
void variogram_trans_cut(Vario *vario, int nh, double ycut)
{
  double variance, sum, cyp, cyy;
  int ih, idisc, ndisc, idir, ipas;
  static double disc = 0.01;

  /* Initializations */

  ndisc = (int) (2. / disc + 1.);

  /* Core allocation */

  VectorDouble ro(ndisc);
  VectorDouble covyp(ndisc);
  for (idisc = 0; idisc < ndisc; idisc++)
    ro[idisc] = disc * idisc - 1.;

  /* Calculate the first normalized Hermite polynomials for ycut */

  VectorDouble psic = hermiteCoefLower(ycut, nh);

  /* Variance */

  variance = 0.;
  for (ih = 1; ih < nh; ih++)
    variance += psic[ih] * psic[ih];

  for (idisc = 0; idisc < ndisc; idisc++)
  {
    sum = 0.;
    for (ih = 1; ih < nh; ih++)
      sum += psic[ih] * psic[ih] * pow(ro[idisc], ih);
    covyp[idisc] = sum;
  }

  /* Loop on the directions */

  for (idir = 0; idir < vario->getDirectionNumber(); idir++)
  {

    /* Loop on the lags */

    for (ipas = 0; ipas < vario->getLagNumber(idir); ipas++)
    {
      cyp = variance - vario->getGg(idir, 0, 0, ipas);
      cyy = st_linear_interpolate(ndisc, covyp.data(), ro.data(), cyp);
      vario->setGg(idir, 0, 0, ipas, MAX(0, 1. - cyy));
    }
  }

  /* Set the variance */

  vario->setVar(1., 0, 0);
}

/****************************************************************************/
/*!
 **  Determine the samples used for a variogram in multilayers framework
 **
 ** \return  Error return code
 **
 ** \param[in]  db     Db description
 ** \param[in]  seltab Number of sample definition (0, 1 or 2)
 ** \param[in]  vario  Vario structure
 **
 ** \param[out]  vorder Vario_Order struct
 ure
 **
 *****************************************************************************/
int variogram_mlayers(Db *db, int *seltab, Vario *vario, Vario_Order *vorder)
{
  int iiech, iech, jjech, jech, ipas, npair;
  SpaceTarget T1(vario->getSpace());
  SpaceTarget T2(vario->getSpace());

  /* Initializations */

  if (db == nullptr) return 1;
  if (vario == nullptr) return 1;

  // Local variables to speed up calculations
  bool hasSel = db->hasLocVariable(ELoc::SEL);
  int nech = db->getSampleNumber();
  double dist = 0.;

  /* Loop on the directions */

  for (int idir = 0; idir < vario->getDirectionNumber(); idir++)
  {
    const DirParam &dirparam = vario->getDirParam(idir);


    /* Loop on the first point */

    for (iech = iiech = 0; iech < nech; iech++)
    {
      if (hasSel && !db->isActive(iech)) continue;
      db->getSampleAsST(iech, T1);

      if (seltab[iech] == 0) continue;
      for (int ifois = 0; ifois < seltab[iech]; ifois++, iiech++)
      {
        for (jech = jjech = 0; jech < nech; jech++)
        {
          if (hasSel && !db->isActive(jech)) continue;
          if (seltab[jech] == 0) continue;
          db->getSampleAsST(jech, T2);

          for (int jfois = 0; jfois < seltab[jech]; jfois++, jjech++)
          {

            // Reject the point as soon as one BiTargetChecker is not correct
            if (! vario->keepPair(idir, T1, T2, &dist)) continue;

            /* Get the rank of the lag */

            ipas = dirparam.getLagRank(dist);
            if (IFFFF(ipas)) continue;

            /* Internal storage */

            vario_order_add(vorder, iiech, jjech, &iech, &jech, ipas, idir,
                            ABS(dist));
          }
        }
      }
    }
  }
  vorder = vario_order_final(vorder, &npair);
  return (0);
}

/*****************************************************************************/
/*!
 **  Calculate the experimental variogram of the Raw starting from the Model
 **  of the Gaussian variable
 **
 ** \return  Error return code
 **
 ** \param[in,out] vario    Experimental variogram
 ** \param[in]  anam        Point anamorphosis
 ** \param[in]  model       Model of the Punctual Gaussian
 **
 ** \remark  At entrance, the input variogram only serves in providing
 ** \remark  the calculation parameters
 **
 *****************************************************************************/
int variogram_y2z(Vario *vario, AAnam *anam, Model *model)
{
  double chh, cov_value;

  /* Preliminary checks */

  if (vario == nullptr) return 1;
  if (anam == (AAnam*) NULL) return 1;
  if (model == nullptr) return 1;
  if (anam->getType() != EAnam::HERMITIAN)
  {
    messerr("This function is restricted to Gaussian Anamorphosis");
    return 1;
  }
  AnamHermite *anam_hermite = dynamic_cast<AnamHermite*>(anam);
  if (anam_hermite->getRCoef() != 1.)
  {
    messerr("This function is restricted to Punctual Anamorphosis");
    return 1;
  }
  if (vario == nullptr) return 1;
  if (vario->getVariableNumber() != 1)
  {
    messerr("This function is restricted to Monovariate Variogram");
    return 1;
  }
  if (model->getVariableNumber() != 1)
  {
    messerr("This function requires a Monovariate Model");
    return 1;
  }
  if (model->getDimensionNumber() != vario->getDimensionNumber())
  {
    messerr("Variogram and Model should share the same Space Dimension");
    return 1;
  }

  /* Initializations */

  int ndim = vario->getDimensionNumber();
  VectorDouble d1(ndim, 0.);

  /* Calculate the theoretical variance of Z */

  double varz = anam_hermite->computeVariance(1.);

  /* Loop on the directions of the variogram */

  for (int idir = 0, ndir = vario->getDirectionNumber(); idir < ndir; idir++)
  {
    /* Loop on the lags */

    for (int ipas = 0, npas = vario->getLagNumber(idir); ipas < npas; ipas++)
    {
      for (int idim = 0; idim < ndim; idim++)
        d1[idim] = (ipas + 1) * vario->getDPas(idir)
                   * vario->getCodir(idir, idim);

      model_calcul_cov(NULL,model, nullptr, 1, 1., d1, &chh);
      if (chh < 0.)
      {
        messerr("Gaussian covariance is negative in direction %d for lag %d",
                idir + 1, ipas + 1);
        messerr("Calculation is impossible");
        return 1;
      }

      cov_value = anam_hermite->computeVariance(chh);
      vario->setGg(idir, 0, 0, ipas, varz - cov_value);
      vario->setHh(idir, 0, 0, ipas, (ipas + 1) * vario->getDPas(idir));
      vario->setSw(idir, 0, 0, ipas, 1.);
    }
  }
  return 0;
}

/****************************************************************************/
/*!
 **  Evaluate the experimental conditional expectation
 **
 ** \param[in]  db1           Db descriptor (for target variable)
 ** \param[in]  db2           Db descriptor (for auxiliary variables)
 ** \param[in]  icol1         Rank of the target variable
 ** \param[in]  icol2         Rank of the explanatory variable
 ** \param[in]  mini          Minimum value for the explanaroty variable
 ** \param[in]  maxi          Maximum value for the explanaroty variable
 ** \param[in]  nclass        Number of classes
 ** \param[in]  verbose       Verbose flag
 **
 *****************************************************************************/
VectorVectorDouble condexp(Db *db1,
                           Db *db2,
                           int icol1,
                           int icol2,
                           double mini,
                           double maxi,
                           int nclass,
                           bool verbose)
{
  VectorVectorDouble xycond(2);
  xycond[0].resize(nclass);
  xycond[1].resize(nclass);
  VectorInt ncond(nclass,0);

  /* Loop on the samples */

  for (int iech = 0; iech < db1->getSampleNumber(); iech++)
  {
    if (!db1->isActive(iech)) continue;
    double val1 = db1->getArray(iech, icol1);
    if (FFFF(val1)) continue;
    double val2 = db2->getArray(iech, icol2);
    if (FFFF(val2)) continue;
    if (val2 < mini || val2 > maxi) continue;

    int rank = int((nclass - 1.) * (val2 - mini) / (maxi - mini));

    xycond[0][rank] += val1;
    xycond[1][rank] += val2;
    ncond[rank]++;
  }

  /* Normation */

  for (int iclass = 0; iclass < nclass; iclass++)
  {
    if (ncond[iclass] <= 0)
    {
      xycond[0][iclass] = TEST;
      xycond[1][iclass] = TEST;
    }
    else
    {
      xycond[0][iclass] /= ncond[iclass];
      xycond[1][iclass] /= ncond[iclass];
    }
  }

  /* Optional printout */

  if (verbose)
  {
    message("Experimental Conditional Expectation\n");
    for (int iclass = 0; iclass < nclass; iclass++)
    {
      if (ncond[iclass] > 0)
        message("Class %2d : V1=%lf V2=%lf\n", iclass + 1, xycond[0][iclass],
                xycond[1][iclass]);
    }
  }
  return xycond;
}

/****************************************************************************/
/*!
 **  Evaluate the experimental variogram cloud
 **
 ** \return  Error return code
 **
 ** \param[in]  db           Db descriptor
 ** \param[in]  varioparam   VarioParam structure
 ** \param[in]  lagmax       Maximum distance
 ** \param[in]  varmax       Maximum Variance value
 ** \param[in]  lagnb        Number of discretization steps along distance axis
 ** \param[in]  varnb        Number of discretization steps along variance axis
 ** \param[in]  namconv      Naming convention
 **
 *****************************************************************************/
DbGrid* db_variogram_cloud(Db *db,
                           const VarioParam *varioparam,
                           double lagmax,
                           double varmax,
                           int lagnb,
                           int varnb,
                           const NamingConvention& namconv)
{
  if (FFFF(lagmax)) lagmax = db->getExtensionDiagonal();
  if (FFFF(varmax)) (void) variogram_cloud_dim(db, varioparam, &varmax);

  // Create a grid as a support for the variogram cloud calculations

  VectorInt nx(2);
  nx[0] = lagnb;
  nx[1] = varnb;
  VectorDouble dx(2);
  dx[0] = lagmax / (double) lagnb;
  dx[1] = varmax / (double) varnb;
  VectorDouble x0(2);
  x0[0] = 0.;
  x0[1] = 0.;
  DbGrid *dbgrid = DbGrid::create(nx, dx, x0);

  // Calling the variogram cloud calculation function

  int error = variogram_cloud(db, varioparam, dbgrid, namconv);

  // In case of error, free the newly created structure

  if (error)
  {
    delete dbgrid;
    dbgrid = nullptr;
  }
  return dbgrid;
}

