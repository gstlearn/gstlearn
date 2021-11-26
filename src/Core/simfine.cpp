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
#include "Basic/Law.hpp"
#include "geoslib_e.h"
#include "geoslib_old_f.h"

static int NDIM, NMULT, FLAG_KS, IX[2][5], IY[2][5], IZ[2][5];
static double XN[2][5], YN[2][5], ZN[2][5], WGT[2][2][5], STDV[2][2];
static VectorInt NX1(3), NX2(3);
static VectorDouble DX1(3), DX2(3), X01(3), X02(3);

/*! \cond */
#define LHS(i,j) (lhs[(i) * neq + (j)])
#define RHS(i)   (rhs[(i)])
/*! \endcond */

/****************************************************************************/
/*!
 **  Define the characteristics from Db1 to Db2
 **
 ** \param[in]  db  Staring grid Db structure
 **
 *****************************************************************************/
static void st_dim_1_to_2(Db *db)

{

  /* Input file */

  NX1[0] = (NDIM >= 1) ? db->getNX(0) :
                         1;
  NX1[1] = (NDIM >= 2) ? db->getNX(1) :
                         1;
  NX1[2] = (NDIM >= 3) ? db->getNX(2) :
                         1;
  DX1[0] = (NDIM >= 1) ? db->getDX(0) :
                         1.;
  DX1[1] = (NDIM >= 2) ? db->getDX(1) :
                         1.;
  DX1[2] = (NDIM >= 3) ? db->getDX(2) :
                         1.;
  X01[0] = (NDIM >= 1) ? db->getX0(0) :
                         0.;
  X01[1] = (NDIM >= 2) ? db->getX0(1) :
                         0.;
  X01[2] = (NDIM >= 3) ? db->getX0(2) :
                         0.;

  /* Output file */

  NX2[0] = NX1[0] * 2 + 1;
  NX2[1] = NX1[1] * 2 + 1;
  NX2[2] = NX1[2];
  DX2[0] = DX1[0] / 2.;
  DX2[1] = DX1[1] / 2.;
  DX2[2] = DX1[2];
  X02[0] = X01[0] - DX2[0];
  X02[1] = X01[1] - DX2[1];
  X02[2] = X01[2];
}

/****************************************************************************/
/*!
 **  Define the characteristics from Db2 to Db1
 **
 ** \param[in]  db  Starting grid Db structure
 **
 *****************************************************************************/
static void st_dim_2_to_1(Db *db)

{

  /* Input file */

  NX2[0] = (NDIM >= 1) ? db->getNX(0) :
                         1;
  NX2[1] = (NDIM >= 2) ? db->getNX(1) :
                         1;
  NX2[2] = (NDIM >= 3) ? db->getNX(2) :
                         1;
  DX2[0] = (NDIM >= 1) ? db->getDX(0) :
                         1.;
  DX2[1] = (NDIM >= 2) ? db->getDX(1) :
                         1.;
  DX2[2] = (NDIM >= 3) ? db->getDX(2) :
                         1.;
  X02[0] = (NDIM >= 1) ? db->getX0(0) :
                         0.;
  X02[1] = (NDIM >= 2) ? db->getX0(1) :
                         0.;
  X02[2] = (NDIM >= 3) ? db->getX0(2) :
                         0.;

  /* Output file */

  NX1[0] = NX2[0] - 2;
  NX1[1] = NX2[1] - 2;
  NX1[2] = NX2[2];
  DX1[0] = DX2[0];
  DX1[1] = DX2[1];
  DX1[2] = DX2[2];
  X01[0] = X02[0] + DX2[0];
  X01[1] = X02[1] + DX2[1];
  X01[2] = X02[2];
}

/****************************************************************************/
/*!
 **  Read a value in a Db
 **
 ** \return  The value read
 **
 ** \param[in]  db     Grid Db structure
 ** \param[in]  iatt   Rank of the attribute to be read
 ** \param[in]  ix0    Index of the target along X
 ** \param[in]  iy0    Index of the target along Y
 ** \param[in]  iz0    Index of the target along Z
 ** \param[in]  idx    Shift along X
 ** \param[in]  idy    Shift along Y
 ** \param[in]  idz    Shift along Z
 **
 *****************************************************************************/
static double st_read(Db *db,
                      int iatt,
                      int ix0,
                      int iy0,
                      int iz0,
                      int idx,
                      int idy,
                      int idz)
{
  int ix, iy, iz, iad, ind[3];
  double value;

  ix = ix0 + idx;
  if (ix < 0 || ix >= db->getNX(0)) ix = ix0 - idx;
  iy = iy0 + idy;
  if (iy < 0 || iy >= db->getNX(1)) iy = iy0 - idy;
  iz = iz0 + idz;
  if (iz < 0 || iz >= db->getNX(2)) iz = iz0 - idz;
  ind[0] = ix;
  ind[1] = iy;
  ind[2] = iz;
  iad = db_index_grid_to_sample(db, ind);
  value = db->getArray(iad, iatt);
  return (value);
}

/****************************************************************************/
/*!
 **  Write a value in a Db
 **
 ** \param[in]  db     Grid Db structure
 ** \param[in]  iatt   Rank of the attribute to be written
 ** \param[in]  ix0    Index of the target along X
 ** \param[in]  iy0    Index of the target along Y
 ** \param[in]  iz0    Index of the target along Z
 ** \param[in]  value  Value to be written
 **
 *****************************************************************************/
static void st_write(Db *db, int iatt, int ix0, int iy0, int iz0, double value)
{
  int iad, ind[3];

  ind[0] = ix0;
  ind[1] = iy0;
  ind[2] = iz0;
  iad = db_index_grid_to_sample(db, ind);
  db->setArray(iad, iatt, value);
}

/****************************************************************************/
/*!
 **  Copy the initial information from db1 into db2
 **
 ** \param[in]  db1     Input grid Db structure
 ** \param[in]  iatt1   Rank of the attribute to be read from db1
 ** \param[in]  db2     Output grid Db structure
 ** \param[in]  iatt2   Rank of the attribute to be written into db2
 **
 *****************************************************************************/
static void st_merge_data(Db *db1, int iatt1, Db *db2, int iatt2)
{
  int ix1, iy1, iz1, ix2, iy2, iz2;
  double value;

  for (ix1 = 0; ix1 < NX1[0]; ix1++)
    for (iy1 = 0; iy1 < NX1[1]; iy1++)
      for (iz1 = 0; iz1 < NX1[2]; iz1++)
      {
        ix2 = 1 + 2 * ix1;
        iy2 = 1 + 2 * iy1;
        iz2 = iz1;
        value = st_read(db1, iatt1, ix1, iy1, iz1, 0, 0, 0);
        st_write(db2, iatt2, ix2, iy2, iz2, value);
      }
}

/****************************************************************************/
/*!
 **  Truncate the resulting information from db2 into db1
 **
 ** \param[in]  db2     Input grid Db structure
 ** \param[in]  iatt2   Rank of the attribute to be read from db2
 ** \param[in]  db1     Output grid Db structure
 ** \param[in]  iatt1   Rank of the attribute to be written into db1
 **
 *****************************************************************************/
static void st_truncate_result(Db *db2, int iatt2, Db *db1, int iatt1)
{
  int ix, iy, iz;
  double value;

  for (ix = 0; ix < NX1[0]; ix++)
    for (iy = 0; iy < NX1[1]; iy++)
      for (iz = 0; iz < NX1[2]; iz++)
      {
        value = st_read(db2, iatt2, ix, iy, iz, 1, 1, 0);
        st_write(db1, iatt1, ix, iy, iz, value);
      }
}

/****************************************************************************/
/*!
 **  Define the location of a neighborhood point
 **
 ** \param[in]  type   Type of kriging
 ** \param[in]  rank   Rank of the neighboring data
 ** \param[in]  idx    Shift along X
 ** \param[in]  idy    Shift along Y
 ** \param[in]  idz    Shift along Z
 **
 *****************************************************************************/
static void st_neigh(int type, int rank, int idx, int idy, int idz)
{
  IX[type][rank] = idx;
  IY[type][rank] = idy;
  IZ[type][rank] = idz;
  XN[type][rank] = idx * DX2[0];
  YN[type][rank] = idy * DX2[1];
  ZN[type][rank] = idz * DX2[2];
}

/****************************************************************************/
/*!
 **  Simulate the target cell
 **
 ** \param[in]  db     Grid Db structure
 ** \param[in]  type   Type of kriging
 ** \param[in]  iatt   Rank of the attribute to be read
 ** \param[in]  ix0    Index of the target along X
 ** \param[in]  iy0    Index of the target along Y
 ** \param[in]  iz0    Index of the target along Z
 **
 *****************************************************************************/
static void st_simulate_target(Db *db,
                               int type,
                               int iatt,
                               int ix0,
                               int iy0,
                               int iz0)
{
  double value;
  int i;

  value = 0.;
  if (iz0 == 0)
  {

    /* Case of the first layer */

    for (i = 0; i < 4; i++)
      value += (WGT[type][0][i]
          * st_read(db, iatt, ix0, iy0, iz0, IX[type][i], IY[type][i],
                    IZ[type][i]));
    value += STDV[type][0] * law_gaussian();
  }
  else
  {

    /* Case of a subsequent layer */

    for (i = 0; i < 5; i++)
      value += (WGT[type][1][i]
          * st_read(db, iatt, ix0, iy0, iz0, IX[type][i], IY[type][i],
                    IZ[type][i]));
    value += STDV[type][1] * law_gaussian();
  }

  st_write(db, iatt, ix0, iy0, iz0, value);
}

/****************************************************************************/
/*!
 **  Solve the kriging system
 **
 ** \return  Error return code
 **
 ** \param[in]  type   Type of kriging
 ** \param[in]  rank   Rank of the neighboring data
 ** \param[in]  nb     Number of equations
 ** \param[in]  model  Model structure
 **
 *****************************************************************************/
static int st_kriging_solve(int type, int rank, int nb, Model *model)
{
  double lhs[36], rhs[6], var[2], variance;
  int i, j, neq;
  VectorDouble d1;
  CovCalcMode mode;

  /* Initializations */

  neq = (FLAG_KS) ? nb :
                    nb + 1;
  d1.resize(3);

  /* Establish the kriging L.H.S. */

  for (i = 0; i < nb; i++)
    for (j = 0; j < nb; j++)
    {
      d1[0] = XN[type][i] - XN[type][j];
      d1[1] = YN[type][i] - YN[type][j];
      d1[2] = ZN[type][i] - ZN[type][j];
      model_calcul_cov(model, mode, 1, 1., d1, &LHS(i, j));
    }

  /* Establish the kriging R.H.S. */

  mode.setMember(ECalcMember::RHS);
  for (i = 0; i < nb; i++)
  {
    d1[0] = XN[type][i];
    d1[1] = YN[type][i];
    d1[2] = ZN[type][i];
    model_calcul_cov(model, mode, 1, 1., d1, &RHS(i));
  }

  /* Add the Universality condition (optional) */

  if (!FLAG_KS)
  {
    for (i = 0; i < nb; i++)
      LHS(i,nb) = LHS(nb,i) = 1.;
    LHS(nb,nb) = 0;
    RHS(nb) = 1.;
  }

  /* Derive the kriging weights */

  if (matrix_invert(lhs, neq, -1))
  {
    messerr("Kriging matrix inversion failed");
    messerr("Check the consistency between the model and the SK/OK option");
    return (1);
  }
  matrix_product(neq, neq, 1, lhs, rhs, WGT[type][rank]);

  /* Calculate the variance */

  mode.setMember(ECalcMember::VAR);
  for (i = 0; i < 3; i++)
    d1[i] = 0.;
  model_calcul_cov(model, mode, 1, 1., d1, &var[0]);
  matrix_product(1, neq, 1, rhs, WGT[type][rank], &var[1]);
  variance = var[0] - var[1];
  STDV[type][rank] = (variance > 0) ? sqrt(variance) :
                                      0.;

  /* Printout of the weights */

  if (debug_query("kriging"))
  {
    message("\nDisplay of the Kriging weights\n");
    for (i = 0; i < nb; i++)
      message("X=%10.3lf Y=%10.3lf Z=%10.3lf W=%10.6lf\n", XN[type][i],
              YN[type][i], ZN[type][i], WGT[type][rank][i]);
    message("Variance of error           = %10.6lf\n", variance);
    message("Standard deviation of error = %10.6lf\n", STDV[type][rank]);
  }

  return (0);
}

/****************************************************************************/
/*!
 **  Establish and solve the different kriging systems
 **
 ** \return  Eror return code
 **
 ** \param[in]  model  Model structure
 **
 *****************************************************************************/
static int st_kriging_define(Model *model)

{
  /* Define the kriging system for the cell centers */

  st_neigh(0, 0, -1, -1, 0);
  st_neigh(0, 1, 1, -1, 0);
  st_neigh(0, 2, 1, 1, 0);
  st_neigh(0, 3, -1, 1, 0);
  st_neigh(0, 4, 0, 0, -1);

  if (st_kriging_solve(0, 0, 4, model)) return (1);
  if (st_kriging_solve(0, 1, 5, model)) return (1);

  /* Define the kriging system for the mid-vertices */

  st_neigh(1, 0, -1, 0, 0);
  st_neigh(1, 1, 0, -1, 0);
  st_neigh(1, 2, 1, 0, 0);
  st_neigh(1, 3, 0, 1, 0);
  st_neigh(1, 4, 0, 0, -1);

  if (st_kriging_solve(1, 0, 4, model)) return (1);
  if (st_kriging_solve(1, 1, 5, model)) return (1);

  return (0);
}

/****************************************************************************/
/*!
 **  Simulate the missing nodes
 **
 ** \param[in]  db     Grid Db structure
 ** \param[in]  iatt   Rank of the column
 **
 *****************************************************************************/
static void st_simulate_nodes(Db *db, int iatt)
{
  int ix, iy, iz;

  /* Perform the cell centers */

  for (iz = 0; iz < NX2[2]; iz++)
    for (ix = 0; ix < NX2[0]; ix++)
      for (iy = 0; iy < NX2[1]; iy++)
        if ((ix % 2 == 0) && (iy % 2 == 0))
          st_simulate_target(db, 0, iatt, ix, iy, iz);

  /* Perform the cell mid-vertices */

  for (iz = 0; iz < NX2[2]; iz++)
    for (ix = 0; ix < NX2[0]; ix++)
      for (iy = 0; iy < NX2[1]; iy++)
        if (((ix % 2 == 0) && (iy % 2 == 1)) || ((ix % 2 == 1) && (iy % 2 == 0)))
          st_simulate_target(db, 1, iatt, ix, iy, iz);
}

/****************************************************************************/
/*!
 **  Refine the simulation
 **
 ** \return  Error return code
 **
 ** \param[in]  dbin       Input grid Db structure
 ** \param[in]  model      Model srtucture
 ** \param[in]  flag_ks    1 for SK and 0 for OK
 ** \param[in]  nmult      Refinement factor
 ** \param[in]  seed       Seed for the random number generator
 **
 ** \param[out] tab        Output array
 **
 *****************************************************************************/
GSTLEARN_EXPORT int simfine_f(Db *dbin,
                              Model *model,
                              int flag_ks,
                              int nmult,
                              int seed,
                              VectorDouble &tab)
{
  int error, imult, iatt1, iatt2, idim;
  double diag;
  Db *db1, *db2;

  /* Initializations */

  error = 1;
  db1 = db2 = nullptr;
  NDIM = dbin->getNDim();
  NMULT = nmult;
  db1 = dbin;
  FLAG_KS = flag_ks;
  law_set_random_seed(seed);

  /* Preliminary check */

  if (!is_grid(dbin))
  {
    messerr("This simulation refinement is dedicated to Grid File");
    goto label_end;
  }
  if (!dbin->isVariableNumberComparedTo(1)) goto label_end;

  /* Patch the model with maximum dimension for OK */

  diag = 0.;
  for (idim = 0; idim < NDIM; idim++)
    diag += dbin->getDX(idim) * dbin->getDX(idim);
  model->setField(sqrt(diag));

  /* Store information from the input grid */

  iatt1 = db_attribute_identify(dbin, ELoc::Z, 0);
  if (iatt1 <= 0) goto label_end;

  /* Loop on the refinement factors */

  for (imult = 0; imult < NMULT; imult++)
  {

    /* Create the output grid */

    st_dim_1_to_2(db1);
    db2 = db_create_grid(0, NDIM, 1, ELoadBy::SAMPLE, 1, NX2, X02, DX2,
                         dbin->getGrid().getRotAngles());
    iatt2 = db2->addFields(1, TEST);
    if (iatt2 <= 0) goto label_end;

    /* Establish the kriging system */

    if (st_kriging_define(model)) goto label_end;

    /* Copy the initial data */

    st_merge_data(db1, iatt1, db2, iatt2);

    /* Perform the simulation */

    st_simulate_nodes(db2, iatt2);

    /* Create the new input file (for next step) */

    if (db1 != dbin) db1 = db_delete(db1);
    st_dim_2_to_1(db2);
    db1 = db_create_grid(0, NDIM, 1, ELoadBy::SAMPLE, 1, NX1, X01, DX1,
                         dbin->getGrid().getRotAngles());
    iatt1 = db1->addFields(1, TEST);
    if (iatt1 <= 0) goto label_end;

    /* Truncate the output grid for next step */

    st_truncate_result(db2, iatt2, db1, iatt1);

    /* Delete the output file */

    db2 = db_delete(db2);
  }

  /* Load the resulting variable in the output array */

  db_vector_get_att(db1, iatt1, tab.data());

  /* Set the error returned code */

  error = 0;

  label_end: if (db1 != dbin) db1 = db_delete(db1);
  return (error);
}

/****************************************************************************/
/*!
 **  Define the characteristics of the refined grid
 **
 ** \return  Error returned code
 **
 ** \param[in]  dbin       Input grid Db structure
 ** \param[in]  nmult      Refinement factor
 **
 ** \param[out] ndim       Space Dimension
 ** \param[out] ntot       Total number of cells in the output grid
 ** \param[out] nx         Array of number of cells
 ** \param[out] x0         Array of grid origin coordinates
 ** \param[out] dx         Array of grid mesh dimensions
 **
 *****************************************************************************/
GSTLEARN_EXPORT int simfine_dim(Db *dbin,
                                int nmult,
                                int *ndim,
                                int *ntot,
                                int *nx,
                                double *x0,
                                double *dx)
{
  double rmult;
  int i;

  if (!is_grid(dbin))
  {
    messerr("This simulation refinement is dedicated to Grid File");
    return (1);
  }
  rmult = pow(2., (double) nmult);
  NDIM = *ndim = dbin->getNDim();
  *ntot = 1;
  for (i = 0; i < NDIM; i++)
  {
    x0[i] = dbin->getX0(i);
    dx[i] = dbin->getDX(i);
    nx[i] = dbin->getNX(i);
    if (i < 2)
    {
      nx[i] = (int) (1 + (nx[i] - 1) * rmult);
      dx[i] = dx[i] / rmult;
    }
    (*ntot) *= nx[i];
  }
  return (0);
}

