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

/*! \cond */
#define POISSON_MEANS(ivar) (VARIO->means[(ivar)])
#define DATES(idate,i)      (vario->dates[2 * (idate) + (i)])
#define IAD(ivar,jvar)      ((ivar) * nvar + (jvar))
#define ADD(ix,iy,iz,nx)    ((iz) + nx[2] * ((iy) + nx[1] * (ix)))
#define OPP(idim,i)         (dims[idim] - i - 1)

/*! \endcond */

static Model *MODEL;
static Db *DBMAP;

static VectorDouble BETA;
static VectorDouble DRFDIAG;
static MatrixRectangular DRFXA;
static MatrixRectangular DRFGX;
static MatrixRectangular DRFTAB;
static MatrixSquareGeneral DRFXGX;
static MatrixSquareGeneral MATDRF;
static int IPTV, IPTW;

/****************************************************************************/
/*!
 **  Get variable order
 **
 ** \return  Rank of the pair of variables (-1 if incorrect arguments)
 **
 ** \param[in]  nvar      Number of variables
 ** \param[in]  ivar0     Rank of the first variable
 ** \param[in]  jvar0     Rank of the second variable
 **
 *****************************************************************************/
static int st_get_variable_order(int nvar, int ivar0, int jvar0)
{
  int rank = 0;
  for (int ivar = 0; ivar < nvar; ivar++)
    for (int jvar = 0; jvar <= ivar; jvar++, rank++)
    {
      if (ivar == ivar0 && jvar == jvar0) return rank;
      if (ivar == jvar0 && jvar == ivar0) return rank;
    }
  return -1;
}

/****************************************************************************/
/*!
 **  Manage the drift removal option
 **
 ** \param[in]  type    Type of action
 **                     0 (initialization)
 **                     +1 (allocation)
 **                     -1 (deallocation)
 ** \param[in]  db      Db structure
 ** \param[in]  model   Model structure
 **
 ** \remark If the drift removal is not called, set argument 'model' to NULL
 **
 *****************************************************************************/
void manage_drift_removal(int type, Db *db, Model *model)
{
  /* Dispatch */

  switch (type)
  {
    case 0:
      MODEL = nullptr;
      break;

    case 1:
      if (model != nullptr)
      {
        int nbfl = model->getDriftNumber();
        int nech = db->getActiveAndDefinedNumber(0);
        MODEL = model;
        BETA.resize(nbfl,0.);
        DRFDIAG.resize(nech, 0.);
        DRFTAB.reset(nech, nbfl, 0.);
        DRFXA.reset(nech, nbfl, 0.);
        DRFGX.reset(nech, nbfl, 0.);
        DRFXGX.reset(nbfl, nbfl, 0.);
        MATDRF.reset(nbfl, nbfl, 0.);
      }
      break;

    case -1:
      MODEL = nullptr;
      break;
  }
}

/****************************************************************************/
/*!
 **  Estimate the coefficients of the global drift
 **
 ** \return Error return code
 **
 ** \param[in]  db      Db descriptor
 ** \param[in]  verbose Verbose flag
 **
 *****************************************************************************/
int estimate_drift_coefficients(Db *db, int verbose)

{
  if (MODEL == nullptr) return (0);

  int iiech;
  int nbfl = MODEL->getDriftNumber();
  int nech = db->getActiveAndDefinedNumber(0);
  VectorDouble b(nbfl, 0.);

  /* Calculate: t(X) %*% X */

  for (int iech = iiech = 0; iech < db->getSampleNumber(); iech++)
  {
    if (!db->isActiveAndDefined(iech, 0)) continue;
    VectorDouble drfloc = MODEL->evalDriftVec(db, iech, ECalcMember::LHS);
    double zval = db->getLocVariable(ELoc::Z, iech, 0);

    for (int il = 0; il < nbfl; il++)
    {
      if (FFFF(drfloc[il]))
      {
        messerr("Drift cannot be calculated: term (%d) is undefined at sample (%d)",
                il + 1, iech + 1);
        return 1;
      }
      DRFTAB.setValue(iiech, il, drfloc[il]);
      b[il] += drfloc[il] * zval;
      for (int jl = 0; jl < nbfl; jl++)
        MATDRF.setValue(il,jl, MATDRF(il,jl) + drfloc[il] * drfloc[jl]);
    }
    iiech++;
  }

  /* Calculate: A = (t(X) %*% X)-1 */

  if (MATDRF.invert()) return 1;

  /* Calculate: A %*% t(X) %*% Y */

  MATDRF.prodMatVecInPlace(b, BETA);

  /* Optional printout */

  if (verbose)
  {
    message("Drift removal initial step\n");
    print_matrix("Drift Coefficients Matrix", 0, 1, nbfl, nbfl, NULL, MATDRF.getValues().data());
  }

  /* Pre-process the vector X %*% A */

  matrix_product_safe(nech, nbfl, nbfl, DRFTAB.getValues().data(),
                      MATDRF.getValues().data(), DRFXA.getValues().data());

  return 0;
}

/****************************************************************************/
/*!
 **  Calculate the a bias term between samples iech and jech
 **
 ** \param[in]  db    Db structure
 ** \param[in]  nbfl  Number of drift terms
 ** \param[in]  iiech Relative rank of the first sample
 ** \param[in]  jjech Relative rank of the second sample
 **
 *****************************************************************************/
double get_bias_value(Db *db, int nbfl, int iiech, int jjech)
{
  double bias0 = 0.;
  for (int il = 0; il < nbfl; il++)
    for (int jl = 0; jl < nbfl; jl++)
      bias0 += DRFXA.getValue(iiech, il) * DRFXGX.getValue(il, jl) * DRFXA.getValue(jjech, jl);

  double bias1 = 0.;
  for (int il = 0; il < nbfl; il++)
    bias1 += DRFXA(iiech, il) * DRFGX.getValue(jjech, il);

  double bias2 = 0.;
  for (int il = 0; il < nbfl; il++)
    bias2 += DRFGX.getValue(iiech, il) * DRFXA.getValue(jjech, il);

  return (bias0 - (bias1 + bias2));
}

/****************************************************************************/
/*!
 **  Calculate the global bias terms
 **
 ** \param[in]  db        Db description
 ** \param[in]  d1        Working vector (Dimension: ndim)
 **
 *****************************************************************************/
void calculateBiasGlobal(Db *db, VectorDouble d1)
{
  double c00, covtab, value;

  /* Initializations */

  int nbfl = MODEL->getDriftNumber();
  int ndim = MODEL->getDimensionNumber();
  int nech = db->getActiveAndDefinedNumber(0);

  /* Calculate the c00 term */

  for (int idim = 0; idim < ndim; idim++) d1[idim] = 0.;
  model_calcul_cov(NULL,MODEL, nullptr, 1, 1., d1, &c00);

  /* Calculate the term: G %*% X */

  int iiech = 0;
  for (int iech = 0; iech < db->getSampleNumber(); iech++)
  {
    if (!db->isActiveAndDefined(iech, 0)) continue;
    for (int il = 0; il < nbfl; il++)
    {
      value = 0;
      int jjech = 0;
      for (int jech = 0; jech < db->getSampleNumber(); jech++)
      {
        if (!db->isActiveAndDefined(jech, 0)) continue;
        for (int idim = 0; idim < ndim; idim++)
          d1[idim] = db->getDistance1D(iech, jech, idim);
        model_calcul_cov(NULL,MODEL, nullptr, 1, 1., d1, &covtab);
        value += (c00 - covtab) * DRFTAB.getValue(jjech, il);
        jjech++;
      }
      DRFGX.setValue(iiech, il, value);
    }
    iiech++;
  }

  /* Calculate the term: t(X) %*% G %*% X */

  for (int il = 0; il < nbfl; il++)
    for (int jl = 0; jl < nbfl; jl++)
    {
      value = 0;
      for (int iech = 0; iech < nech; iech++)
        value += DRFGX.getValue(iech, il) * DRFTAB.getValue(iech, jl);
      DRFXGX.setValue(il,jl,value);
    }

  /* Calculate the term: diag(bias) */

  iiech = 0;
  for (int iech = 0; iech < db->getSampleNumber(); iech++)
  {
    if (!db->isActiveAndDefined(iech, 0)) continue;
    DRFDIAG[iiech] = get_bias_value(db, nbfl, iiech, iiech);
    iiech++;
  }
}

double get_DRFDIAG(int iech)
{
  return DRFDIAG[iech];
}
VectorDouble get_BETA()
{
  return BETA;
}
MatrixSquareGeneral get_DRFXGX()
{
  return DRFXGX;
}

/****************************************************************************/
/*!
 **  Calculate the data value (possibly after removing the global trend)
 **
 ** \return    The data value (or the residual)
 **
 ** \param[in]  db      Db descriptor
 ** \param[in]  iech    Rank of the sample
 ** \param[in]  ivar    Rank of the variable
 **
 ** \remark  The trend removal only applies on the first variable
 ** \remark  Therefore, if applied on any variable rank other than 0,
 ** \remark  TEST is returned
 **
 *****************************************************************************/
static double st_get_IVAR(const Db *db, int iech, int ivar)
{
  double zz = db->getLocVariable(ELoc::Z, iech, ivar);
  if (FFFF(zz)) return (TEST);
  if (MODEL == nullptr) return (zz);
  if (ivar != 0) return (TEST);
  double drfval = model_drift_evaluate(0, MODEL, db, iech, 0, BETA.data());
  if (FFFF(drfval)) return (TEST);
  return (zz - drfval);
}

/****************************************************************************/
/*!
 **  Internal function for setting a VMAP value
 **
 ** \param[in]  nvar        Number of variables
 ** \param[in]  ipas        Rank of the variogram lag
 ** \param[in]  ivar        Index of the first variable
 ** \param[in]  jvar        Index of the second variable
 ** \param[in]  ww          Weight
 ** \param[in]  value       Variogram value
 **
 *****************************************************************************/
static void st_vmap_set(const ECalcVario& /*calcul_type*/,
                        int nvar,
                        int ipas,
                        int ivar,
                        int jvar,
                        int /*orient*/,
                        double ww,
                        double /*dist*/,
                        double value)
{
  int ijvar = st_get_variable_order(nvar, ivar, jvar);

  DBMAP->updArray(ipas, IPTV + ijvar, 0, ww * value);
  DBMAP->updArray(ipas, IPTW + ijvar, 0, ww);
  return;
}

/****************************************************************************/
/*!
 **  Update the variogram values
 **
 ** \param[in]  db             Db descriptor
 ** \param[in]  calcul_type    Type of calculation (ECalcVario)
 ** \param[in]  nvar           Number of variables
 ** \param[in]  iech1          Rank of the first sample
 ** \param[in]  iech2          Rank of the second sample
 ** \param[in]  ipas           Rank of the lag
 ** \param[in]  dist           Distance value
 ** \param[in]  do_asym        When FALSE, do not perform the symmetry
 ** \param[in]  st_generic_set Generic function for setting the variogram
 **
 ** \remarks: The argument 'do_asym' allows performing the double assignment
 ** \remarks: for the asymmetric functions (such as covariance)
 ** \remarks: This is due to the fact that the double is called entirely
 ** \remarks: in the calling function
 **
 *****************************************************************************/
void variogram_evaluate(Db *db,
                        const ECalcVario &calcul_type,
                        int nvar,
                        int iech1,
                        int iech2,
                        int ipas,
                        double dist,
                        int do_asym,
                        void (*st_generic_set)(const ECalcVario &calcul_type,
                                               int nvar,
                                               int iadlag,
                                               int ivar,
                                               int jvar,
                                               int orient,
                                               double ww,
                                               double dist,
                                               double value))
{
  double z11, z12, z21, z22, scale, value;

  double w1 = db->getWeight(iech1);
  double w2 = db->getWeight(iech2);
  if (FFFF(w1) || FFFF(w2)) return;
  int orient = (dist > 0) ? 1 : -1;
  dist = ABS(dist);

  switch (calcul_type.toEnum())
  {
    case ECalcVario::E_VARIOGRAM:
    case ECalcVario::E_TRANS1:
    case ECalcVario::E_TRANS2:
    case ECalcVario::E_BINORMAL:
      scale = w1 * w2;
      for (int ivar = 0; ivar < nvar; ivar++)
        for (int jvar = 0; jvar <= ivar; jvar++)
        {
          z11 = st_get_IVAR(db, iech1, ivar);
          z12 = st_get_IVAR(db, iech2, ivar);
          z21 = st_get_IVAR(db, iech1, jvar);
          z22 = st_get_IVAR(db, iech2, jvar);
          if (!FFFF(z11) && !FFFF(z21) && !FFFF(z12) && !FFFF(z22))
          {
            value = (z12 - z11) * (z22 - z21) / 2.;
            st_generic_set(calcul_type, nvar, ipas, ivar, jvar, 0, scale, dist,
                           value);
          }
        }
      break;

    case ECalcVario::E_MADOGRAM:
      scale = w1 * w2;
      for (int ivar = 0; ivar < nvar; ivar++)
        for (int jvar = 0; jvar <= ivar; jvar++)
        {
          z11 = st_get_IVAR(db, iech1, ivar);
          z12 = st_get_IVAR(db, iech2, ivar);
          z21 = st_get_IVAR(db, iech1, jvar);
          z22 = st_get_IVAR(db, iech2, jvar);
          if (!FFFF(z11) && !FFFF(z21) && !FFFF(z12) && !FFFF(z22))
          {
            value = sqrt(ABS((z12 - z11) * (z22 - z21))) / 2.;
            st_generic_set(calcul_type, nvar, ipas, ivar, jvar, 0, scale, dist,
                           value);
          }
        }
      break;

    case ECalcVario::E_RODOGRAM:
      scale = w1 * w2;
      for (int ivar = 0; ivar < nvar; ivar++)
        for (int jvar = 0; jvar <= ivar; jvar++)
        {
          ;
          z11 = st_get_IVAR(db, iech1, ivar);
          z12 = st_get_IVAR(db, iech2, ivar);
          z21 = st_get_IVAR(db, iech1, jvar);
          z22 = st_get_IVAR(db, iech2, jvar);
          if (!FFFF(z11) && !FFFF(z21) && !FFFF(z12) && !FFFF(z22))
          {
            value = pow(ABS((z12 - z11) * (z22 - z21)), 0.25) / 2.;
            st_generic_set(calcul_type, nvar, ipas, ivar, jvar, 0, scale, dist,
                           value);
          }
        }
      break;

    case ECalcVario::E_POISSON:
      scale = (w1 * w2) / (w1 + w2);
      for (int ivar = 0; ivar < nvar; ivar++)
        for (int jvar = 0; jvar <= ivar; jvar++)
        {
          z11 = st_get_IVAR(db, iech1, ivar);
          z12 = st_get_IVAR(db, iech2, ivar);
          z21 = st_get_IVAR(db, iech1, jvar);
          z22 = st_get_IVAR(db, iech2, jvar);
          if (!FFFF(z11) && !FFFF(z21) && !FFFF(z12) && !FFFF(z22)
              && (w1 > 0. && w2 > 0.))
          {
            value = (z12 - z11) * (z22 - z21) / 2.;
            st_generic_set(calcul_type, nvar, ipas, ivar, jvar, 0, scale, dist,
                           value);
          }
        }
      break;

    case ECalcVario::E_COVARIANCE:
    case ECalcVario::E_COVARIANCE_NC:
      scale = w1 * w2;
      for (int ivar = 0; ivar < nvar; ivar++)
        for (int jvar = 0; jvar <= ivar; jvar++)
        {
          z11 = st_get_IVAR(db, iech1, ivar);
          z12 = st_get_IVAR(db, iech2, ivar);
          z21 = st_get_IVAR(db, iech1, jvar);
          z22 = st_get_IVAR(db, iech2, jvar);
          if (!FFFF(z11) && !FFFF(z22))
          {
            value = z11 * z22;
            st_generic_set(calcul_type, nvar, ipas, ivar, jvar, orient, scale,
                           dist, value);
          }
          if (!FFFF(z12) && !FFFF(z21) && do_asym)
          {
            value = z12 * z21;
            st_generic_set(calcul_type, nvar, ipas, ivar, jvar, -orient, scale,
                           dist, value);
          }
        }
      break;

    case ECalcVario::E_COVARIOGRAM:
      scale = w2;
      for (int ivar = 0; ivar < nvar; ivar++)
        for (int jvar = 0; jvar <= ivar; jvar++)
        {
          z11 = st_get_IVAR(db, iech1, ivar);
          z12 = st_get_IVAR(db, iech2, ivar);
          z21 = st_get_IVAR(db, iech1, jvar);
          z22 = st_get_IVAR(db, iech2, jvar);
          if (!FFFF(z11) && !FFFF(z22))
          {
            value = z11 * z22;
            st_generic_set(calcul_type, nvar, ipas, ivar, jvar, orient, scale,
                           dist, value);
          }
          if (!FFFF(z12) && !FFFF(z21) && do_asym)
          {
            value = z12 * z21;
            st_generic_set(calcul_type, nvar, ipas, ivar, jvar, -orient, scale,
                           dist, value);
          }
        }
      break;

    case ECalcVario::E_ORDER4:
      scale = w1 * w2;
      for (int ivar = 0; ivar < nvar; ivar++)
        for (int jvar = 0; jvar <= ivar; jvar++)
        {
          z11 = st_get_IVAR(db, iech1, ivar);
          z12 = st_get_IVAR(db, iech2, ivar);
          z21 = st_get_IVAR(db, iech1, jvar);
          z22 = st_get_IVAR(db, iech2, jvar);
          if (!FFFF(z11) && !FFFF(z21) && !FFFF(z12) && !FFFF(z22))
          {
            value = (z12 - z11) * (z22 - z21);
            value = value * value / 2.;
            st_generic_set(calcul_type, nvar, ipas, ivar, jvar, 0, scale, dist,
                           value);
          }
        }
      break;

    default:
      messageAbort("st_variogram_evaluate");
      break;
  }

  return;
}

/****************************************************************************/
/*!
 **  Scale the variogram maps
 **
 ** \param[in]  dbmap     Db of Grid type containing the Variogram Maps
 ** \param[in]  nv2       nvar ( nvar + 1) /2
 **
 *****************************************************************************/
static void st_vmap_scale(Db *dbmap, int nv2)
{
  for (int iech = 0; iech < dbmap->getSampleNumber(); iech++)
  {
    for (int ijvar = 0; ijvar < nv2; ijvar++)
    {
      double value = dbmap->getArray(iech, IPTW + ijvar);
      if (value <= 0.)
        dbmap->setArray(iech, IPTV + ijvar, TEST);
      else
        dbmap->updArray(iech, IPTV + ijvar, 3, value);
    }
  }
}

/****************************************************************************/
/*!
 **  Get the absolute index of a grid node, shifted within a neighborhood
 **
 ** \return  Returned absolute index (<0 is not relevant)
 **
 ** \param[in]  dbmap     Db of Grid type containing the Variogram Maps
 ** \param[in]  indg0     Index decomposition for the central cell
 ** \param[in]  neigh     Array describing the neighborhooed
 **                       (Dimension: ndim * nbmax)
 ** \param[in]  rank      Rank of the cell within the neghborhood
 **
 ** \param[out] indg1     Working array for grid indices
 **
 *****************************************************************************/
static int st_find_neigh_cell(DbGrid *dbmap,
                              int *indg0,
                              int *neigh,
                              int rank,
                              int *indg1)
{
  int ndim;

  // Initializations

  ndim = dbmap->getNDim();

  // Get the indices of the neighboring cell

  for (int idim = 0; idim < ndim; idim++)
    indg1[idim] = indg0[idim] + neigh[rank * ndim + idim];

  return (db_index_grid_to_sample(dbmap, indg1));
}

/****************************************************************************/
/*!
 **  Calculate the variogram map when data are isolated points
 **
 ** \return  Error return code
 **
 ** \param[in]  db           Db containing the data
 ** \param[in]  dbmap        Db of Grid type containing the Variogram Maps
 ** \param[in]  calcul_type  Type of calculation (ECalcVario)
 ** \param[in]  radius       Dilation radius (used to smooth the resulting maps)
 ** \param[in]  namconv      Naming convention
 **
 *****************************************************************************/
static int st_vmap_general(Db *db,
                           DbGrid *dbmap,
                           const ECalcVario &calcul_type,
                           int radius,
                           const NamingConvention& namconv)
{
  int error, nvar, nv2, i, idim, flag_out, nbmax;
  int *indg0, *indg1, *ind1, iech0, iech1, iech2, jech1, jech2, nech, ndim;
  double *delta, *mid, *coor, x0;
  VectorInt neigh;

  /* Preliminary checks */

  error = 0;
  indg0 = indg1 = ind1 = nullptr;
  delta = coor = mid = nullptr;
  if (db == nullptr) return (1);
  if (dbmap == nullptr) return (1);

  if (! dbmap->isGrid())
  {
    messerr("This feature requires a Grid Data Base");
    messerr("to store the Variogram Maps");
    return (1);
  }
  if (db->getNDim() != 2 && db->getNDim() != 3)
  {
    messerr("The Variogram Map can only be calculated on a grid data set");
    messerr("with dimension equal to 2 or 3");
    return (1);
  }
  if (dbmap->getNDim() > db->getNDim())
  {
    messerr("The space dimension of the VMAP (%d)", dbmap->getNDim());
    messerr(
        "must not be larger than the space dimension of the input Grid (%d)",
        db->getNDim());
    return (1);
  }

  /* Initializations */

  ndim = dbmap->getNDim();
  nvar = db->getLocNumber(ELoc::Z);
  nech = db->getSampleNumber();
  nv2 = nvar * (nvar + 1) / 2;

  /* Core allocation */

  indg0 = db_indg_alloc(dbmap);
  if (indg0 == nullptr) goto label_end;
  indg1 = db_indg_alloc(dbmap);
  if (indg1 == nullptr) goto label_end;
  ind1 = (int*) mem_alloc(sizeof(int) * nech, 0);
  if (ind1 == nullptr) goto label_end;
  delta = db_sample_alloc(db, ELoc::X);
  if (delta == nullptr) goto label_end;
  mid = db_sample_alloc(db, ELoc::X);
  if (mid == nullptr) goto label_end;
  coor = db_vector_alloc(db);
  if (coor == nullptr) goto label_end;

  /* Create the variables in the Variogram Map file */

  IPTV = dbmap->addColumnsByConstant(nv2, 0.);
  if (IPTV < 0) goto label_end;
  IPTW = dbmap->addColumnsByConstant(nv2, 0.);
  if (IPTW < 0) goto label_end;

  /* Calculate a neighborhood (if radius > 0) */

  neigh = gridcell_neigh(ndim, 1, radius, false, false);
  nbmax = (int) neigh.size() / ndim;

  /* Calculate the VMAP half-extension */

  for (idim = 0; idim < ndim; idim++)
    mid[idim] = dbmap->getNX(idim) * dbmap->getDX(idim) / 2;

  /* Sorting the samples according to their first coordinate */

  if (db_coorvec_get(db, 0, coor)) goto label_end;
  for (i = 0; i < nech; i++)
    ind1[i] = i;
  ut_sort_double(1, nech, ind1, coor);

  /* Loop on the first data */

  for (jech1 = 0; jech1 < nech; jech1++)
  {
    iech1 = ind1[jech1];
    if (!db->isActive(iech1)) continue;
    x0 = db->getCoordinate(iech1, 0);

    /* Loop on the second data */

    for (jech2 = jech1; jech2 < nech; jech2++)
    {
      iech2 = ind1[jech2];
      if (!db->isActive(iech2)) continue;
      delta[0] = db->getCoordinate(iech2, 0) - x0;
      if (delta[0] > mid[0]) break;

      for (idim = 1, flag_out = 0; idim < ndim && flag_out == 0; idim++)
      {
        delta[idim] = db->getDistance1D(iech2, iech1, idim);
        if (delta[idim] > mid[idim]) flag_out = 1;
      }
      if (flag_out) continue;

      /* Evaluate the variogram map */

      DBMAP = dbmap;

      // Apply to the target cell
      if (point_to_grid(dbmap, delta, 0, indg0)) continue;
      for (int in = 0; in < nbmax; in++)
      {
        iech0 = st_find_neigh_cell(dbmap, indg0, neigh.data(), in, indg1);
        if (iech0 < 0) continue;
        variogram_evaluate(db, calcul_type, nvar, iech1, iech2, iech0, TEST,
                              0, st_vmap_set);
      }

      // Avoid symmetry if point is compared to itself
      if (iech1 == iech2) continue;

      // Apply to the opposite target cell
      for (idim = 0; idim < ndim; idim++)
        delta[idim] = -delta[idim];
      if (point_to_grid(dbmap, delta, 0, indg0)) continue;
      for (int in = 0; in < nbmax; in++)
      {
        iech0 = st_find_neigh_cell(dbmap, indg0, neigh.data(), in, indg1);
        if (iech0 < 0) continue;
        variogram_evaluate(db, calcul_type, nvar, iech1, iech2, iech0, TEST,
                              0, st_vmap_set);
      }
    }
  }

  /* Normalization */

  st_vmap_scale(dbmap, nv2);

  /* Set the error return code */

  namconv.setNamesAndLocators(db, VectorString(), ELoc::Z, -1, dbmap, IPTW, "Nb", 1, false);
  namconv.setNamesAndLocators(db, VectorString(), ELoc::Z, -1, dbmap, IPTV, "Var");
  error = 0;

  label_end: indg0 = db_indg_free(indg0);
  indg1 = db_indg_free(indg1);
  ind1 = (int*) mem_free((char* ) ind1);
  delta = db_sample_free(delta);
  mid = db_sample_free(mid);
  coor = db_vector_free(coor);
  return (error);
}

/****************************************************************************/
/*!
 **  Calculate the variogram map when data are defined on a Grid
 **
 ** \return  Error return code
 **
 ** \param[in]  dbgrid       Db of Grid type containing the data
 ** \param[in]  dbmap        Db of Grid type containing the Variogram Maps
 ** \param[in]  calcul_type  Type of calculation (ECalcVario)
 ** \param[in]  namconv      Naming Convention
 **
 *****************************************************************************/
static int st_vmap_grid(DbGrid *dbgrid,
                        DbGrid *dbmap,
                        const ECalcVario &calcul_type,
                        const NamingConvention& namconv)
{
  int error, nvar, nv2, idim, delta;
  int *ind1, *ind2, *ind0, iech0, iech1, iech2, flag_out, ndim;

  /* Preliminary checks */

  error = 0;
  ind0 = ind1 = ind2 = nullptr;
  if (dbgrid == nullptr) return (1);
  if (dbmap == nullptr) return (1);

  if (! dbgrid->isGrid())
  {
    messerr("This Variogram Map is defined for Grid Data Base only");
    return (1);
  }
  if (! dbmap->isGrid())
  {
    messerr("This feature requires a Grid Data Base");
    messerr("to store the Variogram Maps");
    return (1);
  }
  if (dbgrid->getNDim() != 2 && dbgrid->getNDim() != 3)
  {
    messerr("The Variogram Map can only be calculated on a grid data set");
    messerr("with dimension equal to 2 or 3");
    return (1);
  }
  if (dbmap->getNDim() > dbgrid->getNDim())
  {
    messerr("The space dimension of the VMAP (%d)", dbmap->getNDim());
    messerr(
        "must not be larger than the space dimension of the input Grid (%d)",
        dbgrid->getNDim());
    return (1);
  }
  for (idim = 0; idim < dbmap->getNDim(); idim++)
  {
    if (ABS(dbmap->getDX(idim) - dbgrid->getDX(idim)) > 1.e-03)
    {
      messerr("The grid mesh in the direction %d (dx=%lf)", idim + 1,
              dbgrid->getDX(idim));
      messerr("must match the mesh of the Variogram Map grid (dx=%lf)",
              dbmap->getDX(idim));
      return (1);
    }
  }

  /* Initializations */

  ndim = dbmap->getNDim();
  nvar = dbgrid->getLocNumber(ELoc::Z);
  nv2 = nvar * (nvar + 1) / 2;

  /* Core allocation */

  ind0 = db_indg_alloc(dbmap);
  if (ind0 == nullptr) goto label_end;
  ind1 = db_indg_alloc(dbgrid);
  if (ind1 == nullptr) goto label_end;
  ind2 = db_indg_alloc(dbgrid);
  if (ind2 == nullptr) goto label_end;

  /* Create the variables in the Variogram Map file */

  IPTV = dbmap->addColumnsByConstant(nv2, 0.);
  if (IPTV < 0) goto label_end;
  IPTW = dbmap->addColumnsByConstant(nv2, 0.);
  if (IPTW < 0) goto label_end;

  /* Loop on the first data */

  for (iech1 = 0; iech1 < dbgrid->getSampleNumber(); iech1++)
  {
    if (!dbgrid->isActive(iech1)) continue;
    db_index_sample_to_grid(dbgrid, iech1, ind1);

    /* Loop on the second data */

    for (iech2 = 0; iech2 < dbgrid->getSampleNumber(); iech2++)
    {
      if (!dbgrid->isActive(iech2)) continue;
      db_index_sample_to_grid(dbgrid, iech2, ind2);

      for (idim = flag_out = 0; idim < ndim && flag_out == 0; idim++)
      {
        delta = ind1[idim] - ind2[idim];
        int moitie = (dbmap->getNX(idim) - 1) / 2;
        if (delta < -moitie || delta > moitie) flag_out = 1;
        ind0[idim] = delta + moitie;
      }
      if (flag_out) continue;

      /* Evaluate the variogram map */

      DBMAP = dbmap;
      iech0 = db_index_grid_to_sample(dbmap, ind0);
      variogram_evaluate(dbgrid, calcul_type, nvar, iech1, iech2, iech0, TEST,
                            0, st_vmap_set);
    }
  }

  /* Normalization */

  st_vmap_scale(dbmap, nv2);

  /* Set the error return code */

  namconv.setNamesAndLocators(dbgrid, VectorString(), ELoc::Z, -1, dbmap, IPTW, "Nb", 1, false);
  namconv.setNamesAndLocators(dbgrid, VectorString(), ELoc::Z, -1, dbmap, IPTV, "Var");
  error = 0;

  label_end: ind0 = db_indg_free(ind0);
  ind1 = db_indg_free(ind1);
  ind2 = db_indg_free(ind2);
  return (error);
}

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
    z1 = st_get_IVAR(db, iech, 0);
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
      z2 = st_get_IVAR(db, jech, 0);
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
  double *ids, dist, z1, z2, value;
  int *indg, *rank, nech, iech, jech, igrid, idir, ideb;
  VectorDouble coor;
  SpaceTarget T1(vario->getSpace());
  SpaceTarget T2(vario->getSpace());

  /* Initializations */

  indg = rank = nullptr;
  ids = nullptr;
  const VarioParam &varioparam = vario->getVarioParam();

  // Local variables to speed up calculations
  bool hasSel = db->hasLocVariable(ELoc::SEL);
  bool hasWeight = db->hasLocVariable(ELoc::W);
  double w1 = 1.;
  double w2 = 1.;

  /* Core allocation */

  nech = db->getSampleNumber();
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
    z1 = st_get_IVAR(db, iech, 0);
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
      z2 = st_get_IVAR(db, jech, 0);
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
  double value, w1, w2, z1, z2;
  int nech, iech, jech, ideb;
  SpaceTarget T1(varioparam->getSpace());
  SpaceTarget T2(varioparam->getSpace());

  /* Preliminary calculations */

  const DirParam &dirparam = varioparam->getDirParam(idir);
  nech = db->getSampleNumber();

  // Creating a local Vario structure (to constitute the BiTargetCheck list
  Vario* vario = Vario::create(*varioparam);
  vario->setDb(db);
  if (vario->prepare()) return;

  // Local variables to speed up calculations
  bool hasSel = db->hasLocVariable(ELoc::SEL);
  bool hasWeight = db->hasLocVariable(ELoc::W);
  bool hasDate = varioparam->isDateUsed(db);

  /* Loop on the first point */

  for (iech = 0; iech < nech - 1; iech++)
  {
    if (hasSel && !db->isActive(iech)) continue;
    if (hasWeight && FFFF(db->getWeight(iech))) continue;
    db->getSampleAsST(iech, T1);

    ideb = (hasDate) ? 0 : iech + 1;
    for (jech = ideb; jech < nech; jech++)
    {
      if (hasSel && !db->isActive(jech)) continue;
      if (hasWeight && FFFF(db->getWeight(jech))) continue;
      db->getSampleAsST(jech, T2);

      // Reject the point as soon as one BiTargetChecker is not correct
      if (! vario->keepPair(idir, T1, T2, &dist)) continue;

      if (floor(dist / dirparam.getDPas() + 0.5) >= dirparam.getLagNumber())
        continue;

      w1 = db->getWeight(iech);
      w2 = db->getWeight(jech);
      if (FFFF(w1) || FFFF(w2)) continue;
      z1 = st_get_IVAR(db, iech, 0);
      z2 = st_get_IVAR(db, jech, 0);
      if (FFFF(z1) || FFFF(z2)) continue;
      value = w1 * w2 * (z2 - z1) * (z2 - z1) / 2.;
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
  int idir, iptr;

  /* Initializations */

  if (db == nullptr) return (1);
  if (dbgrid == nullptr) return (1);
  if (varioparam == (VarioParam*) NULL) return (1);
  manage_drift_removal(0, NULL, NULL);

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
  iptr = dbgrid->addColumnsByConstant(ndir, 0.);
  if (iptr < 0) return (1);

  /* Loop on the directions to evaluate */

  for (idir = 0; idir < ndir; idir++)
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
  int idir;

  /* Initializations */

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
  for (idir = 0; idir < varioparam->getDirectionNumber(); idir++)
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
 **  Product of a vector by its conjugate
 **
 ** \param[in]  size      Dimension of the vectors
 ** \param[in]  coef      Multiplicative coefficient
 ** \param[in]  tab1      First complex array
 ** \param[in]  tab2      Second complex array
 **
 ** \param[out]  tab      Output complex array
 **
 *****************************************************************************/
static void st_product_conj(int size,
                            double coef,
                            double *tab1[2],
                            double *tab2[2],
                            double *tab[2])
{
  int i;

  for (i = 0; i < size; i++)
  {
    tab[0][i] += coef * (tab1[0][i] * tab2[0][i] + tab1[1][i] * tab2[1][i]);
    tab[1][i] += coef * (tab1[0][i] * tab2[1][i] - tab1[1][i] * tab2[0][i]);
  }
}

/****************************************************************************/
/*!
 **  Load the resulting array in the output VMAP grid
 **
 ** \param[in]  dbmap     Db containing the VMAP
 ** \param[in]  tab       Array to be stored (in FFT space)
 ** \param[in]  iptr      Pointer for storage
 **
 *****************************************************************************/
static void st_vmap_store(DbGrid *dbmap, double *tab, int iptr)
{
  int ix, iy, iz, ecr, iech, indice[3];
  VectorDouble dims(3);
  int ndim = dbmap->getNDim();

  for (int idim = 0; idim < 3; idim++)
  {
    if (idim < ndim)
      dims[idim] = dbmap->getNX(idim);
    else
      dims[idim] = 1;
  }

  /* Loop on the sample (supposedly ordered) */

  ecr = 0;
  for (ix = 0; ix < dims[0]; ix++)
    for (iy = 0; iy < dims[1]; iy++)
      for (iz = 0; iz < dims[2]; iz++, ecr++)
      {
        indice[0] = ix;
        indice[1] = iy;
        indice[2] = iz;
        iech = db_index_grid_to_sample(dbmap, indice);
        dbmap->setArray(iech, iptr, tab[ecr]);
      }

  return;
}

/****************************************************************************/
/*!
 **  Blank the FFT arrays
 **
 ** \param[in] size    Dimension of the vectors
 ** \param[in] tab     Complex array to be blanked out
 **
 *****************************************************************************/
static void st_vmap_blank(int size, double *tab[2])
{
  int i, ic;

  for (ic = 0; ic < 2; ic++)
    for (i = 0; i < size; i++)
      tab[ic][i] = 0.;
  return;
}

/****************************************************************************/
/*!
 **  Scale the FFT arrays
 **
 ** \param[in] size     Dimension of the vectors
 ** \param[in] scale    Scaling factor
 ** \param[in,out] tab1 First complex array
 ** \param[in] tab2     Second complex array
 **
 *****************************************************************************/
static void st_vmap_rescale(int size, double scale, double *tab1, double *tab2)
{
  int i;
  double value;

  for (i = 0; i < size; i++)
  {
    value = tab2[i];
    if (value > EPSILON8) tab1[i] /= (scale * value);
  }
  return;
}

/****************************************************************************/
/*!
 **  Shift the product of means for the FFT arrays
 **
 ** \param[in] size    Dimension of the vectors
 ** \param[in,out] tab Input/Output complex array
 ** \param[in] tabm1   Complex array for mean of variable 1
 ** \param[in] tabm2   Complex array for mean of variable 2
 **
 *****************************************************************************/
static void st_vmap_shift(int size, double *tab, double *tabm1, double *tabm2)
{
  int i;

  for (i = 0; i < size; i++)
    tab[i] -= tabm1[i] * tabm2[i];

  return;
}

/****************************************************************************/
/*!
 **  Load the data into the FFT arrays
 **
 ** \return  Error returned code
 **
 ** \param[in] dbgrid    Db structure containing the input grid
 ** \param[in] ndim      Space dimension
 ** \param[in] sizetot   Dimension of the vectors
 ** \param[in] dims      Array of dimensions of the extended images
 ** \param[in] dinv      Array of dimensions of the extended images (inverted)
 ** \param[in] ivar      Rank of the first variable
 ** \param[in] jvar      Rank of the second variable
 **
 ** \param[out] i1i1    Pointer to the complex array 1(Z1)
 ** \param[out] z1i1    Pointer to the complex array Z1
 ** \param[out] i2i2    Pointer to the complex array 1(Z2)
 ** \param[out] z2i2    Pointer to the complex array Z2
 ** \param[out] i1i2    Pointer to the complex array 1(Z1) * 1(Z2)
 ** \param[out] z1i2    Pointer on the complex array Z1 * 1(Z2)
 ** \param[out] z2i1    Pointer on the complex array Z2 * 1(Z1)
 ** \param[out] z1z2    Pointer on the complex array Z1 * Z2
 **
 ** \remark The arrays are evaluated only if the input pointer is defined
 **
 *****************************************************************************/
static int st_vmap_load(DbGrid *dbgrid,
                        int ndim,
                        int sizetot,
                        int *dims,
                        int *dinv,
                        int ivar,
                        int jvar,
                        double *i1i1[2],
                        double *z1i1[2],
                        double *i2i2[2],
                        double *z2i2[2],
                        double *i1i2[2],
                        double *z1i2[2],
                        double *z2i1[2],
                        double *z1z2[2])
{
  int ix, iy, iz, iech, ecr, ind1, ind2, indice[3];
  double val1, val2;

  /* Initialize the complex array */

  if (i1i1 != NULL) st_vmap_blank(sizetot, i1i1);
  if (z1i1 != NULL) st_vmap_blank(sizetot, z1i1);
  if (i2i2 != NULL) st_vmap_blank(sizetot, i2i2);
  if (z2i2 != NULL) st_vmap_blank(sizetot, z2i2);
  if (i1i2 != NULL) st_vmap_blank(sizetot, i1i2);
  if (z1i2 != NULL) st_vmap_blank(sizetot, z1i2);
  if (z2i1 != NULL) st_vmap_blank(sizetot, z2i1);
  if (z1z2 != NULL) st_vmap_blank(sizetot, z1z2);

  /* Loop on the grid cells */

  for (ix = ecr = 0; ix < dims[0]; ix++)
    for (iy = 0; iy < dims[1]; iy++)
      for (iz = 0; iz < dims[2]; iz++, ecr++)
      {
        if (ndim >= 1 && ix >= dbgrid->getNX(0)) continue;
        if (ndim >= 2 && iy >= dbgrid->getNX(1)) continue;
        if (ndim >= 3 && iz >= dbgrid->getNX(2)) continue;
        indice[0] = ix;
        indice[1] = iy;
        indice[2] = iz;
        iech = db_index_grid_to_sample(dbgrid, indice);
        if (!dbgrid->getSelection(iech)) continue;
        val1 = st_get_IVAR(dbgrid, iech, jvar);
        val2 = st_get_IVAR(dbgrid, iech, ivar);
        ind1 = (!FFFF(val1));
        ind2 = (!FFFF(val2));

        if (i1i1 != NULL) i1i1[0][ecr] = (ind1);
        if (z1i1 != NULL) z1i1[0][ecr] = (ind1) ? val1 :
                                                  0.;
        if (i2i2 != NULL) i2i2[0][ecr] = (ind2);
        if (z2i2 != NULL) z2i2[0][ecr] = (ind2) ? val2 :
                                                  0.;
        if (i1i2 != NULL) i1i2[0][ecr] = (ind1 && ind2);
        if (z1i2 != NULL) z1i2[0][ecr] = (ind1 && ind2) ? val1 :
                                                          0.;
        if (z2i1 != NULL) z2i1[0][ecr] = (ind1 && ind2) ? val2 :
                                                          0.;
        if (z1z2 != NULL) z1z2[0][ecr] = (ind1 && ind2) ? val1 * val2 :
                                                          0.;
      }

  /* Perform the Fourier transform of the complex arrays defined */

  if (i1i1 != NULL && fftn(ndim, dinv, i1i1[0], i1i1[1], 1, 1.)) return (1);
  if (z1i1 != NULL && fftn(ndim, dinv, z1i1[0], z1i1[1], 1, 1.)) return (1);
  if (i2i2 != NULL && fftn(ndim, dinv, i2i2[0], i2i2[1], 1, 1.)) return (1);
  if (z2i2 != NULL && fftn(ndim, dinv, z2i2[0], z2i2[1], 1, 1.)) return (1);
  if (i1i2 != NULL && fftn(ndim, dinv, i1i2[0], i1i2[1], 1, 1.)) return (1);
  if (z1i2 != NULL && fftn(ndim, dinv, z1i2[0], z1i2[1], 1, 1.)) return (1);
  if (z2i1 != NULL && fftn(ndim, dinv, z2i1[0], z2i1[1], 1, 1.)) return (1);
  if (z1z2 != NULL && fftn(ndim, dinv, z1z2[0], z1z2[1], 1, 1.)) return (1);
  return (0);
}

/*****************************************************************to be updated***********/
/*!
 **  Allocate an array of complex values
 **
 ** \return  Error returned code
 **
 ** \param[in]  size    Dimension of the complex array
 **
 ** \param[out] tab     Complex array to be allocated
 **
 *****************************************************************************/
static int st_complex_array_alloc(int size, double *tab[2])
{
  int ic;

  for (ic = 0; ic < 2; ic++)
    tab[ic] = nullptr;

  for (ic = 0; ic < 2; ic++)
  {
    tab[ic] = (double*) mem_alloc(sizeof(double) * size, 0);
    if (tab[ic] == nullptr) return (1);
  }
  return (0);
}

/****************************************************************************/
/*!
 **  Free an array of complex values
 **
 ** \param[out] tab     Complex array to be freed
 **
 *****************************************************************************/
static void st_complex_array_free(double *tab[2])
{
  int ic;

  if (tab == NULL) return;
  for (ic = 0; ic < 2; ic++)
  {
    if (tab[ic] != nullptr) tab[ic] = (double*) mem_free((char* ) tab[ic]);
  }
  return;
}

/****************************************************************************/
/*!
 **  Load the resulting array in the output VMAP grid
 **
 ** \param[in]  nxmap   Dimensions of the Vmap Grid
 ** \param[in]  nxgrid  Dimensions of the Db Grid
 ** \param[in]  dims    Array of dimensions of the extended images
 ** \param[in]  tabin   Array to be stored (in FFT space)
 **
 ** \param[out] tabout  Array containing the resulting VMAP
 **
 *****************************************************************************/
static void st_vmap_extract(int *nxmap,
                            int *nxgrid,
                            int *dims,
                            double *tabin,
                            double *tabout)
{
  int ix, iy, iz, jx, jy, jz, nxs2, nys2, nzs2, nxloc, nyloc, nzloc;

  /* Initializations */

  nxs2 = nxmap[0] / 2;
  nys2 = nxmap[1] / 2;
  nzs2 = nxmap[2] / 2;
  nxloc = MIN(1 + nxs2, nxgrid[0]);
  nyloc = MIN(1 + nys2, nxgrid[1]);
  nzloc = MIN(1 + nzs2, nxgrid[2]);

  /* Fill the array (IX+,IY+,IZ+) */

  for (ix = 0; ix < nxloc; ix++)
    for (iy = 0; iy < nyloc; iy++)
      for (iz = 0; iz < nzloc; iz++)
      {
        jx = nxs2 + ix;
        jy = nys2 + iy;
        jz = nzs2 + iz;
        tabout[ADD(jx, jy, jz, nxmap)] = tabin[ADD(ix, iy, iz, dims)];
      }

  /* Fill the array (IX-,IY+,IZ+) */

  for (ix = 0; ix < nxloc - 1; ix++)
    for (iy = 0; iy < nyloc; iy++)
      for (iz = 0; iz < nzloc; iz++)
      {
        jx = nxs2 - ix - 1;
        jy = nys2 + iy;
        jz = nzs2 + iz;
        tabout[ADD(jx, jy, jz, nxmap)] = tabin[ADD(OPP(0,ix), iy, iz, dims)];
      }

  /* Fill the array (IX+,IY-,IZ+) */

  for (ix = 0; ix < nxloc; ix++)
    for (iy = 0; iy < nyloc - 1; iy++)
      for (iz = 0; iz < nzloc; iz++)
      {
        jx = nxs2 + ix;
        jy = nys2 - iy - 1;
        jz = nzs2 + iz;
        tabout[ADD(jx, jy, jz, nxmap)] = tabin[ADD(ix, OPP(1,iy), iz, dims)];
      }

  /* Fill the array (IX-,IY-,IZ+) */

  for (ix = 0; ix < nxloc - 1; ix++)
    for (iy = 0; iy < nyloc - 1; iy++)
      for (iz = 0; iz < nzloc; iz++)
      {
        jx = nxs2 - ix - 1;
        jy = nys2 - iy - 1;
        jz = nzs2 + iz;
        tabout[ADD(jx, jy, jz, nxmap)] = tabin[ADD(OPP(0,ix), OPP(1,iy), iz,
                                                   dims)];
      }

  /* Fill the array (IX+,IY+,IZ-) */

  for (ix = 0; ix < nxloc; ix++)
    for (iy = 0; iy < nyloc; iy++)
      for (iz = 0; iz < nzloc - 1; iz++)
      {
        jx = nxs2 + ix;
        jy = nys2 + iy;
        jz = nzs2 - iz - 1;
        tabout[ADD(jx, jy, jz, nxmap)] = tabin[ADD(ix, iy, OPP(2,iz), dims)];
      }

  /* Fill the array (IX-,IY+,IZ-) */

  for (ix = 0; ix < nxloc - 1; ix++)
    for (iy = 0; iy < nyloc; iy++)
      for (iz = 0; iz < nzloc - 1; iz++)
      {
        jx = nxs2 - ix - 1;
        jy = nys2 + iy;
        jz = nzs2 - iz - 1;
        tabout[ADD(jx, jy, jz, nxmap)] = tabin[ADD(OPP(0,ix), iy, OPP(2,iz),
                                                   dims)];
      }

  /* Fill the array (IX+,IY-,IZ-) */

  for (ix = 0; ix < nxloc; ix++)
    for (iy = 0; iy < nyloc - 1; iy++)
      for (iz = 0; iz < nzloc - 1; iz++)
      {
        jx = nxs2 + ix;
        jy = nys2 - iy - 1;
        jz = nzs2 - iz - 1;
        tabout[ADD(jx, jy, jz, nxmap)] = tabin[ADD(ix, OPP(1,iy), OPP(2,iz),
                                                   dims)];
      }

  /* Fill the array (IX-,IY-,IZ-) */

  for (ix = 0; ix < nxloc - 1; ix++)
    for (iy = 0; iy < nyloc - 1; iy++)
      for (iz = 0; iz < nzloc - 1; iz++)
      {
        jx = nxs2 - ix - 1;
        jy = nys2 - iy - 1;
        jz = nzs2 - iz - 1;
        tabout[ADD(jx, jy, jz, nxmap)] = tabin[ADD(OPP(0,ix), OPP(1,iy),
                                                   OPP(2,iz), dims)];
      }

  return;
}

/*****************************************************************************/
/*!
 **  Calculates the 2-D variogram map on a grid using FFT
 **
 ** \return  Error return code
 **
 ** \param[in]  dbgrid       Db of Grid type containing the data
 ** \param[in]  dbmap        Db of Grid type containing the Variogram Maps
 ** \param[in]  calcul_type  Type of calculation (ECalcVario)
 ** \param[in]  namconv      Naming convention
 **
 *****************************************************************************/
static int st_vmap_grid_fft(DbGrid *dbgrid,
                            DbGrid *dbmap,
                            const ECalcVario &calcul_type,
                            const NamingConvention& namconv)
{
  int dims[3], dinv[3], nxmap[3], nxgrid[3], sizetot, sizemap, sizegrid;
  int i, ndim, ic, error, idim, ivar, jvar, ijvar, nvar, nvs2;
  double *i1i1[2], *z1i1[2], *i2i2[2], *z2i2[2], *i1i2[2];
  double *z1i2[2], *z2i1[2], *z1z2[2], *ztab[2];
  double *res_nn, *res_gg, *res_m1, *res_m2;
  static int verbose = 0;

  /* Initializations */

  error = 1;
  sizetot = 0;
  res_nn = res_gg = res_m1 = res_m2 = nullptr;
  for (ic = 0; ic < 2; ic++)
  {
    i1i1[ic] = z1i1[ic] = i2i2[ic] = z2i2[ic] = nullptr;
    i1i2[ic] = z1i2[ic] = z2i1[ic] = z1z2[ic] = nullptr;
    ztab[ic] = nullptr;
  }

  /* Preliminary checks */

  if (dbgrid == nullptr) return (1);
  if (dbmap == nullptr) return (1);

  if (calcul_type != ECalcVario::VARIOGRAM && calcul_type
      != ECalcVario::COVARIOGRAM
      && calcul_type != ECalcVario::COVARIANCE
      && calcul_type != ECalcVario::COVARIANCE_NC)
  {
    messerr("This function is limited to the calculation of");
    messerr("Variogram, Covariance (centered or not) or Covariogram");
    return (1);
  }
  if (! dbgrid->isGrid())
  {
    messerr("This Variogram Map is defined for Grid Data Base only");
    return (1);
  }
  if (! dbmap->isGrid())
  {
    messerr("This feature requires a Grid Data Base");
    messerr("to store the Variogram Maps");
    return (1);
  }
  if (dbgrid->getNDim() != 2 && dbgrid->getNDim() != 3)
  {
    messerr("The Variogram Map can only be calculated on a grid data set");
    messerr("with dimension equal to 2 or 3");
    return (1);
  }
  if (dbmap->getNDim() > dbgrid->getNDim())
  {
    messerr("The space dimension of the VMAP (%d)", dbmap->getNDim());
    messerr(
        "must not be larger than the space dimension of the input Grid (%d)",
        dbgrid->getNDim());
    return (1);
  }
  for (idim = 0; idim < dbmap->getNDim(); idim++)
  {
    if (ABS(dbmap->getDX(idim) - dbgrid->getDX(idim)) > 1.e-03)
    {
      messerr("The grid mesh in the direction %d (dx=%lf)", idim + 1,
              dbgrid->getDX(idim));
      messerr("must match the mesh of the Variogram Map grid (dx=%lf)",
              dbmap->getDX(idim));
      return (1);
    }
  }

  for (idim = 0; idim < 3; idim++)
    nxgrid[idim] = nxmap[idim] = 1;
  for (idim = 0; idim < dbgrid->getNDim(); idim++)
    nxgrid[idim] = dbgrid->getNX(idim);
  for (idim = 0; idim < dbmap->getNDim(); idim++)
    nxmap[idim] = dbmap->getNX(idim);

  /* Preliminary calculations */

  nvar = dbgrid->getLocNumber(ELoc::Z);
  nvs2 = nvar * (nvar + 1) / 2;
  ndim = 0;
  sizetot = sizemap = sizegrid = 1;
  for (i = 0; i < 3; i++)
  {
    dinv[i] = 1;
    if (nxgrid[i] <= 1)
    {
      dims[i] = 1;
    }
    else
    {
      dims[i] = (int) ceil((double) (nxgrid[i] + nxmap[i] - 1) / 8.) * 8;
      sizegrid *= nxgrid[i];
      sizemap *= nxmap[i];
      sizetot *= dims[i];
      ndim++;
    }
  }
  for (i = 0; i < ndim; i++)
    dinv[i] = dims[ndim - i - 1];
  if (verbose)
  {
    mestitle(0, "Simulation of a grid using FFT");
    message("Grid dimension:");
    for (idim = 0; idim < ndim; idim++)
      message(" %4d", nxgrid[idim]);
    message("\n");
    message("Variogram Map :");
    for (idim = 0; idim < ndim; idim++)
      message(" %4d", nxmap[idim]);
    message("\n");
    message("Working array :");
    for (idim = 0; idim < ndim; idim++)
      message(" %4d", dinv[idim]);
    message("\n");
  }

  /* Create the variables in the Variogram Map file */

  IPTV = dbmap->addColumnsByConstant(nvs2, TEST);
  if (IPTV < 0) goto label_end;
  IPTW = dbmap->addColumnsByConstant(nvs2, TEST);
  if (IPTW < 0) goto label_end;

  /* Core allocation */

  if (st_complex_array_alloc(sizetot, ztab)) goto label_end;
  if (calcul_type == ECalcVario::VARIOGRAM)
  {
    if (st_complex_array_alloc(sizetot, i1i2)) goto label_end;
    if (st_complex_array_alloc(sizetot, z1i2)) goto label_end;
    if (st_complex_array_alloc(sizetot, z2i1)) goto label_end;
    if (st_complex_array_alloc(sizetot, z1z2)) goto label_end;
  }
  else
  {
    if (st_complex_array_alloc(sizetot, i1i1)) goto label_end;
    if (st_complex_array_alloc(sizetot, z1i1)) goto label_end;
    if (st_complex_array_alloc(sizetot, i2i2)) goto label_end;
    if (st_complex_array_alloc(sizetot, z2i2)) goto label_end;
  }

  res_nn = (double*) mem_alloc(sizeof(double) * sizemap, 0);
  if (res_nn == nullptr) goto label_end;
  res_gg = (double*) mem_alloc(sizeof(double) * sizemap, 0);
  if (res_gg == nullptr) goto label_end;
  for (i = 0; i < sizemap; i++)
    res_nn[i] = res_gg[i] = 0.;
  if (calcul_type == ECalcVario::COVARIANCE || calcul_type
      == ECalcVario::COVARIOGRAM)
  {
    res_m1 = (double*) mem_alloc(sizeof(double) * sizemap, 0);
    if (res_m1 == nullptr) goto label_end;
    res_m2 = (double*) mem_alloc(sizeof(double) * sizemap, 0);
    if (res_m2 == nullptr) goto label_end;
    for (i = 0; i < sizemap; i++)
      res_m1[i] = res_m2[i] = 0.;
  }

  /* Loop on the variables */

  for (ivar = ijvar = 0; ivar < nvar; ivar++)
    for (jvar = 0; jvar <= ivar; jvar++, ijvar++)
    {

      /* Calculate the structural function */

      if (calcul_type == ECalcVario::VARIOGRAM)
      {
        if (st_vmap_load(dbgrid, ndim, sizetot, dims, dinv, ivar, jvar,
        NULL, NULL, NULL, NULL, i1i2, z1i2, z2i1, z1z2)) continue;

        /* Calculate the number of pairs */
        st_vmap_blank(sizetot, ztab);
        st_product_conj(sizetot, 1., i1i2, i1i2, ztab);
        if (fftn(ndim, dinv, ztab[0], ztab[1], -1, -1.)) continue;
        st_vmap_extract(nxmap, nxgrid, dims, ztab[0], res_nn);

        /* Structural function */
        st_vmap_blank(sizetot, ztab);
        st_product_conj(sizetot, 1., z1z2, i1i2, ztab);
        st_product_conj(sizetot, 1., i1i2, z1z2, ztab);
        st_product_conj(sizetot, -1., z1i2, z2i1, ztab);
        st_product_conj(sizetot, -1., z2i1, z1i2, ztab);
        if (fftn(ndim, dinv, ztab[0], ztab[1], -1, -1.)) continue;
        st_vmap_extract(nxmap, nxgrid, dims, ztab[0], res_gg);
        st_vmap_rescale(sizemap, 2., res_gg, res_nn);
      }
      else
      {
        if (st_vmap_load(dbgrid, ndim, sizetot, dims, dinv, ivar, jvar, i1i1,
                         z1i1, i2i2, z2i2, NULL, NULL, NULL, NULL)) continue;

        /* Calculate the number of pairs */
        st_vmap_blank(sizetot, ztab);
        st_product_conj(sizetot, 1., i1i1, i2i2, ztab);
        if (fftn(ndim, dinv, ztab[0], ztab[1], -1, -1.)) continue;
        st_vmap_extract(nxmap, nxgrid, dims, ztab[0], res_nn);

        if (calcul_type == ECalcVario::COVARIANCE || calcul_type
            == ECalcVario::COVARIOGRAM)
        {
          /* Calculate the means */
          st_vmap_blank(sizetot, ztab);
          st_product_conj(sizetot, 1., z1i1, i2i2, ztab);
          if (fftn(ndim, dinv, ztab[0], ztab[1], -1, -1.)) continue;
          st_vmap_extract(nxmap, nxgrid, dims, ztab[0], res_m1);
          st_vmap_rescale(sizemap, 1., res_m1, res_nn);

          st_vmap_blank(sizetot, ztab);
          st_product_conj(sizetot, 1., i1i1, z2i2, ztab);
          if (fftn(ndim, dinv, ztab[0], ztab[1], -1, -1.)) continue;
          st_vmap_extract(nxmap, nxgrid, dims, ztab[0], res_m2);
          st_vmap_rescale(sizemap, 1., res_m2, res_nn);
        }

        /* Structural function */
        st_vmap_blank(sizetot, ztab);
        st_product_conj(sizetot, 1., z1i1, z2i2, ztab);
        if (fftn(ndim, dinv, ztab[0], ztab[1], -1, -1.)) continue;
        st_vmap_extract(nxmap, nxgrid, dims, ztab[0], res_gg);
        st_vmap_rescale(sizemap, 1., res_gg, res_nn);

        if (calcul_type == ECalcVario::COVARIANCE || calcul_type
            == ECalcVario::COVARIOGRAM)
          st_vmap_shift(sizemap, res_gg, res_m1, res_m2);
      }

      /* Store the results */
      st_vmap_store(dbmap, res_nn, IPTW + ijvar);
      st_vmap_store(dbmap, res_gg, IPTV + ijvar);
    }

  /* Set the error return code */

  namconv.setNamesAndLocators(dbgrid, VectorString(), ELoc::Z, -1, dbmap, IPTW, "Nb", 1, false);
  namconv.setNamesAndLocators(dbgrid, VectorString(), ELoc::Z, -1, dbmap, IPTV, "Var", 1);
  error = 0;

  label_end: st_complex_array_free(ztab);
  if (calcul_type == ECalcVario::VARIOGRAM)
  {
    st_complex_array_free(i1i2);
    st_complex_array_free(z1i2);
    st_complex_array_free(z2i1);
    st_complex_array_free(z1z2);
  }
  else
  {
    st_complex_array_free(i1i1);
    st_complex_array_free(z1i1);
    st_complex_array_free(i2i2);
    st_complex_array_free(z2i2);
  }
  res_nn = (double*) mem_free((char* ) res_nn);
  res_gg = (double*) mem_free((char* ) res_gg);
  res_m1 = (double*) mem_free((char* ) res_m1);
  res_m1 = (double*) mem_free((char* ) res_m2);
  return (error);
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

/****************************************************************************/
/*!
 **  Calculate the variogram map
 **
 ** \return  Error return code
 **
 ** \param[in]  db          Db containing the data
 ** \param[in]  dbmap       VMAP grid structure
 ** \param[in]  calcul_type Type of calculation (ECalcVario)
 ** \param[in]  radius      Dilation radius (mooth resulting maps) only on points
 ** \param[in]  flag_FFT    Use FFT method (only valid on grid)
 ** \param[in]  namconv     Naming convention
 **
 *****************************************************************************/
int vmap_compute(Db *db,
                 DbGrid *dbmap,
                 const ECalcVario &calcul_type,
                 int radius,
                 bool flag_FFT,
                 const NamingConvention& namconv)
{
  int error = 0;

  // Calculating the variogram map in different ways

  if (db->isGrid())
  {
    DbGrid* dbgrid = dynamic_cast<DbGrid*>(db);

    // Case where Data are on a regular grid

    if (flag_FFT)
      error = st_vmap_grid_fft(dbgrid, dbmap, calcul_type, namconv);
    else
      error = st_vmap_grid(dbgrid, dbmap, calcul_type, namconv);
  }
  else
  {

    // Case where Data are on a set of points

    error = st_vmap_general(db, dbmap, calcul_type, radius, namconv);
  }

  return error;
}

/****************************************************************************/
/*!
 **  Calculate the variogram map
 **
 ** \return  Error return code
 **
 ** \param[in]  db          Db containing the data
 ** \param[in]  calcul_type Type of calculation (ECalcVario)
 ** \param[in]  nx_arg      Vector of (Half-) number of nodes for Vmap (def:20)
 ** \param[in]  dxx         Vector of mesh for Vmap (see details)
 ** \param[in]  radius      Dilation radius (mooth resulting maps) only on points
 ** \param[in]  flag_FFT    Use FFT method (only valid on grid)
 ** \param[in]  namconv     Naming convention
 **
 ** \remarks For calculating the default values:
 ** \remarks - for nx: it is set to 20 in all directions
 ** \remarks - for dx:
 ** \remarks   . If 'Db' is a grid, the mesh of the grid is used
 ** \remarks   - Otherwise, the mesh is set to the field extension / nx
 **
 *****************************************************************************/
DbGrid* db_vmap_compute(Db *db,
                        const ECalcVario &calcul_type,
                        const VectorInt& nx_arg,
                        const VectorDouble& dxx,
                        int radius,
                        bool flag_FFT,
                        const NamingConvention& namconv)
{
  int error = 0;

  // Creating the output Variogram Map grid

  int ndim = db->getNDim();
  VectorInt nxx = nx_arg;
  if (nxx.empty()) nxx.resize(ndim, 20);
  if (ndim != (int) nxx.size())
  {
    messerr("Argument 'nxx' should have same Space Dimension as 'db'");
    return nullptr;
  }
  if (! dxx.empty() && ndim != (int) dxx.size())
  {
    messerr("Argument 'dxx'  should have same Space Dimension as 'db'");
    return nullptr;
  }
  VectorInt nx_map(ndim);
  VectorDouble dx_map(ndim);
  VectorDouble x0_map(ndim);

  for (int idim = 0; idim<ndim; idim++)
    nx_map[idim] = 2 * nxx[idim] + 1;
  if (db->isGrid())
  {
    DbGrid* dbgrid = dynamic_cast<DbGrid*>(db);
    for (int idim = 0; idim < ndim; idim++)
      dx_map[idim] = dbgrid->getDX(idim);
  }
  else
  {
    for (int idim = 0; idim < ndim; idim++)
      dx_map[idim] = (! dxx.empty() && !FFFF(dxx[idim])) ?
          dxx[idim] : db->getExtension(idim) / (double) nxx[idim];
  }
  for (int idim = 0; idim < ndim; idim++)
    x0_map[idim] = -nxx[idim] * dx_map[idim];

  DbGrid *dbmap = DbGrid::create(nx_map, dx_map, x0_map);

  // Calculating the variogram map in different ways

  error = vmap_compute(db, dbmap, calcul_type, radius, flag_FFT, namconv);

  // In case of error, free the newly created VMAP structure

  if (error)
  {
    delete dbmap;
    dbmap = nullptr;
  }
  return dbmap;
}

