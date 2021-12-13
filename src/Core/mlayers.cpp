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
#include "geoslib_old_f.h"
#include "Variogram/Vario.hpp"
#include "Basic/Utilities.hpp"
#include "Model/Model.hpp"
#include "Neigh/Neigh.hpp"
#include "Db/Db.hpp"

#include <math.h>

/*! \cond */
#define IAD(n,i,j)        ((n) * (i) + (j))
#define A(i,j)            (a[IAD(neq,i,j)])
#define ACOV(i,j)         (acov[IAD(nech,i,j)])
#define GS(i,j)           (gs[IAD(npar,i,j)])
#define C(i,j)            (covtab[IAD(nlayers,i,j)])
#define PHIA(i,ilayer)    (phia[IAD(nlayers,i,ilayer)])
#define PHIB(i,ilayer)    (phib[IAD(nlayers,i,ilayer)])
#define AA(i,ilayer2)     (aa[IAD(nlayer2,i,ilayer2)])
#define A2(n,i,j)         (a2[IAD(n,i,j)])
#define B2(n,i,j)         (b2[IAD(n,i,j)])
#define INVS(npar,i,j)    (invS[IAD(npar,i,j)])
#define FFTAB(ipar,iech)  (fftab[(iech) * npar + (ipar)])
#define POST_S(npar,i,j)  (post_S[IAD(npar,i,j)])
#define ATAB(n,i,j)       (atab[IAD(n,i,j)])
#define VARS(n,i,j)       (vars[IAD(n,i,j)])
/*! \endcond */

typedef struct
{
  int flag_same; /* 1 if input and output files coincide */
  int flag_vel; /* 1 for velocity; 0 for thickness */
  int flag_cumul; /* 1 for cumulating in depth */
  int flag_ext; /* Use external drift */
  int flag_z; /* 1 if output must be converted into depth */
  int colrefd; /* Reference depth map (if >= 0) */
  int colreft; /* Reference time map (if >= 0) */
  int colrefb; /* Bottom map (if >= 0) */
  int match_time; /* 1 if Time provided through External Drift */
  ELoc ptime; /* Pointer to the Time variables */
  int nlayers; /* Number of layers */
  int nbfl; /* Number of drift functions */
  int nech; /* Number of active samples */
  int neq; /* Number of equations */
  int npar; /* Nb. of layers * Nb. of drift functions */
} LMlayers;

/****************************************************************************/
/*!
 **  Free the Multi-Layers internal structure
 **
 ** \return Pointer to the freed structure
 **
 ** \param[in]  lmlayers  Pointer to the LMlayers structure to be freed
 **
 *****************************************************************************/
static LMlayers* lmlayers_free(LMlayers *lmlayers)
{
  if (lmlayers == nullptr) return (lmlayers);
  lmlayers = (LMlayers*) mem_free((char* ) lmlayers);
  return (lmlayers);
}

/****************************************************************************/
/*!
 **  Returns the number of drift functions
 **
 ** \return  Number of drift conditions
 **
 ** \param[in]  irf_rank  Rank of the Intrinsic Random Function (0 or 1)
 ** \param[in]  flag_ext  1 if external drift must be used; 0 otherwise
 **
 *****************************************************************************/
static int st_get_number_drift(int irf_rank, int flag_ext)
{
  switch (irf_rank)
  {
    case -1:
      return (0);
      break;

    case 0:
      if (!flag_ext)
        return (1);
      else
        return (2);
      break;

    case 1:
      if (!flag_ext)
        return (3);
      else
        return (4);
      break;

    default:
      messageAbort("Irf_rank must be -1, 0 or 1");
      break;
  }
  return (0);
}

/****************************************************************************/
/*!
 **  Allocate the Multi-Layers internal structure
 **
 ** \return Pointer to the allocated structure
 **
 ** \param[in]  flag_same  1 if input and output files coincide
 ** \param[in]  flag_vel   1 if work is performed in Velocity, 0 for Depth
 ** \param[in]  flag_cumul 1 if work must be done on Depth; 0 on Thickness
 ** \param[in]  flag_ext   1 if external drift must be used; 0 otherwise
 ** \param[in]  flag_z     1 if the output must be provided in depth
 ** \param[in]  colrefd    Rank of the reference depth variable
 ** \param[in]  colreft    Rank of the reference Time variable in Dbout
 ** \param[in]  colrefb    Rank of the reference Bottom Depth variable
 ** \param[in]  irf_rank   Rank of the Intrinsic Random Function (0 or 1)
 ** \param[in]  match_time Pointer to the Time pointer
 **                        (1 if defined as ELoc::F or 0 for ELoc::TIME)
 ** \param[in]  nlayers    Number of layers
 **
 *****************************************************************************/
static LMlayers* lmlayers_alloc(int flag_same,
                                int flag_vel,
                                int flag_cumul,
                                int flag_ext,
                                int flag_z,
                                int colrefd,
                                int colreft,
                                int colrefb,
                                int irf_rank,
                                int match_time,
                                int nlayers)
{
  LMlayers *lmlayers;

  lmlayers = (LMlayers*) mem_alloc(sizeof(LMlayers), 1);
  lmlayers->flag_same = flag_same;
  lmlayers->flag_vel = flag_vel;
  lmlayers->flag_cumul = flag_cumul;
  lmlayers->flag_ext = flag_ext;
  lmlayers->flag_z = flag_z;
  lmlayers->colrefd = colrefd;
  lmlayers->colreft = colreft;
  lmlayers->colrefb = colrefb;
  lmlayers->match_time = match_time;
  lmlayers->ptime = (match_time) ? ELoc::F :
                                   ELoc::TIME;
  lmlayers->nlayers = nlayers;
  lmlayers->nbfl = st_get_number_drift(irf_rank, flag_ext);
  lmlayers->nech = 0;
  lmlayers->neq = 0;
  lmlayers->npar = lmlayers->nbfl * nlayers;
  return (lmlayers);
}

/****************************************************************************/
/*!
 **  Print the Multi-Layers internal structure
 **
 ** \param[in]  lmlayers  Pointer to the LMlayers structure
 **
 *****************************************************************************/
static void lmlayers_print(LMlayers *lmlayers)

{
  static const char *NOK[] = { "NO", "YES" };

  if (lmlayers == nullptr) return;

  mestitle(0, "Multi-Layers Environments");
  if (lmlayers->flag_vel)
    message("Working in Velocity\n");
  else
    message("Working in Depth\n");
  if (lmlayers->flag_cumul)
    message("Producing estimation in Depth\n");
  else
    message("Producing estimation in Thickness\n");
  if (lmlayers->flag_z) message("Results are converted into Depth\n");

  message("Do the Input and Output Db coincide: %s\n",
          NOK[lmlayers->flag_same]);
  message("Using External Drift functions: %s\n", NOK[lmlayers->flag_ext]);
  message("Is Time used as External Drift: %s\n", NOK[lmlayers->match_time]);
  if (lmlayers->colrefd >= 0)
    message("Rank of the Reference Depth Map = %d\n", lmlayers->colrefd);
  if (lmlayers->colreft >= 0)
    message("Rank of the Reference Time Map = %d\n", lmlayers->colreft);
  if (lmlayers->colrefb >= 0)
    message("Rank of the Bottom Depth Map = %d\n", lmlayers->colrefb);

  message("\n");
  message("Number of layers = %d\n", lmlayers->nlayers);
  message("Number of drift functions (per layer) = %d\n", lmlayers->nbfl);
  message("Number of active samples (including collocated duplicates) = %d\n",
          lmlayers->nech);
  message("\n");
  return;
}

/****************************************************************************/
/*!
 **  Returns the absolute grid node absolute index which is the closest to a
 **  given sample of a Db
 **  In the case of same input and output file, simply return 'iech'
 **
 ** \return 1 if the sample does not belong to the grid; 0 otherwise
 **
 ** \param[in]  lmlayers  LMlayers structure
 ** \param[in]  dbin      Input Db structure
 ** \param[in]  dbout     Output Db structure
 ** \param[in]  iech      Rank in the input Db
 ** \param[out] igrid     Rank of the node in the output Db
 **
 *****************************************************************************/
static int st_locate_sample_in_output(LMlayers *lmlayers,
                                      Db *dbin,
                                      Db *dbout,
                                      int iech,
                                      int *igrid)
{
  int idim, indg[2];
  double coor[2];

  /* In the case the input and output files coincide, simply return 'iech' */
  if (lmlayers->flag_same)
  {
    *igrid = iech;
    return (0);
  }

  /* The files are different */
  for (idim = 0; idim < dbin->getNDim(); idim++)
    coor[idim] = dbin->getCoordinate(iech, idim);
  if (point_to_grid(dbout, coor, 0, indg) != 0) return (1);
  *igrid = db_index_grid_to_sample(dbout, indg);
  return (0);
}

/****************************************************************************/
/*!
 **  Check if the target layer rank is consistent
 **
 ** \param[in]  string    Name of the calling procedure
 ** \param[in]  lmlayers  LMlayers structure
 ** \param[in]  ilayer0   Rank of the target layer (starting from 1)
 **
 ** \remarks If this target layer rank is not correct, mes_abort() is called
 ** \remarks and the program is interrupted as this must never happen.
 **
 *****************************************************************************/
static void st_check_layer(const char *string, LMlayers *lmlayers, int ilayer0)
{
  if (ilayer0 >= 1 && ilayer0 <= lmlayers->nlayers) return;

  messerr("Error when calling function %s", string);
  messerr("- Number of layers         = %d", lmlayers->nlayers);
  messerr("- Rank of the target layer = %d", ilayer0);
  messageAbort("This error should never happen");
}

/****************************************************************************/
/*!
 **  Fill the proportion vector at a output location, up to the target layer
 **
 ** \returns 1 if the proportion vector cannot be defined; 0 otherwise
 **
 ** \param[in]  lmlayers  LMlayers structure
 ** \param[in]  dbout     Output Db structure
 ** \param[in]  iech      Rank of the target sample (in output Db)
 ** \param[in]  ilayer0   Rank of the target layer (starting from 1)
 **
 ** \param[out] props     Working array (Dimension: nlayers)
 **
 *****************************************************************************/
static int st_get_props_result(LMlayers *lmlayers,
                               Db *dbout,
                               int iech,
                               int ilayer0,
                               double *props)
{
  double pval, t0, t1, tlast, tt;
  int ilayer;

  /* Initializations */

  st_check_layer("st_get_props_result", lmlayers, ilayer0);
  for (ilayer = 0; ilayer < lmlayers->nlayers; ilayer++)
    props[ilayer] = 0.;

  /* Dispatch */

  if (lmlayers->flag_vel)
  {

    /* Working in velocities */

    t0 = (lmlayers->colreft >= 0) ? dbout->getArray(iech, lmlayers->colreft) :
                                    0.;
    if (FFFF(t0)) return (1);
    t1 = get_LOCATOR_ITEM(dbout, lmlayers->ptime, ilayer0 - 1, iech);
    if (FFFF(t1)) return (1);
    tlast = t0;

    /* Loop on the layers */

    for (ilayer = 0; ilayer < ilayer0; ilayer++)
    {
      tt = get_LOCATOR_ITEM(dbout, lmlayers->ptime, ilayer, iech);
      if (FFFF(tt)) return (1);
      pval = (tt - tlast) / (t1 - t0);
      if (pval < 0 || pval > 1) return (1);
      tlast = tt;
      props[ilayer] = pval;
    }
  }
  else
  {

    /* Working in depth */

    for (ilayer = 0; ilayer < ilayer0; ilayer++)
      props[ilayer] = 1.;
  }
  return (0);
}

/****************************************************************************/
/*!
 **  Fill the proportion vector at a data location, up to the target layer
 **
 ** \returns 1 if the proportion vector cannot be defined; 0 otherwise
 **
 ** \param[in]  lmlayers  LMlayers structure
 ** \param[in]  dbin      Input Db structure
 ** \param[in]  dbout     Output Db structure
 ** \param[in]  iech      Rank of the target sample
 ** \param[in]  ilayer0   Rank of the target layer (starting from 1)
 **
 ** \param[out] props     Working array (Dimension: nlayers)
 **
 *****************************************************************************/
static int st_get_props_data(LMlayers *lmlayers,
                             Db *dbin,
                             Db *dbout,
                             int iech,
                             int ilayer0,
                             double *props)
{
  int igrid, ilayer;

  /* Initializations */

  st_check_layer("st_get_props_data", lmlayers, ilayer0);
  for (ilayer = 0; ilayer < lmlayers->nlayers; ilayer++)
    props[ilayer] = 0.;

  /* Get the sample rank in the output Db of the sample from the input Db */

  if (st_locate_sample_in_output(lmlayers, dbin, dbout, iech, &igrid))
    return (1);

  /* Evaluate the proportion vector */

  if (st_get_props_result(lmlayers, dbout, igrid, ilayer0, props)) return (1);

  return (0);
}

/****************************************************************************/
/*!
 **  Return the external drift value at an output location for a target layer
 **
 ** \returns  The external drift value of TEST
 **
 ** \param[in]  lmlayers  LMlayers structure
 ** \param[in]  dbout     Output Db structure
 ** \param[in]  iech      Rank of the target sample (in output Db)
 ** \param[in]  ilayer0   Rank of the target layer (Starting from 1)
 **
 *****************************************************************************/
static double st_get_drift_result(LMlayers *lmlayers,
                                  Db *dbout,
                                  int iech,
                                  int ilayer0)
{
  double drift;

  if (!lmlayers->flag_ext) return (TEST);
  st_check_layer("st_get_drift_result", lmlayers, ilayer0);

  drift = dbout->getExternalDrift(iech, ilayer0 - 1);
  return (drift);
}

/****************************************************************************/
/*!
 **  Return the external drift value at an input location for a target layer
 **
 ** \returns  The external drift value of TEST
 **
 ** \param[in]  lmlayers  LMlayers structure
 ** \param[in]  dbin      Input Db structure
 ** \param[in]  dbout     Output Db structure
 ** \param[in]  iech      Rank of the target sample (in output Db)
 ** \param[in]  ilayer0   Rank of the target layer (Starting from 1)
 **
 *****************************************************************************/
static double st_get_drift_data(LMlayers *lmlayers,
                                Db *dbin,
                                Db *dbout,
                                int iech,
                                int ilayer0)
{
  int igrid;
  double drift;

  if (!lmlayers->flag_ext) return (TEST);
  st_check_layer("st_get_drift_data", lmlayers, ilayer0);

  /* Get the sample rank in the output Db of the sample from the input Db*/

  if (st_locate_sample_in_output(lmlayers, dbin, dbout, iech, &igrid))
    return (TEST);

  drift = st_get_drift_result(lmlayers, dbout, igrid, ilayer0);
  return (drift);
}

/****************************************************************************/
/*!
 **  Calculate the array of covariances for zero distance
 **
 ** \param[in]  lmlayers  Pointer to the LMlayers structure to be freed
 ** \param[in]  model     Model
 ** \param[in]  prop1     Working array at first point (Dimension: nlayers)
 **
 ** \param[out] covtab    Working array (Dimension = nlayers * nlayers)
 ** \param[out] c00       Returned array (Dimension = nlayers)
 **
 ** \remarks:  This array depends on the target location through proportions
 **
 *****************************************************************************/
static void st_covariance_c00(LMlayers *lmlayers,
                              Model *model,
                              double *prop1,
                              double *covtab,
                              double *c00)
{
  int nlayers, flag_interrupt;
  double value;
  CovCalcMode mode;

  nlayers = lmlayers->nlayers;
  model_calcul_cov(model, mode, 1, 1., VectorDouble(), covtab);

  if (lmlayers->flag_cumul)
  {
    for (int k = 0; k < nlayers; k++)
    {
      value = 0.;
      flag_interrupt = 0;
      for (int i = 0; i <= k && flag_interrupt == 0; i++)
        for (int j = 0; j <= k && flag_interrupt == 0; j++)
        {
          if (FFFF(prop1[i]) || FFFF(prop1[j]))
            flag_interrupt = 1;
          else
            value += prop1[i] * prop1[j] * C(i, j);
        }
      c00[k] = (flag_interrupt) ? TEST :
                                  value;
    }
  }
  else
  {
    for (int k = 0; k < nlayers; k++)
      c00[k] = C(k, k);
  }
}

/****************************************************************************/
/*!
 **  Calculate the covariance between data and data
 **
 ** \return  The covariance terms or TEST
 **
 ** \param[in]  lmlayers  Pointer to the LMlayers structure to be freed
 ** \param[in]  model     Model
 ** \param[in]  ilayer    Layer of interest (first point). Starting from 1
 ** \param[in]  prop1     Working array at first point (Dimension: nlayers)
 ** \param[in]  jlayer    Layer of interest (second point). Starting from 1
 ** \param[in]  prop2     Working array at second point (Dimension: nlayers)
 ** \param[in]  dd        Distance vector (or NULL for zero-distance)
 **
 ** \param[out] covtab    Working array (Dimension = nlayers * nlayers)
 **
 ** \remarks:  As this function may return TEST, TEST value should be tested
 **
 *****************************************************************************/
static double st_cij(LMlayers *lmlayers,
                     Model *model,
                     int ilayer,
                     double *prop1,
                     int jlayer,
                     double *prop2,
                     double *dd,
                     double *covtab)
{
  double value;
  int i, j, nlayers;
  VectorDouble d1;
  CovCalcMode mode;

  /* Initializations */

  nlayers = lmlayers->nlayers;
  d1.resize(2);
  st_check_layer("st_cij", lmlayers, ilayer);
  st_check_layer("st_cij", lmlayers, jlayer);

  /* Calculate the covariance matrix */

  d1[0] = (dd != nullptr) ? dd[0] :
                            0.;
  d1[1] = (dd != nullptr) ? dd[1] :
                            0.;
  model_calcul_cov(model, mode, 1, 1., d1, covtab);

  /* Evaluate the covariance term */

  value = 0.;
  for (i = 0; i < ilayer; i++)
    for (j = 0; j < jlayer; j++)
    {
      if (FFFF(prop1[i])) return (TEST);
      if (FFFF(prop2[j])) return (TEST);
      value += prop1[i] * prop2[j] * C(i, j);
    }
  return (value);
}

/****************************************************************************/
/*!
 **  Calculate the covariance between data and target
 **
 ** \return  The covariance terms or TEST
 **
 ** \param[in]  lmlayers  LMlayers structure
 ** \param[in]  model     Model
 ** \param[in]  ilayer    Layer of interest (data point). Starting from 1
 ** \param[in]  prop1     Working array at data point (Dimension: nlayers)
 ** \param[in]  jlayer    Layer of interest (target point). Starting from 1
 ** \param[in]  dd        Distance vector (or NULL for zero-distance)
 **
 ** \param[out] covtab    Working array (Dimension = nlayers * nlayers)
 **
 ** \remarks:  As this function may return TEST, TEST value should be tested
 **
 *****************************************************************************/
static double st_ci0(LMlayers *lmlayers,
                     Model *model,
                     int ilayer,
                     double *prop1,
                     int jlayer,
                     double *dd,
                     double *covtab)
{
  double value;
  int i, nlayers;
  VectorDouble d1;
  CovCalcMode mode;

  /* Initializations */

  nlayers = lmlayers->nlayers;
  d1.resize(2);
  st_check_layer("st_ci0", lmlayers, ilayer);
  st_check_layer("st_ci0", lmlayers, jlayer);

  /* Calculate the covariance matrix */

  d1[0] = (dd != nullptr) ? dd[0] :
                            0.;
  d1[1] = (dd != nullptr) ? dd[1] :
                            0.;
  model_calcul_cov(model, mode, 1, 1., d1, covtab);

  /* Evaluate the covariance term */

  value = 0.;
  for (i = 0; i < ilayer; i++)
  {
    if (FFFF(prop1[i])) return (1);
    value += prop1[i] * C(i, jlayer - 1);
  }
  return (value);
}

/****************************************************************************/
/*!
 **  Calculate the drift terms
 **
 ** \return  Error return code
 **
 ** \param[in]  lmlayers    Pointer to the LMlayers structure to be freed
 ** \param[in]  coor        Array of coordinates
 ** \param[in]  propval     Value for the proportion (used if flag_cumul=TRUE)
 ** \param[in]  drext       Value of the external drift
 ** \param[in,out] ipos_loc Address for the first drift term.
 **                         On output, address for the next term after the drift
 **
 ** \param[out] b         Array for storing the drift
 **
 *****************************************************************************/
static int st_drift(LMlayers *lmlayers,
                    double *coor,
                    double propval,
                    double drext,
                    int *ipos_loc,
                    double *b)
{
  int ipos;

  if (lmlayers->flag_ext && FFFF(drext)) return (1);
  ipos = *ipos_loc;
  switch (lmlayers->nbfl)
  {
    case 0:
      break;

    case 1:
      b[ipos++] = propval;
      break;

    case 2:
      b[ipos++] = propval;
      b[ipos++] = propval * drext;
      break;

    case 3:
      b[ipos++] = propval;
      b[ipos++] = propval * coor[0];
      b[ipos++] = propval * coor[1];
      break;

    case 4:
      b[ipos++] = propval;
      b[ipos++] = propval * coor[0];
      b[ipos++] = propval * coor[1];
      b[ipos++] = propval * drext;
      break;
  }
  *ipos_loc = ipos;
  return (0);
}

/****************************************************************************/
/*!
 **  Calculates the L.H.S. for one data
 **
 ** \returns:  1 if the L.H.S. vector cannot be calculated; 0 otherwise
 **
 ** \param[in]  lmlayers  LMlayers structure
 ** \param[in]  dbin      Input Db structure
 ** \param[in]  dbout     Output Db structure
 ** \param[in]  model     Model
 ** \param[in]  seltab    Number of sample definition (0, 1 or 2)
 ** \param[in]  iech0     Rank of the target sample (used for ext. drift)
 ** \param[in]  ilayer0   Rank of the layer of interest (Starting from 1)
 ** \param[in]  coor      Coordinates of the data
 ** \param[in]  prop0     Working array at first point (Dimension: nlayers)
 **
 ** \param[out] prop2     Working array (Dimension: nlayers)
 ** \param[out] covtab    Working array (Dimension = nlayers * nlayers)
 ** \param[out] b         R.H.S. vector (Dimension = neq)
 **
 *****************************************************************************/
static int st_lhs_one(LMlayers *lmlayers,
                      Db *dbin,
                      Db *dbout,
                      Model *model,
                      int *seltab,
                      int iech0,
                      int ilayer0,
                      double *coor,
                      double *prop0,
                      double *prop2,
                      double *covtab,
                      double *b)
{
  int jech, jjech, jfois, jlayer, nlayers, i;
  double drext, coor2[2], d1[2];

  /* Initializations */

  nlayers = lmlayers->nlayers;

  /* Covariance part */

  for (jech = jjech = 0; jech < dbin->getSampleNumber(); jech++)
  {
    if (seltab[jech] == 0) continue;
    coor2[0] = dbin->getCoordinate(jech, 0);
    coor2[1] = dbin->getCoordinate(jech, 1);
    for (jfois = 0; jfois < seltab[jech]; jfois++, jjech++)
    {
      jlayer =
          (jfois == 0) ? (int) get_LOCATOR_ITEM(dbin, ELoc::LAYER, 0, jech) :
                         nlayers;

      /* Evaluate the proportion vector */

      if (st_get_props_data(lmlayers, dbin, dbout, jech, jlayer, prop2))
        return (1);

      /* Calculate the distance vector */

      d1[0] = coor[0] - coor2[0];
      d1[1] = coor[1] - coor2[1];
      b[jjech] = st_cij(lmlayers, model, ilayer0, prop0, jlayer, prop2, d1,
                        covtab);
      if (FFFF(b[jjech])) return (1);
    }
  }

  /* Drift part */

  for (i = 0; i < ilayer0; i++)
  {
    drext = st_get_drift_data(lmlayers, dbin, dbout, iech0, i + 1);
    if (st_drift(lmlayers, coor, prop0[i], drext, &jjech, b)) return (1);
  }
  return (0);
}

/****************************************************************************/
/*!
 **  Calculates the Kriging R.H.S.
 **
 ** \returns:  1 if the R.H.S. has not been calculated; 0 otherwise
 **
 ** \param[in]  lmlayers  LMlayers structure
 ** \param[in]  dbin      Input Db structure
 ** \param[in]  dbout     Output Db structure
 ** \param[in]  model     Model
 ** \param[in]  coor      Coordinates of the target sample
 ** \param[in]  seltab    Number of sample definition (0, 1 or 2)
 ** \param[in]  iechout   Rank of the target sample
 ** \param[in]  ilayer0   Rank of the layer of interest (Starting from 1)
 **
 ** \param[out] prop0     Working array (Dimension: nlayers)
 ** \param[out] prop2     Working array (Dimension: nlayers)
 ** \param[out] covtab    Working array (Dimension = nlayers * nlayers)
 ** \param[out] b         R.H.S. vector (Dimension = neq)
 **
 *****************************************************************************/
static int st_rhs(LMlayers *lmlayers,
                  Db *dbin,
                  Db *dbout,
                  Model *model,
                  double *coor,
                  int *seltab,
                  int iechout,
                  int ilayer0,
                  double *prop0,
                  double *prop2,
                  double *covtab,
                  double *b)
{
  int jech, jjech, i, jlayer, ipos, ifois, nlayers, ideb;
  double drext, d1[2], coor2[2], propval;

  /* Get the coordinates of the target */

  nlayers = lmlayers->nlayers;
  st_check_layer("st_rhs", lmlayers, ilayer0);
  coor[0] = dbout->getCoordinate(iechout, 0);
  coor[1] = dbout->getCoordinate(iechout, 1);

  /* Initialize the vector with zeroes */

  for (i = 0; i < lmlayers->neq; i++)
    b[i] = 0.;

  /* Covariance part */

  for (jech = jjech = 0; jech < dbin->getSampleNumber(); jech++)
  {
    if (seltab[jech] == 0) continue;
    coor2[0] = dbin->getCoordinate(jech, 0);
    coor2[1] = dbin->getCoordinate(jech, 1);
    for (ifois = 0; ifois < seltab[jech]; ifois++, jjech++)
    {
      jlayer =
          (ifois == 0) ? (int) get_LOCATOR_ITEM(dbin, ELoc::LAYER, 0, jech) :
                         nlayers;

      /* Evaluate the proportion vector */

      (void) st_get_props_data(lmlayers, dbin, dbout, jech, jlayer, prop2);

      /* Calculate the distance vector */

      d1[0] = coor2[0] - coor[0];
      d1[1] = coor2[1] - coor[1];
      if (lmlayers->flag_cumul)
        b[jjech] = st_cij(lmlayers, model, ilayer0, prop0, jlayer, prop2, d1,
                          covtab);
      else
        b[jjech] = st_ci0(lmlayers, model, jlayer, prop2, ilayer0, d1, covtab);
      if (FFFF(b[jjech])) return (1);
    }
  }

  /* Drift part */

  ideb = (lmlayers->flag_cumul) ? 0 :
                                  ilayer0 - 1;
  for (i = ideb; i < ilayer0; i++)
  {
    ipos = lmlayers->nech + lmlayers->nbfl * i;
    drext = st_get_drift_result(lmlayers, dbout, iechout, i + 1);
    propval = (lmlayers->flag_cumul) ? prop0[i] :
                                       1.;
    if (st_drift(lmlayers, coor, propval, drext, &ipos, b)) return (1);
  }
  return (0);
}

/****************************************************************************/
/*!
 **  Calculates the Kriging L.H.S.
 **
 ** \returns:  1 if the L.H.S. vector cannot be calculated; 0 otherwise
 **
 ** \param[in]  lmlayers  LMlayers structure
 ** \param[in]  dbin      Input Db structure
 ** \param[in]  dbout     Output Db structure
 ** \param[in]  model     Model
 ** \param[in]  seltab    Number of sample definition (0, 1 or 2)
 **
 ** \param[out] prop1     Working array (Dimension: nlayers)
 ** \param[out] prop2     Working array (Dimension: nlayers)
 ** \param[out] covtab    Working array (Dimension = nlayers * nlayers)
 ** \param[out] a         L.H.S. (square) matrix
 ** \param[out] acov      L.H.S. (square) covariance matrix
 **
 *****************************************************************************/
static int st_lhs(LMlayers *lmlayers,
                  Db *dbin,
                  Db *dbout,
                  Model *model,
                  int *seltab,
                  double *prop1,
                  double *prop2,
                  double *covtab,
                  double *a,
                  double *acov)
{
  int    iiech,jjech;
  double coor[2];

  /* Initialize the matrix with zeroes */

  int nech    = lmlayers->nech;
  int neq     = lmlayers->neq;
  int nlayers = lmlayers->nlayers;
  for (int i=0; i<neq * neq; i++) a[i] = 0.;

  /* Loop on the first sample */

  iiech = 0;
  for (int iech=0; iech<dbin->getSampleNumber(); iech++)
  {
    if (seltab[iech] == 0) continue;
    coor[0] = dbin->getCoordinate(iech,0);
    coor[1] = dbin->getCoordinate(iech,1);
    for (int ifois=0; ifois<seltab[iech]; ifois++, iiech++)
    {
      int ilayer = (ifois == 0) ?
        (int) get_LOCATOR_ITEM(dbin,ELoc::LAYER,0,iech) : nlayers;
      
      /* Evaluate the proportion vector */

      if (st_get_props_data(lmlayers, dbin, dbout, iech, ilayer, prop1))
        return (1);

      /* Loop on the second sample */

      if (st_lhs_one(lmlayers, dbin, dbout, model, seltab, iech, ilayer, coor,
                     prop1, prop2, covtab, &A(iiech, 0))) return (1);
    }
  }

  /* Symmetrization */

  for (iiech = 0; iiech < neq; iiech++)
    for (jjech = 0; jjech <= iiech; jjech++)
      A(iiech,jjech) = A(jjech, iiech);

  /* Extraction of the Covariance matrix */

  for (int iech=0; iech<nech; iech++)
    for (int jech=0; jech<nech; jech++)
      ACOV(iech,jech) = A(iech,jech);
  
  if (get_keypone("Bayes_Debug_Flag",0))
    set_keypair("Mlayers_LHS_Matrix",1,neq,neq,a);

  return(0);
}

/****************************************************************************/
/*!
 **  Establish the vector of data
 **
 ** \param[in]  lmlayers  LMlayers structure
 ** \param[in]  dbin      Input Db structure
 ** \param[in]  dbout     Output Db structure
 ** \param[in]  seltab    Number of sample definition (0, 1 or 2)
 **
 ** \param[out] zval      The data vector (Dimension: neq)
 **
 *****************************************************************************/
static void st_data_vector(LMlayers *lmlayers,
                           Db *dbin,
                           Db *dbout,
                           int *seltab,
                           double *zval)
{
  double value, dtime;
  int i, iech, iiech, igrid, ifois, ilayer, nlayers;

  /* Initialize the vector with zeroes */

  igrid = 0;
  nlayers = lmlayers->nlayers;
  for (i = 0; i < lmlayers->neq; i++)
    zval[i] = 0.;

  /* Loop on the samples */

  for (iech = iiech = 0; iech < dbin->getSampleNumber(); iech++)
  {
    if (seltab[iech] == 0) continue;

    /* Calculate the grid node index (optional) */

    if (lmlayers->colrefd >= 0 || lmlayers->colreft >= 0
        || lmlayers->colrefb >= 0 || lmlayers->flag_vel)
      (void) st_locate_sample_in_output(lmlayers, dbin, dbout, iech, &igrid);

    for (ifois = 0; ifois < seltab[iech]; ifois++, iiech++)
    {
      ilayer =
          (ifois == 0) ? (int) get_LOCATOR_ITEM(dbin, ELoc::LAYER, 0, iech) :
                         nlayers;

      if (ifois == 0)
      {

        /* Depth of the actual sample */

        value = dbin->getVariable(iech, 0);
      }
      else
      {

        /* Depth of the collocated bottom sample */

        value = dbout->getArray(igrid, lmlayers->colrefb);
      }

      /* Centering to the reference Depth surface */

      if (lmlayers->colrefd >= 0)
        value -= dbout->getArray(igrid, lmlayers->colrefd);

      /* Converting into velocities */

      if (lmlayers->flag_vel)
      {
        dtime = get_LOCATOR_ITEM(dbout, lmlayers->ptime, ilayer - 1, igrid);
        if (lmlayers->colreft >= 0)
          dtime -= dbout->getArray(igrid, lmlayers->colreft);
        value /= dtime;
      }
      zval[iiech] = value;
    }
  }

  if (get_keypone("Bayes_Debug_Flag", 0))
    set_keypair("Mlayers_Zval_Matrix", 1, lmlayers->neq, 1, zval);
}

/****************************************************************************/
/*!
 **  Calculate the Drift and subtract it from the Data
 **
 ** \param[in]  lmlayers  LMlayers structure
 ** \param[in]  verbose   1 for a  verbose option
 ** \param[in]  dbin      Input Db structure
 ** \param[in]  dbout     Output Db structure
 ** \param[in]  seltab    Number of sample definition (0, 1 or 2)
 **
 ** \param[out] zval      The data vector (Dimension: neq)
 **
 ** \remarks In the current version, the optimal coefficients of the Drift
 ** \remarks are output using the set_keypair mechanism using the keyword:
 ** \remarks "Optim_Drift_Coeffs" which returns 'ipos' values
 **
 *****************************************************************************/
static int st_subtract_optimal_drift(LMlayers *lmlayers,
                                     int verbose,
                                     Db *dbin,
                                     Db *dbout,
                                     int *seltab,
                                     double *zval)
{
  double *atab, *btab, *drift, *props, *coeff, drext, coor[2];
  int nlayers, error, iech, iiech, ifois, ilayer, nbfl, neq, ipos;
  int flag_subtract = 1;

  /* Initializations */

  error = 1;
  nlayers = lmlayers->nlayers;
  drift = props = atab = btab = coeff = nullptr;
  nbfl = lmlayers->nbfl;
  neq = nbfl * nlayers;

  /* Core allocation */

  coeff = (double*) mem_alloc(sizeof(double) * neq, 0);
  if (coeff == nullptr) goto label_end;
  drift = (double*) mem_alloc(sizeof(double) * neq, 0);
  if (drift == nullptr) goto label_end;
  props = (double*) mem_alloc(sizeof(double) * nlayers, 0);
  if (props == nullptr) goto label_end;
  atab = (double*) mem_alloc(sizeof(double) * neq * neq, 0);
  if (atab == nullptr) goto label_end;
  btab = (double*) mem_alloc(sizeof(double) * neq, 0);
  if (btab == nullptr) goto label_end;
  for (int i = 0; i < neq; i++)
    btab[i] = 0.;
  for (int i = 0; i < neq * neq; i++)
    atab[i] = 0.;

  /* Find the vector of optimal mean values */

  for (iech = iiech = 0; iech < dbin->getSampleNumber(); iech++)
  {
    if (seltab[iech] == 0) continue;

    for (ifois = 0; ifois < seltab[iech]; ifois++, iiech++)
    {
      ilayer =
          (ifois == 0) ? (int) get_LOCATOR_ITEM(dbin, ELoc::LAYER, 0, iech) :
                         nlayers;

      /* Evaluate the proportion vector */

      if (st_get_props_data(lmlayers, dbin, dbout, iech, ilayer, props))
        goto label_end;

      /* Get the coordinates of the samples */

      coor[0] = dbin->getCoordinate(iech, 0);
      coor[1] = dbin->getCoordinate(iech, 1);

      /* Get the drift vector */

      ipos = 0;
      for (int i = 0; i < nlayers; i++)
      {
        drext = st_get_drift_data(lmlayers, dbin, dbout, iech, i + 1);
        if (st_drift(lmlayers, coor, props[i], drext, &ipos, drift)) continue;
      }

      /* Calculate the contribution to the different arrays */

      for (int k1 = 0; k1 < ipos; k1++)
      {
        btab[k1] += zval[iiech] * drift[k1];
        for (int k2 = 0; k2 < ipos; k2++)
          ATAB(neq,k1,k2) += drift[k1] * drift[k2];
      }
    }
  }

  /* Find the optimal drift coefficients */

  if (matrix_invert(atab, neq, -1)) goto label_end;
  matrix_product(neq, neq, 1, atab, btab, coeff);

  /* Optional printout of the result */

  if (verbose)
    print_matrix("Estimated Drift", 0, 1, nlayers, nbfl, NULL, coeff);

  /* Subtract the optimal mean */

  for (iech = iiech = 0; iech < dbin->getSampleNumber(); iech++)
  {
    if (seltab[iech] == 0) continue;

    for (ifois = 0; ifois < seltab[iech]; ifois++, iiech++)
    {
      ilayer =
          (ifois == 0) ? (int) get_LOCATOR_ITEM(dbin, ELoc::LAYER, 0, iech) :
                         nlayers;

      /* Evaluate the proportion vector */

      if (st_get_props_data(lmlayers, dbin, dbout, iech, ilayer, props))
        goto label_end;

      /* Get the coordinates of the samples */

      coor[0] = dbin->getCoordinate(iech, 0);
      coor[1] = dbin->getCoordinate(iech, 1);

      /* Get the drift vector */

      ipos = 0;
      for (int i = 0; i < nlayers; i++)
      {
        drext = st_get_drift_data(lmlayers, dbin, dbout, iech, i + 1);
        if (st_drift(lmlayers, coor, props[i], drext, &ipos, drift)) continue;
      }

      /* Subtract the optimal estimation of the drift */

      if (flag_subtract) for (int k1 = 0; k1 < ipos; k1++)
        zval[iiech] -= coeff[k1] * drift[k1];

      /* Save the results of the optimal drift */

      set_keypair("Optim_Drift_Coeffs", 1, 1, ipos, coeff);

      /* Print the residuals (optional) */

      if (debug_query("variogram"))
        message("Sample %d (Layer %d) - Coor = %lf %lf - Residual = %lf\n",
                iech + 1, ilayer, coor[0], coor[1], zval[iiech]);
    }
  }

  /* Set the error return code */

  error = 0;

  label_end:

  /* Core deallocation */

  atab = (double*) mem_free((char* ) atab);
  btab = (double*) mem_free((char* ) btab);
  props = (double*) mem_free((char* ) props);
  drift = (double*) mem_free((char* ) drift);
  return (error);
}

/****************************************************************************/
/*!
 **  Check if an intercept with the bottom layer is located close enough
 **  to the current sample
 **
 ** \return  1 if a duplicate must be generated; 0 otherwise
 **
 ** \param[in]  lmlayers  LMlayers structure
 ** \param[in]  dbin      Input Db structure
 ** \param[in]  iech0     Rank of sample to be discarded (or -1)
 ** \param[in]  coor      Coordinates of the target
 **
 *****************************************************************************/
static int st_get_close_sample(LMlayers *lmlayers,
                               Db *dbin,
                               int iech0,
                               double *coor)
{
  int iech, ilayer;
  double dx, dy;
  static double EPS = 1.e-05;

  /* Check if a close sample has already been reviewed */

  for (iech = 0; iech < iech0; iech++)
  {
    dx = dbin->getCoordinate(iech, 0) - coor[0];
    if (ABS(dx) > EPS) continue;
    dy = dbin->getCoordinate(iech, 1) - coor[1];
    if (ABS(dy) > EPS) continue;
    return (0);
  }

  /* Check among the subsequent samples if a sample with matching coordinates */
  /* and belonging to the bottom surface exists */

  for (iech = iech0 + 1; iech < dbin->getSampleNumber(); iech++)
  {
    dx = dbin->getCoordinate(iech, 0) - coor[0];
    if (ABS(dx) > EPS) continue;
    dy = dbin->getCoordinate(iech, 1) - coor[1];
    if (ABS(dy) > EPS) continue;
    ilayer = (int) get_LOCATOR_ITEM(dbin, ELoc::LAYER, 0, iech);
    if (ilayer == lmlayers->nlayers) return (0);
  }
  return (1);
}

/****************************************************************************/
/*!
 **  Perform the per-calculation for estimation with collocated option
 **
 ** \param[in]  lmlayers  LMlayers structure
 ** \param[in]  iechout   Rank of the target
 ** \param[in]  coor      Coordinates of the target
 ** \param[in]  dbin      Input Db structure
 ** \param[in]  dbout     Output Db structure
 ** \param[in]  model     Model
 ** \param[in]  seltab    Number of sample definition (0, 1 or 2)
 ** \param[in]  a         L.H.S. (square) inverted matrix
 ** \param[in]  zval      Data vector (extended)
 **
 ** \param[out] prop1     Working array (Dimension: nlayers)
 ** \param[out] prop2     Working array (Dimension: nlayers)
 ** \param[out] covtab    Working array (Dimension = nlayers * nlayers)
 ** \param[out] b2        Working vector (Dimension = neq)
 ** \param[out] baux      Working vector (Dimension = neq)
 ** \param[out] ratio     Ratio value
 **
 *****************************************************************************/
static int st_collocated_prepare(LMlayers *lmlayers,
                                 int iechout,
                                 double *coor,
                                 Db *dbin,
                                 Db *dbout,
                                 Model *model,
                                 int *seltab,
                                 double *a,
                                 double *zval,
                                 double *prop1,
                                 double *prop2,
                                 double *covtab,
                                 double *b2,
                                 double *baux,
                                 double *ratio)
{
  double botval, c0, coefa, coefz;
  int nlayers, neq;

  (*ratio) = 0.;
  neq = lmlayers->neq;
  nlayers = lmlayers->nlayers;

  botval = dbout->getArray(iechout, lmlayers->colrefb);
  if (lmlayers->colrefd >= 0)
    botval -= dbout->getArray(iechout, lmlayers->colrefd);
  if (st_get_props_result(lmlayers, dbout, iechout, nlayers, prop1)) return (1);
  c0 = st_cij(lmlayers, model, nlayers, prop1, nlayers, prop1, NULL, covtab);
  if (FFFF(c0)) return (1);

  if (st_lhs_one(lmlayers, dbin, dbout, model, seltab, iechout, nlayers, coor,
                 prop1, prop2, covtab, baux)) return (1);
  matrix_product(neq, neq, 1, a, baux, b2);
  matrix_product(1, neq, 1, b2, zval, &coefz);
  matrix_product(1, neq, 1, b2, baux, &coefa);
  (*ratio) = (ABS(c0 - coefa) > 1.e-6) ? (botval - coefz) / (c0 - coefa) :
                                         0.;

  return (0);
}

/****************************************************************************/
/*!
 **  Perform the estimation at the grid nodes in regular case
 **
 ** \param[in]  lmlayers  LMlayers structure
 ** \param[in]  flag_std  1 if the estimation error must be calculated
 ** \param[in]  c00       Variance for target
 ** \param[in]  a         L.H.S. (square) inverted matrix
 ** \param[in]  b         Working vector (Dimension = neq)
 ** \param[in]  dual      Dual vector
 ** \param[in]  wgt       Working array (Dimension = neq)
 **
 ** \param[out] estim     Estimated value
 ** \param[out] stdev     Standard deviation of estimation error
 **
 *****************************************************************************/
static void st_estimate_regular(LMlayers *lmlayers,
                                int flag_std,
                                double c00,
                                double *a,
                                double *b,
                                double *dual,
                                double *wgt,
                                double *estim,
                                double *stdev)
{
  double c00val, stdv;
  int neq;

  /* Initializations */

  neq = lmlayers->neq;
  *estim = *stdev = TEST;

  /* Perform the estimation (in Dual form) */

  matrix_product(1, neq, 1, dual, b, estim);

  /* Perform the variance of estimation error */

  if (flag_std)
  {
    c00val = c00;
    if (FFFF(c00val))
      stdv = TEST;
    else
    {
      matrix_product(neq, neq, 1, a, b, wgt);
      matrix_product(1, neq, 1, b, wgt, &stdv);
      stdv = c00val - stdv;
      stdv = (stdv > 0) ? sqrt(stdv) :
                          0.;
    }
    *stdev = stdv;
  }
}

/****************************************************************************/
/*!
 **  Perform the estimation at the grid nodes in the bayesian case
 **
 ** \param[in]  lmlayers  LMlayers structure
 ** \param[in]  flag_std  1 if the estimation error must be calculated
 ** \param[in]  c00       Variance at target
 ** \param[in]  acov      L.H.S. (square) inverted matrix
 ** \param[in]  zval      Data vector
 ** \param[in]  b         Working vector (Dimension = neq)
 ** \param[in]  wgt       Working array (Dimension = neq)
 ** \param[out] post_mean  Array of posterior mean
 ** \param[out] a0         Constant term
 ** \param[out] cc         Output value
 ** \param[out] ss         Output value
 ** \param[out] gs         Output value
 **
 ** \param[out] estim     Estimated value
 ** \param[out] stdev     Standard deviation of estimation error
 **
 *****************************************************************************/
static void st_estimate_bayes(LMlayers *lmlayers,
                              int flag_std,
                              double c00,
                              double *acov,
                              double *zval,
                              double *b,
                              double *wgt,
                              double *post_mean,
                              double *a0,
                              double *cc,
                              double *ss,
                              double *gs,
                              double *estim,
                              double *stdev)
{
  double *rhs, *ff0, *temp, *fsf0, *c2, estim1, estim2, stdv;
  int nech, npar;

  /* Initializations */

  *estim = *stdev = TEST;
  nech = lmlayers->nech;
  npar = lmlayers->npar;
  rhs = &b[0];
  ff0 = &b[nech];

  /* Core allocation */

  temp = (double*) mem_alloc(sizeof(double) * npar, 1);
  fsf0 = (double*) mem_alloc(sizeof(double) * nech, 1);
  c2 = (double*) mem_alloc(sizeof(double) * nech, 1);

  /* Perform the estimation */

  matrix_product(nech, npar, 1, a0, ff0, fsf0);
  for (int iech = 0; iech < nech; iech++)
    c2[iech] = rhs[iech] + fsf0[iech];
  matrix_product(nech, nech, 1, cc, c2, wgt);

  matrix_product(1, nech, 1, wgt, zval, &estim1);
  matrix_product(1, npar, 1, ff0, post_mean, &estim2);
  *estim = estim1 + estim2;

  /* Calculate the standard deviation */

  if (flag_std)
  {
    matrix_product(1, nech, npar, rhs, ss, temp);
    for (int ipar = 0; ipar < npar; ipar++)
      temp[ipar] -= ff0[ipar];

    stdv = c00;
    for (int iech = 0; iech < nech; iech++)
      for (int jech = 0; jech < nech; jech++)
        stdv -= rhs[iech] * ACOV(iech, jech) * rhs[jech];

    for (int ipar = 0; ipar < npar; ipar++)
      for (int jpar = 0; jpar < npar; jpar++)
        stdv += temp[ipar] * GS(ipar, jpar) * temp[jpar];

    stdv = (stdv > 0) ? sqrt(stdv) :
                        0.;
    *stdev = stdv;
  }

  /* Core deallocation */

  temp = (double*) mem_free((char* ) temp);
  fsf0 = (double*) mem_free((char* ) fsf0);
  c2 = (double*) mem_free((char* ) c2);
}

/****************************************************************************/
/*!
**  Perform the estimation at the grid nodes
**
** \param[in]  lmlayers   LMlayers structure
** \param[in]  dbin       Input Db structure
** \param[in]  dbout      Output Db structure
** \param[in]  model      Model
** \param[in]  seltab     Number of sample definition (0, 1 or 2)
** \param[in]  flag_bayes 1 if the Bayesian hypothesis is used on drift coeffs
** \param[in]  flag_std   1 if the estimation error must be calculated
** \param[in]  a          L.H.S. (square) inverted matrix
** \param[in]  zval       Data vector (extended)
** \param[in]  dual       Dual vector
**
** \param[out] prop1      Working array (Dimension: nlayers) 
** \param[out] prop2      Working array (Dimension: nlayers) 
** \param[out] covtab     Working array (Dimension = nlayers * nlayers)
** \param[out] b          Working vector (Dimension = neq)
** \param[out] b2         Working vector (Dimension = neq)
** \param[out] baux       Working vector (Dimension = neq)
** \param[out] wgt        Working array (Dimension = neq)
** \param[out] c00        Working array (Dimension = nlayers)
** \param[out] a0         Constant term
** \param[out] cc         Output value
** \param[out] ss         Output value
** \param[out] gs         Output value
** \param[out] post_mean  Array of posterior mean
**
*****************************************************************************/
static void st_estimate(LMlayers *lmlayers,
                        Db *dbin,
                        Db *dbout,
                        Model *model,
                        int *seltab,
                        int flag_bayes,
                        int flag_std,
                        double *a,
                        double *zval,
                        double *dual,
                        double *prop1,
                        double *prop2,
                        double *covtab,
                        double *b,
                        double *b2,
                        double *baux,
                        double *wgt,
                        double *c00,
                        double * /*fftab*/,
                        double *a0,
                        double *cc,
                        double *ss,
                        double *gs,
                        double * /*prior_mean*/,
                        double *post_mean)
{
  double estim, cx, coor[2], coefb, botval, ratio, stdv;
  int iechout, ilayer, flag_correc, nlayers, neq;

  /* Loop on the grid nodes */

  nlayers = lmlayers->nlayers;
  neq = lmlayers->neq;
  coefb = ratio = botval = stdv = 0.;
  if (flag_std && !lmlayers->flag_cumul)
    st_covariance_c00(lmlayers, model, NULL, covtab, c00);

  for (iechout = 0; iechout < dbout->getSampleNumber(); iechout++)
  {
    debug_index(iechout + 1);
    if (!dbout->isActive(iechout)) continue;
    coor[0] = dbout->getCoordinate(iechout, 0);
    coor[1] = dbout->getCoordinate(iechout, 1);
    if (debug_query("kriging") || debug_query("nbgh") || debug_query("results"))
    {
      mestitle(1, "Target location");
      db_sample_print(dbout, iechout, 1, 0, 0);
    }

    /* Correction in the case of collocation of the bottom surface */

    flag_correc = 0;
    if (lmlayers->colrefb >= 0)
    {
      flag_correc = st_get_close_sample(lmlayers, dbin, -1, coor);
      if (flag_correc)
      {
        if (st_collocated_prepare(lmlayers, iechout, coor, dbin, dbout, model,
                                  seltab, a, zval, prop1, prop2, covtab, b2,
                                  baux, &ratio)) continue;
      }
    }

    /* Loop on the layers */

    for (ilayer = 0; ilayer < nlayers; ilayer++)
    {
      if (debug_query("kriging") || debug_query("nbgh"))
        mestitle(2, "Layer #%d", ilayer + 1);

      /* Find the proportions for the target if flag_cumul=TRUE */

      if (lmlayers->flag_cumul)
      {
        if (st_get_props_result(lmlayers, dbout, iechout, ilayer + 1, prop1))
          continue;
      }

      /* Find the C00 terms */

      if (flag_std && lmlayers->flag_cumul)
        st_covariance_c00(lmlayers, model, prop1, covtab, c00);

      /* Establish the R.H.S. */

      if (st_rhs(lmlayers, dbin, dbout, model, coor, seltab, iechout,
                 ilayer + 1, prop1, prop2, covtab, b)) continue;
      if (debug_query("kriging"))
        krige_rhs_print(1, lmlayers->nech, neq, neq, NULL, b);

      /* Perform estimation */

      if (flag_bayes)
        st_estimate_bayes(lmlayers, flag_std, c00[ilayer], a, zval, b, wgt,
                          post_mean, a0, cc, ss, gs, &estim, &stdv);
      else
        st_estimate_regular(lmlayers, flag_std, c00[ilayer], a, b, dual, wgt,
                            &estim, &stdv);

      /* Perform the correction (in case of collocated bottom) */

      if (flag_correc)
      {
        cx = st_ci0(lmlayers, model, nlayers, prop1, ilayer + 1, NULL, covtab);
        if (FFFF(cx)) continue;
        matrix_product(1, neq, 1, b2, b, &coefb);
        estim += (cx - coefb) * ratio;
      }

      /* Store the result */

      dbout->setVariable(iechout, ilayer, estim);
      if (flag_std) dbout->setVariable(iechout, nlayers + ilayer, stdv);
      if (debug_query("results"))
      {
        message("Estimate = %lf", ilayer + 1, estim);
        if (flag_std) message(" - Variance = %lf", stdv * stdv);
        message("\n");
      }
    }
  }
  debug_index(0);
}

/****************************************************************************/
/*!
 **  Check all the auxiliary variables
 **
 ** \return The number of active samples
 **
 ** \param[in]  lmlayers  LMlayers structure
 ** \param[in]  dbin      Input Db structure
 ** \param[in]  dbout     Output Db structure
 **
 ** \param[in,out]  seltab    Number of sample definition (0, 1 or 2)
 **
 *****************************************************************************/
static int st_check_auxiliary_variables(LMlayers *lmlayers,
                                        Db *dbin,
                                        Db *dbout,
                                        int *seltab)
{
  int iech, ilayer, igrid, newval, nechtot;
  double drift, value, coor[2];

  nechtot = 0;
  for (iech = 0; iech < dbin->getSampleNumber(); iech++)
  {
    if (seltab[iech] == 0) continue;
    coor[0] = dbin->getCoordinate(iech, 0);
    coor[1] = dbin->getCoordinate(iech, 1);
    ilayer = (int) get_LOCATOR_ITEM(dbin, ELoc::LAYER, 0, iech);
    if (st_locate_sample_in_output(lmlayers, dbin, dbout, iech, &igrid))
      goto label_suppress;

    /* Case of an external drift */

    if (lmlayers->flag_ext)
    {
      drift = st_get_drift_data(lmlayers, dbin, dbout, iech, ilayer);
      if (FFFF(drift)) goto label_suppress;
    }

    /* Case of a Depth Reference variable */

    if (lmlayers->colrefd >= 0)
    {
      value = dbout->getArray(igrid, lmlayers->colrefd);
      if (FFFF(value)) goto label_suppress;
    }

    /* Case of a Bottom Depth Reference variable */

    newval = 1;
    if (lmlayers->colrefb >= 0)
    {
      value = dbout->getArray(igrid, lmlayers->colrefb);
      if (FFFF(value)) goto label_suppress;

      /* Check if a duplicate sample must be added:       */
      /* - the sample does not belong to the bottom layer */
      /* - an analoguous sample does not already exist    */

      if (ilayer < lmlayers->nlayers && st_get_close_sample(lmlayers, dbin,
                                                            iech, coor))
        newval = 2;
    }

    /* The sample is valid */

    seltab[iech] = newval;
    nechtot += newval;
    continue;

    label_suppress: seltab[iech] = 0;
    continue;
  }

  return (nechtot);
}

/****************************************************************************/
/*!
 **  Convert the results in the Depth scale
 **
 ** \param[in]  lmlayers  LMlayers structure
 ** \param[in]  dbout     Output Db structure
 ** \param[in]  flag_std  1 if the estimation error must be calculated
 **
 ** \remarks The conversion is performed:
 ** \remarks - if the calculations have been performed in Velocity or Thickness
 ** \remarks - if the calculations have been performed in cumulative or not
 ** \remarks The standard deviation is also transformed
 **
 *****************************************************************************/
static void st_convert_results(LMlayers *lmlayers, Db *dbout, int flag_std)
{
  double depth0, depth, value, stdv, time0, time, depth_prev, time_prev, delta;
  int iechout, ilayer, nlayers;

  /* Initializations */

  nlayers = lmlayers->nlayers;
  time = depth = stdv = 0.;

  /* If Depth converion is not required, nothing to be done */

  if (!lmlayers->flag_z) return;

  /* Loop on the target points */

  for (iechout = 0; iechout < dbout->getSampleNumber(); iechout++)
  {

    /* Identify the reference surface */

    depth0 =
        (lmlayers->colrefd >= 0) ? dbout->getArray(iechout, lmlayers->colrefd) :
                                   0.;

    /* Identify the reference time (for velocity) */

    time0 =
        (lmlayers->colreft >= 0) ? dbout->getArray(iechout, lmlayers->colreft) :
                                   0.;

    depth_prev = depth0;
    time_prev = time0;

    /* Loop on the layers */

    for (ilayer = 0; ilayer < lmlayers->nlayers; ilayer++)
    {

      /* Read the estimated value */

      value = dbout->getVariable(iechout, ilayer);
      if (flag_std) stdv = dbout->getVariable(iechout, nlayers + ilayer);

      if (lmlayers->flag_cumul)
      {

        /* Case when calculations have been processed in cumulative way */

        if (lmlayers->flag_vel)
        {
          time = get_LOCATOR_ITEM(dbout, lmlayers->ptime, ilayer, iechout);
          delta = time - time0;
          depth = depth0 + value * delta;
          if (flag_std) stdv *= delta;
        }
        else
        {
          depth = depth0 + value;
        }
      }
      else
      {

        /* Case when calculations have been processed in individual way */

        if (lmlayers->flag_vel)
        {
          time = get_LOCATOR_ITEM(dbout, lmlayers->ptime, ilayer, iechout);
          delta = time - time_prev;
          depth = depth_prev + value * delta;
          if (flag_std) stdv *= delta;
        }
        else
        {
          depth = depth_prev + value;
        }
        time_prev = time;
        depth_prev = depth;
      }

      /* Store the transformed results */

      dbout->setVariable(iechout, ilayer, depth);
      if (flag_std) dbout->setVariable(iechout, nlayers + ilayer, stdv);
    }
  }
}

/****************************************************************************/
/*!
 **  Fill the array of drift values at data points
 **
 ** \returns 1 Error return code
 **
 ** \param[in]  lmlayers  LMlayers structure
 ** \param[in]  dbin      Input Db structure
 ** \param[in]  dbout     Output Db structure
 ** \param[in]  seltab    Number of sample definition (0, 1 or 2)
 ** \param[in]  prop1     Working array (Dimension: nlayers)
 **
 ** \param[out] fftab     Drift array (Dimension: npar[nrow] * nech[ncol])
 **
 *****************************************************************************/
static int st_drift_data(LMlayers *lmlayers,
                         Db *dbin,
                         Db *dbout,
                         int *seltab,
                         double *prop1,
                         double *fftab)
{
  int npar, nech, iech, iiech, ilayer, ipos;
  double coor[2], drext;

  /* Initializations */

  npar = lmlayers->npar;
  nech = lmlayers->nech;
  for (int i = 0; i < npar * nech; i++)
    fftab[i] = 0.;

  for (iech = iiech = 0; iech < dbin->getSampleNumber(); iech++)
  {
    if (seltab[iech] == 0) continue;
    coor[0] = dbin->getCoordinate(iech, 0);
    coor[1] = dbin->getCoordinate(iech, 1);
    for (int ifois = 0; ifois < seltab[iech]; ifois++, iiech++)
    {
      ilayer =
          (ifois == 0) ? (int) get_LOCATOR_ITEM(dbin, ELoc::LAYER, 0, iech) :
                         lmlayers->nlayers;

      /* Evaluate the proportion vector */

      if (st_get_props_data(lmlayers, dbin, dbout, iech, ilayer, prop1))
        return (1);

      /* Loop on the second sample */

      ipos = iech * lmlayers->npar;
      for (int i = 0; i < ilayer; i++)
      {
        drext = st_get_drift_data(lmlayers, dbin, dbout, iech, i + 1);
        if (st_drift(lmlayers, coor, prop1[i], drext, &ipos, fftab)) return (1);
      }
    }
  }

  if (get_keypone("Bayes_Debug_Flag", 0))
    set_keypair("Mlayers_Drift_Matrix", 1, nech, npar, fftab);
  return (0);
}

/****************************************************************************/
/*!
 **  Calculate the posterior mean and variances in Bayesian case
 **
 ** \return  Error return code
 **
 ** \remarks At the end of this function,
 ** \remarks invS   contains the inverse of prior_S
 ** \remarks post_S contains the inverse of post_S
 **
 *****************************************************************************/
static int st_drift_bayes(LMlayers *lmlayers,
                          int verbose,
                          double *prior_mean,
                          double *prior_vars,
                          double *acov,
                          double *zval,
                          double *fftab,
                          double *a0,
                          double *cc,
                          double *ss,
                          double *gs,
                          double *post_mean,
                          double *post_S)
{
  double *ffc, *fft, *fm1z, *gg, *invH, *invS;
  int error, npar, nech, npar2, nech2;

  /* Initializations */

  error = 1;
  nech = lmlayers->nech;
  npar = lmlayers->npar;
  npar2 = npar * npar;
  nech2 = nech * nech;
  fft = ffc = fm1z = gg = invH = invS = nullptr;

  /* Core allocation */

  fft = (double*) mem_alloc(sizeof(double) * npar * nech, 0);
  if (fft == nullptr) goto label_end;
  ffc = (double*) mem_alloc(sizeof(double) * npar * nech, 0);
  if (ffc == nullptr) goto label_end;
  fm1z = (double*) mem_alloc(sizeof(double) * npar, 0);
  if (fm1z == nullptr) goto label_end;
  gg = (double*) mem_alloc(sizeof(double) * npar * npar, 0);
  if (gg == nullptr) goto label_end;
  invH = (double*) mem_alloc(sizeof(double) * npar * npar, 0);
  if (invH == nullptr) goto label_end;
  invS = (double*) mem_alloc(sizeof(double) * npar * npar, 0);
  if (invS == nullptr) goto label_end;

  /* Constitute the prior Variance-Covariance matrix */

  for (int ipar = 0; ipar < npar; ipar++)
    for (int jpar = 0; jpar < npar; jpar++)
      INVS(npar,ipar,jpar) = (ipar == jpar) ? prior_vars[ipar] :
                                              0.;

  /* Optional printout */

  if (verbose)
  {
    print_matrix("Prior Mean", 0, 1, lmlayers->nlayers, lmlayers->nbfl,
    NULL,
                 prior_mean);
    print_matrix("Prior Variance", 0, 1, lmlayers->npar, lmlayers->npar,
    NULL,
                 invS);
  }

  /* Invert the Data Variance-Covariance matrix */

  if (get_keypone("Bayes_Debug_Flag", 0))
    set_keypair("Bayes_ACov_Matrix", 1, nech, nech, acov);
  if (matrix_invert(acov, nech, -1)) goto label_end;

  /* Invert the prior Variance-Covariance matrix */

  if (get_keypone("Bayes_Debug_Flag", 0))
    set_keypair("Bayes_S_Matrix", 1, npar, npar, invS);
  if (matrix_invert(invS, npar, -1)) goto label_end;

  /* Auxiliary calculations */

  matrix_product(npar, nech, nech, fftab, acov, ffc);
  matrix_transpose(npar, nech, fftab, fft);
  matrix_product(npar, nech, npar, ffc, fft, invH);
  matrix_product(npar, nech, 1, ffc, zval, fm1z);
  if (get_keypone("Bayes_Debug_Flag", 0))
    set_keypair("Bayes_InvH_Matrix", 1, npar, npar, invH);
  if (get_keypone("Bayes_Debug_Flag", 0))
    set_keypair("Bayes_Fm1z_Matrix", 1, npar, 1, fm1z);

  /* Calculate the Posterior Variance-Covariance matrix */

  for (int i = 0; i < npar2; i++)
    post_S[i] = invS[i] + invH[i];
  if (matrix_invert(post_S, npar, -1)) goto label_end;
  if (get_keypone("Bayes_Debug_Flag", 0))
    set_keypair("Bayes_Post_S_Matrix", 1, npar, npar, post_S);

  /* Calculate the Posterior Mean vector */

  matrix_product(npar, npar, 1, invS, prior_mean, post_mean);
  for (int i = 0; i < npar; i++)
    fm1z[i] += post_mean[i];
  matrix_product(npar, npar, 1, post_S, fm1z, post_mean);
  if (get_keypone("Bayes_Debug_Flag", 0))
    set_keypair("Bayes_Post_Mean_Matrix", 1, npar, 1, post_mean);

  /* Optional printout */

  if (verbose)
  {
    print_matrix("Posterior Mean", 0, 1, lmlayers->nlayers, lmlayers->nbfl,
    NULL,
                 post_mean);
    print_matrix("Posterior Variance", 0, 1, lmlayers->npar, lmlayers->npar,
    NULL,
                 post_S);
  }

  /* Modify the Data vector */

  for (int iech = 0; iech < nech; iech++)
    for (int ipar = 0; ipar < npar; ipar++)
      zval[iech] -= FFTAB(ipar,iech) * post_mean[ipar];
  if (get_keypone("Bayes_Debug_Flag", 0))
    set_keypair("Mlayers_Zval2_Matrix", 1, lmlayers->neq, 1, zval);

  /* Auxiliary arrays prepared for estimation */

  matrix_product(nech, npar, npar, fft, post_S, a0);
  matrix_product(nech, nech, npar, acov, fft, ss);
  for (int i = 0; i < npar2; i++)
    invS[i] = post_S[i];
  if (matrix_invert(invS, npar, -1)) goto label_end;
  for (int i = 0; i < npar2; i++)
    gg[i] = invH[i] + invS[i];
  if (matrix_invert(gg, npar, -1)) goto label_end;
  if (get_keypone("Bayes_Debug_Flag", 0))
    set_keypair("Mlayers_GG_Matrix", 1, npar, npar, gg);

  matrix_prod_norme(1, nech, npar, ss, gg, cc);
  for (int i = 0; i < nech2; i++)
    cc[i] = acov[i] - cc[i];
  for (int i = 0; i < npar2; i++)
    gs[i] = invH[i] + invS[i];
  if (matrix_invert(gs, npar, -1)) goto label_end;
  if (get_keypone("Bayes_Debug_Flag", 0))
    set_keypair("Mlayers_CC_Matrix", 1, nech, nech, cc);

  /* Set the error return code */

  error = 0;

  label_end: ffc = (double*) mem_free((char* ) ffc);
  fft = (double*) mem_free((char* ) fft);
  fm1z = (double*) mem_free((char* ) fm1z);
  gg = (double*) mem_free((char* ) gg);
  invH = (double*) mem_free((char* ) invH);
  invS = (double*) mem_free((char* ) invS);
  return (error);
}

/****************************************************************************/
/*!
 **  Multi-layers architecture estimation
 **
 ** \return  Error return code
 **
 ** \param[in]  dbin       Input Db structure
 ** \param[in]  dbout      Output Db structure
 ** \param[in]  model      Model structure
 ** \param[in]  neigh      Neigh structrue
 ** \param[in]  flag_same  1 if input and output files coincide
 ** \param[in]  flag_z     1 if the output must be converted back into depth
 ** \param[in]  flag_vel   1 if work is performed in Velocity, 0 for Depth
 ** \param[in]  flag_cumul 1 if work is performed in Depth; 0 in Thickness
 ** \param[in]  flag_ext   1 if external drift must be used; 0 otherwise
 ** \param[in]  flag_std   1 if the estimation error must be calculated
 ** \param[in]  flag_bayes 1 if the Bayesian hypothesis is used on drift coeffs
 ** \param[in]  irf_rank   Rank of the Intrinsic Random Function (0 or 1)
 ** \param[in]  match_time 1 if external drift matches time; 0 otherwise
 ** \param[in]  dim_prior  Dimension of the prior information (for verification)
 ** \param[in]  prior_mean Vector of prior means for drift coefficients
 ** \param[in]  prior_vars Vector of prior variances for drift coefficients
 ** \param[in]  colrefd    Rank of the reference Depth variable in Dbout
 ** \param[in]  colreft    Rank of the reference Time variable in Dbout
 ** \param[in]  colrefb    Rank of the Bottom Depth variable in Dbout (or -1)
 ** \param[in]  verbose    Verbose option
 **
 *****************************************************************************/
int multilayers_kriging(Db *dbin,
                                        Db *dbout,
                                        Model *model,
                                        Neigh *neigh,
                                        int flag_same,
                                        int flag_z,
                                        int flag_vel,
                                        int flag_cumul,
                                        int flag_ext,
                                        int flag_std,
                                        int flag_bayes,
                                        int irf_rank,
                                        int match_time,
                                        int dim_prior,
                                        double *prior_mean,
                                        double *prior_vars,
                                        int colrefd,
                                        int colreft,
                                        int colrefb,
                                        int verbose)
{
  int *seltab, iptr, nlayers, ilayer, nechmax, nech, iech, neq, nvar, npar,
      error;
  double *a, *b, *b2, *baux, *zval, *dual, *covtab, *prop1, *prop2, *c00, *wgt;
  double *acov, *atot;
  double *fftab, *a0, *cc, *ss, *gs, *post_mean, *post_S;
  ELoc ptime;
  LMlayers *lmlayers;

  /* Preliminary checks */

  error = 1;

  iptr = -1;
  seltab = nullptr;
  covtab = nullptr;
  a = b = b2 = baux = zval = prop1 = prop2 = dual = nullptr;
  c00 = wgt = nullptr;
  acov = atot = nullptr;
  fftab = a0 = cc = ss = gs = post_mean = post_S = nullptr;
  lmlayers = nullptr;
  nlayers = model->getVariableNumber();
  nechmax = dbin->getSampleNumber();
  ptime = (match_time) ? ELoc::F :
                         ELoc::TIME;
  if (krige_koption_manage(1, 1, EKrigOpt::PONCTUAL, 1, VectorInt()))
    goto label_end;
  if (dbin->getNDim() != 2)
  {
    messerr("The input Db must be defined in 2-D");
    goto label_end;
  }
  if (dbout->getNDim() != 2)
  {
    messerr("The output Db must be defined in 2-D");
    goto label_end;
  }
  if (!dbin->isVariableNumberComparedTo(1)) goto label_end;
  if (!flag_same && !is_grid(dbout))
  {
    messerr("If Input and Output are different, Output should be a Grid Db");
    goto label_end;
  }
  if (!exist_LOCATOR(dbin, ELoc::LAYER))
  {
    messerr("The input Db must contain a LAYER locator");
    goto label_end;
  }
  if (flag_ext && nlayers != dbout->getExternalDriftNumber())
  {
    messerr("Inconsistency between:");
    messerr("- the number of variables in the Model (%d)", nlayers);
    messerr("- the number of external drifts in the Output Db File (%d)",
            dbout->getExternalDriftNumber());
    goto label_end;
  }
  if (flag_vel && nlayers != get_LOCATOR_NITEM(dbout, ptime))
  {
    messerr("Inconsistency between:");
    messerr("- the number of variables in the Model (%d)", nlayers);
    messerr("- the number of time variables in the Output Db File (%d)",
            get_LOCATOR_NITEM(dbout, ptime));
    goto label_end;
  }
  if (neigh->getType() != ENeigh::UNIQUE)
  {
    messerr("This procedure is only available in Unique Neighborhood");
    goto label_end;
  }
  if (flag_std && colrefb >= 0)
  {
    messerr("Calculation of the standard deviation of the estimation error");
    messerr("has not been programmed yet in collocation case");
    goto label_end;
  }
  if (flag_bayes && colrefb >= 0)
  {
    messerr("Use of Bayesian hypothesis has not been programmed yet");
    messerr("in collocation case");
    goto label_end;
  }
  if (flag_cumul && colrefb >= 0)
  {
    messerr("Collocation option is not coded when the results are expected");
    messerr("directly expressed in Depth (rather than Thickness)");
    goto label_end;
  }
  if (prior_mean == nullptr || prior_vars == nullptr) flag_bayes = 0;
  if (flag_bayes && dim_prior
      != st_get_number_drift(irf_rank, flag_ext) * nlayers)
  {
    messerr("The dimension of the Prior information (%d)", dim_prior);
    messerr("must be equal to %d (nlayers) x %d (nbfl)", nlayers,
            st_get_number_drift(irf_rank, flag_ext));
  }
  if (manage_external_info(1, ELoc::F, dbin, dbout, &iptr)) goto label_end;

  /* Allocating the output variables */

  nvar = nlayers;
  if (flag_std) nvar += nlayers;
  iptr = dbout->addFieldsByConstant(nvar, TEST, String(), ELoc::Z);

  /* Fill the Multi-Layers internal structure */

  lmlayers = lmlayers_alloc(flag_same, flag_vel, flag_cumul, flag_ext, flag_z,
                            colrefd, colreft, colrefb, irf_rank, match_time,
                            nlayers);

  /* Core allocation */

  seltab = (int*) mem_alloc(sizeof(int) * nechmax, 1);
  prop1 = (double*) mem_alloc(sizeof(double) * nlayers, 1);
  prop2 = (double*) mem_alloc(sizeof(double) * nlayers, 1);

  /* Calculate the number of active samples */

  for (iech = 0; iech < nechmax; iech++)
  {
    seltab[iech] = 0;
    ilayer = (int) get_LOCATOR_ITEM(dbin, ELoc::LAYER, 0, iech);
    if (ilayer < 1 || ilayer > nlayers) continue;
    if (st_get_props_data(lmlayers, dbin, dbout, iech, ilayer, prop1)) continue;
    seltab[iech] = 1;
  }

  /* Check the definition of all auxiliary variables defined on output file */
  /* Count the number of active samples (including the duplicates) */

  nech = st_check_auxiliary_variables(lmlayers, dbin, dbout, seltab);
  lmlayers->nech = nech;

  /* Allocation */

  npar = lmlayers->npar;
  neq = lmlayers->nech + npar;
  atot = (double*) mem_alloc(sizeof(double) * neq * neq, 1);
  acov = (double*) mem_alloc(sizeof(double) * nech * nech, 1);
  b = (double*) mem_alloc(sizeof(double) * neq, 1);
  b2 = (double*) mem_alloc(sizeof(double) * neq, 1);
  baux = (double*) mem_alloc(sizeof(double) * neq, 1);
  zval = (double*) mem_alloc(sizeof(double) * neq, 1);
  dual = (double*) mem_alloc(sizeof(double) * neq, 1);
  wgt = (double*) mem_alloc(sizeof(double) * neq, 1);
  covtab = (double*) mem_alloc(sizeof(double) * nlayers * nlayers, 1);
  c00 = (double*) mem_alloc(sizeof(double) * nlayers, 1);
  if (flag_bayes)
  {
    fftab = (double*) mem_alloc(sizeof(double) * nech * npar, 1);
    a0 = (double*) mem_alloc(sizeof(double) * nech * npar, 1);
    cc = (double*) mem_alloc(sizeof(double) * nech * nech, 1);
    ss = (double*) mem_alloc(sizeof(double) * nech * npar, 1);
    gs = (double*) mem_alloc(sizeof(double) * npar * npar, 1);
    post_S = (double*) mem_alloc(sizeof(double) * npar * npar, 1);
    post_mean = (double*) mem_alloc(sizeof(double) * npar, 1);
  }
  lmlayers->neq = neq;
  if (verbose) lmlayers_print(lmlayers);

  /* Establish the kriging matrix */

  st_lhs(lmlayers, dbin, dbout, model, seltab, prop1, prop2, covtab, atot,
         acov);
  if (is_debug_reference_defined() || debug_query("kriging"))
    krige_lhs_print(nech, neq, neq, NULL, atot);

  /* Establish the data vector */

  st_data_vector(lmlayers, dbin, dbout, seltab, zval);
  if (is_debug_reference_defined() || debug_query("kriging"))
  {
    mestitle(0, "Data Vector");
    message("Number of active samples  = %d\n", nech);
    message("Total number of equations = %d\n", neq);
    print_matrix("Data", 0, 1, 1, nech, NULL, zval);
  }

  /* Assign the Variance-Covariance matrix */

  a = (flag_bayes) ? acov :
                     atot;

  /* Calculate the Posterior in the Bayesian case */

  if (flag_bayes)
  {
    if (st_drift_data(lmlayers, dbin, dbout, seltab, prop1, fftab))
      goto label_end;
    if (st_drift_bayes(lmlayers, verbose, prior_mean, prior_vars, a, zval,
                       fftab, a0, cc, ss, gs, post_mean, post_S))
      goto label_end;
  }
  else
  {
    if (matrix_invert(a, neq, -1)) goto label_end;
    matrix_product(neq, neq, 1, a, zval, dual);
  }

  /* Perform the estimation over the grid nodes */

  st_estimate(lmlayers, dbin, dbout, model, seltab, flag_bayes, flag_std, a,
              zval, dual, prop1, prop2, covtab, b, b2, baux, wgt, c00, fftab,
              a0, cc, ss, gs, prior_mean, post_mean);

  /* Reconstitute the surfaces (optional) */

  st_convert_results(lmlayers, dbout, flag_std);

  /* Set the error return code */

  error = 0;

  label_end: (void) krige_koption_manage(-1, 1, EKrigOpt::PONCTUAL, 1,
                                         VectorInt());
  (void) manage_external_info(-1, ELoc::F, dbin, dbout, &iptr);
  seltab = (int*) mem_free((char* ) seltab);
  prop1 = (double*) mem_free((char* ) prop1);
  prop2 = (double*) mem_free((char* ) prop2);
  covtab = (double*) mem_free((char* ) covtab);
  zval = (double*) mem_free((char* ) zval);
  dual = (double*) mem_free((char* ) dual);
  atot = (double*) mem_free((char* ) atot);
  acov = (double*) mem_free((char* ) acov);
  b = (double*) mem_free((char* ) b);
  b2 = (double*) mem_free((char* ) b2);
  baux = (double*) mem_free((char* ) baux);
  c00 = (double*) mem_free((char* ) c00);
  wgt = (double*) mem_free((char* ) wgt);
  fftab = (double*) mem_free((char* ) fftab);
  a0 = (double*) mem_free((char* ) a0);
  cc = (double*) mem_free((char* ) cc);
  ss = (double*) mem_free((char* ) ss);
  gs = (double*) mem_free((char* ) gs);
  post_S = (double*) mem_free((char* ) post_S);
  post_mean = (double*) mem_free((char* ) post_mean);
  lmlayers = lmlayers_free(lmlayers);
  return (error);
}

/****************************************************************************/
/*!
 **  Evaluate the sill matrix for a given lag
 **
 ** \return  Error return code (proportions not calculatable)
 **
 ** \param[in]  lmlayers   LMlayers structure
 ** \param[in]  dbin       Input Db structure
 ** \param[in]  dbout      Output Db structure
 ** \param[in]  vorder     Vario_Order structure
 ** \param[in]  nlayers    Number of layers
 ** \param[in]  ifirst     Index of the first pair (included)
 ** \param[in]  ilast      Index of the last pair (excluded)
 ** \param[in]  zval       Array containing the sample values
 **
 ** \param[out] nval       Number of relevant pairs
 ** \param[out] distsum    Average distance
 ** \param[out] stat       Working array (Dimension: nlayers * nlayers)
 ** \param[out] phia       Working array for proportions (Dimension: nlayers)
 ** \param[out] phib       Working array for proportions (Dimension: nlayers)
 ** \param[out] atab       Working array (Dimension: nhalf * nhalf)
 ** \param[out] btab       Working array (Dimension: nhalf)
 **
 *****************************************************************************/
static int st_evaluate_lag(LMlayers *lmlayers,
                           Db *dbin,
                           Db *dbout,
                           Vario_Order *vorder,
                           int nlayers,
                           int ifirst,
                           int ilast,
                           double *zval,
                           int *nval,
                           double *distsum,
                           int *stat,
                           double *phia,
                           double *phib,
                           double *atab,
                           double *btab)
{
  int iech, jech, iiech, jjech, ilayer, jlayer, ecr1, ecr2, nhalf;
  double z1, z2, dist, fact1, fact2;

  /* Local initializations */

  (*nval) = 0;
  (*distsum) = 0.;
  nhalf = nlayers * (nlayers + 1) / 2;
  for (int i = 0; i < nhalf; i++)
    btab[i] = 0.;
  for (int i = 0; i < nhalf * nhalf; i++)
    atab[i] = 0.;
  for (int i = 0; i < nlayers * nlayers; i++)
    stat[i] = 0;

  /* Loop on the pairs contributing to the lag */

  for (int ipair = ifirst; ipair < ilast; ipair++)
  {
    vario_order_get_indices(vorder, ipair, &iiech, &jjech, &dist);
    vario_order_get_auxiliary(vorder, ipair, (char*) &iech, (char*) &jech);
    z1 = zval[iiech];
    z2 = zval[jjech];
    (*distsum) += dist;

    ilayer = (int) get_LOCATOR_ITEM(dbin, ELoc::LAYER, 0, iech);
    if (st_get_props_data(lmlayers, dbin, dbout, iech, ilayer, phia))
      return (1);

    jlayer = (int) get_LOCATOR_ITEM(dbin, ELoc::LAYER, 0, jech);
    if (st_get_props_data(lmlayers, dbin, dbout, jech, jlayer, phib))
      return (1);

    ecr1 = 0;
    stat[(ilayer - 1) * nlayers + (jlayer - 1)] += 1;

    for (int il1 = 0; il1 < nlayers; il1++)
      for (int jl1 = 0; jl1 <= il1; jl1++, ecr1++)
      {
        fact1 = phia[il1] * phib[jl1];
        if (il1 != jl1) fact1 += phia[jl1] * phib[il1];
        btab[ecr1] += fact1 * z1 * z2;

        ecr2 = 0;
        for (int il2 = 0; il2 < nlayers; il2++)
          for (int jl2 = 0; jl2 <= il2; jl2++, ecr2++)
          {
            fact2 = phia[il2] * phib[jl2];
            if (il2 != jl2) fact2 += phia[jl2] * phib[il2];
            ATAB(nhalf,ecr1,ecr2) += fact1 * fact2;
          }
      }

    (*nval)++;
  }
  (*distsum) /= (*nval);
  return (0);
}

/****************************************************************************/
/*!
 **  Calculate the variogram for each lag per direction
 **
 ** \return  Error return code
 **
 ** \param[in]  lmlayers   LMlayers structure
 ** \param[in]  verbose    1 for a verbose option
 ** \param[in]  dbin       Input Db structure
 ** \param[in]  dbout      Output Db structure
 ** \param[in]  vorder     Vario_Order structure
 ** \param[in]  zval       Data vector
 ** \param[in]  idir       Rank of the Direction
 **
 ** \param[out] vario      Vario structure
 **
 *****************************************************************************/
static int st_varioexp_chh(LMlayers *lmlayers,
                           int verbose,
                           Db *dbin,
                           Db *dbout,
                           Vario_Order *vorder,
                           double *zval,
                           int idir,
                           Vario *vario)
{
  double *atab, *btab, *phia, *phib, *sill, distsum;
  int *stat, error, nlayers, iadlag, nhalf, nhalf2, nval;
  int ipas, number, ifirst, ilast, ilayer, jlayer, ijl;

  /* Initializations */

  error = 1;
  sill = atab = btab = phia = phib = nullptr;
  stat = nullptr;
  nlayers = lmlayers->nlayers;
  nhalf = nlayers * (nlayers + 1) / 2;
  nhalf2 = nhalf * nhalf;

  /* Core allocation */

  phia = (double*) mem_alloc(sizeof(double) * nlayers, 0);
  if (phia == nullptr) goto label_end;
  phib = (double*) mem_alloc(sizeof(double) * nlayers, 0);
  if (phib == nullptr) goto label_end;
  btab = (double*) mem_alloc(sizeof(double) * nhalf, 0);
  if (btab == nullptr) goto label_end;
  atab = (double*) mem_alloc(sizeof(double) * nhalf2, 0);
  if (atab == nullptr) goto label_end;
  sill = (double*) mem_alloc(sizeof(double) * nhalf, 0);
  if (sill == nullptr) goto label_end;
  stat = (int*) mem_alloc(sizeof(int) * nlayers * nlayers, 0);
  if (stat == nullptr) goto label_end;

  /* Loop on the lags */

  for (ipas = 0; ipas < vario->getLagNumber(idir); ipas++)
  {
    vario_order_get_bounds(vorder, idir, ipas, &ifirst, &ilast);
    number = ilast - ifirst;
    if (number <= 0) continue;

    /* Loop on the pairs contributing to the lag */

    if (st_evaluate_lag(lmlayers, dbin, dbout, vorder, nlayers, ifirst, ilast,
                        zval, &nval, &distsum, stat, phia, phib, atab, btab))
      goto label_end;

    if (debug_query("variogram"))
    {
      message("Lag %d\n", ipas + 1);
      print_matrix("L.H.S.", 0, 1, nhalf, nhalf, NULL, atab);
      print_matrix("R.H.S.", 0, 1, 1, nhalf, NULL, btab);
    }

    if (matrix_invert(atab, nhalf, -2))
    {
      messerr("--> Inversion problem for lag %d", ipas + 1);
      if (verbose)
      {
        /* Matrix must be evaluated (as it has been destroyed by inversion) */
        (void) st_evaluate_lag(lmlayers, dbin, dbout, vorder, nlayers, ifirst,
                               ilast, zval, &nval, &distsum, stat, phia, phib,
                               atab, btab);
        messerr("Number of pairs  = %d", nval);
        messerr("Average distance = %lf", distsum);
        print_imatrix("Number of samples per layer", 0, 1, nlayers, nlayers,
        NULL,
                      stat);
        print_matrix("L.H.S.", 0, 1, nhalf, nhalf, NULL, atab);
        print_matrix("R.H.S.", 0, 1, 1, nhalf, NULL, btab);
      }
      continue;
    }
    matrix_product(1, nhalf, nhalf, btab, atab, sill);

    /* Optional printout */

    if (debug_query("variogram")) print_trimat("C(h)", 2, nlayers, sill);

    /* Store the covariance values */

    ijl = 0;
    for (ilayer = 0; ilayer < nlayers; ilayer++)
      for (jlayer = 0; jlayer <= ilayer; jlayer++, ijl++)
      {
        iadlag = vario->getDirAddress(idir, ilayer, jlayer, ipas, false, 1);
        vario->setGgByRank(idir, iadlag, sill[ijl]);
        vario->setHhByRank(idir, iadlag, distsum);
        vario->setSwByRank(idir, iadlag, nval);
        iadlag = vario->getDirAddress(idir, ilayer, jlayer, ipas, false, -1);
        vario->setGgByRank(idir, iadlag, sill[ijl]);
        vario->setHhByRank(idir, iadlag, -distsum);
        vario->setSwByRank(idir, iadlag, nval);
      }
  }

  /* Set the error return code */

  error = 0;

  label_end: phia = (double*) mem_free((char* ) phia);
  phib = (double*) mem_free((char* ) phib);
  atab = (double*) mem_free((char* ) atab);
  btab = (double*) mem_free((char* ) btab);
  sill = (double*) mem_free((char* ) sill);
  stat = (int*) mem_free((char* ) stat);
  return (error);
}

/****************************************************************************/
/*!
 **  Multi-layers architecture experimental variogram
 **
 ** \return  Error return code
 **
 ** \param[in]  dbin       Input Db structure
 ** \param[in]  dbout      Output Db structure
 ** \param[in]  vario      Vario structure
 ** \param[in]  nlayers    Number of layers
 ** \param[in]  flag_vel   1 if work is performed in Velocity, 0 for Depth
 ** \param[in]  flag_ext   1 if external drift must be used; 0 otherwise
 ** \param[in]  irf_rank   Rank of the Intrinsic Random Function (0 or 1)
 ** \param[in]  match_time 1 if external drift matches time; 0 otherwise
 ** \param[in]  colrefd    Rank of the reference Depth variable in Dbout
 ** \param[in]  colreft    Rank of the reference Time variable in Dbout
 ** \param[in]  verbose    1 for a  verbose option
 **
 *****************************************************************************/
int multilayers_vario(Db *dbin,
                                      Db *dbout,
                                      Vario *vario,
                                      int nlayers,
                                      int flag_vel,
                                      int flag_ext,
                                      int irf_rank,
                                      int match_time,
                                      int colrefd,
                                      int colreft,
                                      int verbose)
{
  int *seltab, error, ilayer, nechmax, nech, iech, idir, iptr;
  double *prop1, *zval;
  ELoc ptime;
  LMlayers *lmlayers;
  Vario_Order *vorder;

  /* Preliminary checks */

  error = 1;
  seltab = nullptr;
  prop1 = zval = nullptr;
  lmlayers = nullptr;
  vorder = nullptr;
  nechmax = dbin->getSampleNumber();
  ptime = (match_time) ? ELoc::F :
                         ELoc::TIME;
  if (dbin->getNDim() != 2)
  {
    messerr("The input Db must be defined in 2-D");
    goto label_end;
  }
  if (flag_vel && dbout->getNDim() != 2)
  {
    messerr("The output Db must be defined in 2-D");
    goto label_end;
  }
  if (!dbin->isVariableNumberComparedTo(1)) goto label_end;
  if (!exist_LOCATOR(dbin, ELoc::LAYER))
  {
    messerr("The input Db must contain a LAYER locator");
    goto label_end;
  }
  if (flag_ext && nlayers != dbout->getExternalDriftNumber())
  {
    messerr("Inconsistency between:");
    messerr("- the number of variables in the Model (%d)", nlayers);
    messerr("- the number of external drifts in the Output Db File (%d)",
            dbout->getExternalDriftNumber());
    goto label_end;
  }
  if (flag_vel && nlayers != get_LOCATOR_NITEM(dbout, ptime))
  {
    messerr("Inconsistency between:");
    messerr("- the number of layers (%d)", nlayers);
    messerr("- the number of time variables in the Output Db File (%d)",
            get_LOCATOR_NITEM(dbout, ptime));
    goto label_end;
  }
  if (manage_external_info(1, ELoc::F, dbin, dbout, &iptr)) goto label_end;

  /* Fill the Multi-Layers internal structure */

  lmlayers = lmlayers_alloc(0, flag_vel, 0, flag_ext, 0, colrefd, colreft, -1,
                            irf_rank, match_time, nlayers);
  if (verbose) lmlayers_print(lmlayers);

  /* Core allocation */

  seltab = (int*) mem_alloc(sizeof(int) * nechmax, 1);
  prop1 = (double*) mem_alloc(sizeof(double) * nlayers, 1);
  zval = (double*) mem_alloc(sizeof(double) * nechmax, 1);

  /* Calculate the number of active samples */

  for (iech = 0; iech < nechmax; iech++)
  {
    seltab[iech] = 0;
    ilayer = (int) get_LOCATOR_ITEM(dbin, ELoc::LAYER, 0, iech);
    if (ilayer < 1 || ilayer > nlayers) continue;
    if (st_get_props_data(lmlayers, dbin, dbout, iech, ilayer, prop1)) continue;
    seltab[iech] = 1;
  }

  /* Check the definition of all auxiliary variables defined on output file */
  /* Count the number of active samples (including the duplicates) */

  nech = st_check_auxiliary_variables(lmlayers, dbin, dbout, seltab);
  lmlayers->nech = nech;
  lmlayers->neq = nech + lmlayers->npar;

  /* Establish the data vector */

  st_data_vector(lmlayers, dbin, dbout, seltab, zval);

  /* Subtract the optimal average or drift */

  if (st_subtract_optimal_drift(lmlayers, verbose, dbin, dbout, seltab, zval))
    goto label_end;

  /* Evaluate the Geometry */

  vorder = vario_order_manage(1, 1, sizeof(int), NULL);
  if (variogram_mlayers(dbin, seltab, vario, vorder)) goto label_end;

  /* Evaluate the variogram */

  for (idir = 0; idir < vario->getDirectionNumber(); idir++)
  {
    if (st_varioexp_chh(lmlayers, verbose, dbin, dbout, vorder, zval, idir,
                        vario)) goto label_end;
  }

  /* Set the error return code */

  error = 0;

  label_end: (void) manage_external_info(-1, ELoc::F, dbin, dbout, &iptr);
  seltab = (int*) mem_free((char* ) seltab);
  prop1 = (double*) mem_free((char* ) prop1);
  zval = (double*) mem_free((char* ) zval);
  vorder = vario_order_manage(-1, 1, sizeof(int), vorder);
  lmlayers = lmlayers_free(lmlayers);
  return (error);
}

/****************************************************************************/
/*!
 **  Determine the mean and variance of drift coefficients
 **
 ** \return  Error return code
 **
 ** \param[in]  nech       Number of samples
 ** \param[in]  npar       Number of drift coefficients
 ** \param[in]  zval       The data vector (Dimension: neq)
 ** \param[in]  fftab      Drift array (Dimension: npar[nrow] * nech[ncol])
 **
 ** \param[out] mean       Array of means
 ** \param[out] vars       Array of variances
 **
 *****************************************************************************/
static int st_get_prior(int nech,
                        int npar,
                        double *zval,
                        double *fftab,
                        double *mean,
                        double *vars)
{
  double *atab, *btab, *atab0, *btab0, *result;
  int error, pivot, size, ecr;

  /* Initializations */

  error = 1;
  atab = btab = atab0 = btab0 = result = nullptr;
  size = npar * (npar + 1) / 2;

  /* Core allocation */

  atab = (double*) mem_alloc(sizeof(double) * size, 0);
  if (atab == nullptr) goto label_end;
  atab0 = (double*) mem_alloc(sizeof(double) * size, 0);
  if (atab0 == nullptr) goto label_end;
  btab = (double*) mem_alloc(sizeof(double) * npar, 0);
  if (btab == nullptr) goto label_end;
  btab0 = (double*) mem_alloc(sizeof(double) * npar, 0);
  if (btab0 == nullptr) goto label_end;
  result = (double*) mem_alloc(sizeof(double) * npar, 0);
  if (result == nullptr) goto label_end;
  for (int i = 0; i < npar; i++)
    btab0[i] = 0.;
  for (int i = 0; i < size; i++)
    atab0[i] = 0.;
  for (int i = 0; i < npar; i++)
    mean[i] = 0.;
  for (int i = 0; i < npar * npar; i++)
    vars[i] = 0.;

  /* Loop on the data */

  for (int iech = 0; iech < nech; iech++)
    for (int ipar = ecr = 0; ipar < npar; ipar++)
    {
      btab0[ipar] += zval[iech] * FFTAB(ipar, iech);
      for (int jpar = 0; jpar <= ipar; jpar++, ecr++)
        atab0[ecr] += FFTAB(ipar,iech) * FFTAB(jpar, iech);
    }

  /* Optional printout */

  if (get_keypone("Bayes_Debug_Flag", 0))
  {
    set_keypair("Bayes_Get_Prior_ATAB0", 1, size, 1, atab0);
    set_keypair("Bayes_Get_Prior_BTAB0", 1, npar, 1, btab0);
  }

  /* Bootstrap for the variance-covariance */

  for (int iech = 0; iech < nech; iech++)
  {
    for (int i = 0; i < npar; i++)
      btab[i] = btab0[i];
    for (int i = 0; i < size; i++)
      atab[i] = atab0[i];

    /* Update the arrays by suppressing the current data */

    for (int ipar = ecr = 0; ipar < npar; ipar++)
    {
      btab[ipar] -= zval[iech] * FFTAB(ipar, iech);
      for (int jpar = 0; jpar <= ipar; jpar++, ecr++)
        atab[ecr] -= FFTAB(ipar,iech) * FFTAB(jpar, iech);
    }

    /* Solve the system */

    if (matrix_solve(0, atab, btab, result, npar, 1, &pivot)) goto label_end;

    /* Update the statistics */

    for (int ipar = 0; ipar < npar; ipar++)
    {
      mean[ipar] += result[ipar];
      for (int jpar = 0; jpar < npar; jpar++)
        VARS(npar,ipar,jpar) += result[ipar] * result[jpar];
    }
  }

  /* Normalize the results */

  for (int ipar = 0; ipar < npar; ipar++)
    mean[ipar] /= nech;
  for (int ipar = 0; ipar < npar; ipar++)
    for (int jpar = 0; jpar < npar; jpar++)
      VARS(npar,ipar,jpar) = VARS(npar,ipar,jpar) / nech
          - mean[ipar] * mean[jpar];

  /* Set the error return code */

  error = 0;

  label_end: atab = (double*) mem_free((char* ) atab);
  atab0 = (double*) mem_free((char* ) atab0);
  btab = (double*) mem_free((char* ) btab);
  btab0 = (double*) mem_free((char* ) btab0);
  result = (double*) mem_free((char* ) result);
  return (error);
}

/****************************************************************************/
/*!
 **  Multi-layers get the mean and prior matrices for Bayesian prior
 **
 ** \return  Error return code
 **
 ** \param[in]  dbin       Input Db structure
 ** \param[in]  dbout      Output Db structure
 ** \param[in]  model      Model structure
 ** \param[in]  flag_same  1 if input and output files coincide
 ** \param[in]  flag_vel   1 if work is performed in Velocity, 0 for Depth
 ** \param[in]  flag_ext   1 if external drift must be used; 0 otherwise
 ** \param[in]  irf_rank   Rank of the Intrinsic Random Function (0 or 1)
 ** \param[in]  match_time 1 if external drift matches time; 0 otherwise
 ** \param[in]  colrefd    Rank of the reference Depth variable in Dbout
 ** \param[in]  colreft    Rank of the reference Time variable in Dbout
 ** \param[in]  colrefb    Rank of the Bottom Depth variable in Dbout (or -1)
 ** \param[in]  verbose    Verbose option
 **
 ** \param[out] npar_arg   Number of drift terms
 ** \param[out] mean       Array of means
 ** \param[out] vars       Array of variances
 **
 *****************************************************************************/
int multilayers_get_prior(Db *dbin,
                                          Db *dbout,
                                          Model *model,
                                          int flag_same,
                                          int flag_vel,
                                          int flag_ext,
                                          int irf_rank,
                                          int match_time,
                                          int colrefd,
                                          int colreft,
                                          int colrefb,
                                          int verbose,
                                          int *npar_arg,
                                          double **mean,
                                          double **vars)
{
  int nlayers, ilayer, nechmax, nech, iech, npar, error, iptr, neq;
  int *seltab;
  double *zval, *props, *fftab;
  ELoc ptime;
  LMlayers *lmlayers;

  /* Preliminary checks */

  error = 1;
  iptr = -1;
  seltab = nullptr;
  fftab = zval = props = nullptr;
  lmlayers = nullptr;
  nlayers = model->getVariableNumber();
  nechmax = dbin->getSampleNumber();
  ptime = (match_time) ? ELoc::F :
                         ELoc::TIME;
  if (krige_koption_manage(1, 1, EKrigOpt::PONCTUAL, 1, VectorInt()))
    goto label_end;
  if (dbin->getNDim() != 2)
  {
    messerr("The input Db must be defined in 2-D");
    goto label_end;
  }
  if (dbout->getNDim() != 2)
  {
    messerr("The output Db must be defined in 2-D");
    goto label_end;
  }
  if (!dbin->isVariableNumberComparedTo(1)) goto label_end;
  if (!flag_same && !is_grid(dbout))
  {
    messerr("If Input and Output are different, Output should be a Grid Db");
    goto label_end;
  }
  if (!exist_LOCATOR(dbin, ELoc::LAYER))
  {
    messerr("The input Db must contain a LAYER locator");
    goto label_end;
  }
  if (flag_ext && nlayers != dbout->getExternalDriftNumber())
  {
    messerr("Inconsistency between:");
    messerr("- the number of variables in the Model (%d)", nlayers);
    messerr("- the number of external drifts in the Output Db File (%d)",
            dbout->getExternalDriftNumber());
    goto label_end;
  }
  if (flag_vel && nlayers != get_LOCATOR_NITEM(dbout, ptime))
  {
    messerr("Inconsistency between:");
    messerr("- the number of variables in the Model (%d)", nlayers);
    messerr("- the number of time variables in the Output Db File (%d)",
            get_LOCATOR_NITEM(dbout, ptime));
    goto label_end;
  }
  if (manage_external_info(1, ELoc::F, dbin, dbout, &iptr)) goto label_end;

  /* Fill the Multi-Layers internal structure */

  lmlayers = lmlayers_alloc(flag_same, flag_vel, 0, flag_ext, 1, colrefd,
                            colreft, colrefb, irf_rank, match_time, nlayers);

  /* Core allocation */

  seltab = (int*) mem_alloc(sizeof(int) * nechmax, 1);
  props = (double*) mem_alloc(sizeof(double) * nlayers, 1);

  /* Calculate the number of active samples */

  for (iech = 0; iech < nechmax; iech++)
  {
    seltab[iech] = 0;
    ilayer = (int) get_LOCATOR_ITEM(dbin, ELoc::LAYER, 0, iech);
    if (ilayer < 1 || ilayer > nlayers) continue;
    if (st_get_props_data(lmlayers, dbin, dbout, iech, ilayer, props)) continue;
    seltab[iech] = 1;
  }

  /* Check the definition of all auxiliary variables defined on output file */
  /* Count the number of active samples (including the duplicates) */

  nech = st_check_auxiliary_variables(lmlayers, dbin, dbout, seltab);
  lmlayers->nech = nech;
  lmlayers->neq = nech + lmlayers->npar;
  if (verbose) lmlayers_print(lmlayers);

  /* Allocation */

  npar = lmlayers->npar;
  neq = lmlayers->nech + npar;
  zval = (double*) mem_alloc(sizeof(double) * neq, 1);
  fftab = (double*) mem_alloc(sizeof(double) * nech * npar, 1);
  *mean = (double*) mem_alloc(sizeof(double) * npar, 1);
  *vars = (double*) mem_alloc(sizeof(double) * npar * npar, 1);

  /* Establish the data vector */

  st_data_vector(lmlayers, dbin, dbout, seltab, zval);

  /* Establish the drift matrix */

  if (st_drift_data(lmlayers, dbin, dbout, seltab, props, fftab))
    goto label_end;

  /* Estimate the optimal drift matrices */

  if (st_get_prior(nech, npar, zval, fftab, *mean, *vars)) goto label_end;

  /* Set the error return code */

  *npar_arg = npar;
  error = 0;

  label_end: (void) krige_koption_manage(-1, 1, EKrigOpt::PONCTUAL, 1,
                                         VectorInt());
  (void) manage_external_info(-1, ELoc::F, dbin, dbout, &iptr);
  seltab = (int*) mem_free((char* ) seltab);
  props = (double*) mem_free((char* ) props);
  fftab = (double*) mem_free((char* ) fftab);
  zval = (double*) mem_free((char* ) zval);
  lmlayers = lmlayers_free(lmlayers);
  if (error)
  {
    *mean = (double*) mem_free((char* ) *mean);
    *vars = (double*) mem_free((char* ) *vars);
  }
  return (error);
}
