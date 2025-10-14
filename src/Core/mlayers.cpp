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
#include "Basic/OptDbg.hpp"
#include "Basic/Utilities.hpp"
#include "Calculators/CalcMigrate.hpp"
#include "Core/Keypair.hpp"
#include "Db/Db.hpp"
#include "Db/DbGrid.hpp"
#include "Model/Model.hpp"
#include "Neigh/ANeigh.hpp"
#include "geoslib_f.h"
#include "geoslib_old_f.h"
#include <cmath>

/*! \cond */
#define IAD(n, i, j)       ((n) * (i) + (j))
#define A(i, j)            (a[IAD(neq, i, j)])
#define ACOV(i, j)         (acov[IAD(nech, i, j)])
#define GS(i, j)           (gs[IAD(npar, i, j)])
#define PHIA(i, ilayer)    (phia[IAD(nlayers, i, ilayer)])
#define PHIB(i, ilayer)    (phib[IAD(nlayers, i, ilayer)])
#define AA(i, ilayer2)     (aa[IAD(nlayer2, i, ilayer2)])
#define A2(n, i, j)        (a2[IAD(n, i, j)])
#define B2(n, i, j)        (b2[IAD(n, i, j)])
#define INVS(npar, i, j)   (invS[IAD(npar, i, j)])
#define FFTAB(ipar, iech)  (fftab[(iech) * npar + (ipar)])
#define POST_S(npar, i, j) (post_S[IAD(npar, i, j)])
#define ATAB(n, i, j)      (atab[IAD(n, i, j)])
#define VARS(n, i, j)      (vars[IAD(n, i, j)])
#define IAD(n, i, j)       ((n) * (i) + (j))
#define A(i, j)            (a[IAD(neq, i, j)])
#define GS(i, j)           (gs[IAD(npar, i, j)])
#define PHIA(i, ilayer)    (phia[IAD(nlayers, i, ilayer)])
#define PHIB(i, ilayer)    (phib[IAD(nlayers, i, ilayer)])
#define AA(i, ilayer2)     (aa[IAD(nlayer2, i, ilayer2)])
#define A2(n, i, j)        (a2[IAD(n, i, j)])
#define B2(n, i, j)        (b2[IAD(n, i, j)])
#define INVS(npar, i, j)   (invS[IAD(npar, i, j)])
#define FFTAB(ipar, iech)  (fftab[(iech) * npar + (ipar)])
#define POST_S(npar, i, j) (post_S[IAD(npar, i, j)])
#define ATAB(n, i, j)      (atab[IAD(n, i, j)])
#define VARS(n, i, j)      (vars[IAD(n, i, j)])
/*! \endcond */

namespace gstlrn
{
typedef struct
{
  bool flag_same;  /* True if input and output files coincide */
  bool flag_vel;   /* True for velocity; False for thickness */
  bool flag_cumul; /* True for cumulating in depth */
  bool flag_ext;   /* Use external drift */
  bool flag_z;     /* True if output must be converted into depth */
  Id colrefd;      /* Reference depth map (if >= 0) */
  Id colreft;      /* Reference time map (if >= 0) */
  Id colrefb;      /* Bottom map (if >= 0) */
  Id match_time;   /* True if Time provided through External Drift */
  ELoc ptime;      /* Pointer to the Time variables */
  Id nlayers;      /* Number of layers */
  Id nbfl;         /* Number of drift functions */
  Id nech;         /* Number of active samples */
  Id neq;          /* Number of equations */
  Id npar;         /* Nb. of layers * Nb. of drift functions */
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
static LMlayers* lmlayers_free(LMlayers* lmlayers)
{
  if (lmlayers == nullptr) return (lmlayers);
  delete lmlayers;
  lmlayers = nullptr;
  return (lmlayers);
}

/****************************************************************************/
/*!
 **  Returns the number of drift functions
 **
 ** \return  Number of drift conditions
 **
 ** \param[in]  irf_rank  Rank of the Intrinsic Random Function (0 or 1)
 ** \param[in]  flag_ext  True if external drift must be used; False otherwise
 **
 *****************************************************************************/
static Id st_get_number_drift(Id irf_rank, bool flag_ext)
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

static LMlayers* lmlayers_alloc(const Db* dbout,
                                bool flag_same,
                                bool flag_vel,
                                bool flag_cumul,
                                bool flag_ext,
                                bool flag_z,
                                const String& namerefd,
                                const String& namereft,
                                const String& namerefb,
                                Id irf_rank,
                                bool match_time,
                                Id nlayers)
{
  Id colrefd           = (namerefd.empty()) ? -1 : dbout->getColIdx(namerefd);
  Id colreft           = (namereft.empty()) ? -1 : dbout->getColIdx(namereft);
  Id colrefb           = (namerefb.empty()) ? -1 : dbout->getColIdx(namerefb);
  auto* lmlayers       = new LMlayers;
  lmlayers->flag_same  = flag_same;
  lmlayers->flag_vel   = flag_vel;
  lmlayers->flag_cumul = flag_cumul;
  lmlayers->flag_ext   = flag_ext;
  lmlayers->flag_z     = flag_z;
  lmlayers->colrefd    = colrefd;
  lmlayers->colreft    = colreft;
  lmlayers->colrefb    = colrefb;
  lmlayers->match_time = match_time;
  lmlayers->ptime      = (match_time) ? ELoc::F : ELoc::TIME;
  lmlayers->nlayers    = nlayers;
  lmlayers->nbfl       = st_get_number_drift(irf_rank, flag_ext);
  lmlayers->nech       = 0;
  lmlayers->neq        = 0;
  lmlayers->npar       = lmlayers->nbfl * nlayers;
  return (lmlayers);
}

/****************************************************************************/
/*!
 **  Print the Multi-Layers internal structure
 **
 ** \param[in]  lmlayers  Pointer to the LMlayers structure
 **
 *****************************************************************************/
static void lmlayers_print(LMlayers* lmlayers)

{
  static const char* NOK[] = {"NO", "YES"};

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
}

/****************************************************************************/
/*!
 **  Returns the absolute grid node absolute index which is the closest to a
 **  given sample of a Db
 **  In the case of same input and output file, simply return 'iech'
 **
 ** \return The rank or ITEST
 **
 ** \param[in]  lmlayers  LMlayers structure
 ** \param[in]  dbin      Input Db structure
 ** \param[in]  dbout     Output Db structure
 ** \param[in]  iech      Rank in the input Db
 **
 *****************************************************************************/
static Id st_locate_sample_in_output(LMlayers* lmlayers,
                                     Db* dbin,
                                     DbGrid* dbout,
                                     Id iech)
{
  /* In the case the input and output files coincide, simply return 'iech' */
  if (lmlayers->flag_same) return iech;

  Id ndim = dbout->getNDim();
  VectorInt indg(ndim, 0);
  VectorDouble coor(ndim);

  /* The files are different */
  for (Id idim = 0; idim < dbin->getNDim(); idim++)
    coor[idim] = dbin->getCoordinate(iech, idim);
  if (point_to_grid(dbout, coor.data(), 0, indg.data()) != 0) return ITEST;
  return dbout->indiceToRank(indg);
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
static void st_check_layer(const char* string, LMlayers* lmlayers, Id ilayer0)
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
static Id st_get_props_result(LMlayers* lmlayers,
                              Db* dbout,
                              Id iech,
                              Id ilayer0,
                              VectorDouble& props)
{
  double pval, t0, t1, tlast, tt;
  Id ilayer;

  /* Initializations */

  st_check_layer("st_get_props_result", lmlayers, ilayer0);
  for (ilayer = 0; ilayer < lmlayers->nlayers; ilayer++)
    props[ilayer] = 0.;

  /* Dispatch */

  if (lmlayers->flag_vel)
  {

    /* Working in velocities */

    t0 = (lmlayers->colreft >= 0) ? dbout->getArray(iech, lmlayers->colreft) : 0.;
    if (FFFF(t0)) return (1);
    t1 = dbout->getFromLocator(lmlayers->ptime, iech, ilayer0 - 1);
    if (FFFF(t1)) return (1);
    tlast = t0;

    /* Loop on the layers */

    for (ilayer = 0; ilayer < ilayer0; ilayer++)
    {
      tt = dbout->getFromLocator(lmlayers->ptime, iech, ilayer);
      if (FFFF(tt)) return (1);
      pval = (tt - tlast) / (t1 - t0);
      if (pval < 0 || pval > 1) return (1);
      tlast         = tt;
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
static Id st_get_props_data(LMlayers* lmlayers,
                            Db* dbin,
                            DbGrid* dbout,
                            Id iech,
                            Id ilayer0,
                            VectorDouble& props)
{
  Id igrid, ilayer;

  /* Initializations */

  st_check_layer("st_get_props_data", lmlayers, ilayer0);
  for (ilayer = 0; ilayer < lmlayers->nlayers; ilayer++)
    props[ilayer] = 0.;

  /* Get the sample rank in the output Db of the sample from the input Db */

  igrid = st_locate_sample_in_output(lmlayers, dbin, dbout, iech);
  if (IFFFF(igrid)) return (1);

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
static double st_get_drift_result(LMlayers* lmlayers,
                                  Db* dbout,
                                  Id iech,
                                  Id ilayer0)
{
  double drift;

  if (!lmlayers->flag_ext) return (TEST);
  st_check_layer("st_get_drift_result", lmlayers, ilayer0);

  drift = dbout->getLocVariable(ELoc::F, iech, ilayer0 - 1);
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
static double st_get_drift_data(LMlayers* lmlayers,
                                Db* dbin,
                                DbGrid* dbout,
                                Id iech,
                                Id ilayer0)
{
  Id igrid;
  double drift;

  if (!lmlayers->flag_ext) return (TEST);
  st_check_layer("st_get_drift_data", lmlayers, ilayer0);

  /* Get the sample rank in the output Db of the sample from the input Db*/

  igrid = st_locate_sample_in_output(lmlayers, dbin, dbout, iech);
  if (IFFFF(igrid)) return (TEST);

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
static void st_covariance_c00(LMlayers* lmlayers,
                              Model* model,
                              const VectorDouble& prop1,
                              MatrixSquare& covtab,
                              VectorDouble& c00)
{
  Id nlayers, flag_interrupt;
  double value;

  nlayers = lmlayers->nlayers;
  model->evaluateMatInPlace(nullptr, VectorDouble(), covtab, true);

  if (lmlayers->flag_cumul)
  {
    for (Id k = 0; k < nlayers; k++)
    {
      value          = 0.;
      flag_interrupt = 0;
      for (Id i = 0; i <= k && flag_interrupt == 0; i++)
        for (Id j = 0; j <= k && flag_interrupt == 0; j++)
        {
          if (FFFF(prop1[i]) || FFFF(prop1[j]))
            flag_interrupt = 1;
          else
            value += prop1[i] * prop1[j] * covtab.getValue(i, j);
        }
      c00[k] = (flag_interrupt) ? TEST : value;
    }
  }
  else
  {
    for (Id k = 0; k < nlayers; k++)
      c00[k] = covtab.getValue(k, k);
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
static double st_cij(LMlayers* lmlayers,
                     Model* model,
                     Id ilayer,
                     const VectorDouble& prop1,
                     Id jlayer,
                     const VectorDouble& prop2,
                     const VectorDouble& dd,
                     MatrixSquare& covtab)
{
  VectorDouble d1(2);
  st_check_layer("st_cij", lmlayers, ilayer);
  st_check_layer("st_cij", lmlayers, jlayer);

  /* Calculate the covariance matrix */

  d1[0] = (!dd.empty()) ? dd[0] : 0.;
  d1[1] = (!dd.empty()) ? dd[1] : 0.;
  model->evaluateMatInPlace(nullptr, d1, covtab, true);

  /* Evaluate the covariance term */

  double value = 0.;
  for (Id i = 0; i < ilayer; i++)
    for (Id j = 0; j < jlayer; j++)
    {
      if (FFFF(prop1[i])) return (TEST);
      if (FFFF(prop2[j])) return (TEST);
      value += prop1[i] * prop2[j] * covtab.getValue(i, j);
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
static double st_ci0(LMlayers* lmlayers,
                     Model* model,
                     Id ilayer,
                     const VectorDouble& prop1,
                     Id jlayer,
                     const VectorDouble& dd,
                     MatrixSquare& covtab)
{
  VectorDouble d1(2);
  st_check_layer("st_ci0", lmlayers, ilayer);
  st_check_layer("st_ci0", lmlayers, jlayer);

  /* Calculate the covariance matrix */

  d1[0] = (!dd.empty()) ? dd[0] : 0.;
  d1[1] = (!dd.empty()) ? dd[1] : 0.;
  model->evaluateMatInPlace(nullptr, d1, covtab, true);

  /* Evaluate the covariance term */

  double value = 0.;
  for (Id i = 0; i < ilayer; i++)
  {
    if (FFFF(prop1[i])) return (1);
    value += prop1[i] * covtab.getValue(i, jlayer - 1);
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
 ** \param[in]  propval     Value for the proportion (used if flag_cumul)
 ** \param[in]  drext       Value of the external drift
 ** \param[in,out] ipos_loc Address for the first drift term.
 **                         On output, address for the next term after the drift
 **
 ** \param[out] b         Array for storing the drift
 **
 *****************************************************************************/
static Id st_drift(LMlayers* lmlayers,
                   const VectorDouble& coor,
                   double propval,
                   double drext,
                   Id* ipos_loc,
                   VectorDouble& b)
{
  Id ipos;

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
static Id st_lhs_one(LMlayers* lmlayers,
                     Db* dbin,
                     DbGrid* dbout,
                     Model* model,
                     VectorInt& seltab,
                     Id iech0,
                     Id ilayer0,
                     VectorDouble& coor,
                     VectorDouble& prop0,
                     VectorDouble& prop2,
                     MatrixSquare& covtab,
                     VectorDouble& b)
{
  Id jech, jjech, jfois, jlayer, nlayers, i;
  double drext, coor2[2];
  VectorDouble d1(2);

  /* Initializations */

  nlayers = lmlayers->nlayers;
  b.fill(0);

  /* Covariance part */

  for (jech = jjech = 0; jech < dbin->getNSample(); jech++)
  {
    if (seltab[jech] == 0) continue;
    coor2[0] = dbin->getCoordinate(jech, 0);
    coor2[1] = dbin->getCoordinate(jech, 1);
    for (jfois = 0; jfois < seltab[jech]; jfois++, jjech++)
    {
      jlayer = (jfois == 0) ? static_cast<Id>(dbin->getFromLocator(ELoc::LAYER, jech)) : nlayers;

      /* Evaluate the proportion vector */

      if (st_get_props_data(lmlayers, dbin, dbout, jech, jlayer, prop2)) return (1);

      /* Calculate the distance vector */

      d1[0]    = coor[0] - coor2[0];
      d1[1]    = coor[1] - coor2[1];
      b[jjech] = st_cij(lmlayers, model, ilayer0, prop0, jlayer, prop2, d1, covtab);
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
static Id st_rhs(LMlayers* lmlayers,
                 Db* dbin,
                 DbGrid* dbout,
                 Model* model,
                 VectorDouble& coor,
                 VectorInt& seltab,
                 Id iechout,
                 Id ilayer0,
                 VectorDouble& prop0,
                 VectorDouble& prop2,
                 MatrixSquare& covtab,
                 VectorDouble& b)
{
  Id jech, jjech, i, jlayer, ipos, ifois, nlayers, ideb;
  double drext, coor2[2], propval;
  VectorDouble d1(2);

  /* Get the coordinates of the target */

  nlayers = lmlayers->nlayers;
  st_check_layer("st_rhs", lmlayers, ilayer0);
  coor[0] = dbout->getCoordinate(iechout, 0);
  coor[1] = dbout->getCoordinate(iechout, 1);

  /* Initialize the vector with zeroes */

  for (i = 0; i < lmlayers->neq; i++)
    b[i] = 0.;

  /* Covariance part */

  for (jech = jjech = 0; jech < dbin->getNSample(); jech++)
  {
    if (seltab[jech] == 0) continue;
    coor2[0] = dbin->getCoordinate(jech, 0);
    coor2[1] = dbin->getCoordinate(jech, 1);
    for (ifois = 0; ifois < seltab[jech]; ifois++, jjech++)
    {
      jlayer = (ifois == 0) ? static_cast<Id>(dbin->getFromLocator(ELoc::LAYER, jech)) : nlayers;

      /* Evaluate the proportion vector */

      (void)st_get_props_data(lmlayers, dbin, dbout, jech, jlayer, prop2);

      /* Calculate the distance vector */

      d1[0] = coor2[0] - coor[0];
      d1[1] = coor2[1] - coor[1];
      if (lmlayers->flag_cumul)
        b[jjech] = st_cij(lmlayers, model, ilayer0, prop0, jlayer, prop2, d1, covtab);
      else
        b[jjech] = st_ci0(lmlayers, model, jlayer, prop2, ilayer0, d1, covtab);
      if (FFFF(b[jjech])) return (1);
    }
  }

  /* Drift part */

  ideb = (lmlayers->flag_cumul) ? 0 : ilayer0 - 1;
  for (i = ideb; i < ilayer0; i++)
  {
    ipos    = lmlayers->nech + lmlayers->nbfl * i;
    drext   = st_get_drift_result(lmlayers, dbout, iechout, i + 1);
    propval = (lmlayers->flag_cumul) ? prop0[i] : 1.;
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
static Id st_lhs(LMlayers* lmlayers,
                 Db* dbin,
                 DbGrid* dbout,
                 Model* model,
                 VectorInt& seltab,
                 VectorDouble& prop1,
                 VectorDouble& prop2,
                 MatrixSquare& covtab,
                 MatrixSquare& a,
                 MatrixSquare& acov)
{
  Id iiech, jjech;
  VectorDouble coor(2);

  /* Initialize the matrix with zeroes */

  Id nech    = lmlayers->nech;
  Id neq     = lmlayers->neq;
  Id nlayers = lmlayers->nlayers;
  a.fill(0);
  VectorDouble b(neq);

  /* Loop on the first sample */

  iiech = 0;
  for (Id iech = 0; iech < dbin->getNSample(); iech++)
  {
    if (seltab[iech] == 0) continue;
    coor[0] = dbin->getCoordinate(iech, 0);
    coor[1] = dbin->getCoordinate(iech, 1);
    for (Id ifois = 0; ifois < seltab[iech]; ifois++, iiech++)
    {
      Id ilayer = (ifois == 0) ? static_cast<Id>(dbin->getFromLocator(ELoc::LAYER, iech)) : nlayers;

      /* Evaluate the proportion vector */

      if (st_get_props_data(lmlayers, dbin, dbout, iech, ilayer, prop1))
        return (1);

      /* Loop on the second sample */

      if (st_lhs_one(lmlayers, dbin, dbout, model, seltab, iech, ilayer, coor,
                     prop1, prop2, covtab, b)) return (1);

      for (Id i = 0; i < neq; i++)
        a.setValue(iiech, i, b[i]);
    }
  }

  /* Symmetrization */
  for (iiech = 0; iiech < neq; iiech++)
    for (jjech = 0; jjech <= iiech; jjech++)
      a.setValue(iiech, jjech, a.getValue(jjech, iiech));

  /* Extraction of the Covariance matrix */
  for (Id iech = 0; iech < nech; iech++)
    for (Id jech = 0; jech < nech; jech++)
      acov.setValue(iech, jech, a.getValue(iech, jech));

  if (get_keypone("Bayes_Debug_Flag", 0))
    set_keypair("Mlayers_LHS_Matrix", 1, neq, neq, a.getValues().data());

  return (0);
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
static void st_data_vector(LMlayers* lmlayers,
                           Db* dbin,
                           DbGrid* dbout,
                           VectorInt& seltab,
                           VectorDouble& zval)
{
  double value, dtime;
  Id iech, iiech, igrid, ifois, ilayer, nlayers;

  /* Initialize the vector with zeroes */

  igrid   = 0;
  nlayers = lmlayers->nlayers;
  zval.fill(0.);

  /* Loop on the samples */

  for (iech = iiech = 0; iech < dbin->getNSample(); iech++)
  {
    if (seltab[iech] == 0) continue;

    /* Calculate the grid node index (optional) */

    if (lmlayers->colrefd >= 0 || lmlayers->colreft >= 0 ||
        lmlayers->colrefb >= 0 || lmlayers->flag_vel)
      igrid = st_locate_sample_in_output(lmlayers, dbin, dbout, iech);

    for (ifois = 0; ifois < seltab[iech]; ifois++, iiech++)
    {
      ilayer = (ifois == 0) ? static_cast<Id>(dbin->getFromLocator(ELoc::LAYER, iech)) : nlayers;

      if (ifois == 0)
      {

        /* Depth of the actual sample */

        value = dbin->getZVariable(iech, 0);
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
        dtime = dbout->getFromLocator(lmlayers->ptime, igrid, ilayer - 1);
        if (lmlayers->colreft >= 0)
          dtime -= dbout->getArray(igrid, lmlayers->colreft);
        value /= dtime;
      }
      zval[iiech] = value;
    }
  }

  if (get_keypone("Bayes_Debug_Flag", 0))
    set_keypair("Mlayers_Zval_Matrix", 1, lmlayers->neq, 1, zval.data());
}

/****************************************************************************/
/*!
 **  Calculate the Drift and subtract it from the Data
 **
 ** \param[in]  lmlayers  LMlayers structure
 ** \param[in]  verbose   True for a  verbose option
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
static Id st_subtract_optimal_drift(LMlayers* lmlayers,
                                    bool verbose,
                                    Db* dbin,
                                    DbGrid* dbout,
                                    VectorInt& seltab,
                                    VectorDouble& zval)
{
  double drext;
  Id iech, iiech, ifois, ilayer, ipos;
  Id flag_subtract = 1;
  VectorDouble coor(2);

  /* Initializations */

  Id error   = 1;
  Id nlayers = lmlayers->nlayers;
  Id nbfl    = lmlayers->nbfl;
  Id neq     = nbfl * nlayers;

  /* Core allocation */

  VectorDouble props(nlayers);
  VectorDouble drift(neq);
  VectorDouble coeff(neq);
  VectorDouble atab(neq * neq);
  VectorDouble btab(neq);

  /* Find the vector of optimal mean values */

  for (iech = iiech = 0; iech < dbin->getNSample(); iech++)
  {
    if (seltab[iech] == 0) continue;

    for (ifois = 0; ifois < seltab[iech]; ifois++, iiech++)
    {
      ilayer = (ifois == 0) ? static_cast<Id>(dbin->getFromLocator(ELoc::LAYER, iech)) : nlayers;

      /* Evaluate the proportion vector */

      if (st_get_props_data(lmlayers, dbin, dbout, iech, ilayer, props))
        goto label_end;

      /* Get the coordinates of the samples */

      coor[0] = dbin->getCoordinate(iech, 0);
      coor[1] = dbin->getCoordinate(iech, 1);

      /* Get the drift vector */

      ipos = 0;
      for (Id i = 0; i < nlayers; i++)
      {
        drext = st_get_drift_data(lmlayers, dbin, dbout, iech, i + 1);
        if (st_drift(lmlayers, coor, props[i], drext, &ipos, drift)) continue;
      }

      /* Calculate the contribution to the different arrays */

      for (Id k1 = 0; k1 < ipos; k1++)
      {
        btab[k1] += zval[iiech] * drift[k1];
        for (Id k2 = 0; k2 < ipos; k2++)
          ATAB(neq, k1, k2) += drift[k1] * drift[k2];
      }
    }
  }

  /* Find the optimal drift coefficients */

  if (matrix_invert(atab.data(), neq, -1)) goto label_end;
  matrix_product_safe(neq, neq, 1, atab.data(), btab.data(), coeff.data());

  /* Optional printout of the result */

  if (verbose)
    print_matrix("Estimated Drift", 0, 1, nlayers, nbfl, NULL, coeff.data());

  /* Subtract the optimal mean */

  for (iech = iiech = 0; iech < dbin->getNSample(); iech++)
  {
    if (seltab[iech] == 0) continue;

    for (ifois = 0; ifois < seltab[iech]; ifois++, iiech++)
    {
      ilayer = (ifois == 0) ? static_cast<Id>(dbin->getFromLocator(ELoc::LAYER, iech)) : nlayers;

      /* Evaluate the proportion vector */

      if (st_get_props_data(lmlayers, dbin, dbout, iech, ilayer, props))
        goto label_end;

      /* Get the coordinates of the samples */

      coor[0] = dbin->getCoordinate(iech, 0);
      coor[1] = dbin->getCoordinate(iech, 1);

      /* Get the drift vector */

      ipos = 0;
      for (Id i = 0; i < nlayers; i++)
      {
        drext = st_get_drift_data(lmlayers, dbin, dbout, iech, i + 1);
        if (st_drift(lmlayers, coor, props[i], drext, &ipos, drift)) continue;
      }

      /* Subtract the optimal estimation of the drift */

      if (flag_subtract)
        for (Id k1 = 0; k1 < ipos; k1++)
          zval[iiech] -= coeff[k1] * drift[k1];

      /* Save the results of the optimal drift */

      set_keypair("Optim_Drift_Coeffs", 1, 1, ipos, coeff.data());

      /* Print the residuals (optional) */

      if (OptDbg::query(EDbg::VARIOGRAM))
        message("Sample %d (Layer %d) - Coor = %lf %lf - Residual = %lf\n",
                iech + 1, ilayer, coor[0], coor[1], zval[iiech]);
    }
  }

  /* Set the error return code */

  error = 0;

label_end:
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
static Id st_get_close_sample(LMlayers* lmlayers,
                              Db* dbin,
                              Id iech0,
                              const VectorDouble& coor)
{
  Id iech, ilayer;
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

  for (iech = iech0 + 1; iech < dbin->getNSample(); iech++)
  {
    dx = dbin->getCoordinate(iech, 0) - coor[0];
    if (ABS(dx) > EPS) continue;
    dy = dbin->getCoordinate(iech, 1) - coor[1];
    if (ABS(dy) > EPS) continue;
    ilayer = static_cast<Id>(dbin->getFromLocator(ELoc::LAYER, iech));
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
static Id st_collocated_prepare(LMlayers* lmlayers,
                                Id iechout,
                                VectorDouble& coor,
                                Db* dbin,
                                DbGrid* dbout,
                                Model* model,
                                VectorInt& seltab,
                                MatrixSquare* a,
                                VectorDouble& zval,
                                VectorDouble& prop1,
                                VectorDouble& prop2,
                                MatrixSquare& covtab,
                                VectorDouble& b2,
                                VectorDouble& baux,
                                double* ratio)
{
  double coefa, coefz;

  (*ratio)   = 0.;
  Id neq     = lmlayers->neq;
  Id nlayers = lmlayers->nlayers;

  double botval = dbout->getArray(iechout, lmlayers->colrefb);
  if (lmlayers->colrefd >= 0)
    botval -= dbout->getArray(iechout, lmlayers->colrefd);
  if (st_get_props_result(lmlayers, dbout, iechout, nlayers, prop1)) return (1);
  double c0 = st_cij(lmlayers, model, nlayers, prop1, nlayers, prop1,
                     VectorDouble(), covtab);
  if (FFFF(c0)) return (1);

  if (st_lhs_one(lmlayers, dbin, dbout, model, seltab, iechout, nlayers, coor,
                 prop1, prop2, covtab, baux)) return (1);
  matrix_product_safe(neq, neq, 1, a->getValues().data(), baux.data(), b2.data());
  matrix_product_safe(1, neq, 1, b2.data(), zval.data(), &coefz);
  matrix_product_safe(1, neq, 1, b2.data(), baux.data(), &coefa);
  (*ratio) = (ABS(c0 - coefa) > 1.e-6) ? (botval - coefz) / (c0 - coefa) : 0.;

  return (0);
}

/****************************************************************************/
/*!
 **  Perform the estimation at the grid nodes in regular case
 **
 ** \param[in]  lmlayers  LMlayers structure
 ** \param[in]  flag_std  True if the estimation error must be calculated
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
static void st_estimate_regular(LMlayers* lmlayers,
                                bool flag_std,
                                double c00,
                                MatrixSquare* a,
                                VectorDouble& b,
                                VectorDouble& dual,
                                VectorDouble& wgt,
                                double* estim,
                                double* stdev)
{
  double c00val, stdv;

  /* Initializations */

  Id neq = lmlayers->neq;

  /* Perform the estimation (in Dual form) */

  *estim = VH::innerProduct(dual, b, neq);

  /* Perform the variance of estimation error */

  if (flag_std)
  {
    c00val = c00;
    if (FFFF(c00val))
      stdv = TEST;
    else
    {
      matrix_product_safe(neq, neq, 1, a->getValues().data(), b.data(), wgt.data());
      matrix_product_safe(1, neq, 1, b.data(), wgt.data(), &stdv);
      stdv = c00val - stdv;
      stdv = (stdv > 0) ? sqrt(stdv) : 0.;
    }
    *stdev = stdv;
  }
}

/****************************************************************************/
/*!
 **  Perform the estimation at the grid nodes in the bayesian case
 **
 ** \param[in]  lmlayers  LMlayers structure
 ** \param[in]  flag_std  True if the estimation error must be calculated
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
static void st_estimate_bayes(LMlayers* lmlayers,
                              bool flag_std,
                              double c00,
                              const MatrixSquare* acov,
                              VectorDouble& zval,
                              VectorDouble& b,
                              VectorDouble& wgt,
                              VectorDouble& post_mean,
                              VectorDouble& a0,
                              VectorDouble& cc,
                              VectorDouble& ss,
                              const VectorDouble& gs,
                              double* estim,
                              double* stdev)
{
  double *rhs, estim1, estim2, stdv;
  Id nech, npar;
  VectorDouble temp;
  VectorDouble fsf0;
  VectorDouble c2;

  /* Initializations */

  nech = lmlayers->nech;
  npar = lmlayers->npar;
  rhs  = b.data();
  VectorDouble ff0(b.begin() + nech, b.end());

  /* Core allocation */

  temp.resize(npar);
  fsf0.resize(nech);
  c2.resize(nech);

  /* Perform the estimation */

  matrix_product_safe(nech, npar, 1, a0.data(), ff0.data(), fsf0.data());
  for (Id iech = 0; iech < nech; iech++)
    c2[iech] = rhs[iech] + fsf0[iech];
  matrix_product_safe(nech, nech, 1, cc.data(), c2.data(), wgt.data());

  estim1 = VH::innerProduct(wgt, zval);
  estim2 = VH::innerProduct(ff0, post_mean);
  *estim = estim1 + estim2;

  /* Calculate the standard deviation */

  if (flag_std)
  {
    matrix_product_safe(1, nech, npar, rhs, ss.data(), temp.data());
    for (Id ipar = 0; ipar < npar; ipar++)
      temp[ipar] -= ff0[ipar];

    stdv = c00;
    for (Id iech = 0; iech < nech; iech++)
      for (Id jech = 0; jech < nech; jech++)
        stdv -= rhs[iech] * acov->getValue(iech, jech) * rhs[jech];

    for (Id ipar = 0; ipar < npar; ipar++)
      for (Id jpar = 0; jpar < npar; jpar++)
        stdv += temp[ipar] * GS(ipar, jpar) * temp[jpar];

    stdv   = (stdv > 0) ? sqrt(stdv) : 0.;
    *stdev = stdv;
  }
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
** \param[in]  flag_bayes True if the Bayesian hypothesis is used on drift coeffs
** \param[in]  flag_std   True if the estimation error must be calculated
** \param[in]  a          L.H.S. (square) inverted matrix
** \param[in]  zval       Data vector (extended)
** \param[in]  dual       Dual vector
** \param[in] prior_mean  Array of prior means
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
** \param[out] post_mean  Array of posterior means
**
*****************************************************************************/
static void st_estimate(LMlayers* lmlayers,
                        Db* dbin,
                        DbGrid* dbout,
                        Model* model,
                        VectorInt& seltab,
                        bool flag_bayes,
                        bool flag_std,
                        MatrixSquare* a,
                        VectorDouble& zval,
                        VectorDouble& dual,
                        const VectorDouble& prior_mean,
                        VectorDouble& prop1,
                        VectorDouble& prop2,
                        MatrixSquare& covtab,
                        VectorDouble& b,
                        VectorDouble& b2,
                        VectorDouble& baux,
                        VectorDouble& wgt,
                        VectorDouble& c00,
                        VectorDouble& /*fftab*/,
                        VectorDouble& a0,
                        VectorDouble& cc,
                        VectorDouble& ss,
                        VectorDouble& gs,
                        VectorDouble& post_mean)
{
  DECLARE_UNUSED(prior_mean);
  double estim, cx, coefb, botval, ratio, stdv;
  Id iechout, ilayer, flag_correc, nlayers, neq;
  VectorDouble coor(2);

  /* Loop on the grid nodes */

  nlayers = lmlayers->nlayers;
  neq     = lmlayers->neq;
  coefb = ratio = botval = stdv = 0.;
  if (flag_std && !lmlayers->flag_cumul)
    st_covariance_c00(lmlayers, model, VectorDouble(), covtab, c00);

  for (iechout = 0; iechout < dbout->getNSample(); iechout++)
  {
    OptDbg::setCurrentIndex(iechout + 1);
    if (!dbout->isActive(iechout)) continue;
    coor[0] = dbout->getCoordinate(iechout, 0);
    coor[1] = dbout->getCoordinate(iechout, 1);
    if (OptDbg::query(EDbg::KRIGING) || OptDbg::query(EDbg::NBGH) || OptDbg::query(EDbg::RESULTS))
    {
      mestitle(1, "Target location");
      db_sample_print(dbout, iechout, 1, 0, 0, 0);
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
      if (OptDbg::query(EDbg::KRIGING) || OptDbg::query(EDbg::NBGH))
        mestitle(2, "Layer #%d", ilayer + 1);

      /* Find the proportions for the target if flag_cumul */

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
      if (OptDbg::query(EDbg::KRIGING))
        krige_rhs_print(1, lmlayers->nech, neq, neq, NULL, b.data());

      /* Perform estimation */

      estim = stdv = TEST;
      if (flag_bayes)
        st_estimate_bayes(lmlayers, flag_std, c00[ilayer], a, zval, b, wgt,
                          post_mean, a0, cc, ss, gs, &estim, &stdv);
      else
        st_estimate_regular(lmlayers, flag_std, c00[ilayer], a, b, dual, wgt,
                            &estim, &stdv);

      /* Perform the correction (in case of collocated bottom) */

      if (flag_correc)
      {
        cx = st_ci0(lmlayers, model, nlayers, prop1, ilayer + 1,
                    VectorDouble(), covtab);
        if (FFFF(cx)) continue;
        matrix_product_safe(1, neq, 1, b2.data(), b.data(), &coefb);
        estim += (cx - coefb) * ratio;
      }

      /* Store the result */

      dbout->setLocVariable(ELoc::Z, iechout, ilayer, estim);
      if (flag_std) dbout->setLocVariable(ELoc::Z, iechout, nlayers + ilayer, stdv);
      if (OptDbg::query(EDbg::RESULTS))
      {
        message("Estimate = %lf", ilayer + 1, estim);
        if (flag_std) message(" - Variance = %lf", stdv * stdv);
        message("\n");
      }
    }
  }
  OptDbg::setCurrentIndex(0);
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
static Id st_check_auxiliary_variables(LMlayers* lmlayers,
                                       Db* dbin,
                                       DbGrid* dbout,
                                       VectorInt& seltab)
{
  Id iech, ilayer, igrid, newval, nechtot;
  double drift, value;
  VectorDouble coor(2);

  nechtot = 0;
  for (iech = 0; iech < dbin->getNSample(); iech++)
  {
    if (seltab[iech] == 0) continue;
    coor[0] = dbin->getCoordinate(iech, 0);
    coor[1] = dbin->getCoordinate(iech, 1);
    ilayer  = static_cast<Id>(dbin->getFromLocator(ELoc::LAYER, iech));
    igrid   = st_locate_sample_in_output(lmlayers, dbin, dbout, iech);
    if (IFFFF(igrid)) goto label_suppress;

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

      if (ilayer < lmlayers->nlayers && st_get_close_sample(lmlayers, dbin, iech, coor))
        newval = 2;
    }

    /* The sample is valid */

    seltab[iech] = newval;
    nechtot += newval;
    continue;

  label_suppress:
    seltab[iech] = 0;
  }

  return (nechtot);
}

/****************************************************************************/
/*!
 **  Convert the results in the Depth scale
 **
 ** \param[in]  lmlayers  LMlayers structure
 ** \param[in]  dbout     Output Db structure
 ** \param[in]  flag_std  True if the estimation error must be calculated
 **
 ** \remarks The conversion is performed:
 ** \remarks - if the calculations have been performed in Velocity or Thickness
 ** \remarks - if the calculations have been performed in cumulative or not
 ** \remarks The standard deviation is also transformed
 **
 *****************************************************************************/
static void st_convert_results(LMlayers* lmlayers,
                               Db* dbout,
                               bool flag_std)
{
  double depth0, depth, value, stdv, time0, time, depth_prev, time_prev, delta;
  Id iechout, ilayer, nlayers;

  /* Initializations */

  nlayers = lmlayers->nlayers;
  time = stdv = 0.;

  /* If Depth converion is not required, nothing to be done */

  if (!lmlayers->flag_z) return;

  /* Loop on the target points */

  for (iechout = 0; iechout < dbout->getNSample(); iechout++)
  {

    /* Identify the reference surface */

    depth0 =
      (lmlayers->colrefd >= 0) ? dbout->getArray(iechout, lmlayers->colrefd) : 0.;

    /* Identify the reference time (for velocity) */

    time0 =
      (lmlayers->colreft >= 0) ? dbout->getArray(iechout, lmlayers->colreft) : 0.;

    depth_prev = depth0;
    time_prev  = time0;

    /* Loop on the layers */

    for (ilayer = 0; ilayer < lmlayers->nlayers; ilayer++)
    {

      /* Read the estimated value */

      value = dbout->getZVariable(iechout, ilayer);
      if (flag_std) stdv = dbout->getZVariable(iechout, nlayers + ilayer);

      if (lmlayers->flag_cumul)
      {

        /* Case when calculations have been processed in cumulative way */

        if (lmlayers->flag_vel)
        {
          time  = dbout->getFromLocator(lmlayers->ptime, iechout, ilayer);
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
          time  = dbout->getFromLocator(lmlayers->ptime, iechout, ilayer);
          delta = time - time_prev;
          depth = depth_prev + value * delta;
          if (flag_std) stdv *= delta;
        }
        else
        {
          depth = depth_prev + value;
        }
        time_prev  = time;
        depth_prev = depth;
      }

      /* Store the transformed results */

      dbout->setLocVariable(ELoc::Z, iechout, ilayer, depth);
      if (flag_std) dbout->setLocVariable(ELoc::Z, iechout, nlayers + ilayer, stdv);
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
static Id st_drift_data(LMlayers* lmlayers,
                        Db* dbin,
                        DbGrid* dbout,
                        VectorInt& seltab,
                        VectorDouble& prop1,
                        VectorDouble& fftab)
{
  Id iiech, ilayer, ipos;
  double drext;
  VectorDouble coor(2);

  /* Initializations */

  Id npar = lmlayers->npar;
  Id nech = lmlayers->nech;
  for (Id i = 0; i < npar * nech; i++)
    fftab[i] = 0.;

  for (Id iech = iiech = 0; iech < dbin->getNSample(); iech++)
  {
    if (seltab[iech] == 0) continue;
    coor[0] = dbin->getCoordinate(iech, 0);
    coor[1] = dbin->getCoordinate(iech, 1);
    for (Id ifois = 0; ifois < seltab[iech]; ifois++, iiech++)
    {
      ilayer = (ifois == 0) ? static_cast<Id>(dbin->getFromLocator(ELoc::LAYER, iech)) : lmlayers->nlayers;

      /* Evaluate the proportion vector */

      if (st_get_props_data(lmlayers, dbin, dbout, iech, ilayer, prop1))
        return (1);

      /* Loop on the second sample */

      ipos = iech * lmlayers->npar;
      for (Id i = 0; i < ilayer; i++)
      {
        drext = st_get_drift_data(lmlayers, dbin, dbout, iech, i + 1);
        if (st_drift(lmlayers, coor, prop1[i], drext, &ipos, fftab)) return (1);
      }
    }
  }

  if (get_keypone("Bayes_Debug_Flag", 0))
    set_keypair("Mlayers_Drift_Matrix", 1, nech, npar, fftab.data());
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
static Id st_drift_bayes(LMlayers* lmlayers,
                         bool verbose,
                         const VectorDouble& prior_mean,
                         const VectorDouble& prior_vars,
                         MatrixSquare* acov,
                         VectorDouble& zval,
                         VectorDouble& fftab,
                         VectorDouble& a0,
                         VectorDouble& cc,
                         VectorDouble& ss,
                         VectorDouble& gs,
                         VectorDouble& post_mean,
                         VectorDouble& post_S)
{
  Id error, npar, nech, npar2, nech2;
  VectorDouble ffc;
  VectorDouble fm1z;
  VectorDouble gg;
  VectorDouble invH;
  VectorDouble invS;

  /* Initializations */

  error = 1;
  nech  = lmlayers->nech;
  npar  = lmlayers->npar;
  npar2 = npar * npar;
  nech2 = nech * nech;

  /* Core allocation */

  VectorDouble fft(npar * nech, 0);
  ffc.resize(npar * nech);
  fm1z.resize(npar);
  gg.resize(npar * npar);
  invH.resize(npar * npar);
  invS.resize(npar * npar);

  /* Constitute the prior Variance-Covariance matrix */

  for (Id ipar = 0; ipar < npar; ipar++)
    for (Id jpar = 0; jpar < npar; jpar++)
      INVS(npar, ipar, jpar) = (ipar == jpar) ? prior_vars[ipar] : 0.;

  /* Optional printout */

  if (verbose)
  {
    print_matrix("Prior Mean", 0, 1, lmlayers->nlayers, lmlayers->nbfl, NULL, prior_mean.data());
    print_matrix("Prior Variance", 0, 1, lmlayers->npar, lmlayers->npar, NULL, invS.data());
  }

  /* Invert the Data Variance-Covariance matrix */

  if (get_keypone("Bayes_Debug_Flag", 0))
    set_keypair("Bayes_ACov_Matrix", 1, nech, nech, acov->getValues().data());
  if (matrix_invert(acov->getValues().data(), nech, -1)) goto label_end;

  /* Invert the prior Variance-Covariance matrix */

  if (get_keypone("Bayes_Debug_Flag", 0))
    set_keypair("Bayes_S_Matrix", 1, npar, npar, invS.data());
  if (matrix_invert(invS.data(), npar, -1)) goto label_end;

  /* Auxiliary calculations */

  matrix_product_safe(npar, nech, nech, fftab.data(), acov->getValues().data(), ffc.data());
  matrix_transpose(npar, nech, fftab, fft);
  matrix_product_safe(npar, nech, npar, ffc.data(), fft.data(), invH.data());
  matrix_product_safe(npar, nech, 1, ffc.data(), zval.data(), fm1z.data());
  if (get_keypone("Bayes_Debug_Flag", 0))
    set_keypair("Bayes_InvH_Matrix", 1, npar, npar, invH.data());
  if (get_keypone("Bayes_Debug_Flag", 0))
    set_keypair("Bayes_Fm1z_Matrix", 1, npar, 1, fm1z.data());

  /* Calculate the Posterior Variance-Covariance matrix */

  for (Id i = 0; i < npar2; i++)
    post_S[i] = invS[i] + invH[i];
  if (matrix_invert(post_S.data(), npar, -1)) goto label_end;
  if (get_keypone("Bayes_Debug_Flag", 0))
    set_keypair("Bayes_Post_S_Matrix", 1, npar, npar, post_S.data());

  /* Calculate the Posterior Mean vector */

  matrix_product_safe(npar, npar, 1, invS.data(), prior_mean.data(),
                      post_mean.data());
  for (Id i = 0; i < npar; i++)
    fm1z[i] += post_mean[i];
  matrix_product_safe(npar, npar, 1, post_S.data(), fm1z.data(),
                      post_mean.data());
  if (get_keypone("Bayes_Debug_Flag", 0))
    set_keypair("Bayes_Post_Mean_Matrix", 1, npar, 1, post_mean.data());

  /* Optional printout */

  if (verbose)
  {
    print_matrix("Posterior Mean", 0, 1,
                 lmlayers->nlayers, lmlayers->nbfl, NULL, post_mean.data());
    print_matrix("Posterior Variance", 0, 1,
                 lmlayers->npar, lmlayers->npar, NULL, post_S.data());
  }

  /* Modify the Data vector */

  for (Id iech = 0; iech < nech; iech++)
    for (Id ipar = 0; ipar < npar; ipar++)
      zval[iech] -= FFTAB(ipar, iech) * post_mean[ipar];
  if (get_keypone("Bayes_Debug_Flag", 0))
    set_keypair("Mlayers_Zval2_Matrix", 1, lmlayers->neq, 1, zval.data());

  /* Auxiliary arrays prepared for estimation */

  matrix_product_safe(nech, npar, npar, fft.data(), post_S.data(), a0.data());
  matrix_product_safe(nech, nech, npar, acov->getValues().data(), fft.data(), ss.data());
  for (Id i = 0; i < npar2; i++)
    invS[i] = post_S[i];
  if (matrix_invert(invS.data(), npar, -1)) goto label_end;
  for (Id i = 0; i < npar2; i++)
    gg[i] = invH[i] + invS[i];
  if (matrix_invert(gg.data(), npar, -1)) goto label_end;
  if (get_keypone("Bayes_Debug_Flag", 0))
    set_keypair("Mlayers_GG_Matrix", 1, npar, npar, gg.data());

  matrix_prod_norme(1, nech, npar, ss.data(), gg.data(), cc.data());
  for (Id i = 0; i < nech2; i++)
    cc[i] = acov->getValues().data()[i] - cc[i];
  for (Id i = 0; i < npar2; i++)
    gs[i] = invH[i] + invS[i];
  if (matrix_invert(gs.data(), npar, -1)) goto label_end;
  if (get_keypone("Bayes_Debug_Flag", 0))
    set_keypair("Mlayers_CC_Matrix", 1, nech, nech, cc.data());

  /* Set the error return code */

  error = 0;

label_end:
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
 ** \param[in]  neigh      ANeigh structure
 ** \param[in]  flag_same  True if input and output files coincide
 ** \param[in]  flag_z     True if the output must be converted back into depth
 ** \param[in]  flag_vel   True if work is performed in Velocity, False for Depth
 ** \param[in]  flag_cumul True if work is performed in Depth; False in Thickness
 ** \param[in]  flag_ext   True if external drift must be used; False otherwise
 ** \param[in]  flag_std   True if the estimation error must be calculated
 ** \param[in]  flag_bayes True if the Bayesian hypothesis is used on drift coeffs
 ** \param[in]  irf_rank   Rank of the Intrinsic Random Function (0 or 1)
 ** \param[in]  match_time True if external drift matches time; False otherwise
 ** \param[in]  dim_prior  Dimension of the prior information (for verification)
 ** \param[in]  prior_mean Vector of prior means for drift coefficients
 ** \param[in]  prior_vars Vector of prior variances for drift coefficients
 ** \param[in]  namerefd   Name of the reference Depth variable in Dbout
 ** \param[in]  namereft   Name of the reference Time variable in Dbout
 ** \param[in]  namerefb   Name of the Bottom Depth variable in Dbout (or -1)
 ** \param[in]  verbose    Verbose option
 ** \param[in]  namconv    Naming convention for output
 **
 *****************************************************************************/
Id multilayers_kriging(Db* dbin,
                       DbGrid* dbout,
                       Model* model,
                       ANeigh* neigh,
                       bool flag_same,
                       bool flag_z,
                       bool flag_vel,
                       bool flag_cumul,
                       bool flag_ext,
                       bool flag_std,
                       bool flag_bayes,
                       Id irf_rank,
                       bool match_time,
                       Id dim_prior,
                       const VectorDouble& prior_mean,
                       const VectorDouble& prior_vars,
                       const String& namerefd,
                       const String& namereft,
                       const String& namerefb,
                       bool verbose,
                       const NamingConvention& namconv)
{
  Id nlayers, ilayer, nechmax, nech, iech, neq, nvar, npar, error, iptr;

  MatrixSquare covtab;
  bool flag_created;
  ELoc ptime;
  VectorInt seltab;
  VectorDouble zval;
  VectorDouble prop1;
  VectorDouble prop2;
  VectorDouble fftab;
  VectorDouble baux;
  VectorDouble b;
  VectorDouble b2;
  VectorDouble dual;
  VectorDouble wgt;
  VectorDouble c00;
  VectorDouble a0;
  VectorDouble cc;
  VectorDouble ss;
  VectorDouble gs;
  VectorDouble post_S;
  VectorDouble post_mean;

  LMlayers* lmlayers;
  MatrixSquare acov;
  MatrixSquare atot;
  MatrixSquare* a;

  /* Preliminary checks */

  error        = 1;
  flag_created = false;
  a            = nullptr;
  lmlayers     = nullptr;
  nlayers      = model->getNVar();
  nechmax      = dbin->getNSample();
  ptime        = (match_time) ? ELoc::F : ELoc::TIME;
  global_init(dbin, dbout);
  if (krige_koption_manage(1, 1, EKrigOpt::POINT, 1, VectorInt()))
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
  if (!dbin->isNVarComparedTo(1)) goto label_end;
  if (!flag_same && !dbout->isGrid())
  {
    messerr("If Input and Output are different, Output should be a Grid Db");
    goto label_end;
  }
  if (!dbin->hasLocator(ELoc::LAYER))
  {
    messerr("The input Db must contain a LAYER locator");
    goto label_end;
  }
  if (flag_ext && nlayers != dbout->getNLoc(ELoc::F))
  {
    messerr("Inconsistency between:");
    messerr("- the number of variables in the Model (%d)", nlayers);
    messerr("- the number of external drifts in the Output Db File (%d)",
            dbout->getNLoc(ELoc::F));
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
  if (flag_std && !namerefb.empty())
  {
    messerr("Calculation of the standard deviation of the estimation error");
    messerr("has not been programmed yet in collocation case");
    goto label_end;
  }
  if (flag_bayes && !namerefb.empty())
  {
    messerr("Use of Bayesian hypothesis has not been programmed yet");
    messerr("in collocation case");
    goto label_end;
  }
  if (flag_cumul && !namerefb.empty())
  {
    messerr("Collocation option is not coded when the results are expected");
    messerr("directly expressed in Depth (rather than Thickness)");
    goto label_end;
  }
  if (prior_mean.empty() || prior_vars.empty()) flag_bayes = false;
  if (flag_bayes && dim_prior != st_get_number_drift(irf_rank, flag_ext) * nlayers)
  {
    messerr("The dimension of the Prior information (%d)", dim_prior);
    messerr("must be equal to %d (nlayers) x %d (nbfl)", nlayers,
            st_get_number_drift(irf_rank, flag_ext));
  }
  if (manageExternalInformation(1, ELoc::F, dbin, dbout, &flag_created)) goto label_end;

  /* Allocating the output variables */

  nvar = nlayers;
  if (flag_std) nvar += nlayers;
  iptr = dbout->addColumnsByConstant(nvar, TEST, String(), ELoc::Z);

  /* Fill the Multi-Layers internal structure */

  lmlayers = lmlayers_alloc(dbout, flag_same, flag_vel, flag_cumul, flag_ext, flag_z,
                            namerefd, namereft, namerefb, irf_rank, match_time,
                            nlayers);

  /* Core allocation */

  seltab.resize(nechmax);
  prop1.resize(nlayers);
  prop2.resize(nlayers);

  /* Calculate the number of active samples */

  for (iech = 0; iech < nechmax; iech++)
  {
    seltab[iech] = 0;
    ilayer       = static_cast<Id>(dbin->getFromLocator(ELoc::LAYER, iech));
    if (ilayer < 1 || ilayer > nlayers) continue;
    if (st_get_props_data(lmlayers, dbin, dbout, iech, ilayer, prop1)) continue;
    seltab[iech] = 1;
  }

  /* Check the definition of all auxiliary variables defined on output file */
  /* Count the number of active samples (including the duplicates) */

  nech           = st_check_auxiliary_variables(lmlayers, dbin, dbout, seltab);
  lmlayers->nech = nech;

  /* Allocation */

  npar = lmlayers->npar;
  neq  = lmlayers->nech + npar;
  atot.reset(neq, neq);
  acov.reset(nech, nech);
  b2.resize(neq);
  b.resize(neq);
  baux.resize(neq);
  zval.resize(neq);
  dual.resize(neq);
  wgt.resize(neq);
  covtab = MatrixSquare(nlayers);
  c00.resize(nlayers);
  if (flag_bayes)
  {
    fftab.resize(nech * npar, 0);
    a0.resize(nech * npar);
    cc.resize(nech * nech);
    ss.resize(nech * npar);
    gs.resize(npar * npar);
    post_S.resize(npar * npar);
    post_mean.resize(npar);
  }
  lmlayers->neq = neq;
  if (verbose) lmlayers_print(lmlayers);

  /* Establish the kriging matrix */

  st_lhs(lmlayers, dbin, dbout, model, seltab, prop1, prop2, covtab, atot, acov);
  if (OptDbg::isReferenceDefined() || OptDbg::query(EDbg::KRIGING))
    krige_lhs_print(nech, neq, neq, NULL, atot.getValues().data());

  /* Establish the data vector */

  st_data_vector(lmlayers, dbin, dbout, seltab, zval);
  if (OptDbg::isReferenceDefined() || OptDbg::query(EDbg::KRIGING))
  {
    mestitle(0, "Data Vector");
    message("Number of active samples  = %d\n", nech);
    message("Total number of equations = %d\n", neq);
    print_matrix("Data", 0, 1, 1, nech, NULL, zval.data());
  }

  /* Assign the Variance-Covariance matrix */

  a = (flag_bayes) ? &acov : &atot;

  /* Calculate the Posterior in the Bayesian case */

  if (flag_bayes)
  {
    if (st_drift_data(lmlayers, dbin, dbout, seltab, prop1, fftab)) goto label_end;
    if (st_drift_bayes(lmlayers, verbose, prior_mean, prior_vars, a, zval,
                       fftab, a0, cc, ss, gs, post_mean, post_S)) goto label_end;
  }
  else
  {
    if (matrix_invertFromMatrixSquare(a, neq, -1)) goto label_end;
    matrix_product_safe(neq, neq, 1, a->getValues().data(), zval.data(), dual.data());
  }

  /* Perform the estimation over the grid nodes */

  st_estimate(lmlayers, dbin, dbout, model, seltab, flag_bayes, flag_std, a,
              zval, dual, prior_mean, prop1, prop2, covtab, b, b2,
              baux, wgt, c00, fftab, a0, cc, ss, gs, post_mean);

  /* Reconstitute the surfaces (optional) */

  st_convert_results(lmlayers, dbout, flag_std);

  // Rename the output variables
  namconv.setNamesAndLocators(dbout, iptr, "estim", nvar);

  /* Set the error return code */

  error = 0;

label_end:
  (void)krige_koption_manage(-1, 1, EKrigOpt::POINT, 1, VectorInt());
  (void)manageExternalInformation(-1, ELoc::F, dbin, dbout, &flag_created);
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
static Id st_evaluate_lag(LMlayers* lmlayers,
                          Db* dbin,
                          DbGrid* dbout,
                          Vario_Order* vorder,
                          Id nlayers,
                          Id ifirst,
                          Id ilast,
                          VectorDouble& zval,
                          Id* nval,
                          double* distsum,
                          VectorInt& stat,
                          VectorDouble& phia,
                          VectorDouble& phib,
                          VectorDouble& atab,
                          VectorDouble& btab)
{
  Id iech, jech, iiech, jjech, ilayer, jlayer, ecr1, ecr2, nhalf;
  double z1, z2, dist, fact1, fact2;

  /* Local initializations */

  (*nval)    = 0;
  (*distsum) = 0.;
  nhalf      = nlayers * (nlayers + 1) / 2;
  for (Id i = 0; i < nhalf; i++)
    btab[i] = 0.;
  for (Id i = 0; i < nhalf * nhalf; i++)
    atab[i] = 0.;
  for (Id i = 0; i < nlayers * nlayers; i++)
    stat[i] = 0;

  /* Loop on the pairs contributing to the lag */

  for (Id ipair = ifirst; ipair < ilast; ipair++)
  {
    vario_order_get_indices(vorder, ipair, &iiech, &jjech, &dist);
    vario_order_get_auxiliary(vorder, ipair,
                              reinterpret_cast<char*>(&iech),
                              reinterpret_cast<char*>(&jech));
    z1 = zval[iiech];
    z2 = zval[jjech];
    (*distsum) += dist;

    ilayer = static_cast<Id>(dbin->getFromLocator(ELoc::LAYER, iech));
    if (st_get_props_data(lmlayers, dbin, dbout, iech, ilayer, phia))
      return (1);

    jlayer = static_cast<Id>(dbin->getFromLocator(ELoc::LAYER, jech));
    if (st_get_props_data(lmlayers, dbin, dbout, jech, jlayer, phib))
      return (1);

    ecr1 = 0;
    stat[(ilayer - 1) * nlayers + (jlayer - 1)] += 1;

    for (Id il1 = 0; il1 < nlayers; il1++)
      for (Id jl1 = 0; jl1 <= il1; jl1++, ecr1++)
      {
        fact1 = phia[il1] * phib[jl1];
        if (il1 != jl1) fact1 += phia[jl1] * phib[il1];
        btab[ecr1] += fact1 * z1 * z2;

        ecr2 = 0;
        for (Id il2 = 0; il2 < nlayers; il2++)
          for (Id jl2 = 0; jl2 <= il2; jl2++, ecr2++)
          {
            fact2 = phia[il2] * phib[jl2];
            if (il2 != jl2) fact2 += phia[jl2] * phib[il2];
            ATAB(nhalf, ecr1, ecr2) += fact1 * fact2;
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
 ** \param[in]  verbose    True for a verbose option
 ** \param[in]  dbin       Input Db structure
 ** \param[in]  dbout      Output Db structure
 ** \param[in]  vorder     Vario_Order structure
 ** \param[in]  zval       Data vector
 ** \param[in]  idir       Rank of the Direction
 **
 ** \param[out] vario      Vario structure
 **
 *****************************************************************************/
static Id st_varioexp_chh(LMlayers* lmlayers,
                          bool verbose,
                          Db* dbin,
                          DbGrid* dbout,
                          Vario_Order* vorder,
                          VectorDouble& zval,
                          Id idir,
                          Vario* vario)
{
  double distsum;
  Id error, nlayers, iadlag, nhalf, nhalf2, nval;
  Id ilag, number, ifirst, ilast, ilayer, jlayer, ijl;
  VectorDouble phia;
  VectorDouble phib;
  VectorDouble btab;
  VectorDouble atab;
  VectorDouble sill;
  VectorInt stat;

  /* Initializations */

  error   = 1;
  nlayers = lmlayers->nlayers;
  nhalf   = nlayers * (nlayers + 1) / 2;
  nhalf2  = nhalf * nhalf;

  /* Core allocation */

  phia.resize(nlayers);
  phib.resize(nlayers);
  btab.resize(nhalf);
  atab.resize(nhalf2);
  sill.resize(nhalf);
  stat.resize(nlayers * nlayers);

  /* Loop on the lags */

  for (ilag = 0; ilag < vario->getNLag(idir); ilag++)
  {
    vario_order_get_bounds(vorder, idir, ilag, &ifirst, &ilast);
    number = ilast - ifirst;
    if (number <= 0) continue;

    /* Loop on the pairs contributing to the lag */

    if (st_evaluate_lag(lmlayers, dbin, dbout, vorder, nlayers, ifirst, ilast,
                        zval, &nval, &distsum, stat, phia, phib,
                        atab, btab))
      goto label_end;

    if (OptDbg::query(EDbg::VARIOGRAM))
    {
      message("Lag %d\n", ilag + 1);
      print_matrix("L.H.S.", 0, 1, nhalf, nhalf, NULL, atab.data());
      print_matrix("R.H.S.", 0, 1, 1, nhalf, NULL, btab.data());
    }

    if (matrix_invert(atab.data(), nhalf, -2))
    {
      messerr("--> Inversion problem for lag %d", ilag + 1);
      if (verbose)
      {
        /* Matrix must be evaluated (as it has been destroyed by inversion) */
        (void)st_evaluate_lag(lmlayers, dbin, dbout, vorder, nlayers, ifirst,
                              ilast, zval, &nval, &distsum, stat, phia, phib,
                              atab, btab);
        messerr("Number of pairs  = %d", nval);
        messerr("Average distance = %lf", distsum);
        print_imatrix("Number of samples per layer", 0, 1, nlayers, nlayers,
                      NULL, stat.data());
        print_matrix("L.H.S.", 0, 1, nhalf, nhalf, NULL, atab.data());
        print_matrix("R.H.S.", 0, 1, 1, nhalf, NULL, btab.data());
      }
      continue;
    }
    matrix_product_safe(1, nhalf, nhalf, btab.data(), atab.data(), sill.data());

    /* Optional printout */

    if (OptDbg::query(EDbg::VARIOGRAM))
      print_trimat("C(h)", 2, nlayers, sill.data());

    /* Store the covariance values */

    ijl = 0;
    for (ilayer = 0; ilayer < nlayers; ilayer++)
      for (jlayer = 0; jlayer <= ilayer; jlayer++, ijl++)
      {
        iadlag = vario->getDirAddress(idir, ilayer, jlayer, ilag, false, 1);
        vario->setGgByIndex(idir, iadlag, sill[ijl]);
        vario->setHhByIndex(idir, iadlag, distsum);
        vario->setSwByIndex(idir, iadlag, nval);
        iadlag = vario->getDirAddress(idir, ilayer, jlayer, ilag, false, -1);
        vario->setGgByIndex(idir, iadlag, sill[ijl]);
        vario->setHhByIndex(idir, iadlag, -distsum);
        vario->setSwByIndex(idir, iadlag, nval);
      }
  }

  /* Set the error return code */

  error = 0;

label_end:
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
 ** \param[in]  flag_vel   True if work is performed in Velocity, False for Depth
 ** \param[in]  flag_ext   True if external drift must be used; False otherwise
 ** \param[in]  irf_rank   Rank of the Intrinsic Random Function (0 or 1)
 ** \param[in]  match_time True if external drift matches time; False otherwise
 ** \param[in]  namerefd   Name of the reference Depth variable in Dbout
 ** \param[in]  namereft   Name of the reference Time variable in Dbout
 ** \param[in]  verbose    True for a  verbose option
 **
 *****************************************************************************/
Id multilayers_vario(Db* dbin,
                     DbGrid* dbout,
                     Vario* vario,
                     Id nlayers,
                     bool flag_vel,
                     bool flag_ext,
                     Id irf_rank,
                     bool match_time,
                     const String& namerefd,
                     const String& namereft,
                     bool verbose)
{
  Id error, ilayer, nechmax, nech, iech, idir;
  bool flag_created;
  ELoc ptime;
  LMlayers* lmlayers;
  VectorInt seltab;
  VectorDouble prop1;
  VectorDouble zval;
  Vario_Order* vorder;

  /* Preliminary checks */

  error        = 1;
  flag_created = false;
  lmlayers     = nullptr;
  vorder       = nullptr;
  nechmax      = dbin->getNSample();
  ptime        = (match_time) ? ELoc::F : ELoc::TIME;
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
  if (!dbin->isNVarComparedTo(1)) goto label_end;
  if (!dbin->hasLocator(ELoc::LAYER))
  {
    messerr("The input Db must contain a LAYER locator");
    goto label_end;
  }
  if (flag_ext && nlayers != dbout->getNLoc(ELoc::F))
  {
    messerr("Inconsistency between:");
    messerr("- the number of variables in the Model (%d)", nlayers);
    messerr("- the number of external drifts in the Output Db File (%d)",
            dbout->getNLoc(ELoc::F));
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
  if (manageExternalInformation(1, ELoc::F, dbin, dbout, &flag_created)) goto label_end;

  /* Fill the Multi-Layers internal structure */

  lmlayers = lmlayers_alloc(dbout, 0, flag_vel, 0, flag_ext, 0, namerefd, namereft, String(),
                            irf_rank, match_time, nlayers);
  if (verbose) lmlayers_print(lmlayers);

  /* Core allocation */

  seltab.resize(nechmax);
  prop1.resize(nlayers);
  zval.resize(nechmax);

  /* Calculate the number of active samples */

  for (iech = 0; iech < nechmax; iech++)
  {
    seltab[iech] = 0;
    ilayer       = static_cast<Id>(dbin->getFromLocator(ELoc::LAYER, iech));
    if (ilayer < 1 || ilayer > nlayers) continue;
    if (st_get_props_data(lmlayers, dbin, dbout, iech, ilayer, prop1)) continue;
    seltab[iech] = 1;
  }

  /* Check the definition of all auxiliary variables defined on output file */
  /* Count the number of active samples (including the duplicates) */

  nech           = st_check_auxiliary_variables(lmlayers, dbin, dbout, seltab);
  lmlayers->nech = nech;
  lmlayers->neq  = nech + lmlayers->npar;

  /* Establish the data vector */

  st_data_vector(lmlayers, dbin, dbout, seltab, zval);

  /* Subtract the optimal average or drift */

  if (st_subtract_optimal_drift(lmlayers, verbose, dbin, dbout, seltab, zval))
    goto label_end;

  /* Evaluate the Geometry */

  vorder = vario_order_manage(1, 1, sizeof(Id), NULL);
  if (vario->computeGeometryMLayers(dbin, seltab, vorder)) goto label_end;

  /* Evaluate the variogram */

  for (idir = 0; idir < vario->getNDir(); idir++)
  {
    if (st_varioexp_chh(lmlayers, verbose, dbin, dbout, vorder, zval, idir,
                        vario)) goto label_end;
  }

  /* Set the error return code */

  error = 0;

label_end:
  (void)manageExternalInformation(-1, ELoc::F, dbin, dbout, &flag_created);
  vario_order_manage(-1, 1, sizeof(Id), vorder);
  lmlayers_free(lmlayers);
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
static Id st_get_prior(Id nech,
                       Id npar,
                       VectorDouble& zval,
                       VectorDouble& fftab,
                       VectorDouble& mean,
                       VectorDouble& vars)
{
  MatrixSymmetric atab(npar);
  MatrixSymmetric atab0(npar);
  VectorDouble btab(npar);
  VectorDouble btab0(npar);
  VectorDouble result(npar);

  atab0.fill(0.);
  btab0.fill(0.);
  for (Id i = 0; i < npar; i++)
    mean[i] = 0.;
  for (Id i = 0; i < npar * npar; i++)
    vars[i] = 0.;

  /* Loop on the data */

  for (Id iech = 0; iech < nech; iech++)
    for (Id ipar = 0; ipar < npar; ipar++)
    {
      btab0[ipar] += zval[iech] * FFTAB(ipar, iech);
      for (Id jpar = 0; jpar <= ipar; jpar++)
        atab0.updValue(ipar, jpar, EOperator::ADD, FFTAB(ipar, iech) * FFTAB(jpar, iech));
    }

  /* Optional printout */

  if (get_keypone("Bayes_Debug_Flag", 0))
  {
    set_keypair("Bayes_Get_Prior_ATAB0", 1, npar * npar, 1, atab0.getValues().data());
    set_keypair("Bayes_Get_Prior_BTAB0", 1, npar, 1, btab0.data());
  }

  /* Bootstrap for the variance-covariance */

  for (Id iech = 0; iech < nech; iech++)
  {
    btab = btab0;
    atab = atab0;

    /* Update the arrays by suppressing the current data */

    for (Id ipar = 0; ipar < npar; ipar++)
    {
      btab[ipar] -= zval[iech] * FFTAB(ipar, iech);
      for (Id jpar = 0; jpar <= ipar; jpar++)
        atab.setValue(ipar, jpar, atab.getValue(ipar, jpar) - FFTAB(ipar, iech) * FFTAB(jpar, iech));
    }

    /* Solve the system */

    if (atab.solve(btab, result)) return 1;

    /* Update the statistics */

    for (Id ipar = 0; ipar < npar; ipar++)
    {
      mean[ipar] += result[ipar];
      for (Id jpar = 0; jpar < npar; jpar++)
        VARS(npar, ipar, jpar) += result[ipar] * result[jpar];
    }
  }

  /* Normalize the results */

  for (Id ipar = 0; ipar < npar; ipar++)
    mean[ipar] /= nech;
  for (Id ipar = 0; ipar < npar; ipar++)
    for (Id jpar = 0; jpar < npar; jpar++)
      VARS(npar, ipar, jpar) = VARS(npar, ipar, jpar) / nech - mean[ipar] * mean[jpar];

  return 0;
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
 ** \param[in]  flag_same  True if input and output files coincide
 ** \param[in]  flag_vel   True if work is performed in Velocity, False for Depth
 ** \param[in]  flag_ext   True if external drift must be used; False otherwise
 ** \param[in]  irf_rank   Rank of the Intrinsic Random Function (0 or 1)
 ** \param[in]  match_time True if external drift matches time; False otherwise
 ** \param[in]  namerefd   Name of the reference Depth variable in Dbout
 ** \param[in]  namereft   Name of the reference Time variable in Dbout
 ** \param[in]  namerefb   Name of the Bottom Depth variable in Dbout (or -1)
 ** \param[in]  verbose    Verbose option
 **
 *****************************************************************************/
Id multilayers_get_prior(Db* dbin,
                         DbGrid* dbout,
                         Model* model,
                         bool flag_same,
                         bool flag_vel,
                         bool flag_ext,
                         Id irf_rank,
                         bool match_time,
                         const String& namerefd,
                         const String& namereft,
                         const String& namerefb,
                         bool verbose)
{
  Id nlayers, ilayer, nechmax, nech, iech, npar, error, neq;
  bool flag_created;
  VectorInt seltab;
  VectorDouble zval;
  VectorDouble props;
  VectorDouble fftab;
  ELoc ptime;
  LMlayers* lmlayers;
  VectorDouble mean;
  VectorDouble vars;

  /* Preliminary checks */

  error        = 1;
  flag_created = false;
  lmlayers     = nullptr;
  nlayers      = model->getNVar();
  nechmax      = dbin->getNSample();
  ptime        = (match_time) ? ELoc::F : ELoc::TIME;
  if (krige_koption_manage(1, 1, EKrigOpt::POINT, 1, VectorInt())) goto label_end;
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
  if (!dbin->isNVarComparedTo(1)) goto label_end;
  if (!flag_same && !dbout->isGrid())
  {
    messerr("If Input and Output are different, Output should be a Grid Db");
    goto label_end;
  }
  if (!dbin->hasLocator(ELoc::LAYER))
  {
    messerr("The input Db must contain a LAYER locator");
    goto label_end;
  }
  if (flag_ext && nlayers != dbout->getNLoc(ELoc::F))
  {
    messerr("Inconsistency between:");
    messerr("- the number of variables in the Model (%d)", nlayers);
    messerr("- the number of external drifts in the Output Db File (%d)",
            dbout->getNLoc(ELoc::F));
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
  if (manageExternalInformation(1, ELoc::F, dbin, dbout, &flag_created)) goto label_end;

  /* Fill the Multi-Layers internal structure */

  lmlayers = lmlayers_alloc(dbout, flag_same, flag_vel, 0, flag_ext, 1, namerefd,
                            namereft, namerefb, irf_rank, match_time, nlayers);

  /* Core allocation */

  seltab.resize(nechmax);
  props.resize(nlayers);

  /* Calculate the number of active samples */

  for (iech = 0; iech < nechmax; iech++)
  {
    seltab[iech] = 0;
    ilayer       = static_cast<Id>(dbin->getFromLocator(ELoc::LAYER, iech));
    if (ilayer < 1 || ilayer > nlayers) continue;
    if (st_get_props_data(lmlayers, dbin, dbout, iech, ilayer, props)) continue;
    seltab[iech] = 1;
  }

  /* Check the definition of all auxiliary variables defined on output file */
  /* Count the number of active samples (including the duplicates) */

  nech           = st_check_auxiliary_variables(lmlayers, dbin, dbout, seltab);
  lmlayers->nech = nech;
  lmlayers->neq  = nech + lmlayers->npar;
  if (verbose) lmlayers_print(lmlayers);

  /* Allocation */

  npar = lmlayers->npar;
  neq  = lmlayers->nech + npar;
  zval.resize(neq);
  fftab.resize(nech * npar);
  mean.resize(npar);
  vars.resize(npar * npar);

  /* Establish the data vector */

  st_data_vector(lmlayers, dbin, dbout, seltab, zval);

  /* Establish the drift matrix */

  if (st_drift_data(lmlayers, dbin, dbout, seltab, props, fftab)) goto label_end;

  /* Estimate the optimal drift matrices */

  if (st_get_prior(nech, npar, zval, fftab, mean, vars)) goto label_end;

  /* Print the resulting values */
  message("Number of parameters = %d\n", npar);
  VH::dump("Means", mean);
  VH::dump("Variances", vars);

  /* Set the error return code */

  error = 0;

label_end:
  (void)krige_koption_manage(-1, 1, EKrigOpt::POINT, 1, VectorInt());
  (void)manageExternalInformation(-1, ELoc::F, dbin, dbout, &flag_created);
  lmlayers = lmlayers_free(lmlayers);
  return (error);
}
} // namespace gstlrn