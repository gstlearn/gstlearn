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
#include "Model/Constraints.hpp"
#include "Basic/AException.hpp"
#include "Basic/File.hpp"
#include "Basic/Utilities.hpp"
#include "Covariances/CovAniso.hpp"
#include "Model/Option_AutoFit.hpp"
#include "Model/Model.hpp"
#include "Model/ConsItem.hpp"
#include "Basic/EJustify.hpp"
#include "Basic/String.hpp"
#include "Db/Db.hpp"
#include "Variogram/Vario.hpp"

#include <math.h>

/*! \cond */
#define TAKE_ROT       (( optvar.getLockSamerot() && first_covrot < 0) ||  \
                        ! optvar.getLockSamerot())
#define DEFINE_AIC     (! optvar.getFlagGoulardUsed())
#define DEFINE_THIRD   (flag_param)
#define DEFINE_RANGE   (flag_range > 0)
#define DEFINE_ANICOEF (flag_range != 0 && optvar.getAuthAniso())
#define DEFINE_ANIROT  (flag_range != 0 && optvar.getAuthAniso() &&  \
                        optvar.getAuthRotation())
#define DEFINE_T_RANGE (model->getModTransMode() == EModelProperty::TAPE)
#define FLAG_COMPRESS(imod,icov) (flag_compress[(imod) * ncova + (icov)])

#define COMP_INDEX(i,j)          ((i) * ((i)+1) / 2 + (j))
#define IJDIR(ijvar,ipadir)      ((ijvar) * npadir + (ipadir))
#define SILL(icov,ijvar)         sill[(ijvar) + (icov) * nvs2]
#define GG(ijvar,ipadir)         gg[IJDIR(ijvar,ipadir)]
#define GG2(ijvar,ipadir)        gg2[IJDIR(ijvar,ipadir)]
#define GE1(ijvar,ipadir)        ge1[IJDIR(ijvar,ipadir)]
#define WT(ijvar,ipadir)         wt[IJDIR(ijvar,ipadir)]
#define WT2(ijvar,ipadir)        wt2[IJDIR(ijvar,ipadir)]
#define MP(ijvar,ipadir)         mp[IJDIR(ijvar,ipadir)]
#define GE(icov,ijvar,ipadir)    ge[IJDIR(ijvar,ipadir)  + (icov)*nvs2*npadir]
#define GE2(icov,ijvar,ipadir)   ge2[IJDIR(ijvar,ipadir) + (icov)*nvs2*npadir]
#define FK(icov,ijvar,ipadir)    fk[IJDIR(ijvar,ipadir)  + (icov)*nvs2*npadir]
#define DD(idim,ijvar,ipadir)    dd[IJDIR(ijvar,ipadir)  + (idim)*nvs2*npadir]
#define TAB(ijvar,ipadir)        tabin[IJDIR(ijvar,ipadir)]

#define AD(ivar,jvar)            (ivar) + nvar * (jvar)
#define VARS(ivar,jvar)          vario->vars[(ivar) * vario->getNVar() + (jvar)]
#define VECPRO(ivar,jvar)        vecpro[AD(ivar,jvar)]
#define VARCHOL(ivar,jvar)       varchol[(ivar) * nvar + (jvar)]

#define AA(icov,jcov)            aa[(icov) * ncova + (jcov)]
#define CC(ivar,jvar)            cc[AD(ivar,jvar)]
#define AUX(ivar,jvar)           aux[AD(ivar,jvar)]

#define AIC(icov,ijvar)          aic[(icov)*nvs2 + (ijvar)]
#define ALPHAK(icov,ijvar)       alphak[(icov)*nvs2 + (ijvar)]

#define CORRECT(idir,k)         (vario->getHhByRank(idir,k) != 0. && ! FFFF(vario->getHhByRank(idir,k)) && \
                                 vario->getSwByRank(idir,k) != 0. && ! FFFF(vario->getSwByRank(idir,k)) && \
                                 ! FFFF(vario->getGgByRank(idir,k)))
#define VALPRO(ivar)             valpro[(ivar)]
#define MATCOR(icov,ivar,jvar)   matcor[(icov)*nvar*nvar  + AD(ivar,jvar)]
#define MATCORU(icov,ivar,jvar)  matcoru[(icov)*nvar*nvar  + AD(ivar,jvar)]
#define AIRKV(ipadir,ivar)       Airk_v[(ipadir)*nvar + (ivar)]
#define BIRKV(ipadir,ivar)       Birk_v[(ipadir)*nvar + (ivar)]
#define ALPHA(icov,ivar,jvar)    alpha[(icov)*nvar*nvar + AD(ivar,jvar)]
#define MUOLD(ivar,jvar)         muold[AD(ivar,jvar)]

typedef struct
{
  int npadir;
  VectorDouble gg;
  VectorDouble ge;
  VectorDouble wt;
  VectorDouble sill;
  VectorDouble covtab;
  VectorDouble wtc;
  VectorDouble ggc;
  VectorDouble dd;
  VectorDouble wt2;
  VectorDouble ge1;
  VectorDouble ge2;
  VectorDouble gg2;
  VectorDouble alphau;
  VectorDouble sill1;
} Recint;

typedef struct
{
  int ivar;
  int jvar;
  VectorDouble dd;
} StrExp;

typedef struct
{
  int flag_regularize;
  int ndim;
  int ndisc[3];
  double support[3];
} Regularize;

/*! \endcond */

// TODO : rename this (and remove static string below)
static char string[STRING_LENGTH];
static char cov_name[STRING_LENGTH];

static int CONGRUENCY = 50;
static double EpsFit = 1.e-12;

static Regularize REGULARIZE;
static std::vector<StrExp> STREXPS;
static StrMod *STRMOD = nullptr;
static Option_AutoFit MAUTO;
static int *INDG1;
static int *INDG2;
static const Db *DBMAP;
static void (*ST_PREPAR_GOULARD)(int imod);
static Recint RECINT;

/****************************************************************************/
/*!
 **  Compute the name of the range
 **
 ** \param[in]  ivar Rank of the variable
 **
 *****************************************************************************/
static void st_name_range(int ivar)
{
  if (ivar == 0)
    (void) gslStrcpy(string, "Range U");
  else if (ivar == 1)
    (void) gslStrcpy(string, "Range V");
  else if (ivar == 2)
    (void) gslStrcpy(string, "Range W");
  else
    (void) gslSPrintf(string, "Range in direction %d", ivar);
}

/****************************************************************************/
/*!
 **  Compute the name of the scale factor
 **
 ** \param[in]  ivar Rank of the variable
 **
 *****************************************************************************/
static void st_name_scale(int ivar)
{
  if (ivar == 0)
    (void) gslStrcpy(string, "Scale U");
  else if (ivar == 1)
    (void) gslStrcpy(string, "Scale V");
  else if (ivar == 2)
    (void) gslStrcpy(string, "Scale W");
  else
    (void) gslSPrintf(string, "Scale in direction %d", ivar);
}

/****************************************************************************/
/*!
 **  Compute the name of the rotation
 **
 ** \param[in]  rank  Rank of the angle
 **
 *****************************************************************************/
static void st_name_rotation(int rank)
{
  if (rank == 0)
    (void) gslStrcpy(string, "Anisotropy Rotation Angle around Oz");
  else if (rank == 1)
    (void) gslStrcpy(string, "Anisotropy Rotation Angle around Oy");
  else if (rank == 2)
    (void) gslStrcpy(string, "Anisotropy Rotation Angle around Ox");
  else
    (void) gslSPrintf(string, "Anisotropy Rotation Angle %d", rank);
}

/****************************************************************************/
/*!
 **  Decode the parameter identificator
 **
 ** \param[in]   parid     Parameter identificator
 **
 ** \param[out]  imod      Rank of the model
 ** \param[out]  icov      Rank of the covariance
 ** \param[out]  icons     Type of the constraint (EConsElem)
 ** \param[out]  ivar      Rank of the first index
 ** \param[out]  jvar      Rank of the second index
 **
 *****************************************************************************/
static void st_parid_decode(int parid,
                            int *imod,
                            int *icov,
                            EConsElem *icons,
                            int *ivar,
                            int *jvar)
{
  int value, divide, iic;

  value = parid;
  divide = value / CONGRUENCY;
  *jvar = value - divide * CONGRUENCY;
  value = divide;
  divide = value / CONGRUENCY;
  *ivar = value - divide * CONGRUENCY;
  value = divide;
  divide = value / CONGRUENCY;
  iic = value - divide * CONGRUENCY;
  value = divide;
  divide = value / CONGRUENCY;
  *icov = value - divide * CONGRUENCY;
  value = divide;
  divide = value / CONGRUENCY;
  *imod = value - divide * CONGRUENCY;
  *icons = EConsElem::fromValue(iic);
  return;
}

/****************************************************************************/
/*!
 **  Blank out an array of values
 **
 ** \param[in]  tab   Array to be initialized
 ** \param[in]  ntab  Size of the Array
 **
 *****************************************************************************/
static void st_blank(VectorDouble &tab, int ntab)
{
  int i;

  for (i = 0; i < ntab; i++)
    tab[i] = TEST;
}

/****************************************************************************/
/*!
 **  Encode the parameter identificator
 **
 ** \return  Encoded parameter identificator
 **
 ** \param[in]  imod      Rank of the model
 ** \param[in]  icov      Rank of the covariance
 ** \param[in]  icons     Type of the constraint (EConsElem)
 ** \param[in]  ivar      Rank of the first index
 ** \param[in]  jvar      Rank of the second index
 **
 *****************************************************************************/
static int st_parid_encode(int imod,
                           int icov,
                           const EConsElem &icons,
                           int ivar,
                           int jvar)
{
  int value;

  value = imod;
  value = value * CONGRUENCY + icov;
  value = value * CONGRUENCY + icons.getValue();
  value = value * CONGRUENCY + ivar;
  value = value * CONGRUENCY + jvar;

  return (value);
}

/****************************************************************************/
/*!
 **  Define the list of parameters identificators
 **
 ** \return Number of parameters
 **
 ** \param[in]  strmod    StrMod structure
 ** \param[in]  npar0     Initial number of parameters to be inferred
 **
 *****************************************************************************/
static int st_parid_alloc(StrMod *strmod, int npar0)
{
  int ntot, jcov, flag_range, flag_param, flag_aniso, flag_rotation, ivar, jvar;
  int min_order, max_ndim, flag_int_1d, flag_int_2d, idim, ndim, nvar;
  int first_covrot, imod;
  double scalfac, parmax;
  Option_VarioFit optvar;
  Model *model;

  /* Initializations */

  optvar = strmod->optvar;
  ndim = strmod->models[0]->getDimensionNumber();
  nvar = strmod->models[0]->getVariableNumber();
  first_covrot = -1;

  /* Core allocation */

  strmod->parid.resize(npar0);

  /* Loop on the models */

  ntot = 0;
  for (imod = 0; imod < strmod->nmodel; imod++)
  {
    model = strmod->models[imod];

    /* Loop on the basic structures */

    for (jcov = 0; jcov < model->getCovaNumber(); jcov++)
    {
      model_cova_characteristics(model->getCovaType(jcov), cov_name,
                                 &flag_range, &flag_param, &min_order,
                                 &max_ndim, &flag_int_1d, &flag_int_2d,
                                 &flag_aniso, &flag_rotation, &scalfac,
                                 &parmax);

      /* AIC coefficients -> Sill */
      if (DEFINE_AIC)
        for (ivar = 0; ivar < nvar; ivar++)
          for (jvar = 0; jvar <= ivar; jvar++)
            strmod->parid[ntot++] = st_parid_encode(imod, jcov, EConsElem::SILL,
                                                    ivar, jvar);

      /* Third parameter */
      if (DEFINE_THIRD)
      {
        strmod->parid[ntot++] = st_parid_encode(imod, jcov, EConsElem::PARAM, 0,
                                                0);
      }

      /* Range */
      if (DEFINE_RANGE)
      {
        strmod->parid[ntot++] = st_parid_encode(imod, jcov, EConsElem::RANGE, 0,
                                                0);
      }

      /* Anisotropy coefficients */
      if (DEFINE_ANICOEF)
      {
        if (ndim == 2)
        {
          strmod->parid[ntot++] = st_parid_encode(imod, jcov, EConsElem::RANGE,
                                                  1, 0);
        }
        else if (ndim == 3)
        {
          if (!optvar.getLockIso2d())
          {
            strmod->parid[ntot++] = st_parid_encode(imod, jcov,
                                                    EConsElem::RANGE, 1, 0);
          }
          if (!optvar.getLockNo3d())
          {
            strmod->parid[ntot++] = st_parid_encode(imod, jcov,
                                                    EConsElem::RANGE, 2, 0);
          }
        }
        else
        {
          for (idim = 1; idim < ndim; idim++)
          {
            strmod->parid[ntot++] = st_parid_encode(imod, jcov,
                                                    EConsElem::RANGE, idim, 0);
          }
        }
      }

      /* Anisotropy angles */
      if (DEFINE_ANIROT)
      {
        if (ndim == 2)
        {
          if (TAKE_ROT)
          {
            first_covrot = jcov;
            strmod->parid[ntot++] = st_parid_encode(imod, jcov,
                                                    EConsElem::ANGLE, 0, 0);
          }
        }
        else if (ndim == 3 && optvar.getLockRot2d())
        {
          if (TAKE_ROT)
          {
            first_covrot = jcov;
            strmod->parid[ntot++] = st_parid_encode(imod, jcov,
                                                    EConsElem::ANGLE, 0, 0);
          }
        }
        else
        {
          if (TAKE_ROT)
          {
            first_covrot = jcov;
            for (idim = 0; idim < ndim; idim++)
            {
              strmod->parid[ntot++] = st_parid_encode(imod, jcov,
                                                      EConsElem::ANGLE, idim,
                                                      0);
            }
          }
        }
      }

      /* Tapering Range */
      if (DEFINE_T_RANGE)
      {
        strmod->parid[ntot++] = st_parid_encode(imod, 0, EConsElem::T_RANGE, 0,
                                                0);
      }
    }
  }

  return (ntot);
}

/****************************************************************************/
/*!
 **  Free the pointer on the StrMod structure
 **
 ** \return  Newly freed pointer
 **
 ** \param[in]  strmod  Pointer to the StrMod to be freed
 **
 *****************************************************************************/
static StrMod* st_model_auto_strmod_free(StrMod *strmod)

{
  if (strmod == nullptr) return (strmod);

  /* Check that the user_data area has been freed */
  if (strmod->user_data != NULL)
  {
    messerr("The User_Data area of the StrMod structure has not been freed");
    messerr("Before the StrMod structure is released");
  }

  delete strmod;
  strmod = nullptr;
  return (strmod);
}

/****************************************************************************/
/*!
 **  Allocate the pointers on the StrMod structure
 **
 ** \return  Pointer to the newly allocated StrMod structure
 **
 ** \param[in]  model1    Model structure
 ** \param[in]  model2    Model structure
 ** \param[in]  npar0     Initial number of parameters
 ** \param[in]  norder    Order of the Generalized Variogram
 ** \param[in]  hmax      Maximum experimental distance
 ** \param[in]  angles    Reference angles (coming from the Variogram)
 ** \param[in]  optvar    Opt_Vario structure
 **
 ** \param[out] npar      Final number of parameters
 **
 *****************************************************************************/
static StrMod* st_model_auto_strmod_alloc(Model *model1,
                                          Model *model2,
                                          int npar0,
                                          int norder,
                                          double hmax,
                                          VectorDouble &angles,
                                          const Option_VarioFit &optvar,
                                          int *npar)
{
  StrMod *strmod;
  Model *model;
  int i, ncovmax, icov, ivar, jvar, nvar, error, nmodel;

  /* Initialization */

  error = 1;
  strmod = nullptr;

  /* Core allocation */

  strmod = new StrMod;

  /* Load the structure */

  strmod->npar_init = npar0;
  strmod->norder = norder;
  strmod->models[0] = model1;
  strmod->models[1] = model2;
  strmod->optvar = optvar;
  strmod->user_data = NULL;
  strmod->parid = VectorInt();
  strmod->covtab = VectorDouble();

  /* Count the number of models defined */

  ncovmax = nvar = 0;
  for (i = nmodel = 0; i < 2; i++)
  {
    model = strmod->models[i];
    if (model == nullptr) break;
    nmodel++;
    nvar = model->getVariableNumber();
    if (ncovmax < model->getCovaNumber()) ncovmax = model->getCovaNumber();

    /* Set the default value for the range */
    /* For models where range is not asked, as it is redundant with sill */

    model->setField(hmax);
    for (icov = 0; icov < model->getCovaNumber(); icov++)
    {
      // Set the default range

      CovAniso *cova = model->getCova(icov);
      cova->setRange(hmax);

      // Set the default values for the sill matrix

      for (ivar = 0; ivar < nvar; ivar++)
        for (jvar = 0; jvar < nvar; jvar++)
        {
          double sill = (ivar == jvar) ? 1. :
                                         0.;
          cova->setSill(ivar, jvar, sill);
          cova->setSill(jvar, ivar, sill);
        }

      // Set the rotation angle to the one given by the variogram
      // Used in the case where the rotation angle is not inferred

      cova->setAnisoAngles(angles);
    }
  }
  if (nmodel <= 0)
  {
    messerr("The number of models must be strictly positive");
    goto label_end;
  }
  strmod->nmodel = nmodel;

  /* Core allocation */

  strmod->covtab.resize(nvar * nvar);
  *npar = st_parid_alloc(strmod, npar0);

  /* Set the error return code */

  error = 0;

  label_end: if (error) strmod = st_model_auto_strmod_free(strmod);
  return (strmod);
}

/****************************************************************************/
/*!
 **  Return the number of experimental conditions
 **
 ** \return  Error return code if the returned argument are null
 **
 ** \param[in]  vario     Vario structure
 **
 ** \param[out] nbexp_ret  Total number of experimental values
 **                        (number of valid directions/lags/variables)
 ** \param[out] npadir_ret Total number of lags for all directions
 **
 *****************************************************************************/
static int st_get_vario_dimension(const Vario *vario,
                                  int *nbexp_ret,
                                  int *npadir_ret)

{
  int idir, nbexp, ipas, i, ivar, jvar, nvar, npadir;

  /* Initializations */

  nbexp = npadir = 0;
  nvar = vario->getVariableNumber();

  /* Calculate the total number of lags */

  for (idir = 0; idir < vario->getDirectionNumber(); idir++)
  {
    npadir += vario->getLagTotalNumber(idir);
    for (ipas = 0; ipas < vario->getLagNumber(idir); ipas++)
      for (ivar = 0; ivar < nvar; ivar++)
        for (jvar = 0; jvar <= ivar; jvar++)
        {
          i = vario->getDirAddress(idir, ivar, jvar, ipas, false, 1);
          if (CORRECT(idir, i)) nbexp++;
        }
  }

  /* Setting the return arguments */

  *nbexp_ret = nbexp;
  *npadir_ret = npadir;

  if (nbexp <= 0)
  {
    messerr("No active experimental variogram");
    return (1);
  }
  return (0);
}

/****************************************************************************/
/*!
 **  Return the number of experimental conditions
 **
 ** \return  Error return code if the returned argument are null
 **
 ** \param[in]  dbmap     Db structure
 ** \param[in]  nvar      Number of variables
 **
 ** \param[out] nbexp_ret  Total number of experimental values
 **                        (number of valid directions/lags/variables)
 ** \param[out] npadir_ret Maximum number of lags for all directions
 **
 *****************************************************************************/
static int st_get_vmap_dimension(const Db *dbmap,
                                 int nvar,
                                 int *nbexp_ret,
                                 int *npadir_ret)
{
  int nbexp, nvs2, ndef, npadir, nech, iech, ijvar;

  /* Initializations */

  nbexp = npadir = 0;
  nvs2 = nvar * (nvar + 1) / 2;
  nech = dbmap->getSampleNumber();

  /* Calculate the total number of lags */

  for (iech = 0; iech < nech; iech++)
  {
    ndef = 0;
    for (ijvar = 0; ijvar < nvs2; ijvar++)
      if (!FFFF(dbmap->getVariable(iech, ijvar))) ndef++;
    nbexp += ndef;
    if (ndef > 0) npadir++;
  }

  /* Setting the return arguments */

  *nbexp_ret = nbexp;
  *npadir_ret = npadir;

  if (nbexp <= 0)
  {
    messerr("No active experimental variogram map samples");
    return (1);
  }
  return (0);
}

/****************************************************************************/
/*!
 **  Update the parameters of Option_AutoFit structure according to the
 **  parameters of the Experimental Variogram
 **
 ** \param[in]  nvar      Number of variables
 ** \param[in]  varchol   Cholesky decomposition of the variance matrix
 ** \param[in]  mauto     Option_AutoFit
 **
 *****************************************************************************/
static void st_mauto_rescale(int nvar,
                             VectorDouble &varchol,
                             Option_AutoFit &mauto)
{
  double total;

  // Compute the mean variance over the components

  total = 0.;
  for (int ivar = 0; ivar < nvar; ivar++)
    total += VARCHOL(ivar,ivar)* VARCHOL(ivar,ivar);
  mauto.setTolred(mauto.getTolstop() * total / nvar);
  return;
}

/****************************************************************************/
/*!
 **  Manage the verbose option for iterative Goulard when called from Foxleg
 **
 ** \param[in]  mode      0 (before) or 1 (after) Goulard
 ** \param[in]  mauto     Option_AutoFit
 **
 *****************************************************************************/
static void st_goulard_verbose(int mode, Option_AutoFit &mauto)
{
  static int verbose;
  static int flag_converge;

  /* Dispatch */

  if (mode == 0)
  {
    verbose = mauto.getVerbose();
    mauto.setVerbose(0);
    flag_converge = debug_query("converge");
    debug_define("converge", 0);
  }
  else
  {
    mauto.setVerbose(verbose);
    debug_define("converge", flag_converge);
  }
}

/****************************************************************************/
/*!
 **  Manage the array of StrExp structures
 **
 ** \return  The vector of StrExp structures
 **
 ** \param[in]  nbexp     Number of items
 ** \param[in]  ndim      Space dimension
 **
 *****************************************************************************/
static std::vector<StrExp> st_strexp_manage(int nbexp, int ndim)
{
  std::vector<StrExp> strexps = std::vector<StrExp>(nbexp);
  for (int i = 0; i < nbexp; i++)
  {
    strexps[i].ivar = -1;
    strexps[i].jvar = -1;
    strexps[i].dd = VectorDouble(ndim);
  }
  return strexps;
}

/****************************************************************************/
/*!
 **  Compress the weights for the experimental variograms
 **
 ** \param[in]  vario     Vario structure
 ** \param[in]  npadir    Total number of lags
 ** \param[in]  tabin     Uncompressed array
 **
 ** \param[out] tabout    Compressed array
 **
 *****************************************************************************/
static void st_compress_array(const Vario *vario,
                              int npadir,
                              VectorDouble &tabin,
                              VectorDouble &tabout)
{
  int nvar = vario->getVariableNumber();

  int ecr = 0;
  int ipadir = 0;
  for (int idir = 0; idir < vario->getDirectionNumber(); idir++)
    for (int ipas = 0; ipas < vario->getLagNumber(idir); ipas++, ipadir++)
    {
      int ijvar = 0;
      for (int ivar = 0; ivar < nvar; ivar++)
        for (int jvar = 0; jvar <= ivar; jvar++, ijvar++)
        {
          if (FFFF(TAB(ijvar, ipadir))) continue;
          tabout[ecr++] = TAB(ijvar, ipadir);
        }
    }
}

/****************************************************************************/
/*!
 **  Find a possible value for C00 (for asymetric case only)
 **
 ** \return The value of C00
 **
 ** \param[in]  vario    Vario structure
 ** \param[in]  idir     Direction structure
 ** \param[in]  ivar     First variable
 ** \param[in]  jvar     Second variable
 **
 ** \remark This function is meant to provide a reasonable value for C00.
 ** \remark By default it returns the value of the central point
 ** \remark If this value if zero (with sw=0), it returns the closest non-zero
 ** \remark value
 **
 *****************************************************************************/
static double st_get_c00(const Vario *vario, int idir, int ivar, int jvar)
{
  int ipas, iad, iad0;

  iad = iad0 = vario->getDirAddress(idir, ivar, jvar, 0, false, 0);
  if (vario->getGgByRank(idir, iad) != 0. || vario->getSwByRank(idir, iad) > 0)
    goto label_end;

  for (ipas = 0; ipas < vario->getLagNumber(idir); ipas++)
  {
    iad = vario->getDirAddress(idir, ivar, jvar, ipas, false, 1);
    if (vario->getGgByRank(idir, iad) != 0) goto label_end;
    iad = vario->getDirAddress(idir, ivar, jvar, ipas, false, -1);
    if (vario->getGgByRank(idir, iad) != 0) goto label_end;
  }
  iad = iad0;

  label_end: return (vario->getGgByRank(idir, iad));
}

/****************************************************************************/
/*!
 **  Fill the array of pointers on the experimental conditions
 **
 ** \param[in]  vario     Vario structure
 ** \param[in]  npadir    Total number of lags
 **
 ** \param[out] strexps  Allocated array of StrExp (optional)
 ** \param[out] gg       Allocated array of experimental values
 **
 *****************************************************************************/
static void st_load_gg(const Vario *vario,
                       int npadir,
                       std::vector<StrExp> &strexps,
                       VectorDouble &gg)
{
  int idir, ecr, ipas, iad, jad, ivar, jvar, ijvar, nvar, idim, ipadir;
  double n1, n2, g1, g2, c00, dist;

  /* Initializations */

  nvar = vario->getVariableNumber();

  /* Load the Experimental conditions structure */

  ecr = 0;
  for (idir = ipadir = 0; idir < vario->getDirectionNumber(); idir++)
  {
    for (ipas = 0; ipas < vario->getLagNumber(idir); ipas++, ipadir++)
      for (ivar = ijvar = 0; ivar < nvar; ivar++)
        for (jvar = 0; jvar <= ivar; jvar++, ijvar++)
        {

          /* Calculate the variogram value */

          GG(ijvar,ipadir)= TEST;
          if (vario->getFlagAsym())
          {
            iad = vario->getDirAddress(idir,ivar,jvar,ipas,false,1);
            jad = vario->getDirAddress(idir,ivar,jvar,ipas,false,-1);
            c00 = st_get_c00(vario,idir,ivar,jvar);
            n1 = vario->getSwByRank(idir,iad);
            n2 = vario->getSwByRank(idir,jad);
            if (n1 + n2 <= 0) continue;
            g1 = vario->getGgByRank(idir,iad);
            g2 = vario->getGgByRank(idir,jad);
            if (! CORRECT(idir,iad) || ! CORRECT(idir,jad)) continue;
            GG(ijvar,ipadir) = c00 - (n1 * g1 + n2 * g2) / (n1 + n2);
            dist = (ABS(vario->getHhByRank(idir,iad)) + ABS(vario->getHhByRank(idir,jad))) / 2.;
          }
          else
          {
            iad = vario->getDirAddress(idir,ivar,jvar,ipas,false,1);
            if (! CORRECT(idir,iad)) continue;
            GG(ijvar,ipadir) = vario->getGgByRank(idir,iad);
            dist = ABS(vario->getHhByRank(idir,iad));
          }

          /* Define the item of the StrExp array (if defined) */

          if (! strexps.empty())
          {
            strexps[ecr].ivar = ivar;
            strexps[ecr].jvar = jvar;

            for (idim=0; idim<vario->getDimensionNumber(); idim++)
            strexps[ecr].dd[idim] = dist * vario->getCodir(idir,idim);
            ecr++;
          }
        }
      }

  return;
}

/****************************************************************************/
/*!
 **  Prepare the array for Goulard's algorithm
 **  in the case of Variogram calculation
 **
 ** \param[in]  imod  Rank of the model
 **
 *****************************************************************************/
static void st_prepar_goulard_vario(int imod)

{
  Model *model = STRMOD->models[imod];
  int npadir = RECINT.npadir;
  int ndim = model->getDimensionNumber();
  int nvar = model->getVariableNumber();
  int nvs2 = nvar * (nvar + 1) / 2;
  VectorDouble &dd = RECINT.dd;
  VectorDouble &ge = RECINT.ge;
  VectorDouble d0(ndim);
  VectorDouble tab(nvar * nvar);
  CovCalcMode mode;
  mode.update(0, ITEST, ECalcMember::LHS, -1, 0, 0);
  mode.setOrderVario(STRMOD->norder);

  /* Loop on the basic structures */

  for (int icov = 0; icov < model->getCovaNumber(); icov++)
  {
    mode.setKeepOnlyCovIdx(icov);

    /* Loop on the experiments */

    for (int ipadir = 0; ipadir < npadir; ipadir++)
    {
      int ijvar = 0;
      for (int ivar = 0; ivar < nvar; ivar++)
        for (int jvar = 0; jvar <= ivar; jvar++, ijvar++)
        {
          int flag_test = 0;
          for (int idim = 0; idim < ndim && flag_test == 0; idim++)
          {
            d0[idim] = DD(idim, ijvar, ipadir);
            if (FFFF(d0[idim])) flag_test = 1;
          }
          if (flag_test)
          {
            GE(icov,ijvar,ipadir)= TEST;
          }
          else
          {
            GE(icov,ijvar,ipadir) = model_calcul_cov_ij(model, mode, ivar, jvar, d0);
          }
        }
      }
    }

  return;
}

/*****************************************************************************/
/*!
 **  Calculates the values of a generic covariance model corresponding
 **  to the lags of an experimental variogram
 **
 ** \param[in]  vario   Vario structure
 ** \param[in]  model   Model structure
 ** \param[in]  npadir  Total number of lags
 **
 ** \param[out] dd      Array of distances (optional)
 ** \param[out] ge      Array of generic covariance values (optional)
 **
 *****************************************************************************/
static void st_load_ge(const Vario *vario,
                       Model *model,
                       int npadir,
                       VectorDouble &dd,
                       VectorDouble &ge)
{
  int ndim = model->getDimensionNumber();
  int ndir = vario->getDirectionNumber();
  int nvar = vario->getVariableNumber();
  int nvs2 = nvar * (nvar + 1) / 2;
  int norder = 0;
  if (vario->getCalculType() == ECalcVario::GENERAL1) norder = 1;
  if (vario->getCalculType() == ECalcVario::GENERAL2) norder = 2;
  if (vario->getCalculType() == ECalcVario::GENERAL3) norder = 3;
  VectorDouble d1(ndim);
  CovCalcMode mode;
  mode.setUnitary(true);
  mode.setAsVario(true);
  if (norder > 0) mode.setOrderVario(norder);

  /* Loop on the basic structures */

  for (int icov = 0; icov < model->getCovaNumber(); icov++)
  {
    ACov *cova = model->getCova(icov);
    for (int idim = 0; idim < ndim; idim++)
      d1[idim] = 0.;

    /* Loop on the experiments */

    int ipadir = 0;
    for (int idir = 0; idir < ndir; idir++)
    {
      for (int ipas = 0; ipas < vario->getLagNumber(idir); ipas++, ipadir++)
      {
        int ijvar = 0;
        for (int ivar = 0; ivar < nvar; ivar++)
          for (int jvar = 0; jvar <= ivar; jvar++, ijvar++)
          {
            int shift = ijvar * vario->getLagTotalNumber(idir);
            if (!ge.empty()) GE(icov,ijvar,ipadir)= 0.;

            double dist = 0.;
            if (vario->getFlagAsym())
            {
              int iad = shift + vario->getLagNumber(idir) + ipas + 1;
              int jad = shift + vario->getLagNumber(idir) - ipas - 1;
              if (!CORRECT(idir, iad) || !CORRECT(idir, jad)) continue;
              dist = (ABS(vario->getHhByRank(idir,iad)) + ABS(vario->getHhByRank(idir, jad)))
                  / 2.;
            }
            else
            {
              int iad = shift + ipas;
              if (!CORRECT(idir, iad)) continue;
              dist = ABS(vario->getHhByRank(idir, iad));
            }
            for (int idim = 0; idim < ndim; idim++)
              d1[idim] = dist * vario->getCodir(idir, idim);
            if (!ge.empty())
            GE(icov,ijvar,ipadir)= cova->eval(ivar,jvar,1.,d1,VectorDouble(),mode);

            if (!dd.empty()) for (int idim = 0; idim < ndim; idim++)
              DD(idim,ijvar,ipadir)= dist * vario->getCodir(idir,idim);
          }
        }
      }
    }
  return;
}

/*****************************************************************************/
/*!
 **  Calculates the weighting factors for each experimental variogram
 **  value of each directional variogram
 **
 ** \param[in]  vario     Vario structure
 ** \param[in]  wmode     Type of the weighting procedure
 ** \param[in]  npadir    Total number of lags
 **
 ** \param[out] wt        Array of weights attached to variogram lags
 **
 *****************************************************************************/
static void st_load_wt(const Vario *vario,
                       int wmode,
                       int npadir,
                       VectorDouble &wt)
{
  double d1, d2, dd, n1, n2, nn, total, ratio;
  int ivar, jvar, ijvar, ipas, idir, ndir, iad, jad, shift, nvs2, nvar, ipadir;

  /* Initializations */

  ndir = vario->getDirectionNumber();
  nvar = vario->getVariableNumber();
  nvs2 = nvar * (nvar + 1) / 2;
  VectorDouble flag(ndir);

  /* Determine the count of significant directions */

  for (idir = 0; idir < ndir; idir++)
  {
    flag[idir] = 0.;
    for (ipas = 0; ipas < vario->getLagNumber(idir); ipas++)
      for (ijvar = 0; ijvar < nvs2; ijvar++)
      {
        shift = ijvar * vario->getLagTotalNumber(idir);
        if (vario->getFlagAsym())
        {
          iad = shift + vario->getLagNumber(idir) + ipas + 1;
          jad = shift + vario->getLagNumber(idir) - ipas - 1;
          n1 = vario->getSwByRank(idir, iad);
          n2 = vario->getSwByRank(idir, jad);
          if (CORRECT(idir, iad)) flag[idir] += n1;
          if (CORRECT(idir, jad)) flag[idir] += n2;
        }
        else
        {
          iad = shift + ipas;
          nn = vario->getSwByRank(idir, iad);
          if (CORRECT(idir, iad)) flag[idir] += nn;
        }
      }
  }

  switch (wmode)
  {
    case 1:
      for (idir = ipadir = 0; idir < ndir; idir++)
      {
        for (ipas = 0; ipas < vario->getLagNumber(idir); ipas++, ipadir++)
        {
          if (flag[idir] == 0.) continue;
          for (ijvar = 0; ijvar < nvs2; ijvar++)
          {
            shift = ijvar * vario->getLagTotalNumber(idir);
            if (vario->getFlagAsym())
            {
              iad = shift + vario->getLagNumber(idir) + ipas + 1;
              jad = shift + vario->getLagNumber(idir) - ipas - 1;
              if (CORRECT(idir,iad) && CORRECT(idir, jad))
              WT(ijvar,ipadir)= flag[idir];
            }
            else
            {
              iad = shift + ipas;
              if (CORRECT(idir,iad))
              WT(ijvar,ipadir) = flag[idir];
            }
          }
        }
      }
      break;

      case 2:
      for (idir=ipadir=0; idir<ndir; idir++)
      {
        for (ipas=0; ipas<vario->getLagNumber(idir); ipas++,ipadir++)
        {
          if (flag[idir] == 0.) continue;
          for (ijvar=0; ijvar<nvs2; ijvar++)
          {
            shift = ijvar * vario->getLagTotalNumber(idir);
            if (vario->getFlagAsym())
            {
              iad = shift + vario->getLagNumber(idir) + ipas + 1;
              jad = shift + vario->getLagNumber(idir) - ipas - 1;
              if (! (CORRECT(idir,iad) && CORRECT(idir,jad))) continue;
              n1 = vario->getSwByRank(idir,iad);
              n2 = vario->getSwByRank(idir,jad);
              d1 = ABS(vario->getHhByRank(idir,iad));
              d2 = ABS(vario->getHhByRank(idir,jad));
              if (d1 > 0 && d2 > 0)
              WT(ijvar,ipadir) = sqrt((n1+n2) * (n1+n2) / (n1*d1+n2*d2) / 2.);
            }
            else
            {
              iad = shift + ipas;
              if (! CORRECT(idir,iad)) continue;
              nn = vario->getSwByRank(idir,iad);
              dd = ABS(vario->getHhByRank(idir,iad));
              if (dd > 0)
              WT(ijvar,ipadir) = nn / dd;
            }
          }
        }
      }
      break;

      case 3:
      for (idir=ipadir=0; idir<ndir; idir++)
      {
        for (ipas=0; ipas<vario->getLagNumber(idir); ipas++,ipadir++)
        {
          if (flag[idir] == 0.) continue;
          for (ijvar=0; ijvar<nvs2; ijvar++)
          {
            shift = ijvar * vario->getLagTotalNumber(idir);
            if (vario->getFlagAsym())
            {
              iad = shift + vario->getLagNumber(idir) + ipas + 1;
              jad = shift + vario->getLagNumber(idir) - ipas - 1;
              if (CORRECT(idir,iad) && CORRECT(idir,jad))
              WT(ijvar,ipadir) = 1. / vario->getLagNumber(idir);
            }
            else
            {
              iad = shift + ipas;
              if (CORRECT(idir,iad))
              WT(ijvar,ipadir) = 1. / vario->getLagNumber(idir);
            }
          }
        }
      }
      break;

      default:
      for (idir=ipadir=0; idir<ndir; idir++)
      {
        for (ipas=0; ipas<vario->getLagNumber(idir); ipas++,ipadir++)
        {
          if (flag[idir] == 0.) continue;
          for (ijvar=0; ijvar<nvs2; ijvar++)
          {
            shift = ijvar * vario->getLagTotalNumber(idir);
            if (vario->getFlagAsym())
            {
              iad = shift + vario->getLagNumber(idir) + ipas + 1;
              jad = shift + vario->getLagNumber(idir) - ipas - 1;
              if (CORRECT(idir,iad) && CORRECT(idir,jad))
              WT(ijvar,ipadir) = 1.;
            }
            else
            {
              iad = shift + ipas;
              if (CORRECT(idir,iad))
              WT(ijvar,ipadir) = 1.;
            }
          }
        }
      }
      break;
    }

    /* Scaling by direction and by variable */

  for (ijvar = 0; ijvar < nvs2; ijvar++)
  {
    for (idir = ipadir = 0; idir < ndir; idir++)
    {
      total = 0.;
      for (ipas = 0; ipas < vario->getLagNumber(idir); ipas++, ipadir++)
      {
        if (flag[idir] == 0.) continue;
        if (WT(ijvar,ipadir)> 0 && ! FFFF(WT(ijvar,ipadir)))
        total += WT(ijvar,ipadir);
      }
      if (total == 0.) continue;
      ipadir -= vario->getLagNumber(idir);
      for (ipas = 0; ipas < vario->getLagNumber(idir); ipas++, ipadir++)
      {
        if (flag[idir] == 0.) continue;
        if (WT(ijvar,ipadir)> 0 && ! FFFF(WT(ijvar,ipadir)))
        WT(ijvar,ipadir) /= total;
      }
    }
  }

  /* Scaling by variable variances */

  for (ivar = ijvar = 0; ivar < nvar; ivar++)
    for (jvar = 0; jvar <= ivar; jvar++, ijvar++)
    {
      ratio =
          (vario->getVarBivar(ivar, jvar) > 0 && vario->getVarBivar(jvar, ivar) > 0) ? sqrt(
                                                                                   vario->getVarBivar(
                                                                                       ivar,
                                                                                       jvar)
                                                                                   * vario->getVarBivar(
                                                                                       jvar,
                                                                                       ivar)) :
                                                                               1.;
      for (idir = ipadir = 0; idir < ndir; idir++)
      {
        for (ipas = 0; ipas < vario->getLagNumber(idir); ipas++, ipadir++)
          if (!FFFF(WT(ijvar, ipadir))) WT(ijvar,ipadir)/= ratio;
        }
      }

  return;
}

/****************************************************************************/
/*!
 **  Display the Goulard final score
 **
 ** \param[in]  mode        0 without; 1 with constraints
 ** \param[in]  mauto       Option_AutoFit structure
 ** \param[in]  ncova       Number of basic structures
 ** \param[in]  niter       Number of iterations
 ** \param[in]  crit        Convergence criterion
 **
 *****************************************************************************/
static void st_goulard_score(const Option_AutoFit &mauto,
                             int mode,
                             int ncova,
                             int niter,
                             double crit)
{
  if (mauto.getVerbose() > 0)
  {
    if (mode == 0)
      mestitle(1, "Statistics for Goulard algorithm");
    else
      mestitle(1, "Statistics for Goulard algorithm (with sill constraints)");
    message("Number of sills fitted   = %d\n", ncova);
    message("Number of iterations     = %d/%d\n", niter, mauto.getMaxiter());
    message("Conv. criterion          = %f\n", mauto.getTolstop());
    message("Conv. criterion (scaled) = %f\n", mauto.getTolred());
    message("Convergence score        = %f\n", ABS(crit));
  }
}

/****************************************************************************/
/*!
 **  Display the title for the Goulard algorithm
 **
 ** \param[in]  nvar        Number of variables
 ** \param[in]  ncova       Number of basic structures
 **
 *****************************************************************************/
static void st_goulard_debug_title(int nvar, int ncova)
{
  int icov, ivar, jvar;
  static char loc_string[20];

  if (!debug_query("converge")) return;
  mestitle(1, "Trajectory of parameters in Goulard Algorithm");
  message("(Sti(V1-V2) : Sill for structure 'i' for variables 'V1' and 'V2'\n");
  tab_prints(NULL, 1, EJustify::RIGHT, "Iteration");
  tab_prints(NULL, 1, EJustify::RIGHT, "Score");
  for (icov = 0; icov < ncova; icov++)
    for (ivar = 0; ivar < nvar; ivar++)
      for (jvar = 0; jvar <= ivar; jvar++)
      {
        (void) gslSPrintf(loc_string, "St%d(%d-%d)", icov + 1, ivar + 1,
                          jvar + 1);
        tab_prints(NULL, 1, EJustify::RIGHT, loc_string);
      }
  message("\n");
}

/****************************************************************************/
/*!
 **  Display the current status for the Goulard algorithm
 **
 ** \param[in]  nvar        Number of variables
 ** \param[in]  ncova       Number of basic structures
 ** \param[in]  iter        Rank of the iteration
 ** \param[in]  sill        Matrix of sills
 ** \param[in]  crit        Valeu for the convergence criterion
 **
 *****************************************************************************/
static void st_goulard_debug_current(int nvar,
                                     int ncova,
                                     int iter,
                                     VectorDouble &sill,
                                     double crit)
{
  int icov, ivar, jvar, ijvar, nvs2;

  if (!debug_query("converge")) return;
  nvs2 = nvar * (nvar + 1) / 2;
  tab_printi(NULL, 1, EJustify::RIGHT, iter + 1);
  if (FFFF(crit))
    tab_prints(NULL, 1, EJustify::RIGHT, "     ");
  else
    tab_printd(NULL, 1, EJustify::RIGHT, crit);

  for (icov = 0; icov < ncova; icov++)
    for (ivar = ijvar = 0; ivar < nvar; ivar++)
      for (jvar = 0; jvar <= ivar; jvar++, ijvar++)
        tab_printg(NULL, 1, EJustify::RIGHT, SILL(icov, ijvar));
  message("\n");
}

/****************************************************************************/
/*!
 **  Save the sill matrices in the keypair mechanism
 **
 ** \param[in]  mode      1 for setting; -1 for deleting
 ** \param[in]  model     Model structure containing the basic structures
 **
 *****************************************************************************/
static void st_keypair_sill(int mode, Model *model)
{
  int ncova, nvar;
  char loc_string[100];

  if (model == nullptr) return;
  ncova = model->getCovaNumber();
  nvar = model->getVariableNumber();

  if (mode < 0)
  {
    del_keypair("Fitted_Sill", 0);
  }
  else
  {
    for (int icova = 0; icova < ncova; icova++)
    {
      (void) gslSPrintf(loc_string, "Fitted_Sill_%d", icova + 1);
      set_keypair(loc_string, 1, nvar, nvar,
                  model->getSill(icova).getValues().data());
    }
  }
}

/****************************************************************************/
/*!
 **  Store the results using the keypair mechanism
 **
 ** \param[in]  mode      1 for setting; -1 for deleting
 ** \param[in]  icov      Rank of the covariance
 ** \param[in]  nvar      Number of variables
 ** \param[in]  valpro    Array of eigen values
 ** \param[in]  vecpro    Array of eigen vectors
 **
 *****************************************************************************/
static void st_keypair_results(int mode,
                               int icov,
                               int nvar,
                               double *valpro,
                               double *vecpro)
{
  char loc_string[50];

  if (mode < 0)
  {
    del_keypair("Model_Auto_Eigen_Values", 0);
    del_keypair("Model_Auto_Eigen_Vector", 0);
  }
  else
  {
    (void) gslSPrintf(loc_string, "Model_Auto_Eigen_Values_%d", icov + 1);
    set_keypair(loc_string, 1, 1, nvar, valpro);
    (void) gslSPrintf(loc_string, "Model_Auto_Eigen_Vector_%d", icov + 1);
    set_keypair(loc_string, 1, nvar, nvar, vecpro);
  }
}

/****************************************************************************/
/*!
 **  Reset the array of sills
 **
 ** \param[in]  nvar        Number of variables
 ** \param[in]  ncova       Number of covariances
 **
 ** \param[out] sill        Array of resulting sills
 **
 *****************************************************************************/
static void st_sill_reset(int nvar, int ncova, VectorDouble &sill)
{
  int ijvar, nvs2;

  nvs2 = nvar * (nvar + 1) / 2;
  for (int icov = 0; icov < ncova; icov++)
  {
    for (int ivar = ijvar = 0; ivar < nvar; ivar++)
      for (int jvar = 0; jvar <= ivar; jvar++, ijvar++)
        SILL(icov,ijvar)= (ivar == jvar);
  }
}

/****************************************************************************/
/*!
 **  Routine for fitting a model using an experimental variogram
 **
 ** \return  Error return code
 **
 ** \param[in]  mauto       Option_AutoFit structure
 ** \param[in]  nvar        Number of variables
 ** \param[in]  ncova       Number of covariances
 ** \param[in]  npadir      Maximum number of lags for all directions
 ** \param[in]  wt          Array of weights (Dimension: npadir)
 ** \param[in]  gg          Array of experimental values (Dimension: npadir)
 ** \param[in]  ge          Array of model values (Dimension: npadir)
 **
 ** \param[out] sill        Array of resulting sills
 ** \param[out] crit_arg    Convergence criterion
 **
 ** \remark  Internal arrays:
 ** \remark  MP : Contains the current Model (ijvar,ipadir)
 **
 *****************************************************************************/
static int st_goulard_without_constraint(const Option_AutoFit &mauto,
                                         int nvar,
                                         int ncova,
                                         int npadir,
                                         VectorDouble &wt,
                                         VectorDouble &gg,
                                         VectorDouble &ge,
                                         VectorDouble &sill,
                                         double *crit_arg)
{
  int icov, ivar, jvar, kvar, ijvar, iter, ipadir, allpos, nvs2;
  double sum, sum1, sum2, temp, crit, crit_mem, value;

  /*******************/
  /* Initializations */
  /*******************/

  sum = sum1 = sum2 = 0.;
  nvs2 = nvar * (nvar + 1) / 2;

  /* Core allocation */

  VectorDouble valpro(nvar);
  VectorDouble vecpro(nvar * nvar);
  VectorDouble aic(ncova * nvs2);
  VectorDouble cc(nvar * nvar);
  VectorDouble mp(npadir * nvs2);
  VectorDouble fk(npadir * nvs2 * ncova);
  VectorDouble alphak(ncova * nvs2);
  for (int i = 0; i < npadir * nvs2 * ncova; i++)
    fk[i] = 0.;

  /********************/
  /* Pre-calculations */
  /********************/

  for (ivar = ijvar = 0; ivar < nvar; ivar++)
    for (jvar = 0; jvar <= ivar; jvar++, ijvar++)
      for (ipadir = 0; ipadir < npadir; ipadir++)
      {
        MP(ijvar,ipadir)= 0.;
        for (icov=0; icov<ncova; icov++)
        MP(ijvar,ipadir) += SILL(icov,ijvar) * GE(icov,ijvar,ipadir);
      }

  for (icov = 0; icov < ncova; icov++)
    for (ivar = ijvar = 0; ivar < nvar; ivar++)
      for (jvar = 0; jvar <= ivar; jvar++, ijvar++)
      {
        sum1 = sum2 = AIC(icov,ijvar)= 0;
        for (ipadir=0; ipadir<npadir; ipadir++)
        {
          if (FFFF(WT(ijvar,ipadir))) continue;
          temp = WT(ijvar,ipadir) * GE(icov,ijvar,ipadir);
          FK(icov,ijvar,ipadir) = temp;
          sum1 += temp * GG(ijvar,ipadir);
          sum2 += temp * GE(icov,ijvar,ipadir);
        }
        ALPHAK(icov,ijvar) = 1. / sum2;
        AIC(icov,ijvar) = sum1 * ALPHAK(icov,ijvar);
      }

  crit = crit_mem = 0.;
  for (ivar = ijvar = 0; ivar < nvar; ivar++)
    for (jvar = 0; jvar <= ivar; jvar++, ijvar++)
      for (ipadir = 0; ipadir < npadir; ipadir++)
      {
        if (FFFF(WT(ijvar, ipadir))) continue;
        temp = GG(ijvar,ipadir)- MP(ijvar,ipadir);
        value = (ivar != jvar) ? 2. :
                                 1.;
        crit += value * WT(ijvar, ipadir) * temp * temp;
      }

  /***********************/
  /* Iterative procedure */
  /***********************/

  for (iter = 0; iter < mauto.getMaxiter(); iter++)
  {

    /* Loop on the elementary structures */

    for (icov = 0; icov < ncova; icov++)
    {

      /* Establish the coregionalization matrix */

      for (ivar = ijvar = 0; ivar < nvar; ivar++)
      {
        for (jvar = 0; jvar <= ivar; jvar++, ijvar++)
        {
          sum = 0;
          for (ipadir = 0; ipadir < npadir; ipadir++)
          {
            MP(ijvar,ipadir)-= SILL(icov,ijvar) * GE(icov,ijvar,ipadir);
            sum += FK(icov,ijvar,ipadir) * MP(ijvar,ipadir);
          }
          CC(ivar,jvar)= AIC(icov,ijvar) - ALPHAK(icov,ijvar) * sum;
          CC(jvar,ivar)= CC(ivar,jvar);
        }
      }
      /* Computing and sorting the eigen values and eigen vectors */

      if (matrix_eigen(cc.data(),nvar,valpro.data(),vecpro.data())) return 1;

      /* Store the values using the keypair mechanism */

      st_keypair_results(1,icov,nvar,valpro.data(),vecpro.data());

      ivar=0;
      allpos=1;
      while((ivar<nvar) && allpos)
      {
        if (valpro[ivar++] < 0) allpos=0;
      }

      /* Calculate the new coregionalization matrix */

      for (ivar=ijvar=0; ivar<nvar; ivar++)
      {
        for (jvar=0; jvar<=ivar; jvar++, ijvar++)
        {
          if (allpos)
          {
            SILL(icov,ijvar) = CC(ivar,jvar);
          }
          else
          {
            sum = 0.;
            for (kvar=0; kvar<nvar; kvar++)
            sum += (MAX(VALPRO(kvar),0.) *
                VECPRO(ivar,kvar) * VECPRO(jvar,kvar));
            SILL(icov,ijvar) = sum;
          }
          for (ipadir=0; ipadir<npadir; ipadir++)
          MP(ijvar,ipadir) += SILL(icov,ijvar) * GE(icov,ijvar,ipadir);
        }
      }
    }

    /* Update the global criterion */

    crit_mem = crit;
    crit = 0.;
    for (ipadir=0; ipadir<npadir; ipadir++)
    for (ivar=ijvar=0; ivar<nvar; ivar++)
    for (jvar=0; jvar<=ivar; jvar++, ijvar++)
    {
      if (FFFF(WT(ijvar,ipadir))) continue;
      temp = GG(ijvar,ipadir) - MP(ijvar,ipadir);
      value = (ivar != jvar) ? 2. : 1.;
      crit += value * WT(ijvar,ipadir) * temp * temp;
    }

    /* Optional printout */

    st_goulard_debug_current(nvar,ncova,iter,sill,crit);

    /* Stopping criterion */

    if (ABS(crit) < mauto.getTolred() ||
        ABS(crit-crit_mem) / ABS(crit) < mauto.getTolred()) break;
  }

  st_goulard_score(mauto, 0, ncova, iter, crit);

  *crit_arg = crit;
  return (0);
}

/****************************************************************************/
/*!
 **  Affect values for default, lower and upper parameters
 **
 ** \param[in]  rank       Rank of the parameter
 ** \param[in]  def_val    Default distance
 ** \param[in]  lower_val  Lowest possible value
 ** \param[in]  upper_val  Highest possible value
 ** \param[in]  param      Array of current values
 ** \param[in]  lower      Array of lower values
 ** \param[in]  upper      Array of upper values
 **
 ****************************************************************************/
static void st_affect(int rank,
                      double def_val,
                      double lower_val,
                      double upper_val,
                      VectorDouble &param,
                      VectorDouble &lower,
                      VectorDouble &upper)
{

  /* Lower bound */

  if (FFFF(lower[rank]))
    lower[rank] = lower_val;
  else if (!FFFF(lower_val)) lower[rank] = MAX(lower_val, lower[rank]);

  /* Upper bound */

  if (FFFF(upper[rank]))
    upper[rank] = upper_val;
  else if (!FFFF(upper_val)) upper[rank] = MIN(upper_val, upper[rank]);

  /* Initial parameter */
  if (FFFF(param[rank])) param[rank] = (!FFFF(def_val)) ? def_val :
                                                          0.;

  if (!FFFF(lower[rank]) && !FFFF(upper[rank]))
  {
    if (param[rank] < lower[rank] || param[rank] > upper[rank])
      param[rank] = (lower[rank] + upper[rank]) / 2.;
  }
  else if (!FFFF(lower[rank]))
  {
    if (param[rank] < lower[rank]) param[rank] = lower[rank] + 1;
  }
  else if (!FFFF(upper[rank]))
  {
    if (param[rank] > upper[rank]) param[rank] = upper[rank] - 1;
  }

  return;
}

/****************************************************************************/
/*!
 **  Compress the default, lower and upper arrays
 **
 ** \return  nb_tot      : Number of parameters (input)
 **
 ** \param[in]  n_init       Initial number of parameters
 ** \param[in,out] parid     Array of parameters identificators
 ** \param[in,out] param     Array of current values
 ** \param[in,out] lower     Array of lower values
 ** \param[in,out] upper     Array of upper values
 **
 ****************************************************************************/
static int st_compress_parid(int n_init,
                             VectorInt &parid,
                             VectorDouble &param,
                             VectorDouble &lower,
                             VectorDouble &upper)
{
  int i, n_tot;

  for (i = n_tot = 0; i < n_init; i++)
  {
    if (FFFF(param[i])) continue;
    parid[n_tot] = parid[i];
    param[n_tot] = param[i];
    upper[n_tot] = upper[i];
    lower[n_tot] = lower[i];
    n_tot++;
  }

  return (n_tot);
}

/****************************************************************************/
/*!
 **  Prints a parameter of the basic structure
 **
 ** \param[in]  name      Name of the parameter
 ** \param[in]  flag_end  1 to add the end_of_line symbol - 0 otherwise
 ** \param[in]  rank      Rank of the parameter
 ** \param[in]  param     Array of current values
 ** \param[in]  lower     Array of lower values
 ** \param[in]  upper     Array of upper values
 **
 ****************************************************************************/
static void st_print(const char *name,
                     int flag_end,
                     int rank,
                     VectorDouble &param,
                     VectorDouble &lower,
                     VectorDouble &upper)
{
  /* Title */

  message("%2d - %s : ", rank + 1, name);

  /* Current value */

  if (FFFF(param[rank]))
    message("NA ");
  else
    message("%lf ", param[rank]);

  /* Lower bound */

  if (FFFF(lower[rank]))
    message("]NA,");
  else
    message("[%lf,", lower[rank]);

  /* Upper bound */

  if (FFFF(upper[rank]))
    message("NA[");
  else
    message("%lf]", upper[rank]);

  /* Add the end_of_line */

  if (flag_end) message("\n");

  return;
}

/****************************************************************************/
/*!
 **  Print the resulting parameters of the Model
 **
 ** \param[in]  flag_title  1 if the conditions must be printed
 ** \param[in]  strmod    StrMod structure
 ** \param[in]  mauto     Option_AutoFit structure
 ** \param[in]  param     Current values for parameters
 ** \param[in]  lower     Array of lower values
 ** \param[in]  upper     Array of upper values
 ** \param[in]  npar      Number of parameters to be inferred
 ** \param[in]  nbexp     Number of experiments
 **
 *****************************************************************************/
static void st_model_auto_strmod_print(int flag_title,
                                       StrMod *strmod,
                                       const Option_AutoFit &mauto,
                                       VectorDouble &param,
                                       VectorDouble &lower,
                                       VectorDouble &upper,
                                       int npar,
                                       int nbexp)
{
  int ntot, icov, ivar, jvar, ndim, nvar, imod, imod_mem, icov_mem;
  Option_VarioFit optvar;
  EConsElem icons;
  static const char *NOK[] = { "OFF", "ON" };

  /* Initializations */

  if (!(mauto.getVerbose() > 0 || debug_query("converge"))) return;
  optvar = strmod->optvar;
  ndim = strmod->models[0]->getDimensionNumber();
  nvar = strmod->models[0]->getVariableNumber();
  imod_mem = icov_mem = -1;

  /* Title */

  if (flag_title)
  {
    mestitle(0, "%s", "Optimization Conditions");
    message("- Number of variables       %d  \n", nvar);
    message("- Space dimension           %d  \n", ndim);
    message("- Number of experiments     %d  \n", nbexp);
    message("- Number of parameters      %d  \n", npar);
    message("- Constrained Minimization  %s\n",
            NOK[!FFFF(mauto.getConstantSillValue())]);
    message("- Intrinsic option          %s\n", NOK[mauto.getFlagIntrinsic()]);
    messageFlush(optvar.toString());
  }

  /* Loop on the models */

  for (ntot = 0; ntot < npar; ntot++)
  {
    st_parid_decode(strmod->parid[ntot], &imod, &icov, &icons, &ivar, &jvar);
    if (imod != imod_mem || icov != icov_mem)
    {
      if (imod != imod_mem)
      {
        if (strmod->nmodel > 1)
          mestitle(1, "Model %d", imod + 1);
        else
          mestitle(1, "Model");
      }
      message("Structure : %s\n",
              strmod->models[imod]->getCovName(icov).c_str());
    }
    imod_mem = imod;
    icov_mem = icov;

    switch (icons.toEnum())
    {
      case EConsElem::E_SILL:
        st_print("AIC", 1, ntot, param, lower, upper);
        break;

      case EConsElem::E_PARAM:
        st_print("Parameter", 1, ntot, param, lower, upper);
        break;

      case EConsElem::E_RANGE:
        st_name_range(ivar);
        st_print(string, 1, ntot, param, lower, upper);
        break;

      case EConsElem::E_SCALE:
        st_name_scale(ivar);
        st_print(string, 1, ntot, param, lower, upper);
        break;

      case EConsElem::E_ANGLE:
        st_name_rotation(ivar);
        st_print(string, 1, ntot, param, lower, upper);
        break;

      case EConsElem::E_T_RANGE:
        st_print("Tapering Range", 1, ntot, param, lower, upper);
        break;

      case EConsElem::E_TENSOR:
        st_print("Anisotropy Matrix", 1, ntot, param, lower, upper);
        break;

      default:
        messerr("Unknown constraint!\n");
        break;
    }
  }
  return;
}

/****************************************************************************/
/*!
 **  Manage the scaling factors of the parameters
 **
 ** \param[in]  strmod          StrMod structure
 ** \param[in]  npar            Number of parameters to be inferred
 ** \param[in]  hmax            Maximum experimental distance
 ** \param[in]  varchol         Cholesky decomposition of the variance matrix
 **
 ** \param[out]  scale          Array of scales
 **
 *****************************************************************************/
static void st_model_auto_scldef(StrMod *strmod,
                                 int npar,
                                 double hmax,
                                 VectorDouble &varchol,
                                 VectorDouble &scale)
{
  int icov, ivar, jvar, imod, lec;
  int flag_range, flag_param, flag_aniso, flag_rotation;
  int min_order, max_ndim, flag_int_1d, flag_int_2d, ntot;
  double scalfac, parmax, dunit, dvar;
  Model *model;
  EConsElem icons;

  /* Loop on the models */

  if (npar <= 0) return;
  for (ntot = 0; ntot < npar; ntot++)
  {
    st_parid_decode(strmod->parid[ntot], &imod, &icov, &icons, &ivar, &jvar);
    model = strmod->models[imod];
    model_cova_characteristics(model->getCovaType(icov), cov_name, &flag_range,
                               &flag_param, &min_order, &max_ndim, &flag_int_1d,
                               &flag_int_2d, &flag_aniso, &flag_rotation,
                               &scalfac, &parmax);
    switch (icons.toEnum())
    {
      case EConsElem::E_SILL:
        lec = ivar * (ivar + 1) / 2 + jvar;
        dvar = ABS(varchol[lec]) / sqrt(model->getCovaNumber());
        scale[ntot] = dvar;
        break;

      case EConsElem::E_PARAM:
        if (parmax < 0 || FFFF(parmax)) parmax = 1.;
        scale[ntot] = parmax;
        break;

      case EConsElem::E_RANGE:
        dunit = hmax / model_get_nonugget_cova(model) / 2.;
        scale[ntot] = dunit;
        break;

      case EConsElem::E_ANGLE:
        scale[ntot] = 1800.;
        break;

      case EConsElem::E_T_RANGE:
        dunit = hmax / 10.;
        scale[ntot] = dunit;
        break;

      default:
        break;
    }
  }
  return;
}

/****************************************************************************/
/*!
 **  Update default values, lower and upper bounds
 **
 ** \param[in]  strmod          StrMod structure
 ** \param[in]  npar            Number of parameters to be inferred
 ** \param[in]  constraints     Constraints structure
 **
 ** \param[out]  param          Current values for parameters
 ** \param[out]  lower          Lower values for parameters
 ** \param[out]  upper          Upper values for parameters
 **
 *****************************************************************************/
static void st_model_auto_constraints_apply(StrMod *strmod,
                                            int npar,
                                            const Constraints &constraints,
                                            VectorDouble &param,
                                            VectorDouble &lower,
                                            VectorDouble &upper)
{
  int ipar, icov, ivar, jvar, imod;
  double param_loc, lower_loc, upper_loc;
  EConsElem icons;

  /* Loop on the models */

  if (npar <= 0) return;
  for (ipar = 0; ipar < npar; ipar++)
  {
    st_parid_decode(strmod->parid[ipar], &imod, &icov, &icons, &ivar, &jvar);
    param_loc = constraints_get(constraints, EConsType::DEFAULT, imod, icov,
                                icons, ivar, jvar);
    lower_loc = constraints_get(constraints, EConsType::LOWER, imod, icov,
                                icons, ivar, jvar);
    upper_loc = constraints_get(constraints, EConsType::UPPER, imod, icov,
                                icons, ivar, jvar);
    st_affect(ipar, param_loc, lower_loc, upper_loc, param, lower, upper);
  }
  return;
}

/****************************************************************************/
/*!
 **  Manage the default values of the parameters
 **
 ** \param[in]  strmod          StrMod structure
 ** \param[in]  npar            Number of parameters to be inferred
 ** \param[in]  hmax            Maximum experimental distance
 ** \param[in]  varchol         Cholesky decomposition of the variance matrix
 ** \param[in]  angles          Reference angles
 **
 ** \param[out]  param          Current values for parameters
 ** \param[out]  lower          Lower values for parameters
 ** \param[out]  upper          Upper values for parameters
 **
 *****************************************************************************/
static void st_model_auto_pardef(StrMod *strmod,
                                 int npar,
                                 double hmax,
                                 VectorDouble &varchol,
                                 VectorDouble &angles,
                                 VectorDouble &param,
                                 VectorDouble &lower,
                                 VectorDouble &upper)
{
  int icov, icovm, ivar, jvar, imod, lec;
  int flag_range, flag_param, flag_aniso, flag_rotation;
  int min_order, max_ndim, flag_int_1d, flag_int_2d, ntot;
  double scalfac, parmax, dist, dunit, dmin, dvar, valdef;
  Model *model;
  ECov type;
  EConsElem icons;

  /* Loop on the models */

  icovm = 0;
  if (npar <= 0) return;
  for (ntot = 0; ntot < npar; ntot++)
  {
    st_parid_decode(strmod->parid[ntot], &imod, &icov, &icons, &ivar, &jvar);
    model = strmod->models[imod];
    type = model->getCovaType(icov);
    model_cova_characteristics(type, cov_name, &flag_range, &flag_param,
                               &min_order, &max_ndim, &flag_int_1d,
                               &flag_int_2d, &flag_aniso, &flag_rotation,
                               &scalfac, &parmax);
    if (type == ECov::NUGGET && ivar == 0 && jvar == 0) icovm++;
    switch (icons.toEnum())
    {
      case EConsElem::E_SILL:
        lec = ivar * (ivar + 1) / 2 + jvar;
        dvar = varchol[lec] / sqrt(model->getCovaNumber());
        st_affect(ntot, dvar, TEST, TEST, param, lower, upper);
        break;

      case EConsElem::E_PARAM:
        if (parmax < 0) parmax = TEST;
        valdef = 1.;
        if (type == ECov::COSEXP) valdef = hmax / 3.;
        st_affect(ntot, valdef, 0.001, parmax, param, lower, upper);
        break;

      case EConsElem::E_RANGE:
        dunit = hmax / model_get_nonugget_cova(model) / 2.;
        dmin = hmax / 1.e6;
        dist = dunit * (icov + 1 - icovm);
        st_affect(ntot, dist, dmin, TEST, param, lower, upper);
        break;

      case EConsElem::E_ANGLE:
        st_affect(ntot, angles[ivar], TEST, TEST, param, lower, upper);
        break;

      case EConsElem::E_T_RANGE:
        dmin = hmax / 1.e6;
        dist = hmax / 10.;
        st_affect(ntot, dist, dmin, TEST, param, lower, upper);
        break;

      default:
        break;
    }
  }
  return;
}

/****************************************************************************/
/*!
 **  Load the parameters in the Model structure
 **
 ** \param[in]  strmod  StrMod structure
 ** \param[in]  npar    Number of parameters
 ** \param[in]  param   Current values for parameters
 **
 *****************************************************************************/
static void st_model_auto_strmod_define(StrMod *strmod,
                                        int npar,
                                        VectorDouble &param)
{
  int icov, nvar, ntot, imod, imod_mem, icov_mem, ivar, jvar, size;
  int flag_rot, flag_aic, found, ipos, ndim;
  Model *model;
  CovAniso *cova;
  CovAniso *cova1;
  Tapering *tape;
  Option_VarioFit optvar;
  EConsElem icons;

  /* Initializations */

  nvar = strmod->models[0]->getVariableNumber();
  ndim = strmod->models[0]->getDimensionNumber();
  optvar = strmod->optvar;
  size = nvar * (nvar + 1) / 2;
  VectorDouble ranges(ndim, 0.);
  VectorDouble angles(ndim, 0.);
  VectorDouble tritab(size);

  /* Loop on the parameters */

  flag_rot = flag_aic = 0;
  imod_mem = icov_mem = -1;
  for (ntot = 0; ntot < npar; ntot++)
  {
    st_parid_decode(strmod->parid[ntot], &imod, &icov, &icons, &ivar, &jvar);

    // Store the global parameters for the previous structure

    if ((imod != imod_mem || icov != icov_mem) && (imod_mem >= 0
        && icov_mem >= 0))
    {
      cova = strmod->models[imod_mem]->getCova(icov_mem);
      if (optvar.getAuthAniso())
        cova->setRanges(ranges);
      else
        cova->setRange(ranges[0]);
      if (flag_rot) cova->setAnisoAngles(angles);
      if (flag_aic)
      {
        VectorDouble sill = matrix_produit_lu_VD(nvar, tritab.data());
        MatrixSquareGeneral mat(nvar);
        mat.setValues(sill);
        cova->setSill(mat);
      }
      flag_rot = flag_aic = 0;
    }

    model = strmod->models[imod];
    cova = model->getCova(icov);

    // Load the parameters of the current model / structure

    if (imod_mem != imod || icov_mem != icov)
    {
      ranges = cova->getRanges();
      angles = cova->getAnisoAngles();
    }
    imod_mem = imod;
    icov_mem = icov;

    switch (icons.toEnum())
    {
      case EConsElem::E_SILL:
        ipos = ivar * (ivar + 1) / 2 + jvar;
        tritab[ipos] = param[ntot];
        flag_aic = 1;
        break;

      case EConsElem::E_PARAM:
        cova->setParam(param[ntot]);
        break;

      case EConsElem::E_RANGE:
        if (ivar == 0) ut_vector_fill(ranges, param[ntot]);
        if (ivar < ndim) ranges[ivar] = param[ntot];
        break;

      case EConsElem::E_ANGLE:
        if (ivar < ndim) angles[ivar] = param[ntot];
        flag_rot = 1;
        break;

      case EConsElem::E_T_RANGE:
        tape = model->getModTrans().getTape();
        tape->setRange(param[ntot]);
        break;

      default:
        break;
    }
  }

  // If previous parameter was a rotation, set the rotation matrix

  if (imod_mem >= 0 && icov_mem >= 0)
  {
    cova = strmod->models[imod_mem]->getCova(icov_mem);
    if (optvar.getAuthAniso())
      cova->setRanges(ranges);
    else
      cova->setRange(ranges[0]);
    if (flag_rot) cova->setAnisoAngles(angles);
    if (flag_aic)
    {
      VectorDouble sill = matrix_produit_lu_VD(nvar, tritab.data());
      MatrixSquareGeneral mat(nvar);
      mat.setValues(sill);
      cova->setSill(mat);
    }
    flag_rot = flag_aic = 0;
  }

  // If 'lock_samerot' is ON, copy the rotation to all structures

  if (strmod->optvar.getLockSamerot())
  {
    for (imod = 0; imod < strmod->nmodel; imod++)
    {
      model = strmod->models[imod];

      // Look for the first basic structure with a rotation defined

      for (icov = 0, found = -1; icov < model->getCovaNumber() && found < 0;
          icov++)
      {
        if (model->getCova(icov)->hasRange()) found = icov;
      }
      if (found < 0) continue;
      cova1 = model->getCova(found);

      for (icov = 1; icov < model->getCovaNumber(); icov++)
      {
        if (icov == found) continue;
        cova = model->getCova(icov);
        if (!cova->getAnisoRotMatVec().empty())
          cova->setAnisoAngles(cova1->getAnisoAngles());
      }
    }
  }

  // Update the internal function (which may have changed during the fit 

  for (imod = 0; imod < strmod->nmodel; imod++)
    model_setup(strmod->models[imod]);
  return;
}

/****************************************************************************/
/*!
 **  Check if the basic structure can be compressed
 **
 ** \return  1 if it can be removed; 0 otherwise
 **
 ** \param[in]  strmod      StrMod structure
 ** \param[in]  imod        Rank of the Model
 ** \param[in]  icov        Rank of the basic structure
 ** \param[in]  hmax        Maximum distance
 ** \param[in]  gmax        Maximum variance
 ** \param[in]  tolsigma    Percentage of the variance
 **
 *****************************************************************************/
static int st_structure_reduce(StrMod *strmod,
                               int imod,
                               int icov,
                               double hmax,
                               double gmax,
                               double tolsigma)
{
  Model *model = strmod->models[imod];
  int nvar = model->getVariableNumber();
  int ndim = model->getDimensionNumber();
  VectorDouble d1(ndim, hmax);
  VectorDouble tab(nvar * nvar);
  CovCalcMode mode;

  mode.update(0, 0, ECalcMember::LHS, icov, 0, 0);
  mode.setOrderVario(STRMOD->norder);
  model_calcul_cov(model, mode, 1, 1., d1, tab.data());

  for (int ivar = 0; ivar < nvar; ivar++)
  {
    if (tab[ivar + nvar * ivar] > tolsigma * gmax / 100.) return (0);
  }

  return (1);
}

/*****************************************************************************/
/*!
 **  Calculates the values of the model corresponding
 **  to the experimentations stored in strexp
 **
 ** \param[in]  nbexp    Number of experimental conditions
 ** \param[in]  imod     Rank of the model
 ** \param[in]  strexps  Array of StrExp structures to be freed
 ** \param[in]  strmod   StrMod structure
 **
 ** \param[out] tabge    Array of generic covariance values
 **
 *****************************************************************************/
static void st_evaluate_vario(int imod,
                              int nbexp,
                              std::vector<StrExp> &strexps,
                              StrMod *strmod,
                              VectorDouble &tabge)
{
  Model *model = strmod->models[imod];
  int nvar = model->getVariableNumber();
  int ndim = strmod->models[0]->getDimensionNumber();
  VectorDouble d0(ndim);
  VectorDouble tab(nvar * nvar);
  CovCalcMode mode;

  mode.update(0, 0, ECalcMember::LHS, -1, 0, 0);
  mode.setOrderVario(strmod->norder);

  /* Loop on the experimental conditions */

  for (int i = 0; i < nbexp; i++)
  {
    int ivar = strexps[i].ivar;
    int jvar = strexps[i].jvar;

    for (int idim = 0; idim < ndim; idim++)
      d0[idim] = strexps[i].dd[idim];
    tabge[i] = model_calcul_cov_ij(model, mode, ivar, jvar, d0);
  }
  return;
}

/*****************************************************************************/
/*!
 **  Calculates the values of the model corresponding
 **  to the experimentations stored in vmap
 **
 ** \param[in]  imod     Rank of the model
 ** \param[in]  strmod   StrMod structure
 **
 ** \param[out] tabge    Array of generic covariance values
 **
 *****************************************************************************/
static void st_evaluate_vmap(int imod, StrMod *strmod, VectorDouble &tabge)
{
  Model *model = strmod->models[imod];
  int ndim = model->getDimensionNumber();
  int nvar = model->getVariableNumber();
  int nech = DBMAP->getSampleNumber();
  VectorDouble d0(ndim);
  VectorDouble tab(nvar * nvar);
  db_index_sample_to_grid(DBMAP, nech / 2, INDG1);

  CovCalcMode mode;
  mode.update(0, 0, ECalcMember::LHS, -1, 0, 0);
  mode.setOrderVario(strmod->norder);

  /* Loop on the experimental conditions */

  int ecr = 0;
  for (int iech = 0; iech < nech; iech++)
  {
    db_index_sample_to_grid(DBMAP, iech, INDG2);
    for (int idim = 0; idim < ndim; idim++)
      d0[idim] = (INDG2[idim] - INDG1[idim]) * DBMAP->getDX(idim);

    int ijvar = 0;
    for (int ivar = 0; ivar < nvar; ivar++)
      for (int jvar = 0; jvar <= ivar; jvar++, ijvar++)
      {
        if (FFFF(DBMAP->getVariable(iech, ijvar))) continue;
        tabge[ecr++] = model_calcul_cov_ij(model, mode, ivar, jvar, d0);
      }
  }
  return;
}

/****************************************************************************/
/*!
 **  Find the rank of the parid which matches several criteria
 **
 ** \return  Rank of the parid or -1 if not found
 **
 ** \param[in]  strmod      StrMod structure
 ** \param[in]  npar        Number of parid
 ** \param[in]  imod0       Rank of the target model (or -1)
 ** \param[in]  icov0       Rank of the target covariance (or -1)
 ** \param[in]  icons0      Target type of the constraint (EConsElem)
 ** \param[in]  ivar0       Target first variable (or -1)
 ** \param[in]  jvar0       Target second variable (or -1)
 **
 *****************************************************************************/
static int st_parid_match(StrMod *strmod,
                          int npar,
                          int imod0,
                          int icov0,
                          const EConsElem &icons0,
                          int ivar0,
                          int jvar0)
{
  int ntot, imod, icov, ivar, jvar;
  EConsElem icons;

  for (ntot = 0; ntot < npar; ntot++)
  {
    st_parid_decode(strmod->parid[ntot], &imod, &icov, &icons, &ivar, &jvar);

    /* Check the answer */

    if (imod0 >= 0 && imod != imod0) continue;
    if (icov0 >= 0 && icov != icov0) continue;
    if (icons != icons0) continue;
    if (ivar0 >= 0 && ivar != ivar0) continue;
    if (jvar0 >= 0 && jvar != jvar0) continue;
    return (ntot);
  }
  return (-1);
}

/****************************************************************************/
/*!
 **  Check if the resulting Model is definite positive
 **
 ** \return  1 if the model is not definite positive
 **
 ** \param[in]  model  Model structure
 **
 *****************************************************************************/
static int st_check_definite_positive(Model *model)
{
  int nvar = model->getVariableNumber();
  VectorDouble cc(nvar * nvar);
  VectorDouble vecpro(nvar * nvar);
  VectorDouble valpro(nvar);

  /* Loop on the basic structures */

  for (int icov = 0; icov < model->getCovaNumber(); icov++)
  {
    message("\nCheck the Sill Matrix for structure %d\n", icov + 1);

    /* Load the matrix of sills */

    for (int ivar = 0; ivar < nvar; ivar++)
      for (int jvar = 0; jvar < nvar; jvar++)
        CC(ivar,jvar)= model->getSill(icov,ivar,jvar);

        /* Check definite positiveness */

    if (!is_matrix_definite_positive(nvar, cc.data(), valpro.data(),
                                     vecpro.data(), 1)) return 1;
  }

  return (0);
}

/****************************************************************************/
/*!
 **  Cut negative eigen values
 **
 ** \return  1 if all eigen values are strictly positive
 **
 ** \param[in]  nvar   Number of variables
 ** \param[in]  icov0  Rank of the target basic structure
 ** \param[in]  matcor Matrix of sills
 **
 ** \param[out]  matcoru Matrix of sills (can coincide with input matcor)
 **
 *****************************************************************************/
static int st_truncate_negative_eigen(int nvar,
                                      int icov0,
                                      VectorDouble &matcor,
                                      VectorDouble &matcoru)

{
  VectorDouble cc(nvar * nvar);
  VectorDouble vecpro(nvar * nvar);
  VectorDouble valpro(nvar);

  /* Load into temporary arrays */

  for (int ivar = 0; ivar < nvar; ivar++)
    for (int jvar = 0; jvar < nvar; jvar++)
      CC(ivar,jvar)= MATCOR(icov0,ivar,jvar);

  if (matrix_eigen(cc.data(), nvar, valpro.data(), vecpro.data()))
    messageAbort("st_truncate_negative_eigen");

  /* Check positiveness */

  int flag_positive = 1;
  for (int ivar = 0; ivar < nvar; ivar++)
    if (VALPRO(ivar)<= 0) flag_positive = 0;
  if (!flag_positive)
  {

    /* Calculate the new coregionalization matrix */

    for (int ivar = 0; ivar < nvar; ivar++)
      for (int jvar = 0; jvar < nvar; jvar++)
      {
        double sum = 0.;
        for (int kvar = 0; kvar < nvar; kvar++)
          sum += MAX(VALPRO(kvar),0.) * VECPRO(ivar, kvar) * VECPRO(jvar, kvar);
        MATCORU(icov0,ivar,jvar)= sum;
      }
    }
    else
    {
      for (int ivar=0; ivar<nvar; ivar++)
      for (int jvar=0; jvar<nvar; jvar++)
      MATCORU(icov0,ivar,jvar) = MATCOR(icov0,ivar,jvar);
    }

  return (flag_positive);
}

/****************************************************************************/
/*!
 **  Sum the sill of all the stuctures
 **
 ** \return The sum of the sills over all the sill matrices for a given variable
 **
 ** \param[in]  ivar0     Index of the variable
 ** \param[in]  ncova     Number of basic structures
 ** \param[in]  nvar      Number of variables
 ** \param[in]  alpha     The coregionalisation matrices (Dim: nvar^2 x ncova )
 **
 **
 *****************************************************************************/
static double st_sum_sills(int ivar0, int ncova, int nvar, VectorDouble &alpha)
{
  int icov;
  double Sr;

  /* Initializations */

  Sr = 0;

  for (icov = 0; icov < ncova; icov++)
    Sr += ALPHA(icov, ivar0, ivar0);

  return Sr;
}

/****************************************************************************/
/*!
 **  Evaluate the score
 **
 ** \return  Value for the score
 **
 ** \param[in]  nvar        Number of variables
 ** \param[in]  nvs2        Dimension of the triangular matrix (Dimension: nvar)
 ** \param[in]  ncova       Number of basic structures
 ** \param[in]  npadir      Total number of lags
 ** \param[in]  wt          Array of weights attached to variogram lags
 ** \param[in]  gg          Array of experimental values
 ** \param[in]  ge          Array of generic covariance values
 ** \param[in]  matcor      Matrix of sills
 **
 *****************************************************************************/
static double st_score(int nvar,
                       int nvs2,
                       int ncova,
                       int npadir,
                       VectorDouble &wt,
                       VectorDouble &gg,
                       VectorDouble &ge,
                       VectorDouble &matcor)
{
  double score, dd, coeff;
  int ivar, jvar, ijvar, icov, ipadir;

  /* Loop on the variables */

  score = 0.;
  for (ivar = ijvar = 0; ivar < nvar; ivar++)
    for (jvar = 0; jvar <= ivar; jvar++, ijvar++)
    {
      coeff = (ivar == jvar) ? 1 :
                               2;
      for (ipadir = 0; ipadir < npadir; ipadir++)
      {
        dd = GG(ijvar, ipadir);
        if (FFFF(dd)) continue;
        for (icov = 0; icov < ncova; icov++)
          dd -= MATCOR(icov,ivar,jvar)* GE(icov,ijvar,ipadir);
        score += coeff * WT(ijvar, ipadir) * dd * dd;
      }
    }

  return (score);
}

/****************************************************************************/
/*!
 ** Return the rank of the variable pair, given the two variable ranks
 **
 ** \param[in] ivar0    Index of the first variable
 ** \param[in] jvar0    Index of the second variable
 **
 *****************************************************************************/
static int st_combineVariables(int ivar0, int jvar0)
{
  if (ivar0 > jvar0)
    return (ivar0 + jvar0 * (jvar0 + 1) / 2);
  else
    return (jvar0 + ivar0 * (ivar0 + 1) / 2);
}

/****************************************************************************/
/*! 
 ** Minimization of the order-4 polynomial
 **
 ** \return Value for which the fourth order polynomial is minimum.
 **
 ** \param[in] icov0    Index of the target basic structure
 ** \param[in] ivar0    Index of the target variable
 ** \param[in] ncova    Number of basic structures
 ** \param[in] nvar     Number of variables
 ** \param[in] npadir   Maximum number of lags for all directions
 ** \param[in] xrmax    Maximum value for the solution
 ** \param[in] xr       Current vector of sqrt(constraint/(sum of the sills))
 ** \param[in] alpha    Current auxiliary matrices alpha
 ** \param[in] wt       Array of weights attached to variogram lags
 ** \param[in] gg       Array of experimental values
 ** \param[in] ge       Array of generic covariance values
 ** \param[in] consSill Vector of constant required sills (optional)
 **
 *****************************************************************************/
static double st_minimize_P4(int icov0,
                             int ivar0,
                             int ncova,
                             int nvar,
                             int npadir,
                             double xrmax,
                             VectorDouble &xr,
                             VectorDouble &alpha,
                             VectorDouble &wt,
                             VectorDouble &gg,
                             VectorDouble &ge,
                             const VectorDouble &consSill)
{
  int number, nin, ivar, k, irr, irl, jcov, l, icov, nvs2;
  double retval, value;
  double a, c, d, s, x[3];

  /* Core allocation */

  VectorDouble Nir_v(nvar);
  VectorDouble Mrr_v(npadir);
  VectorDouble Crr_v(npadir);
  VectorDouble Airk_v(npadir * nvar);
  VectorDouble Birk_v(npadir * nvar);
  VectorDouble xx(2);
  VectorDouble xt(2);
  VectorDouble xest(2);

  irr = st_combineVariables(ivar0, ivar0);
  nvs2 = nvar * (nvar + 1) / 2;

  for (ivar = 0; ivar < nvar; ivar++)
  {
    irl = st_combineVariables(ivar0, ivar);
    Nir_v[ivar] = 0.;
    for (k = 0; k < npadir; k++)
      Nir_v[ivar] += WT(irl,k)* GE(icov0,0,k) * GE(icov0,0,k);
    }

  for (k = 0; k < npadir; k++)
  {
    Mrr_v[k] = 0.;
    for (jcov = 0; jcov < ncova; jcov++)
    {
      if (jcov == icov0) continue;
      Mrr_v[k] += ALPHA(jcov,ivar0,ivar0)* (GE(jcov,0,k) - GE(icov0,0,k));
    }
  }

  for (k = 0; k < npadir; k++)
    for (ivar = 0; ivar < nvar; ivar++)
    {
      irl = st_combineVariables(ivar0, ivar);
      value = 0.;
      for (l = 0; l < npadir; l++)
        value += WT(irl,l)* GG(irl,l) * GE(icov0,0,l);
      value *= GE(icov0,0,k)/ Nir_v[ivar];
      AIRKV(k,ivar)= GG(irl,k) - value;
    }

  for (k = 0; k < npadir; k++)
    for (ivar = 0; ivar < nvar; ivar++)
    {
      irl = st_combineVariables(ivar0, ivar);
      BIRKV(k,ivar)= 0.;
      for (icov = 0; icov < ncova; icov++)
      {
        if (icov == icov0) continue;
        value = 0.;
        for (l = 0; l < npadir; l++)
          value += WT(irl,l)* GE(icov,0,l) * GE(icov,0,l);
        BIRKV(k,ivar)+= ALPHA(icov,ivar0,ivar) *
        (GE(icov,0,k) - value * GE(icov0,0,k))/Nir_v[ivar];
      }
    }

  for (k = 0; k < npadir; k++)
    Crr_v[k] = GG(irr,k)- consSill[ivar0] * GE(icov0,0,k);

  a = 0.;
  for (k = 0; k < npadir; k++)
    a += WT(irr,k)* Mrr_v[k] * Mrr_v[k];

  c = 0.;
  for (k = 0; k < npadir; k++)
    for (ivar = 0; ivar < nvar; ivar++)
    {
      if (ivar != ivar0)
      {
        irl = st_combineVariables(ivar0, ivar);
        s = xr[ivar] * BIRKV(k, ivar);
        c += WT(irl,k)* s * s;
      }
      else
      {
        c -= WT(irr,k) * Mrr_v[k] * Crr_v[k];
      }
    }
  c *= 2.;

  d = 0.;
  for (k = 0; k < npadir; k++)
    for (ivar = 0; ivar < nvar; ivar++)
    {
      if (ivar == ivar0) continue;
      irl = st_combineVariables(ivar0, ivar);
      d -= WT(irl,k)* AIRKV(k,ivar) * xr[ivar] * BIRKV(k,ivar);
    }

  d *= 4.;

  number = solve_P3(4. * a, 0., 2. * c, d, x);

  switch (number)
  {
    case 1:
      retval = MAX(0., MIN(xrmax, x[0]));
      break;

    case 3:
    {
      nin = 0;
      xx[0] = x[0];
      xx[1] = x[2];
      for (k = 0; k < 2; k++)
      {
        xt[k] = MAX(0., MIN(xrmax, xx[k]));
        xest[k] = (a * xt[k] * xt[k] * xt[k] * xt[k] + c * xt[k] * xt[k]
                   + d * xt[k])
                  / 2.;
        if (xt[k] == xx[k]) nin++;
      }
      if (nin == 1)
      {
        retval = (xt[0] == xx[0]) ? xx[0] :
                                    xx[1];
      }
      else
      {
        retval = (xest[0] < xest[1]) ? xt[0] :
                                       xt[1];
      }
    }
      break;

    default:
      retval = xr[ivar0];
      //mes_abort("st_minimize_P4: Number of extrema is impossible");
      break;
  }

  if (retval < 1.e-9) retval = 0.;

  return retval;
}

/****************************************************************************/
/*!
 ** Update alpha from xr (diagonal only)
 **
 ** \param[in]     icov0    Target basic structure
 ** \param[in]     ivar0    Target variable
 ** \param[in]     ncova    Number of basic structures
 ** \param[in]     nvar     Number of variables
 ** \param[in]     xr       Current vector of sqrt(constraint/(sum of the sills))
 ** \param[in,out] alpha    Current auxiliary matrices alpha
 ** \param[in]     consSill Vector of constant Sill (optional)
 **
 *****************************************************************************/
void st_updateAlphaDiag(int icov0,
                        int ivar0,
                        int ncova,
                        int nvar,
                        VectorDouble &xr,
                        VectorDouble &alpha,
                        const VectorDouble &consSill)
{
  double srm, value;

  srm = st_sum_sills(ivar0, ncova, nvar, alpha) - ALPHA(icov0, ivar0, ivar0);
  value = consSill[ivar0] / (xr[ivar0] * xr[ivar0]) - srm;
  ALPHA(icov0,ivar0,ivar0)= MAX(0.,value);
}

/*****************************************************************************/
/*!
 ** Update 'sills' for the structures other than the current one
 **
 ** \param[in]  icov0    Target basic structure
 ** \param[in]  ivar0    Target variable
 ** \param[in]  ncova    Number of basic structures
 ** \param[in]  nvar     Number of variables
 ** \param[in]  alpha    Current auxiliary matrices alpha
 ** \param[in]  xr       Current vector of sqrt(constraint/(sum of the sills))
 ** \param[out] matcor   Current sills matrices
 **
 ******************************************************************************/
static void st_updateOtherSills(int icov0,
                                int ivar0,
                                int ncova,
                                int nvar,
                                VectorDouble &alpha,
                                VectorDouble &xr,
                                VectorDouble &matcor)
{
  int jcov, ivar;

  for (jcov = 0; jcov < ncova; jcov++)
  {
    if (jcov == icov0) continue;
    for (ivar = 0; ivar < nvar; ivar++)
      MATCOR(jcov,ivar0, ivar)= MATCOR(jcov,ivar, ivar0) =
      ALPHA(jcov,ivar0, ivar) * xr[ivar0] * xr[ivar];
    }
  }

  /*****************************************************************************/
  /*!
   ** Update the sill matrix for the current structure
   ** (except diagonal in the constrained case)
   **
   ** \param[in]  icov0    Target basic structure
   ** \param[in]  ivar0    Target variable
   ** \param[in]  npadir   Maximum number of lags for all directions
   ** \param[in]  nvar     Number of variables
   ** \param[in]  ncova    Number of basic structures
   ** \param[in]  wt       Array of weights attached to variogram lags
   ** \param[in]  ge       Array of generic covariance values
   ** \param[in]  gg       Array of experimental values
   ** \param[in]  consSill Vector of optional constrained sills
   ** \param[out] matcor   Matrices of sills
   **
   *****************************************************************************/
static void st_updateCurrentSillGoulard(int icov0,
                                        int ivar0,
                                        int npadir,
                                        int nvar,
                                        int ncova,
                                        VectorDouble &wt,
                                        VectorDouble &ge,
                                        VectorDouble &gg,
                                        const VectorDouble &consSill,
                                        VectorDouble &matcor)
{
  int ivar, ilagdir, icov, ivs2, nvs2;
  double tot1, tot2, wtloc, ggloc, geloc;

  VectorDouble mv(npadir);
  nvs2 = nvar * (nvar + 1) / 2;

  // Loop on the variables

  for (ivar = 0; ivar < nvar; ivar++)
  {
    if (ivar == ivar0 && !FFFF(consSill[ivar0])) continue;

    for (ilagdir = 0; ilagdir < npadir; ilagdir++)
    {
      mv[ilagdir] = 0.;
      for (icov = 0; icov < ncova; icov++)
      {
        if (icov == icov0) continue;
        mv[ilagdir] += MATCOR(icov,ivar0,ivar)* GE(icov,0,ilagdir);
      }
    }

    tot1 = 0.;
    tot2 = 0.;
    ivs2 = st_combineVariables(ivar0, ivar);

    for (ilagdir = 0; ilagdir < npadir; ilagdir++)
    {
      wtloc = WT(ivs2, ilagdir);
      ggloc = GG(ivs2, ilagdir);
      geloc = GE(icov0, 0, ilagdir);

      if (!FFFF(ggloc))
      {
        tot1 += wtloc * geloc * (ggloc - mv[ilagdir]);
        tot2 += wtloc * geloc * geloc;
      }
    }

    MATCOR(icov0,ivar0, ivar)= MATCOR(icov0,ivar, ivar0) = tot1 / tot2;
  }

  return;
}

/*****************************************************************************/
/*!
 ** Update 'sills' for the structures other than the current one
 **
 ** \param[in] icov0    Target basic structure
 ** \param[in] ivar0    Index of the variable
 ** \param[in] nvar     Number of variables
 ** \param[in] xr       Current vector of sqrt(constraint/(sum of the sills))
 ** \param[in] consSill Vector of optional constant sills
 ** \param[in] matcor   Matrices of sills
 ** \param[out] alpha   Current auxiliary matrices alpha
 **
 ******************************************************************************/
static void st_updateAlphaNoDiag(int icov0,
                                 int ivar0,
                                 int nvar,
                                 VectorDouble &xr,
                                 const VectorDouble &consSill,
                                 VectorDouble &matcor,
                                 VectorDouble &alpha)
{
  int ivar;
  double value;
  for (ivar = 0; ivar < nvar; ivar++)
  {
    if (ivar == ivar0 && !FFFF(consSill[ivar0])) continue;
    value = MATCOR(icov0,ivar0,ivar)/ (xr[ivar0] * xr[ivar]);
    ALPHA(icov0,ivar0,ivar)= ALPHA(icov0,ivar,ivar0) = value;
  }
}

/****************************************************************************/
/*!
 * Update the sill matrix for the current structure 'icov0' (diagonal only)
 *
 ** \param[in] icov0    Target basic structure
 ** \param[in] ivar0    Index of the variable
 ** \param[in] nvar     Number of variables
 ** \param[in] alpha    Current auxiliary matrices alpha  
 ** \param[in] xr       Current vector of sqrt(constraint/(sum of the sills))
 **
 ** \param[out] matcor   Matrices of sills     
 **
 ****************************************************************************/
static void st_updateCurrentSillDiag(int icov0,
                                     int ivar0,
                                     int nvar,
                                     VectorDouble &alpha,
                                     VectorDouble &xr,
                                     VectorDouble &matcor)
{
  double value;
  value = xr[ivar0] * xr[ivar0] * ALPHA(icov0, ivar0, ivar0);
  if (value < 0.) value = 0.;
  MATCOR(icov0,ivar0,ivar0)= value;
}

/*****************************************************************************/
/*!
 ** Make sure the current matrix of sills is definite positive
 ** (diagonal unchanged)
 **
 ** \param[in]      icov0    Index of the target basic structure
 ** \param[in]      nvar     Number of variables
 ** \param[in]      consSill Vector of optional constant sills
 ** \param[in,out]  matcor   Matrices of sills
 **
 *****************************************************************************/
static int st_makeDefinitePositive(int icov0,
                                   int nvar,
                                   const VectorDouble &consSill,
                                   VectorDouble &matcor)
{
  VectorDouble muold(nvar);
  VectorDouble norme1(nvar);

  for (int ivar = 0; ivar < nvar; ivar++)
    muold[ivar] = MATCOR(icov0, ivar, ivar);

  int flag_positive = st_truncate_negative_eigen(nvar, icov0, matcor, matcor);

  if (flag_positive) goto label_end;

  for (int ivar = 0; ivar < nvar; ivar++)
  {
    if (FFFF(consSill[ivar]))
      norme1[ivar] = 1.;
    else
    {
      if (ABS(MATCOR(icov0,ivar, ivar)) > EpsFit)
        norme1[ivar] = sqrt(muold[ivar] / MATCOR(icov0, ivar, ivar));
      else
        norme1[ivar] = (ABS(muold[ivar]) < EpsFit) ? 1. :
                                                     0.;
    }
  }

  for (int ivar = 0; ivar < nvar; ivar++)
    for (int jvar = 0; jvar < nvar; jvar++)
      MATCOR(icov0,ivar, jvar)*= norme1[ivar] * norme1[jvar];

  label_end: return flag_positive;
}

/*****************************************************************************/
/*!
 **  Optimization under constraints
 **
 ** \return The number of iteration required to reach convergence
 **
 ** \param[in]  nvar        Number of variables
 ** \param[in]  nvs2        Dimension of the triangular matrix (Dimension: nvar)
 ** \param[in]  ncova       Number of basic structures
 ** \param[in]  npadir      Total number of lags
 ** \param[in]  mauto       Option_AutoFit structure
 ** \param[in]  wt          Array of weights attached to variogram lags
 ** \param[in]  gg          Array of experimental values
 ** \param[in]  ge          Array of generic covariance values
 ** \param[in]  matcor      Matrix of sills
 **
 ** \param[out] score       The convergence score
 **
 *****************************************************************************/
static int st_optimize_under_constraints(int nvar,
                                         int nvs2,
                                         int ncova,
                                         int npadir,
                                         const Option_AutoFit &mauto,
                                         VectorDouble &wt,
                                         VectorDouble &gg,
                                         VectorDouble &ge,
                                         VectorDouble &matcor,
                                         double *score)
{
  int iter, icov, icov0, ivar, ivar0, jcov, jvar;
  double score_old, score_new, xrmax, denom;

  /* Core allocation */

  VectorDouble xr(nvar);
  VectorDouble alpha(nvar * nvar * ncova);

  /* Calculate the initial score */

  score_new = st_score(nvar, nvs2, ncova, npadir, wt, gg, ge, matcor);

  for (icov = 0; icov < ncova; icov++)
    for (ivar = 0; ivar < nvar; ivar++)
      for (jvar = 0; jvar < nvar; jvar++)
        ALPHA(icov,ivar,jvar)= MATCOR(icov,ivar,jvar);

        /***********************/
        /* Iterative procedure */
        /***********************/

        /* Optional printout */

  st_goulard_debug_title(nvar, ncova);
  st_goulard_debug_current(nvar, ncova, 0, matcor, TEST);

  for (ivar = 0; ivar < nvar; ivar++)
  {
    if (FFFF(mauto.getConstantSills(ivar)))
      xr[ivar] = 1.;
    else
      xr[ivar] = sqrt(
          mauto.getConstantSills(ivar) / st_sum_sills(ivar, ncova, nvar,
                                                      alpha));
  }

  for (iter = 0; iter < mauto.getMaxiter(); iter++)
  {
    for (icov0 = 0; icov0 < ncova; icov0++)
    {
      for (ivar0 = 0; ivar0 < nvar; ivar0++)
      {
        denom = st_sum_sills(ivar0, ncova, nvar, alpha)-
        ALPHA(icov0,ivar0,ivar0);
        if (!FFFF(mauto.getConstantSills(ivar0)) && denom < 1e-30) continue;
        if (!FFFF(mauto.getConstantSills(ivar0)))
        {
          xrmax = sqrt(mauto.getConstantSills(ivar0) / denom);
          xr[ivar0] = st_minimize_P4(icov0, ivar0, ncova, nvar, npadir, xrmax,
                                     xr, alpha, wt, gg, ge,
                                     mauto.getConstantSills());

          if (xr[ivar0] == 0)
          {
            xr[ivar0] = 1.;
            for (jcov = 0; jcov < ncova; jcov++)
            {
              if (jcov == icov0) continue;
              for (jvar = 0; jvar < nvar; jvar++)
                ALPHA(jcov,ivar0,jvar)= ALPHA(jcov,jvar,ivar0) = 0.;
              }
            }
          }

          /* Update 'alpha' (diagonal only) */

        if (!FFFF(mauto.getConstantSills(ivar0)))
          st_updateAlphaDiag(icov0, ivar0, ncova, nvar, xr, alpha,
                             mauto.getConstantSills());

        /* Update 'sills' for the structures other than the current one */

        st_updateOtherSills(icov0, ivar0, ncova, nvar, alpha, xr, matcor);

        /* Update the sill matrix for the current structure */
        /* (except diagonal in the constrained case)        */

        st_updateCurrentSillGoulard(icov0, ivar0, npadir, nvar, ncova, wt, ge,
                                    gg, mauto.getConstantSills(), matcor);

        /* Update sill matrix for the current structure (for diagonal) */

        if (!FFFF(mauto.getConstantSills(ivar0)))
          st_updateCurrentSillDiag(icov0, ivar0, nvar, alpha, xr, matcor);

        /* Make sure the current matrix of sills if definite positive */
        /* (diagonal unchanged)                                       */

        (void) st_makeDefinitePositive(icov0, nvar, mauto.getConstantSills(),
                                       matcor);

        /* Update 'alpha' for the current structure */

        for (ivar = 0; ivar < nvar; ivar++)
          st_updateAlphaNoDiag(icov0, ivar, nvar, xr, mauto.getConstantSills(),
                               matcor, alpha);
      }
    }

    /* Update the score */

    score_old = score_new;
    score_new = st_score(nvar, nvs2, ncova, npadir, wt, gg, ge, matcor);
    if (ABS(score_new - score_old) / score_old < mauto.getTolred()) break;
  }

  /* Optional printout */

  st_goulard_debug_current(nvar, ncova, iter, matcor, score_new);
  *score = score_new;
  return (iter);
}

/****************************************************************************/
/*!
 **  Initialize the system for Goulard algorithm
 **
 ** \return Error return code
 **
 ** \param[in]  nvar        Number of variables
 ** \param[in]  nvs2        Dimension of the triangular matrix (Dimension: nvar)
 ** \param[in]  ncova       Number of basic structures
 ** \param[in]  npadir      Total number of lags
 ** \param[in]  wt          Array of weights attached to variogram lags
 ** \param[in]  gg          Array of experimental values
 ** \param[in]  ge          Array of generic covariance values
 ** \param[in]  consSill    Array containing the constraints (optional)
 **
 ** \param[out] matcor     Matrix of sills (Dimension: ncova * nvar * nvar)
 **
 *****************************************************************************/
static int st_initialize_goulard(int nvar,
                                 int nvs2,
                                 int ncova,
                                 int npadir,
                                 VectorDouble &wt,
                                 VectorDouble &gg,
                                 VectorDouble &ge,
                                 const VectorDouble &consSill,
                                 VectorDouble &matcor)
{
  int ivar, jvar, ijvar, icov, jcov, ipadir, error, nae, nai;
  double wtloc, ggloc, geloc1, geloc2, temp, be;

  /* Initializations */

  error = 1;

  /* Core allocation */

  VectorDouble aa(ncova * ncova);
  VectorDouble bb(ncova);
  VectorDouble res(ncova);
  VectorDouble Ae(ncova);
  VectorDouble Ai(ncova * ncova);
  VectorDouble bi(ncova);

  /* Initialize the constraints matrices */

  for (icov = 0; icov < ncova; icov++)
  {
    Ae[icov] = 1.;
    bi[icov] = 0.;
    for (jcov = 0; jcov < ncova; jcov++)
      Ai[icov * ncova + jcov] = (icov == jcov);
  }

  /* Loop on the variables */

  for (ivar = ijvar = 0; ivar < nvar; ivar++)
    for (jvar = 0; jvar <= ivar; jvar++, ijvar++)
    {
      /* Reset the arrays */

      for (icov = 0; icov < ncova; icov++)
      {
        bb[icov] = 0.;
        for (jcov = 0; jcov < ncova; jcov++)
          AA(icov,jcov)= 0.;
        }

        for (ipadir=0; ipadir<npadir; ipadir++)
        {
          wtloc = WT(ijvar,ipadir);
          ggloc = GG(ijvar,ipadir);
          if (FFFF(wtloc) || FFFF(ggloc)) continue;
          for (icov=0; icov<ncova; icov++)
          {
            geloc1 = GE(icov,ijvar,ipadir);
            if (FFFF(geloc1)) continue;
            bb[icov] += wtloc * ggloc * geloc1;
            for (jcov=0; jcov<ncova; jcov++)
            {
              geloc2 = GE(jcov,ijvar,ipadir);
              if (FFFF(geloc2)) continue;
              AA(icov,jcov) += wtloc * geloc1 * geloc2;
            }
          }
        }

        nae = 0;
        nai = 0;

        if (ivar == jvar && ! consSill.empty() && !FFFF(consSill[ivar]))
        {
          be = consSill[ivar];
          nae = 1;
          nai = ncova;
        }

        /* Update (taking into account possible constraints) */

        if (matrix_qoci(ncova,aa.data(),bb.data(),nae,Ae.data(),&be,nai,
                Ai.data(),bi.data(),res.data()))
        {
          for (icov=0; icov<ncova; icov++)
          res[icov] = (ivar == jvar);
          if (ivar == jvar && ! consSill.empty() && !FFFF(consSill[ivar]))
          {
            temp = consSill[ivar] / ncova;
            for(icov=0; icov<ncova; icov++)
            res[icov] = temp;
          }
        }

        /* Store in the output matrix */

        for (icov=0; icov<ncova; icov++)
        MATCOR(icov,ivar,jvar) = MATCOR(icov,jvar,ivar) = res[icov];
      }

      /* Set the error return code */

  error = 0;
  return (error);
}

/****************************************************************************/
/*!
 **  Copy the resulting sill into the Model
 **
 ** \param[in]  nvar        Number of variables
 ** \param[in]  ncova       Number of covariances
 ** \param[in]  sill        Array of sills
 **
 ** \param[out] model       Fitted Model structure
 **
 *****************************************************************************/
static void st_goulard_sill_to_model(int nvar,
                                     int ncova,
                                     VectorDouble &sill,
                                     Model *model)
{
  int ijvar, nvs2;

  nvs2 = nvar * (nvar + 1) / 2;
  for (int icov = 0; icov < ncova; icov++)
    for (int ivar = ijvar = 0; ivar < nvar; ivar++)
      for (int jvar = 0; jvar <= ivar; jvar++, ijvar++)
      {
        model->setSill(icov, ivar, jvar, SILL(icov, ijvar));
        model->setSill(icov, jvar, ivar, SILL(icov, ijvar));
      }
}

/****************************************************************************/
/*!
 **  Internal function for Goulard under constraints
 **
 ** \return  Error returned code
 **
 ** \param[in]  mauto       Option_AutoFit structure
 ** \param[in]  nvar        Number of variables
 ** \param[in]  ncova       Number of covariances
 ** \param[in]  npadir      Maximum number of lags for all directions
 ** \param[in]  wt          Array of weights (Dimension: npadir)
 ** \param[in]  gg          Array of experimental values (Dimension: npadir)
 ** \param[in]  ge          Array of model values (Dimension: npadir)
 **
 ** \param[out] sill        Array of resulting sills
 **
 *****************************************************************************/
static int st_goulard_with_constraints(const Option_AutoFit &mauto,
                                       int nvar,
                                       int ncova,
                                       int npadir,
                                       VectorDouble &wt,
                                       VectorDouble &gg,
                                       VectorDouble &ge,
                                       VectorDouble &sill)
{
  double crit;
  int nvs2, icov, ivar, jvar, ijvar, flag_positive, iter;

  /* Initializations */

  nvs2 = nvar * (nvar + 1) / 2;

  /* Core allocation */

  VectorDouble matcor(nvar * nvar * ncova);

  /* Initialize the Goulard system */

  st_initialize_goulard(nvar, nvs2, ncova, npadir, wt, gg, ge,
                        mauto.getConstantSills(), matcor);

  /* Update according to the eigen values */

  flag_positive = 1;
  for (icov = 0; icov < ncova; icov++)
    if (!st_makeDefinitePositive(icov, nvar, mauto.getConstantSills(), matcor))
      flag_positive = 0;

  if (!flag_positive)
  {
    for (icov = 0; icov < ncova; icov++)
      for (ivar = 0; ivar < nvar; ivar++)
      {
        if (!FFFF(mauto.getConstantSills(ivar)))
          MATCOR(icov,ivar,ivar)= mauto.getConstantSills(ivar) / ncova;
          else
          MATCOR(icov,ivar,ivar)=1.;

          for(jvar=0;jvar<ivar;jvar++)
          MATCOR(icov,ivar,jvar) = MATCOR(icov,jvar,ivar)=0.;
        }

        /* Perform the optimization under constraints */

        iter = st_optimize_under_constraints(nvar,nvs2,ncova,npadir,mauto,
            wt,gg,ge,matcor,&crit);

        /* Optional printout */

        st_goulard_score(mauto,1,ncova,iter,crit);
      }

      /* Load the parameters in the final model */

  for (icov = 0; icov < ncova; icov++)
    for (ivar = ijvar = 0; ivar < nvar; ivar++)
      for (jvar = 0; jvar < nvar; jvar++, ijvar++)
        SILL(icov,ijvar)= MATCOR(icov,ivar,jvar);

  return (0);
}

/****************************************************************************/
/*!
 **  Routine for fitting a model using an experimental variogram
 **  for the Intrinsic case
 **
 ** \return  Error return code
 **
 ** \param[in]  model       Model to be fitted
 ** \param[in]  mauto       Option_AutoFit structure
 ** \param[in]  npadir      Maximum number of lags for all directions
 ** \param[in]  wt          Array of weights (Dimension: npadir)
 ** \param[in]  gg          Array of experimental values (Dimension: npadir)
 ** \param[in]  ge          Array of model values (Dimension: npadir)
 ** \param[in]  wt2         Array of weights (Dimension: npadir)
 ** \param[in]  ge1         Array of model values (Dimension: npadir)
 ** \param[in]  ge2         Array of model values (Dimension: npadir)
 ** \param[in]  gg2         Array of experimental values (Dimension: npadir)
 **
 ** \param[out] alphau      Array
 ** \param[out] sill1       Array of resulting sills
 **
 ** \remark  Internal arrays:
 ** \remark  MP : Contains the current Model (ijvar,ipadir)
 **
 *****************************************************************************/
static int st_sill_fitting_int(Model *model,
                               const Option_AutoFit& mauto,
                               int npadir,
                               VectorDouble &wt,
                               VectorDouble &gg,
                               VectorDouble &ge,
                               VectorDouble &wt2,
                               VectorDouble &ge1,
                               VectorDouble &ge2,
                               VectorDouble &gg2,
                               VectorDouble &alphau,
                               VectorDouble &sill1)
{
  double sum, pivot, newval, crit, crit_mem;
  int error, icov, nvs2, ivar, jvar, ijvar, ipadir, nvar, ncova, iter;

  /* Initializations */

  error = 1;
  nvar = model->getVariableNumber();
  ncova = model->getCovaNumber();
  nvs2 = nvar * (nvar + 1) / 2;
  crit_mem = 1.e30;
  for (icov = 0; icov < ncova; icov++)
    alphau[icov] = 1. / (double) ncova;

  /* Iterative procedure */

  Option_AutoFit mauto_new(mauto);
  mauto_new.setMaxiter(1);
  for (iter = 0; iter < mauto.getMaxiter(); iter++)
  {

    /* Initialize the arrays for the first pass */

    for (ipadir = 0; ipadir < npadir; ipadir++)
    {
      for (ivar = ijvar = 0; ivar < nvar; ivar++)
      {
        for (jvar = 0; jvar <= ivar; jvar++, ijvar++)
        {
          sum = 0.;
          for (icov = 0; icov < ncova; icov++)
            sum += alphau[icov] * GE(icov, ijvar, ipadir);
          GE1(ijvar,ipadir)= sum;
        }
      }
    }

    /* Call Goulard with 1 structure (no constraint) */

    st_sill_reset(nvar,1,sill1);
    if (st_goulard_without_constraint(mauto_new,nvar,1,npadir,
            wt,gg,ge1,sill1,&crit)) goto label_end;

    /* Initialize the arrays for the second pass */

    for (ivar=ijvar=0; ivar<nvar; ivar++)
    for (jvar=0; jvar<=ivar; jvar++, ijvar++)
    {
      pivot = sill1[ijvar];
      for (ipadir=0; ipadir<npadir; ipadir++)
      {
        GG2(ijvar,ipadir) = (pivot == 0) ? 0. : GG(ijvar,ipadir) / pivot;
        WT2(ijvar,ipadir) = WT(ijvar,ipadir) * pivot * pivot;
        for (icov=0; icov<ncova; icov++)
        GE2(icov,ijvar,ipadir) = GE(icov,ijvar,ipadir);
      }
    }

    /* Call Goulard with 1 variable (no constraint) */

    if (st_goulard_without_constraint(mauto_new,1,ncova,npadir*nvs2,
            wt2,gg2,ge2,alphau,&crit)) goto label_end;

    /* Stopping criterion */

    if (ABS(crit) < mauto_new.getTolred() ||
        ABS(crit-crit_mem) / ABS(crit) < mauto_new.getTolred()) break;
    crit_mem = crit;
  }

  /* Patch the final model */

  for (icov = 0; icov < ncova; icov++)
    for (ivar = ijvar = 0; ivar < nvar; ivar++)
      for (jvar = 0; jvar <= ivar; jvar++, ijvar++)
      {
        newval = alphau[icov] * sill1[ijvar];
        model->setSill(icov, ivar, jvar, newval);
        model->setSill(icov, jvar, ivar, newval);
      }

  /* Error return code */

  error = 0;

  label_end:

  return (error);
}

/****************************************************************************/
/*!
 **  General Routine for fitting a model using an experimental variogram
 **
 ** \return  Error return code
 **
 ** \param[in]  flag_reset  1 to reset the array of sill before usage
 ** \param[in]  flag_title  1 to print the title
 ** \param[in]  model       Model to be fitted
 ** \param[in]  mauto       Option_AutoFit structure
 **
 *****************************************************************************/
static int st_goulard_fitting(int flag_reset,
                              int flag_title,
                              Model *model,
                              const Option_AutoFit& mauto)
{
  int status;
  double crit;

  /* Initialize the array of sills */

  if (flag_reset)
    st_sill_reset(model->getVariableNumber(), model->getCovaNumber(),
                  RECINT.sill);

  /* Print the debug title (optional) */

  if (flag_title)
    st_goulard_debug_title(model->getVariableNumber(), model->getCovaNumber());

  /* Dispatch */

  if (!mauto.getFlagIntrinsic())
  {

    /* No intrinsic hypothesis */

    if (FFFF(mauto.getConstantSillValue()))
    {
      /* Without constraint on the sill */

      status = st_goulard_without_constraint(mauto, model->getVariableNumber(),
                                             model->getCovaNumber(),
                                             RECINT.npadir, RECINT.wt,
                                             RECINT.gg, RECINT.ge, RECINT.sill,
                                             &crit);
    }
    else
    {

      /* With constraint on the sill */

      status = st_goulard_with_constraints(mauto, model->getVariableNumber(),
                                           model->getCovaNumber(),
                                           RECINT.npadir, RECINT.wt, RECINT.gg,
                                           RECINT.ge, RECINT.sill);
    }

    /* Copy the array 'sill' in the Model */

    st_goulard_sill_to_model(model->getVariableNumber(), model->getCovaNumber(),
                             RECINT.sill, model);
  }
  else
  {
    status = st_sill_fitting_int(model, mauto, RECINT.npadir, RECINT.wt,
                                 RECINT.gg, RECINT.ge, RECINT.wt2, RECINT.ge1,
                                 RECINT.ge2, RECINT.gg2, RECINT.alphau,
                                 RECINT.sill1);
  }

  return (status);
}

/****************************************************************************/
/*!
 **  Check if the model has (at least) one intrinsic structure
 **
 ** \return 1 if model contains at least one intrinsic structure; 0 otherwise
 **
 ** \param[in]  model    Model structure containing the basic structures
 ** \param[in]  filter   Array specifying if a basic structure is filtered (1)
 **                      or not (0). This array is optional
 **
 *****************************************************************************/
static int st_model_has_intrinsic(Model *model, int *filter)
{
  int flag_range, flag_param, flag_aniso, flag_rotation, icov, n_int;
  int min_order, max_ndim, flag_int_1d, flag_int_2d;
  double scalfac, parmax;

  /* Loop on the basic structures */

  n_int = 0;
  for (icov = 0; icov < model->getCovaNumber(); icov++)
  {
    if (filter != nullptr && filter[icov]) continue;
    model_cova_characteristics(model->getCovaType(icov), cov_name, &flag_range,
                               &flag_param, &min_order, &max_ndim, &flag_int_1d,
                               &flag_int_2d, &flag_aniso, &flag_rotation,
                               &scalfac, &parmax);
    if (min_order == 0) n_int++;
  }
  return (n_int > 0);
}

/****************************************************************************/
/*!
 **  Reduce the Model by discarding unnecessary basic structure(s)
 **
 ** \return  1 if the Model has been reduced; 0 otherwise
 **
 ** \param[in]  strmod      StrMod structure
 ** \param[in]  npar        Number of parameters
 ** \param[in]  hmax        Maximum distance value
 ** \param[in]  gmax        Maximum variogram value
 ** \param[in]  param       Current values for parameters
 ** \param[in]  lower       Lower values for parameters
 ** \param[in]  upper       Upper values for parameters
 ** \param[in]  mauto       Option_AutoFit structure
 **
 *****************************************************************************/
static int st_model_auto_strmod_reduce(StrMod *strmod,
                                       int *npar,
                                       double hmax,
                                       double gmax,
                                       VectorDouble &param,
                                       VectorDouble &lower,
                                       VectorDouble &upper,
                                       Option_AutoFit &mauto)
{
  int ntot, nparloc, icov, jcov, ncova, ivar, jvar, jmod, kcov, ncovleft;
  int flag_modified, imod, nmodel;
  int flag_range, flag_param, min_order, max_ndim, flag_int_1d;
  int flag_int_2d, flag_aniso, flag_rotation, rank;
  int lost_rank, lost_imod, lost_icov;
  double scalfac, parmax;
  VectorInt flag_compress;
  Option_VarioFit optvar;
  Model *model;
  EConsElem icons;

  /* Initializations */

  nparloc = *npar;
  optvar = strmod->optvar;

  /* Load the parameters in the Model */

  st_model_auto_strmod_define(strmod, nparloc, param);

  /* Run the last Goulard algorithm (if necessary) */

  st_goulard_verbose(0, mauto);
  if (optvar.getFlagGoulardUsed()) for (imod = 0; imod < strmod->nmodel; imod++)
  {
    ST_PREPAR_GOULARD(imod);
    (void) st_goulard_fitting(1, 1, STRMOD->models[imod], mauto);
  }
  st_goulard_verbose(1, mauto);

  /* Initializations */

  if (optvar.getFlagNoreduce()) return (0);
  nmodel = strmod->nmodel;
  ncova = 0;
  for (imod = 0; imod < nmodel; imod++)
    ncova = MAX(ncova, strmod->models[imod]->getCovaNumber());
  if (ncova <= 1) return (0);
  flag_modified = 0;
  flag_compress.resize(ncova * nmodel);

  /* Check if the basic structure must be discarded */

  lost_rank = lost_imod = lost_icov = -1;
  ncovleft = 0;
  for (imod = 0; imod < strmod->nmodel; imod++)
  {
    model = strmod->models[imod];
    for (icov = 0; icov < model->getCovaNumber(); icov++)
    {
      FLAG_COMPRESS(imod,icov) = st_structure_reduce(strmod, imod, icov, hmax,
                                                     gmax, mauto.getTolsigma());

      // Check that at least one intrinsic structure is kept.
      // Otherwise, do not suppress the current structure
      if (optvar.getKeepIntstr() && FLAG_COMPRESS(imod, icov))
      {
        if (!st_model_has_intrinsic(model, &FLAG_COMPRESS(imod, 0)))
        FLAG_COMPRESS(imod,icov) = 0;
      }
      if (FLAG_COMPRESS(imod, icov))
      {

        flag_modified++;
        if (mauto.getVerbose() > 0 || debug_query("converge"))
        {
          if (flag_modified == 1)
            mestitle(0, "Suppressing the unnecessary basic structures");
          message("Structure '%s' in model #%d is suppressed\n",
                  model->getCovName(icov).c_str(), imod + 1);
        }

        /* Check if the structure to be deleted contains the only rotation */
        /* (this may be TRUE only if lock_samerot is ON */

        if (optvar.getLockSamerot())
        {
          rank = st_parid_match(strmod, nparloc, imod, icov, EConsElem::ANGLE,
                                -1, -1);
          if (rank >= 0 && lost_rank < 0)
          {
            lost_rank = rank;
            st_parid_decode(strmod->parid[lost_rank], &lost_imod, &lost_icov,
                            &icons, &ivar, &jvar);
            if (mauto.getVerbose() > 0 || debug_query("converge"))
            {
              message("Note: This structure contains rotation parameters.\n");
              message("As the fitting method considers a shared rotation\n");
              message("This rotation must be swapped to another structure\n");
            }
          }
        }
      }
      else
        ncovleft++;
    }
  }

  if (ncovleft <= 0)
  {
    message(
        "Due to the tolerance (%lf (percent)), no basic structure would be left\n",
        mauto.getTolsigma());
    message("No structure is discarded\n");
    flag_modified = 0;
    goto label_end;
  }

  /* Loop on the basic structures */

  for (ntot = 0; ntot < nparloc; ntot++)
  {
    st_parid_decode(strmod->parid[ntot], &imod, &icov, &icons, &ivar, &jvar);
    if (imod == lost_imod && icov == lost_icov && icons == EConsElem::ANGLE)
      continue;
    if (FLAG_COMPRESS(imod, icov)) param[ntot] = TEST;
  }

  /* If the rotation is "lost", set it elsewhere */

  if (lost_rank >= 0)
  {
    for (imod = 0; imod < strmod->nmodel; imod++)
    {
      model = strmod->models[imod];
      for (icov = 0; icov < model->getCovaNumber(); icov++)
      {
        if (FLAG_COMPRESS(imod, icov)) continue;
        model_cova_characteristics(model->getCovaType(icov), cov_name,
                                   &flag_range, &flag_param, &min_order,
                                   &max_ndim, &flag_int_1d, &flag_int_2d,
                                   &flag_aniso, &flag_rotation, &scalfac,
                                   &parmax);
        if (! DEFINE_ANIROT) continue;

        /* This non-masked component can be assigned the lost rotation */

        if (mauto.getVerbose() > 0 || debug_query("converge"))
        {
          message("The Rotation is swapped to Structure '%s' in model #%d\n",
                  model->getCovName(icov).c_str(), imod + 1);
        }

        ivar = 0;
        while (1)
        {
          rank = st_parid_match(strmod, nparloc, lost_imod, lost_icov,
                                EConsElem::ANGLE, ivar, -1);
          if (rank < 0) goto label_compress;
          strmod->parid[rank] = st_parid_encode(imod, icov, EConsElem::ANGLE,
                                                ivar, 0);
          ivar++;
        }
      }
    }
  }

  /* Compress the vector of parameters and bounds */

  label_compress: nparloc = st_compress_parid(ntot, strmod->parid, param, lower,
                                              upper);

  /* Modifying the covariance ranks in parid */

  for (imod = 0; imod < strmod->nmodel; imod++)
    for (icov = jcov = 0; icov < strmod->models[imod]->getCovaNumber(); icov++)
    {
      if (FLAG_COMPRESS(imod, icov)) continue;

      /* Shift icov parameter */

      for (ntot = 0; ntot < nparloc; ntot++)
      {
        st_parid_decode(strmod->parid[ntot], &jmod, &kcov, &icons, &ivar,
                        &jvar);
        if (jmod == imod && kcov == icov)
          strmod->parid[ntot] = st_parid_encode(jmod, jcov, icons, ivar, jvar);
      }
      jcov++;
    }

  // Suppress the basic structures from the model (
  // Warning: We start from the end in order to avoid having to compress
  // FLAG_COMPRESS consequently

  for (imod = strmod->nmodel - 1; imod >= 0; imod--)
  {
    ncova = strmod->models[imod]->getCovaNumber();
    for (icov = ncova - 1; icov >= 0; icov--)
    {
      if (!FLAG_COMPRESS(imod, icov)) continue;
      strmod->models[imod]->delCova(icov);
    }
  }

  label_end: *npar = nparloc;
  return (flag_modified);
}

/****************************************************************************/
/*!
 **  Define the options of the model
 **
 ** \return Error return code
 **
 ** \param[in]  model    Model structure containing the basic structures
 ** \param[in]  optvar   Opt_Vario structure
 **
 *****************************************************************************/
static int st_model_define(Model *model, const Option_VarioFit &optvar)
{
  int flag_range, flag_param, flag_aniso, flag_rotation, jcov;
  int min_order, max_ndim, flag_int_1d, flag_int_2d;
  double scalfac, parmax;
  CovAniso *cova;

  /* Loop on the basic structures */

  for (jcov = 0; jcov < model->getCovaNumber(); jcov++)
  {
    cova = model->getCova(jcov);
    model_cova_characteristics(cova->getType(), cov_name, &flag_range,
                               &flag_param, &min_order, &max_ndim, &flag_int_1d,
                               &flag_int_2d, &flag_aniso, &flag_rotation,
                               &scalfac, &parmax);

    /// TODO [Cova] : To be restored? Why this?
    /*
     cova->setFlagAniso((DEFINE_ANICOEF && flag_aniso));
     // Rotation is not suppressed, even if the rotation flag is
     // switched off, as a rotation angle may still be defined
     cova->setFlagRotation(flag_aniso);
     */
  }

  /* Check that the model contains at least one intrinsic structure */

  if (optvar.getKeepIntstr())
  {
    if (!st_model_has_intrinsic(model, nullptr))
    {
      messerr("Automatic Fitting must keep one Intrinsic Basic Structure");
      messerr("No such structure is provided");
      return (1);
    }
  }
  return (0);
}

/****************************************************************************/
/*!
 **  Define the options of the structure Opt_Vario
 **
 ** \return Error return code
 **
 ** \param[in]  vario       Vario structure containing the exp. variogram
 ** \param[in]  model       Model structure
 ** \param[in]  constraints Constraints structure
 **
 ** \param[out]  optvar  Opt_Vario structure
 **
 *****************************************************************************/
static int st_alter_model_optvar(const Vario *vario,
                                 Model *model,
                                 Constraints &constraints,
                                 Option_VarioFit &optvar)
{
  int ndim, ndir, idir, n_2d, n_3d;

  /* Initializations */

  ndim = model->getDimensionNumber();
  ndir = vario->getDirectionNumber();
  n_2d = n_3d = 0;

  /* Calculate the number of directions in 2-D and 3-D */

  if (ndim == 2)
  {
    n_2d = ndir;
    n_3d = 0;
  }

  /* 3-D case */

  if (ndim == 3)
  {
    for (idir = 0; idir < ndir; idir++)
    {
      if (vario->getCodir(idir, 2) == 0)
        n_2d++;
      else
        n_3d++;
    }
    optvar.setLockNo3d(n_3d <= 0);
    optvar.setLockIso2d(n_2d <= 0);
  }

  /* Clever setting of options */

  if (ndir <= ndim) optvar.setAuthRotation(0);
  if (ndir <= 1 || ndim <= 1) optvar.setAuthAniso(0);
  if (ndir <= 1 || ndim <= 1) optvar.setAuthRotation(0);

  if (n_3d <= 0) optvar.setLockNo3d(1);
  if (n_2d <= 1) optvar.setLockIso2d(1);
  if (optvar.getLockIso2d()) optvar.setAuthRotation(0);
  if (optvar.getLockNo3d()) optvar.setLockRot2d(1);

  /* Consequences of no anisotropy */

  if (!optvar.getAuthAniso())
  {
    optvar.setAuthRotation(0);
    optvar.setLockSamerot(0);
    optvar.setLockRot2d(0);
    optvar.setLockNo3d(0);
    optvar.setLockIso2d(0);
  }

  /* Case when properties are defined: Goulard is switch off */

  if (model->getModTransMode() == EModelProperty::ANAM && model->getModTrans().getAnam()->getType()
      != EAnam::HERMITIAN
      && optvar.getFlagGoulardUsed())
  {
    message("Goulard option is switched OFF");
    message("due to presence of ANAM Properties (type != EAnam::HERMITIAN)\n");
    optvar.setFlagGoulardUsed(0);
  }

  /* Case when constraints involve sill(s) */

  if (constraints.isDefinedForSill() && optvar.getFlagGoulardUsed())
  {
    message(
        "Goulard option is switched OFF due to presence of Sills in Constraints\n");
    if (modify_constraints_on_sill(constraints)) return (1);
    optvar.setFlagGoulardUsed(0);
  }

  /* Return an error if Goulard is not used in multivariate case */

  if (model->getVariableNumber() > 1 && !optvar.getFlagGoulardUsed())
  {
    messerr("In Multivariate case, Goulard option is mandatory");
    messerr("It seems that it has been switched OFF. This is an error");
    return (1);
  }
  return (0);
}

/****************************************************************************/
/*!
 **  Define the options of the structure Opt_Vario
 **
 ** \param[in]  dbmap       Db Grid structure containing the Vmap
 ** \param[in]  model       Model structure
 ** \param[in]  constraints Constraints structure
 **
 ** \param[out]  optvar  Opt_Vario structure
 **
 *****************************************************************************/
static int st_alter_vmap_optvar(const Db *dbmap,
                                Model *model,
                                Constraints &constraints,
                                Option_VarioFit &optvar)
{
  /* Clever setting of options */

  optvar.setAuthAniso(1);
  optvar.setAuthRotation(1);
  optvar.setLockNo3d(dbmap->getNDim() <= 2);

  /* Case when properties are defined: Goulard is switch off */

  if (model->getModTransMode() == EModelProperty::ANAM && model->getModTrans().getAnam()->getType()
      != EAnam::HERMITIAN
      && optvar.getFlagGoulardUsed())
  {
    message("Goulard option is switched OFF");
    message("due to presence of ANAM Properties (type != EAnam::HERMITIAN)\n");
    optvar.setFlagGoulardUsed(0);
  }

  /* Case when constraints involve sill(s) */

  if (constraints.isDefinedForSill() && optvar.getFlagGoulardUsed())
  {
    message(
        "Goulard option is switched OFF due to presence of Sills in Constraints\n");
    optvar.setFlagGoulardUsed(0);
    if (modify_constraints_on_sill(constraints)) return (1);
  }

  /* Return an error if Goulard is not used in multivariate case */

  if (model->getVariableNumber() > 1 && !optvar.getFlagGoulardUsed())
  {
    messerr("In Multivariate case, Goulard option is mandatory");
    messerr("It seems that it has been switched OFF. This is an error");
    return (1);
  }

  return (0);
}

/****************************************************************************/
/*!
 **  Count the number of parameters
 **
 ** \return  Number of parameters
 ** \return  -1 if one of the model components is not valid or
 ** \return  -2 if Goulard is not used althrough in multivariate
 **
 ** \param[in]  vario       Vario structure containing the exp. variogram
 ** \param[in]  model1      Model first structure
 ** \param[in]  model2      Model second structure
 ** \param[in]  constraints Constraints structure
 ** \param[in]  optvar      Opt_Vario structure
 **
 ** \param[out] param Array giving the default parameter value
 ** \param[out] lower Array giving the minimum parameter value
 ** \param[out] upper Array giving the maximum parameter value
 **
 ** \remarks The arrays param, lower and upper are allocated by this function
 ** \remarks (Dimension: Number of Parameters)
 ** \remarks They should be freed by the user
 **
 *****************************************************************************/
static int st_model_auto_count(const Vario *vario,
                               Model *model1,
                               Model *model2,
                               Constraints &constraints,
                               Option_VarioFit &optvar,
                               VectorDouble &param,
                               VectorDouble &lower,
                               VectorDouble &upper)
{
  int ntot, jcov, idim, ndim, flag_range, flag_param, flag_aniso, flag_rotation;
  int min_order, max_ndim, flag_int_1d, flag_int_2d, first_covrot, nvar;
  int imod, iparam;
  double scalfac, parmax;
  Model *model;

  /* Loop on the two candidate models */

  ntot = 0;
  for (imod = 0; imod < 2; imod++)
  {
    model = (imod == 0) ? model1 :
                          model2;
    if (model == nullptr) continue;

    /* Initializations */

    ndim = model->getDimensionNumber();
    nvar = model->getVariableNumber();
    first_covrot = -1;
    if (st_alter_model_optvar(vario, model, constraints, optvar)) return (-2);

    /* Define the model */

    if (st_model_define(model, optvar)) return (-1);

    /* Count the number of parameters */

    for (jcov = 0; jcov < model->getCovaNumber(); jcov++)
    {
      model_cova_characteristics(model->getCovaType(jcov), cov_name,
                                 &flag_range, &flag_param, &min_order,
                                 &max_ndim, &flag_int_1d, &flag_int_2d,
                                 &flag_aniso, &flag_rotation, &scalfac,
                                 &parmax);
      if (max_ndim > 0 && ndim > max_ndim)
      {
        messerr("The structure '%s' is limited to dimension (%d)", cov_name,
                max_ndim);
        messerr("The current study is carried out in dimension (%d)", ndim);
        return (-1);
      }

      /* AIC coefficient -> Sill */
      if (DEFINE_AIC) ntot += nvar * (nvar + 1) / 2;

      /* Third parameter */
      if (DEFINE_THIRD) ntot++;

      /* Range or scale factor */
      if (DEFINE_RANGE) ntot++;

      /* Anisotropy coefficients */
      if (DEFINE_ANICOEF)
      {
        if (ndim == 2)
          ntot++;
        else if (ndim == 3)
        {
          if (!optvar.getLockIso2d()) ntot++;
          if (!optvar.getLockNo3d()) ntot++;
        }
        else
        {
          for (idim = 1; idim < ndim; idim++)
            ntot++;
        }
      }

      /* Anisotropy angles */
      if (DEFINE_ANIROT)
      {
        if (ndim == 2)
        {
          if (TAKE_ROT)
          {
            first_covrot = jcov;
            ntot++;
          }
        }
        else if (ndim == 3 && optvar.getLockRot2d())
        {
          if (TAKE_ROT)
          {
            first_covrot = jcov;
            ntot++;
          }
        }
        else
        {
          if (TAKE_ROT)
          {
            first_covrot = jcov;
            for (idim = 0; idim < ndim; idim++)
              ntot++;
          }
        }
      }

      /* Tapering range */

      if (DEFINE_T_RANGE) ntot++;
    }
  }

  /* Allocation */

  param.resize(ntot);
  lower.resize(ntot);
  upper.resize(ntot);
  for (iparam = 0; iparam < ntot; iparam++)
    param[iparam] = lower[iparam] = upper[iparam] = TEST;

  return (ntot);
}

/****************************************************************************/
/*!
 **  Evaluate the model for an experiment
 **
 ** \return  Value of the model
 **
 ** \param[in]  nbexp        Number of experimental conditions
 ** \param[in]  npar         Number of parameters
 ** \param[in]  param        Current values for parameters
 **
 ** \param[out] tabge        Array of resulting values
 **
 *****************************************************************************/
static void st_strmod_vario_evaluate(int nbexp,
                                     int npar,
                                     VectorDouble &param,
                                     VectorDouble &tabge)
{
  int imod;

  /* Define the current values of the parameters */

  st_model_auto_strmod_define(STRMOD, npar, param);

  /* Run the Goulard algorithm (if necessary) */

  st_goulard_verbose(0, MAUTO);
  if (STRMOD->optvar.getFlagGoulardUsed())
    for (imod = 0; imod < STRMOD->nmodel; imod++)
    {
      ST_PREPAR_GOULARD(imod);
      (void) st_goulard_fitting(1, 0, STRMOD->models[imod], MAUTO);
    }
  st_goulard_verbose(1, MAUTO);

  /* Calculate the array of model values */

  st_evaluate_vario(0, nbexp, STREXPS, STRMOD, tabge);

  return;
}

/****************************************************************************/
/*!
 **  Prepare the array for Goulard's algorithm
 **  in the case of VarioMap calculation
 **
 ** \param[in]  imod      Rank of the model
 **
 *****************************************************************************/
static void st_prepar_goulard_vmap(int imod)

{
  Model *model = STRMOD->models[imod];
  VectorDouble &ge = RECINT.ge;
  int npadir = RECINT.npadir;
  int ndim = model->getDimensionNumber();
  int nvar = model->getVariableNumber();
  int ncova = model->getCovaNumber();
  int nvs2 = nvar * (nvar + 1) / 2;
  int nech = DBMAP->getSampleNumber();
  VectorDouble d0(ndim);
  VectorDouble tab(nvar * nvar);
  db_index_sample_to_grid(DBMAP, nech / 2, INDG1);
  CovCalcMode mode;

  /* Loop on the basic structures */

  for (int icov = 0; icov < ncova; icov++)
  {
    mode.update(0, ITEST, ECalcMember::LHS, icov, 0, 0);
    mode.setOrderVario(STRMOD->norder);

    /* Loop on the experiments */

    for (int ipadir = 0; ipadir < RECINT.npadir; ipadir++)
    {
      db_index_sample_to_grid(DBMAP, ipadir, INDG2);
      for (int idim = 0; idim < ndim; idim++)
        d0[idim] = (INDG2[idim] - INDG1[idim]) * DBMAP->getDX(idim);
      model_calcul_cov(model, mode, 1, 1., d0, tab.data());

      /* Loop on the variables */

      int ijvar = 0;
      for (int ivar = 0; ivar < nvar; ivar++)
        for (int jvar = 0; jvar <= ivar; jvar++, ijvar++)
          GE(icov,ijvar,ipadir)= tab[ijvar];
        }
      }
  return;
}

/****************************************************************************/
/*!
 **  Evaluate the model for an experiment
 **
 ** \return  Value of the model
 **
 ** \param[in]  npar         Number of parameters
 ** \param[in]  param        Current values for parameters
 **
 ** \param[out] tabge        Array of resulting values
 **
 *****************************************************************************/
static void st_strmod_vmap_evaluate(int /*nbexp*/,
                                    int npar,
                                    VectorDouble &param,
                                    VectorDouble &tabge)
{
  int imod;

  /* Define the current values of the parameters */

  st_model_auto_strmod_define(STRMOD, npar, param);

  /* Run the Goulard algorithm (if necessary) */

  st_goulard_verbose(0, MAUTO);
  if (STRMOD->optvar.getFlagGoulardUsed())
    for (imod = 0; imod < STRMOD->nmodel; imod++)
    {
      ST_PREPAR_GOULARD(imod);
      (void) st_goulard_fitting(1, 0, STRMOD->models[imod], MAUTO);
    }
  st_goulard_verbose(1, MAUTO);

  /* Calculate the array of model values */

  for (imod = 0; imod < STRMOD->nmodel; imod++)
    st_evaluate_vmap(imod, STRMOD, tabge);

  return;
}

/****************************************************************************/
/*!
 **  Define Cholesky decomposition of the variance-covariance matrix
 **
 ** \param[in]  vario    Vario structure
 ** \param[in]  model    Model structure
 **
 ** \param[out] varchol  Cholesky array
 **
 ** \remark  In the case of Automatic Model Fitting with Properties,
 ** \remark  we consider that the variance (of the Variogram) cannot
 ** \remark  serve as initial value for the sill of the Model
 **
 *****************************************************************************/
static void st_vario_varchol_manage(const Vario *vario,
                                    Model *model,
                                    VectorDouble &varchol)
{
  int nvar, size, nvar2, i, ivar, jvar;
  VectorDouble aux;
  Model *model_nugget;
  CovCalcMode mode;

  /* Initializations */

  nvar = vario->getVariableNumber();
  size = nvar * (nvar + 1) / 2;
  nvar2 = nvar * nvar;
  model_nugget = nullptr;

  /* Allocation */

  aux.resize(nvar2);
  varchol.resize(size);

  /* Particular case of Properties attached to the Model */

  if (model->getModTransMode() != EModelProperty::NONE)
  {
    model_nugget = model_default(model->getDimensionNumber(),
                                 model->getVariableNumber());
    model_calcul_cov(model, mode, 1, 1., VectorDouble(), aux.data());
    for (i = 0; i < nvar2; i++)
      aux[i] = vario->getVarIJ(i) / aux[i];
  }
  else
  {
    for (i = 0; i < nvar2; i++)
      aux[i] = vario->getVarIJ(i);
  }

  if (matrix_cholesky_decompose(aux.data(), varchol.data(), nvar))
  {

    /* The matrix is filled arbitrarily */
    for (ivar = 0; ivar < nvar; ivar++)
      for (jvar = 0; jvar < nvar; jvar++)
        AUX(ivar,jvar)= (ivar == jvar);
        (void) matrix_cholesky_decompose(aux.data(),varchol.data(),nvar);
      }
  model_nugget = model_free(model_nugget);
  return;
}

/****************************************************************************/
/*!
 **  Define Cholesky decomposition of the variance-covariance matrix
 **
 ** \param[in]  dbmap    Db structure
 **
 ** \param[out] varchol  Cholesky array
 **
 *****************************************************************************/
static void st_vmap_varchol_manage(const Db *dbmap, VectorDouble &varchol)
{
  int nvar, size, nvar2, i, iloc, ivar;
  double mini, maxi, gmax;
  VectorDouble aux;

  /* Initializations */

  nvar = dbmap->getVariableNumber();
  size = nvar * (nvar + 1) / 2;
  nvar2 = nvar * nvar;

  /* Allocation */

  aux.resize(nvar2);
  for (i = 0; i < nvar2; i++)
    aux[i] = 0.;
  for (ivar = 0; ivar < nvar; ivar++)
  {
    iloc = db_attribute_identify(dbmap, ELoc::Z, ivar);
    (void) db_attribute_range(dbmap, iloc, &mini, &maxi, &gmax);
    AUX(ivar,ivar)= gmax;
  }
  varchol.resize(size);
  if (matrix_cholesky_decompose(aux.data(), varchol.data(), nvar))
    messageAbort("Error in the Cholesky decomposition of the variance matrix");
  return;
}

/****************************************************************************/
/*!
 **  Update of the resulting model after the automatic fit
 **
 ** \param[in,out]  strmod          StrMod structure
 ** \param[in]      optvar          Opt_Vario structure
 **
 *****************************************************************************/
static void st_model_post_update(StrMod *strmod, const Option_VarioFit &optvar)
{
  for (int imod = 0; imod < strmod->nmodel; imod++)
  {
    Model *model = strmod->models[imod];

    for (int icov = 0; icov < model->getCovaNumber(); icov++)
    {
      CovAniso *cova = model->getCova(icov);
      if (!cova->hasRange()) continue;
      if (cova->getAnisoCoeffs().empty()) continue;

      // Set the isotropy (by application)

      if (!optvar.getAuthAniso())
      {
        if (!cova->isIsotrop())
        my_throw("Posterior Check: The covariance should be isotropic");
      }
    }
  }
}

/****************************************************************************/
/*!
 **  Manage memory for variogram fitting
 **
 ** \return  Error returned code
 **
 ** \param[in]  mauto     Option_AutoFit structure
 ** \param[in]  flag_exp  1 for experimental variogram
 ** \param[in]  ndim      Space dimension
 ** \param[in]  nvar      Number of variables
 ** \param[in]  nbexp     Number of experimental variogram values
 ** \param[in]  ncova     Number of covariances
 ** \param[in]  npadir    Total number of lags
 **
 *****************************************************************************/
static int st_manage_recint(const Option_AutoFit &mauto,
                            int flag_exp,
                            int ndim,
                            int nvar,
                            int nbexp,
                            int ncova,
                            int npadir)
{
  int nvs2, nv2;

  /* Initializations */

  nv2 = nvar * nvar;
  nvs2 = nvar * (nvar + 1) / 2;

  RECINT.npadir = npadir;
  RECINT.wt.resize(npadir * nvs2);
  st_blank(RECINT.wt, npadir * nvs2);
  RECINT.gg.resize(npadir * nvs2);
  st_blank(RECINT.gg, npadir * nvs2);
  RECINT.ge.resize(npadir * nvs2 * ncova);
  st_blank(RECINT.ge, npadir * nvs2 * ncova);
  RECINT.sill.resize(nvs2 * ncova);
  st_blank(RECINT.sill, nvs2 * ncova);
  RECINT.covtab.resize(nv2);

  if (flag_exp)
  {
    RECINT.wtc.resize(nbexp);
    st_blank(RECINT.wtc, nbexp);
    RECINT.ggc.resize(nbexp);
    st_blank(RECINT.ggc, nbexp);
    RECINT.dd.resize(npadir * nvs2 * ndim);
    st_blank(RECINT.dd, npadir * nvs2 * ndim);
  }

  if (mauto.getFlagIntrinsic())
  {
    RECINT.alphau.resize(ncova);
    st_blank(RECINT.alphau, ncova);
    RECINT.sill1.resize(nvs2);
    st_blank(RECINT.sill1, nvs2);
    RECINT.ge1.resize(nvs2 * npadir);
    st_blank(RECINT.ge1, nvs2 * npadir);
    RECINT.ge2.resize(nvs2 * npadir * ncova);
    st_blank(RECINT.ge2, nvs2 * npadir * ncova);
    RECINT.wt2.resize(nvs2 * npadir);
    st_blank(RECINT.wt2, nvs2 * npadir);
    RECINT.gg2.resize(nvs2 * npadir);
    st_blank(RECINT.gg2, nvs2 * npadir);
  }
  return (0);
}

/****************************************************************************/
/*!
 **  Initialize the internal structure for regularization
 **
 *****************************************************************************/
static void st_regularize_init()
{
  REGULARIZE.flag_regularize = 0;
  REGULARIZE.ndim = 0;
  for (int idim = 0; idim < 3; idim++)
  {
    REGULARIZE.ndisc[idim] = 1;
    REGULARIZE.support[idim] = 0.;
  }
}

/****************************************************************************/
/*!
 **  Automatic model fitting
 **
 ** \return  Error returned code
 **
 ** \param[in]  vario       Vario structure containing the exp. variogram
 ** \param[in]  model       Model structure containing the basic structures
 ** \param[in]  verbose     Verbose flag
 ** \param[in]  mauto_arg   Option_AutoFit structure
 ** \param[in]  cons_arg    Constraints structure
 ** \param[in]  optvar_arg  Opt_Vario structure
 **
 *****************************************************************************/
int model_auto_fit(const Vario *vario,
                   Model *model,
                   bool verbose,
                   const Option_AutoFit &mauto_arg,
                   const Constraints &cons_arg,
                   const Option_VarioFit &optvar_arg)
{
  int i, error, status, nbexp, norder, npar, npadir, npar0;
  int flag_hneg, flag_gneg, flag_reduce, nvar, ncova, ndim, flag_regular;
  double c0, hmin, hmax, gmin, gmax;
  std::vector<StrExp> strexps;
  StrMod *strmod;
  VectorDouble varchol, scale, param, lower, upper;
  static int flag_check_result = 0;

  // Getting local copy of const references

  Option_AutoFit mauto = mauto_arg;
  Option_VarioFit optvar = optvar_arg;
  Constraints constraints = cons_arg;

  /* Initializations */

  nbexp = status = npadir = 0;
  strmod = nullptr;
  ncova = model->getCovaNumber();
  ndim = model->getDimensionNumber();
  nvar = vario->getVariableNumber();
  VectorDouble angles;
  st_regularize_init();
  mauto.setVerbose(verbose);

  /* Preliminary checks */

  error = 1;
  norder = 0;
  if (vario->getCalculType() == ECalcVario::GENERAL1) norder = 1;
  if (vario->getCalculType() == ECalcVario::GENERAL2) norder = 2;
  if (vario->getCalculType() == ECalcVario::GENERAL3) norder = 3;
  if (model->getDimensionNumber() > 3)
  {
    messerr("Procedure cannot be used for space dimension (%d) larger than 3",
            model->getDimensionNumber());
    goto label_end;
  }
  if (vario->getCalculType() == ECalcVario::MADOGRAM || vario->getCalculType()
      == ECalcVario::RODOGRAM
      || vario->getCalculType() == ECalcVario::GENERAL1
      || vario->getCalculType() == ECalcVario::GENERAL2
      || vario->getCalculType() == ECalcVario::GENERAL3)
  {
    messerr("Procedure is designed only for symmetric covariance");
    return (1);
  }
  if (!FFFF(mauto.getConstantSillValue()))
  {
    if (!optvar.getFlagGoulardUsed())
    {
      messerr("When Constraints on the sum of Sills are defined");
      messerr("The Goulard option must be switched ON");
      return (1);
    }
    mauto.setConstantSills(nvar);
  }

  // Define regularizing constraints (temporarily) using "keypair" mechanism

  flag_regular = (int) get_keypone("Data_Discretization", 0.);
  if (flag_regular)
  {
    REGULARIZE.flag_regularize = 1;
    REGULARIZE.ndim = ndim;
    REGULARIZE.ndisc[0] = (int) get_keypone("Data_Discretization_NX", 1.);
    REGULARIZE.ndisc[1] = (int) get_keypone("Data_Discretization_NY", 1.);
    REGULARIZE.ndisc[2] = (int) get_keypone("Data_Discretization_NZ", 1.);
    REGULARIZE.support[0] = get_keypone("Data_Discretization_DX", 0.);
    REGULARIZE.support[1] = get_keypone("Data_Discretization_DY", 0.);
    REGULARIZE.support[2] = get_keypone("Data_Discretization_DZ", 0.);
  }

  /* Free the "keypair" mechanism strings */

  st_keypair_sill(-1, model);
  st_keypair_results(-1, 0, 0, NULL, NULL);

  /* Calculate the variogram extension */

  variogram_extension(vario, 0, 0, -1, 0, 1, TEST, TEST, TEST, TEST, &flag_hneg,
                      &flag_gneg, &c0, &hmin, &hmax, &gmin, &gmax);
  angles.resize(ndim);
  (void) ut_angles_from_codir(vario->getDimensionNumber(), 1,
                              vario->getCodir(0), angles);
  st_vario_varchol_manage(vario, model, varchol);

  /* Scale the parameters in the mauto structure */

  st_mauto_rescale(nvar, varchol, mauto);

  /* Create the experimental structures */

  if (st_get_vario_dimension(vario, &nbexp, &npadir)) goto label_end;
  strexps = st_strexp_manage(nbexp, ndim);

  /* Fill the weight and experimental tabulated arrays */

  if (st_manage_recint(mauto, 1, ndim, nvar, nbexp, ncova, npadir))
    goto label_end;

  /* Generate the default values */

  npar0 = st_model_auto_count(vario, model, nullptr, constraints, optvar, param,
                              lower, upper);

  /* Create the Model structures */

  strmod = st_model_auto_strmod_alloc(model, NULL, npar0, norder, hmax, angles,
                                      optvar, &npar);
  if (strmod == nullptr) goto label_end;

  /* Load the arrays */

  st_load_wt(vario, mauto.getWmode(), npadir, RECINT.wt);
  st_compress_array(vario, npadir, RECINT.wt, RECINT.wtc);
  st_load_gg(vario, npadir, strexps, RECINT.gg);
  st_compress_array(vario, npadir, RECINT.gg, RECINT.ggc);
  st_load_ge(vario, model, npadir, RECINT.dd, RECINT.ge);

  if (npar > 0) scale.resize(npar);

  /* Set the default values and bounds */

  st_model_auto_pardef(strmod, npar, hmax, varchol, angles, param, lower,
                       upper);
  st_model_auto_scldef(strmod, npar, hmax, varchol, scale);
  st_model_auto_constraints_apply(strmod, npar, constraints, param, lower,
                                  upper);

  /* Minimization algorithm */

  STREXPS = strexps;
  STRMOD = strmod;
  MAUTO = mauto;
  ST_PREPAR_GOULARD = st_prepar_goulard_vario;
  do
  {
    st_model_auto_strmod_print(1, strmod, mauto, param, lower, upper, npar,
                               nbexp);
    if (npar > 0)
    {
      status = foxleg_f(nbexp, npar, 0, VectorDouble(), param, lower, upper,
                        scale, mauto, 0, st_strmod_vario_evaluate, RECINT.ggc,
                        RECINT.wtc);
    }
    else
    {
      status = st_goulard_fitting(1, 1, model, mauto);
    }
    if (status > 0) goto label_end;

    st_model_auto_strmod_print(0, strmod, mauto, param, lower, upper, npar,
                               nbexp);
    flag_reduce = st_model_auto_strmod_reduce(strmod, &npar, hmax, gmax, param,
                                              lower, upper, mauto);

    /* Reset the parameters to their default values */

    if (status < 0)
    {
      for (i = 0; i < npar; i++)
        param[i] = lower[i] = upper[i] = TEST;
      st_model_auto_pardef(strmod, npar, hmax, varchol, angles, param, lower,
                           upper);
    }
    st_model_auto_scldef(strmod, npar, hmax, varchol, scale);
  }
  while (flag_reduce && npar > 0);

  /* Perform the last cosmetic updates */

  st_model_post_update(strmod, optvar);

  /* Set the returned error code */

  if (mauto.getVerbose() >= 0 && status < 0)
    messerr("\nConvergence not reached after %d iterations (%d parameters)",
            mauto.getMaxiter(), npar);

  /* Check the quality of the result */

  if (flag_check_result)
  {
    if (st_check_definite_positive(model)) goto label_end;
  }

  /* Store the sills in the keypair mechanism */

  st_keypair_sill(1, model);

  /* Set the error return code */

  error = 0;

  label_end: strmod = st_model_auto_strmod_free(strmod);

  return (error);
}

/****************************************************************************/
/*!
 **  Fitting a model using an experimental variogram
 **
 ** \return  Error return code
 **
 ** \param[in]     vario     Vario structure
 ** \param[in,out] model     Model to be fitted
 ** \param[in]     mauto     Option_AutoFit structure
 **
 *****************************************************************************/
int model_fitting_sills(const Vario *vario, Model *model, const Option_AutoFit &mauto)
{
  int nvar, ncova, ndir, nbexp, npadir, ndim;
  std::vector<StrExp> strexps;

  /*******************/
  /* Initializations */
  /*******************/

  if (model == nullptr) return (1);
  if (vario == nullptr) return (1);
  ndir = vario->getDirectionNumber();
  ndim = model->getDimensionNumber();
  nvar = model->getVariableNumber();
  ncova = model->getCovaNumber();

  /* Reset the coregionalization matrix */

  if (st_get_vario_dimension(vario, &nbexp, &npadir) || nvar <= 0 || ndir <= 0
      || ncova <= 0)
  {
    messerr("The automatic fitting tool does not function as :");
    messerr("- the number of variables is zero");
    messerr("- the number of basic structures is zero");
    messerr("- the number of directions in which variograms are");
    messerr("  calculated is zero");
    return (1);
  }

  /* Core allocation */

  if (st_manage_recint(mauto, 0, ndim, nvar, 0, ncova, npadir)) return 1;

  /* Free the keypair mechanism strings */

  st_keypair_sill(-1, model);

  /* Load the arrays */

  st_load_wt(vario, mauto.getWmode(), npadir, RECINT.wt);

  st_load_gg(vario, npadir, strexps, RECINT.gg);
  st_load_ge(vario, model, npadir, RECINT.dd, RECINT.ge);

  /* Automatic Sill Fitting procedure */

  if (st_goulard_fitting(1, 1, model, mauto)) return 1;

  /* Store the sills in the keypair mechanism */

  st_keypair_sill(1, model);

  return (0);
}

/****************************************************************************/
/*!
 **  Count the number of parameters
 **
 ** \return  Number of parameters
 ** \return  -1 if one of the model components is not valid
 ** \return  -2 if Goulard is not used althrough in multivariate
 **
 ** \param[in]  dbmap       Db grid structure containing the Vmap
 ** \param[in]  model       Model structure
 ** \param[in]  constraints Constraint structure
 ** \param[in]  optvar      Opt_Vario structure
 **
 ** \param[out] param Array giving the default parameter value
 ** \param[out] lower Array giving the minimum parameter value
 ** \param[out] upper Array giving the maximum parameter value
 **
 *****************************************************************************/
static int st_vmap_auto_count(const Db *dbmap,
                              Model *model,
                              Constraints &constraints,
                              Option_VarioFit &optvar,
                              VectorDouble &param,
                              VectorDouble &lower,
                              VectorDouble &upper)
{
  int ntot, jcov, idim, ndim, flag_range, flag_param, flag_aniso, flag_rotation;
  int min_order, max_ndim, flag_int_1d, flag_int_2d, first_covrot, nvar, iparam;
  double scalfac, parmax;

  /* Initializations */

  ndim = model->getDimensionNumber();
  nvar = model->getVariableNumber();
  first_covrot = -1;

  /* Check the number of Variogram Maps */

  if (nvar * (nvar + 1) / 2 != dbmap->getVariableNumber())
  {
    messerr("The number of items in the Db Grid for Variogram maps (%d)",
            dbmap->getVariableNumber());
    messerr("is not compatible with the number of variables in the Model (%d)",
            nvar);
    return (-1);
  }

  if (st_alter_vmap_optvar(dbmap, model, constraints, optvar)) return (-2);

  /* Define the model */

  if (st_model_define(model, optvar)) return (-1);

  /* Count the number of parameters */

  for (jcov = ntot = 0; jcov < model->getCovaNumber(); jcov++)
  {
    model_cova_characteristics(model->getCovaType(jcov), cov_name, &flag_range,
                               &flag_param, &min_order, &max_ndim, &flag_int_1d,
                               &flag_int_2d, &flag_aniso, &flag_rotation,
                               &scalfac, &parmax);
    if (max_ndim > 0 && ndim > max_ndim)
    {
      messerr("The structure '%s' is limited to dimension (%d)", cov_name,
              max_ndim);
      messerr("The current study is carried out in dimension (%d)", ndim);
      return (-1);
    }

    /* AIC coefficient -> Sill */
    if (DEFINE_AIC) ntot += nvar * (nvar + 1) / 2;

    /* Third parameter */
    if (DEFINE_THIRD) ntot++;

    /* Range or scale factor */
    if (DEFINE_RANGE) ntot++;

    /* Anisotropy coefficients */
    if (DEFINE_ANICOEF)
    {
      if (ndim == 2)
        ntot++;
      else if (ndim == 3)
      {
        if (!optvar.getLockIso2d()) ntot++;
        if (!optvar.getLockNo3d()) ntot++;
      }
      else
      {
        for (idim = 1; idim < ndim; idim++)
          ntot++;
      }
    }

    /* Anisotropy angles */
    if (DEFINE_ANIROT)
    {
      if (ndim == 2)
      {
        if (TAKE_ROT)
        {
          first_covrot = jcov;
          ntot++;
        }
      }
      else if (ndim == 3 && optvar.getLockRot2d())
      {
        if (TAKE_ROT)
        {
          first_covrot = jcov;
          ntot++;
        }
      }
      else
      {
        if (TAKE_ROT)
        {
          first_covrot = jcov;
          for (idim = 0; idim < ndim; idim++)
            ntot++;
        }
      }
    }
  }

  /* Allocation */

  param.resize(ntot);
  lower.resize(ntot);
  upper.resize(ntot);
  for (iparam = 0; iparam < ntot; iparam++)
    param[iparam] = lower[iparam] = upper[iparam] = TEST;

  return (ntot);
}

/****************************************************************************/
/*!
 **  Fill the array of pointers on the experimental conditions
 **
 ** \param[in]  npadir  Maximum number of lags for all directions
 **
 ** \param[out] gg      Allocated array of experimental uncompressed values
 ** \param[out] wt      Allocated array of experimental uncompressed weights
 **
 *****************************************************************************/
static void st_load_vmap(int npadir, VectorDouble &gg, VectorDouble &wt)
{
  int ijvar, nvar, nech, iech, nvs2, ipadir, ntest;
  double value, dist, wgt;

  /* Initializations */

  nech = DBMAP->getSampleNumber();
  nvar = DBMAP->getVariableNumber();
  nvs2 = nvar * (nvar + 1) / 2;
  db_index_sample_to_grid(DBMAP, nech / 2, INDG1);

  /* Load the Experimental conditions structure */

  ipadir = 0;
  for (iech = 0; iech < nech; iech++)
  {
    db_index_sample_to_grid(DBMAP, iech, INDG2);
    dist = distance_intra(DBMAP, nech / 2, iech, NULL);
    wgt = (dist > 0) ? 1. / dist :
                       0.;

    /* Check samples containing only undefined values */

    ntest = 0;
    for (ijvar = 0; ijvar < nvs2; ijvar++)
      if (!FFFF(DBMAP->getVariable(iech, ijvar))) ntest++;
    if (ntest <= 0) continue;

    for (ijvar = 0; ijvar < nvs2; ijvar++)
    {
      WT(ijvar,ipadir)= 0.;
      GG(ijvar,ipadir) = 0.;

      value = DBMAP->getVariable(iech,ijvar);
      if (FFFF(value)) continue;

      WT(ijvar,ipadir) = wgt;
      GG(ijvar,ipadir) = value;
    }
    ipadir++;
  }

  return;
}

/****************************************************************************/
/*!
 **  Automatic model fitting
 **
 ** \return  Error returned code
 **
 ** \param[in]  dbmap       Db Grid structure containing the Vmap
 ** \param[in]  model       Model structure containing the basic structures
 ** \param[in]  verbose     Verbose flag
 ** \param[in]  mauto_arg   Option_AutoFit structure
 ** \param[in]  cons_arg    Constraints structure
 ** \param[in]  optvar_arg  Opt_Vario structure
 **
 *****************************************************************************/
int vmap_auto_fit(const Db *dbmap,
                  Model *model,
                  bool verbose,
                  const Option_AutoFit &mauto_arg,
                  const Constraints &cons_arg,
                  const Option_VarioFit &optvar_arg)
{
  int i, error, status, nbexp, norder, npar0, npar, npadir, ndim;
  int flag_reduce, ncova, nvar, idim;
  double hmax, gmax;
  StrMod *strmod;
  VectorDouble varchol, scale, param, lower, upper;

  // Copy of const reference classes

  Option_AutoFit mauto = mauto_arg;
  Constraints constraints = cons_arg;
  Option_VarioFit optvar = optvar_arg;

  /* Initializations */

  nbexp = status = norder = npadir = 0;
  hmax = gmax = 0.;
  strmod = nullptr;
  ncova = model->getCovaNumber();
  nvar = model->getVariableNumber();
  ndim = model->getDimensionNumber();
  VectorDouble angles;
  mauto.setVerbose(verbose);

  /* Preliminary checks */

  error = 1;
  if (ndim > 3)
  {
    messerr("Procedure cannot be used for space dimension (%d) larger than 3",
            ndim);
    goto label_end;
  }
  nvar = model->getVariableNumber();
  if (nvar != dbmap->getVariableNumber())
  {
    messerr("Number of variables in Db (%d) must match the one in Model (%d)",
            model->getVariableNumber(), dbmap->getVariableNumber());
    goto label_end;
  }
  if (!FFFF(mauto.getConstantSillValue()))
  {
    if (!optvar.getFlagGoulardUsed())
    {
      messerr("When Constraints on the sum of Sills are defined");
      messerr("The Goulard option must be switched ON");
      goto label_end;
    }
    mauto.setConstantSills(nvar);
  }
  if (st_get_vmap_dimension(dbmap, nvar, &npadir, &nbexp)) goto label_end;
  angles.resize(ndim);
  for (idim = 0; idim < ndim; idim++)
    angles[idim] = 0.;
  if (db_extension_diag(dbmap, &hmax)) goto label_end;
  st_vmap_varchol_manage(dbmap, varchol);

  /* Scale the parameters in the Option_AutoFit structure */

  st_mauto_rescale(nvar, varchol, mauto);

  /* Free the keypair mechanism strings */

  st_keypair_sill(-1, model);
  st_keypair_results(-1, 0, 0, NULL, NULL);

  /* Create the experimental structures */

  hmax /= 2.;
  DBMAP = dbmap;
  INDG1 = db_indg_alloc(dbmap);
  if (INDG1 == nullptr) goto label_end;
  INDG2 = db_indg_alloc(dbmap);
  if (INDG2 == nullptr) goto label_end;

  /* Core allocation */

  if (st_manage_recint(mauto, 0, ndim, nvar, 0, ncova, npadir)) goto label_end;
  st_load_vmap(npadir, RECINT.gg, RECINT.wt);

  /* Generate the default values */

  npar0 = st_vmap_auto_count(dbmap, model, constraints, optvar, param, lower,
                             upper);

  /* Create the Model structures */

  strmod = st_model_auto_strmod_alloc(model, NULL, npar0, norder, hmax, angles,
                                      optvar, &npar);
  if (strmod == nullptr) goto label_end;
  if (npar == 0)
  {
    messerr("The VMAP Automatic Fitting procedure");
    messerr("does not allow using Goulard algorithm only");
    goto label_end;
  }
  scale.resize(npar);

  /* Set the default values and bounds */

  st_model_auto_pardef(strmod, npar, hmax, varchol, angles, param, lower,
                       upper);
  st_model_auto_scldef(strmod, npar, hmax, varchol, scale);
  st_model_auto_constraints_apply(strmod, npar, constraints, param, lower,
                                  upper);

  /* Minimization algorithm */

  STRMOD = strmod;
  MAUTO = mauto;
  ST_PREPAR_GOULARD = st_prepar_goulard_vmap;
  do
  {
    st_model_auto_strmod_print(1, strmod, mauto, param, lower, upper, npar,
                               nbexp);

    status = foxleg_f(nbexp, npar, 0, VectorDouble(), param, lower, upper,
                      scale, mauto, 0, st_strmod_vmap_evaluate, RECINT.gg,
                      RECINT.wt);
    if (status > 0) goto label_end;

    st_model_auto_strmod_print(0, strmod, mauto, param, lower, upper, npar,
                               nbexp);
    flag_reduce = st_model_auto_strmod_reduce(strmod, &npar, hmax, gmax, param,
                                              lower, upper, mauto);

    /* Reset the parameters to their default values */

    if (status < 0)
    {
      for (i = 0; i < npar; i++)
        param[i] = lower[i] = upper[i] = TEST;
      st_model_auto_pardef(strmod, npar, hmax, varchol, angles, param, lower,
                           upper);
    }
    st_model_auto_scldef(strmod, npar, hmax, varchol, scale);
  }
  while (flag_reduce && npar > 0);

  /* Perform the last cosmetic updates */

  st_model_post_update(strmod, optvar);

  /* Set the returned error code */

  if (mauto.getVerbose() >= 0 && status < 0)
    messerr("Convergence not reached after %d iterations (%d parameters)",
            mauto.getMaxiter(), npar);

  /* Store the sills in the keypair mechanism */

  st_keypair_sill(1, model);

  error = 0;

  label_end: INDG1 = db_indg_free(INDG1);
  INDG2 = db_indg_free(INDG2);
  strmod = st_model_auto_strmod_free(strmod);
  return (error);
}

/****************************************************************************/
/*!
 **  Print the Auto Fitting Constraints Structure
 **
 ** \param[in]  constraints  Constraints structure
 **
 *****************************************************************************/
void constraints_print(const Constraints &constraints)
{
  constraints.display();
}

/****************************************************************************/
/*!
 **  If a constraint concerns a sill, take its square root
 **  as it corresponds to a constraints on AIC (not on a sill directly)
 **  due to the fact that it will be processed in FOXLEG (not in GOULARD)
 **  This transform only makes sense for MONOVARIATE case (the test should
 **  have been performed beforehand)
 **
 ** \return Error code (if the sill constraint is negative)
 **
 ** \param[in]  constraints  Constraints structure
 **
 *****************************************************************************/
int modify_constraints_on_sill(Constraints &constraints)

{
  for (int i = 0; i < (int) constraints.getConsItemNumber(); i++)
  {
    const ConsItem *consitem = constraints.getConsItems(i);
    if (consitem->getType() != EConsElem::SILL) continue;
    if (consitem->getValue() < 0) return (1);
    constraints.setValue(i, sqrt(consitem->getValue()));
  }
  return (0);
}

/****************************************************************************/
/*!
 **  Return the constraint value (if defined) or TEST
 **
 ** \return Returned value or TEST
 **
 ** \param[in,out]  constraints  Constraints structure
 ** \param[in]      icase        Parameter type (EConsType)
 ** \param[in]      igrf         Rank of the Gaussian Random Function
 ** \param[in]      icov         Rank of the structure (starting from 0)
 ** \param[in]      icons        Type of the constraint (EConsElem)
 ** \param[in]      iv1          Rank of the first variable
 ** \param[in]      iv2          Rank of the second variable
 **
 *****************************************************************************/
double constraints_get(const Constraints &constraints,
                       const EConsType &icase,
                       int igrf,
                       int icov,
                       const EConsElem &icons,
                       int iv1,
                       int iv2)
{
  if (!constraints.isDefined()) return (TEST);

  for (int i = 0; i < (int) constraints.getConsItemNumber(); i++)
  {
    const ConsItem *item = constraints.getConsItems(i);
    if (item->getIGrf() != igrf || item->getICov() != icov
        || item->getType() != icons || item->getIV1() != iv1) continue;
    if (icons == EConsElem::SILL && item->getIV2() != iv2) continue;

    if (item->getIcase() == EConsType::EQUAL)
    {
      if (icase == EConsType::LOWER || icase == EConsType::UPPER)
        return (item->getValue());
    }
    else
    {
      if (icase == item->getIcase()) return (item->getValue());
    }
  }
  return (TEST);
}
