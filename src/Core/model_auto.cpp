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

#include "Enum/EAnam.hpp"
#include "Enum/EJustify.hpp"

#include "Basic/AException.hpp"
#include "Basic/File.hpp"
#include "Basic/Utilities.hpp"
#include "Basic/String.hpp"
#include "Basic/OptDbg.hpp"
#include "Basic/VectorHelper.hpp"
#include "Covariances/CovAniso.hpp"
#include "Covariances/CovLMCTapering.hpp"
#include "Covariances/CovLMCConvolution.hpp"
#include "Covariances/CovLMCAnamorphosis.hpp"
#include "Model/Option_AutoFit.hpp"
#include "Model/Model.hpp"
#include "Model/Constraints.hpp"
#include "Model/ConsItem.hpp"
#include "Db/Db.hpp"
#include "Db/DbGrid.hpp"
#include "Variogram/Vario.hpp"
#include "Geometry/GeometryHelper.hpp"
#include "Matrix/MatrixSquareGeneral.hpp"

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
#define DEFINE_T_RANGE (model->getCovMode() == EModelProperty::TAPE)
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
#define VARCHOL(ivar,jvar)       varchol[COMP_INDEX(ivar,jvar)]

#define AA(icov,jcov)            aa[(icov) * ncova + (jcov)]
#define AUX(ivar,jvar)           aux[AD(ivar,jvar)]

#define AIC(icov,ijvar)          aic[(icov)*nvs2 + (ijvar)]
#define ALPHAK(icov,ijvar)       alphak[(icov)*nvs2 + (ijvar)]

#define CORRECT(idir,k)         (! isZero(vario->getHhByIndex(idir,k)) && ! FFFF(vario->getHhByIndex(idir,k)) && \
                                 ! isZero(vario->getSwByIndex(idir,k)) && ! FFFF(vario->getSwByIndex(idir,k)) && \
                                 ! FFFF(vario->getGgByIndex(idir,k)))
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
static Constraints CONSTRAINTS;
static int *INDG1;
static int *INDG2;
static const DbGrid *DBMAP;
static void (*ST_PREPAR_GOULARD)(int imod);
static Recint RECINT;

static void st_modify_optvar_for_anam(Model* model, Option_VarioFit &optvar)
{
  const CovLMCAnamorphosis* covanam = dynamic_cast<const CovLMCAnamorphosis*>(model->getCovAnisoList());
  if (covanam != nullptr)
  {
    const EAnam anamtype = covanam->getAnamType();
    if (anamtype != EAnam::HERMITIAN && optvar.getFlagGoulardUsed())
      optvar.setFlagGoulardUsed(0);
  }
}

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
  int divide, iic;

  int value = parid;
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
  int value = imod;
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
  int flag_range, flag_param, flag_aniso, flag_rotation;
  int min_order, max_ndim, flag_int_1d, flag_int_2d;
  double scalfac, parmax;

  /* Initializations */

  Option_VarioFit optvar = strmod->optvar;
  int ndim = strmod->models[0]->getDimensionNumber();
  int nvar = strmod->models[0]->getVariableNumber();
  int first_covrot = -1;

  /* Core allocation */

  strmod->parid.resize(npar0);

  /* Loop on the models */

  int ntot = 0;
  for (int imod = 0; imod < strmod->nmodel; imod++)
  {
    Model* model = strmod->models[imod];

    /* Loop on the basic structures */

    for (int jcov = 0; jcov < model->getCovaNumber(); jcov++)
    {
      model_cova_characteristics(model->getCovaType(jcov), cov_name,
                                 &flag_range, &flag_param, &min_order,
                                 &max_ndim, &flag_int_1d, &flag_int_2d,
                                 &flag_aniso, &flag_rotation, &scalfac,
                                 &parmax);

      /* AIC coefficients -> Sill */
      if (DEFINE_AIC)
        for (int ivar = 0; ivar < nvar; ivar++)
          for (int jvar = 0; jvar <= ivar; jvar++)
            strmod->parid[ntot++] = st_parid_encode(imod, jcov, EConsElem::SILL, ivar, jvar);

      /* Third parameter */
      if (DEFINE_THIRD)
      {
        strmod->parid[ntot++] = st_parid_encode(imod, jcov, EConsElem::PARAM, 0, 0);
      }

      /* Range */
      if (DEFINE_RANGE)
      {
        strmod->parid[ntot++] = st_parid_encode(imod, jcov, EConsElem::RANGE, 0, 0);
      }

      /* Anisotropy coefficients */
      if (DEFINE_ANICOEF)
      {
        if (ndim == 2)
        {
          strmod->parid[ntot++] = st_parid_encode(imod, jcov, EConsElem::RANGE, 1, 0);
        }
        else if (ndim == 3)
        {
          if (!optvar.getLockIso2d())
          {
            strmod->parid[ntot++] = st_parid_encode(imod, jcov, EConsElem::RANGE, 1, 0);
          }
          if (!optvar.getLockNo3d())
          {
            strmod->parid[ntot++] = st_parid_encode(imod, jcov, EConsElem::RANGE, 2, 0);
          }
        }
        else
        {
          for (int idim = 1; idim < ndim; idim++)
          {
            strmod->parid[ntot++] = st_parid_encode(imod, jcov, EConsElem::RANGE, idim, 0);
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
            strmod->parid[ntot++] = st_parid_encode(imod, jcov, EConsElem::ANGLE, 0, 0);
          }
        }
        else if (ndim == 3 && optvar.getLockRot2d())
        {
          if (TAKE_ROT)
          {
            first_covrot = jcov;
            strmod->parid[ntot++] = st_parid_encode(imod, jcov, EConsElem::ANGLE, 0, 0);
          }
        }
        else
        {
          if (TAKE_ROT)
          {
            first_covrot = jcov;
            for (int idim = 0; idim < ndim; idim++)
            {
              strmod->parid[ntot++] = st_parid_encode(imod, jcov, EConsElem::ANGLE, idim, 0);
            }
          }
        }
      }

      /* Tapering Range */
      if (DEFINE_T_RANGE)
      {
        strmod->parid[ntot++] = st_parid_encode(imod, 0, EConsElem::T_RANGE, 0, 0);
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
  return nullptr;
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
  int error = 1;

  /* Core allocation */

  StrMod* strmod = new StrMod;
  strmod->npar_init = npar0;
  strmod->norder = norder;
  strmod->models[0] = model1;
  strmod->models[1] = model2;
  strmod->optvar = optvar;
  strmod->user_data = NULL;
  strmod->parid = VectorInt();
  strmod->covtab = VectorDouble();

  /* Count the number of models defined */

  int ncovmax = 0;
  int nvar = 0;
  int nmodel = 0;
  for (int i = 0; i < 2; i++)
  {
    Model* model = strmod->models[i];
    if (model == nullptr) break;
    nmodel++;
    nvar = model->getVariableNumber();
    if (ncovmax < model->getCovaNumber()) ncovmax = model->getCovaNumber();

    /* Set the default value for the range */
    /* For models where range is not asked, as it is redundant with sill */

    model->setField(hmax);
    for (int icov = 0; icov < model->getCovaNumber(); icov++)
    {
      // Set the default range

      CovAniso *cova = model->getCova(icov);
      cova->setRangeIsotropic(hmax);

      // Set the default values for the sill matrix

      for (int ivar = 0; ivar < nvar; ivar++)
        for (int jvar = 0; jvar < nvar; jvar++)
        {
          double sill = (ivar == jvar) ? 1. : 0.;
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

  label_end:
  if (error) strmod = st_model_auto_strmod_free(strmod);
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
static int st_get_vario_dimension(Vario *vario,
                                  int *nbexp_ret,
                                  int *npadir_ret)

{
  int nbexp = 0;
  int npadir = 0;
  int nvar = vario->getVariableNumber();

  // Possibly update the distance for first lag
  // if equal to 0 but corresponds to lots of pairs attached
  // This patch is not performed for asymetrical case as the h=0 is only conventional.
  for (int idir = 0; idir < vario->getDirectionNumber(); idir++)
  {
    for (int ivar = 0; ivar < nvar; ivar++)
      for (int jvar = 0; jvar <= ivar; jvar++)
      {
        int iad0 = vario->getCenter(ivar, jvar, idir);
        double sw0 = vario->getSwByIndex(idir, iad0);
        double hh0 = vario->getHhByIndex(idir, iad0);
        // The test on the number of pairs avoids hacking in the case
        // of a conventional construction where the number of pairs
        // for the first lag is arbitrarily set to 1.
        if (isZero(hh0) && sw0 > 1.)
        {
          int iad = vario->getNext(idir, ivar, jvar);
          double sw1 = vario->getSwByIndex(idir, iad);
          double hh1 = vario->getHhByIndex(idir, iad);

          if (! vario->getFlagAsym())
          {
            hh0 = hh1 * sw0 / sw1;
            vario->setHhByIndex(idir, iad0, hh0);
          }
        }
      }
  }

  /* Calculate the total number of lags */

  for (int idir = 0; idir < vario->getDirectionNumber(); idir++)
  {
    npadir += vario->getLagTotalNumber(idir);
    for (int ipas = 0; ipas < vario->getLagNumber(idir); ipas++)
      for (int ivar = 0; ivar < nvar; ivar++)
        for (int jvar = 0; jvar <= ivar; jvar++)
        {
          int i = vario->getDirAddress(idir, ivar, jvar, ipas, false, 1);
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
  int nbexp = 0;
  int npadir = 0;
  int nvs2 = nvar * (nvar + 1) / 2;
  int nech = dbmap->getSampleNumber();

  /* Calculate the total number of lags */

  for (int iech = 0; iech < nech; iech++)
  {
    int ndef = 0;
    for (int ijvar = 0; ijvar < nvs2; ijvar++)
      if (!FFFF(dbmap->getLocVariable(ELoc::Z,iech, ijvar))) ndef++;
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
  // Compute the mean variance over the components
  double total = 0.;
  for (int ivar = 0; ivar < nvar; ivar++)
    total += VARCHOL(ivar,ivar) * VARCHOL(ivar,ivar);
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
  static bool local_verbose;
  static bool local_converge;

  /* Dispatch */

  if (mode == 0)
  {
    local_verbose = mauto.getVerbose();
    local_converge = OptDbg::query(EDbg::CONVERGE);
    mauto.setVerbose(false);
    OptDbg::undefine(EDbg::CONVERGE);
  }
  else
  {
    mauto.setVerbose(local_verbose);
    if (local_converge)
      OptDbg::define(EDbg::CONVERGE);
    else
      OptDbg::undefine(EDbg::CONVERGE);
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
  for (int idir = 0, ndir = vario->getDirectionNumber(); idir < ndir; idir++)
    for (int ipas = 0, npas = vario->getLagNumber(idir); ipas < npas; ipas++, ipadir++)
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
  int iad0 = vario->getDirAddress(idir, ivar, jvar, 0, false, 0);
  int iad = iad0;
  if (! isZero(vario->getGgByIndex(idir, iad)) || vario->getSwByIndex(idir, iad) > 0)
    goto label_end;

  for (int ipas = 0, npas = vario->getLagNumber(idir); ipas < npas; ipas++)
  {
    iad = vario->getDirAddress(idir, ivar, jvar, ipas, false, 1);
    if (! isZero(vario->getGgByIndex(idir, iad))) goto label_end;
    iad = vario->getDirAddress(idir, ivar, jvar, ipas, false, -1);
    if (! isZero(vario->getGgByIndex(idir, iad))) goto label_end;
  }
  iad = iad0;

  label_end:
  return (vario->getGgByIndex(idir, iad));
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
  int nvar = vario->getVariableNumber();
  int ndim = vario->getDimensionNumber();

  /* Load the Experimental conditions structure */

  int ecr = 0;
  int ipadir = 0;
  for (int idir = 0, ndir = vario->getDirectionNumber(); idir < ndir; idir++)
  {
    for (int ipas = 0, npas = vario->getLagNumber(idir); ipas < npas; ipas++, ipadir++)
    {
      int ijvar = 0;
      for (int ivar = ijvar = 0; ivar < nvar; ivar++)
        for (int jvar = 0; jvar <= ivar; jvar++, ijvar++)
        {

          /* Calculate the variogram value */

          double dist = 0.;
          GG(ijvar,ipadir)= TEST;
          if (vario->getFlagAsym())
          {
            int iad = vario->getDirAddress(idir,ivar,jvar,ipas,false,1);
            int jad = vario->getDirAddress(idir,ivar,jvar,ipas,false,-1);
            double c00 = st_get_c00(vario,idir,ivar,jvar);
            double n1 = vario->getSwByIndex(idir,iad);
            double n2 = vario->getSwByIndex(idir,jad);
            if (n1 + n2 > 0)
            {
              double g1 = vario->getGgByIndex(idir,iad);
              double g2 = vario->getGgByIndex(idir,jad);
              if (CORRECT(idir,iad) && CORRECT(idir,jad))
              {
                GG(ijvar,ipadir) = c00 - (n1 * g1 + n2 * g2) / (n1 + n2);
                dist = (ABS(vario->getHhByIndex(idir,iad)) +
                        ABS(vario->getHhByIndex(idir,jad))) / 2.;
              }
            }
          }
          else
          {
            int iad = vario->getDirAddress(idir,ivar,jvar,ipas,false,1);
            if (CORRECT(idir,iad))
            {
              GG(ijvar,ipadir) = vario->getGgByIndex(idir,iad);
              dist = ABS(vario->getHhByIndex(idir,iad));
            }
          }

          /* Define the item of the StrExp array (if defined) */

          if (! strexps.empty())
          {
            int i = vario->getDirAddress(idir, ivar, jvar, ipas, false, 1);
            if (! CORRECT(idir, i)) continue;

            strexps[ecr].ivar = ivar;
            strexps[ecr].jvar = jvar;

            for (int idim=0; idim<ndim; idim++)
              strexps[ecr].dd[idim] = dist * vario->getCodir(idir,idim);
            ecr++;
          }
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
  CovCalcMode mode(ECalcMember::LHS);
  mode.setAsVario(true);
  mode.setUnitary(true);
  mode.setOrderVario(STRMOD->norder);

  /* Loop on the basic structures */

  for (int icov = 0, ncov = model->getCovaNumber(); icov < ncov; icov++)
  {
    mode.setActiveCovListFromOne(icov);

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
            GE(icov,ijvar,ipadir) = model->evalIvarIpas(1., d0, ivar, jvar, &mode);
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
  if (vario->getCalcul() == ECalcVario::GENERAL1) norder = 1;
  if (vario->getCalcul() == ECalcVario::GENERAL2) norder = 2;
  if (vario->getCalcul() == ECalcVario::GENERAL3) norder = 3;
  VectorDouble d1(ndim);
  CovCalcMode mode = CovCalcMode(ECalcMember::LHS);
  mode.setAsVario(true);
  mode.setUnitary(true);
  mode.setOrderVario(norder);

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
      for (int ipas = 0, npas = vario->getLagNumber(idir); ipas < npas; ipas++, ipadir++)
      {
        int ijvar = 0;
        for (int ivar = 0; ivar < nvar; ivar++)
          for (int jvar = 0; jvar <= ivar; jvar++, ijvar++)
          {
            int shift = ijvar * vario->getLagTotalNumber(idir);
            if (!ge.empty()) GE(icov,ijvar,ipadir) = 0.;

            double dist = 0.;
            if (vario->getFlagAsym())
            {
              int iad = shift + vario->getLagNumber(idir) + ipas + 1;
              int jad = shift + vario->getLagNumber(idir) - ipas - 1;
              if (!CORRECT(idir, iad) || !CORRECT(idir, jad)) continue;
              dist = (ABS(vario->getHhByIndex(idir,iad)) + ABS(vario->getHhByIndex(idir,jad))) / 2.;
            }
            else
            {
              int iad = shift + ipas;
              if (!CORRECT(idir, iad)) continue;
              dist = ABS(vario->getHhByIndex(idir, iad));
            }
            for (int idim = 0; idim < ndim; idim++)
              d1[idim] = dist * vario->getCodir(idir, idim);

            if (!ge.empty())
              GE(icov,ijvar,ipadir) = cova->evalIvarIpas(1.,d1,ivar,jvar,&mode);

            if (!dd.empty())
              for (int idim = 0; idim < ndim; idim++)
                DD(idim,ijvar,ipadir) = d1[idim];
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
  int ipadir;

  /* Initializations */

  int ndir = vario->getDirectionNumber();
  int nvar = vario->getVariableNumber();
  int nvs2 = nvar * (nvar + 1) / 2;
  VectorDouble flag(ndir);

  /* Determine the count of significant directions */

  for (int idir = 0; idir < ndir; idir++)
  {
    flag[idir] = 0.;
    for (int ipas = 0, npas = vario->getLagNumber(idir); ipas < npas; ipas++)
      for (int ijvar = 0; ijvar < nvs2; ijvar++)
      {
        int shift = ijvar * vario->getLagTotalNumber(idir);
        if (vario->getFlagAsym())
        {
          int iad = shift + vario->getLagNumber(idir) + ipas + 1;
          int jad = shift + vario->getLagNumber(idir) - ipas - 1;
          double n1 = vario->getSwByIndex(idir, iad);
          double n2 = vario->getSwByIndex(idir, jad);
          if (CORRECT(idir, iad)) flag[idir] += n1;
          if (CORRECT(idir, jad)) flag[idir] += n2;
        }
        else
        {
          int iad = shift + ipas;
          double nn = vario->getSwByIndex(idir, iad);
          if (CORRECT(idir, iad)) flag[idir] += nn;
        }
      }
  }

  switch (wmode)
  {
    case 1:
      ipadir = 0;
      for (int idir = 0; idir < ndir; idir++)
      {
        for (int ipas = 0, npas = vario->getLagNumber(idir); ipas < npas; ipas++, ipadir++)
        {
          if (isZero(flag[idir])) continue;
          for (int ijvar = 0; ijvar < nvs2; ijvar++)
          {
            int shift = ijvar * vario->getLagTotalNumber(idir);
            if (vario->getFlagAsym())
            {
              int iad = shift + vario->getLagNumber(idir) + ipas + 1;
              int jad = shift + vario->getLagNumber(idir) - ipas - 1;
              if (CORRECT(idir,iad) && CORRECT(idir, jad))
              WT(ijvar,ipadir)= flag[idir];
            }
            else
            {
              int iad = shift + ipas;
              if (CORRECT(idir,iad))
              WT(ijvar,ipadir) = flag[idir];
            }
          }
        }
      }
      break;

      case 2:
      ipadir = 0;
      for (int idir=0; idir<ndir; idir++)
      {
        for (int ipas=0, npas = vario->getLagNumber(idir); ipas < npas; ipas++,ipadir++)
        {
          if (isZero(flag[idir])) continue;
          for (int ijvar=0; ijvar<nvs2; ijvar++)
          {
            int shift = ijvar * vario->getLagTotalNumber(idir);
            if (vario->getFlagAsym())
            {
              int iad = shift + vario->getLagNumber(idir) + ipas + 1;
              int jad = shift + vario->getLagNumber(idir) - ipas - 1;
              if (! (CORRECT(idir,iad) && CORRECT(idir,jad))) continue;
              double n1 = vario->getSwByIndex(idir,iad);
              double n2 = vario->getSwByIndex(idir,jad);
              double d1 = ABS(vario->getHhByIndex(idir,iad));
              double d2 = ABS(vario->getHhByIndex(idir,jad));
              if (d1 > 0 && d2 > 0)
              WT(ijvar,ipadir) = sqrt((n1+n2) * (n1+n2) / (n1*d1+n2*d2) / 2.);
            }
            else
            {
              int iad = shift + ipas;
              if (! CORRECT(idir,iad)) continue;
              double nn = vario->getSwByIndex(idir,iad);
              double dd = ABS(vario->getHhByIndex(idir,iad));
              if (dd > 0)
              WT(ijvar,ipadir) = nn / dd;
            }
          }
        }
      }
      break;

      case 3:
      ipadir = 0;
      for (int idir=0; idir<ndir; idir++)
      {
        for (int ipas=0, npas=vario->getLagNumber(idir); ipas < npas; ipas++,ipadir++)
        {
          if (isZero(flag[idir])) continue;
          for (int ijvar=0; ijvar<nvs2; ijvar++)
          {
            int shift = ijvar * vario->getLagTotalNumber(idir);
            if (vario->getFlagAsym())
            {
              int iad = shift + vario->getLagNumber(idir) + ipas + 1;
              int jad = shift + vario->getLagNumber(idir) - ipas - 1;
              if (CORRECT(idir,iad) && CORRECT(idir,jad))
              WT(ijvar,ipadir) = 1. / vario->getLagNumber(idir);
            }
            else
            {
              int iad = shift + ipas;
              if (CORRECT(idir,iad))
              WT(ijvar,ipadir) = 1. / vario->getLagNumber(idir);
            }
          }
        }
      }
      break;

      default:
      ipadir = 0;
      for (int idir=0; idir<ndir; idir++)
      {
        for (int ipas=0, npas=vario->getLagNumber(idir); ipas < npas; ipas++,ipadir++)
        {
          if (isZero(flag[idir])) continue;
          for (int ijvar=0; ijvar<nvs2; ijvar++)
          {
            int shift = ijvar * vario->getLagTotalNumber(idir);
            if (vario->getFlagAsym())
            {
              int iad = shift + vario->getLagNumber(idir) + ipas + 1;
              int jad = shift + vario->getLagNumber(idir) - ipas - 1;
              if (CORRECT(idir,iad) && CORRECT(idir,jad))
              WT(ijvar,ipadir) = 1.;
            }
            else
            {
              int iad = shift + ipas;
              if (CORRECT(idir,iad))
              WT(ijvar,ipadir) = 1.;
            }
          }
        }
      }
      break;
    }

    /* Scaling by direction and by variable */

  for (int ijvar = 0; ijvar < nvs2; ijvar++)
  {
    ipadir = 0;
    for (int idir = 0; idir < ndir; idir++)
    {
      double total = 0.;
      for (int ipas = 0, npas = vario->getLagNumber(idir); ipas < npas; ipas++, ipadir++)
      {
        if (isZero(flag[idir])) continue;
        if (WT(ijvar,ipadir)> 0 && ! FFFF(WT(ijvar,ipadir)))
        total += WT(ijvar,ipadir);
      }
      if (isZero(total)) continue;
      ipadir -= vario->getLagNumber(idir);
      for (int ipas = 0, npas = vario->getLagNumber(idir); ipas < npas; ipas++, ipadir++)
      {
        if (isZero(flag[idir])) continue;
        if (WT(ijvar,ipadir)> 0 && ! FFFF(WT(ijvar,ipadir)))
        WT(ijvar,ipadir) /= total;
      }
    }
  }

  /* Scaling by variable variances */

  int ijvar0 = 0;
  for (int ivar = 0; ivar < nvar; ivar++)
    for (int jvar = 0; jvar <= ivar; jvar++, ijvar0++)
    {
      double ratio = (vario->getVar(ivar, jvar) > 0 && vario->getVar(jvar, ivar) > 0) ?
          sqrt(vario->getVar(ivar,jvar) * vario->getVar(jvar,ivar)) : 1.;
      ipadir = 0;
      for (int idir = 0; idir < ndir; idir++)
      {
        for (int ipas = 0, npas = vario->getLagNumber(idir); ipas < npas; ipas++, ipadir++)
          if (!FFFF(WT(ijvar0, ipadir))) WT(ijvar0,ipadir)/= ratio;
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
 ** \param[in]  iter        Number of iterations
 ** \param[in]  crit        Convergence criterion
 **
 *****************************************************************************/
static void st_goulard_score(const Option_AutoFit &mauto,
                             int mode,
                             int ncova,
                             int iter,
                             double crit)
{
  if (mauto.getVerbose())
  {
    if (mode == 0)
      mestitle(1, "Statistics for Goulard algorithm");
    else
      mestitle(1, "Statistics for Goulard algorithm (with sill constraints)");
    message("Number of sills fitted   = %d\n", ncova);
    message("Number of iterations     = %d/%d\n", iter, mauto.getMaxiter());
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
  static char loc_string[20];

  if (!OptDbg::query(EDbg::CONVERGE)) return;
  mestitle(1, "Trajectory of parameters in Goulard Algorithm");
  message("(Sti(V1-V2) : Sill for structure 'i' for variables 'V1' and 'V2'\n");
  tab_prints(NULL, "Iteration");
  tab_prints(NULL, "Score");
  for (int icov = 0; icov < ncova; icov++)
    for (int ivar = 0; ivar < nvar; ivar++)
      for (int jvar = 0; jvar <= ivar; jvar++)
      {
        (void) gslSPrintf(loc_string, "St%d(%d-%d)", icov + 1, ivar + 1, jvar + 1);
        tab_prints(NULL, loc_string);
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
  if (!OptDbg::query(EDbg::CONVERGE)) return;
  int nvs2 = nvar * (nvar + 1) / 2;
  tab_printi(NULL, iter + 1);
  if (FFFF(crit))
    tab_prints(NULL, "     ");
  else
    tab_printd(NULL, crit);

  for (int icov = 0; icov < ncova; icov++)
  {
    int ijvar = 0;
    for (int ivar = 0; ivar < nvar; ivar++)
      for (int jvar = 0; jvar <= ivar; jvar++, ijvar++)
        tab_printg(NULL, SILL(icov, ijvar));
  }
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
  char loc_string[100];

  if (model == nullptr) return;
  int ncova = model->getCovaNumber();
  int nvar  = model->getVariableNumber();

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
                  model->getSillValues(icova).getValues().data());
    }
  }
}

/****************************************************************************/
/*!
 ** Store the results using the keypair mechanism
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
  int nvs2 = nvar * (nvar + 1) / 2;
  for (int icov = 0; icov < ncova; icov++)
  {
    int ijvar = 0;
    for (int ivar = 0; ivar < nvar; ivar++)
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
  int allpos;
  double temp, crit, crit_mem, value;
  VectorDouble valpro;
  const MatrixSquareGeneral* vecpro;

  /*******************/
  /* Initializations */
  /*******************/

  double sum  = 0.;
  double sum1 = 0.;
  double sum2 = 0.;
  int nvs2 = nvar * (nvar + 1) / 2;

  /* Core allocation */

  VectorDouble aic(ncova * nvs2);
  VectorDouble mp(npadir * nvs2);
  VectorDouble fk(npadir * nvs2 * ncova, 0.);
  VectorDouble alphak(ncova * nvs2);
  MatrixSquareSymmetric cc(nvar);

  /********************/
  /* Pre-calculations */
  /********************/

  int ijvar = 0;
  for (int ivar = 0; ivar < nvar; ivar++)
    for (int jvar = 0; jvar <= ivar; jvar++, ijvar++)
      for (int ipadir = 0; ipadir < npadir; ipadir++)
      {
        MP(ijvar,ipadir)= 0.;
        for (int icov=0; icov<ncova; icov++)
        MP(ijvar,ipadir) += SILL(icov,ijvar) * GE(icov,ijvar,ipadir);
      }

  for (int icov = 0; icov < ncova; icov++)
  {
    ijvar = 0;
    for (int ivar = 0; ivar < nvar; ivar++)
      for (int jvar = 0; jvar <= ivar; jvar++, ijvar++)
      {
        sum1 = sum2 = AIC(icov,ijvar)= 0;
        for (int ipadir=0; ipadir<npadir; ipadir++)
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
  }

  crit = crit_mem = 0.;
  ijvar = 0;
  for (int ivar = 0; ivar < nvar; ivar++)
    for (int jvar = 0; jvar <= ivar; jvar++, ijvar++)
      for (int ipadir = 0; ipadir < npadir; ipadir++)
      {
        if (FFFF(WT(ijvar, ipadir))) continue;
        temp = GG(ijvar,ipadir) - MP(ijvar,ipadir);
        value = (ivar != jvar) ? 2. : 1.;
        crit += value * WT(ijvar, ipadir) * temp * temp;
      }

  /***********************/
  /* Iterative procedure */
  /***********************/

  int iter;
  for (iter = 0; iter < mauto.getMaxiter(); iter++)
  {

    /* Loop on the elementary structures */

    for (int icov = 0; icov < ncova; icov++)
    {

      /* Establish the coregionalization matrix */

      ijvar = 0;
      for (int ivar = 0; ivar < nvar; ivar++)
      {
        for (int jvar = 0; jvar <= ivar; jvar++, ijvar++)
        {
          sum = 0;
          for (int ipadir = 0; ipadir < npadir; ipadir++)
          {
            MP(ijvar,ipadir)-= SILL(icov,ijvar) * GE(icov,ijvar,ipadir);
            sum += FK(icov,ijvar,ipadir) * MP(ijvar,ipadir);
          }
          value = AIC(icov,ijvar) - ALPHAK(icov,ijvar) * sum;
          cc.setValue(ivar, jvar, value);
          cc.setValue(jvar, ivar, value);
        }
      }
      /* Computing and sorting the eigen values and eigen vectors */

      if (cc.computeEigen()) return 1;
      valpro = cc.getEigenValues();
      vecpro = cc.getEigenVectors();

      int kvar = 0;
      allpos = 1;
      while ((kvar < nvar) && allpos)
      {
        if (valpro[kvar++] < 0) allpos=0;
      }

      /* Calculate the new coregionalization matrix */

      ijvar = 0;
      for (int ivar=0; ivar<nvar; ivar++)
      {
        for (int jvar=0; jvar<=ivar; jvar++, ijvar++)
        {
          if (allpos)
          {
            SILL(icov,ijvar) = cc.getValue(ivar,jvar);
          }
          else
          {
            sum = 0.;
            for (kvar=0; kvar<nvar; kvar++)
              sum += (MAX(valpro[kvar],0.) * vecpro->getValue(ivar,kvar) * vecpro->getValue(jvar,kvar));
            SILL(icov,ijvar) = sum;
          }
          for (int ipadir=0; ipadir<npadir; ipadir++)
            MP(ijvar,ipadir) += SILL(icov,ijvar) * GE(icov,ijvar,ipadir);
        }
      }
    }

    /* Update the global criterion */

    crit_mem = crit;
    crit = 0.;
    for (int ipadir=0; ipadir<npadir; ipadir++)
    {
      ijvar = 0;
      for (int ivar=0; ivar<nvar; ivar++)
        for (int jvar=0; jvar<=ivar; jvar++, ijvar++)
        {
          if (FFFF(WT(ijvar,ipadir))) continue;
          temp = GG(ijvar,ipadir) - MP(ijvar,ipadir);
          value = (ivar != jvar) ? 2. : 1.;
          crit += value * WT(ijvar,ipadir) * temp * temp;
       }
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
  else
  {
    if (!FFFF(lower_val)) lower[rank] = MAX(lower_val, lower[rank]);
  }

  /* Upper bound */

  if (FFFF(upper[rank]))
    upper[rank] = upper_val;
  else
  {
    if (!FFFF(upper_val)) upper[rank] = MIN(upper_val, upper[rank]);
  }

  /* Initial parameter */
  if (FFFF(param[rank])) param[rank] = (!FFFF(def_val)) ? def_val : 0.;

  if (!FFFF(lower[rank]) && !FFFF(upper[rank]))
  {
    if (param[rank] < lower[rank] || param[rank] > upper[rank])
    {
      if (lower[rank] > 0)
        param[rank] = (lower[rank] + upper[rank]) / 2.;
      else
        // This only happens in the case of lower bound on AIC. This avoids the mid-interval
        // which will cause a zero-gradient in model_auto.
        param[rank] = (upper[rank]) / 2.;
    }
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
  int ntot = 0;
  for (int i = 0; i < n_init; i++)
  {
    if (FFFF(param[i])) continue;
    parid[ntot] = parid[i];
    param[ntot] = param[i];
    upper[ntot] = upper[i];
    lower[ntot] = lower[i];
    ntot++;
  }
  return (ntot);
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
 ** \param[in]  strmod      StrMod structure
 ** \param[in]  constraints Constraints structure
 ** \param[in]  mauto       Option_AutoFit structure
 ** \param[in]  param       Current values for parameters
 ** \param[in]  lower       Array of lower values
 ** \param[in]  upper       Array of upper values
 ** \param[in]  npar        Number of parameters to be inferred
 ** \param[in]  nbexp      Number of experiments
 **
 *****************************************************************************/
static void st_model_auto_strmod_print(int flag_title,
                                       StrMod *strmod,
                                       const Constraints& constraints,
                                       const Option_AutoFit &mauto,
                                       VectorDouble &param,
                                       VectorDouble &lower,
                                       VectorDouble &upper,
                                       int npar,
                                       int nbexp)
{
  int icov, ivar, jvar, imod;
  EConsElem icons;
  static const char *NOK[] = { "OFF", "ON" };

  /* Initializations */

  bool skip = false;
  if (! mauto.getVerbose()) skip = true;
  if (! OptDbg::query(EDbg::CONVERGE)) skip = true;
  if (skip) return;

  Option_VarioFit optvar = strmod->optvar;
  int ndim = strmod->models[0]->getDimensionNumber();
  int nvar = strmod->models[0]->getVariableNumber();
  int imod_mem = -1;
  int icov_mem = -1;

  /* Title */

  if (flag_title)
  {
    mestitle(0, "%s", "Optimization Conditions");
    message("- Number of variables       %d  \n", nvar);
    message("- Space dimension           %d  \n", ndim);
    message("- Number of experiments     %d  \n", nbexp);
    message("- Number of parameters      %d  \n", npar);
    message("- Constrained Minimization  %s\n",
            NOK[!FFFF(constraints.getConstantSillValue())]);
    message("- Intrinsic option          %s\n",
            NOK[mauto.getFlagIntrinsic()]);
    messageFlush(optvar.toString());
  }

  /* Loop on the models */

  for (int ntot = 0; ntot < npar; ntot++)
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
  int icov, ivar, jvar, imod;
  int flag_range, flag_param, flag_aniso, flag_rotation;
  int min_order, max_ndim, flag_int_1d, flag_int_2d;
  double scalfac, parmax;
  EConsElem icons;

  /* Loop on the models */

  if (npar <= 0) return;
  for (int ntot = 0; ntot < npar; ntot++)
  {
    st_parid_decode(strmod->parid[ntot], &imod, &icov, &icons, &ivar, &jvar);
    Model* model = strmod->models[imod];
    model_cova_characteristics(model->getCovaType(icov), cov_name, &flag_range,
                               &flag_param, &min_order, &max_ndim, &flag_int_1d,
                               &flag_int_2d, &flag_aniso, &flag_rotation,
                               &scalfac, &parmax);
    switch (icons.toEnum())
    {
      case EConsElem::E_SILL:
      {
        int lec = ivar * (ivar + 1) / 2 + jvar;
        scale[ntot] = ABS(varchol[lec]) / sqrt(model->getCovaNumber());
        break;
      }

      case EConsElem::E_PARAM:
      {
        if (parmax < 0 || FFFF(parmax)) parmax = 1.;
        scale[ntot] = parmax;
        break;
      }

      case EConsElem::E_RANGE:
      {
        scale[ntot] = hmax / model->getCovaNumber(true) / 2.;
        break;
      }

      case EConsElem::E_ANGLE:
      {
        scale[ntot] = 1800.;
        break;
      }

      case EConsElem::E_T_RANGE:
      {
        scale[ntot] = hmax / 10.;
        break;
      }

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
  int icov, ivar, jvar, imod;
  double param_loc, lower_loc, upper_loc;
  EConsElem icons;

  /* Loop on the models */

  if (npar <= 0) return;
  for (int ipar = 0; ipar < npar; ipar++)
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
  int icov, ivar, jvar, imod;
  int flag_range, flag_param, flag_aniso, flag_rotation;
  int min_order, max_ndim, flag_int_1d, flag_int_2d;
  double scalfac, parmax;
  ECov type;
  EConsElem icons;

  /* Loop on the models */

  int icovm = 0;
  if (npar <= 0) return;
  for (int ntot = 0; ntot < npar; ntot++)
  {
    st_parid_decode(strmod->parid[ntot], &imod, &icov, &icons, &ivar, &jvar);
    Model* model = strmod->models[imod];
    type = model->getCovaType(icov);
    model_cova_characteristics(type, cov_name, &flag_range, &flag_param,
                               &min_order, &max_ndim, &flag_int_1d,
                               &flag_int_2d, &flag_aniso, &flag_rotation,
                               &scalfac, &parmax);
    if (type == ECov::NUGGET && ivar == 0 && jvar == 0) icovm++;
    switch (icons.toEnum())
    {
      case EConsElem::E_SILL:
      {
        int lec = ivar * (ivar + 1) / 2 + jvar;
        double dvar = varchol[lec] / sqrt(model->getCovaNumber());
        st_affect(ntot, dvar, TEST, TEST, param, lower, upper);
        break;
      }

      case EConsElem::E_PARAM:
      {
        if (parmax < 0) parmax = TEST;
        double valdef = 1.;
        if (type == ECov::COSEXP) valdef = hmax / 3.;
        st_affect(ntot, valdef, 0.001, parmax, param, lower, upper);
        break;
      }

      case EConsElem::E_RANGE:
      {
        double dunit = hmax / model->getCovaNumber(true) / 2.;
        double dmin = hmax / 1.e6;
        double dist = dunit * (icov + 1 - icovm);
        st_affect(ntot, dist, dmin, TEST, param, lower, upper);
        break;
      }

      case EConsElem::E_ANGLE:
      {
        st_affect(ntot, angles[ivar], TEST, TEST, param, lower, upper);
        break;
      }

      case EConsElem::E_T_RANGE:
      {
        double dmin = hmax / 1.e6;
        double dist = hmax / 10.;
        st_affect(ntot, dist, dmin, TEST, param, lower, upper);
        break;
      }

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
  int icov, ntot, imod, ivar, jvar;
  Model *model;
  CovAniso *cova;
  CovAniso *cova1;
  EConsElem icons;

  /* Initializations */

  int nvar = strmod->models[0]->getVariableNumber();
  int ndim = strmod->models[0]->getDimensionNumber();
  Option_VarioFit optvar = strmod->optvar;
  int size = nvar * (nvar + 1) / 2;
  VectorDouble ranges(ndim, 0.);
  VectorDouble angles(ndim, 0.);
  VectorDouble tritab(size);

  /* Loop on the parameters */

  int flag_rot = 0;
  int flag_aic = 0;
  int imod_mem = -1;
  int icov_mem = -1;
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
        cova->setRangeIsotropic(ranges[0]);
      if (flag_rot) cova->setAnisoAngles(angles);
      if (flag_aic)
      {
        MatrixSquareSymmetric* mat = MatrixSquareSymmetric::createFromTLTU(nvar, tritab);
        cova->setSill(*mat);
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
      {
        int ipos = ivar * (ivar + 1) / 2 + jvar;
        tritab[ipos] = param[ntot];
        flag_aic = 1;
        break;
      }

      case EConsElem::E_PARAM:
      {
        cova->setParam(param[ntot]);
        break;
      }

      case EConsElem::E_RANGE:
      {
        if (ivar == 0) VH::fill(ranges, param[ntot]);
        if (ivar < ndim) ranges[ivar] = param[ntot];
        break;
      }

      case EConsElem::E_ANGLE:
      {
        if (ivar < ndim) angles[ivar] = param[ntot];
        flag_rot = 1;
        break;
      }

      case EConsElem::E_T_RANGE:
      {
        model->setTapeRange(param[ntot]);
        break;
      }

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
      cova->setRangeIsotropic(ranges[0]);
    if (flag_rot) cova->setAnisoAngles(angles);
    if (flag_aic)
    {
      MatrixSquareSymmetric* mat = MatrixSquareSymmetric::createFromTLTU(nvar, tritab);
      cova->setSill(*mat);
    }
    flag_rot = flag_aic = 0;
  }

  // If 'lock_samerot' is ON, copy the rotation to all structures

  if (strmod->optvar.getLockSamerot())
  {
    for (int jmod = 0; jmod < strmod->nmodel; jmod++)
    {
      model = strmod->models[jmod];

      // Look for the first basic structure with a rotation defined

      int found = -1;
      for (int jcov = 0; jcov < model->getCovaNumber() && found < 0; jcov++)
      {
        if (model->getCova(jcov)->hasRange()) found = jcov;
      }
      if (found < 0) continue;
      cova1 = model->getCova(found);

      for (int jcov = 1; jcov < model->getCovaNumber(); jcov++)
      {
        if (jcov == found) continue;
        cova = model->getCova(jcov);
        if (!cova->getAnisoRotMat().empty())
          cova->setAnisoAngles(cova1->getAnisoAngles());
      }
    }
  }
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
  MatrixSquareGeneral tab(nvar);
  CovCalcMode mode(ECalcMember::LHS);
  mode.setAsVario(true);
  mode.setActiveCovListFromOne(icov);
  mode.setOrderVario(STRMOD->norder);
  model->evaluateMatInPlace(nullptr, d1, tab, true, 1., &mode);

  for (int ivar = 0; ivar < nvar; ivar++)
  {
    if (tab.getValue(ivar, ivar) > tolsigma * gmax / 100.) return (0);
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
  CovCalcMode mode(ECalcMember::LHS);
  mode.setAsVario(true);
  mode.setOrderVario(strmod->norder);

  /* Loop on the experimental conditions */

  for (int i = 0; i < nbexp; i++)
  {
    int ivar = strexps[i].ivar;
    int jvar = strexps[i].jvar;
    VectorDouble d0 = strexps[i].dd;
    tabge[i] = model->evalIvarIpas(1., d0, ivar, jvar, &mode);
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

  CovCalcMode mode(ECalcMember::LHS);
  mode.setAsVario(true);
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
        if (FFFF(DBMAP->getLocVariable(ELoc::Z,iech, ijvar))) continue;
        tabge[ecr++] = model->evalIvarIpas(1., d0, ivar, jvar, &mode);
      }
  }
  return;
}

/****************************************************************************/
/*!
 **  Find the rank of the 'parid' which matches several criteria
 **
 ** \return  Rank of the 'parid' or -1 if not found
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
  int imod, icov, ivar, jvar;
  EConsElem icons;

  for (int ntot = 0; ntot < npar; ntot++)
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
  MatrixSquareSymmetric cc(nvar);

  /* Loop on the basic structures */

  for (int icov = 0, ncov = model->getCovaNumber(); icov < ncov; icov++)
  {
    message("\nCheck the Sill Matrix for structure %d\n", icov + 1);

    /* Load the matrix of sills */

    for (int ivar = 0; ivar < nvar; ivar++)
      for (int jvar = 0; jvar < nvar; jvar++)
        cc.setValue(ivar,jvar,model->getSill(icov,ivar,jvar));

        /* Check definite positiveness */

    if (! cc.isDefinitePositive()) return 1;
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
  MatrixSquareSymmetric cc(nvar);
  for (int ivar = 0; ivar < nvar; ivar++)
    for (int jvar = 0; jvar < nvar; jvar++)
      cc.setValue(ivar,jvar,MATCOR(icov0,ivar,jvar));

  if (cc.computeEigen())
    messageAbort("st_truncate_negative_eigen");

  VectorDouble valpro = cc.getEigenValues();
  const MatrixSquareGeneral* vecpro = cc.getEigenVectors();

  /* Check positiveness */

  int flag_positive = 1;
  for (int ivar = 0; ivar < nvar; ivar++)
    if (valpro[ivar]<= 0) flag_positive = 0;
  if (!flag_positive)
  {

    /* Calculate the new coregionalization matrix */

    for (int ivar = 0; ivar < nvar; ivar++)
      for (int jvar = 0; jvar < nvar; jvar++)
      {
        double sum = 0.;
        for (int kvar = 0; kvar < nvar; kvar++)
          sum += MAX(valpro[kvar],0.) * vecpro->getValue(ivar, kvar) * vecpro->getValue(jvar, kvar);
        MATCORU(icov0,ivar,jvar)= sum;
      }
    }
    else
    {
      for (int ivar=0; ivar<nvar; ivar++)
        for (int jvar=0; jvar<nvar; jvar++)
          MATCORU(icov0,ivar,jvar) = MATCOR(icov0,ivar,jvar);
    }

  return flag_positive;
}

/****************************************************************************/
/*!
 **  Sum the sill of all the structures
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
  double Sr = 0;
  for (int icov = 0; icov < ncova; icov++)
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
  double score = 0.;
  int ijvar = 0;
  for (int ivar = 0; ivar < nvar; ivar++)
    for (int jvar = 0; jvar <= ivar; jvar++, ijvar++)
    {
      double coeff = (ivar == jvar) ? 1 : 2;
      for (int ipadir = 0; ipadir < npadir; ipadir++)
      {
        double dd = GG(ijvar, ipadir);
        if (FFFF(dd)) continue;
        for (int icov = 0; icov < ncova; icov++)
          dd -= MATCOR(icov,ivar,jvar) * GE(icov,ijvar,ipadir);
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
  double retval, a, c, d, s, x[3];

  /* Core allocation */

  VectorDouble Nir_v(nvar);
  VectorDouble Mrr_v(npadir);
  VectorDouble Crr_v(npadir);
  VectorDouble Airk_v(npadir * nvar);
  VectorDouble Birk_v(npadir * nvar);
  VectorDouble xx(2);
  VectorDouble xt(2);
  VectorDouble xest(2);

  int irr = st_combineVariables(ivar0, ivar0);
  int nvs2 = nvar * (nvar + 1) / 2;

  for (int ivar = 0; ivar < nvar; ivar++)
  {
    int irl = st_combineVariables(ivar0, ivar);
    Nir_v[ivar] = 0.;
    for (int k = 0; k < npadir; k++)
      Nir_v[ivar] += WT(irl,k)* GE(icov0,0,k) * GE(icov0,0,k);
    }

  for (int k = 0; k < npadir; k++)
  {
    Mrr_v[k] = 0.;
    for (int jcov = 0; jcov < ncova; jcov++)
    {
      if (jcov == icov0) continue;
      Mrr_v[k] += ALPHA(jcov,ivar0,ivar0)* (GE(jcov,0,k) - GE(icov0,0,k));
    }
  }

  for (int k = 0; k < npadir; k++)
    for (int ivar = 0; ivar < nvar; ivar++)
    {
      int irl = st_combineVariables(ivar0, ivar);
      double value = 0.;
      for (int l = 0; l < npadir; l++)
        value += WT(irl,l)* GG(irl,l) * GE(icov0,0,l);
      value *= GE(icov0,0,k)/ Nir_v[ivar];
      AIRKV(k,ivar)= GG(irl,k) - value;
    }

  for (int k = 0; k < npadir; k++)
    for (int ivar = 0; ivar < nvar; ivar++)
    {
      int irl = st_combineVariables(ivar0, ivar);
      BIRKV(k,ivar)= 0.;
      for (int icov = 0; icov < ncova; icov++)
      {
        if (icov == icov0) continue;
        double value = 0.;
        for (int l = 0; l < npadir; l++)
          value += WT(irl,l) * GE(icov,0,l) * GE(icov,0,l);
        BIRKV(k,ivar)+= ALPHA(icov,ivar0,ivar) *
        (GE(icov,0,k) - value * GE(icov0,0,k))/Nir_v[ivar];
      }
    }

  for (int k = 0; k < npadir; k++)
    Crr_v[k] = GG(irr,k)- consSill[ivar0] * GE(icov0,0,k);

  a = 0.;
  for (int k = 0; k < npadir; k++)
    a += WT(irr,k)* Mrr_v[k] * Mrr_v[k];

  c = 0.;
  for (int k = 0; k < npadir; k++)
    for (int ivar = 0; ivar < nvar; ivar++)
    {
      if (ivar != ivar0)
      {
        int irl = st_combineVariables(ivar0, ivar);
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
  for (int k = 0; k < npadir; k++)
    for (int ivar = 0; ivar < nvar; ivar++)
    {
      if (ivar == ivar0) continue;
      int irl = st_combineVariables(ivar0, ivar);
      d -= WT(irl,k)* AIRKV(k,ivar) * xr[ivar] * BIRKV(k,ivar);
    }

  d *= 4.;

  int number = solve_P3(4. * a, 0., 2. * c, d, x);

  switch (number)
  {
    case 1:
      retval = MAX(0., MIN(xrmax, x[0]));
      break;

    case 3:
    {
      int nin = 0;
      xx[0] = x[0];
      xx[1] = x[2];
      for (int k = 0; k < 2; k++)
      {
        xt[k] = MAX(0., MIN(xrmax, xx[k]));
        xest[k] = (a * xt[k] * xt[k] * xt[k] * xt[k] + c * xt[k] * xt[k]
                   + d * xt[k]) / 2.;
        if (areEqual(xt[k], xx[k])) nin++;
      }
      if (nin == 1)
      {
        retval = (areEqual(xt[0], xx[0])) ? xx[0] : xx[1];
      }
      else
      {
        retval = (xest[0] < xest[1]) ? xt[0] : xt[1];
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
  double srm = st_sum_sills(ivar0, ncova, nvar, alpha) - ALPHA(icov0, ivar0, ivar0);
  double value = consSill[ivar0] / (xr[ivar0] * xr[ivar0]) - srm;
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
  for (int jcov = 0; jcov < ncova; jcov++)
  {
    if (jcov == icov0) continue;
    for (int ivar = 0; ivar < nvar; ivar++)
    {
      MATCOR(jcov,ivar0, ivar)= MATCOR(jcov,ivar, ivar0) =
          ALPHA(jcov,ivar0, ivar) * xr[ivar0] * xr[ivar];
    }
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
  VectorDouble mv(npadir);
  int nvs2 = nvar * (nvar + 1) / 2;

  // Loop on the variables

  for (int ivar = 0; ivar < nvar; ivar++)
  {
    if (ivar == ivar0 && !FFFF(consSill[ivar0])) continue;

    for (int ilagdir = 0; ilagdir < npadir; ilagdir++)
    {
      mv[ilagdir] = 0.;
      for (int icov = 0; icov < ncova; icov++)
      {
        if (icov == icov0) continue;
        mv[ilagdir] += MATCOR(icov,ivar0,ivar) * GE(icov,0,ilagdir);
      }
    }

    double tot1 = 0.;
    double tot2 = 0.;
    int ivs2 = st_combineVariables(ivar0, ivar);

    for (int ilagdir = 0; ilagdir < npadir; ilagdir++)
    {
      double wtloc = WT(ivs2, ilagdir);
      double ggloc = GG(ivs2, ilagdir);
      double geloc = GE(icov0, 0, ilagdir);

      if (!FFFF(ggloc))
      {
        tot1 += wtloc * geloc * (ggloc - mv[ilagdir]);
        tot2 += wtloc * geloc * geloc;
      }
    }

    MATCOR(icov0,ivar0, ivar)= MATCOR(icov0,ivar, ivar0) = tot1 / tot2;
  }
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
  for (int ivar = 0; ivar < nvar; ivar++)
  {
    if (ivar == ivar0 && !FFFF(consSill[ivar0])) continue;
    double value = MATCOR(icov0,ivar0,ivar)/ (xr[ivar0] * xr[ivar]);
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
  double value = xr[ivar0] * xr[ivar0] * ALPHA(icov0, ivar0, ivar0);
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

  if (flag_positive) return flag_positive;

  for (int ivar = 0; ivar < nvar; ivar++)
  {
    if (FFFF(consSill[ivar]))
      norme1[ivar] = 1.;
    else
    {
      if (ABS(MATCOR(icov0,ivar, ivar)) > EpsFit)
        norme1[ivar] = sqrt(muold[ivar] / MATCOR(icov0, ivar, ivar));
      else
        norme1[ivar] = (ABS(muold[ivar]) < EpsFit) ? 1. : 0.;
    }
  }

  for (int ivar = 0; ivar < nvar; ivar++)
    for (int jvar = 0; jvar < nvar; jvar++)
      MATCOR(icov0,ivar, jvar)*= norme1[ivar] * norme1[jvar];

  return flag_positive;
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
 ** \param[in]  constraints Constraints structure
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
                                         const Constraints& constraints,
                                         const Option_AutoFit &mauto,
                                         VectorDouble &wt,
                                         VectorDouble &gg,
                                         VectorDouble &ge,
                                         VectorDouble &matcor,
                                         double *score)
{
  double score_old, xrmax;

  /* Core allocation */

  VectorDouble xr(nvar);
  VectorDouble alpha(nvar * nvar * ncova);
  int iter = 0;

  /* Calculate the initial score */

  double score_new = st_score(nvar, nvs2, ncova, npadir, wt, gg, ge, matcor);

  for (int icov = 0; icov < ncova; icov++)
    for (int ivar = 0; ivar < nvar; ivar++)
      for (int jvar = 0; jvar < nvar; jvar++)
        ALPHA(icov,ivar,jvar)= MATCOR(icov,ivar,jvar);

        /***********************/
        /* Iterative procedure */
        /***********************/

        /* Optional printout */

  st_goulard_debug_title(nvar, ncova);
  st_goulard_debug_current(nvar, ncova, 0, matcor, TEST);

  for (int ivar = 0; ivar < nvar; ivar++)
  {
    if (FFFF(constraints.getConstantSills(ivar)))
      xr[ivar] = 1.;
    else
      xr[ivar] = sqrt(
          constraints.getConstantSills(ivar) / st_sum_sills(ivar, ncova, nvar, alpha));
  }

  for (iter = 0; iter < mauto.getMaxiter(); iter++)
  {
    for (int icov0 = 0; icov0 < ncova; icov0++)
    {
      for (int ivar0 = 0; ivar0 < nvar; ivar0++)
      {
        double denom = st_sum_sills(ivar0, ncova, nvar, alpha)-
        ALPHA(icov0,ivar0,ivar0);
        if (!FFFF(constraints.getConstantSills(ivar0)) && denom < 1e-30) continue;
        if (!FFFF(constraints.getConstantSills(ivar0)))
        {
          xrmax = sqrt(constraints.getConstantSills(ivar0) / denom);
          xr[ivar0] = st_minimize_P4(icov0, ivar0, ncova, nvar, npadir, xrmax,
                                     xr, alpha, wt, gg, ge,
                                     constraints.getConstantSills());

          if (isZero(xr[ivar0]))
          {
            xr[ivar0] = 1.;
            for (int jcov = 0; jcov < ncova; jcov++)
            {
              if (jcov == icov0) continue;
              for (int jvar = 0; jvar < nvar; jvar++)
                ALPHA(jcov,ivar0,jvar)= ALPHA(jcov,jvar,ivar0) = 0.;
              }
            }
          }

          /* Update 'alpha' (diagonal only) */

        if (!FFFF(constraints.getConstantSills(ivar0)))
          st_updateAlphaDiag(icov0, ivar0, ncova, nvar, xr, alpha,
                             constraints.getConstantSills());

        /* Update 'sills' for the structures other than the current one */

        st_updateOtherSills(icov0, ivar0, ncova, nvar, alpha, xr, matcor);

        /* Update the sill matrix for the current structure */
        /* (except diagonal in the constrained case)        */

        st_updateCurrentSillGoulard(icov0, ivar0, npadir, nvar, ncova, wt, ge,
                                    gg, constraints.getConstantSills(), matcor);

        /* Update sill matrix for the current structure (for diagonal) */

        if (!FFFF(constraints.getConstantSills(ivar0)))
          st_updateCurrentSillDiag(icov0, ivar0, nvar, alpha, xr, matcor);

        /* Make sure the current matrix of sills if definite positive */
        /* (diagonal unchanged)                                       */

        (void) st_makeDefinitePositive(icov0, nvar, constraints.getConstantSills(),
                                       matcor);

        /* Update 'alpha' for the current structure */

        for (int ivar = 0; ivar < nvar; ivar++)
          st_updateAlphaNoDiag(icov0, ivar, nvar, xr, constraints.getConstantSills(),
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
  MatrixSquareSymmetric aa(ncova);
  VectorDouble bb(ncova);

  MatrixRectangular Ae(ncova, 1);
  VectorDouble be(1);
  MatrixRectangular Ai(ncova, ncova);
  VectorDouble bi(ncova);
  VectorDouble res(ncova);

  /* Initialize the constraints matrices */

  Ae.fill(1.);
  bi.fill(0.);
  for (int icov = 0; icov < ncova; icov++)
    for (int jcov = 0; jcov < ncova; jcov++)
      Ai.setValue(icov, jcov, (icov == jcov));

  /* Loop on the variables */

  int ijvar = 0;
  for (int ivar = 0; ivar < nvar; ivar++)
    for (int jvar = 0; jvar <= ivar; jvar++, ijvar++)
    {
      /* Reset the arrays */

      bb.fill(0.);
      aa.fill(0.);

      for (int ipadir=0; ipadir<npadir; ipadir++)
      {
        double wtloc = WT(ijvar,ipadir);
        double ggloc = GG(ijvar,ipadir);
        if (FFFF(wtloc) || FFFF(ggloc)) continue;
        for (int icov=0; icov<ncova; icov++)
        {
          double geloc1 = GE(icov,ijvar,ipadir);
          if (FFFF(geloc1)) continue;
          bb[icov] += wtloc * ggloc * geloc1;
          for (int jcov=0; jcov<=icov; jcov++)
          {
            double geloc2 = GE(jcov,ijvar,ipadir);
            if (FFFF(geloc2)) continue;
            aa.updValue(icov,jcov,EOperator::ADD,wtloc * geloc1 * geloc2);
          }
        }
      }

      int retcode = 0;
      if (ivar == jvar && ! consSill.empty() && !FFFF(consSill[ivar]))
      {
        be[0] = consSill[ivar];
        retcode = aa.minimizeWithConstraintsInPlace(bb, Ae, be, Ai, bi, res);
      }
      else
      {
        retcode = aa.minimizeWithConstraintsInPlace(bb,
                                                    MatrixRectangular(), VectorDouble(),
                                                    MatrixRectangular(), VectorDouble(),
                                                    res);
      }

      /* Update (taking into account possible constraints) */

      if (retcode)
      {
        for (int icov=0; icov<ncova; icov++)
          res[icov] = (ivar == jvar);
        if (ivar == jvar && ! consSill.empty() && !FFFF(consSill[ivar]))
        {
          double temp = consSill[ivar] / ncova;
          for(int icov=0; icov<ncova; icov++)
            res[icov] = temp;
        }
      }

      /* Store in the output matrix */

      for (int icov=0; icov<ncova; icov++)
        MATCOR(icov,ivar,jvar) = MATCOR(icov,jvar,ivar) = res[icov];
    }
  return 0;
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
  int nvs2 = nvar * (nvar + 1) / 2;
  for (int icov = 0; icov < ncova; icov++)
  {
    int ijvar = 0;
    for (int ivar = ijvar = 0; ivar < nvar; ivar++)
      for (int jvar = 0; jvar <= ivar; jvar++, ijvar++)
      {
        model->setSill(icov, ivar, jvar, SILL(icov, ijvar));
      }
  }
}

/****************************************************************************/
/*!
 **  Internal function for Goulard under constraints
 **
 ** \return  Error returned code
 **
 ** \param[in]  constraints Constraints structure
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
static int st_goulard_with_constraints(const Constraints& constraints,
                                       const Option_AutoFit &mauto,
                                       int nvar,
                                       int ncova,
                                       int npadir,
                                       VectorDouble &wt,
                                       VectorDouble &gg,
                                       VectorDouble &ge,
                                       VectorDouble &sill)
{
  double crit;
  int flag_positive, iter;

  /* Initializations */

  int nvs2 = nvar * (nvar + 1) / 2;

  /* Core allocation */

  VectorDouble matcor(nvar * nvar * ncova);

  /* Initialize the Goulard system */

  st_initialize_goulard(nvar, nvs2, ncova, npadir, wt, gg, ge,
                        constraints.getConstantSills(), matcor);

  /* Update according to the eigen values */

  flag_positive = 1;
  for (int icov = 0; icov < ncova; icov++)
    if (!st_makeDefinitePositive(icov, nvar, constraints.getConstantSills(), matcor))
      flag_positive = 0;

  if (!flag_positive)
  {
    for (int icov = 0; icov < ncova; icov++)
      for (int ivar = 0; ivar < nvar; ivar++)
      {
        if (!FFFF(constraints.getConstantSills(ivar)))
          MATCOR(icov,ivar,ivar)= constraints.getConstantSills(ivar) / ncova;
        else
          MATCOR(icov,ivar,ivar)=1.;

        for (int jvar=0;jvar<ivar;jvar++)
          MATCOR(icov,ivar,jvar) = MATCOR(icov,jvar,ivar)=0.;
      }

        /* Perform the optimization under constraints */

        iter = st_optimize_under_constraints(nvar,nvs2,ncova,npadir,constraints,mauto,
            wt,gg,ge,matcor,&crit);

        /* Optional printout */

        st_goulard_score(mauto,1,ncova,iter,crit);
      }

      /* Load the parameters in the final model */

  for (int icov = 0; icov < ncova; icov++)
  {
    int ijvar = 0;
    for (int ivar = 0; ivar < nvar; ivar++)
      for (int jvar = 0; jvar <= ivar; jvar++, ijvar++)
        SILL(icov,ijvar)= MATCOR(icov,ivar,jvar);
  }
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
  double newval, crit;

  /* Initializations */

  int error = 1;
  int nvar = model->getVariableNumber();
  int ncova = model->getCovaNumber();
  int nvs2 = nvar * (nvar + 1) / 2;
  double crit_mem = 1.e30;
  for (int icov = 0; icov < ncova; icov++)
    alphau[icov] = 1. / (double) ncova;

  /* Iterative procedure */

  Option_AutoFit mauto_new(mauto);
  mauto_new.setMaxiter(1);
  for (int iter = 0; iter < mauto.getMaxiter(); iter++)
  {

    /* Initialize the arrays for the first pass */

    for (int ipadir = 0; ipadir < npadir; ipadir++)
    {
      int ijvar = 0;
      for (int ivar = 0; ivar < nvar; ivar++)
      {
        for (int jvar = 0; jvar <= ivar; jvar++, ijvar++)
        {
          double sum = 0.;
          for (int icov = 0; icov < ncova; icov++)
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

    int ijvar = 0;
    for (int ivar=0; ivar<nvar; ivar++)
      for (int jvar=0; jvar<=ivar; jvar++, ijvar++)
      {
        double pivot = sill1[ijvar];
        for (int ipadir=0; ipadir<npadir; ipadir++)
        {
          GG2(ijvar,ipadir) = (isZero(pivot)) ? 0. : GG(ijvar,ipadir) / pivot;
          WT2(ijvar,ipadir) = WT(ijvar,ipadir) * pivot * pivot;
          for (int icov=0; icov<ncova; icov++)
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

  for (int icov = 0; icov < ncova; icov++)
  {
    int ijvar = 0;
    for (int ivar = 0; ivar < nvar; ivar++)
      for (int jvar = 0; jvar <= ivar; jvar++, ijvar++)
      {
        newval = alphau[icov] * sill1[ijvar];
        model->setSill(icov, ivar, jvar, newval);
        model->setSill(icov, jvar, ivar, newval);
      }
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
 ** \param[in]  constraints Constraints structure
 ** \param[in]  mauto       Option_AutoFit structure
 **
 *****************************************************************************/
static int st_goulard_fitting(int flag_reset,
                              int flag_title,
                              Model *model,
                              const Constraints& constraints,
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

    if (FFFF(constraints.getConstantSillValue()))
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

      status = st_goulard_with_constraints(constraints, mauto,
                                           model->getVariableNumber(),
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
  int flag_range, flag_param, flag_aniso, flag_rotation;
  int min_order, max_ndim, flag_int_1d, flag_int_2d;
  double scalfac, parmax;

  /* Loop on the basic structures */

  int n_int = 0;
  for (int icov = 0; icov < model->getCovaNumber(); icov++)
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
 ** \param[in]  constraints Constraints structure
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
                                       Constraints& constraints,
                                       Option_AutoFit &mauto)
{
  int icov, ncova, ivar, jvar, jmod, kcov, imod, nmodel, ncovleft, rank;
  int flag_modified, flag_range, flag_param, min_order, max_ndim, flag_int_1d;
  int flag_int_2d, flag_aniso, flag_rotation;
  int lost_rank, lost_imod, lost_icov;
  double scalfac, parmax;
  VectorInt flag_compress;
  Model *model;
  EConsElem icons;

  /* Initializations */

  int nparloc = *npar;
  int ntot = 0;
  Option_VarioFit optvar = strmod->optvar;

  /* Load the parameters in the Model */

  st_model_auto_strmod_define(strmod, nparloc, param);

  /* Run the last Goulard algorithm (if necessary) */

  st_goulard_verbose(0, mauto);
  if (optvar.getFlagGoulardUsed())
    for (imod = 0; imod < strmod->nmodel; imod++)
    {
      ST_PREPAR_GOULARD(imod);
      (void) st_goulard_fitting(1, 1, STRMOD->models[imod], constraints, mauto);
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
        if (mauto.getVerbose() || OptDbg::query(EDbg::CONVERGE))
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
          rank = st_parid_match(strmod, nparloc, imod, icov, EConsElem::ANGLE, -1, -1);
          if (rank >= 0 && lost_rank < 0)
          {
            lost_rank = rank;
            st_parid_decode(strmod->parid[lost_rank], &lost_imod, &lost_icov,
                            &icons, &ivar, &jvar);
            if (mauto.getVerbose() || OptDbg::query(EDbg::CONVERGE))
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

        if (mauto.getVerbose() || OptDbg::query(EDbg::CONVERGE))
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

  label_compress:
  nparloc = st_compress_parid(ntot, strmod->parid, param, lower, upper);

  /* Modifying the covariance ranks in parid */

  for (imod = 0; imod < strmod->nmodel; imod++)
  {
    int jcov = 0;
    for (icov = 0; icov < strmod->models[imod]->getCovaNumber(); icov++)
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

  label_end:
  *npar = nparloc;
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
  int flag_range, flag_param, flag_aniso, flag_rotation;
  int min_order, max_ndim, flag_int_1d, flag_int_2d;
  double scalfac, parmax;

  /* Loop on the basic structures */

  for (int jcov = 0; jcov < model->getCovaNumber(); jcov++)
  {
    CovAniso* cova = model->getCova(jcov);
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
  int ndim = model->getDimensionNumber();
  int ndir = vario->getDirectionNumber();
  int n_2d = 0;
  int n_3d = 0;

  /* 2-D case */
  if (ndim == 2)
  {
    n_2d = ndir;
    n_3d = 0;
  }

  /* 3-D case */
  if (ndim == 3)
  {
    for (int idir = 0; idir < ndir; idir++)
    {
      if (isZero(vario->getCodir(idir, 2)))
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

  if (ndim == 3 && n_3d <= 0) optvar.setLockNo3d(1);
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

  st_modify_optvar_for_anam(model,optvar);

  /* Case when constraints involve sill(s) */

  if (constraints.isDefinedForSill() && optvar.getFlagGoulardUsed())
  {
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

  st_modify_optvar_for_anam(model,optvar);

  /* Case when constraints involve sill(s) */

  if (constraints.isDefinedForSill() && optvar.getFlagGoulardUsed())
  {
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
  int flag_range, flag_param, flag_aniso, flag_rotation;
  int min_order, max_ndim, flag_int_1d, flag_int_2d;
  double scalfac, parmax;

  /* Loop on the two candidate models */

  int ntot = 0;
  for (int imod = 0; imod < 2; imod++)
  {
    Model* model = (imod == 0) ? model1 : model2;
    if (model == nullptr) continue;

    /* Initializations */

    int ndim = model->getDimensionNumber();
    int nvar = model->getVariableNumber();
    int first_covrot = -1;
    if (st_alter_model_optvar(vario, model, constraints, optvar)) return (-2);

    /* Define the model */

    if (st_model_define(model, optvar)) return (-1);

    /* Count the number of parameters */

    for (int jcov = 0; jcov < model->getCovaNumber(); jcov++)
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
          for (int idim = 1; idim < ndim; idim++)
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
            for (int idim = 0; idim < ndim; idim++)
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
  for (int iparam = 0; iparam < ntot; iparam++)
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
  /* Define the current values of the parameters */

  st_model_auto_strmod_define(STRMOD, npar, param);

  /* Run the Goulard algorithm (if necessary) */

  st_goulard_verbose(0, MAUTO);
  if (STRMOD->optvar.getFlagGoulardUsed())
    for (int imod = 0; imod < STRMOD->nmodel; imod++)
    {
      ST_PREPAR_GOULARD(imod);
      (void) st_goulard_fitting(1, 0, STRMOD->models[imod], CONSTRAINTS, MAUTO);
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
  MatrixSquareGeneral tab(nvar);
  db_index_sample_to_grid(DBMAP, nech / 2, INDG1);
  CovCalcMode mode(ECalcMember::LHS);
  mode.setAsVario(true);
  mode.setUnitary(true);
  mode.setOrderVario(STRMOD->norder);

  /* Loop on the basic structures */

  for (int icov = 0; icov < ncova; icov++)
  {
    mode.setActiveCovListFromOne(icov);

    /* Loop on the experiments */

    for (int ipadir = 0; ipadir < RECINT.npadir; ipadir++)
    {
      db_index_sample_to_grid(DBMAP, ipadir, INDG2);
      for (int idim = 0; idim < ndim; idim++)
        d0[idim] = (INDG2[idim] - INDG1[idim]) * DBMAP->getDX(idim);
      model->evaluateMatInPlace(nullptr, d0, tab, true, 1., &mode);

      /* Loop on the variables */

      int ijvar = 0;
      for (int ivar = 0; ivar < nvar; ivar++)
        for (int jvar = 0; jvar <= ivar; jvar++, ijvar++)
          GE(icov,ijvar,ipadir) = tab.getValue(ivar, jvar);
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
  /* Define the current values of the parameters */

  st_model_auto_strmod_define(STRMOD, npar, param);

  /* Run the Goulard algorithm (if necessary) */

  st_goulard_verbose(0, MAUTO);
  if (STRMOD->optvar.getFlagGoulardUsed())
    for (int imod = 0; imod < STRMOD->nmodel; imod++)
    {
      ST_PREPAR_GOULARD(imod);
      (void) st_goulard_fitting(1, 0, STRMOD->models[imod], CONSTRAINTS, MAUTO);
    }
  st_goulard_verbose(1, MAUTO);

  /* Calculate the array of model values */

  for (int imod = 0; imod < STRMOD->nmodel; imod++)
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
  int ndim = model->getDimensionNumber();
  int nvar = vario->getVariableNumber();

  MatrixSquareSymmetric mat(nvar);

  /* Particular case of Properties attached to the Model */

  if (model->getCovMode() != EModelProperty::NONE)
  {
    Model* model_nugget = Model::createNugget(ndim, nvar);
    MatrixSquareGeneral aux(nvar);
    model->evaluateMatInPlace(nullptr, VectorDouble(), aux);

    for (int ivar = 0; ivar < nvar; ivar++)
      for (int jvar = 0; jvar <= ivar; jvar++)
      {
        double auxval = aux.getValue(ivar, jvar);
        double value = (ABS(auxval) > 0.) ? vario->getVar(ivar, jvar) / auxval : 0.;
        mat.setValue(ivar, jvar, value);
      }
    delete model_nugget;
  }
  else
  {
    for (int ivar = 0; ivar < nvar; ivar++)
      for (int jvar = 0; jvar <= ivar; jvar++)
        mat.setValue(ivar, jvar, vario->getVar(ivar, jvar));
  }
  if (mat.choleskyDecompose())
  {
    /* The matrix is filled arbitrarily */
    for (int ivar = 0; ivar < nvar; ivar++)
      for (int jvar = 0; jvar < nvar; jvar++)
        mat.setValue(ivar, jvar, (ivar == jvar));
    if (mat.choleskyDecompose())
      messageAbort("Error in st_vario_varchol_manage(): This should never happen");
  }
  varchol = mat.getCholeskyTL();
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
  double mini, maxi, gmax;

  /* Initializations */

  int nvar = dbmap->getLocNumber(ELoc::Z);
  int size = nvar * (nvar + 1) / 2;
  int nvar2 = nvar * nvar;

  /* Allocation */

  VectorDouble aux(nvar2,0.);
  for (int ivar = 0; ivar < nvar; ivar++)
  {
    int iloc = db_attribute_identify(dbmap, ELoc::Z, ivar);
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
        if (!cova->isIsotropic())
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
  int nv2 = nvar * nvar;
  int nvs2 = nvar * (nvar + 1) / 2;

  RECINT.npadir = npadir;
  RECINT.wt.fill(TEST, npadir * nvs2);
  RECINT.gg.fill(TEST, npadir * nvs2);
  RECINT.ge.fill(TEST, npadir * nvs2 * ncova);
  RECINT.sill.fill(TEST, nvs2 * ncova);
  RECINT.covtab.resize(nv2);

  if (flag_exp)
  {
    RECINT.wtc.fill(TEST, nbexp);
    RECINT.ggc.fill(TEST, nbexp);
    RECINT.dd.fill(TEST, npadir * nvs2 * ndim);
  }

  if (mauto.getFlagIntrinsic())
  {
    RECINT.alphau.fill(TEST, ncova);
    RECINT.sill1.fill(TEST, nvs2);
    RECINT.ge1.fill(TEST, nvs2 * npadir);
    RECINT.ge2.fill(TEST, nvs2 * npadir * ncova);
    RECINT.wt2.fill(TEST, nvs2 * npadir);
    RECINT.gg2.fill(TEST, nvs2 * npadir);
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
int model_auto_fit(Vario *vario,
                   Model *model,
                   bool verbose,
                   const Option_AutoFit &mauto_arg,
                   const Constraints &cons_arg,
                   const Option_VarioFit &optvar_arg)
{
  int npar, npar0, flag_hneg, flag_gneg, flag_reduce, flag_regular;
  double c0, hmin, hmax, gmin, gmax;
  std::vector<StrExp> strexps;
  VectorDouble varchol, scale, param, lower, upper;
  static int flag_check_result = 0;

  // Getting local copy of const references

  Option_AutoFit mauto = mauto_arg;
  Option_VarioFit optvar = optvar_arg;
  Constraints constraints = cons_arg;

  /* Initializations */

  int nbexp = 0;
  int status = 0;
  int npadir = 0;
  int ncova = model->getCovaNumber();
  int ndim = model->getDimensionNumber();
  int nvar = vario->getVariableNumber();
  VectorDouble angles;
  st_regularize_init();
  mauto.setVerbose(verbose);
  StrMod* strmod = nullptr;

  /* Preliminary checks */

  int error = 1;
  int norder = 0;
  if (vario->getCalcul() == ECalcVario::GENERAL1) norder = 1;
  if (vario->getCalcul() == ECalcVario::GENERAL2) norder = 2;
  if (vario->getCalcul() == ECalcVario::GENERAL3) norder = 3;
  if (vario->getCalcul() == ECalcVario::MADOGRAM ||
      vario->getCalcul() == ECalcVario::RODOGRAM ||
      vario->getCalcul() == ECalcVario::GENERAL1 ||
      vario->getCalcul() == ECalcVario::GENERAL2 ||
      vario->getCalcul() == ECalcVario::GENERAL3)
  {
    messerr("Procedure is designed only for symmetric covariance");
    return (1);
  }
  if (constraints.isConstraintSillDefined())
  {
    if (!optvar.getFlagGoulardUsed())
     {
       messerr("When Constraints on the sum of Sills are defined");
       messerr("The Goulard option must be switched ON");
       return (1);
     }
    if (!FFFF(constraints.getConstantSillValue()))
     constraints.expandConstantSill(nvar);
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

  vario->getExtension(0, 0, -1, 0, 1, TEST, TEST, TEST, TEST, &flag_hneg,
                      &flag_gneg, &c0, &hmin, &hmax, &gmin, &gmax);
  angles.resize(ndim);
  (void) GH::rotationGetAnglesFromCodirInPlace(vario->getCodirs(0), angles);
  st_vario_varchol_manage(vario, model, varchol);

  /* Scale the parameters in the mauto structure */

  st_mauto_rescale(nvar, varchol, mauto);

  /* Create the experimental structures */

  if (st_get_vario_dimension(vario, &nbexp, &npadir)) goto label_end;
  strexps = st_strexp_manage(nbexp, ndim);

  /* Fill the weight and experimental tabulated arrays */

  if (st_manage_recint(mauto, 1, ndim, nvar, nbexp, ncova, npadir)) goto label_end;

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

  st_model_auto_pardef(strmod, npar, hmax, varchol, angles, param, lower, upper);
  st_model_auto_scldef(strmod, npar, hmax, varchol, scale);
  st_model_auto_constraints_apply(strmod, npar, constraints, param, lower, upper);

  /* Minimization algorithm */

  STREXPS = strexps;
  STRMOD = strmod;
  MAUTO = mauto;
  CONSTRAINTS = constraints;
  ST_PREPAR_GOULARD = st_prepar_goulard_vario;
  do
  {
    st_model_auto_strmod_print(1, strmod, constraints, mauto, param, lower,
                               upper, npar, nbexp);
    if (npar > 0)
    {
      status = foxleg_f(nbexp, npar, 0, MatrixRectangular(), param, lower, upper,
                        scale, mauto, 0, st_strmod_vario_evaluate, RECINT.ggc, RECINT.wtc);
    }
    else
    {
      status = st_goulard_fitting(1, 1, model, constraints, mauto);
    }
    if (status > 0) goto label_end;

    st_model_auto_strmod_print(0, strmod, constraints, mauto, param, lower,
                               upper, npar, nbexp);
    flag_reduce = st_model_auto_strmod_reduce(strmod, &npar, hmax, gmax, param,
                                              lower, upper, constraints, mauto);

    /* Reset the parameters to their default values */

    if (status < 0)
    {
      for (int i = 0; i < npar; i++)
        param[i] = lower[i] = upper[i] = TEST;
      st_model_auto_pardef(strmod, npar, hmax, varchol, angles, param, lower, upper);
    }
    st_model_auto_scldef(strmod, npar, hmax, varchol, scale);
  }
  while (flag_reduce && npar > 0);

  /* Perform the last cosmetic updates */

  st_model_post_update(strmod, optvar);

  /* Set the returned error code */

  if (mauto.getVerbose() && status < 0)
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

  label_end:
  strmod = st_model_auto_strmod_free(strmod);
  return (error);
}

/****************************************************************************/
/*!
 **  Fitting a model using an experimental variogram
 **
 ** \return  Error return code
 **
 ** \param[in]     vario       Vario structure
 ** \param[in,out] model       Model to be fitted
 ** \param[in]     constraints Constraints structure
 ** \param[in]     mauto       Option_AutoFit structure
 **
 *****************************************************************************/
int model_fitting_sills(Vario *vario,
                        Model *model,
                        const Constraints& constraints,
                        const Option_AutoFit& mauto)
{
  int nbexp, npadir;
  std::vector<StrExp> strexps;

  /*******************/
  /* Initializations */
  /*******************/

  if (model == nullptr) return (1);
  if (vario == nullptr) return (1);
  int ndir = vario->getDirectionNumber();
  int ndim = model->getDimensionNumber();
  int nvar = model->getVariableNumber();
  int ncova = model->getCovaNumber();

  /* Reset the coregionalization matrix */

  if (st_get_vario_dimension(vario, &nbexp, &npadir) ||
      nvar <= 0 || ndir <= 0 || ncova <= 0)
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

  if (st_goulard_fitting(1, 1, model, constraints, mauto)) return 1;

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
  int flag_range, flag_param, flag_aniso, flag_rotation;
  int min_order, max_ndim, flag_int_1d, flag_int_2d;
  double scalfac, parmax;

  /* Initializations */

  int ndim = model->getDimensionNumber();
  int nvar = model->getVariableNumber();
  int first_covrot = -1;

  /* Check the number of Variogram Maps */

  if (nvar * (nvar + 1) / 2 != dbmap->getLocNumber(ELoc::Z))
  {
    messerr("The number of items in the Db Grid for Variogram maps (%d)",
            dbmap->getLocNumber(ELoc::Z));
    messerr("is not compatible with the number of variables in the Model (%d)",
            nvar);
    return (-1);
  }

  if (st_alter_vmap_optvar(dbmap, model, constraints, optvar)) return (-2);

  /* Define the model */

  if (st_model_define(model, optvar)) return (-1);

  /* Count the number of parameters */

  int ntot = 0;
  for (int jcov = 0; jcov < model->getCovaNumber(); jcov++)
  {
    model_cova_characteristics(model->getCovaType(jcov), cov_name, &flag_range,
                               &flag_param, &min_order, &max_ndim, &flag_int_1d,
                               &flag_int_2d, &flag_aniso, &flag_rotation,
                               &scalfac, &parmax);
    if (max_ndim > 0 && ndim > max_ndim)
    {
      messerr("The structure '%s' is limited to dimension (%d)", cov_name, max_ndim);
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
        for (int idim = 1; idim < ndim; idim++)
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
          for (int idim = 0; idim < ndim; idim++)
            ntot++;
        }
      }
    }
  }

  /* Allocation */

  param.resize(ntot);
  lower.resize(ntot);
  upper.resize(ntot);
  for (int iparam = 0; iparam < ntot; iparam++)
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
  int nech = DBMAP->getSampleNumber();
  int nvar = DBMAP->getLocNumber(ELoc::Z);
  int nvs2 = nvar * (nvar + 1) / 2;
  db_index_sample_to_grid(DBMAP, nech / 2, INDG1);

  /* Load the Experimental conditions structure */

  int ipadir = 0;
  for (int iech = 0; iech < nech; iech++)
  {
    db_index_sample_to_grid(DBMAP, iech, INDG2);
    double dist = distance_intra(DBMAP, nech / 2, iech, NULL);
    double wgt = (dist > 0) ? 1. / dist : 0.;

    /* Check samples containing only undefined values */

    int ntest = 0;
    for (int ijvar = 0; ijvar < nvs2; ijvar++)
      if (!FFFF(DBMAP->getLocVariable(ELoc::Z,iech, ijvar))) ntest++;
    if (ntest <= 0) continue;

    for (int ijvar = 0; ijvar < nvs2; ijvar++)
    {
      WT(ijvar,ipadir)= 0.;
      GG(ijvar,ipadir) = 0.;

      double value = DBMAP->getLocVariable(ELoc::Z,iech,ijvar);
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
int vmap_auto_fit(const DbGrid *dbmap,
                  Model *model,
                  bool verbose,
                  const Option_AutoFit &mauto_arg,
                  const Constraints &cons_arg,
                  const Option_VarioFit &optvar_arg)
{
  int npar0, npar, flag_reduce;
  VectorDouble varchol, scale, param, lower, upper;

  // Copy of const referencse into local classes

  Option_AutoFit mauto = mauto_arg;
  Constraints constraints = cons_arg;
  Option_VarioFit optvar = optvar_arg;

  /* Initializations */

  int error = 1;
  int nbexp = 0;
  int status = 0;
  int norder = 0;
  int npadir = 0;
  int ncova = model->getCovaNumber();
  int nvar = model->getVariableNumber();
  int ndim = model->getDimensionNumber();
  double hmax = 0.;
  double gmax = 0.;
  StrMod* strmod = nullptr;
  VectorDouble angles(ndim, 0.);
  mauto.setVerbose(verbose);

  /* Preliminary checks */

  if (nvar != dbmap->getLocNumber(ELoc::Z))
  {
    messerr("Number of variables in Db (%d) must match the one in Model (%d)",
            model->getVariableNumber(), dbmap->getLocNumber(ELoc::Z));
    goto label_end;
  }
  if (constraints.isConstraintSillDefined())
  {
    if (!optvar.getFlagGoulardUsed())
     {
       messerr("When Constraints on the sum of Sills are defined");
       messerr("The Goulard option must be switched ON");
       goto label_end;
     }
    if (!FFFF(constraints.getConstantSillValue()))
     constraints.expandConstantSill(nvar);
  }
  if (st_get_vmap_dimension(dbmap, nvar, &npadir, &nbexp)) goto label_end;
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

  npar0 = st_vmap_auto_count(dbmap, model, constraints, optvar, param, lower, upper);

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

  st_model_auto_pardef(strmod, npar, hmax, varchol, angles, param, lower, upper);
  st_model_auto_scldef(strmod, npar, hmax, varchol, scale);
  st_model_auto_constraints_apply(strmod, npar, constraints, param, lower, upper);

  /* Minimization algorithm */

  STRMOD = strmod;
  MAUTO = mauto;
  CONSTRAINTS = constraints;
  ST_PREPAR_GOULARD = st_prepar_goulard_vmap;
  do
  {
    st_model_auto_strmod_print(1, strmod, constraints, mauto, param, lower,
                               upper, npar, nbexp);
    status = foxleg_f(nbexp, npar, 0, MatrixRectangular(), param, lower, upper,
                      scale, mauto, 0, st_strmod_vmap_evaluate, RECINT.gg, RECINT.wt);
    if (status > 0) goto label_end;

    st_model_auto_strmod_print(0, strmod, constraints, mauto, param, lower,
                               upper, npar, nbexp);
    flag_reduce = st_model_auto_strmod_reduce(strmod, &npar, hmax, gmax, param,
                                              lower, upper, constraints, mauto);

    /* Reset the parameters to their default values */

    if (status < 0)
    {
      for (int i = 0; i < npar; i++)
        param[i] = lower[i] = upper[i] = TEST;
      st_model_auto_pardef(strmod, npar, hmax, varchol, angles, param, lower, upper);
    }
    st_model_auto_scldef(strmod, npar, hmax, varchol, scale);
  }
  while (flag_reduce && npar > 0);

  /* Perform the last cosmetic updates */

  st_model_post_update(strmod, optvar);

  /* Set the returned error code */

  if (mauto.getVerbose() && status < 0)
    messerr("Convergence not reached after %d iterations (%d parameters)",
            mauto.getMaxiter(), npar);

  /* Store the sills in the keypair mechanism */

  st_keypair_sill(1, model);

  error = 0;

  label_end:
  INDG1 = db_indg_free(INDG1);
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
  int ncons = (int) constraints.getConsItemNumber();
  for (int i = 0; i < ncons; i++)
  {
    const ConsItem *consitem = constraints.getConsItems(i);
    if (consitem->getType() != EConsElem::SILL) continue;
    if (consitem->getValue() < 0) return (1);
    constraints.setValue(i, sqrt(consitem->getValue()));

    // For constraints on the Sill in monovariate case,
    // Add a constraints on AIC for lower bound
    if (consitem->getIV1() == 0 &&
        consitem->getIV2() == 0 &&
        consitem->getIcase() == EConsType::UPPER)
    {
      ConsItem* consjtem = new ConsItem(*consitem);
      consjtem->setValue(-consjtem->getValue());
      consjtem->setIcase(EConsType::LOWER);
      constraints.addItem(consjtem);
    }
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
