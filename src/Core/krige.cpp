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
#include "geoslib_f.h"
#include "geoslib_old_f.h"
#include "geoslib_define.h"

#include "Polynomials/Hermite.hpp"
#include "Db/ELoadBy.hpp"
#include "Db/Db.hpp"
#include "Db/DbStringFormat.hpp"
#include "Model/Model.hpp"
#include "Model/CovInternal.hpp"
#include "Neigh/ANeighParam.hpp"
#include "Neigh/NeighMoving.hpp"
#include "Neigh/NeighImage.hpp"
#include "Neigh/NeighUnique.hpp"
#include "Neigh/NeighWork.hpp"
#include "Anamorphosis/AnamDiscreteDD.hpp"
#include "Anamorphosis/AnamDiscreteIR.hpp"
#include "Anamorphosis/AnamHermite.hpp"
#include "Basic/String.hpp"
#include "Basic/NamingConvention.hpp"
#include "Basic/Utilities.hpp"
#include "Basic/Law.hpp"
#include "Basic/File.hpp"
#include "Basic/OptDbg.hpp"
#include "Covariances/CovLMCAnamorphosis.hpp"
#include "Covariances/CovContext.hpp"
#include "Covariances/ECalcMember.hpp"
#include "Estimation/KrigingSystem.hpp"
#include "Anamorphosis/EAnam.hpp"

#include <math.h>
#include <string.h>
#include <stdio.h>
#include <limits.h>

/*! \cond */
#define NBYPAS 5
#define IND(iech,ivar)    ((iech) + (ivar) * nech)
#define FLAG(iech,ivar)   (flag[IND(iech,ivar)])
#define COVTAB(ivar,jvar) (covtab[(jvar) * nvar_m + (ivar)])
#define COVAUX(ivar,jvar) (covaux[(jvar) * nvar_m + (ivar)])
#define RHS(i,iv,jv)      (rhs [IND(i,iv) + neq * (jv)])
#define FF0(ib,iv)        (ff0 [(ib) + nfeq * (iv)])
#define VAR0(iv,jv)       (var0[(jv) + nvar * (iv)])
#define VARB(iv,jv)       (varb[(jv) + nvar * (iv)])
#define LHS(i,iv,j,jv)    (lhs [IND(i,iv) + neq * IND(j,jv)])
#define LHS_C(i,j)        (lhs [(i) + nred * (j)])
#define LHS_B(i,j)        (lhs_b[(i) + nred * (j)])
#define LHS_LOC(i,j)      (lhs_loc [(i) + (nloc+nfeq) * (j)])
#define RHS_C(i,iv)       (rhs [(i) + nred * (iv)])
#define WGT(i,iv)         (wgt [(i) + nred * (iv)])
#define ZAM1(i)           (zam1[(i)])
#define DISC1(i,idim)     (KOPTION->disc1[(idim) * KOPTION->ntot + (i)])
#define DISC2(i,idim)     (KOPTION->disc2[(idim) * KOPTION->ntot + (i)])
#define FS(ib,il)         (fs [(il) * shift + (ib)])
#define FSF(ib,jb)        (fsf[(jb) * shift + (ib)])
#define FF(ib,il)         (ff [(il) * shift + (ib)])
#define MU(ib)            (mu[(ib)])
#define SMU(ib)           (smu[(ib)])
#define FCA(ib,il)        (fca[(il) * shift + (ib)])
#define FSFPCM1(ib,jb)    (fsfpcm1[(jb)    * shift + (ib)])
#define SIGMA(ib,jb)      (sigma[(jb)    * shift + (ib)])
#define DCOV(il,jl)       (dcov[(jl)   * nfeq  + (il)])
#define RCOV(il,jl)       (rcov[(jl)   * nfeq  + (il)])
#define LHS_EXP(i,j)      (lhs[(i) * neq + (j)])
#define RHS_EXP(i)        (rhs[(i)])
#define COV_REF(iz)       (cov_ref[cov_radius + (iz)])
#define SMEAN(i,isimu)    (smean[(isimu) * nfeq + (i)])
#define IAD(ix,iy,iz,nn,ss) (((iz) + nn[2]) + ss[2] * (((iy) + nn[1]) + ss[1] * ((ix) + nn[0])))
#define COV_RES(ix,iy,iz) cov_res[IAD(ix,iy,iz,cov_nn,cov_ss)]
#define COV_TOT(ix,iy,iz) cov_tot[IAD(ix,iy,iz,cov_nn,cov_ss)]
#define NUM_TOT(ix,iy,iz) num_tot[IAD(ix,iy,iz,cov_nn,cov_ss)]
#define NEI_CUR(ix,iy,iz) nei_cur[IAD(ix,iy,iz,nei_nn,nei_ss)]
#define NEI_REF(ix,iy,iz) nei_ref[IAD(ix,iy,iz,nei_nn,nei_ss)]
#define LTERM(ivar,iz)    lterm[(iz) * nvarin + ivar]
#define LBACK(ivar,iz)    lback[(iz) * nvarin + ivar]
#define CC(iz,jz)         (cc[(jz)*((jz) + 1) / 2 +(iz)])
#define PROPTAB(ivar,iz)  (proptab[(iz) * nvarin + (ivar)])
#define UTAB(i,j)         (utab[(i) + ndat * (j)])
#define TUTIL(i,j)        (tutil[(i) + nutil * (j)])
#define SPART(i,j)        (spart[(i) + npart * (j)])
#define MATCL(i,j)        (matCL[i][j])
#define COVSS(is,js)      (covss  [(js) + ns * (is)])
#define COVGEN(i1,i2)     (covgen [(i2) + n2 * (i1)])
#define COVPP(ip,jp)      (covpp  [(jp) + np * (ip)])
#define COVGP(ig,ip)      (covgp  [(ip) + np * (ig)])
#define DISTGEN(i,is)     (distgen[(is) + ns * (i)])
#define DISTPS(ip,is)     (distps [(is) + ns * (ip)])
#define DISTGS(ig,is)     (distgs [(is) + ns * (ig)])
#define PRODGEN(i,is)     (prodgen[(is) + ns * (i)])
#define PRODPS(ip,is)     (prodps [(is) + ns * (ip)])
#define PRODGS(ig,is)     (prodgs [(is) + ns * (ig)])
#define DRFTAB(ip,il)     (drftab [(il) + nbfl * (ip)])
#define YMAT(ip,il)       (ymat   [(il) + nbfl * (ip)])
/*! \endcond */

// TODO : remove all these static stuffs !
static double *covtab, *covaux, *drftab, *fs, *fsf, *d1_1, *d1_2, *var0;
static VectorDouble d1, d1_t;
static double *lhs, *lhs_b, *rhs, *wgt, *zext, *zam1, *ff0, *varb;
static int *flag;
static int KRIGE_INIT = 0;
static int MODEL_INIT = 0;
static int IECH_OUT   = -1;
static int RAND_INDEX = -1;
static int FLAG_EST, FLAG_STD, FLAG_WGT, FLAG_COLK, FLAG_SIMU, FLAG_LTERM;
static int FLAG_BAYES, FLAG_PROF, FLAG_VARZ, FLAG_DGM;
static int IPTR_EST, IPTR_STD, IPTR_VARZ, IPTR_NBGH;
static int *RANK_COLCOK;
static Db *DBIN, *DBOUT;
static Koption *KOPTION;
static double R_COEFF;
static int INH_FLAG_VERBOSE = 0;
static int INH_FLAG_LIMIT = 1;
static char string[100];

static CovInternal COVINT;

typedef struct
{
  int ndtot;
  int rank1;
  int rank2;
  Model *model;
  int nugget_opt;
  int nostd;
  ECalcMember member;
  int icov_r;
  double weight;
} Disc_Structure;

/****************************************************************************/
/*!
 **  Management of internal array (double)
 **
 ** \return  Pointer to the newly allocated array
 **
 ** \param[in]  nli   Number of lines
 ** \param[in]  nco   Number of columns
 **
 *****************************************************************************/
static double* st_core(int nli, int nco)
{
  double *tab, rsize;
  int size, i;

  /* Initialization */

  tab = nullptr;
  rsize = (double) nli * (double) nco;
  if (rsize < 0 || rsize > INT_MAX)
  {
    messerr("Core allocation problem: Size (%d x %d) too big", nli, nco);
    return (tab);
  }
  size = nli * nco;

  /* Allocation */

  tab = (double*) mem_alloc(sizeof(double) * size, 0);
  if (tab == nullptr)
  {
    messerr("Core allocation problem: Size (%d) too big", size);
    return (tab);
  }

  for (i = 0; i < size; i++)
    tab[i] = 0.;
  return (tab);
}

/****************************************************************************/
/*!
 **  Management of internal array (integer)
 **
 ** \return  Pointer to the newly allocated array
 **
 ** \param[in]  nli   Number of lines
 ** \param[in]  nco   Number of columns
 **
 *****************************************************************************/
static int* st_icore(int nli, int nco)
{
  int *tab, size, i;

  /* Initialization */

  tab = nullptr;
  size = nli * nco;

  /* Allocation */

  tab = (int*) mem_alloc(sizeof(int) * size, 0);
  if (tab != nullptr) for (i = 0; i < size; i++)
    tab[i] = 0;

  return (tab);
}

/****************************************************************************/
/*!
 **  Manage the relative position array
 **
 ** \param[in]  mode    1 for creating, -1 for deleting
 ** \param[in]  neq     Number of kriging equations
 ** \param[in]  rel_arg Relative position array (used for deletion)
 **
 *****************************************************************************/
static int* st_relative_position_array(int mode, int neq, int *rel_arg)
{
  int *rel, i, j;

  /* Dispatch */

  if (mode > 0)
  {

    /* Creation */

    rel = (int*) st_icore(neq, 1);
    if (rel == nullptr) return (rel);
    for (i = j = 0; i < neq; i++)
    {
      if (flag != NULL && flag[i])
        rel[j++] = i + 1;
      else
        rel[j++] = i + 1;
    }
  }
  else
  {

    /* Deletion */

    rel = rel_arg;
    rel = (int*) mem_free((char* ) rel);
  }
  return (rel);
}

/****************************************************************************/
/*!
 **  Initialize the static global variables
 **
 ** \param[in]  dbin   input Db structure
 ** \param[in]  dbout  output Db structure
 **
 *****************************************************************************/
static void st_global_init(Db *dbin, Db *dbout)
{
  FLAG_STD = FLAG_EST = FLAG_WGT = FLAG_LTERM = FLAG_VARZ = 0;
  FLAG_COLK = FLAG_BAYES = FLAG_PROF = FLAG_SIMU = FLAG_DGM = 0;
  IPTR_EST = IPTR_STD = IPTR_VARZ = IPTR_NBGH = 0;
  IECH_OUT = 0;

  /* Set the global variables */

  DBIN = dbin;
  DBOUT = dbout;

  /* Change of support coefficient for DGM */

  R_COEFF = 1.;
  COVINT = CovInternal();

  return;
}

/****************************************************************************/
/*!
 **  Returns the coordinate of the data (at rank if rank >= 0)
 **  or of the target (at IECH_OUT if rank < 0)
 **
 ** \param[in]  loc_rank   Rank of the sample
 ** \param[in]  idim   Rank of the coordinate
 **
 *****************************************************************************/
static double st_get_idim(int loc_rank, int idim)
{
  double value;

  if (loc_rank >= 0)
  {
    value = DBIN->getCoordinate(loc_rank, idim);
  }
  else
  {
    value = DBOUT->getCoordinate(IECH_OUT, idim);
  }
  return (value);
}

/****************************************************************************/
/*!
 **  Calculate the covariance between two samples from two Db
 **
 ** \param[in]  model        Model structure
 ** \param[in]  flag_init    Initialize the array beforehand
 ** \param[in]  nugget_opt   Option for the nugget effect basic structure
 ** \li                       0 : no particular option
 ** \li                       1 : discard the nugget effect
 ** \li                      -1 : only consider the nugget effect
 ** \param[in]  nostd        0 standard; +-1 special; ITEST normalized
 ** \param[in]  member       Member of the Kriging System (ECalcMember)
 ** \param[in]  icov_r       rank of the target covariance or -1 for all
 ** \param[in]  weight       Weight attached to this calculation
 ** \param[in]  rank1        Rank of the first sample
 ** \param[in]  rank2        Rank of the second sample
 **
 ** \param[out] d1loc        Working array
 ** \param[out] covtab_loc   Output covariance array
 **
 *****************************************************************************/
static void st_cov(Model *model,
                   int flag_init,
                   int nugget_opt,
                   int nostd,
                   const ECalcMember &member,
                   int icov_r,
                   double weight,
                   int rank1,
                   int rank2,
                   VectorDouble d1loc,
                   double *covtab_loc)
{

  /* Initializations */

  if (rank1 >= 0)
  {
    COVINT.setDb1(DBIN);
    COVINT.setIcas1(1);
    COVINT.setIech1(rank1);
  }
  else
  {
    COVINT.setDb1(DBOUT);
    COVINT.setIcas1(2);
    COVINT.setIech1(IECH_OUT);
  }

  if (rank2 >= 0)
  {
    COVINT.setDb2(DBIN);
    COVINT.setIcas2(1);
    COVINT.setIech2(rank2);
  }
  else
  {
    COVINT.setDb2(DBOUT);
    COVINT.setIcas2(2);
    COVINT.setIech2(IECH_OUT);
  }

  CovCalcMode mode(member);
  mode.update(nugget_opt, nostd, member, icov_r, 0, 1);
  model_calcul_cov(&COVINT, model, mode, flag_init, weight, d1loc, covtab_loc);
}

/****************************************************************************/
/*!
 **  Internal recursive function for calculating covariance between data
 **  and data, when data are discretized
 **
 ** \param[in]  idim   Space dimension for current iteration (first point)
 ** \param[in]  jdim   Space dimension for current iteration (second point)
 ** \param[in]  it     Pointer to the Internal Disc_Structure
 **
 *****************************************************************************/
static void st_data_discretize_dd(int idim, int jdim, Disc_Structure *it)
{
  double exts2, dsize, decal;

  // Initialization

  if (idim < it->model->getDimensionNumber() - 1)
  {
    idim = idim + 1;

    // Loop in the current dimension

    exts2 = DBIN->getBlockExtension(it->rank1, idim) / 2.;
    dsize = KOPTION->dsize[idim];

    if (exts2 <= 0. || dsize <= 0.)
    {

      /* Punctual support */

      d1_1[idim] = 0.;
      st_data_discretize_dd(idim, jdim, it);
    }
    else
    {

      /* Implicit loop until reaching the edge of the data */

      decal = -exts2 + dsize / 2.;
      do
      {
        d1_1[idim] = decal;
        st_data_discretize_dd(idim, jdim, it);
        decal = decal + dsize;
      }
      while (decal < exts2);
    }
  }
  else if (jdim < it->model->getDimensionNumber() - 1)
  {
    jdim = jdim + 1;

    // Loop in the current dimension

    exts2 = DBIN->getBlockExtension(it->rank2, jdim) / 2.;
    dsize = KOPTION->dsize[jdim];

    if (exts2 <= 0 || dsize <= 0.)
    {
      /* Punctual support */

      d1_2[jdim] = 0.;
      st_data_discretize_dd(idim, jdim, it);
    }
    else
    {

      /* Implicit loop until reaching the edge of the data */

      decal = -exts2 + dsize / 2.;
      do
      {
        d1_2[jdim] = decal;
        st_data_discretize_dd(idim, jdim, it);
        decal = decal + dsize;
      }
      while (decal < exts2);
    }
  }
  else
  {

    // End of implicit loop on dimensions

    it->ndtot++;
    for (int i = 0; i < it->model->getDimensionNumber(); i++)
      d1_t[i] = d1[i] + d1_1[i] + d1_2[i];
    st_cov(it->model, 0, it->nugget_opt, it->nostd, it->member, it->icov_r,
           it->weight, it->rank1, it->rank2, d1_t, covaux);
  }
}

int is_flag_data_disc_defined(void)
{
  return KOPTION->flag_data_disc;
}

/****************************************************************************/
/*!
 **  Calculate the covariance between two Data samples
 **  This covariance takes possible Data discretization into account
 **
 ** \param[in]  model  Model structure
 ** \param[in]  nugget_opt   Option for the nugget effect basic structure
 ** \li                       0 : no particular option
 ** \li                       1 : discard the nugget effect
 ** \li                      -1 : only consider the nugget effect
 ** \param[in]  nostd        0 standard; +-1 special; ITEST normalized
 ** \param[in]  icov_r       rank of the target covariance or -1 for all
 ** \param[in]  weight       Weight attached to this calculation
 ** \param[in]  rank1        Rank of the first sample
 ** \param[in]  rank2        Rank of the second sample
 **
 ** \param[out] d1           Working array
 ** \param[out] covtab       Output covariance array
 **
 *****************************************************************************/
static void st_cov_dd(Model *model,
                      int nugget_opt,
                      int nostd,
                      int icov_r,
                      double weight,
                      int rank1,
                      int rank2,
                      VectorDouble d1,
                      double *covtab)
{
  int nvar_m;
  double scale;
  Disc_Structure int_disc;

  if (!KOPTION->flag_data_disc)
  {

    // Data is considered as punctual

    st_cov(model, 0, nugget_opt, nostd, ECalcMember::LHS, icov_r, weight, rank1,
           rank2, d1, covtab);
  }
  else
  {
    // Data have a support

    nvar_m = model->getVariableNumber();
    model_covtab_init(1, model, covaux);

    // Implicit discretization loop for first datum

    int_disc.rank1 = rank1;
    int_disc.rank2 = rank2;
    int_disc.model = model;
    int_disc.nugget_opt = nugget_opt;
    int_disc.nostd = nostd;
    int_disc.member = ECalcMember::LHS;
    int_disc.icov_r = icov_r;
    int_disc.weight = weight;
    int_disc.ndtot = 0;

    st_data_discretize_dd(-1, -1, &int_disc);

    // Normalization

    scale = (double) int_disc.ndtot;
    for (int ivar_m = 0; ivar_m < nvar_m; ivar_m++)
      for (int jvar_m = 0; jvar_m < nvar_m; jvar_m++)
        COVTAB(ivar_m,jvar_m) += COVAUX(ivar_m,jvar_m) / scale;
  }
}

/****************************************************************************/
/*!
 **  Internal recursive function for calculating covariance between data
 **  and target, when data is discretized
 **
 ** \param[in]  idim   Space dimension for current iteration
 ** \param[in]  it     Pointer to the Internal Disc_Structure
 **
 *****************************************************************************/
static void st_data_discretize_dg(int idim, Disc_Structure *it)
{
  double exts2, dsize, decal;

  // Initialization

  if (idim < it->model->getDimensionNumber() - 1)
  {
    idim = idim + 1;

    // Loop in the current dimension

    exts2 = DBIN->getBlockExtension(it->rank1, idim) / 2.;
    dsize = KOPTION->dsize[idim];

    if (exts2 <= 0. || dsize <= 0.)
    {

      /* Punctual support */

      d1_1[idim] = 0.;
      st_data_discretize_dg(idim, it);
    }
    else
    {

      /* Implicit loop until reaching the edge of the data */

      decal = -exts2 + dsize / 2.;
      do
      {

        d1_1[idim] = decal;
        st_data_discretize_dg(idim, it);
        decal = decal + dsize;
      }
      while (decal < exts2);
    }
  }
  else
  {

    // End of implicit loop on dimensions

    it->ndtot++;
    for (int i = 0; i < it->model->getDimensionNumber(); i++)
      d1_t[i] = d1[i] + d1_1[i];
    st_cov(it->model, 0, it->nugget_opt, it->nostd, it->member, it->icov_r,
           it->weight, it->rank1, it->rank2, d1_t, covaux);
  }
}

/****************************************************************************/
/*!
 **  Calculate the covariance between a Data sample and a Target sample
 **  This covariance takes possible Data discretization into account
 **
 ** \param[in]  model        Model structure
 ** \param[in]  nugget_opt   Option for the nugget effect basic structure
 ** \li                       0 : no particular option
 ** \li                       1 : discard the nugget effect
 ** \li                      -1 : only consider the nugget effect
 ** \param[in]  nostd        0 standard; +-1 special; ITEST normalized
 ** \param[in]  member       Member of the Kriging System (ECalcMember)
 ** \param[in]  icov_r       rank of the target covariance or -1 for all
 ** \param[in]  weight       Weight attached to this calculation
 ** \param[in]  rank1        Rank of the first sample
 ** \param[in]  rank2        Rank of the second sample
 **
 ** \param[out] d1           Working array
 ** \param[out] covtab       Output covariance array
 **
 ** \remarks: Note the initialization of 'covtab' must be performed
 ** \remarks; in the calling function (if necessary)
 **
 *****************************************************************************/
static void st_cov_dg(Model *model,
                      int nugget_opt,
                      int nostd,
                      const ECalcMember &member,
                      int icov_r,
                      double weight,
                      int rank1,
                      int rank2,
                      VectorDouble d1,
                      double *covtab)
{
  int nvar_m;
  double scale;
  Disc_Structure int_disc;

  if (!KOPTION->flag_data_disc)
  {

    // Data is considered as ponctual

    st_cov(model, 0, nugget_opt, nostd, member, icov_r, weight, rank1, rank2,
           d1, covtab);
  }
  else
  {
    // Data have a support

    nvar_m = model->getVariableNumber();
    model_covtab_init(1, model, covaux);

    // Implicit discretization loop on the data

    int_disc.rank1 = rank1;
    int_disc.rank2 = rank2;
    int_disc.model = model;
    int_disc.nugget_opt = nugget_opt;
    int_disc.nostd = nostd;
    int_disc.member = member;
    int_disc.icov_r = icov_r;
    int_disc.weight = weight;
    int_disc.ndtot = 0;

    for (int i = 0; i < model->getDimensionNumber(); i++)
      d1_1[i] = d1[i];
    st_data_discretize_dg(-1, &int_disc);

    // Normalization

    scale = (double) int_disc.ndtot;
    for (int ivar_m = 0; ivar_m < nvar_m; ivar_m++)
      for (int jvar_m = 0; jvar_m < nvar_m; jvar_m++)
        COVTAB(ivar_m,jvar_m) += COVAUX(ivar_m,jvar_m) / scale;
  }
}

/****************************************************************************/
/*!
 **  Returns the value of the variable (at rank if rank >= 0)
 **  or of the target (at IECH_OUT if rank < 0)
 **
 ** \param[in]  rank   Rank of the sample
 ** \param[in]  ivar   Rank of the variable
 **
 ** \remarks   In case of simulation, the variable of the first simulation
 ** \remarks   is systematically returned. This has no influence on the rest
 ** \remarks   of the calculations
 **
 *****************************************************************************/
static double st_get_ivar(int rank, int ivar)
{
  double value;
  int jvar;

  if (rank >= 0)
  {

    // Variable in the Input file

    if (!FLAG_SIMU)

      // Particular case of simulations

      value = DBIN->getVariable(rank, ivar);
    else

      // Case of the traditional kriging based on Z-variables

      value = DBIN->getSimvar(ELoc::SIMU, rank, 0, ivar, 0, 1, 0);
  }
  else
  {

    // Variable in the Output file: colocated case

    jvar = RANK_COLCOK[ivar];
    if (jvar < 0)
      value = TEST;
    else
      value = DBOUT->getArray(IECH_OUT, jvar);
  }

  return (value);
}

/****************************************************************************/
/*!
 **  Returns the value of the measurement error (at rank if rank >= 0)
 **  or of the target (at IECH_OUT if rank < 0)
 **
 ** \param[in]  rank   Rank of the sample
 ** \param[in]  ivar   Rank of the variable
 **
 *****************************************************************************/
static double st_get_verr(int rank, int ivar)
{
  double value;

  if (rank >= 0)
  {
    value = DBIN->getVarianceError(rank, ivar);
  }
  else
  {
    value = DBOUT->getVarianceError(IECH_OUT, ivar);
  }
  return (value);
}

/****************************************************************************/
/*!
 **  Returns the value of the external drift "rank" (if rank >= 0)
 **  or of the target (at IECH_OUT if rank < 0)
 **
 ** \param[in]  rank   Rank of the sample
 ** \param[in]  ibfl   Rank of the external drift
 **
 *****************************************************************************/
static double st_get_fext(int rank, int ibfl)
{
  double value;

  if (rank >= 0)
  {
    value = DBIN->getExternalDrift(rank, ibfl);
  }
  else
  {
    value = DBOUT->getExternalDrift(IECH_OUT, ibfl);
  }
  return (value);
}

/****************************************************************************/
/*!
 **  Returns the value of the array (at rank if rank >= 0)
 **  or of the target (at IECH_OUT if rank < 0)
 **
 ** \param[in]  rank   Rank of the sample
 ** \param[in]  isimu  Rank of the simulation
 ** \param[in]  ivar   Rank of the variable
 ** \param[in]  icase  Rank of the PGS and GRF (or -1)
 ** \param[in]  nbsimu Number of simulations
 ** \param[in]  nvar   Number of variables
 **
 *****************************************************************************/
static double st_get_array(int rank,
                           int isimu,
                           int ivar,
                           int icase,
                           int nbsimu,
                           int nvar)
{
  double value;

  if (rank >= 0)
    value = DBIN->getSimvar(ELoc::SIMU, rank, isimu, ivar, icase, nbsimu, nvar);
  else
    value = DBOUT->getSimvar(ELoc::SIMU, IECH_OUT, isimu, ivar, icase, nbsimu,
                             nvar);
  return (value);
}

/****************************************************************************/
/*!
 **  Checks the kriging environment
 **
 ** \return  Error return code
 **
 ** \param[in]  flag_in    1 if the Input Db is used
 ** \param[in]  flag_out   1 if the Output Db is used
 ** \param[in]  model      Model structure (optional)
 ** \param[in]  neighparam ANeighParam structure (optional)
 **
 ** \remarks The address of the argument 'neigh' is memorized in a local
 ** \remarks static variable
 **
 *****************************************************************************/
static int st_check_environment(int flag_in,
                                int flag_out,
                                Model *model,
                                ANeighParam *neighparam)
{
  double *dbin_mini, *dbin_maxi, *dbout_mini, *dbout_maxi;
  int error, ndim, nvar, nfex;

  /* Initializations */

  error = 1;
  nvar = ndim = nfex = 0;
  dbin_mini = dbin_maxi = dbout_mini = dbout_maxi = nullptr;

  /*********************************/
  /* Compatibility between two Dbs */
  /*********************************/

  ndim = 0;
  if (flag_in && ndim == 0) ndim = DBIN->getNDim();
  if (flag_out && ndim == 0) ndim = DBOUT->getNDim();
  if (flag_in && flag_out && !DBIN->hasSameDimension(DBOUT)) goto label_end;

  /**********************/
  /* Checking the model */
  /**********************/

  if (model != nullptr)
  {
    nvar = model->getVariableNumber();
    if (nvar <= 0)
    {
      messerr("The number of variables must be positive = %d",
              model->getVariableNumber());
      goto label_end;
    }
    // The following test is avoided in the case of simulations
    // as there may be no Z-variable defined as this stage (Gibbs)
    if (flag_in && !FLAG_SIMU && DBIN->getVariableNumber() != nvar)
    {
      messerr("The number of variables of the Data (%d)",
              DBIN->getVariableNumber());
      messerr("does not match the number of variables of the Model (%d)", nvar);
      goto label_end;
    }
    if (model->getCovaNumber() <= 0)
    {
      messerr("The number of covariance must be positive");
      goto label_end;
    }
    if (model->getDimensionNumber() <= 0)
    {
      messerr("The Space Dimension must be positive = %d",
              model->getDimensionNumber());
      goto label_end;
    }
    if (model->getDimensionNumber() != ndim)
    {
      messerr("The Space Dimension of the Db structure (%d)", ndim);
      messerr("Does not correspond to the Space Dimension of the model (%d)",
              model->getDimensionNumber());
      goto label_end;
    }

    // External drifts
    nfex = model_nfex(model);
    if (nfex > 0)
    {
      if (flag_out && DBOUT->getExternalDriftNumber() != nfex)
      {
        messerr("The Model requires %d external drift(s)", model_nfex(model));
        messerr("but the output Db refers to %d external drift variables",
                DBOUT->getExternalDriftNumber());
        goto label_end;
      }

      if (flag_in && DBIN->getExternalDriftNumber() != nfex)
      {
        if (!(flag_out && is_grid(DBOUT)))
        {
          messerr("The Model requires %d external drift(s)", model_nfex(model));
          messerr("but the input Db refers to %d external drift variables",
                  DBIN->getExternalDriftNumber());
          goto label_end;
        }
      }
    }
  }

  /*********************************/
  /* Calculate the field extension */
  /*********************************/

  if (model != nullptr)
  {
    VectorDouble db_mini;
    VectorDouble db_maxi;
    db_mini.resize(ndim,TEST);
    db_maxi.resize(ndim,TEST);

    /* Input Db structure */

    if (flag_in)
      db_extension(DBIN, db_mini, db_maxi, true);

    /* Output Db structure */

    if (flag_out)
      db_extension(DBOUT, db_mini, db_maxi, true);

    model->setField(ut_vector_extension_diagonal(db_mini, db_maxi));
  }

  /*****************************/
  /* Checking the Neighborhood */
  /*****************************/

  if (neighparam != nullptr)
  {
    if (neighparam->getNDim() != ndim)
    {
      messerr("The Space Dimension of the Neighborhood (%d)", neighparam->getNDim());
      messerr("does not correspond to the Space Dimension of the first Db (%d)",
              ndim);
      goto label_end;
    }
    if (neighparam->getType() == ENeigh::IMAGE && (!flag_out || !is_grid(DBOUT)))
    {
      messerr(
          "The Image neighborhood can only be used when the output Db is a grid");
      goto label_end;
    }
  }

  /* Set the error return code */

  error = 0;

  label_end: return (error);
}

/****************************************************************************/
/*!
 **  Management of internal arrays used by cov and drift functions
 **
 ** \return  Error return code
 **
 ** \param[in]  mode   1 for allocation; -1 for deallocation
 ** \param[in]  model  Model structure
 **
 ** \remarks  This function manages covariance internal arrays with dimension
 ** \remarks  equal to the number of variables in the Model
 **
 *****************************************************************************/
static int st_model_manage(int mode, Model *model)

{
  int nvar, nbfl;

  /* Initializations */

  nvar = model->getVariableNumber();
  nbfl = model->getDriftNumber();

  /* Dispatch */

  if (mode == 1)
  {

    /* Allocation */

    if (MODEL_INIT) return (1);
    d1.resize(DBIN->getNDim());
    d1_1 = db_sample_alloc(DBIN, ELoc::X);
    if (d1_1 == nullptr) return (1);
    d1_2 = db_sample_alloc(DBIN, ELoc::X);
    if (d1_2 == nullptr) return (1);
    d1_t.resize(DBIN->getNDim());
    covtab = st_core(nvar, nvar);
    if (covtab == nullptr) return (1);
    covaux = st_core(nvar, nvar);
    if (covaux == nullptr) return (1);
    if (nbfl > 0)
    {
      drftab = st_core(nbfl, 1);
      if (drftab == nullptr) return (1);
    }
    MODEL_INIT = 1;
  }
  else
  {
    if (!MODEL_INIT) return (1);
    d1_1 = db_sample_free(d1_1);
    d1_2 = db_sample_free(d1_2);
    covtab = (double*) mem_free((char* ) covtab);
    covaux = (double*) mem_free((char* ) covaux);
    drftab = (double*) mem_free((char* ) drftab);
    MODEL_INIT = 0;
  }
  return (0);
}

/****************************************************************************/
/*!
 **  Management of internal arrays used by kriging procedure
 **
 ** \return  Error return code
 **
 ** \param[in]  mode   1 for allocation; -1 for deallocation
 ** \param[in]  nech   Number of samples in the Input Db
 **                    (only used for the neighborhood search, if any)
 ** \param[in]  nmax   Maximum number of samples per neighborhood
 ** \param[in]  nvar   Number of variables (to be calculated)
 ** \param[in]  nfeq   Number of drift equations
 **
 *****************************************************************************/
static int st_krige_manage_basic(int mode,
                                 int nech,
                                 int nmax,
                                 int nvar,
                                 int nfeq)
{
  int neqmax, ncmax;

  /* Initializations */

  ncmax = nmax * nvar;
  neqmax = ncmax + nfeq;
  if (FLAG_COLK) neqmax += nvar;
  if (FLAG_COLK) nech += 1;

  /* Dispatch */

  if (mode == 1)
  {

    /* Allocation */

    if (KRIGE_INIT) return (1);
    flag = st_icore(neqmax, 1);
    if (flag == nullptr) return (1);
    lhs = st_core(neqmax, neqmax);
    if (lhs == nullptr) return (1);
    lhs_b = st_core(neqmax, neqmax);
    if (lhs_b == nullptr) return (1);
    rhs = st_core(neqmax, nvar);
    if (rhs == nullptr) return (1);
    zext = st_core(neqmax, 1);
    if (zext == nullptr) return (1);
    zam1 = st_core(neqmax, 1);
    if (zam1 == nullptr) return (1);
    wgt = st_core(neqmax, nvar);
    if (wgt == nullptr) return (1);
    var0 = st_core(nvar, nvar);
    if (var0 == nullptr) return (1);
    if (FLAG_BAYES)
    {
      fsf = st_core(ncmax, ncmax);
      if (fsf == nullptr) return (1);
      varb = st_core(nvar, nvar);
      if (varb == nullptr) return (1);
      if (nfeq > 0)
      {
        fs = st_core(ncmax, nfeq);
        if (fs == nullptr) return (1);
        ff0 = st_core(nfeq, nvar);
        if (ff0 == nullptr) return (1);
      }
    }
    KRIGE_INIT = 1;
  }
  else
  {

    /* Deallocation */

    if (!KRIGE_INIT) return (1);
    flag = (int*) mem_free((char* ) flag);
    lhs = (double*) mem_free((char* ) lhs);
    lhs = (double*) mem_free((char* ) lhs);
    lhs_b = (double*) mem_free((char* ) lhs_b);
    rhs = (double*) mem_free((char* ) rhs);
    zext = (double*) mem_free((char* ) zext);
    zam1 = (double*) mem_free((char* ) zam1);
    wgt = (double*) mem_free((char* ) wgt);
    var0 = (double*) mem_free((char* ) var0);
    if (FLAG_BAYES)
    {
      fs = (double*) mem_free((char* ) fs);
      fsf = (double*) mem_free((char* ) fsf);
      ff0 = (double*) mem_free((char* ) ff0);
    }
    KRIGE_INIT = 0;
  }
  return (0);
}

/****************************************************************************/
/*!
 **  Returns the maximum number of points per neighborhood
 **
 ** \return  Maximum number of points per neighborhood
 **
 ** \param[in]  neighparam ANeighParam structure
 **
 *****************************************************************************/
static int st_get_nmax(ANeighParam *neighparam)
{
  return neighparam->getMaxSampleNumber(DBIN);
}

/****************************************************************************/
/*!
 **  Management of internal arrays used by kriging procedure
 **
 ** \return  Error return code
 **
 ** \param[in]  mode       1 for allocation; -1 for deallocation
 ** \param[in]  nvar       Number of variables to be calculated
 ** \param[in]  model      Model structure
 ** \param[in]  neighparam ANeighParam structure
 **
 ** \remarks  The number of variables corresponds to the number of variables
 ** \remarks  to be calculated. It is not necessarily equal to the number of
 ** \remarks  variables contained in Model (when kriging a linear combination
 ** \remarks  of variables for example): hence the use of the 'nvar' passed
 ** \remarks  as an argument
 **
 *****************************************************************************/
static int st_krige_manage(int mode,
                           int nvar,
                           Model *model,
                           ANeighParam *neighparam)
{
  int nech, nfeq, nmax;

  /* Initializations */

  nvar = model->getVariableNumber();
  nfeq = model->getDriftEquationNumber();
  nech = DBIN->getSampleNumber();
  nmax = st_get_nmax(neighparam);

  return (st_krige_manage_basic(mode, nech, nmax, nvar, nfeq));
}

/****************************************************************************/
/*!
 **  Allocate the Target discretization
 **
 ** \param[in]  ndim       Space dimension
 ** \param[in]  ndisc      Discretization parameters (or NULL)
 **
 *****************************************************************************/
static int st_block_discretize_alloc(int ndim, VectorInt ndisc)
{
  int ntot;

  ntot = 1;
  for (int idim = 0; idim < ndim; idim++)
    ntot *= ndisc[idim];
  if (ntot <= 0) return (1);
  KOPTION->ntot = ntot;

  KOPTION->ndisc = st_icore(ndim, 1);
  if (KOPTION->ndisc == nullptr) return (1);
  KOPTION->disc1 = st_core(ndim, ntot);
  if (KOPTION->disc1 == nullptr) return (1);
  KOPTION->disc2 = st_core(ndim, ntot);
  if (KOPTION->disc2 == nullptr) return (1);
  for (int idim = 0; idim < ndim; idim++)
    KOPTION->ndisc[idim] = ndisc[idim];
  return (0);
}

/****************************************************************************/
/*!
 **  Allocate the Data discretization
 **
 ** \param[in] ndim    Space dimension
 **
 *****************************************************************************/
static void st_data_discretize_alloc(int ndim)

{
  int nrow, ncol;

  KOPTION->flag_data_disc = 0;
  if (DBIN->getBlockExtensionNumber() > 0)
  {
    if (!get_keypair("Data_Discretization", &nrow, &ncol, &KOPTION->dsize))
    {
      if (nrow * ncol != ndim)
      {
        messerr("Data discretization is defined using set_keypair mechanism");
        messerr("with keyword 'Data_Discretization'");
        messerr("But its dimension should be %d (instead of %d x %d)", ndim,
                nrow, ncol);
      }
      else
      {
        KOPTION->flag_data_disc = 1;
      }
    }
    else
    {
      if (DBIN->getBlockExtensionNumber() > 0)
      {
        message("\n");
        message("Your Input Data File contains 'dblk' locator(s)\n");
        message("defining a non-ponctual support to the data\n");
        message("This feature can be taken into account during Kriging\n");
        message("ONLY if you specify the discretization steps\n");
        message("for each space dimension, using\n");
        message("       set.keypair('Data_Discretization',c(hx,hy,...))\n");
        message("before calling the kriging() function\n");
        message("\n");
        message("Currently:\n");
        message("- the support is disregarded\n");
        message("- data are considered as ponctual\n");
        message("\n");
      }
    }
  }
}

/****************************************************************************/
/*!
 **  Discretize a block
 **
 ** \param[in]  mode      0 if the block extension is read from grid
 **                       1 if the block extension is read from variable
 ** \param[in]  flag_rand 0 if the second discretization is regular
 **                       1 if the second point must be randomized
 ** \param[in]  iech      rank of the variable (used when mode=1)
 **
 *****************************************************************************/
static void st_block_discretize(int mode, int flag_rand, int iech)
{
  int i, j, jech, ntot, nd, nval, ndim, idim, memo;
  double taille;

  /* Initializations */

  memo = law_get_random_seed();
  ntot = KOPTION->ntot;
  ndim = KOPTION->ndim;
  law_set_random_seed(1234546);
  DbGrid* dbgrid = dynamic_cast<DbGrid*>(DBOUT);

  /* Loop on the discretization points */

  for (i = 0; i < ntot; i++)
  {
    jech = i;
    nval = ntot;
    for (idim = ndim - 1; idim >= 0; idim--)
    {
      taille = (mode == 0) ? dbgrid->getDX(idim) : DBOUT->getBlockExtension(iech, idim);
      nd = KOPTION->ndisc[idim];
      nval /= nd;
      j = jech / nval;
      jech -= j * nval;
      DISC1(i,idim) = taille * ((j + 0.5) / nd - 0.5);
      DISC2(i,idim) = DISC1(i, idim);
      if (flag_rand)
      DISC2(i,idim) += taille * law_uniform(-0.5, 0.5) / (double) nd;
    }
  }
  law_set_random_seed(memo);
}

/****************************************************************************/
/*!
 **  Management of Kriging option
 **
 ** \return  Error return code
 **
 ** \param[in]  mode        1 for allocation; -1 for deallocation
 ** \param[in]  flag_check  1 if the file should be checked
 ** \param[in]  calcul      Type of calculation (EKrigOpt)
 ** \param[in]  flag_rand   0 if the second discretization is regular
 **                         1 if the second point must be randomized
 ** \param[in]  ndisc       Discretization parameters (or NULL)
 **
 ** \remark  This function manages the global structure KOPTION
 **
 *****************************************************************************/
int krige_koption_manage(int mode,
                         int flag_check,
                         const EKrigOpt &calcul,
                         int flag_rand,
                         VectorInt ndisc)
{
  int ndim, error;

  /* Initializations */

  error = 1;
  ndim = DBOUT->getNDim();

  /* Dispatch */

  if (mode == 1)
  {

    /* Allocation of the structure */

    KOPTION = new Koption();
    KOPTION->calcul = calcul;

    // Target discretization
    KOPTION->ndim = ndim;
    KOPTION->ntot = 0;
    KOPTION->disc1 = nullptr;
    KOPTION->disc2 = nullptr;
    KOPTION->ndisc = nullptr;

    // Data discretization
    KOPTION->flag_data_disc = 0;
    KOPTION->dsize = nullptr;

    /* Data discretization case (optional) */

    st_data_discretize_alloc(ndim);

    /* Block discretization case */

    switch (KOPTION->calcul.toEnum())
    {
      case EKrigOpt::E_PONCTUAL:
      case EKrigOpt::E_DRIFT:
        break;

      case EKrigOpt::E_BLOCK:

        /* Preliminary checks */

        if (flag_check && !is_grid(DBOUT))
        {
          messerr("Discretization is not allowed if the Target is not a Grid");
          goto label_dealloc;
        }
        if (ndisc.empty())
        {
          messerr("For block estimation, Discretization must be provided");
          goto label_dealloc;
        }

        if (st_block_discretize_alloc(ndim, ndisc)) goto label_dealloc;

        st_block_discretize(0, flag_rand, 0);

        break;
    }
    error = 0;
  }
  else
  {
    error = 0;

    /* Deallocation procedure */

    label_dealloc: if (KOPTION != nullptr)
    {
      KOPTION->ndisc = (int*) mem_free((char* ) KOPTION->ndisc);
      KOPTION->disc1 = (double*) mem_free((char* ) KOPTION->disc1);
      KOPTION->disc2 = (double*) mem_free((char* ) KOPTION->disc2);
      KOPTION->dsize = (double*) mem_free((char* ) KOPTION->dsize);
      delete KOPTION;
    }
  }

  return (error);
}

/****************************************************************************/
/*!
 **  Define the array flag to convert from isotropic to heterotopic case
 **
 ** \param[in]  model      Model structure
 ** \param[in]  nbgh_ranks Vector of selected samples
 ** \param[in]  neq        Number of equations
 **
 ** \param[out]  nred  Reduced number of equations
 **
 *****************************************************************************/
static void st_flag_define(Model *model,
                           const VectorInt& nbgh_ranks,
                           int neq,
                           int *nred)
{
  int i, iech, ibfl, ib, il, ivar, idim, valid, count, nvar;

  /* Initializations */

  int nech = (int) nbgh_ranks.size();
  *nred = ivar = 0;
  nvar = model->getVariableNumber();
  for (i = 0; i < neq; i++)
    flag[i] = 1;

  /* Check on the coordinates */

  for (iech = 0; iech < nech; iech++)
  {
    valid = 1;
    for (idim = 0; idim < DBIN->getNDim(); idim++)
      if (FFFF(st_get_idim(nbgh_ranks[iech], idim))) valid = 0;
    if (!valid) for (ivar = 0; ivar < DBIN->getVariableNumber(); ivar++)
      FLAG(iech,ivar) = 0;
  }

  /* Check on the data values */

  for (iech = 0; iech < nech; iech++)
    for (ivar = 0; ivar < nvar; ivar++)
      if (FFFF(st_get_ivar(nbgh_ranks[iech], ivar))) FLAG(iech,ivar) = 0;

  /* Check on the external drifts */

  for (iech = 0; iech < nech; iech++)
    for (ibfl = 0; ibfl < model_nfex(model); ibfl++)
      if (FFFF(st_get_fext(nbgh_ranks[iech], ibfl)))
        for (ivar = 0; ivar < DBIN->getVariableNumber(); ivar++)
          FLAG(iech,ivar) = 0;

  /* Check on the drift */

  for (ib = 0; ib < model->getDriftEquationNumber(); ib++)
  {
    valid = 0;
    for (il = 0; il < model->getDriftNumber(); il++)
      for (ivar = 0; ivar < nvar; ivar++)
      {
        if (model->getCoefDrift(ivar, il, ib) == 0.) continue;
        for (iech = 0; iech < nech; iech++)
          if (!FFFF(st_get_ivar(nbgh_ranks[iech], ivar))) valid++;
      }
    FLAG(nech+ib,DBIN->getVariableNumber()-1) = (valid > 0);
  }

  /* Calculate the new number of equations */

  for (i = count = 0; i < neq; i++)
  {
    if (flag[i] != 0) count++;
  }

  /* Returning arguments */

  *nred = count;
  return;
}

/*****************************************************************************/
/*!
 **  Checks if the number of samples is compatible with the number of
 **  drift equations
 **
 ** \return  Error returned code: 1 if an error is found; 0 otherwise
 **
 ** \param[in]  model  Model structure
 ** \param[in]  nbgh_ranks Vector of selected samples
 **
 *****************************************************************************/
static int st_authorize(Model *model, const VectorInt& nbgh_ranks)
{
  int i, nvar, nfeq, n_cov, n_drf, error, nech;

  /* Initializations */

  error = 1;
  nech = (int) nbgh_ranks.size();
  nvar = model->getVariableNumber();
  nfeq = model->getDriftEquationNumber();

  /* Preliminary check */

  if (nech * nvar < nfeq) goto label_end;

  /* Check that enough information is present */

  n_cov = n_drf = 0;
  for (i = 0; i < nvar * nech; i++)
    n_cov += flag[i];
  for (i = 0; i < nfeq; i++)
    n_drf += flag[i + nvar * nech];
  if (n_cov <= 0 || n_cov < n_drf) goto label_end;

  /* Set the error return code */

  error = 0;

  label_end: return (error);
}

/****************************************************************************/
/*!
 **  Establish the constant term for the variance calculation
 **
 ** \param[in]  model    Model structure
 ** \param[in]  nvar     Number of output variables
 ** \param[in]  matCL    Matrix of linear combination (or NULL)
 **                      (Dimension: nvarCL * model->getNVar())
 **
 ** \remarks When 'matCL' is provided, 'nvar' stands for the first dimension of
 ** \remarks the matrix 'matCL' (its second dimension is equal to model->getNVar()).
 ** \remarks Otherwise nvar designates model->getNVar()
 **
 *****************************************************************************/
static void st_variance0(Model *model, int nvar, VectorVectorDouble matCL)
{
  double value;
  int nvar_m;

  /* Initializations */

  nvar_m = model->getVariableNumber();

  /* In the non-stationary case, the calculation is postponed. */
  /* The value is temporarily set to TEST                      */

  if (model->isNoStat())
  {
    for (int ivar = 0; ivar < nvar; ivar++)
      for (int jvar = 0; jvar < nvar; jvar++)
        VAR0(ivar,jvar) = TEST;
  }
  else
  {
    COVINT.setIcas1(2);
    COVINT.setIcas2(2);
    COVINT.setIech1(0);
    COVINT.setIech2(0);
    COVINT.setNdim(model->getDimensionNumber());
    model_variance0(model, KOPTION, covtab, var0);

    // If 'matCL' is provided, some extra-calculations are needed

    if (!matCL.empty())
    {
      for (int ivar = 0; ivar < nvar; ivar++)
        for (int jvar = 0; jvar < nvar; jvar++)
        {
          value = 0.;
          for (int ivar_m = 0; ivar_m < model->getVariableNumber(); ivar_m++)
            for (int jvar_m = 0; jvar_m < model->getVariableNumber(); jvar_m++)
              value += (MATCL(ivar, ivar_m) * COVTAB(ivar_m, jvar_m)
                        * MATCL(jvar, jvar_m));
          VAR0(ivar,jvar) = value;
        }
    }
  }

  return;
}

/****************************************************************************/
/*!
 **  Establish the calculation of variance or standard deviation
 **
 ** \param[in]  model   Model structure
 ** \param[in]  ivar    Rank of the target variable
 ** \param[in]  jvar    Rank of the auxiliary variable
 ** \param[in]  nred    Reduced number of equations
 **
 *****************************************************************************/
static double st_variance(Model *model, int ivar, int jvar, int nred)
{
  double var;
  int i, nvar;

  nvar = model->getVariableNumber();

  // In the stationary case, var0 has already been calculated
  // Otherwise if must be derived here

  if (model->isNoStat() || model->hasExternalCov())
  {
    COVINT.setIcas1(2);
    COVINT.setIcas2(2);
    COVINT.setIech1(IECH_OUT);
    COVINT.setIech2(IECH_OUT);
    COVINT.setNdim(model->getDimensionNumber());
    model_variance0_nostat(model, KOPTION, &COVINT, covtab, var0);
  }

  var = VAR0(ivar, jvar);
  if (FLAG_BAYES) var += VARB(ivar, ivar);
  for (i = 0; i < nred; i++)
    var -= RHS_C(i,jvar) * WGT(i, ivar);

  return (var);
}

/****************************************************************************/
/*!
 **  Establish the variance of the estimator
 **
 ** \param[in]  ivar    Rank of the target variable
 ** \param[in]  jvar    Rank of the auxiliary variable
 ** \param[in]  nfeq    Number of drift equations
 ** \param[in]  nred    Reduced number of equations
 **
 *****************************************************************************/
static double st_varestimate(int ivar, int jvar, int nfeq, int nred)
{
  double var, signe;
  int i, cumflag;

  cumflag = nred - nfeq;

  var = 0.;
  for (i = 0; i < nred; i++)
  {
    signe = (i < cumflag) ? 1. : -1.;
    var += signe * RHS_C(i, jvar) * WGT(i, ivar);
  }
  return (var);
}

/****************************************************************************/
/*!
 **  Returns the additional variance for continuous moving neighborhood
 **
 ** \return  Additional variance or 0
 **
 ** \param[in]  neighparam ANeighParam structure
 ** \param[in]  db1        First Db
 ** \param[in]  rank1      Rank of the sample in the first Db
 ** \param[in]  db2        Second Db
 ** \param[in]  rank2      Rank of the sample in the second Db
 **
 ** \remarks In the case of a neighborhood which is not MOVING (or undefined),
 ** \remarks this function systematically returns 0.
 **
 *****************************************************************************/
static double st_continuous_variance(ANeighParam *neighparam,
                                     Db *db1,
                                     int rank1,
                                     Db *db2,
                                     int rank2)
{
  static double eps_dist = EPSILON4;

  /* Initializations */

  if (neighparam == nullptr) return (0.);
  if (neighparam->getType() != ENeigh::MOVING) return (0.);
  NeighMoving* neighM = dynamic_cast<NeighMoving*>(neighparam);
  int ndim = neighparam->getNDim();
  VectorDouble dd(ndim);

  /* Calculate the distance increment */

  for (int idim = 0; idim < ndim; idim++)
    dd[idim] = db1->getCoordinate(rank1, idim)
        - db2->getCoordinate(rank2, idim);

  /* Anisotropic neighborhood */

  if (neighM->getFlagAniso())
  {

    /* Rotated anisotropy ellipsoid */

    if (neighM->getFlagRotation())
      matrix_product_safe(1, ndim, ndim, dd.data(), neighM->getAnisoRotMats().data(),
                          dd.data());
    for (int idim = 0; idim < ndim; idim++)
      dd[idim] /= neighM->getAnisoCoeff(idim);
  }

  /* Calculate the distance */

  double dist;
  matrix_product(1, ndim, 1, dd.data(), dd.data(), &dist);
  dist = sqrt(dist) / neighM->getRadius();
  double var = 0.;
  if (dist > neighM->getDistCont())
  {
    if (ABS(1. - dist) < eps_dist) dist = 1. - eps_dist;
    var = (dist - neighM->getDistCont()) / (1. - dist);
    var = var * var;
  }

  return (var);
}

/****************************************************************************/
/*!
 **  Establish the kriging L.H.S.
 **
 ** \param[in]  model      Model structure
 ** \param[in]  neighparam ANeighParam structure
 ** \param[in]  nbgh_ranks Vector of selected samples
 ** \param[in]  neq        Number of equations
 **
 *****************************************************************************/
static void st_lhs(Model *model,
                   ANeighParam *neighparam,
                   const VectorInt& nbgh_ranks,
                   int neq)
{
  int i, iech, jech, idim, ivar, jvar, ib, il, nvar_m, nfeq, nbfl, code1, code2, nech;
  double verr, cref, value;

  /* Initializations */

  nech = (int) nbgh_ranks.size();
  nvar_m = model->getVariableNumber();
  nfeq = model->getDriftEquationNumber();
  nbfl = model->getDriftNumber();
  for (i = 0; i < neq * neq; i++)
    lhs[i] = 0.;

  /* Establish the covariance part */

  for (iech = 0; iech < nech; iech++)
    for (jech = 0; jech < nech; jech++)
    {
      model_covtab_init(1, model, covtab);
      for (idim = 0; idim < DBIN->getNDim(); idim++)
        d1[idim] = (st_get_idim(nbgh_ranks[jech], idim)
            - st_get_idim(nbgh_ranks[iech], idim));
      st_cov_dd(model, 0, 0, -1, 1., nbgh_ranks[iech], nbgh_ranks[jech], d1, covtab);

      for (ivar = 0; ivar < nvar_m; ivar++)
        for (jvar = 0; jvar < nvar_m; jvar++)
        {
          LHS(iech,ivar,jech,jvar) = COVTAB(ivar, jvar);

          /* Correction due to measurement errors */

          verr = 0.;
          if (FLAG_PROF)
          {
            code1 = (int) DBIN->getCode(nbgh_ranks[iech]);
            code2 = (int) DBIN->getCode(nbgh_ranks[jech]);
            if (code1 != 0 && code2 != 0 && code1 == code2)
              verr = DBIN->getVarianceError(nbgh_ranks[iech], 0);
          }
          else
          {
            if (iech == jech)
            {
              verr = DBIN->getVarianceError(nbgh_ranks[iech], ivar);

              if (neighparam->getFlagContinuous())
              {
                // In the case of continuous Kriging, we must update the LHS
                // by considering the distance between data and target

                cref = LHS(iech, ivar, jech, jvar);
                verr = cref
                    * st_continuous_variance(neighparam, DBIN, nbgh_ranks[iech], DBOUT,
                                             IECH_OUT);
              }
            }
          }
          if (!FFFF(verr) && verr > 0) LHS(iech,ivar,jech,jvar) += verr;

          // Correction in the DGM case

          if (FLAG_DGM)
          {
            if (iech != jech)
            LHS(iech,ivar,jech,jvar) *= (R_COEFF * R_COEFF);
          }
        }
    }

  /* Establish the drift part */

  if (nfeq <= 0 || nbfl <= 0) return;
  for (iech = 0; iech < nech; iech++)
  {
    if (nbgh_ranks[iech] >= 0)
      model_calcul_drift(model, ECalcMember::LHS, DBIN, nbgh_ranks[iech], drftab);
    else
      model_calcul_drift(model, ECalcMember::LHS, DBOUT, IECH_OUT, drftab);

    for (ivar = 0; ivar < nvar_m; ivar++)
      for (ib = 0; ib < nfeq; ib++)
      {
        value = 0.;
        for (il = 0; il < nbfl; il++)
          value += drftab[il] * model->getCoefDrift(ivar, il, ib);
        LHS(iech,ivar,ib,nvar_m) = value;
        LHS(ib,nvar_m,iech,ivar) = value;
      }
  }
  return;
}

/*****************************************************************************/
/*!
 **  Converts from isotopic to heterotopic the L.H.S
 **
 ** \param[in]  neq  Number of equations
 **
 *****************************************************************************/
static void st_lhs_iso2hetero(int neq)

{
  int lec_lhs = 0;
  int ecr_lhs = 0;
  for (int i = 0; i < neq; i++)
    for (int j = 0; j < neq; j++)
    {
      if (flag[i] != 0 && flag[j] != 0) lhs[ecr_lhs++] = lhs[lec_lhs];
      lec_lhs++;
    }
  return;
}

/****************************************************************************/
/*!
 **  Print the L.H.S. matrix
 **
 ** \param[in]  nech    Number of active points (optional)
 ** \param[in]  neq     Number of equations
 ** \param[in]  nred    Reduced number of equations
 ** \param[in]  flagloc Flag array (optional)
 ** \param[in]  lhs     Kriging L.H.S
 **
 *****************************************************************************/
void krige_lhs_print(int nech,
                     int neq,
                     int nred,
                     int *flagloc,
                     double *lhs)
{
  int *rel, i, j, ipass, npass, ideb, ifin;

  /* Initializations */

  rel = nullptr;
  rel = st_relative_position_array(1, neq, rel);
  npass = (nred - 1) / NBYPAS + 1;

  /* General Header */

  mestitle(0, "LHS of Kriging matrix (compressed)");
  if (nech > 0) message("Number of active samples    = %d\n", nech);
  message("Total number of equations   = %d\n", neq);
  message("Reduced number of equations = %d\n", nred);

  /* Loop on the passes */

  for (ipass = 0; ipass < npass; ipass++)
  {
    ideb = ipass * NBYPAS;
    ifin = MIN(nred, ideb + NBYPAS);
    message("\n");

    /* Header line */

    tab_prints(NULL, "Rank");
    tab_prints(NULL, "    ");
    for (j = ideb; j < ifin; j++)
      tab_printi(NULL, j + 1);
    message("\n");

    /* Flag line */

    if (flagloc != NULL)
    {
      tab_prints(NULL, "    ");
      tab_prints(NULL, "Flag");
      for (j = ideb; j < ifin; j++)
        tab_printi(NULL, rel[j]);
      message("\n");
    }

    /* Matrix lines */

    for (i = 0; i < nred; i++)
    {
      tab_printi(NULL, i + 1);
      tab_printi(NULL, rel[i]);
      for (j = ideb; j < ifin; j++)
        tab_printg(NULL, LHS_C(i, j));
      message("\n");
    }
  }

  rel = st_relative_position_array(-1, neq, rel);
  return;
}

/****************************************************************************/
/*!
 **  Extract the valid data
 **  Operate the product by the inverse covariance matrix
 **
 ** \param[in]  model  Model structure
 ** \param[in]  rmean  Array giving the posterior means for the drift terms
 ** \param[in]  nbgh_ranks Vector of selected samples
 ** \param[in]  nred   Reduced number of equations
 **
 ** \param[out] lterm  Product Z*C-1*Z
 **                    (only produced if FLAG_LTERM is true)
 **
 *****************************************************************************/
static void st_data_dual(Model *model,
                         double *rmean,
                         const VectorInt& nbgh_ranks,
                         int nred,
                         double *lterm)
{
  int i, iech, ivar, ecr, nvar, nfeq, nech;
  double mean;

  /* Initializations */

  nech = (int) nbgh_ranks.size();
  nvar = model->getVariableNumber();
  nfeq = model->getDriftEquationNumber();

  /* Set the whole array to 0 */

  for (i = 0; i < nred; i++)
    zext[i] = 0.;

  /* Extract the data */

  for (ivar = ecr = 0; ivar < nvar; ivar++)
  {
    for (iech = 0; iech < nech; iech++)
    {
      if (!FLAG(iech, ivar)) continue;
      mean = 0.;
      if (nfeq <= 0) mean = model->getContext().getMean(ivar);
      if (FLAG_BAYES)
        mean = model_drift_evaluate(1, model, DBIN, nbgh_ranks[iech],
                                    ivar, rmean, drftab);
      zext[ecr++] = st_get_ivar(nbgh_ranks[iech], ivar) - mean;
    }
  }

  /* Operate the product : Z * A-1 */

  matrix_product(nred, nred, 1, lhs, zext, zam1);

  /* Operate the product : Z * A-1 * Z */

  if (FLAG_LTERM) matrix_product(1, nred, 1, zam1, zext, lterm);

  return;
}

/****************************************************************************/
/*!
 **  Define the array flag[] and the kriging L.H.S.
 **
 ** \param[in]  model      Model structure
 ** \param[in]  neighparam ANeighParam structure
 ** \param[in]  nbgh_ranks Vector of selected indices
 **
 ** \param[out]  status   Kriging error status
 ** \param[out]  nred_r   Reduced number of active points
 ** \param[out]  neq_r    Number of kriging equations
 **
 ** \remark If flag_inv is not switched ON, the inversion is not performed
 ** \remark and the Dual product neither
 ** \remark i.e.: LHS  contains the direct Kriging Matrix
 ** \remark       ZAM1 does not contain any relevant information
 **
 *****************************************************************************/
static void st_prepar(Model *model,
                      ANeighParam *neighparam,
                      const VectorInt& nbgh_ranks,
                      int *status,
                      int *nred_r,
                      int *neq_r)
{
  int nred, neq;

  /* Initializations */

  int nech = (int) nbgh_ranks.size();
  *nred_r = 0;
  *neq_r = 0;
  *status = 1;
  neq = nech * model->getVariableNumber() + model->getDriftEquationNumber();

  /* Define the array flag */

  st_flag_define(model, nbgh_ranks, neq, &nred);

  /* Check if the number of points is compatible with the model */

  if (st_authorize(model, nbgh_ranks)) return;

  /* Establish the Kriging L.H.S. */

  st_lhs(model, neighparam, nbgh_ranks, neq);
  st_lhs_iso2hetero(neq);

  if (OptDbg::query(EDbg::KRIGING)) krige_lhs_print(nech, neq, nred, flag, lhs);

  /* Backup the Kriging matrix before inversion */

  (void) memcpy((char*) lhs_b, (char*) lhs, nred * nred * sizeof(double));

  /* Invert the L.H.S. matrix */

  if (matrix_invert(lhs, nred, IECH_OUT))
  {
    messerr("The Kriging Matrix (%d,%d) is singular", nred, nred);
    messerr("One of the usual reason is the presence of duplicates");
    messerr("To display the Kriging Matrix, run the Kriging procedure again");
    messerr("typing the following command first:");
    messerr("  set.keypair('Dump_Kriging_Matrix',1)");
    if (get_keypone("Dump_Kriging_Matrix", 0))
      print_matrix("Singular matrix", 0, 1, nred, nred, NULL, lhs);
    return;
  }

  /* Returning arguments */

  *status = 0;
  *nred_r = nred;
  *neq_r = neq;

  return;
}

/****************************************************************************/
/*!
 **  Establish the kriging R.H.S
 **
 ** \param[in]  model    Model structure
 ** \param[in]  nbgh_ranks Vector of selected samples
 ** \param[in]  neq      Number of equations
 ** \param[in]  nvar     Number of output variables
 ** \param[in]  matCL    Matrix of linear combinaison (or NULL)
 **                      (Dimension: nvar * model->getNVar())
 **
 ** \param[out]  status  Kriging error status
 **
 ** \remarks When 'matCL' is provided, 'nvar' stands for the first dimension of
 ** \remarks the matrix 'matCL' (its second dimension is equal to model->getNVar()).
 ** \remarks Otherwise nvar designates nvar_m = model->getVariableNumber()
 **
 *****************************************************************************/
static void st_rhs(Model *model,
                   const VectorInt& nbgh_ranks,
                   int neq,
                   int nvar,
                   VectorVectorDouble matCL,
                   int *status)
{
  int    i,iech,ib,nvar_m,nbfl,nfeq,idim,nscale,nech;
  double value,ratio;

  /* Initializations */

  nscale = 1;
  nech = (int) nbgh_ranks.size();
  nvar_m = model->getVariableNumber();
  nbfl = model->getDriftNumber();
  nfeq = model->getDriftEquationNumber();

  /* Establish the covariance part */

  for (iech = 0; iech < nech; iech++)
  {
    switch (KOPTION->calcul.toEnum())
    {
      case EKrigOpt::E_PONCTUAL:
        nscale = 1;
        model_covtab_init(1, model, covtab);
        for (idim = 0; idim < DBIN->getNDim(); idim++)
        {
          d1[idim] = (DBOUT->getCoordinate(IECH_OUT, idim)
              - st_get_idim(nbgh_ranks[iech], idim));
          // The next option is plugged for the case of target randomization
          // for the case of Point-Block Model
          if (RAND_INDEX >= 0 && KOPTION->disc1 != nullptr)
            d1[idim] += DISC1(RAND_INDEX, idim);
        }
        st_cov_dg(model, 0, 0, ECalcMember::RHS, -1, 1., nbgh_ranks[iech], -1, d1,
                  covtab);
        break;

      case EKrigOpt::E_BLOCK:
        nscale = KOPTION->ntot;
        model_covtab_init(1, model, covtab);
        for (i = 0; i < nscale; i++)
        {
          for (idim = 0; idim < DBIN->getNDim(); idim++)
            d1[idim] = (DBOUT->getCoordinate(IECH_OUT, idim)
                - st_get_idim(nbgh_ranks[iech], idim)
                        + DISC1(i, idim));
          st_cov_dg(model, 0, 0, ECalcMember::RHS, -1, 1., nbgh_ranks[iech], -1, d1,
                    covtab);
        }
        break;

      case EKrigOpt::E_DRIFT:
        nscale = 1;
        for (int ivar_m = 0; ivar_m < nvar_m; ivar_m++)
          for (int jvar_m = 0; jvar_m < nvar_m; jvar_m++)
            COVTAB(ivar_m,jvar_m) = 0;
        break;
    }

    /* Normation */

    ratio = 1. / (double) nscale;
    if (FLAG_DGM) ratio *= R_COEFF;
    for (int ivar_m = 0; ivar_m < nvar_m; ivar_m++)
      for (int jvar_m = 0; jvar_m < nvar_m; jvar_m++)
        COVTAB(ivar_m,jvar_m) *= ratio;

    if (matCL.empty())
    {
      for (int jvar = 0; jvar < nvar; jvar++)
        for (int ivar = 0; ivar < nvar; ivar++)
          RHS(iech,ivar,jvar) = COVTAB(ivar, jvar);
    }
    else
    {
      for (int jvar = 0; jvar < nvar; jvar++)
        for (int ivar_m = 0; ivar_m < nvar_m; ivar_m++)
        {
          value = 0.;
          for (int jvar_m = 0; jvar_m < nvar_m; jvar_m++)
            value += MATCL(jvar,jvar_m) * COVTAB(ivar_m, jvar_m);
          RHS(iech,ivar_m,jvar) = value;
        }
    }
  }

  /* Establish the drift part */

  if (nfeq <= 0) return;

  model_calcul_drift(model,ECalcMember::RHS,DBOUT,IECH_OUT,drftab);
  for (int il=0; il<nbfl; il++)
    if (FFFF(drftab[il]))
    {
      *status = 1;
      return;
    }

  if (! matCL.empty())
  {
    if (model->isFlagLinked())
      messageAbort(
          "Kriging of a Linear combination is incompatible with linked drifts");
    for (int ivar = 0; ivar < nvar; ivar++)
      for (int jvar_m = ib = 0; jvar_m < nvar_m; jvar_m++)
        for (int jl = 0; jl < nbfl; jl++, ib++)
        {
          value = 0.;
          for (int il = 0; il < nbfl; il++)
            value += drftab[il] * model->getCoefDrift(jvar_m, il, ib);
          value *= MATCL(ivar, jvar_m);
          RHS(ib,nvar_m,ivar) = value; /* nvar_m is used to calculate shift */
        }
  }
  else
  {
    for (int ivar = 0; ivar < nvar; ivar++)
      for (ib = 0; ib < nfeq; ib++)
      {
        value = 0.;
        for (int il = 0; il < nbfl; il++)
          value += drftab[il] * model->getCoefDrift(ivar, il, ib);
        RHS(ib,nvar_m,ivar) = value; /* nvar_m is used to calculate shift */
      }
  }

  return;
}

/****************************************************************************/
/*!
 **  Establish the drift part of the kriging R.H.S
 **
 ** \param[in]  model     Model structure
 **
 ** \param[out]  status   Kriging error status
 **
 *****************************************************************************/
static void st_ff0(Model *model, int *status)
{
  int ivar, ib, il, nvar, nbfl, nfeq;
  double value;

  /* Initializations */

  nvar = model->getVariableNumber();
  nbfl = model->getDriftNumber();
  nfeq = model->getDriftEquationNumber();

  /* Establish the drift part */

  if (nbfl <= 0 || nfeq <= 0) return;

  model_calcul_drift(model, ECalcMember::RHS, DBOUT, IECH_OUT, drftab);
  for (il = 0; il < nbfl; il++)
    if (FFFF(drftab[il]))
    {
      *status = 1;
      return;
    }

  for (ivar = 0; ivar < nvar; ivar++)
    for (ib = 0; ib < nfeq; ib++)
    {
      value = 0.;
      for (il = 0; il < nbfl; il++)
        value += drftab[il] * model->getCoefDrift(ivar, il, ib);
      FF0(ib,ivar) = value;
    }

  return;
}

/****************************************************************************/
/*!
 **  Converts from isotopic to heterotopic the R.H.S.
 **
 ** \param[in]  neq  Number of equations
 ** \param[in]  nvar Number of variables
 **
 *****************************************************************************/
static void st_rhs_iso2hetero(int neq, int nvar)
{
  int i, ivar, lec_rhs, ecr_rhs;

  /* Loop on the elements */

  lec_rhs = ecr_rhs = 0;
  for (ivar = 0; ivar < nvar; ivar++)
    for (i = 0; i < neq; i++, lec_rhs++)
    {
      if (flag[i] == 0) continue;
      rhs[ecr_rhs++] = rhs[lec_rhs];
    }
  return;
}

/****************************************************************************/
/*!
 **  Print the R.H.S. matrix
 **
 ** \param[in]  nvar     Number of variables
 ** \param[in]  nech     Number of active points (optional)
 ** \param[in]  neq      Number of equations
 ** \param[in]  nred     Reduced number of equations
 ** \param[in]  flag     Flag array (optional)
 ** \param[in]  rhs      Kriging R.H.S. matrix
 **
 *****************************************************************************/
void krige_rhs_print(int nvar,
                     int nech,
                     int neq,
                     int nred,
                     int *flag,
                     double *rhs)
{
  int *rel, i, ivar, idim;

  /* Initializations */

  rel = nullptr;
  rel = st_relative_position_array(1, neq, NULL);

  /* General Header */

  mestitle(0, "RHS of Kriging matrix (compressed)");
  if (nech > 0) message("Number of active samples    = %d\n", nech);
  message("Total number of equations   = %d\n", neq);
  message("Reduced number of equations = %d\n", nred);
  message("Number of right-hand sides  = %d\n", nvar);

  /* Kriging option */

  if (KOPTION != nullptr)
  {
    switch (KOPTION->calcul.toEnum())
    {
      case EKrigOpt::E_PONCTUAL:
        message("Punctual Estimation\n");
        break;

      case EKrigOpt::E_BLOCK:
        message("Block Estimation : Discretization = ");
        for (idim = 0; idim < KOPTION->ndim; idim++)
        {
          if (idim != 0) message(" x ");
          message("%d", KOPTION->ndisc[idim]);
        }
        message("\n");
        break;

      case EKrigOpt::E_DRIFT:
        message("Drift Estimation\n");
        break;
    }
  }
  message("\n");

  /* Header line */

  tab_prints(NULL, "Rank");
  if (flag != nullptr) tab_prints(NULL, "Flag");
  for (ivar = 0; ivar < nvar; ivar++)
    tab_printi(NULL, ivar + 1);
  message("\n");

  /* Matrix lines */

  for (i = 0; i < nred; i++)
  {
    tab_printi(NULL, i + 1);
    if (flag != nullptr) tab_printi(NULL, rel[i]);
    for (ivar = 0; ivar < nvar; ivar++)
      tab_printg(NULL, RHS_C(i, ivar));
    message("\n");
  }

  rel = st_relative_position_array(-1, neq, rel);
  return;
}

/****************************************************************************/
/*!
 **  Print the Dual matrix
 **
 ** \param[in]  nech     Number of active points (optional)
 ** \param[in]  neq      Number of equations
 ** \param[in]  nred     Reduced number of equations
 ** \param[in]  flag     Flag array (optional)
 ** \param[in]  dual     Kriging Dual matrix
 **
 *****************************************************************************/
void krige_dual_print(int nech, int neq, int nred, int *flag, double *dual)
{
  int *rel, i;

  /* Initializations */

  rel = nullptr;
  rel = st_relative_position_array(1, neq, NULL);

  /* General Header */

  mestitle(0, "Dual Vector (completed with zeroes and compressed)");
  if (nech > 0) message("Number of active samples    = %d\n", nech);
  message("Total number of equations   = %d\n", neq);
  message("Reduced number of equations = %d\n", nred);

  /* Header line */

  tab_prints(NULL, "Rank");
  if (flag != nullptr) tab_prints(NULL, "Flag");
  message("\n");

  /* Matrix lines */

  for (i = 0; i < nred; i++)
  {
    tab_printi(NULL, i + 1);
    if (flag != nullptr) tab_printi(NULL, rel[i]);
    tab_printg(NULL, dual[i]);
    message("\n");
  }

  rel = st_relative_position_array(-1, neq, rel);
  return;
}

/****************************************************************************/
/*!
 **  Calculate the final estimation and storage
 **
 ** \param[in]  model        Model structure
 ** \param[in]  rmean        Array giving the posterior means for drift terms
 ** \param[in]  flag_xvalid  when cross-validation option is switched ON
 **                          1: Z*-Z and (Z*-Z)/S*
 **                          2: Z* and S*
 **                          > 0 for ONE Point out
 **                          < 0 for excluding information with same code
 ** \param[in]  status       Kriging error status
 ** \param[in]  nvar         Number of output variables
 ** \param[in]  nred         Reduced number of equations
 **
 *****************************************************************************/
static void st_estimate(Model  *model,
                        double *rmean,
                        int status,
                        int flag_xvalid,
                        int nvar,
                        int nred)
{
  int i, ivar, nfeq;
  double estim, stdv, var, valdat;

  /* Initializations */

  nfeq = model->getDriftEquationNumber();

  /* Estimation */

  if (FLAG_EST != 0)
  {
    for (ivar = 0; ivar < nvar; ivar++)
    {
      estim = 0.;
      if (nfeq <= 0) estim = model->getContext().getMean(ivar);

      if (status == 0 && (nred > 0 || nfeq <= 0 || FLAG_BAYES))
      {
        if (FLAG_BAYES)
          estim = model_drift_evaluate(0, model, DBOUT, IECH_OUT, ivar, rmean, drftab);
        for (i = 0; i < nred; i++)
          estim += RHS_C(i,ivar) * ZAM1(i);

        if (flag_xvalid && FLAG_EST > 0)
        {
          valdat = DBIN->getVariable(IECH_OUT, ivar);
          estim = (FFFF(valdat)) ? TEST : estim - valdat;
        }
      }
      else
      {
        // In case of failure with KS, set the result to mean
        if (nfeq > 0) estim = TEST;
      }
      DBOUT->setArray(IECH_OUT, IPTR_EST + ivar, estim);
    }
  }

  /* Variance of the estimation error */

  if (FLAG_STD != 0)
  {
    for (ivar = 0; ivar < nvar; ivar++)
    {
      if (status == 0 && (nred > 0 || nfeq <= 0 || FLAG_BAYES))
      {
        stdv = st_variance(model, ivar, ivar, nred);
        if (stdv < 0) stdv = 0.;

        stdv = sqrt(stdv);

        if (flag_xvalid && FLAG_STD > 0)
        {
          estim = DBOUT->getArray(IECH_OUT, IPTR_EST + ivar);
          stdv = (FFFF(estim) || stdv <= 0.) ? TEST : estim / stdv;
        }
      }
      else
      {
        stdv = TEST;
      }
      DBOUT->setArray(IECH_OUT, IPTR_STD + ivar, stdv);
    }
  }

  /* Variance of the estimator */

  if (FLAG_VARZ != 0)
  {
    for (ivar = 0; ivar < nvar; ivar++)
    {
      if (status == 0 && (nred > 0 || nfeq <= 0))
        var = st_varestimate(ivar, ivar, nfeq, nred);
      else
        var = TEST;
      DBOUT->setArray(IECH_OUT, IPTR_VARZ + ivar, var);
    }
  }

  return;
}

/****************************************************************************/
/*!
 **  Calculate the final conditional simulation
 **
 ** \param[in]  model     Model structure
 ** \param[out] smean     Array for simulated posterior mean for the drift means
 ** \param[in]  status    Kriging error status
 ** \param[in]  icase     Rank of the PGS and GRF (or -1)
 ** \param[in]  nbsimu    Number of simulations
 ** \param[in]  nbgh_ranks Vector of selected samples
 ** \param[in]  nred      Reduced number of equations
 **
 ** \remark  KS and BAYES are incompatible: we can use mean in both cases
 **
 *****************************************************************************/
static void st_simulate(Model *model,
                        double *smean,
                        int status,
                        int icase,
                        int nbsimu,
                        const VectorInt& nbgh_ranks,
                        int nred)
{
  int isimu, iech, jech, ivar, jvar, lec, ecr, nvar, nfeq, nech;
  double simu, mean, data, value;

  /* Initializations */

  nech = (int) nbgh_ranks.size();
  nvar = model->getVariableNumber();
  nfeq = model->getDriftEquationNumber();

  /* Simulation */

  for (isimu = ecr = 0; isimu < nbsimu; isimu++)
    for (ivar = 0; ivar < nvar; ivar++, ecr++)
    {
      simu = 0.;
      if (nfeq <= 0) simu = model->getContext().getMean(ivar);

      if (status == 0)
      {
        if (FLAG_BAYES)
          simu = model_drift_evaluate(0, model, DBOUT, IECH_OUT, ivar,
                                      &SMEAN(0, isimu), drftab);

        lec = ivar * nred;
        for (jvar = 0; jvar < nvar; jvar++)
          for (iech = 0; iech < nech; iech++)
          {
            if (!FLAG(iech, jvar)) continue;
            jech = nbgh_ranks[iech];

            mean = 0.;
            if (nfeq <= 0) mean = model->getMean(jvar);
            if (FLAG_BAYES)
              mean = model_drift_evaluate(1, model, DBIN, jech, jvar,
                                          &SMEAN(0, isimu), drftab);
            data = st_get_array(jech, isimu, jvar, icase, nbsimu, nvar);
            simu -= wgt[lec++] * (data + mean);
          }

        if (OptDbg::query(EDbg::KRIGING))
        {
          value = DBOUT->getArray(IECH_OUT, IPTR_EST + ecr);
          message("Non-conditional simulation #%d = %lf\n", isimu + 1, value);
          message("Kriged difference = %lf\n", -simu);
          message("Conditional simulation #%d = %lf\n", isimu + 1,
                  value + simu);
        }
      }
      else
      {
        // In case of failure with KS, set the contidioning to mean
        if (nfeq > 0) simu = TEST;
      }

      /* Add the conditioning kriging to the NC simulation at target */
      DBOUT->updSimvar(ELoc::SIMU, IECH_OUT, isimu, ivar, icase, nbsimu, nvar,
                       0, simu);
    }

  return;
}

/****************************************************************************/
/*!
 **  Print the kriging weights
 **
 ** \param[in]  status  Kriging error status
 ** \param[in]  nvar    Number of variables (output)
 ** \param[in]  nvar_m  Number of variables in the Model
 ** \param[in]  nfeq    Number of drift equations
 ** \param[in]  nbgh_ranks Vector of selected samples
 ** \param[in]  nred    Reduced number of equations
 ** \param[in]  icase   Rank of the PGS or GRF
 ** \param[in]  flag    Flag array
 ** \param[in]  wgt     Array of Kriging weights
 **
 ** \remark In the case of simulations (icase>=0), the data vector is not
 ** \remark printed as it changes for every sample, per simulation
 **
 *****************************************************************************/
static void krige_wgt_print(int status,
                            int nvar,
                            int nvar_m,
                            int nfeq,
                            const VectorInt& nbgh_ranks,
                            int nred,
                            int icase,
                            int *flag,
                            double *wgt)
{
  double *sum, value;
  int iwgt, ivar, jvar_m, ivar_m, iech, lec, cumflag, idim, ndim, ib, number,
      flag_value, nech;

  /* Initializations */

  nech = (int) nbgh_ranks.size();
  ndim = DBIN->getNDim();
  sum = (double*) st_core(nvar_m, 1);
  if (sum == nullptr) return;

  /* Header */

  mestitle(0, "(Co-) Kriging weights");

  /* First line */

  tab_prints(NULL, "Rank");
  for (idim = 0; idim < ndim; idim++)
  {
    String strloc = getLocatorName(ELoc::X, idim);
    tab_prints(NULL, strloc.c_str());
  }
  if (DBIN->hasCode()) tab_prints(NULL, "Code");
  if (DBIN->getVarianceErrorNumber() > 0)
    tab_prints(NULL, "Err.");
  if (KOPTION->flag_data_disc) for (idim = 0; idim < ndim; idim++)
  {
    (void) gslSPrintf(string, "Size%d", idim + 1);
    tab_prints(NULL, string);
  }
  tab_prints(NULL, "Data");
  for (ivar = 0; ivar < nvar; ivar++)
  {
    (void) gslSPrintf(string, "Z%d*", ivar + 1);
    tab_prints(NULL, string);
  }
  message("\n");

  /* Display the information and the weights */

  for (jvar_m = lec = cumflag = 0; jvar_m < nvar_m; jvar_m++)
  {
    if (nvar > 1) message("Using variable Z%-2d\n", jvar_m + 1);

    /* Loop on the samples */

    for (ivar_m = 0; ivar_m < nvar_m; ivar_m++)
      sum[ivar_m] = 0.;
    for (iech = 0; iech < nech; iech++, lec++)
    {
      flag_value = (flag != nullptr) ? flag[lec] :
                                       1;
      tab_printi(NULL, iech + 1);
      for (idim = 0; idim < ndim; idim++)
        tab_printg(NULL, st_get_idim(nbgh_ranks[iech], idim));
      if (DBIN->hasCode())
        tab_printg(NULL, DBIN->getCode(nbgh_ranks[iech]));
      if (DBIN->getVarianceErrorNumber() > 0)
        tab_printg(NULL, st_get_verr(nbgh_ranks[iech], (FLAG_PROF) ? 0 : jvar_m));
      if (KOPTION->flag_data_disc)
      {
        for (idim = 0; idim < ndim; idim++)
          tab_printg(NULL, DBIN->getBlockExtension(nbgh_ranks[iech], idim));
      }
      if (icase < 0)
        tab_printg(NULL, st_get_ivar(nbgh_ranks[iech], jvar_m));
      else
        tab_prints(NULL, "   ");

      for (ivar = 0; ivar < nvar; ivar++)
      {
        iwgt = nred * ivar + cumflag;
        value = (wgt != nullptr && status == 0 && flag_value) ? wgt[iwgt] : TEST;
        if (!FFFF(value)) sum[ivar] += value;
        tab_printg(NULL, value);
      }
      if (flag_value) cumflag++;
      message("\n");
    }

    number = 1 + ndim + 1;
    if (DBIN->getVarianceErrorNumber() > 0) number++;
    if (KOPTION->flag_data_disc) number += ndim;
    tab_prints(NULL, "Sum of weights", number, EJustify::LEFT);
    for (ivar = 0; ivar < nvar; ivar++)
    {
      value = (status == 0) ? sum[ivar] : TEST;
      tab_printg(NULL, value);
    }
    message("\n");
  }

  sum = (double*) mem_free((char* ) sum);
  if (nfeq <= 0 || wgt == nullptr) return;

  /* Header */

  mestitle(0, "Drift coefficients");

  /* First line */

  tab_prints(NULL, "Rank");
  tab_prints(NULL, "Lagrange");
  tab_prints(NULL, "Coeff");
  message("\n");

  /* Loop on the drift coefficients */

  cumflag = nred - nfeq;
  for (ib = 0; ib < nfeq; ib++)
  {
    iwgt = ib + cumflag;
    tab_printi(NULL, ib + 1);
    value = (status == 0) ? wgt[iwgt] : TEST;
    tab_printg(NULL, value);
    value = (status == 0) ? zam1[iwgt] : TEST;
    tab_printg(NULL, value);

    message("\n");
  }

  return;
}

/****************************************************************************/
/*!
 **  Store the neighborhood parameters
 **
 ** \param[in]  status    Kriging error status
 ** \param[in]  ntab      Number of neighborhood parameters
 ** \param[in]  tab       Array containing the neighborhood parameters
 **
 *****************************************************************************/
static void st_store_nbgh(int status, int ntab, double *tab)
{
  double value;
  int i;

  /* Loop on the parameters */

  for (i = 0; i < ntab; i++)
  {

    /* Store the parameter */

    value = (status == 0) ? tab[i] : TEST;
    DBOUT->setArray(IECH_OUT, IPTR_NBGH + i, value);
  }

  return;
}

/****************************************************************************/
/*!
 **  Print the results
 **
 ** \param[in]  flag_xvalid  when cross-validation option is switched ON
 **                          1: Z*-Z and (Z*-Z)/S*
 **                          2: Z* and S*
 **                          > 0 for ONE Point out
 **                          < 0 for excluding information with same code
 ** \param[in]  nvar         Number of variables
 ** \param[in]  status       Kriging error status
 **
 *****************************************************************************/
static void st_result_kriging_print(int flag_xvalid, int nvar, int status)
{
  int ivar;
  double value, estim, estval, esterr, sigma, trueval, sterr, stdev;

  /* Header */

  if (flag_xvalid != 0)
    mestitle(0, "Cross-validation results");
  else
    mestitle(0, "(Co-) Kriging results");
  message("Target Sample = %d\n", IECH_OUT + 1);

  /* Loop on the results */

  for (ivar = 0; ivar < nvar; ivar++)
  {
    if (flag_xvalid != 0)
    {
      message("Variable Z%-2d\n", ivar + 1);
      if (FLAG_EST != 0)
      {
        trueval = (status == 0) ? DBIN->getVariable(IECH_OUT, ivar) : TEST;
        estim = (status == 0) ? DBOUT->getArray(IECH_OUT, IPTR_EST + ivar) : TEST;

        if (FLAG_EST > 0)
        {
          estval = (status == 0) ? estim + trueval : TEST;
          esterr = (status == 0) ? estim : TEST;
        }
        else
        {
          estval = (status == 0) ? estim : TEST;
          esterr = (status == 0) ? estim - trueval : TEST;
        }

        tab_printg(" - True value        = ", trueval);
        message("\n");
        tab_printg(" - Estimated value   = ", estval);
        message("\n");
        tab_printg(" - Estimation Error  = ", esterr);
        message("\n");

        if (FLAG_STD != 0)
        {
          stdev = (status == 0) ? DBOUT->getArray(IECH_OUT, IPTR_STD + ivar) : TEST;

          if (FLAG_STD > 0)
          {
            sterr = stdev;
            sigma = (status == 0) ? esterr / stdev : TEST;
          }
          else
          {
            sigma = stdev;
            sterr = (status == 0) ? esterr / stdev : TEST;
          }

          tab_printg(" - Std. deviation    = ", sigma);
          message("\n");
          tab_printg(" - Normalized Error  = ", sterr);
          message("\n");
        }
      }
    }
    else
    {
      message("Variable Z%-2d\n", ivar + 1);
      if (FLAG_EST != 0)
      {
        value = (status == 0) ? DBOUT->getArray(IECH_OUT, IPTR_EST + ivar) : TEST;
        tab_printg(" - Estimate  = ", value);
        message("\n");
      }
      if (FLAG_STD != 0)
      {
        value = (status == 0) ? DBOUT->getArray(IECH_OUT, IPTR_STD + ivar) : TEST;
        tab_printg(" - Std. Dev. = ", value);
        value = (status == 0) ? VAR0(ivar, ivar) : TEST;
        message("\n");
        tab_printg(" - Cov(h=0)  = ", value);
        message("\n");
      }
      if (FLAG_VARZ != 0)
      {
        value = (status == 0) ? DBOUT->getArray(IECH_OUT, IPTR_VARZ + ivar) : TEST;
        tab_printg(" - Var(Z*)   = ", value);
        message("\n");
      }
    }
  }
  return;
}

/****************************************************************************/
/*!
 **  Print the simulation results
 **
 ** \param[in]  nbsimu       Number of simulations
 ** \param[in]  nvar         Number of variables
 ** \param[in]  status       Kriging error status
 **
 *****************************************************************************/
static void st_result_simulate_print(int nbsimu, int nvar, int status)
{
  int ivar, isimu, ecr;
  double value;

  /* Header */

  mestitle(0, "Simulation results");

  /* Loop on the results */

  for (isimu = ecr = 0; isimu < nbsimu; isimu++)
    for (ivar = 0; ivar < nvar; ivar++, ecr++)
    {
      message("Simulation #%d of Z%-2d : ", isimu + 1, ivar + 1);
      value = (status == 0) ? DBOUT->getArray(IECH_OUT, IPTR_EST + ecr) :
                              TEST;
      tab_printg(" = ", value);
      message("\n");
    }
  return;
}

/****************************************************************************/
/*!
 **  Print the neighborhood parameters
 **
 ** \param[in]  status  Kriging error status
 ** \param[in]  tab     Array of neighborhood parameters
 **
 *****************************************************************************/
static void st_res_nbgh_print(int status, double *tab)
{
  if (status != 0) return;

  /* Header */

  mestitle(0, "Neighborhood Parameters");

  message("Number of selected samples          = %d\n", (int) tab[0]);
  message("Maximum neighborhood distance       = %lf\n", tab[1]);
  message("Minimum neighborhood distance       = %lf\n", tab[2]);
  message("Number of non-empty sectors         = %d\n", (int) tab[3]);
  message("Number of consecutive empty sectors = %d\n", (int) tab[4]);

  return;
}

/****************************************************************************/
/*!
 **  Check the consistency of the Colocation specification
 **
 ** \return  Error return code
 **
 ** \param[in]  dbin          Input Db structure
 ** \param[in]  dbout         Output Db structure
 ** \param[in]  rank_colcok   Array of ranks of colocated variables
 **
 ** \remarks The array 'rank_colcok' (if present) must be dimensioned
 ** \remarks to the number of variables in Dbin.
 ** \remarks Each element gives the rank of the colocated variable within Dbout
 ** \remarks or -1 if not colocated
 ** \remarks If the array 'rank_colcok' is absent, colocation option is OFF.
 **
 ** \remarks In input, the numbering in ; rank_colcok' starts from 1
 ** \remarks In output, the numbering starts from 0
 **
 ** \remarks FLAG_COLK is defined in this function
 **
 *****************************************************************************/
static int st_check_colcok(Db *dbin, Db *dbout, int *rank_colcok)
{
  int ivar, jvar;

  FLAG_COLK = 0;
  if (rank_colcok == nullptr) return (0);

  /* Loop on the ranks of the colocated variables */

  for (ivar = 0; ivar < dbin->getVariableNumber(); ivar++)
  {
    jvar = rank_colcok[ivar];
    if (IFFFF(jvar)) jvar = 0;
    if (jvar > dbout->getColumnNumber())
    {
      messerr("Error in the Colocation array:");
      messerr("Input variable (#%d): rank of the colocated variable is %d",
              ivar + 1, jvar);
      messerr("But the Output file only contains %d attributes(s)",
              dbout->getColumnNumber());
      return (1);
    }
    rank_colcok[ivar] = jvar - 1;
  }

  // Assign the array of ranks as a global variable
  RANK_COLCOK = rank_colcok;
  FLAG_COLK = 1;
  return (0);
}

/****************************************************************************/
/*!
 **  Save the (Co-) Kriging weights using the keypair mechanism
 **
 ** \param[in]  status   Kriging error status
 ** \param[in]  iech_out Rank of the output sample
 ** \param[in]  nvar     Number of variables
 ** \param[in]  nbgh_ranks Vector of selected samples
 ** \param[in]  nred     Reduced number of equations
 ** \param[in]  flag     Flag array
 ** \param[in]  wgt      Array of Kriging weights
 **
 *****************************************************************************/
static void st_save_keypair_weights(int status,
                                    int iech_out,
                                    int nvar,
                                    const VectorInt& nbgh_ranks,
                                    int nred,
                                    int *flag,
                                    double *wgt)
{
  double wgtloc, values[5];
  int lec, flag_value, iwgt, cumflag;

  /* Initializations */

  int nech = (int) nbgh_ranks.size();
  if (status != 0) return;
  values[0] = iech_out;

  /* Loop on the output variables */

  for (int jvar = lec = cumflag = 0; jvar < nvar; jvar++)
  {
    values[1] = jvar;

    /* Loop on the input samples */

    for (int iech = 0; iech < nech; iech++, lec++)
    {
      flag_value = (flag != nullptr) ? flag[lec] : 1;
      if (flag_value)
      {
        values[2] = nbgh_ranks[iech];

        /* Loop on the input variables */

        for (int ivar = 0; ivar < nvar; ivar++)
        {
          iwgt = nred * ivar + cumflag;
          wgtloc = (wgt != nullptr && flag_value) ? wgt[iwgt] : TEST;
          if (!FFFF(wgtloc))
          {
            values[3] = ivar;
            values[4] = wgtloc;
            app_keypair("KrigingWeights", 1, 1, 5, values);
          }
        }
        if (flag_value) cumflag++;
      }
    }
  }
  return;
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
 ** \param[in]  neighparam  ANeighParam structure
 ** \param[in]  calcul      Kriging calculation option (EKrigOpt)
 ** \param[in]  ndisc       Array giving the discretization counts
 ** \param[in]  flag_est    Option for storing the estimation
 ** \param[in]  flag_std    Option for storing the standard deviation
 ** \param[in]  flag_varz   Option for storing the variance of the estimator
 ** \param[in]  rank_colcok Option for running Collocated Cokriging
 ** \param[in]  matCL       Matrix of linear combination (or NULL)
 **                         (Dimension: nvarCL * model->getNVar())
 ** \param[in]  namconv     Naming convention
 **
 *****************************************************************************/
int kriging(Db *dbin,
            Db *dbout,
            Model *model,
            ANeighParam *neighparam,
            const EKrigOpt &calcul,
            int flag_est,
            int flag_std,
            int flag_varz,
            VectorInt ndisc,
            VectorInt rank_colcok,
            VectorVectorDouble matCL,
            const NamingConvention& namconv)
{
  // Preliminary checks

  if (neighparam->getType() == ENeigh::IMAGE)
  {
    messerr("This tool cannot function with an IMAGE neighborhood");
    return 1;
  }

  // Initializations

  int iptr_est  = -1;
  int iptr_std  = -1;
  int iptr_varz = -1;
  int nvar = matCL.empty() ? model->getVariableNumber() : matCL.size();

  /* Add the attributes for storing the results */

  if (flag_est != 0)
  {
    iptr_est = dbout->addColumnsByConstant(nvar, 0.);
    if (iptr_est < 0) return 1;
  }
  if (flag_std != 0)
  {
    iptr_std = dbout->addColumnsByConstant(nvar, 0.);
    if (iptr_std < 0) return 1;
  }
  if (flag_varz != 0)
  {
    iptr_varz = dbout->addColumnsByConstant(nvar, 0.);
    if (iptr_varz < 0) return 1;
  }

  /* Setting options */

  KrigingSystem ksys(dbin, dbout, model, neighparam);
  if (ksys.setKrigOptEstim(iptr_est, iptr_std, iptr_varz)) return 1;
  if (ksys.setKrigOptCalcul(calcul, ndisc)) return 1;
  if (ksys.setKrigOptColCok(rank_colcok)) return 1;
  if (ksys.setKrigOptMatCL(matCL)) return 1;
  if (! ksys.isReady()) return 1;

  /* Loop on the targets to be processed */

  for (int iech_out = 0; iech_out < dbout->getSampleNumber(); iech_out++)
  {
    mes_process("Kriging sample", dbout->getSampleNumber(), iech_out);
    if (ksys.estimate(iech_out)) return 1;
  }

  /* Set the error return flag */

  namconv.setNamesAndLocators(dbin, ELoc::Z, nvar, dbout, iptr_varz, "varz", 1,
                              false);
  namconv.setNamesAndLocators(dbin, ELoc::Z, nvar, dbout, iptr_std, "stdev", 1,
                              false);
  namconv.setNamesAndLocators(dbin, ELoc::Z, nvar, dbout, iptr_est, "estim");

  return 0;
}

/****************************************************************************/
/*!
 **  Standard Cross-Validation
 **
 ** \return  Error return code
 **
 ** \param[in]  db          Db structure
 ** \param[in]  model       Model structure
 ** \param[in]  neighparam  ANeighParam structure
 ** \param[in]  flag_kfold  1 if a code (K-FOLD) is used
 ** \param[in]  flag_est    Option for storing the estimation
 **                         1: Z*-Z; -1: Z*
 ** \param[in]  flag_std    Option for storing the standard deviation
 **                         1: (Z*-Z)/S; -1: S
 ** \param[in]  flag_varz   Option for storing the variance of estimator
 ** \param[in]  rank_colcok Option for running Collocated Cokriging
 ** \param[in]  namconv     Naming Convention
 **
 *****************************************************************************/
int xvalid(Db *db,
           Model *model,
           ANeighParam *neighparam,
           int flag_kfold,
           int flag_est,
           int flag_std,
           int flag_varz,
           VectorInt rank_colcok,
           const NamingConvention& namconv)
{
  // Preliminary checks

  if (neighparam->getType() == ENeigh::IMAGE)
  {
    messerr("This tool cannot function with an IMAGE neighborhood");
    return 1;
  }

  // Initializations

  int iptr_est  = -1;
  int iptr_std  = -1;
  int iptr_varz = -1;
  int nvar = model->getVariableNumber();

  /* Add the attributes for storing the results */

  if (flag_est != 0)
  {
    iptr_est = db->addColumnsByConstant(nvar, 0.);
    if (iptr_est < 0) return 1;
  }
  if (flag_std != 0)
  {
    iptr_std = db->addColumnsByConstant(nvar, 0.);
    if (iptr_std < 0) return 1;
  }
  if (flag_varz != 0)
  {
    iptr_varz = db->addColumnsByConstant(nvar, 0.);
    if (iptr_varz < 0) return 1;
  }

  /* Setting options */

  KrigingSystem ksys(db, db, model, neighparam);
  if (ksys.setKrigOptEstim(iptr_est, iptr_std, iptr_varz)) return 1;
  if (ksys.setKrigOptXValid(true, flag_kfold, flag_est > 0, flag_std > 0)) return 1;
  if (ksys.setKrigOptColCok(rank_colcok)) return 1;
  if (! ksys.isReady()) return 1;

  /* Loop on the targets to be processed */

  for (int iech_out = 0; iech_out < db->getSampleNumber(); iech_out++)
  {
    mes_process("Cross_validating sample", db->getSampleNumber(), iech_out);
    if (ksys.estimate(iech_out)) return 1;
  }

  /* Set the error return flag */

  namconv.setNamesAndLocators(db, ELoc::Z, nvar, db, iptr_varz, "varz", 1,
                              false);
  if (flag_est > 0)
    namconv.setNamesAndLocators(db, ELoc::Z, nvar, db, iptr_std, "stderr", 1,
                                false);
  else
    namconv.setNamesAndLocators(db, ELoc::Z, nvar, db, iptr_std, "stdev", 1,
                                false);
  if (flag_std > 0)
    namconv.setNamesAndLocators(db, ELoc::Z, nvar, db, iptr_est, "esterr");
  else
    namconv.setNamesAndLocators(db, ELoc::Z, nvar, db, iptr_est, "estim");

  return 0;
}

/****************************************************************************/
/*!
 **  Kriging in the Gaussian Discrete Model
 **
 ** \return  Error return code
 **
 ** \param[in]  dbin        Input Db structure
 ** \param[in]  dbout       Output Db structure
 ** \param[in]  model       Model structure
 ** \param[in]  neighparam  ANeighParam structure
 ** \param[in]  flag_est    Option for storing the estimation
 ** \param[in]  flag_std    Option for storing the standard deviation
 ** \param[in]  flag_varz   Option for storing the variance of the estimator
 ** \param[in]  rval        Change of support coefficient
 ** \param[in]  namconv     Naming convention
 **
 *****************************************************************************/
int krigdgm(Db *dbin,
            Db *dbout,
            Model *model,
            ANeighParam *neighparam,
            int flag_est,
            int flag_std,
            int flag_varz,
            double rval,
            const NamingConvention& namconv)
 {
  // Preliminary checks

  if (neighparam->getType() == ENeigh::IMAGE)
  {
    messerr("This tool cannot function with an IMAGE neighborhood");
    return 1;
  }

  // Initializations

  int iptr_est  = -1;
  int iptr_std  = -1;
  int iptr_varz = -1;
  int nvar = model->getVariableNumber();

  /* Add the attributes for storing the results */

  if (flag_est != 0)
  {
    iptr_est = dbout->addColumnsByConstant(nvar, 0.);
    if (iptr_est < 0) return 1;
  }
  if (flag_std != 0)
  {
    iptr_std = dbout->addColumnsByConstant(nvar, 0.);
    if (iptr_std < 0) return 1;
  }
  if (flag_varz != 0)
  {
    iptr_varz = dbout->addColumnsByConstant(nvar, 0.);
    if (iptr_varz < 0) return 1;
  }

  /* Setting options */

  KrigingSystem ksys(dbin, dbout, model, neighparam);
  if (ksys.setKrigOptEstim(iptr_est, iptr_std, iptr_varz)) return 1;
  if (ksys.setKrigOptDGM(true, rval)) return 1;
  if (! ksys.isReady()) return 1;

  /* Loop on the targets to be processed */

  for (int iech_out = 0; iech_out < dbout->getSampleNumber(); iech_out++)
  {
    mes_process("Kriging sample", dbout->getSampleNumber(), iech_out);
    if (ksys.estimate(iech_out)) return 1;
  }

  /* Set the error return flag */

  namconv.setNamesAndLocators(dbin, ELoc::Z, nvar, dbout, iptr_varz, "varz", 1,
                              false);
  namconv.setNamesAndLocators(dbin, ELoc::Z, nvar, dbout, iptr_std, "stdev", 1,
                              false);
  namconv.setNamesAndLocators(dbin, ELoc::Z, nvar, dbout, iptr_est, "estim");

  return 0;
}

/****************************************************************************/
/*!
 **  Punctual Kriging based on profiles
 **
 ** \return  Error return code
 **
 ** \param[in]  dbin       input Db structure
 ** \param[in]  dbout      output Db structure
 ** \param[in]  model      Model structure
 ** \param[in]  neighparam ANeighParam structure
 ** \param[in]  flag_est   Option for the storing the estimation
 ** \param[in]  flag_std   Option for the storing the standard deviation
 ** \param[in]  namconv     Naming convention
 **
 *****************************************************************************/
int krigprof(Db *dbin,
             Db *dbout,
             Model *model,
             ANeighParam *neighparam,
             int flag_est,
             int flag_std,
             const NamingConvention& namconv)
{
  // Preliminary checks

  if (neighparam->getType() == ENeigh::IMAGE)
  {
    messerr("This tool cannot function with an IMAGE neighborhood");
    return 1;
  }

  // Initializations

  int iptr_est = -1;
  int iptr_std = -1;
  int nvar = model->getVariableNumber();

  /* Add the attributes for storing the results */

  if (flag_est != 0)
  {
    iptr_est = dbout->addColumnsByConstant(nvar, 0.);
    if (iptr_est < 0) return 1;
  }
  if (flag_std != 0)
  {
    iptr_std = dbout->addColumnsByConstant(nvar, 0.);
    if (iptr_std < 0) return 1;
  }

  /* Setting options */

  KrigingSystem ksys(dbin, dbout, model, neighparam);
  if (ksys.setKrigOptEstim(iptr_est, iptr_std, -1)) return 1;
  if (ksys.setKrigoptCode(true)) return 1;
  if (!ksys.isReady()) return 1;

  /* Loop on the targets to be processed */

  for (int iech_out = 0; iech_out < dbout->getSampleNumber(); iech_out++)
  {
    mes_process("Kriging sample", dbout->getSampleNumber(), iech_out);
    if (ksys.estimate(iech_out)) return 1;
  }

  /* Set the error return flag */

  namconv.setNamesAndLocators(dbin, ELoc::Z, nvar, dbout, iptr_std, "stdev", 1,
                              false);
  namconv.setNamesAndLocators(dbin, ELoc::Z, nvar, dbout, iptr_est, "estim");

  return 0;
}

/****************************************************************************/
/*!
 **  Management of internal arrays used by bayesian procedure
 **
 ** \return  Error return code
 **
 ** \param[in]  mode      1 for allocation; -1 for deallocation
 ** \param[in]  nbsimu    Number of simulation (0 for kriging)
 ** \param[in]  model     Model structure
 **
 ** \param[out] rmean     Array giving the posterior means for the drift terms
 ** \param[out] rcov      Array containing the posterior covariance matrix
 **                       for the drift terms
 ** \param[out] smean     Array for simulated posterior mean for the drift means
 **                       (only if nbsimu > 0)
 **
 *****************************************************************************/
static int bayes_manage(int mode,
                        int nbsimu,
                        Model *model,
                        double **rmean,
                        double **rcov,
                        double **smean)
{
  int nfeq;

  /* Initializations */

  nfeq = model->getDriftEquationNumber();

  /* Dispatch */

  if (mode == 1)
  {

    /* Allocation */

    *rmean = st_core(nfeq, 1);
    if (*rmean == nullptr) return (1);
    *rcov = st_core(nfeq, nfeq);
    if (*rcov == nullptr) return (1);
    if (nbsimu > 0)
    {
      *smean = st_core(nfeq, nbsimu);
      if (*smean == nullptr) return (1);
    }
  }
  else
  {

    /* Deallocation */

    *rmean = (double*) mem_free((char* ) (*rmean));
    *rcov = (double*) mem_free((char* ) (*rcov));
    if (nbsimu > 0) *smean = (double*) mem_free((char* ) (*smean));
  }
  return (0);
}

/****************************************************************************/
/*!
 **  Perform the Bayesian estimation of the Drift Coefficients
 **
 ** \return  Error return code
 **
 ** \param[in]  model      Model structure
 ** \param[in]  neighparam ANeighParam structure
 ** \param[in]  nbghw      NeighWork structure
 ** \param[in]  dmean      Array giving the prior means for the drift terms
 ** \param[in]  dcov       Array containing the prior covariance matrix
 **                        for the drift terms
 **
 ** \param[out] rmean      Array giving the posterior means for the drift terms
 ** \param[out] rcov       Array containing the  posterior covariance matrix
 **                        for the drift terms
 **
 *****************************************************************************/
static int bayes_precalc(Model *model,
                         ANeighParam *neighparam,
                         NeighWork& nbghw,
                         double *dmean,
                         double *dcov,
                         double *rmean,
                         double *rcov)
{
  int nfeq,error,status,nech,nred,neq,shift,ib,jb,il,jl,flag_fix;
  double *ff,*smu,*sigma,*vars;
  VectorInt nbgh_ranks;

  /* Initializations */

  error = 1;
  nfeq = model->getDriftEquationNumber();
  ff = smu = sigma = vars = nullptr;

  /* Preliminary Checks */

  if (neighparam->getType() != ENeigh::UNIQUE)
  {
    messerr("The Bayesian Estimation of the Drift Coefficients");
    messerr("is only available in Unique Neighborhood");
    goto label_end;
  }

  /* Core allocation */

  IECH_OUT = DBIN->getSampleNumber() / 2;

  /* Check that the variance-covariance matrix is symmetric */

  flag_fix = is_matrix_null(nfeq, nfeq, dcov, 0);
  if (!is_matrix_symmetric(nfeq, dcov, 1)) goto label_end;

  /* Prepare the Kriging matrix (without correction) */

  nbgh_ranks = nbghw.select(DBOUT, IECH_OUT);
  nech = (int) nbgh_ranks.size();
  status = nbgh_ranks.empty();

  FLAG_BAYES = 0;
  st_prepar(model, neighparam, nbgh_ranks, &status, &nred, &neq);
  FLAG_BAYES = 1;
  if (status) goto label_end;
  shift = nred - nfeq;

  /* Complementary core allocation */

  ff = st_core(shift, nfeq);
  if (ff == nullptr) goto label_end;
  smu = st_core(shift, 1);
  if (smu == nullptr) goto label_end;
  sigma = st_core(shift, shift);
  if (sigma == nullptr) goto label_end;
  vars = st_core(shift, 1);
  if (vars == nullptr) goto label_end;

  // Create the array of variables

  ib = 0;
  for (int iech = 0; iech < DBIN->getSampleNumber(); iech++)
  {
    if (!DBIN->isActive(iech)) continue;
    for (int ivar = 0; ivar < DBIN->getVariableNumber(); ivar++)
    {
      double value = DBIN->getVariable(nbgh_ranks[iech], ivar);
      if (FFFF(value)) continue;
      vars[ib++] = value;
    }
  }

  /* Copy DCOV into S and DMEAN into RMEAN */

  (void) memcpy((char*) rcov, (char*) dcov, sizeof(double) * nfeq * nfeq);
  (void) memcpy((char*) rmean, (char*) dmean, sizeof(double) * nfeq);
  if (flag_fix) goto label_print;

  /* Establish the drift array FF */

  for (il = 0; il < nfeq; il++)
    for (ib = 0; ib < shift; ib++)
      FF(ib,il) = LHS_B(ib, shift + il);

  /* Calculate S-1 */

  if (matrix_invert(rcov, nfeq, -1)) goto label_end;

  /* Calculate: SMU = S-1 * MEAN */

  matrix_product(nfeq, nfeq, 1, rcov, dmean, smu);

  /* Covariance matrix SIGMA */

  for (ib = 0; ib < shift; ib++)
    for (jb = 0; jb < shift; jb++)
      SIGMA(ib,jb) = LHS_B(ib, jb);

  /* Calculate SIGMA-1 */

  if (matrix_invert(sigma, shift, -1)) goto label_end;

  /* Inverse of posterior covariance matrix: SC-1 = FFt * SIGMA-1 * FF + S-1 */

  for (il = 0; il < nfeq; il++)
    for (jl = 0; jl < nfeq; jl++)
    {
      double value = 0.;
      for (ib=0; ib<shift; ib++)
        for (jb=0; jb<shift; jb++)
          value += FF(ib,il) * SIGMA(ib,jb) * FF(jb,jl);
      RCOV(il,jl) += value;
    }

  /* Calculating: SMU = FFt * SIGMA-1 * Z + S-1 * MU */

  for (il = 0; il < nfeq; il++)
  {
    double value = 0.;
    for (ib=0; ib<shift; ib++)
      for (jb=0; jb<shift; jb++)
        value += FF(ib,il) * SIGMA(ib,jb) * vars[jb];
    SMU(il) += value;
  }

  /* Posterior mean: RMEAN = SC * SMU */

  if (matrix_invert(rcov, nfeq, -1)) goto label_end;
  matrix_product(nfeq, nfeq, 1, rcov, smu, rmean);

  label_print: if (OptDbg::query(EDbg::BAYES))
  {
    mestitle(0, "Bayesian Drift coefficients");
    print_matrix("Prior Mean", 0, 1, nfeq, 1, NULL, dmean);
    print_matrix("Prior Variance-Covariance", 0, 1, nfeq, nfeq, NULL, dcov);

    print_matrix("Posterior Mean", 0, 1, nfeq, 1, NULL, rmean);
    print_matrix("Posterior Variance-Covariance", 0, 1, nfeq, nfeq, NULL, rcov);
    message("\n");
  }

  /* Set the error return code */

  error = 0;

  label_end:

  /* Core deallocation */

  ff = (double*) mem_free((char* ) ff);
  smu = (double*) mem_free((char* ) smu);
  sigma = (double*) mem_free((char* ) sigma);
  vars = (double*) mem_free((char* ) vars);

  return (error);
}

/****************************************************************************/
/*!
 **  Correct the arrays RHS and VARB in Bayesian case
 **
 ** \param[in]  model  Model structure
 ** \param[in]  rcov   Array containing the posterior covariance matrix
 **                    for the drift terms
 **
 ** \param[out] status Returned status
 **
 *****************************************************************************/
static void st_bayes_correct(Model *model, double *rcov, int *status)
{
  int ivar, jvar, il, jl, nvar, nfeq;

  /* Initializations */

  nvar = model->getVariableNumber();
  nfeq = model->getDriftEquationNumber();

  /* Establish the Drift matrix */

  st_ff0(model, status);
  if (*status) return;

  /* Correct the arrays */

  for (ivar = 0; ivar < nvar; ivar++)
    for (jvar = 0; jvar < nvar; jvar++)
    {
      VARB(ivar,jvar) = 0.;
      for (il = 0; il < nfeq; il++)
        for (jl = 0; jl < nfeq; jl++)
          VARB(ivar,jvar) += FF0(il,ivar) * RCOV(il, jl) * FF0(jl, jvar);
    }
  return;
}

/****************************************************************************/
/*!
 **  Estimation with Bayesian Drift
 **
 ** \return  Error return code
 **
 ** \param[in]  dbin       input Db structure
 ** \param[in]  dbout      output Db structure
 ** \param[in]  model      Model structure
 ** \param[in]  neighparam ANeighParam structure
 ** \param[in]  prior_mean Array giving the prior means for the drift terms
 ** \param[in]  prior_cov  Array containing the prior covariance matrix
 **                        for the drift terms
 ** \param[in]  flag_est   Pointer for the storing the estimation
 ** \param[in]  flag_std   Pointer for the storing the standard deviation
 ** \param[in]  namconv     Naming convention
 **
 *****************************************************************************/
int kribayes(Db *dbin,
             Db *dbout,
             Model *model,
             ANeighParam *neighparam,
             const VectorDouble& prior_mean,
             const VectorDouble& prior_cov,
             int flag_est,
             int flag_std,
             const NamingConvention& namconv)
{
  // Preliminary checks

  if (neighparam->getType() == ENeigh::IMAGE)
  {
    messerr("This tool cannot function with an IMAGE neighborhood");
    return 1;
  }

  // Initializations

  int iptr_est  = -1;
  int iptr_std  = -1;
  int nvar = model->getVariableNumber();

  /* Add the attributes for storing the results */

  if (flag_est != 0)
  {
    iptr_est = dbout->addColumnsByConstant(nvar, 0.);
    if (iptr_est < 0) return 1;
  }
  if (flag_std != 0)
  {
    iptr_std = dbout->addColumnsByConstant(nvar, 0.);
    if (iptr_std < 0) return 1;
  }

  /* Setting options */

  KrigingSystem ksys(dbin, dbout, model, neighparam);
  if (ksys.setKrigOptEstim(iptr_est, iptr_std, -1)) return 1;
  if (ksys.setKrigOptBayes(true, prior_mean, prior_cov)) return 1;
  if (! ksys.isReady()) return 1;

  /* Loop on the targets to be processed */

  for (int iech_out = 0; iech_out < dbout->getSampleNumber(); iech_out++)
  {
    mes_process("Kriging sample", dbout->getSampleNumber(), iech_out);
    if (ksys.estimate(iech_out)) return 1;
  }

  /* Set the error return flag */

  namconv.setNamesAndLocators(dbin, ELoc::Z, nvar, dbout, iptr_std, "stdev", 1,
                              false);
  namconv.setNamesAndLocators(dbin, ELoc::Z, nvar, dbout, iptr_est, "estim");

  return 0;
}

/****************************************************************************/
/*!
 **  Check the Neighborhood
 **
 ** \return  Error return code
 **
 ** \param[in]  dbin       input Db structure
 ** \param[in]  dbout      output Db structure
 ** \param[in]  model      Model structure (optional)
 ** \param[in]  neighparam ANeighParam structure
 ** \param[in]  namconv    Naming Convention
 **
 ** \remark This procedure creates the following arrays:
 ** \remark 1 - The number of selected samples
 ** \remark 2 - The maximum neighborhood distance
 ** \remark 3 - The minimum neighborhood distance
 ** \remark 4 - The number of non-empty sectors
 ** \remark 5 - The number of consecutive empty sectors
 **
 *****************************************************************************/
int test_neigh(Db *dbin,
               Db *dbout,
               Model *model,
               ANeighParam *neighparam,
               const NamingConvention &namconv)
{
  int error, status, nech, ntab, iext;
  NeighWork nbghw;
  VectorInt nbgh_ranks;

  /* Preliminary checks */

  error = 1;
  iext = -1;
  ntab = 5;
  VectorDouble tab(ntab,0.);
  st_global_init(dbin, dbout);
  if (st_check_environment(1, 1, model, neighparam)) goto label_end;
  if (manage_external_info(1, ELoc::F, DBIN, DBOUT, &iext)) goto label_end;
  if (manage_nostat_info(1, model, DBIN, DBOUT)) goto label_end;

  /* Add the attributes for storing the results */

  IPTR_NBGH = dbout->addColumnsByConstant(ntab, 0.);
  if (IPTR_NBGH < 0) goto label_end;

  /* Pre-calculations */

  nbghw.initialize(DBIN, neighparam);
  nbghw.setFlagSimu(FLAG_SIMU);
  if (st_model_manage(1, model)) goto label_end;
  if (st_krige_manage(1, model->getVariableNumber(), model, neighparam))
    goto label_end;

  /* Loop on the targets to be processed */

  status = 0;
  for (IECH_OUT = 0; IECH_OUT < DBOUT->getSampleNumber(); IECH_OUT++)
  {
    mes_process("Neighborhood Test", DBOUT->getSampleNumber(), IECH_OUT);
    OptDbg::setIndex(IECH_OUT + 1);
    if (!dbout->isActive(IECH_OUT)) continue;
    if (OptDbg::query(EDbg::KRIGING) || OptDbg::query(EDbg::NBGH) || OptDbg::query(EDbg::RESULTS))
    {
      mestitle(1, "Target location");
      db_sample_print(dbout, IECH_OUT, 1, 0, 0);
    }

    /* Select the Neighborhood */

    nbgh_ranks = nbghw.select(DBOUT, IECH_OUT);
    nech = (int) nbgh_ranks.size();
    status = nbgh_ranks.empty();

    /* Retrieve the neighborhood parameters */

    tab = nbghw.summary(DBOUT, IECH_OUT);

    /* Store the neighborhood parameters */

    st_store_nbgh(status, ntab, tab.data());
    if (OptDbg::query(EDbg::NBGH)) st_res_nbgh_print(status, tab.data());
  }

  /* Set the error return flag */

  error = 0;
  namconv.setNamesAndLocators(NULL, ELoc::UNKNOWN, 1, dbout, IPTR_NBGH,
                              "Number");
  namconv.setNamesAndLocators(NULL, ELoc::UNKNOWN, 1, dbout, IPTR_NBGH + 1,
                              "MaxDist");
  namconv.setNamesAndLocators(NULL, ELoc::UNKNOWN, 1, dbout, IPTR_NBGH + 2,
                              "MinDist");
  namconv.setNamesAndLocators(NULL, ELoc::UNKNOWN, 1, dbout, IPTR_NBGH + 3,
                              "NbNESect");
  namconv.setNamesAndLocators(NULL, ELoc::UNKNOWN, 1, dbout, IPTR_NBGH + 4,
                              "NbCESect");
  namconv.setLocators(dbout, IPTR_NBGH, ntab);

  label_end: OptDbg::setIndex(0);
  (void) st_model_manage(-1, model);
  (void) st_krige_manage(-1, model->getVariableNumber(), model, neighparam);
  (void) manage_external_info(-1, ELoc::F, DBIN, DBOUT, &iext);
  (void) manage_nostat_info(-1, model, DBIN, DBOUT);
  return (error);
}

/****************************************************************************/
/*!
 **  Conditioning Kriging
 **
 ** \return  Error return code
 **
 ** \param[in]  strloc     Message used for process
 ** \param[in]  dbin       input Db structure
 ** \param[in]  dbout      output Db structure
 ** \param[in]  model      Model structure
 ** \param[in]  neighparam ANeighParam structure
 ** \param[in]  dmean      Array giving the prior means for the drift terms
 ** \param[in]  dcov       Array containing the prior covariance matrix
 **                        for the drift terms
 ** \param[in]  icase      Case for PGS and GRF (or -1)
 ** \param[in]  nbsimu     Number of simulations
 ** \param[in]  flag_dgm   1 if the DGM version of kriging should be used
 ** \param[in]  rval       Change of support coefficient
 **
 *****************************************************************************/
static int st_krigsim_old(const char *strloc,
                          Db *dbin,
                          Db *dbout,
                          Model *model,
                          ANeighParam *neighparam,
                          double *dmean,
                          double *dcov,
                          int icase,
                          int nbsimu,
                          int flag_dgm,
                          double rval)
{
  int error, status, nech, neq, nred, nvar, iext, nfeq;
  double *rmean, *rcov, *smean;
  Model *model_sk;
  NeighWork nbghw;
  VectorInt nbgh_ranks;

  /* Preliminary checks */

  error = 1;
  iext = -1;
  nvar = neq = nred = 0;
  model_sk = nullptr;
  rmean = smean = rcov = nullptr;
  st_global_init(dbin, dbout);
  FLAG_BAYES = (dmean != nullptr);
  FLAG_EST = 1;
  FLAG_STD = 0;
  FLAG_WGT = 1;
  FLAG_SIMU = 1;
  FLAG_DGM = flag_dgm;
  R_COEFF = rval;
  IPTR_EST = dbout->getColIdxByLocator(ELoc::SIMU, 0);
  if (st_check_environment(1, 1, model, neighparam)) goto label_end;
  if (manage_external_info(1, ELoc::F, DBIN, DBOUT, &iext)) goto label_end;
  if (manage_nostat_info(1, model, DBIN, DBOUT)) goto label_end;
  nvar = model->getVariableNumber();
  nfeq = model->getDriftEquationNumber();

  /* Core allocation */

  if (FLAG_BAYES)
  {
    if (bayes_manage(1, nbsimu, model, &rmean, &rcov, &smean)) goto label_end;
  }

  /* Pre-calculations */

  nbghw.initialize(DBIN, neighparam);
  nbghw.setFlagSimu(FLAG_SIMU);
  if (st_model_manage(1, model)) goto label_end;
  if (st_krige_manage(1, nvar, model, neighparam)) goto label_end;
  if (krige_koption_manage(1, 1, EKrigOpt::PONCTUAL, 1, VectorInt()))
    goto label_end;
  if (FLAG_STD != 0) st_variance0(model, nvar, VectorVectorDouble());

  /* Solve the Bayesian estimation of the Drift coefficients */

  if (FLAG_BAYES)
  {
    if (bayes_precalc(model, neighparam, nbghw, dmean, dcov, rmean, rcov))
      goto label_end;
  }

  /* Simulate the drift coefficients from the posterior distributions */

  if (FLAG_BAYES)
  {
    if (bayes_simulate(model, nbsimu, rmean, rcov, smean)) goto label_end;
  }

  /* Duplicate the model, suppressing the Drift terms */

  if (FLAG_BAYES)
  {
    model_sk = model_duplicate(model, 0., -1);
    if (model_sk == nullptr) goto label_end;
  }
  else
  {
    model_sk = model;
  }

  /* Loop on the targets to be processed */

  status = 0;
  for (IECH_OUT = 0; IECH_OUT < DBOUT->getSampleNumber(); IECH_OUT++)
  {
    mes_process(strloc, DBOUT->getSampleNumber(), IECH_OUT);
    OptDbg::setIndex(IECH_OUT + 1);
    if (!dbout->isActive(IECH_OUT)) continue;
    if (OptDbg::query(EDbg::KRIGING) || OptDbg::query(EDbg::NBGH) || OptDbg::query(EDbg::RESULTS))
    {
      mestitle(1, "Target location");
      db_sample_print(dbout, IECH_OUT, 1, 0, 0);
    }

    /* Select the Neighborhood */

    nbgh_ranks = nbghw.select(DBOUT, IECH_OUT);
    nech = (int) nbgh_ranks.size();
    status = nbgh_ranks.empty();
    if (status) goto label_store;

    /* Establish the kriging L.H.S. */

    if (! nbghw.isUnchanged() || neighparam->getFlagContinuous() || OptDbg::force())
    {
      st_prepar(model_sk, neighparam, nbgh_ranks, &status, &nred, &neq);
      if (status) goto label_store;
    }

    /* Establish the kriging R.H.S. */

    st_rhs(model_sk, nbgh_ranks, neq, nvar, VectorVectorDouble(), &status);
    if (status) goto label_store;
    st_rhs_iso2hetero(neq, nvar);

    /* Modify the arrays in the Bayesian case */

    if (FLAG_BAYES)
    {
      st_bayes_correct(model, dcov, &status);
      if (status) goto label_store;
    }
    if (OptDbg::query(EDbg::KRIGING))
      krige_rhs_print(nvar, nech, neq, nred, flag, rhs);

    /* Derive the kriging weights */

    if (FLAG_WGT)
    {
      matrix_product(nred, nred, nvar, lhs, rhs, wgt);
      if (OptDbg::query(EDbg::KRIGING))
        krige_wgt_print(status, nvar, nvar, nfeq, nbgh_ranks, nred, icase, flag, wgt);
    }

    /* Perform the simulation */

    label_store:
    st_simulate(model, smean, status, icase, nbsimu, nbgh_ranks, nred);
    if (OptDbg::query(EDbg::RESULTS))
      st_result_simulate_print(nbsimu, nvar, status);
  }

  /* Set the error return flag */

  error = 0;

  label_end: OptDbg::setIndex(0);
  if (FLAG_BAYES)
  {
    model_sk = model_free(model_sk);
    (void) bayes_manage(-1, nbsimu, model, &rmean, &rcov, &smean);
  }
  (void) st_model_manage(-1, model);
  (void) st_krige_manage(-1, nvar, model, neighparam);
  (void) krige_koption_manage(-1, 1, EKrigOpt::PONCTUAL, 1, VectorInt());
  (void) manage_external_info(-1, ELoc::F, DBIN, DBOUT, &iext);
  (void) manage_nostat_info(-1, model, DBIN, DBOUT);
  return (error);
}

static int st_krigsim_new(const char *strloc,
                          Db *dbin,
                          Db *dbout,
                          Model *model,
                          ANeighParam *neighparam,
                          double *dmean,
                          double *dcov,
                          int icase,
                          int nbsimu,
                          int flag_dgm,
                          double rval)
{
  // Preliminary checks

  if (neighparam->getType() == ENeigh::IMAGE)
  {
    messerr("This tool cannot function with an IMAGE neighborhood");
    return 1;
  }

  // Initializations

  int iptr_est  = -1;
  int nvar = model->getVariableNumber();
  int nfeq = model->getDriftNumber();

  /* Add the attributes for storing the results */

  iptr_est = dbout->getColIdxByLocator(ELoc::SIMU, 0);
  if (iptr_est < 0) return 1;

  /* Bayesian temporary code */

  bool flag_bayes = false;
  VectorDouble prior_mean;
  VectorDouble prior_cov;
  if (dmean != nullptr && dcov != nullptr && nfeq > 0)
  {
    flag_bayes = true;
    prior_mean = ut_vector_set(dmean, nfeq);
    prior_cov  = ut_vector_set(dcov, nfeq * nfeq);
  }

  /* Setting options */

  KrigingSystem ksys(dbin, dbout, model, neighparam);
  if (ksys.setKrigOptFlagSimu(true, nbsimu, icase)) return 1;
  if (ksys.setKrigOptEstim(iptr_est, -1, -1)) return 1;
  if (ksys.setKrigOptBayes(flag_bayes, prior_mean, prior_cov)) return 1;
  if (ksys.setKrigOptDGM(flag_dgm, rval)) return 1;
  if (! ksys.isReady()) return 1;

  /* Loop on the targets to be processed */

  for (int iech_out = 0; iech_out < dbout->getSampleNumber(); iech_out++)
  {
    mes_process(strloc, dbout->getSampleNumber(), iech_out);
    if (ksys.estimate(iech_out)) return 1;
  }

  return 0;
}

int _krigsim(const char *strloc,
             Db *dbin,
             Db *dbout,
             Model *model,
             ANeighParam *neighparam,
             double *dmean,
             double *dcov,
             int icase,
             int nbsimu,
             int flag_dgm,
             double rval)
{
  double value;

  int test = (int) get_keypone("krigsim", 0);
  if (test == 0)
    return st_krigsim_old(strloc, dbin, dbout, model, neighparam, dmean, dcov, icase, nbsimu, flag_dgm, rval);
  else
    return st_krigsim_new(strloc, dbin, dbout, model, neighparam, dmean, dcov, icase, nbsimu, flag_dgm, rval);
}

/****************************************************************************/
/*!
 **  Global estimation over a territory using arithmetic average
 **
 ** \return A Global_Res structure
 **
 ** \param[in]  dbin          Db structure
 ** \param[in]  dbgrid        Db grid structure
 ** \param[in]  model         Model structure
 ** \param[in]  ivar0         Rank of the target variable
 ** \param[in]  flag_verbose  1 for a verbose output
 **
 *****************************************************************************/
Global_Res global_arithmetic(Db *dbin,
                             DbGrid *dbgrid,
                             Model *model,
                             int ivar0,
                             bool flag_verbose)
{
  Global_Res gres;

  /* Preliminary checks */

  if (ivar0 < 0 || ivar0 >= dbin->getVariableNumber())
  {
    messerr("The target variable (%d) must lie between 1 and the number of variables (%d)",
            ivar0 + 1, dbin->getVariableNumber());
    return gres;
  }

  /* Preliminary assignments */

  int ntot = dbin->getSampleNumber(false);
  int np = dbin->getSampleNumber(true);
  int ng = dbgrid->getSampleNumber(true);
  double surface = ng * db_grid_maille(dbgrid);

  /* Average covariance over the data */

  double cxx = model_cxx(model, dbin, dbin, ivar0, ivar0, 0, 0.);

  /* Average covariance between the data and the territory */

  double cxv = model_cxx(model, dbin, dbgrid, ivar0, ivar0, 0, 0.);

  /* Average covariance over the territory */

  double cvv = model_cxx(model, dbgrid, dbgrid, ivar0, ivar0, 0,
                         db_epsilon_distance(dbgrid));

  /* Calculating basic statistics */

  int iatt = db_attribute_identify(dbin, ELoc::Z, ivar0);
  double wtot;
  double ave;
  double var;
  double mini;
  double maxi;
  db_monostat(dbin, iatt, &wtot, &ave, &var, &mini, &maxi);

  /* Filling the resulting structure */

  double sse = cvv - 2. * cxv + cxx;
  sse = (sse > 0) ? sqrt(sse) : 0.;
  double cvsam = (ave != 0.) ? sqrt(var) / ave : TEST;
  double cviid = cvsam / sqrt(np);
  double cvgeo = (ave != 0.) ? sse / ave : TEST;

  /* Filling the output structure */

  gres.ntot = ntot;
  gres.np = np;
  gres.ng = ng;
  gres.surface = surface;
  gres.zest = ave;
  gres.sse  = sse;
  gres.cvgeo = cvgeo;
  gres.cvv = cvv;
  gres.weights.resize(np, 1./np);

  if (flag_verbose)
  {
    mestitle(1,"Global estimation by arithmetic average");
    message("Total number of data             = %d\n", ntot);
    message("Number of active data            = %d\n", np);
    message("Sample variance                  = %lf\n", var);
    message("CVsample                         = %lf\n", cvsam);
    message("CViid                            = %lf\n", cviid);
    message("Cxx                              = %lf\n", cxx);
    message("Cxv                              = %lf\n", cxv);
    message("Cvv                              = %lf\n", cvv);
    if (FFFF(ave))
      message("Estimation by arithmetic average = NA\n");
    else
      message("Estimation by arithmetic average = %lf\n", ave);
    message("Estimation St. dev. of the mean  = %lf\n", sse);
    if (FFFF(cvgeo))
      message("CVgeo                            = NA\n");
    else
      message("CVgeo                            = %lf\n", cvgeo);
    message("Surface                          = %lf\n", surface);
    if (FFFF(ave))
      message("Q (Estimation * Surface)         = NA\n");
    else
      message("Q (Estimation * Surface)         = %lf\n", ave * surface);
    message("\n");
  }

  return gres;
}

/****************************************************************************/
/*!
 **  Global estimation over a territory using kriging
 **
 ** \return A Global_Res structure
 **
 ** \param[in]  dbin          Db structure
 ** \param[in]  dbout         Db grid structure
 ** \param[in]  model         Model structure
 ** \param[in]  ivar0         Rank of the target variable
 ** \param[in]  flag_verbose  1 for a verbose output
 **
 *****************************************************************************/
Global_Res global_kriging(Db *dbin,
                          Db *dbout,
                          Model *model,
                          int ivar0,
                          bool flag_verbose)
{
  NeighUnique neighU;
  Global_Res gres;
  VectorDouble rhsCum;

  /* Preliminary tests */

  if (ivar0 < 0 || ivar0 >= dbin->getVariableNumber())
  {
    messerr("The target variable (%d) must lie between 1 and the number of variables (%d)",
        ivar0 + 1, dbin->getVariableNumber());
    return gres;
  }

  // Initializations

  int iptr_est  = -1;
  int iptr_std  = -1;
  int ndim = dbin->getNDim();
  int nvar = model->getVariableNumber();
  neighU = NeighUnique(ndim, false);

  /* Setting options */

  KrigingSystem ksys(dbin, dbout, model, &neighU);
  if (ksys.setKrigOptFlagGlobal(true)) return gres;
  if (! ksys.isReady()) return gres;

  /* Loop on the targets to be processed */

  int ng = 0;
  for (int iech_out = 0; iech_out < dbout->getSampleNumber(); iech_out++)
  {
    mes_process("Kriging sample", dbout->getSampleNumber(), iech_out);
    if (! dbout->isActive(iech_out)) continue;
    if (ksys.estimate(iech_out)) return gres;

    // Cumulate the R.H.S.

    VectorDouble rhs = ksys.getRHSC(ivar0);
    if (rhsCum.empty()) rhsCum.resize(rhs.size(),0.);
    ut_vector_add_inplace(rhsCum, rhs);
    ng++;
  }

  /* Preliminary checks */

  int ntot = dbin->getSampleNumber(false);
  int np   = dbin->getSampleNumber(true);
  double surface = ng * db_grid_maille(dbout);

  /* Average covariance over the territory */

  double cvv = model_cxx(model, dbout, dbout, ivar0, ivar0, 0,
                         db_epsilon_distance(dbin));

  /* Load the scaled cumulated R.H.S. in the array rhs */

  ut_vector_divide_inplace(rhsCum, (double) ng);

  /* Derive the kriging weights */

  int nred = ksys.getNRed();
  VectorDouble lhsinv = ksys.getLHSInv();
  VectorDouble zam = ksys.getZam();
  VectorDouble wgt(nred);
  matrix_product(nred, nred, nvar, lhsinv.data(), rhsCum.data(), wgt.data());

  /* Perform the estimation */

  double estim = ut_vector_inner_product(rhsCum, zam);
  double stdv = cvv - ut_vector_inner_product(rhsCum, wgt);
  stdv = (stdv > 0) ? sqrt(stdv) : 0.;
  double cvgeo = (estim == 0. || FFFF(estim)) ? TEST : stdv / estim;

  /* Store the results in the output Global_Res struture */

  gres.ntot = ntot;
  gres.np = np;
  gres.ng = ng;
  gres.surface = surface;
  gres.zest = estim;
  gres.sse = stdv;
  gres.cvgeo = cvgeo;
  gres.cvv = cvv;
  gres.weights = wgt;

  /* Printout */

  if (flag_verbose)
  {
    mestitle(1,"Global estimation kriging");
    message("Total number of data             = %d\n", ntot);
    message("Number of active data            = %d\n", np);
    message("Number of variables              = %d\n", nvar);
    message("Cvv                              = %lf\n", cvv);
    if (FFFF(estim))
      message("Estimation by kriging            = NA\n");
    else
    {
      message("Estimation by kriging            = %lf\n", estim);
    }
    message("Estimation St. Dev. of the mean  = %lf\n", stdv);
    if (FFFF(cvgeo))
      message("CVgeo                            = NA\n");
    else
      message("CVgeo                            = %lf\n", cvgeo);
    message("Surface                          = %lf\n", surface);
    if (FFFF(estim))
      message("Q (Estimation * Surface)         = NA\n");
    else
      message("Q (Estimation * Surface)         = %lf\n", estim * surface);
    message("\n");
  }
  return gres;
}

/****************************************************************************/
/*!
 **  Estimation of the variance by transitive method
 **
 ** \return  Error return code
 **
 ** \param[in]  dbgrid       Db structure containing the dicretization grid
 ** \param[in]  model        Model structure
 ** \param[in]  flag_verbose 1 for a verbose output
 ** \param[in]  flag_regular 1 for regular; 0 for stratified
 ** \param[in]  ndisc        Number of Discretization steps
 **
 ** \param[out]  abundance   Global estimated abundance
 ** \param[out]  sse         Global standard deviation
 ** \param[out]  cvtrans     CV transitive
 **
 *****************************************************************************/
int global_transitive(DbGrid *dbgrid,
                      Model *model,
                      int flag_verbose,
                      int flag_regular,
                      int ndisc,
                      double *abundance,
                      double *sse,
                      double *cvtrans)
{
  int idim, ndim, i, ix, iy, ix1, ix2, iy1, iy2, nx, ny, error, flag_value;
  double c00, cvv, dx, dy, dsum, gint, dsse, wtot;
  double value;
  VectorDouble d1;
  CovCalcMode mode;

  /* Initializations */

  error = 1;
  cvv = wtot = dsse = gint = dsum = 0.;
  flag_value = 0;
  st_global_init(dbgrid, dbgrid);
  if (st_check_environment(0, 1, model, NULL)) goto label_end;
  ndim = dbgrid->getNDim();
  d1.resize(2);

  if (ndim < 1 || ndim > 2)
  {
    messerr(
        "The transitive global estimation is implemented for 1 and 2 space only");
    goto label_end;
  }
  if (model->getVariableNumber() != 1)
  {
    messerr(
        "The transitive global estimation is implemented for 1 variable only");
    goto label_end;
  }
  if (!is_grid(dbgrid))
  {
    messerr(
        "The transitive global estimation requires a Db organized as a grid");
    goto label_end;
  }

  /* Core allocation */

  for (idim = 0; idim < ndim; idim++) d1[idim] = 0.;
  model_calcul_cov(NULL,model, mode, 1, 1., d1, &c00);

  /* Abundance estimation */

  flag_value = 0;
  if (dbgrid->getVariableNumber() == 1)
  {
    for (i = 0; i < dbgrid->getSampleNumber(); i++)
    {
      value = dbgrid->getVariable(i, 0);
      if (!FFFF(value)) dsum += value;
    }
    flag_value = 1;
  }

  /* 2-D case */

  if (ndim == 2)
  {
    dx = dbgrid->getDX(0);
    dy = dbgrid->getDX(1);
    nx = dbgrid->getNX(0);
    ny = dbgrid->getNX(1);
    if (flag_value) dsum *= dx * dy;

    /* Estimation */

    if (flag_regular)
    {

      /* Regular case */

      for (ix = -nx + 1; ix <= nx; ix++)
        for (iy = -ny + 1; iy <= ny; iy++)
        {
          d1[0] = dx * ix;
          d1[1] = dy * iy;
          model_calcul_cov(NULL,model, mode, 0, 1., d1, &dsse);
        }
      dsse *= dx * dy;
      // TODO : appeler model_integral
      // if (model_integral(model,ndisc,&gint)) goto label_end;
      *sse = dsse - gint;
    }
    else
    {

      /* Stratified case */

      for (ix1 = 0; ix1 < ndisc; ix1++)
        for (iy1 = 0; iy1 < ndisc; iy1++)
          for (ix2 = 0; ix2 < ndisc; ix2++)
            for (iy2 = 0; iy2 < ndisc; iy2++)
            {
              d1[0] = dx * (ix2 - ix1) / ndisc;
              d1[1] = dy * (iy2 - iy1) / ndisc;
              model_calcul_cov(NULL,model, mode, 0, 1., d1, &cvv);
              wtot += 1.;
            }
      cvv /= wtot;
      *sse = dx * dy * (c00 - cvv);
    }
  }
  else
  {

    /* 1-D case */

    dx = dbgrid->getDX(0);
    nx = dbgrid->getNX(0);
    if (flag_value) dsum *= dx;

    if (flag_regular)
    {

      /* Regular case */

      for (ix = -nx + 1; ix <= nx; ix++)
      {
        d1[0] = dx * ix;
        model_calcul_cov(NULL,model, mode, 0, 1., d1, &dsse);
      }
      dsse *= dx;
      // TODO: appeler model_integral
      // if (model_integral(model,ndisc,&gint)) goto label_end;
      *sse = dsse - gint;
    }
    else
    {

      /* Stratified case */

      for (ix1 = 0; ix1 < ndisc; ix1++)
        for (ix2 = 0; ix2 < ndisc; ix2++)
        {
          d1[0] = dx * (ix2 - ix1) / ndisc;
          model_calcul_cov(NULL,model, mode, 0, 1., d1, &cvv);
          wtot += 1.;
        }
      cvv /= wtot;
      *sse = dx * (c00 - cvv);
    }
  }

  if (flag_value)
  {
    *abundance = dsum;
    *cvtrans = ((*sse) <= 0.) ? TEST :
                                dsum / (*sse);
  }
  else
  {
    *abundance = *cvtrans = TEST;
  }
  (*sse) = (*sse > 0) ? sqrt(*sse) :
                        0.;

  /* Optional printout */

  if (flag_verbose)
  {
    if (flag_regular)
    {
      message("Transitive estimation (Regular case)\n");
      message("====================================\n");
      message("Space dimension           = %d \n", ndim);
      message("s * Sum[G(ks)]            = %lf\n", dsse);
      message("Integral[G(h)]            = %lf\n", gint);
    }
    else
    {
      message("Transitive estimation (Stratified case)\n");
      message("=======================================\n");
      message("Space dimension           = %d \n", ndim);
      message("G(0)                      = %lf\n", c00);
      message("G(s,s)                    = %lf\n", cvv);
    }
    message("Estimation St. Dev.       = %lf\n", (*sse));
    if (flag_value)
    {
      message("Global abundance          = %lf\n", (*abundance));
      if (FFFF(*cvtrans))
        message("Coefficient of Variation  = NA\n");
      else
        message("Coefficient of Variation  = %lf\n", (*cvtrans));
    }
  }

  /* Set the error return code */

  error = 0;

  label_end: return (error);
}

/****************************************************************************/
/*!
 **  Inverse distance estimation when Input DB is a grid
 **
 ** \param[in]  dbin        Input Db
 ** \param[in]  dbout       Output Db
 ** \param[in]  exponent    exponent of the inverse distance
 ** \param[in]  flag_expand 1 for expansion option
 **
 ** \param[out] indg        Working array
 ** \param[out] indref      Working array
 ** \param[out] coor        Working array
 ** \param[out] cooref      Working array
 **
 *****************************************************************************/
static void st_grid_invdist(DbGrid* dbin,
                            Db* dbout,
                            int exponent,
                            int flag_expand,
                            int *indg,
                            int *indref,
                            double *coor,
                            double *cooref)
{
  int idim, ndim, maxneigh, incorrect, ind, rank, iech_neigh;
  double result, total, val_neigh, dist, wgt, dmin;

  /* Initializations */

  DBIN = dbin;
  DBOUT = dbout;
  ndim = dbin->getNDim();
  maxneigh = (int) pow(2., (double) ndim);
  (void) db_extension_diag(dbout, &dmin);
  dmin /= 1.e5;

  /* Loop on the targets to be processed */

  for (IECH_OUT = 0; IECH_OUT < dbout->getSampleNumber(); IECH_OUT++)
  {
    mes_process("Estimation by Inverse distance", dbout->getSampleNumber(),
                IECH_OUT);
    if (!dbout->isActive(IECH_OUT)) continue;
    if (OptDbg::query(EDbg::KRIGING) || OptDbg::query(EDbg::NBGH) || OptDbg::query(EDbg::RESULTS))
    {
      mestitle(1, "Target location");
      db_sample_print(dbout, IECH_OUT, 1, 0, 0);
    }

    /* Find the grid index corresponding to the target */

    for (idim = 0; idim < ndim; idim++)
      cooref[idim] = dbout->getCoordinate(IECH_OUT, idim);
    point_to_grid(dbin, cooref, flag_expand, indref);

    /* Loop on the neighbors */

    result = total = 0.;
    for (rank = 0; rank < maxneigh; rank++)
    {
      for (idim = 0; idim < ndim; idim++)
        indg[idim] = indref[idim];

      /* Decompose the neighborhood rank */

      idim = 0;
      ind = rank;
      while (ind > 0)
      {
        if (ind % 2 == 1) indg[idim] += 1;
        ind /= 2;
        idim++;
      }

      /* Check that the neighboring point lies within the grid */

      for (idim = incorrect = 0; idim < ndim && incorrect == 0; idim++)
      {
        if (indg[idim] >= dbin->getNX(idim))
        {
          if (flag_expand)
            indg[idim]--;
          else
            incorrect = 1;
        }
      }

      /* Process the new neighboring point */

      if (incorrect)
      {
        result = TEST;
        break;
      }
      else
      {

        /* Check the value */

        iech_neigh = db_index_grid_to_sample(dbin, indg);
        val_neigh = dbin->getVariable(iech_neigh, 0);
        if (FFFF(val_neigh))
        {
          result = TEST;
          break;
        }
        else
        {

          /* Calculate the distance from neighborhood to target */

          dist = 0.;
          grid_to_point(dbin, indg, NULL, coor);
          dist = ut_distance(ndim, cooref, coor);
          if (dist < dmin)
          {
            result = val_neigh;
            total = 1.;
            break;
          }
          wgt = 1. / pow(dist, exponent);
          result += wgt * val_neigh;
          total += wgt;
        }
      }
    }
    if (!FFFF(result)) result /= total;
    dbout->setArray(IECH_OUT, IPTR_EST, result);
  }

  return;
}

/****************************************************************************/
/*!
 **  Inverse distance estimation when Input DB is a point file
 **
 ** \param[in]  dbin        Input Db
 ** \param[in]  dbout       Output Db
 ** \param[in]  exponent    exponent of the inverse distance
 ** \param[in]  dmax        Maximum search radius (only used for Point Db)
 **
 ** \param[out] coor        Working array
 ** \param[out] cooref      Working array
 **
 *****************************************************************************/
static void st_point_invdist(Db* dbin,
                             Db* dbout,
                             int exponent,
                             double dmax,
                             double *coor,
                             double *cooref)
{
  int idim, iech_in, ndim;
  double result, total, val_neigh, dist, wgt, dmin;

  /* Initializations */

  DBIN = dbin;
  DBOUT = dbout;
  ndim = dbin->getNDim();
  (void) db_extension_diag(dbout, &dmin);
  dmin /= 1.e5;

  /* Loop on the targets to be processed */

  for (IECH_OUT = 0; IECH_OUT < dbout->getSampleNumber(); IECH_OUT++)
  {
    mes_process("Estimation by Inverse distance", dbout->getSampleNumber(),
                IECH_OUT);
    if (!dbout->isActive(IECH_OUT)) continue;
    if (OptDbg::query(EDbg::KRIGING) || OptDbg::query(EDbg::NBGH) || OptDbg::query(EDbg::RESULTS))
    {
      mestitle(1, "Target location");
      db_sample_print(dbout, IECH_OUT, 1, 0, 0);
    }
    for (idim = 0; idim < ndim; idim++)
      cooref[idim] = dbout->getCoordinate(IECH_OUT, idim);

    /* Loop on the data points */

    result = total = 0.;
    for (iech_in = 0; iech_in < dbin->getSampleNumber(); iech_in++)
    {
      if (!dbin->isActive(iech_in)) continue;
      for (idim = 0; idim < ndim; idim++)
        coor[idim] = dbin->getCoordinate(iech_in, idim);
      val_neigh = dbin->getVariable(iech_in, 0);
      if (FFFF(val_neigh)) continue;

      /* Check that the data point is a valid neighbor */

      dist = ut_distance(ndim, coor, cooref);
      if (!FFFF(dmax) && dist > dmax) continue;

      /* Process the new neighboring point */

      if (dist < dmin)
      {
        result = val_neigh;
        total = 1.;
        break;
      }
      wgt = 1. / pow(dist, exponent);
      result += wgt * val_neigh;
      total += wgt;
    }
    if (!FFFF(result)) result /= total;
    dbout->setArray(IECH_OUT, IPTR_EST, result);
  }

  return;
}

/****************************************************************************/
/*!
 **  Inverse distance estimation
 **
 ** \return  Error return code
 **
 ** \param[in]  dbin        Input Db structure
 ** \param[in]  dbout       Output Db structure
 ** \param[in]  exponent    exponent of the inverse distance
 ** \param[in]  flag_expand 1 for expansion option
 ** \param[in]  dmax        Maximum search radius (used only for Points Db)
 **
 *****************************************************************************/
int invdist_f(Db *dbin, Db *dbout, int exponent, int flag_expand, double dmax)
{
  int *indg, *indref, error;
  double *coor, *cooref;
  DbGrid* dbgrid;

  /* Initializations */

  error = 1;
  indg = indref = nullptr;
  coor = cooref = nullptr;
  st_global_init(dbin, dbout);
  if (st_check_environment(1, 1, NULL, NULL)) goto label_end;

  /* Add the attribute for storing the result */

  IPTR_EST = dbout->addColumnsByConstant(1, 0.);
  if (IPTR_EST < 0) goto label_end;
  coor = db_sample_alloc(dbout, ELoc::X);
  if (coor == nullptr) goto label_end;
  cooref = db_sample_alloc(dbout, ELoc::X);
  if (cooref == nullptr) goto label_end;

  if (!is_grid(dbin))
  {
    st_point_invdist(dbin, dbout, exponent, dmax, coor, cooref);
  }
  else
  {
    dbgrid = dynamic_cast<DbGrid*>(dbin);
    indg = db_indg_alloc(dbgrid);
    if (indg == nullptr) goto label_end;
    indref = db_indg_alloc(dbgrid);
    if (indref == nullptr) goto label_end;
    st_grid_invdist(dbgrid, dbout, exponent, flag_expand, indg, indref, coor, cooref);
  }

  /* Set the error return code */

  error = 0;

  label_end: coor = db_sample_free(coor);
  cooref = db_sample_free(cooref);
  if (is_grid(DBIN))
  {
    indg = db_indg_free(indg);
    indref = db_indg_free(indref);
  }
  return (error);
}

/****************************************************************************/
/*!
 **  Returns the limits of the area of interest
 **
 ** \return  Error returned code
 **
 ** \param[in]  db      input Db structure
 ** \param[in]  top     Elevation of the Top variable
 ** \param[in]  bot     Elevation of the bottom variable
 **
 ** \param[out] ideb    Index of the starting sample
 ** \param[out] ifin    Index of the ending sample
 **
 *****************************************************************************/
static int st_get_limits(DbGrid* db, double top, double bot, int *ideb, int *ifin)
{
  int ndim, nz, iad;
  double z0, dz;

  /* Initializations */

  ndim = db->getNDim();
  z0 = db->getX0(ndim - 1);
  nz = db->getNX(ndim - 1);
  dz = db->getDX(ndim - 1);

  /* Preliminary checks */

  if (!FFFF(bot) && !FFFF(top) && top < bot)
  {
    messerr("Error: Top(%lf) must be larger than Bottom (%lf)", top, bot);
    return (1);
  }

  if (FFFF(bot))
    iad = 0;
  else
    iad = (int) ((bot - z0) / dz);
  if (bot > z0 + iad * dz) iad++;
  *ideb = MAX(0, MIN(iad, nz-1));

  if (FFFF(top))
    iad = nz - 1;
  else
    iad = (int) ((top - z0) / dz);
  *ifin = MAX(0, MIN(iad, nz-1));
  return (0);
}

/****************************************************************************/
/*!
 **  Definition of the neighborhood
 **
 ** \return  Error return code: 1 if the target does not belong to the
 ** \return  area of interest
 **
 ** \param[in]  ideb          Index of the starting sample
 ** \param[in]  ifin          Index of the ending sample
 ** \param[in]  neigh_radius  Radius of the Neighborhood
 **
 ** \param[out] status        Neighborhood error status
 ** \param[out] nbefore       Number of samples in neighborhood before target
 ** \param[out] nafter        Number of samples in neighborhood after target
 **
 *****************************************************************************/
static int st_get_neigh(int ideb,
                        int ifin,
                        int neigh_radius,
                        int *status,
                        int *nbefore,
                        int *nafter)
{
  int iad;

  *status = 1;
  if (IECH_OUT < ideb || IECH_OUT > ifin) return (1);

  iad = MAX(IECH_OUT - neigh_radius, ideb);
  *nbefore = IECH_OUT - iad;

  iad = MIN(IECH_OUT + neigh_radius, ifin);
  *nafter = iad - IECH_OUT;

  *status = 0;
  return (0);
}

/****************************************************************************/
/*!
 **  Calculate the discretized covariance
 **
 ** \param[in]  dist          Integer distance
 ** \param[in]  cov           Array of discretized covariances
 ** \param[in]  cov_radius    Radius of the covariance array
 ** \param[in]  flag_sym      1 for symmetrized covariance
 **
 *****************************************************************************/
static double st_cov_exp(int dist, double *cov, int cov_radius, int flag_sym)
{
  double val1, val2, val;

  if (flag_sym)
  {
    val1 = cov[cov_radius - dist];
    val2 = cov[cov_radius + dist];
    val = (val1 + val2) / 2.;
  }
  else
  {
    val = cov[cov_radius + dist];
  }
  return (val);
}

/****************************************************************************/
/*!
 **  Establish the L.H.S. of the Kriging system
 **  in the case of the discretized covariances
 **
 ** \param[in]  covdd         Array of discretized covariance (data-data)
 ** \param[in]  cov_radius    Radius of the covariance array
 ** \param[in]  flag_sym      1 for symmetrized covariance
 ** \param[in]  nfeq          0 or 1 drift function(s)
 ** \param[in]  nbefore       Number of samples in neighborhood before target
 ** \param[in]  nafter        Number of samples in neighborhood after target
 ** \param[in]  neq           Number of kriging equations
 **
 *****************************************************************************/
static void st_lhs_exp(double *covdd,
                       int cov_radius,
                       int flag_sym,
                       int nfeq,
                       int nbefore,
                       int nafter,
                       int neq)
{
  int i, j;

  /* Covariance part */

  for (i = -nbefore; i <= nafter; i++)
    for (j = -nbefore; j <= nafter; j++)
    {
      LHS_EXP(i+nbefore,j+nbefore) = st_cov_exp(i - j, covdd, cov_radius,
                                                flag_sym);
      LHS_EXP(j+nbefore,i+nbefore) = st_cov_exp(j - i, covdd, cov_radius,
                                                flag_sym);
    }

  /* Drift part */

  if (nfeq == 0) return;
  for (i = -nbefore; i <= nafter; i++)
  {
    LHS_EXP(i+nbefore,neq-1) = 1.;
    LHS_EXP(neq-1,i+nbefore) = 1.;
  }
  LHS_EXP(neq-1,neq-1) = 0.;
  return;
}

/****************************************************************************/
/*!
 **  Establish the R.H.S. of the Kriging system in the case
 **  of the discretized covariances
 **
 ** \param[in]  covd0         Array of discretized covariance (data-data)
 ** \param[in]  cov_radius    Radius of the covariance array
 ** \param[in]  flag_sym      1 for symmetrized covariance
 ** \param[in]  nfeq          0 or 1 drift function(s)
 ** \param[in]  nbefore       Number of samples in neighborhood before target
 ** \param[in]  nafter        Number of samples in neighborhood after target
 ** \param[in]  neq           Number of equations
 **
 *****************************************************************************/
static void st_rhs_exp(double *covd0,
                       int cov_radius,
                       int flag_sym,
                       int nfeq,
                       int nbefore,
                       int nafter,
                       int neq)
{
  int i;

  /* Covariance part */

  for (i = -nbefore; i <= nafter; i++)
    RHS_EXP(i+nbefore) = st_cov_exp(i, covd0, cov_radius, flag_sym);

  /* Drift part */

  if (nfeq == 0) return;
  RHS_EXP(neq-1) = 1.;
  return;
}

/****************************************************************************/
/*!
 **  Perform the estimation for the Factorial Kriging Analysis
 **  in the case of the discretized covariances
 **
 ** \param[in]  db            Db structure
 ** \param[in]  wgt           Array containing the kriging weights
 ** \param[in]  nbefore       Number of samples in neighborhood before target
 ** \param[in]  nafter        Number of samples in neighborhood after target
 **
 *****************************************************************************/
static double st_estim_exp(Db *db, double *wgt, int nbefore, int nafter)
{
  int i;
  double result;

  /* Perform the estimation */

  result = 0.;
  for (i = -nbefore; i <= nafter; i++)
    result += wgt[i + nbefore] * db->getVariable(IECH_OUT + i, 0);

  return (result);
}

/****************************************************************************/
/*!
 **  Factorial Kriging analysis on a 1-D grid file using
 **  discretized covariances for total and partial variables
 **
 ** \return  Error return code
 **
 ** \param[in]  db            input Db structure
 ** \param[in]  covdd         Array of discretized cov. for total variable
 ** \param[in]  covd0         Array of discretized cov. for partial variable
 ** \param[in]  top           Elevation of the Top variable
 ** \param[in]  bot           Elevation of the bottom variable
 ** \param[in]  cov_radius    Radius of the Covariance arrays
 ** \param[in]  neigh_radius  Radius of the Neighborhood
 ** \param[in]  flag_sym      1 for symmetrized covariance
 ** \param[in]  nfeq          0 or 1 drift function(s)
 **
 *****************************************************************************/
int anakexp_f(DbGrid *db,
              double *covdd,
              double *covd0,
              double top,
              double bot,
              int cov_radius,
              int neigh_radius,
              int flag_sym,
              int nfeq)
{
  int i, ndim, nvarin, nech, size, error, ideb, ifin, neq, status;
  int nbefore, nafter, nbefore_mem, nafter_mem;
  double result;
  VectorInt ranks;

  /* Initializations */

  error = 1;
  st_global_init(db, db);
  FLAG_EST = 1;
  lhs = rhs = wgt = nullptr;
  ndim = db->getNDim();
  nvarin = db->getVariableNumber();
  nbefore_mem = nafter_mem = -1;
  size = nech = 0;

  /* Prepare the Koption structure */

  if (krige_koption_manage(1, 1, EKrigOpt::PONCTUAL, 1, VectorInt()))
    return (1);

  /* Preliminary checks */

  if (ndim != 1 || !is_grid(db))
  {
    messerr("This procedure is limited to 1-D grid");
    goto label_end;
  }
  if (nvarin != 1)
  {
    messerr("This procedure is limited to the monovariate case");
    goto label_end;
  }
  nech = db->getNX(ndim - 1);
  if (nfeq != 0 && nfeq != 1)
  {
    messerr("This procedure is limited to Stationary or Intrinsic case");
    messerr("The argument 'nfeq' must be 0 or 1");
    goto label_end;
  }
  if (neigh_radius > cov_radius / 2)
  {
    messerr("The radius of the neighborhood (%d) must be smaller or equal",
            neigh_radius);
    messerr("to the radius of the covariance (%d)", cov_radius);
    goto label_end;
  }

  /* Add the attribute for storing the result */

  IPTR_EST = db->addColumnsByConstant(nvarin, 0.);
  if (IPTR_EST < 0) goto label_end;
  DBOUT = db;

  /* Core allocation */

  size = 2 * neigh_radius + 1;
  st_krige_manage_basic(1, size, size, 1, nfeq);
  ranks.resize(nech);
  for (i = 0; i < nech; i++)
  {
    ranks[i] = i;
    flag[i] = 1;
  }

  /* Get the limits of the area to be processed */

  if (st_get_limits(db, top, bot, &ideb, &ifin)) goto label_end;

  /* Loop on the grid nodes */

  status = 0;
  for (IECH_OUT = 0; IECH_OUT < nech; IECH_OUT++)
  {
    mes_process("Factorial Kriging Analysis", nech, IECH_OUT);
    OptDbg::setIndex(IECH_OUT + 1);
    if (!db->isActive(IECH_OUT)) continue;
    if (OptDbg::query(EDbg::KRIGING) || OptDbg::query(EDbg::NBGH) || OptDbg::query(EDbg::RESULTS))
    {
      mestitle(1, "Target location");
      db_sample_print(db, IECH_OUT, 1, 0, 0);
    }

    /* Discard the grid nodes which doe not belong to the processed area */

    DBOUT->setArray(IECH_OUT, IPTR_EST, TEST);

    /* Look for the neighborhood */

    if (st_get_neigh(ideb, ifin, neigh_radius, &status, &nbefore, &nafter))
      continue;

    /* If the neighborhood has changed, establish the kriging system */

    neq = nafter + nbefore + 1;
    if (nfeq == 1) neq++;
    if (nbefore_mem != nbefore || nafter_mem != nafter || OptDbg::force())
    {
      nbefore_mem = nbefore;
      nafter_mem = nafter;

      /* Establish the L.H.S. of the kriging system */

      st_lhs_exp(covdd, cov_radius, flag_sym, nfeq, nbefore, nafter, neq);
      if (OptDbg::query(EDbg::KRIGING)) krige_lhs_print(nech, neq, neq, flag, lhs);

      /* Invert the kriging system */

      if (matrix_invert(lhs, neq, IECH_OUT))
      {
        status = 1;
        continue;
      }

      /* Establish the R.H.S. of the kriging system */

      st_rhs_exp(covd0, cov_radius, flag_sym, nfeq, nbefore, nafter, neq);
      if (OptDbg::query(EDbg::KRIGING))
        krige_rhs_print(nvarin, nech, neq, neq, flag, rhs);

      /* Derive the kriging weights */

      matrix_product(neq, neq, 1, lhs, rhs, wgt);
    }

    /* Calculate the estimation */

    result = st_estim_exp(db, wgt, nbefore, nafter);
    DBOUT->setArray(IECH_OUT, IPTR_EST, result);
    if (OptDbg::query(EDbg::RESULTS)) st_result_kriging_print(0, nvarin, status);
  }

  /* Set the error return flag */

  error = 0;

  label_end: OptDbg::setIndex(0);
  (void) krige_koption_manage(-1, 1, EKrigOpt::PONCTUAL, 1, VectorInt());
  st_krige_manage_basic(-1, size, size, 1, nfeq);
  return (error);
}

/****************************************************************************/
/*!
 **  Calculate the experimental covariance of the residual variable
 **  defined on the grid
 **
 ** \param[in]  db            input Db structure
 ** \param[in]  model         Model describing the horizontal structure
 ** \param[in]  cov_ref       Array of discretized covariance for target variable
 ** \param[in]  cov_radius    Radius of the covariance array
 ** \param[in]  flag_sym      1 for symmetrized covariance
 ** \param[in]  cov_ss        Array of dimensions of the Covariance array
 ** \param[in]  cov_nn        Array of radius of the Covariance array
 **
 ** \param[out] cov_res       Array containing the covariance of the residual
 **                           variable
 **
 *****************************************************************************/
static void st_calculate_covres(DbGrid *db,
                                Model *model,
                                double *cov_ref,
                                int cov_radius,
                                int flag_sym,
                                int cov_ss[3],
                                int cov_nn[3],
                                double *cov_res)
{
  double dx, dy, c00, covtot, covtab, covver;
  int i, ix, iy, iz;
  VectorDouble d1;
  CovCalcMode mode;

  /* Initializations */

  d1.resize(3);
  dx = db->getDX(0);
  dy = db->getDX(1);
  covtot = COV_REF(0);
  for (i = 0; i < 3; i++)
    d1[i] = 0.;
  model_calcul_cov(NULL,model, mode, 1, 1., d1, &c00);

  /* Evaluate the array of experimental covariance of the residual variable */

  for (ix = -cov_nn[0]; ix <= cov_nn[0]; ix++)
    for (iy = -cov_nn[1]; iy <= cov_nn[1]; iy++)
      for (iz = -cov_nn[2]; iz <= cov_nn[2]; iz++)
      {
        if (!flag_sym)
          covver = COV_REF(iz);
        else
          covver = (COV_REF(iz) + COV_REF(-iz)) / 2.;
        d1[0] = dx * ix;
        d1[1] = dy * iy;
        model_calcul_cov(NULL,model, mode, 1, 1., d1, &covtab);
        COV_RES(ix,iy,iz)= covver * (covtab + covtot - c00) / covtot;
      }

  return;
}

/****************************************************************************/
/*!
 **  Calculate the experimental covariance of the total variable
 **  defined on the grid
 **
 ** \param[in]  db            input Db structure
 ** \param[in]  ix0           index of the grid index along X
 ** \param[in]  iy0           index of the grid index along Y
 ** \param[in]  flag_sym      1 for symmetrized covariance
 ** \param[in]  cov_ss        Array of dimensions of the Covariance array
 ** \param[in]  cov_nn        Array of radius of the Covariance array
 **
 ** \param[out] num_tot       Array containing the numb er of pairs
 ** \param[out] cov_tot       Array containing the covariance of the total
 **                           variable
 **
 *****************************************************************************/
static void st_calculate_covtot(DbGrid *db,
                                int ix0,
                                int iy0,
                                int flag_sym,
                                int cov_ss[3],
                                int cov_nn[3],
                                int *num_tot,
                                double *cov_tot)
{
  int ix, iy, iz, ix1, iy1, iz1, jx1, jy1, jz1, jx2, jy2, jz2, indg[3];
  int idx, idy, idz, jdx, jdy, iad, jad;
  double val1, val2, val, ratio;

  /* Initialization */

  for (ix = -cov_nn[0]; ix <= cov_nn[0]; ix++)
    for (iy = -cov_nn[1]; iy <= cov_nn[1]; iy++)
      for (iz = -cov_nn[2]; iz <= cov_nn[2]; iz++)
      {
        COV_TOT(ix,iy,iz)= 0.;
        NUM_TOT(ix,iy,iz) = 0;
      }

      /* Loop on the first point */

  for (iz1 = 0; iz1 < db->getNX(2); iz1++)
    for (iy1 = -cov_nn[1]; iy1 <= cov_nn[1]; iy1++)
      for (ix1 = -cov_nn[0]; ix1 <= cov_nn[0]; ix1++)
      {
        jx1 = ix0 + ix1;
        if (jx1 < 0 || jx1 >= db->getNX(0)) continue;
        jy1 = iy0 + iy1;
        if (jy1 < 0 || jy1 >= db->getNX(1)) continue;
        jz1 = iz1;

        indg[0] = jx1;
        indg[1] = jy1;
        indg[2] = jz1;
        iad = db_index_grid_to_sample(db, indg);
        if (!db->isActive(iad)) continue;
        val1 = db->getVariable(iad, 0);
        if (FFFF(val1)) continue;

        /* Loop on the second point within the covariance array */

        for (idz = -cov_nn[2]; idz <= cov_nn[2]; idz++)
          for (idy = -cov_nn[1]; idy <= cov_nn[1]; idy++)
            for (idx = -cov_nn[0]; idx <= cov_nn[0]; idx++)
            {
              jx2 = jx1 + idx;
              if (jx2 < 0 || jx2 >= db->getNX(0)) continue;
              jy2 = jy1 + idy;
              if (jy2 < 0 || jy2 >= db->getNX(1)) continue;
              jz2 = jz1 + idz;
              if (jz2 < 0 || jz2 >= db->getNX(2)) continue;

              jdx = jx2 - ix0;
              if (jdx < -cov_nn[0] || jdx > cov_nn[0]) continue;
              jdy = jy2 - iy0;
              if (jdy < -cov_nn[1] || jdy > cov_nn[1]) continue;

              indg[0] = jx2;
              indg[1] = jy2;
              indg[2] = jz2;
              jad = db_index_grid_to_sample(db, indg);
              if (!db->isActive(jad)) continue;
              val2 = db->getVariable(jad, 0);
              if (FFFF(val2)) continue;

              /* Update the Covariance */

              COV_TOT(idx,idy,idz)+= val1 * val2;
              NUM_TOT(idx,idy,idz)+= 1;
            }
          }

          /* Scaling */

  ratio = NUM_TOT(0, 0, 0);
  for (ix = -cov_nn[0]; ix <= cov_nn[0]; ix++)
    for (iy = -cov_nn[1]; iy <= cov_nn[1]; iy++)
      for (iz = -cov_nn[2]; iz <= cov_nn[2]; iz++)
      {
        if (NUM_TOT(ix,iy,iz)<= 0.)
        {
          COV_TOT(ix,iy,iz) = TEST;
        }
        else
        {
          COV_TOT(ix,iy,iz) /= ratio;
        }
      }

      /* Symmetry */

  for (ix = -cov_nn[0]; ix < 0; ix++)
    for (iy = -cov_nn[1]; iy <= cov_nn[1]; iy++)
      for (iz = -cov_nn[2]; iz <= cov_nn[2]; iz++)
      {
        val1 = COV_TOT(ix, iy, iz);
        val2 = COV_TOT(-ix, iy, iz);
        val = (FFFF(val1) || FFFF(val2)) ? TEST :
                                           (val1 + val2) / 2.;
        COV_TOT( ix,iy,iz)= COV_TOT(-ix,iy,iz) = val;
      }

  for (ix = -cov_nn[0]; ix <= cov_nn[0]; ix++)
    for (iy = -cov_nn[1]; iy < 0; iy++)
      for (iz = -cov_nn[2]; iz <= cov_nn[2]; iz++)
      {
        val1 = COV_TOT(ix, -iy, iz);
        val2 = COV_TOT(ix, iy, iz);
        val = (FFFF(val1) || FFFF(val2)) ? TEST :
                                           (val1 + val2) / 2.;
        COV_TOT(ix, iy,iz)= COV_TOT(ix,-iy,iz) = val;
      }

  if (flag_sym) for (ix = -cov_nn[0]; ix <= cov_nn[0]; ix++)
    for (iy = -cov_nn[1]; iy <= cov_nn[1]; iy++)
      for (iz = -cov_nn[2]; iz < 0; iz++)
      {
        val1 = COV_TOT(ix, iy, -iz);
        val2 = COV_TOT(ix, iy, iz);
        val = (FFFF(val1) || FFFF(val2)) ? TEST :
                                           (val1 + val2) / 2.;
        COV_TOT(ix,iy,-iz)= COV_TOT(ix,iy, iz) = val;
}

  return;
}

/****************************************************************************/
/*!
 **  Find the neighborhood of a pixel
 **
 ** \param[in]  db            input Db structure
 ** \param[in]  ix0           index of the pixel along X
 ** \param[in]  iy0           index of the pixel along Y
 ** \param[in]  iz0           index of the pixel along Z
 ** \param[in]  nei_ss        Array of dimensions of the Neighborhood
 ** \param[in]  nei_nn        Array of radius of the Neighborhood
 **
 ** \param[out] nei_cur       Array containing the neighborhood
 **
 *****************************************************************************/
static VectorInt st_neigh_find(DbGrid *db,
                               int ix0,
                               int iy0,
                               int iz0,
                               int nei_ss[3],
                               int nei_nn[3],
                               int *nei_cur)
{
  int ix, iy, iz, jx, jy, jz, indg[3], number, locrank;
  VectorInt nbgh_ranks;

  /* Loop on the pixels of the neighborhood */

  number = 0;
  for (ix = -nei_nn[0]; ix <= nei_nn[0]; ix++)
    for (iy = -nei_nn[1]; iy <= nei_nn[1]; iy++)
      for (iz = -nei_nn[2]; iz <= nei_nn[2]; iz++)
      {
        NEI_CUR(ix,iy,iz)= -1;
        jx = ix0 + ix;
        if (jx < 0 || jx >= db->getNX(0)) continue;
        jy = iy0 + iy;
        if (jy < 0 || jy >= db->getNX(1)) continue;
        jz = iz0 + iz;
        if (jz < 0 || jz >= db->getNX(2)) continue;
        indg[0] = jx;
        indg[1] = jy;
        indg[2] = jz;
        locrank = db_index_grid_to_sample(db,indg);
        if (FFFF(db->getVariable(locrank,0))) continue;
        NEI_CUR(ix,iy,iz) = locrank;
        nbgh_ranks.push_back(locrank);
        flag[number] = 1;
        number++;
      }

      /* Define the returned argument */

  return nbgh_ranks;
}

/****************************************************************************/
/*!
 **  Check if two neighborhood patterns are similar
 **
 ** \return  1 if the patterns are different; 0 otherwise
 **
 ** \param[in]  nei_ss        Array of dimensions of the Neighborhood
 ** \param[in]  nei_nn        Array of radius of the Neighborhood
 ** \param[in]  nei_ref       Array containing the reference neighborhood
 ** \param[in]  nei_cur       Array containing the current neighborhood
 **
 *****************************************************************************/
static int st_neigh_diff(int nei_ss[3],
                         int nei_nn[3],
                         int *nei_ref,
                         int *nei_cur)
{
  int ix, iy, iz, flag1, flag2, flag_diff;

  /* Loop on the pixels of the neighborhood */

  flag_diff = 1;
  for (ix = -nei_nn[0]; ix <= nei_nn[0]; ix++)
    for (iy = -nei_nn[1]; iy <= nei_nn[1]; iy++)
      for (iz = -nei_nn[2]; iz <= nei_nn[2]; iz++)
      {
        flag1 = NEI_REF(ix,iy,iz)< 0;
        flag2 = NEI_CUR(ix,iy,iz) < 0;
        if (flag1 != flag2) goto label_end;
      }
  flag_diff = 0;

  label_end:

  /* Copy the current neighborhood into the reference neighborhood */

  for (ix = -nei_nn[0]; ix <= nei_nn[0]; ix++)
    for (iy = -nei_nn[1]; iy <= nei_nn[1]; iy++)
      for (iz = -nei_nn[2]; iz <= nei_nn[2]; iz++)
        NEI_REF(ix,iy,iz)= NEI_CUR(ix,iy,iz);

  return (flag_diff);
}

/****************************************************************************/
/*!
 **  Establish the kriging L.H.S. using discretized covariances
 **
 ** \param[in]  nech          Number of samples in the Neighborhood
 ** \param[in]  nfeq          Number of drift functions
 ** \param[in]  nei_ss        Array of dimensions of the Neighborhood
 ** \param[in]  nei_nn        Array of radius of the Neighborhood
 ** \param[in]  cov_ss        Array of dimensions of the Covariance array
 ** \param[in]  cov_nn        Array of radius of the Covariance array
 ** \param[in]  nei_cur       Array containing the current neighborhood
 ** \param[in]  cov_tot       Array containing the total variable covariance
 ** \param[in]  nugget        Amount of additional Nugget Effect
 **
 *****************************************************************************/
static void st_lhs_exp_3D(int nech,
                          int nfeq,
                          int nei_ss[3],
                          int nei_nn[3],
                          int cov_ss[3],
                          int cov_nn[3],
                          int *nei_cur,
                          double *cov_tot,
                          double nugget)
{
  int ix, iy, iz, jx, jy, jz, i, j, neq;
  double value;

  /* Initializations */

  neq = nech + nfeq;

  /* Covariance part of the L.H.S. */

  i = 0;
  for (ix = -nei_nn[0]; ix <= nei_nn[0]; ix++)
    for (iy = -nei_nn[1]; iy <= nei_nn[1]; iy++)
      for (iz = -nei_nn[2]; iz <= nei_nn[2]; iz++)
      {
        if (NEI_CUR(ix,iy,iz)< 0) continue;

        j = 0;
        for (jx=-nei_nn[0]; jx<=nei_nn[0]; jx++)
        for (jy=-nei_nn[1]; jy<=nei_nn[1]; jy++)
        for (jz=-nei_nn[2]; jz<=nei_nn[2]; jz++)
        {
          if (NEI_CUR(jx,jy,jz) < 0) continue;
          value = COV_TOT(ix-jx,iy-jy,iz-jz);
          LHS_EXP(i,j) = LHS_EXP(j,i) = value;
          if (i == j) LHS_EXP(i,j) += nugget;
          j++;
        }
        i++;
      }

      /* Drift part */

  if (nfeq == 0) return;
  for (i = 0; i < nech; i++)
  {
    LHS_EXP(i,neq-1) = 1.;
    LHS_EXP(neq-1,i) = 1.;
  }
  LHS_EXP(neq-1,neq-1) = 0.;

  return;
}

/****************************************************************************/
/*!
 **  Establish the kriging R.H.S. using discretized covariances
 **
 ** \param[in]  nech          Number of samples in the neighborhood
 ** \param[in]  nfeq          Number of drift functions
 ** \param[in]  nei_ss        Array of dimensions of the Neighborhood
 ** \param[in]  nei_nn        Array of radius of the Neighborhood
 ** \param[in]  cov_ss        Array of dimensions of the Covariance array
 ** \param[in]  cov_nn        Array of radius of the Covariance array
 ** \param[in]  nei_cur       Array containing the current neighborhood
 ** \param[in]  cov_res       Array containing the residual variable covariance
 **
 *****************************************************************************/
static void st_rhs_exp_3D(int nech,
                          int nfeq,
                          int nei_ss[3],
                          int nei_nn[3],
                          int cov_ss[3],
                          int cov_nn[3],
                          int *nei_cur,
                          double *cov_res)
{
  int ix, iy, iz, neq, i;

  /* Initializations */

  neq = nech + nfeq;

  /* Covariance part of the R.H.S. */

  i = 0;
  for (ix = -nei_nn[0]; ix <= nei_nn[0]; ix++)
    for (iy = -nei_nn[1]; iy <= nei_nn[1]; iy++)
      for (iz = -nei_nn[2]; iz <= nei_nn[2]; iz++)
      {
        if (NEI_CUR(ix,iy,iz)< 0) continue;
        RHS_EXP(i) = COV_RES(ix,iy,iz);
        i++;
      }

      /* Drift part */

  if (nfeq == 0) return;
  RHS_EXP(neq-1) = 1.;

  return;
}

/****************************************************************************/
/*!
 **  Evaluate the Factorial Kriging estimate
 **
 ** \return  The estimation result
 **
 ** \param[in]  db            input Db structure
 ** \param[in]  nei_ss        Array of dimensions of the Neighborhood
 ** \param[in]  nei_nn        Array of radius of the Neighborhood
 ** \param[in]  nei_cur       Array containing the current neighborhood
 ** \param[in]  weight        Array of Kriging weights
 **
 *****************************************************************************/
static double st_estim_exp_3D(Db *db,
                              int nei_ss[3],
                              int nei_nn[3],
                              int *nei_cur,
                              double *weight)
{
  int i, ix, iy, iz;
  double result;

  /* Initializations */

  result = 0.;
  i = 0;
  for (ix = -nei_nn[0]; ix <= nei_nn[0]; ix++)
    for (iy = -nei_nn[1]; iy <= nei_nn[1]; iy++)
      for (iz = -nei_nn[2]; iz <= nei_nn[2]; iz++)
      {
        if (NEI_CUR(ix,iy,iz)< 0) continue;
        result += weight[i] * db->getVariable(NEI_CUR(ix,iy,iz),0);
        i++;
      }

  return (result);
}

/****************************************************************************/
/*!
 **  Dump the contents of the covariance maps
 **
 ** \param[in]  file     FILE structure where the dmp must be produced
 ** \param[in]  ix0      Rank of the trace along X (-1 for the reference)
 ** \param[in]  iy0      Rank of the trace along Y (-1 for the reference)
 ** \param[in]  cov_ss   Array of dimensions of the Covariance array
 ** \param[in]  cov_nn   Array of radius of the Covariance array
 ** \param[out] num_tot  Array containing the numb er of pairs
 ** \param[out] cov_tot  Array containing the covariance of the total variable
 **
 *****************************************************************************/
static void st_vario_dump(FILE *file,
                          int ix0,
                          int iy0,
                          int cov_ss[3],
                          int cov_nn[3],
                          int *num_tot,
                          double *cov_tot)
{
  int ix, iy, iz, num;
  double cov;

  fprintf(file, "*%3d %3d\n", ix0, iy0);

  for (ix = -cov_nn[0]; ix <= cov_nn[0]; ix++)
    for (iy = -cov_nn[1]; iy <= cov_nn[1]; iy++)
      for (iz = -cov_nn[2]; iz <= cov_nn[2]; iz++)
      {
        num = (num_tot == nullptr) ? 0 : NUM_TOT(ix, iy, iz);
        cov = COV_TOT(ix, iy, iz);
        fprintf(file, "%3d %3d %3d %3d %lf\n", ix, iy, iz, num, cov);
      }
  return;
}

/****************************************************************************/
/*!
 **  Factorial Kriging analysis on a grid file using discretized
 **  covariances for the target variable.
 **  The discretized covariance of the total variable is calculated on the fly
 **
 ** \return  Error return code
 **
 ** \param[in]  db            input Db structure
 ** \param[in]  cov_ref       Array of discretized covariance for target variable
 ** \param[in]  cov_radius    Radius of the covariance array
 ** \param[in]  neigh_ver     Radius of the Neighborhood along Vertical
 ** \param[in]  neigh_hor     Radius of the Neighborhood along Horizontal
 ** \param[in]  flag_sym      1 for symmetrized covariance
 ** \param[in]  model         Model structure (only used for horizontal)
 ** \param[in]  nugget        Additional Nugget Effect component
 ** \param[in]  nfeq          0 or 1 drift function(s)
 ** \param[in]  dbg_ix        Rank of the trace along X for variogram debug
 ** \param[in]  dbg_iy        Rank of the trace along Y for variogram debug
 **
 ** \remark  The discretized covariance of the target variable is provided
 ** \remark  in 1-D along the vertical. Its extension to the space dimension
 ** \remark  is performed using the theoretical factorized model
 **
 ** \remark  If dbg_ix < -1 || dbg_iy < -1, no variogram debug file is created
 **
 *****************************************************************************/
int anakexp_3D(DbGrid *db,
               double *cov_ref,
               int cov_radius,
               int neigh_ver,
               int neigh_hor,
               int flag_sym,
               Model *model,
               double nugget,
               int nfeq,
               int dbg_ix,
               int dbg_iy)
{
  int i, ix, iy, iz, ndim, nvarin, nech, error, neq, status, ecr;
  int size_cov, size_nei, flag_new, flag_col;
  int cov_ss[3], cov_nn[3], nei_ss[3], nei_nn[3], indg[3];
  int *num_tot, *nei_cur, *nei_ref;
  double *cov_tot, *cov_res, result;
  FILE *fildmp;
  VectorInt nbgh_ranks;

  /* Initializations */

  error = 1;
  st_global_init(db, db);
  FLAG_EST = 1;
  fildmp = nullptr;
  cov_tot = cov_res = nullptr;
  num_tot = nei_cur = nei_ref = nullptr;
  lhs = rhs = wgt = nullptr;
  ndim = db->getNDim();
  nvarin = db->getVariableNumber();
  size_nei = 0;

  /* Prepare the Koption structure */

  if (krige_koption_manage(1, 1, EKrigOpt::PONCTUAL, 1, VectorInt()))
    return (1);

  /* Preliminary checks */

  if (ndim != 3 || !is_grid(db))
  {
    messerr("This procedure is limited to 3-D grid");
    goto label_end;
  }
  if (nvarin != 1)
  {
    messerr("This procedure is limited to the monovariate case");
    goto label_end;
  }
  if (nfeq != 0 && nfeq != 1)
  {
    messerr("This procedure is limited to Stationary or Intrinsic case");
    messerr("The argument 'nfeq' must be 0 or 1");
    goto label_end;
  }
  if (neigh_ver > cov_radius / 2)
  {
    messerr("The radius of the neighborhood (%d) must be smaller or equal",
            neigh_ver);
    messerr("to the radius of the covariance (%d)", cov_radius);
    goto label_end;
  }

  /* Open the Variogram debugging file */

  if (dbg_ix >= -1 && dbg_ix < db->getNX(0) && dbg_iy >= -1
      && dbg_iy < db->getNX(1))
  {
    fildmp = gslFopen("Vario.dat", "w");
    if (fildmp == nullptr) goto label_end;
  }

  /* Add the attribute for storing the result */

  IPTR_EST = db->addColumnsByConstant(nvarin, 0.);
  if (IPTR_EST < 0) goto label_end;
  DBOUT = db;

  /* Define essential variables */

  nei_nn[0] = MIN(db->getNX(0) - 1, neigh_hor);
  nei_nn[1] = MIN(db->getNX(1) - 1, neigh_hor);
  nei_nn[2] = MIN(db->getNX(2) - 1, neigh_ver);
  size_nei = size_cov = 1;
  for (i = 0; i < db->getNDim(); i++)
  {
    nei_ss[i] = 2 * nei_nn[i] + 1;
    cov_nn[i] = 2 * nei_nn[i];
    cov_ss[i] = 2 * cov_nn[i] + 1;
    size_nei *= nei_ss[i];
    size_cov *= cov_ss[i];
  }
  size_nei += nfeq;

  /* Core allocation */

  num_tot = st_icore(size_cov, 1);
  if (num_tot == nullptr) goto label_end;
  nei_cur = st_icore(size_nei, 1);
  if (nei_cur == nullptr) goto label_end;
  nei_ref = st_icore(size_nei, 1);
  if (nei_ref == nullptr) goto label_end;
  cov_tot = st_core(size_cov, 1);
  if (cov_tot == nullptr) goto label_end;
  cov_res = st_core(size_cov, 1);
  if (cov_res == nullptr) goto label_end;
  st_krige_manage_basic(1, size_nei, size_nei, 1, nfeq);
  for (i = 0; i < size_nei; i++)
    nei_ref[i] = -1;

  /* Calculate the discretized covariance of residual variable */

  st_calculate_covres(db, model, cov_ref, cov_radius, flag_sym, cov_ss, cov_nn,
                      cov_res);
  if (dbg_ix == -1 && dbg_iy == -1)
    st_vario_dump(fildmp, -1, -1, cov_ss, cov_nn, nullptr, cov_res);

  /* Loop on the grid nodes */

  status = nech = neq = 0;
  IECH_OUT = ecr = 0;
  for (ix = 0; ix < db->getNX(0); ix++)
    for (iy = 0; iy < db->getNX(1); iy++)
    {
      flag_col = 1;

      /* Calculate the experimental covariance of total variable */

      st_calculate_covtot(db, ix, iy, flag_sym, cov_ss, cov_nn, num_tot,
                          cov_tot);
      if (dbg_ix == ix && dbg_iy == iy)
        st_vario_dump(fildmp, ix, iy, cov_ss, cov_nn, num_tot, cov_tot);

      for (iz = 0; iz < db->getNX(2); iz++, ecr++)
      {
        mes_process("3-D Factorial Kriging Analysis", DBOUT->getSampleNumber(),
                    ecr);
        indg[0] = ix;
        indg[1] = iy;
        indg[2] = iz;
        IECH_OUT = db_index_grid_to_sample(db, indg);
        OptDbg::setIndex(IECH_OUT + 1);

        /* Initialize the result to TEST */

        DBOUT->setArray(IECH_OUT, IPTR_EST, TEST);

        if (FFFF(db->getVariable(IECH_OUT, 0)) || !db->isActive(IECH_OUT))
          continue;
        if (OptDbg::query(EDbg::KRIGING) || OptDbg::query(EDbg::NBGH)
            || OptDbg::query(EDbg::RESULTS))
        {
          mestitle(1, "Target location");
          db_sample_print(db, IECH_OUT, 1, 0, 0);
        }

        /* Look for the neighborhood */

        nbgh_ranks = st_neigh_find(db, ix, iy, iz, nei_ss, nei_nn, nei_cur);
        nech = (int) nbgh_ranks.size();
        if (nech <= 0) continue;
        neq = (nfeq == 0) ? nech : nech + 1;

        /* Check if the neighborhood has changed */

        flag_new = flag_col || st_neigh_diff(nei_ss, nei_nn, nei_ref, nei_cur);

        /* If the neighborhood has changed, establish the kriging system */

        flag_col = 0;
        if (flag_new || OptDbg::force())
        {

          /* Establish the L.H.S. of the kriging system */

          st_lhs_exp_3D(nech, nfeq, nei_ss, nei_nn, cov_ss, cov_nn, nei_cur,
                        cov_tot, nugget);
          if (OptDbg::query(EDbg::KRIGING))
            krige_lhs_print(nech, neq, neq, flag, lhs);

          /* Invert the kriging system */

          if (matrix_invert(lhs, neq, IECH_OUT))
          {
            status = 1;
            continue;
          }

          /* Establish the R.H.S. of the kriging system */

          st_rhs_exp_3D(nech, nfeq, nei_ss, nei_nn, cov_ss, cov_nn, nei_cur, cov_res);
          if (OptDbg::query(EDbg::KRIGING))
            krige_rhs_print(nvarin, nech, neq, neq, flag, rhs);

          /* Derive the kriging weights */

          matrix_product(neq, neq, 1, lhs, rhs, wgt);
        }

        /* Calculate the estimation */

        result = st_estim_exp_3D(db, nei_ss, nei_nn, nei_cur, wgt);
        DBOUT->setArray(IECH_OUT, IPTR_EST, result);
        if (OptDbg::query(EDbg::RESULTS)) st_result_kriging_print(0, nvarin, status);
      }
    }

  /* Set the error return flag */

  error = 0;
  if (fildmp != nullptr) fclose(fildmp);

  label_end: OptDbg::setIndex(0);
  (void) krige_koption_manage(-1, 1, EKrigOpt::PONCTUAL, 1, VectorInt());
  st_krige_manage_basic(-1, size_nei, size_nei, 1, nfeq);
  num_tot = (int*) mem_free((char* ) num_tot);
  nei_cur = (int*) mem_free((char* ) nei_cur);
  nei_ref = (int*) mem_free((char* ) nei_ref);
  cov_tot = (double*) mem_free((char* ) cov_tot);
  cov_res = (double*) mem_free((char* ) cov_res);
  return (error);
}

/****************************************************************************/
/*!
 **  Simulate the drift coefficients from the posterior distributions
 **
 ** \return  Error returned code
 **
 ** \param[in] model      Model structure
 ** \param[in] nbsimu     Number of simulation (0 for kriging)
 ** \param[in] rmean      Array giving the posterior means for the drift terms
 ** \param[in] rcov       Array containing the posterior covariance matrix
 **                       for the drift terms
 **
 ** \param[out] smean     Array for simulated posterior mean for the drift means
 **
 ** \remark The input array rcov is modified by this routine, due to
 ** \remark the use of the routine matrix_cholesky_decompose
 **
 *****************************************************************************/
int bayes_simulate(Model *model,
                   int nbsimu,
                   double *rmean,
                   double *rcov,
                   double *smean)
{
  double *trimat, *rndmat;
  int nfeq, il, isimu, nftri, error, rank, memo;

  /* Initializations */

  error = 1;
  nfeq = model->getDriftEquationNumber();
  nftri = nfeq * (nfeq + 1) / 2;
  trimat = rndmat = nullptr;
  memo = law_get_random_seed();

  /* Core allocation */

  trimat = (double*) mem_alloc(sizeof(double) * nftri, 0);
  if (trimat == nullptr) goto label_end;
  rndmat = (double*) mem_alloc(sizeof(double) * nfeq, 0);
  if (rndmat == nullptr) goto label_end;

  /* Cholesky decomposition */

  rank = matrix_cholesky_decompose(rcov, trimat, nfeq);
  if (rank > 0)
  {
    messerr("Error in the Cholesky Decomposition of the covariance matrix");
    messerr("Rank of the Matrix = %d", rank);
    messerr("The Drift coefficients have been set to their posterior mean");
    for (isimu = 0; isimu < nbsimu; isimu++)
      for (il = 0; il < nfeq; il++)
        SMEAN(il,isimu) = rmean[il];
    goto label_suite;
  }

  /* Loop on the simulations */

  for (isimu = 0; isimu < nbsimu; isimu++)
  {

    /* Draw a vector of gaussian independent values */

    for (il = 0; il < nfeq; il++)
      rndmat[il] = law_gaussian();

    /* Product of the Lower triangular matrix by the random vector */

    matrix_cholesky_product(1, nfeq, 1, trimat, rndmat, &SMEAN(0, isimu));

    /* Add the mean */

    for (il = 0; il < nfeq; il++)
      SMEAN(il,isimu) += rmean[il];
  }

  /* If DEBUG option is switched ON, the values are printed out */

  label_suite: if (OptDbg::query(EDbg::BAYES))
  {
    mestitle(1, "Simulation of Drift Coefficients (for Bayesian Simulation)");
    message("Rank     Drift Coefficients\n");
    for (isimu = 0; isimu < nbsimu; isimu++)
    {
      message(" %3d ", isimu + 1);
      for (il = 0; il < nfeq; il++)
        message(" %lf", SMEAN(il, isimu));
      message("\n");
    }
  }

  /* Set the returned error code */

  error = 0;

  label_end: trimat = (double*) mem_free((char* ) trimat);
  rndmat = (double*) mem_free((char* ) rndmat);
  law_set_random_seed(memo);
  return (error);
}

/****************************************************************************/
/*!
 **  Smooth a regular grid
 **
 ** \return  Error return code
 **
 ** \param[in]  dbgrid    input and output Db grid structure
 ** \param[in]  neighI    Neigh structure
 ** \param[in]  type      1 for Uniform; 2 for Gaussian
 ** \param[in]  range     Range (used for Gaussian only)
 ** \param[in]  namconv   Naming Convention
 **
 *****************************************************************************/
int image_smoother(DbGrid *dbgrid,
                   NeighImage *neighI,
                   int type,
                   double range,
                   const NamingConvention& namconv)
{
  /* Add the attribute for storing the results */

  int iptr_est = dbgrid->addColumnsByConstant(1, 0.);
  if (iptr_est < 0) return 1;

  /* Setting options */
  // Here ALL options are set, even if most of them could keep their default values

  KrigingSystem ksys(dbgrid, dbgrid, nullptr, neighI);
  if (ksys.setKrigOptEstim(1, -1, -1)) return 1;
  if (ksys.setKrigOptImageSmooth(true, type, range)) return 1;
  if (! ksys.isReady()) return 1;

  /* Loop on the targets to be processed */

  for (int iech_out = 0; iech_out < dbgrid->getSampleNumber(); iech_out++)
  {
    mes_process("Image Smoothing", dbgrid->getSampleNumber(), iech_out);
     if (ksys.estimate(iech_out)) return 1;
  }

  /* Set the error return flag */

  namconv.setNamesAndLocators(dbgrid, ELoc::Z, 1, dbgrid, iptr_est, "estim");

  return 0;
}

/****************************************************************************/
/*!
 **  Punctual Multivariate Kriging under a constraint
 **
 ** \return  Error return code
 **
 ** \param[in]  dbin       input Db structure
 ** \param[in]  dbout      output Db structure
 ** \param[in]  model      Model structure (univariate)
 ** \param[in]  neighU     NeighUnique structure
 ** \param[in]  flag_positive  1 for a positive constraints
 **
 ** \remark  All the variables are estimated using the same model
 ** \remark  In this procedure, we assume that:
 ** \remark  - the problem is multivariate ("z" variables)
 ** \remark  - the constraints is stored in "f" (only used in dbout)
 **
 *****************************************************************************/
int krigsum_f(Db *dbin,
              Db *dbout,
              Model *model,
              NeighUnique *neighU,
              int flag_positive)
{
  double *lterm, seisloc, seistot, estim;
  int *icols, *active, error, iptr_mem, correct;
  int nvarmod, nvarin, status, nech, ivar, nred, neq;
  NeighWork nbghw;
  VectorInt nbgh_ranks;

  /* Preliminary checks */

  error = 1;
  neq = nred = 0;
  st_global_init(dbin, dbout);
  icols = active = nullptr;
  lterm = nullptr;
  nvarin = dbin->getVariableNumber();
  nvarmod = model->getVariableNumber();
  FLAG_EST = 1;
  FLAG_LTERM = 1;

  if (nvarmod != 1)
  {
    messerr("This procedure requires a monovariate model");
    goto label_end;
  }
  if (model_nfex(model) != 0)
  {
    messerr("This procedure requires a model with no External Drift");
    goto label_end;
  }
  if (dbout->getExternalDriftNumber() != 1)
  {
    messerr("This procedure requires one External Drift in the Output Db");
    messerr("The number of External Drift is currently equal to %d",
            dbout->getExternalDriftNumber());
    goto label_end;
  }

  /* Core allocation */

  icols = (int*) mem_alloc(sizeof(int) * nvarin, 0);
  if (icols == nullptr) goto label_end;
  active = (int*) mem_alloc(sizeof(int) * nvarin, 0);
  if (active == nullptr) goto label_end;
  lterm = (double*) mem_alloc(sizeof(double) * nvarin, 0);
  if (lterm == nullptr) goto label_end;

  /* Save the columns for variable definitions */

  for (ivar = 0; ivar < nvarin; ivar++)
    icols[ivar] = db_attribute_identify(dbin, ELoc::Z, ivar);
  dbin->clearLocators(ELoc::Z);
  dbin->setLocatorByUID(icols[0], ELoc::Z);
  if (st_check_environment(1, 1, model, neighU)) goto label_end;

  /* Add the attributes for storing the results */

  iptr_mem = dbout->addColumnsByConstant(nvarin, 0);
  if (iptr_mem < 0) goto label_end;

  /* Pre-calculations */

  nbghw.initialize(DBIN, neighU);
  nbghw.setFlagSimu(FLAG_SIMU);
  if (st_model_manage(1, model)) goto label_end;
  if (st_krige_manage(1, nvarin, model, neighU)) goto label_end;
  if (krige_koption_manage(1, 1, EKrigOpt::PONCTUAL, 1, VectorInt()))
    goto label_end;

  /* Loop on the variables */

  status = 0;
  for (ivar = 0; ivar < nvarin; ivar++)
  {
    dbin->clearLocators(ELoc::Z);
    dbin->setLocatorByUID(icols[ivar], ELoc::Z);
    IPTR_EST = iptr_mem + ivar;
    (void) gslSPrintf(string, "Kriging of variable #%d at sample", ivar + 1);

    /* Loop on the targets to be processed */

    for (IECH_OUT = 0; IECH_OUT < DBOUT->getSampleNumber(); IECH_OUT++)
    {
      mes_process(string, DBOUT->getSampleNumber(), IECH_OUT);
      OptDbg::setIndex(IECH_OUT + 1);
      if (!dbout->isActive(IECH_OUT)) continue;
      if (OptDbg::query(EDbg::KRIGING) || OptDbg::query(EDbg::NBGH)
          || OptDbg::query(EDbg::RESULTS))
      {
        mestitle(1, "Target location");
        db_sample_print(dbout, IECH_OUT, 1, 0, 0);
      }

      /* Select the Neighborhood */

      nbgh_ranks = nbghw.select(DBOUT,  IECH_OUT);
      nech = (int) nbgh_ranks.size();
      status = nbgh_ranks.empty();
      if (status) goto label_store;

      /* Establish the kriging L.H.S. */

      if (! nbghw.isUnchanged() || neighU->getFlagContinuous() || OptDbg::force())
      {
        st_prepar(model, neighU, nbgh_ranks, &status, &nred, &neq);
        if (status) goto label_store;
        st_data_dual(model, NULL, nbgh_ranks, nred, &lterm[ivar]);
      }

      /* Establish the kriging R.H.S. */

      st_rhs(model, nbgh_ranks, neq, nvarmod, VectorVectorDouble(), &status);
      if (status) goto label_store;
      st_rhs_iso2hetero(neq, nvarmod);
      if (OptDbg::query(EDbg::KRIGING))
        krige_rhs_print(nvarmod, nech, neq, nred, flag, rhs);

      /* Perform the estimation */

      label_store: st_estimate(model, NULL, status, 0, nvarmod, nred);
      if (OptDbg::query(EDbg::RESULTS))
        st_result_kriging_print(neighU->getFlagXvalid(), nvarmod, status);
    }
  }

  /* Posterior scaling */

  for (IECH_OUT = 0; IECH_OUT < DBOUT->getSampleNumber(); IECH_OUT++)
  {
    correct = 0;
    for (ivar = 0; ivar < nvarin; ivar++)
      active[ivar] = 0;

    /* Implicit loop until the solution is acceptable */

    while (!correct)
    {
      seistot = 0.;
      seisloc = dbout->getExternalDrift(IECH_OUT, 0);
      for (ivar = 0; ivar < nvarin; ivar++)
      {
        if (active[ivar]) continue;
        seistot += lterm[ivar];
        seisloc -= dbout->getArray(IECH_OUT, iptr_mem + ivar);
      }
      if (seistot == 0.)
      {
        messerr("The sum of scaling terms is zero. No correction is possible");
        goto label_end;
      }

      for (ivar = 0; ivar < nvarin; ivar++)
      {
        if (active[ivar])
          estim = 0;
        else
          estim = (dbout->getArray(IECH_OUT, iptr_mem + ivar)
              + lterm[ivar] * seisloc / seistot);
        dbout->setArray(IECH_OUT, iptr_mem + ivar, estim);
      }

      correct = 1;
      for (ivar = 0; ivar < nvarin; ivar++)
      {
        active[ivar] = (dbout->getArray(IECH_OUT, iptr_mem + ivar) < 0);
        if (active[ivar]) correct = 0;
      }
      if (!flag_positive) correct = 1;
    }
  }

  /* Set the error return flag */

  error = 0;

  label_end: OptDbg::setIndex(0);
  (void) st_model_manage(-1, model);
  (void) st_krige_manage(-1, nvarin, model, neighU);
  (void) krige_koption_manage(-1, 1, EKrigOpt::PONCTUAL, 1, VectorInt());
  icols = (int*) mem_free((char* ) icols);
  active = (int*) mem_free((char* ) active);
  lterm = (double*) mem_free((char* ) lterm);
  return (error);
}

/****************************************************************************/
/*!
 **  Check that the proportions are positive
 **
 ** \return  Number of erroneous proportions
 **
 ** \param[in]  db3grid   Db 3-D Grid
 ** \param[in]  ix        Index of the grid node along X
 ** \param[in]  iy        Index of the grid node along Y
 ** \param[in]  nvarin    Total number of facies
 ** \param[in]  iptr_prop Rank of the proportion attribute
 ** \param[in]  proptab   Array of proportions (dimension: nfac * nx * ny * nz)
 **
 ** \param[out]  lterm   The resulting term Z * C-1 * Z
 **
 ** \remark  If a proportion is negative in a cell, the proportion is set to 0
 ** \remark  and the corresponding LTERM is set to 0
 **
 *****************************************************************************/
static int st_check_positivity(DbGrid *db3grid,
                               int ix,
                               int iy,
                               int nvarin,
                               int iptr_prop,
                               double *proptab,
                               double *lterm)
{
  int indg[3], n_wrong, iech, iz, ivar;
  double prop;

  indg[0] = ix;
  indg[1] = iy;
  n_wrong = 0;
  for (iz = 0; iz < db3grid->getNX(2); iz++)
  {
    indg[2] = iz;
    iech = db_index_grid_to_sample(db3grid, indg);
    if (!db3grid->isActive(iech)) continue;

    /* No estimated proportion at this level */

    for (ivar = 0; ivar < nvarin; ivar++)
    {
      prop = PROPTAB(ivar, iz);
      if (prop >= 0) continue;
      db3grid->setArray(iech, iptr_prop + ivar, 0.);
      LTERM(ivar,iz)= 0.;
      n_wrong++;
    }
  }
  return (n_wrong);
}

/****************************************************************************/
/*!
 **  Check the seismic constraint
 **
 ** \return  Error return code
 **
 ** \param[in]  ix       Index of the grid node along X
 ** \param[in]  iy       Index of the grid node along Y
 ** \param[in]  nz       Number of layers along Z
 ** \param[in]  nvarin   Number of facies
 ** \param[in]  fsum     Rank of the proportion facies which average up
 **                      to the seismic (no constraint if negative)
 ** \param[in]  seisval  Seismic value
 ** \param[in]  proptab  Array of proportions (dimension: nfac * nx * ny * nz)
 **
 *****************************************************************************/
static int st_check_constraint_seismic(int ix,
                                       int iy,
                                       int nz,
                                       int nvarin,
                                       int fsum,
                                       double seisval,
                                       double *proptab)
{
  double prop;
  int iz, error;

  error = 0;
  if (FFFF(seisval)) return (0);

  /* Calculate the average proportion corresponding to the seismic */

  prop = 0.;
  for (iz = 0; iz < nz; iz++)
    prop += PROPTAB(fsum, iz);
  prop /= nz;

  /* Compare seismic and resulting average propotion */

  if (ABS(prop - seisval) > EPSILON5)
  {
    messerr(
        "Block (%d,%d,%d) - Mismatch between proportion (%lf) and seismic (%lf)",
        ix + 1, iy + 1, iz + 1, prop, seisval);
    error++;
  }
  return (error);
}

/****************************************************************************/
/*!
 **  Simultaneous Kriging of proportions
 **  - the proportions must sum up to 1 (per cell)
 **  - the average of the proportions of a given facies is equal to
 **    the seismic value read in the "f" variable of the db2grid
 **
 ** \return  Error return code
 **
 ** \param[in]  dbin       input Db structure
 ** \param[in]  db3grid    output 3-D Grid Db structure
 ** \param[in]  db2grid    output 2-D Grid Db structure
 ** \param[in]  fsum       Rank of the proportion facies which average up
 **                        to the seismic (no constraint if negative)
 ** \param[in]  model      Model structure (monovariate)
 ** \param[in]  neighparam ANeighParam structure (Unique or Bench)
 **
 *****************************************************************************/
int krigmvp_f(Db *dbin,
              DbGrid *db3grid,
              DbGrid *db2grid,
              int fsum,
              Model *model,
              ANeighParam *neighparam)
{
  int *icols, indg[3], nloc, error;
  int status, nech, neq, nred, nvarin, nvarmod, nfeq, nz, pivot;
  int ix, iy, iz, ivar, i, iech, jech, iptr_prop, n_wrong, iter, nsize;
  int flag_correc;
  double *lterm, *lback, *proptab, *cc, *xx, *bb;
  double seisval, correc, proploc, lsum;
  NeighWork nbghw;
  VectorInt nbgh_ranks;

  /* Preliminary checks */

  error = 1;
  neq = nred = 0;
  st_global_init(dbin, db3grid);
  icols = nullptr;
  lterm = lback = proptab = cc = xx = bb = nullptr;
  nvarin = dbin->getVariableNumber();
  nvarmod = model->getVariableNumber();
  nfeq = model->getDriftEquationNumber();
  iptr_prop = iter = nsize = 0;
  seisval = TEST;
  FLAG_EST = 1;
  FLAG_WGT = 1;
  FLAG_LTERM = 1;

  if (nvarmod != 1)
  {
    messerr("This procedure requires a monovariate model");
    goto label_end;
  }
  if (model_nfex(model) != 0)
  {
    messerr("This procedure requires a model with no External Drift");
    goto label_end;
  }
  if (neighparam->getType() != ENeigh::UNIQUE && neighparam->getType() != ENeigh::BENCH)
  {
    messerr(
        "This procedure is currently limited to the Unique or Bench Neighborhood");
    goto label_end;
  }
  if (!is_grid(db3grid))
  {
    messerr("This procedure needs a Grid 3-D Db structure");
    goto label_end;
  }
  if (fsum >= 0)
  {
    if (!is_grid(db2grid))
    {
      messerr("This procedure needs a Grid 2-D Db structure");
      goto label_end;
    }
    if (!db_grid_match(db2grid, db3grid)) goto label_end;
    if (db2grid->getExternalDriftNumber() != 1)
    {
      messerr("This procedure requires one External Drift in the 2-D Grid Db");
      messerr("The number of External Drift is currently equal to %d",
              db2grid->getExternalDriftNumber());
      goto label_end;
    }
    if (fsum >= nvarin)
    {
      messerr(
          "Error in the rank of the variable whose average vertical proportion");
      messerr("Should match the seismic of the 2-D Db grid (%d)", fsum + 1);
      messerr("It should lie between 1 and %d", nvarin);
      goto label_end;
    }
  }
  nz = db3grid->getNX(2);

  /* Core allocation */

  lback = (double*) mem_alloc(sizeof(double) * nvarin * nz, 0);
  if (lback == nullptr) goto label_end;
  lterm = (double*) mem_alloc(sizeof(double) * nvarin * nz, 0);
  if (lterm == nullptr) goto label_end;
  for (i = 0; i < nvarin * nz; i++)
    lterm[i] = lback[i] = TEST;
  icols = (int*) mem_alloc(sizeof(int) * nvarin, 0);
  if (icols == nullptr) goto label_end;
  if (fsum >= 0)
  {
    proptab = (double*) mem_alloc(sizeof(double) * nz * nvarin, 0);
    if (proptab == nullptr) goto label_end;
  }

  /* Save the columns for variable definitions */

  for (ivar = 0; ivar < nvarin; ivar++)
    icols[ivar] = db_attribute_identify(dbin, ELoc::Z, ivar);
  dbin->clearLocators(ELoc::Z);
  dbin->setLocatorByUID(icols[0], ELoc::Z);
  if (st_check_environment(1, 1, model, neighparam)) goto label_end;

  /* Add the attributes for storing the results */

  iptr_prop = db3grid->addColumnsByConstant(nvarin, 0);
  if (iptr_prop < 0) goto label_end;

  /* Pre-calculations */

  nbghw.initialize(DBIN, neighparam);
  nbghw.setFlagSimu(FLAG_SIMU);
  if (st_model_manage(1, model)) goto label_end;
  if (st_krige_manage(1, nvarmod, model, neighparam)) goto label_end;
  if (krige_koption_manage(1, 1, EKrigOpt::PONCTUAL, 1, VectorInt()))
    goto label_end;

  /* Loop on the target grid nodes to be processed */

  status = 0;
  for (ivar = 0; ivar < nvarin; ivar++)
  {
    dbin->clearLocators(ELoc::Z);
    dbin->setLocatorByUID(icols[ivar], ELoc::Z);
    IPTR_EST = iptr_prop + ivar;
    (void) gslSPrintf(string, "Kriging of proportion #%d at sample", ivar + 1);

    /* Loop on the target grid nodes */

    for (iz = IECH_OUT = 0; iz < db3grid->getNX(2); iz++)
    {
      for (iy = 0; iy < db3grid->getNX(1); iy++)
        for (ix = 0; ix < db3grid->getNX(0); ix++, IECH_OUT++)
        {
          mes_process(string, DBOUT->getSampleNumber(), IECH_OUT);
          OptDbg::setIndex(IECH_OUT + 1);
          if (!db3grid->isActive(IECH_OUT)) continue;
          if (OptDbg::query(EDbg::KRIGING) || OptDbg::query(EDbg::NBGH)
              || OptDbg::query(EDbg::RESULTS))
          {
            mestitle(1, "Target location");
            db_sample_print(db3grid, IECH_OUT, 1, 0, 0);
          }

          /* Select the Neighborhood */

          nbgh_ranks = nbghw.select(DBOUT,  IECH_OUT);
          nech = (int) nbgh_ranks.size();
          status = nbgh_ranks.empty();
          if (status) goto label_store;

          /* Establish the kriging L.H.S. */

          if (! nbghw.isUnchanged() || neighparam->getFlagContinuous() || OptDbg::force())
          {
            st_prepar(model, neighparam, nbgh_ranks, &status, &nred, &neq);
            if (status) goto label_store;
            st_data_dual(model, NULL, nbgh_ranks, nred, &LBACK(ivar, iz));
          }

          /* Establish the kriging R.H.S. */

          st_rhs(model, nbgh_ranks, neq, nvarmod, VectorVectorDouble(), &status);
          if (status) goto label_store;
          st_rhs_iso2hetero(neq, nvarmod);
          if (OptDbg::query(EDbg::KRIGING))
            krige_rhs_print(nvarmod, nech, neq, nred, flag, rhs);

          /* Derive the kriging weights */

          if (FLAG_WGT)
          {
            matrix_product(nred, nred, nvarmod, lhs, rhs, wgt);
            if (OptDbg::query(EDbg::KRIGING))
              krige_wgt_print(status, nvarmod, nvarmod, nfeq, nbgh_ranks, nred, -1,
                              flag, wgt);
          }

          /* Perform the estimation */

          label_store: st_estimate(model, NULL, status, 0, nvarmod, nred);
          if (OptDbg::query(EDbg::RESULTS))
            st_result_kriging_print(neighparam->getFlagXvalid(), nvarmod, status);
        }
    }
  }

  /* Discard impossible constraints */

  lsum = 0.;
  for (ivar = 0; ivar < nvarin; ivar++)
    for (iz = 0; iz < nz; iz++)
      lsum += LBACK(ivar, iz);
  if (lsum <= EPSILON5)
  {
    error = 0;
    goto label_end;
  }

  /* Core allocation */

  nsize = (fsum >= 0) ? nz + 1 :
                        nz;
  cc = (double*) mem_alloc(sizeof(double) * nsize * (nsize + 1) / 2, 0);
  if (cc == nullptr) goto label_end;
  xx = (double*) mem_alloc(sizeof(double) * nsize, 0);
  if (xx == nullptr) goto label_end;
  bb = (double*) mem_alloc(sizeof(double) * nsize, 0);
  if (bb == nullptr) goto label_end;

  /* Posterior corrections */

  for (iy = 0; iy < db3grid->getNX(1); iy++)
    for (ix = 0; ix < db3grid->getNX(0); ix++)
    {
      iter = 0;
      for (i = 0; i < nvarin * nz; i++)
        lterm[i] = lback[i];

      /* Loop on the iterations for the same block */

      label_suite: iter++;

      indg[0] = ix;
      indg[1] = iy;
      nloc = nz;
      if (fsum >= 0)
      {
        indg[2] = 0;
        jech = db_index_grid_to_sample(db2grid, indg);
        seisval = db2grid->getExternalDrift(jech, 0);
        if (!FFFF(seisval)) nloc++;
      }

      /* Establish the system */

      for (i = 0; i < nloc * nloc; i++)
        cc[i] = 0.;
      for (iz = 0; iz < db3grid->getNX(2); iz++)
      {
        bb[iz] = -1;
        indg[2] = iz;
        iech = db_index_grid_to_sample(db3grid, indg);
        for (ivar = 0; ivar < nvarin; ivar++)
        {
          bb[iz] += db3grid->getArray(iech, iptr_prop + ivar);
          CC(iz,iz) += LTERM(ivar, iz);
        }
      }

      if (nloc > nz)
      {
        bb[nz] = -2 * seisval;
        for (iz = 0; iz < db3grid->getNX(2); iz++)
        {
          indg[2] = iz;
          iech = db_index_grid_to_sample(db3grid, indg);
          bb[nz] += db3grid->getArray(iech, iptr_prop + fsum);
          CC(iz,nz) += LTERM(fsum, iz);
          CC(nz,nz) += LTERM(fsum, iz);
        }
      }

      /* Check if correction is needed */

      flag_correc = 0;
      for (i = 0; i < nloc && flag_correc == 0; i++)
        if (ABS(bb[i]) > EPSILON5) flag_correc = 1;
      if (!flag_correc) continue;

      /* Solve the system */

      if (matrix_solve(0, cc, bb, xx, nloc, 1, &pivot))
        messageAbort("Core problem in matrix_solve");
      if (pivot > 0)
      {
        messerr("Error during the inversion of the constraint matrix");
        messerr("The constraints may be redundant");
        goto label_end;
      }

      /* Perform the correction */

      for (iz = 0; iz < db3grid->getNX(2); iz++)
      {
        indg[2] = iz;
        iech = db_index_grid_to_sample(db3grid, indg);
        for (ivar = 0; ivar < nvarin; ivar++)
        {
          proploc = db3grid->getArray(iech, iptr_prop + ivar);
          correc = xx[iz];
          if (nloc > nz && ivar == fsum) correc += xx[nz];
          PROPTAB(ivar,iz) = proploc - correc * LTERM(ivar, iz);
        }
      }

      /* Check the seismic constraints */

      if (fsum >= 0)
        (void) st_check_constraint_seismic(ix, iy, nz, nvarin, fsum, seisval,
                                           proptab);

      /* Check positivity of the proportions */

      n_wrong = st_check_positivity(db3grid, ix, iy, nvarin, iptr_prop, proptab,
                                    lterm);
      if (n_wrong > 0) goto label_suite;

      /* Final update the proportions */

      for (iz = 0; iz < nz; iz++)
        for (ivar = 0; ivar < nvarin; ivar++)
        {
          indg[2] = iz;
          iech = db_index_grid_to_sample(db3grid, indg);
          db3grid->setArray(iech, iptr_prop + ivar, PROPTAB(ivar, iz));
        }
    }

  /* Set the error return flag */

  error = 0;

  label_end: OptDbg::setIndex(0);
  (void) st_model_manage(-1, model);
  (void) st_krige_manage(-1, nvarmod, model, neighparam);
  (void) krige_koption_manage(-1, 1, EKrigOpt::PONCTUAL, 1, VectorInt());
  icols = (int*) mem_free((char* ) icols);
  lback = (double*) mem_free((char* ) lback);
  lterm = (double*) mem_free((char* ) lterm);
  cc = (double*) mem_free((char* ) cc);
  xx = (double*) mem_free((char* ) xx);
  bb = (double*) mem_free((char* ) bb);
  proptab = (double*) mem_free((char* ) proptab);
  return (error);
}

/****************************************************************************/
/*!
 **  Perform kriging and return the calculation elements
 **
 ** \return  A Krigtest_Res structure
 **
 ** \param[in]  dbin       input Db structure
 ** \param[in]  dbout      output Db structure
 ** \param[in]  model      Model structure
 ** \param[in]  neighparam ANeighParam structure
 ** \param[in]  iech0      Rank of the target sample
 ** \param[in]  calcul     Kriging calculation option (EKrigOpt)
 ** \param[in]  ndisc      Array giving the discretization counts
 **
 *****************************************************************************/
Krigtest_Res krigtest(Db *dbin,
                      Db *dbout,
                      Model *model,
                      ANeighParam *neighparam,
                      int iech0,
                      const EKrigOpt &calcul,
                      VectorInt ndisc)
{
  Krigtest_Res ktest;

  // Preliminary checks

  if (neighparam->getType() == ENeigh::IMAGE)
  {
    messerr("This tool cannot function with an IMAGE neighborhood");
    return ktest;
  }

  // Initializations

  int iptr_est  = -1;
  int iptr_std  = -1;
  int iptr_varz = -1;
  int nvar = model->getVariableNumber();

  /* Add the attributes for storing the results */

  iptr_est = dbout->addColumnsByConstant(nvar, 0.);
  if (iptr_est < 0) return ktest;
  iptr_std = dbout->addColumnsByConstant(nvar, 0.);
  if (iptr_std < 0) return ktest;

  /* Setting options */

  KrigingSystem ksys(dbin, dbout, model, neighparam);
  if (ksys.setKrigOptEstim(iptr_est, iptr_std, iptr_varz)) return ktest;
  if (ksys.setKrigOptCalcul(calcul, ndisc)) return ktest;
  if (! ksys.isReady()) return ktest;

  /* Loop on the targets to be processed */

  if (ksys.estimate(iech0)) return ktest;

  /* Extract relevant information */

  ktest.ndim = ksys.getNDim();
  ktest.nech = ksys.getNRed();
  ktest.nrhs = 1;
  ktest.neq  = ksys.getNeq();
  ktest.nbgh = ksys.getSampleIndices();
  ktest.xyz  = ksys.getSampleCoordinates();
  ktest.data = ksys.getSampleData();
  ktest.zam  = ksys.getZam();
  ktest.lhs  = ksys.getLHS();
  ktest.rhs  = ksys.getRHSC();
  ktest.wgt  = ksys.getWeights();
  ktest.var  = ksys.getVariance();

  /* Delete fields added in Dbout during calculations */

  for (int ivar = 0; ivar < nvar; ivar++)
    dbout->deleteColumnByUID(iptr_est + ivar);
  for (int ivar = 0; ivar < nvar; ivar++)
    dbout->deleteColumnByUID(iptr_std + ivar);

  return ktest;
}

/****************************************************************************/
/*!
 **  Transform the Kriging results from gaussian to raw
 **
 ** \param[in]  anam      AAnam structure
 **
 ** \remark  This procedure is designed for the monovariate case
 ** \remark  It assumes that the kriging estimate and variance are already
 ** \remark  calculated
 **
 *****************************************************************************/
static void st_transform_gaussian_to_raw(AAnam *anam)
{
  if (anam == nullptr) return;
  AnamHermite *anam_hermite = dynamic_cast<AnamHermite*>(anam);

  /* Get the estimation */

  double est = DBOUT->getArray(IECH_OUT, IPTR_EST);

  /* Get the variance of the kriging error */

  double std = sqrt(DBOUT->getArray(IECH_OUT, IPTR_STD));

  /* Calculate the conditional expectation */

  double condexp = hermiteCondExpElement(est, std, anam_hermite->getPsiHn());
  DBOUT->setArray(IECH_OUT, IPTR_EST, condexp);

  /* Calculate the conditional variance */

  double condvar = hermiteEvaluateZ2(est, std, anam_hermite->getPsiHn());
  condvar -= condexp * condexp;
  DBOUT->setArray(IECH_OUT, IPTR_STD, condvar);
}

/****************************************************************************/
/*!
 **  Punctual Kriging in the Anamorphosed Gaussian Model
 **
 ** \return  Error return code
 **
 ** \param[in]  dbin       input Db structure
 ** \param[in]  dbout      output Db structure
 ** \param[in]  anam       AAnam structure
 ** \param[in]  model      Model structure
 ** \param[in]  neighparam ANeighParam structure
 **
 *****************************************************************************/
int kriggam_f(Db *dbin,
              Db *dbout,
              AAnam *anam,
              Model *model,
              ANeighParam *neighparam)
{
  int error, status, nech, neq, nred, nvar, nfeq;
  double ldum;
  NeighWork nbghw;
  VectorInt nbgh_ranks;

  /* Preliminary checks */

  error = 1;
  nvar = nred = neq = 0;
  st_global_init(dbin, dbout);
  FLAG_EST = 1;
  FLAG_STD = 1;
  if (st_check_environment(1, 1, model, neighparam)) goto label_end;
  nvar = model->getVariableNumber();
  nfeq = model->getDriftEquationNumber();

  /* Preliminary check */

  if (nvar != 1)
  {
    messerr("This procedure is limited to the monovariate case");
    goto label_end;
  }

  /* Add the attributes for storing the results */

  if (FLAG_EST != 0)
  {
    IPTR_EST = dbout->addColumnsByConstant(nvar, 0.);
    if (IPTR_EST < 0) goto label_end;
  }
  if (FLAG_STD != 0)
  {
    IPTR_STD = dbout->addColumnsByConstant(nvar, 0.);
    if (IPTR_STD < 0) goto label_end;
  }

  /* Pre-calculations */

  nbghw.initialize(DBIN, neighparam);
  nbghw.setFlagSimu(FLAG_SIMU);
  if (st_model_manage(1, model)) goto label_end;
  if (st_krige_manage(1, nvar, model, neighparam)) goto label_end;
  if (krige_koption_manage(1, 1, EKrigOpt::PONCTUAL, 1, VectorInt()))
    goto label_end;
  if (FLAG_STD != 0) st_variance0(model, nvar, VectorVectorDouble());

  /* Loop on the targets to be processed */

  status = 0;
  for (IECH_OUT = 0; IECH_OUT < DBOUT->getSampleNumber(); IECH_OUT++)
  {
    mes_process("Kriging sample", DBOUT->getSampleNumber(), IECH_OUT);
    OptDbg::setIndex(IECH_OUT + 1);
    if (!dbout->isActive(IECH_OUT)) continue;
    if (OptDbg::query(EDbg::KRIGING) || OptDbg::query(EDbg::NBGH) || OptDbg::query(EDbg::RESULTS))
    {
      mestitle(1, "Target location");
      db_sample_print(dbout, IECH_OUT, 1, 0, 0);
    }

    /* Select the Neighborhood */

    nbgh_ranks = nbghw.select(DBOUT,  IECH_OUT);
    nech = (int) nbgh_ranks.size();
    status = nbgh_ranks.empty();
    if (status) goto label_store;

    /* Establish the kriging L.H.S. */

    if (! nbghw.isUnchanged() || neighparam->getFlagContinuous() || OptDbg::force())
    {
      st_prepar(model, neighparam, nbgh_ranks, &status, &nred, &neq);
      if (status) goto label_store;
      st_data_dual(model, NULL, nbgh_ranks, nred, &ldum);
    }

    /* Establish the kriging R.H.S. */

    st_rhs(model, nbgh_ranks, neq, nvar, VectorVectorDouble(), &status);
    if (status) goto label_store;
    st_rhs_iso2hetero(neq, nvar);
    if (OptDbg::query(EDbg::KRIGING))
      krige_rhs_print(nvar, nech, neq, nred, flag, rhs);

    /* Derive the kriging weights */

    if (FLAG_WGT)
    {
      matrix_product(nred, nred, nvar, lhs, rhs, wgt);
      if (OptDbg::query(EDbg::KRIGING))
        krige_wgt_print(status, nvar, nvar, nfeq, nbgh_ranks, nred, -1, flag, wgt);
    }

    /* Perform the estimation */

    label_store:
    st_estimate(model, NULL, status, neighparam->getFlagXvalid(), nvar,nred);

    /* Transform the gaussian estimates into raw estimates */

    st_transform_gaussian_to_raw(anam);

    if (OptDbg::query(EDbg::RESULTS))
      st_result_kriging_print(neighparam->getFlagXvalid(), nvar, status);
  }

  /* Set the error return flag */

  error = 0;

  label_end: OptDbg::setIndex(0);
  (void) st_model_manage(-1, model);
  (void) st_krige_manage(-1, nvar, model, neighparam);
  (void) krige_koption_manage(-1, 1, EKrigOpt::PONCTUAL, 1, VectorInt());
  return (error);
}

/****************************************************************************/
/*!
 **  Standard Block Kriging with variable cell dimension
 **
 ** \return  Error return code
 **
 ** \param[in]  dbin        Input Db structure
 ** \param[in]  dbout       Output Db structure
 ** \param[in]  model       Model structure
 ** \param[in]  neighparam  ANeighParam structure
 ** \param[in]  ndisc       Array giving the discretization counts
 ** \param[in]  flag_est    Option for the storing the estimation
 ** \param[in]  flag_std    Option for the storing the standard deviation
 ** \param[in]  rank_colcok Option for running Collocated Cokriging
 ** \param[in]  namconv     Naming convention
 **
 *****************************************************************************/
int krigcell(Db *dbin,
             Db *dbout,
             Model *model,
             ANeighParam *neighparam,
             int flag_est,
             int flag_std,
             VectorInt ndisc,
             VectorInt rank_colcok,
             const NamingConvention& namconv)
{
  // Preliminary checks

  if (neighparam->getType() == ENeigh::IMAGE)
  {
    messerr("This tool cannot function with an IMAGE neighborhood");
    return 1;
  }

  // Initializations

  int iptr_est  = -1;
  int iptr_std  = -1;
  int nvar = model->getVariableNumber();

  /* Add the attributes for storing the results */

  if (flag_est != 0)
  {
    iptr_est = dbout->addColumnsByConstant(nvar, 0.);
    if (iptr_est < 0) return 1;
  }
  if (flag_std != 0)
  {
    iptr_std = dbout->addColumnsByConstant(nvar, 0.);
    if (iptr_std < 0) return 1;
  }

  /* Setting options */

  KrigingSystem ksys(dbin, dbout, model, neighparam);
  if (ksys.setKrigOptEstim(iptr_est, iptr_std, -1)) return 1;
  if (ksys.setKrigOptCalcul(EKrigOpt::BLOCK, ndisc, true)) return 1;
  if (ksys.setKrigOptColCok(rank_colcok)) return 1;
  if (! ksys.isReady()) return 1;

  /* Loop on the targets to be processed */

  for (int iech_out = 0; iech_out < dbout->getSampleNumber(); iech_out++)
  {
    mes_process("Kriging sample", dbout->getSampleNumber(), iech_out);
    if (ksys.estimate(iech_out)) return 1;
  }

  /* Set the error return flag */

  namconv.setNamesAndLocators(dbin, ELoc::Z, nvar, dbout, iptr_std, "stdev", 1,
                              false);
  namconv.setNamesAndLocators(dbin, ELoc::Z, nvar, dbout, iptr_est, "estim");

  return 0;
}

/****************************************************************************/
/*!
 **  Calculate the Hermite polynomials at the data samples
 **
 ** \return  Error return code
 **
 ** \param[in]  db        Input Db structure (containing the factors)
 ** \param[in]  nfactor   Number of factors to be estimated (0: all)
 **
 ** \remark At the end, the newly created variables are transformed into
 ** \remark Z variables for future steps
 **
 *****************************************************************************/
static int st_calculate_hermite_factors(Db *db, int nfactor)
{
  VectorDouble hn;
  int iptr, error;

  /* Initializations */

  error = 1;

  /* Create the new variables */

  iptr = db->addColumnsByConstant(nfactor, 0.);
  if (iptr < 0) goto label_end;

  /* Loop on the samples */

  for (int iech = 0; iech < db->getSampleNumber(); iech++)
  {
    if (!db->isActive(iech)) continue;

    /* Calculate the factors */

    hn = hermitePolynomials(db->getVariable(iech, 0), 1., nfactor + 1);

    /* Store the factors */

    for (int ih = 0; ih < nfactor; ih++)
      db->setArray(iech, iptr + ih, hn[ih + 1]);
  }

  /* Set the newly created variables to Z locator */

  db->setLocatorsByUID(nfactor, iptr, ELoc::Z);

  /* Set the error return code */

  error = 0;

  label_end: return (error);
}

/****************************************************************************/
/*!
 **  Disjunctive Kriging
 **
 ** \return  Error return code
 **
 ** \param[in]  dbin       input Db structure (containing the factors)
 ** \param[in]  dbgrid     output Grid Db structure
 ** \param[in]  model      Model structure
 ** \param[in]  neighparam ANeighParam structure
 ** \param[in]  nfactor    Number of factors to be estimated (0: all)
 ** \param[in]  nmult      Array of Multiplicity for Partition
 ** \param[in]  ndisc      Discretization parameters (or NULL)
 ** \param[in]  flag_est   Option for the storing the estimation
 ** \param[in]  flag_std   Option for the storing the standard deviation
 **
 ** \remark In case the Model handles a Punctual Anamophossis, the
 ** \remark estimation of block average quantities is performed. This initiates
 ** \remark a block estimation kriging and therefore requires the definition of
 ** \remark the block discretization ('ndisc')
 **
 ** \remark In case the Model handles a Block Anamorphosis, the estimation
 ** \remark of probabilities per block is performed. This uses a simple point
 ** \remark (no use of 'ndisc'). Then the user may wish either to exhibit Q,T
 ** \remark on a SMU basis ('dbgrid'should then be the grid of SMUs) or to
 ** \remark directly produce Q,T on the panels (a panel is partitioned into
 ** \remark a set of SMU: this partition is defined through the argument 'nmult')
 **
 ** TODO : Check if the following remark is up to date!
 ** \remark The value 'KOPTION->calcul' must be set to:
 ** \remark - EKrigOpt::PONCTUAL for point-block estimation (flag_block = TRUE)
 ** \remark - EKrigOpt::BLOCK for point estimation
 ** \remark Nevertheless, for Point-Block model, if dbgrid is a grid of Panels
 ** \remark 'KOPTION->calcul' is set to EKrigOpt::BLOCK to provoke the discretization
 ** \remark of Panel into SMUs.
 **
 *****************************************************************************/
int dk_f(Db *dbin,
         DbGrid *dbgrid,
         Model *model,
         ANeighParam *neighparam,
         int nfactor,
         const VectorInt &nmult,
         const VectorInt &ndisc,
         int flag_est,
         int flag_std)
{
  DbGrid *dbsmu;
  int error, status, nech, neq, nred, nvar, nvarz, nfeq, iclass, imult, nb_mult;
  int iptr_est_bck, iptr_std_bck, ivar, i;
  int *varloc, flag_block, flag_panel, flag_continue, neqmax;
  double *rhs_cum, ldum;
  CovLMCAnamorphosis* covanam;
  static double perturb = 1.e-8;
  NeighWork nbghw;
  VectorInt nbgh_ranks;

  /* Preliminary checks */

  error = 1;
  iptr_est_bck = iptr_std_bck = -1;
  rhs_cum = nullptr;
  varloc = nullptr;
  dbsmu = nullptr;
  st_global_init(dbin, dbgrid);
  FLAG_EST = flag_est;
  FLAG_STD = flag_std;
  FLAG_WGT = flag_std;

  // The model is not checked against the Data, as the number of variables
  // is not consistent: Model (1) whereas Data (nfactor-1)
  if (st_check_environment(1, 1, NULL, neighparam)) goto label_end;

  if (model->getVariableNumber() > 1)
  {
    messerr("This application is limited to the monovariate Model case");
    goto label_end;
  }
  if (model->getCovMode() != EModelProperty::ANAM)
  {
    messerr("When using Disjunctive Kriging, the Model be incremented");
    messerr("with Properties beforehad");
    messerr("For that sake, use 'model.properties'");
    goto label_end;
  }
  covanam = dynamic_cast<CovLMCAnamorphosis*>(model->getCovAnisoList());
  if (covanam == nullptr)
  {
    messerr("The Covariance does not seem to be a CovLMCAnamorphosis");
    goto label_end;
  }

  if (IFFFF(nfactor)) nfactor = covanam->getAnamNClass();
  if (covanam->getAnamType() == EAnam::HERMITIAN)
  {
    /* In the gaussian case, calculate the 'nfactor-1' factors */

    if (!DBIN->isVariableNumberComparedTo(1))
    {
      messerr("In Gaussian case, Input File must contain a single variable");
      goto label_end;
    }
    if (st_calculate_hermite_factors(DBIN, nfactor - 1)) goto label_end;
  }
  nvarz = DBIN->getVariableNumber();

  if (nfactor - 1 != nvarz)
  {
    messerr("The number of variables in Input Db (%d) does not match", nvarz);
    messerr("the number of factors (%d)-1", nfactor);
    goto label_end;
  }
  flag_block = 0;
  if (covanam->getAnamType() == EAnam::HERMITIAN)
  {
    AnamHermite *anam_hermite = dynamic_cast<AnamHermite*>(covanam->getAnam());
    if (anam_hermite->getRCoef() < 1.) flag_block = 1;
  }
  else if (covanam->getAnamType() == EAnam::DISCRETE_DD)
  {
    AnamDiscreteDD *anam_discrete_DD = dynamic_cast<AnamDiscreteDD*>(covanam->getAnam());
    if (anam_discrete_DD->getSCoef() > 0.) flag_block = 1;
  }
  else
  {
    AnamDiscreteIR *anam_discrete_IR = dynamic_cast<AnamDiscreteIR*>(covanam->getAnam());
    if (anam_discrete_IR->getRCoef() < 1.) flag_block = 1;
  }
  flag_panel = (flag_block && !nmult.empty());

  /* Add the attributes for storing the results */

  if (FLAG_EST != 0)
  {
    iptr_est_bck = DBOUT->addColumnsByConstant(nfactor - 1, TEST);
    if (iptr_est_bck < 0) goto label_end;
  }
  if (FLAG_STD != 0)
  {
    iptr_std_bck = DBOUT->addColumnsByConstant(nfactor - 1, TEST);
    if (iptr_std_bck < 0) goto label_end;
  }

  /* Pre-calculations */

  nbghw.initialize(DBIN, neighparam);
  nbghw.setFlagSimu(FLAG_SIMU);
  if (st_model_manage(1, model)) goto label_end;
  if (st_krige_manage(1, model->getVariableNumber(), model, neighparam)) goto label_end;
  nb_mult = 1;
  if (flag_block)
  {
    if (flag_panel)
    {
      if (krige_koption_manage(1, 1, EKrigOpt::BLOCK, 1, nmult)) goto label_end;
      KOPTION->calcul = EKrigOpt::PONCTUAL;
      nb_mult = KOPTION->ntot;
    }
    else
    {
      if (krige_koption_manage(1, 1, EKrigOpt::PONCTUAL, 1, VectorInt()))
        goto label_end;
    }
  }
  else
  {
    if (krige_koption_manage(1, 1, EKrigOpt::BLOCK, 0, ndisc)) goto label_end;
  }

  /* Centering the data */

  if (flag_block)
  {
    if (flag_panel)
    {
      dbsmu = db_create_grid_divider(dbgrid, nmult, 1);
      if (db_center_point_to_grid(DBIN, dbsmu, perturb)) goto label_end;
      dbsmu = db_delete(dbsmu);
    }
    else
    {
      if (db_center_point_to_grid(DBIN, dbgrid, perturb)) goto label_end;
    }
  }

  /* Core allocation */

  nvar = model->getVariableNumber();
  nfeq = model->getDriftEquationNumber();
  neqmax = 0;
  if (flag_panel)
  {
    neqmax = nvar * st_get_nmax(neighparam) + model->getDriftEquationNumber();
    rhs_cum = st_core(neqmax, 1);
  }
  varloc = (int*) mem_alloc(sizeof(int) * nvarz, 1);
  for (ivar = 0; ivar < nvarz; ivar++)
    varloc[ivar] = DBIN->getColIdxByLocator(ELoc::Z, ivar);

  /* Loop on the targets to be processed */

  status = 0;
  for (IECH_OUT = 0; IECH_OUT < DBOUT->getSampleNumber(); IECH_OUT++)
  {
    mes_process("Disjunctive Kriging for cell", DBOUT->getSampleNumber(),
                IECH_OUT);
    OptDbg::setIndex(IECH_OUT + 1);
    if (!DBOUT->isActive(IECH_OUT)) continue;
    if (OptDbg::query(EDbg::KRIGING) || OptDbg::query(EDbg::NBGH) || OptDbg::query(EDbg::RESULTS))
    {
      mestitle(1, "Target location");
      db_sample_print(DBOUT, IECH_OUT, 1, 0, 0);
    }

    /* Select the Neighborhood */

    DBIN->clearLocators(ELoc::Z);
    DBIN->setLocatorByUID(varloc[0], ELoc::Z);
    nbgh_ranks = nbghw.select(DBOUT,  IECH_OUT);
    nech = (int) nbgh_ranks.size();
    status = nbgh_ranks.empty();
    if (status) continue;

    /* Loop on the effective factors */

    flag_continue = 1;
    for (iclass = 1; iclass < nfactor && flag_continue; iclass++)
    {
      if (FLAG_EST != 0) IPTR_EST = iptr_est_bck + iclass - 1;
      if (FLAG_STD != 0) IPTR_STD = iptr_std_bck + iclass - 1;
      DBIN->clearLocators(ELoc::Z);
      DBIN->setLocatorByUID(varloc[iclass - 1], ELoc::Z);

      /* Set the rank of the current factor in the model */

      if (model_anamorphosis_set_factor(model, iclass)) goto label_end;

      /* Constant Calculation for the variance */

      if (FLAG_STD != 0) st_variance0(model, nvar, VectorVectorDouble());

      /* Establish the kriging L.H.S. (always performed) */

      st_prepar(model, neighparam, nbgh_ranks, &status, &nred, &neq);
      if (status)
      {
        flag_continue = 0;
        continue;
      }
      st_data_dual(model, NULL, nbgh_ranks, nred, &ldum);

      /* Blank out the estimate and the R.H.S. */

      if (flag_panel) for (i = 0; i < neqmax; i++)
        rhs_cum[i] = 0.;

      /* Loop on Panels (*1) or SMU (*nb_mult) */

      for (imult = 0; imult < nb_mult; imult++)
      {

        /* Establish the kriging R.H.S. */

        if (flag_panel) RAND_INDEX = imult;
        st_rhs(model, nbgh_ranks, neq, nvar, VectorVectorDouble(), &status);
        if (status) goto label_store;
        st_rhs_iso2hetero(neq, nvar);

        /* Cumulate the R.H.S. */

        if (flag_panel) for (i = 0; i < neqmax; i++)
          rhs_cum[i] += rhs[i];
      }

      /* Scale the R.H.S. */

      if (flag_panel) for (i = 0; i < neqmax; i++)
        rhs[i] = rhs_cum[i] / (double) nb_mult;
      if (OptDbg::query(EDbg::KRIGING))
        krige_rhs_print(nvar, nech, neq, nred, flag, rhs);

      /* Derive the kriging weights */

      if (FLAG_WGT)
      {
        matrix_product(nred, nred, nvar, lhs, rhs, wgt);
        if (OptDbg::query(EDbg::KRIGING))
          krige_wgt_print(status, nvar, nvar, nfeq, nbgh_ranks, nred, -1, flag, wgt);
      }

      /* Perform the estimation */

      label_store:
      st_estimate(model, NULL, status, neighparam->getFlagXvalid(),nvar, nred);
      if (OptDbg::query(EDbg::RESULTS))
        st_result_kriging_print(neighparam->getFlagXvalid(), nvar, status);
    }
  }

  /* Set the error return flag */

  error = 0;

  label_end: OptDbg::setIndex(0);
  rhs_cum = (double*) mem_free((char* ) rhs_cum);
  varloc = (int*) mem_free((char* ) varloc);
  (void) st_model_manage(-1, model);
  (void) st_krige_manage(-1, model->getVariableNumber(), model, neighparam);
  (void) krige_koption_manage(-1, 1, EKrigOpt::PONCTUAL, 1, ndisc);
  return (error);
}

/****************************************************************************/
/*!
 **  Perform the Neighborhood search
 **
 ** \return  Array of sample indices of the target neighbors
 **
 ** \param[in]  dbin       input Db structure
 ** \param[in]  model      Model structure
 ** \param[in]  neighparam ANeighParam structure
 ** \param[in]  target     Target location
 **
 ** \param[out] nech_out  Number of samples in the neighborhood
 **
 ** \remarks The resulting array must be freed by the calling procedure
 **
 ** \remark The number of variables in the 'dbin' may be different from
 ** \remark the number of variables in the 'model´.
 ** \remark This happens when the monovariate model is applied systematically
 ** \remark to all variables (such as for DK).
 ** \remark Dbin is modified so as to keep only the first Z-locator
 **
 *****************************************************************************/
int* neigh_calc(Db *dbin,
                Model *model,
                ANeighParam *neighparam,
                double *target,
                int *nech_out)
{
  int *neigh_tab, i, error, status, nech, zloc;
  Db *dbout;
  NeighWork nbghw;
  VectorInt nbgh_ranks;

  /* Preliminary checks */

  neigh_tab = nullptr;
  dbout = nullptr;
  *nech_out = 0;
  error = 1;

  /* Create a temporary dummy Db which contains the target */

  if (model == nullptr) goto label_end;
  dbout = db_create_from_target(target, model->getDimensionNumber(), 1);
  if (dbout == nullptr) goto label_end;
  st_global_init(dbin, dbout);

  /* Modification of 'dbin' */

  if (dbin != nullptr && dbin->getVariableNumber() != model->getVariableNumber()
      && model->getVariableNumber() == 1)
  {
    zloc = dbin->getColIdxByLocator(ELoc::Z);
    dbin->clearLocators(ELoc::Z);
    dbin->setLocatorByUID(zloc, ELoc::Z);
  }
  if (st_check_environment(1, 1, model, neighparam)) goto label_end;

  /* Pre-calculations */

  nbghw.initialize(DBIN, neighparam);
  nbghw.setFlagSimu(FLAG_SIMU);
  if (st_model_manage(1, model)) goto label_end;
  if (st_krige_manage(1, model->getVariableNumber(), model, neighparam))
    goto label_end;

  /* Select the Neighborhood */

  IECH_OUT = 0;
  nbgh_ranks = nbghw.select(DBOUT,  IECH_OUT);
  nech = (int) nbgh_ranks.size();
  status = nbgh_ranks.empty();
  if (status != 0)
  {
    messerr("Neighborhood search failed");
    goto label_end;
  }

  /* Store the neighbor indices */

  neigh_tab = (int*) mem_alloc(sizeof(int) * nech, 1);
  for (i = 0; i < nech; i++)
    neigh_tab[i] = nbgh_ranks[i] + 1;
  *nech_out = nech;

  /* Set the error return flag */

  error = 0;

  label_end: if (error) neigh_tab = (int*) mem_free((char* ) neigh_tab);
  dbout = db_delete(dbout);
  (void) st_model_manage(-1, model);
  (void) st_krige_manage(-1, model->getVariableNumber(), model, neighparam);
  return (neigh_tab);
}

/****************************************************************************/
/*!
 **  Allocate a vector of sample ranks excluding already selected pivots
 **
 ** \return  Pointer to the newly create integer vector
 **
 ** \param[in]  nech       Number of samples
 ** \param[in]  nsize1     Number of exact pivots currently selected
 ** \param[in]  ranks1     Ranks of exact pivots
 ** \param[in]  nsize2     Number of ACP pivots currently selected
 ** \param[in]  ranks2     Ranks of ACP pivots
 **
 ** \remarks The output array must be free by the calling function
 **
 *****************************************************************************/
static int* st_ranks_other(int nech,
                           int nsize1,
                           int *ranks1,
                           int nsize2,
                           int *ranks2)
{
  int *rother, i;

  rother = (int*) mem_alloc(sizeof(int) * nech, 0);
  if (rother == nullptr) return (rother);

  for (i = 0; i < nech; i++)
    rother[i] = i;
  for (i = 0; i < nsize1; i++)
    rother[ranks1[i]] = -1;
  for (i = 0; i < nsize2; i++)
    rother[ranks2[i]] = -1;
  return (rother);
}

/****************************************************************************/
/*!
 **  Establishing the kriging system with exact and ACP points
 **
 ** \return  Error retun code
 **
 ** \param[in]  db        Db structure
 ** \param[in]  model     Model structure
 ** \param[in]  beta      Thresholding value
 ** \param[in]  nsize1    Number of samples (exact)
 ** \param[in]  ranks1    Ranks of samples (exact)
 ** \param[in]  nsize2    Number of samples (ACP)
 ** \param[in]  ranks2    Ranks of samples (ACP)
 ** \param[in]  rother    Ranks of the idle samples (modified by routine)
 **
 ** \param[out] ntot_arg   Number of pivots
 ** \param[out] nutil_arg  Number of active samples
 ** \param[out] rutil_arg  Rank of the active samples
 ** \param[out] tutil_arg  Returned array for the U array
 ** \param[out] invsig_arg Returned array for Inverse Sigma
 **
 ** \remarks  The returned array must be freed by the calling function:
 ** \remarks  - rutil_arg  Array of integer - Dimension: nutil_arg
 ** \remarks  - tutil_arg  Array of float   - Dimension: nutil_arg * ntot_arg
 ** \remarks: - invsig_arg Array of float   - Dimension: ntot_arg * ntot_arg
 **
 *****************************************************************************/
static int st_sampling_krige_data(Db *db,
                                  Model *model,
                                  double beta,
                                  int nsize1,
                                  int *ranks1,
                                  int nsize2,
                                  int *ranks2,
                                  int *rother,
                                  int *ntot_arg,
                                  int *nutil_arg,
                                  int **rutil_arg,
                                  double **tutil_arg,
                                  double **invsig_arg)
{
  int *isort, *ralls, *rutil;
  int ndat, error, i, j, ntot, ntri, nother, npart, n1, ecr, nutil, nmax;
  double *utab, *s, *tl, *xl, *c, *sq, *v, *tn1, *tn2, *eigval, *eigvec, *spart,
      *vsort;
  double *tutil, *invsig, sumval;

  /* Initializations */

  error = 1;
  utab = s = tl = xl = c = v = tn1 = tn2 = sq = nullptr;
  eigval = eigvec = spart = vsort = tutil = invsig = nullptr;
  isort = ralls = rutil = nullptr;
  ndat = db->getSampleNumber(true);
  ntot = nsize1 + nsize2;
  ntri = nsize2 * (nsize2 + 1) / 2;
  nother = ndat - ntot;
  npart = ndat - nsize1;
  nutil = 0;

  /* Core allocation */

  utab = (double*) mem_alloc(sizeof(double) * ndat * ntot, 0);
  if (utab == nullptr) goto label_end;
  for (i = 0; i < ndat * ntot; i++)
    utab[i] = 0.;
  ralls = (int*) mem_alloc(sizeof(int) * ndat, 0);
  if (ralls == nullptr) goto label_end;

  /* Defining 'utab' for exact pivots */

  for (i = 0; i < nsize1; i++)
    UTAB(i,i) = 1.;
  ecr = 0;
  for (i = 0; i < nsize1; i++)
    ralls[ecr++] = ranks1[i];
  for (i = 0; i < nsize2; i++)
    ralls[ecr++] = ranks2[i];
  for (i = 0; i < ndat; i++)
    if (rother[i] >= 0) ralls[ecr++] = rother[i];

  /* Defining 'utab' for ACP pivots */

  if (nsize2 > 0)
  {
    tl = (double*) mem_alloc(sizeof(double) * ntri, 0);
    if (tl == nullptr) goto label_end;
    xl = (double*) mem_alloc(sizeof(double) * ntri, 0);
    if (xl == nullptr) goto label_end;
    v = (double*) mem_alloc(sizeof(double) * nother * nsize2, 0);
    if (v == nullptr) goto label_end;
    sq = (double*) mem_alloc(sizeof(double) * nsize2 * nsize2, 0);
    if (sq == nullptr) goto label_end;
    tn1 = (double*) mem_alloc(sizeof(double) * nsize2 * nsize2, 0);
    if (tn1 == nullptr) goto label_end;
    tn2 = (double*) mem_alloc(sizeof(double) * nsize2 * nsize2, 0);
    if (tn2 == nullptr) goto label_end;
    eigval = (double*) mem_alloc(sizeof(double) * nsize2, 0);
    if (eigval == nullptr) goto label_end;
    eigvec = (double*) mem_alloc(sizeof(double) * nsize2 * nsize2, 0);
    if (eigvec == nullptr) goto label_end;
    if (beta > 0.)
    {
      vsort = (double*) mem_alloc(sizeof(double) * npart, 0);
      if (vsort == nullptr) goto label_end;
      isort = (int*) mem_alloc(sizeof(double) * npart, 0);
      if (isort == nullptr) goto label_end;
    }

    s = model_covmat_by_ranks(model, db, nsize2, ranks2, db, nsize2, ranks2, -1,
                              -1, 0, 1);
    if (s == nullptr) goto label_end;
    if (matrix_cholesky_decompose(s, tl, nsize2)) goto label_end;
    matrix_triangle_to_square(0, nsize2, tl, sq);
    matrix_cholesky_invert(nsize2, tl, xl);
    c = model_covmat_by_ranks(model, db, nsize2, ranks2, db, ndat, rother, -1,
                              -1, 0, 1);
    if (c == nullptr) goto label_end;
    matrix_cholesky_product(4, nsize2, nother, xl, c, v);
    matrix_cholesky_norme(1, nsize2, tl, nullptr, tn1);
    if (matrix_prod_norme(-1, nother, nsize2, v, NULL, tn2)) goto label_end;
    matrix_combine(nsize2 * nsize2, 1, tn1, 1, tn2, tn1);
    if (matrix_eigen(tn1, nsize2, eigval, eigvec)) goto label_end;
    matrix_product_by_diag(3, nsize2, eigvec, eigval, eigvec);
    spart = matrix_bind(1, nsize2, nsize2, sq, nother, nsize2, v, &npart, &n1);
    if (spart == nullptr) goto label_end;
    matrix_product(npart, nsize2, nsize2, spart, eigvec, spart);

    if (beta > 0.)
    {
      for (i = 0; i < npart; i++)
      {
        sumval = 0.;
        for (j = 0; j < nsize2; j++)
          sumval = MAX(sumval, ABS(SPART(i,j)));
        vsort[i] = sumval;
        isort[i] = i;
      }
      ut_sort_double(1, npart, isort, vsort);
      nmax = MIN(npart, (int ) (beta * (double ) npart));
      for (i = 0; i < nmax; i++)
        for (j = 0; j < nsize2; j++)
          SPART(isort[i],j) = 0.;
    }

    for (i = 0; i < npart; i++)
      for (j = 0; j < nsize2; j++)
        UTAB(i+nsize1,j+nsize1) = -SPART(i, j);

    /* Core deallocation */

    tl = (double*) mem_free((char* ) tl);
    xl = (double*) mem_free((char* ) xl);
    v = (double*) mem_free((char* ) v);
    s = (double*) mem_free((char* ) s);
    c = (double*) mem_free((char* ) c);
    sq = (double*) mem_free((char* ) sq);
    tn1 = (double*) mem_free((char* ) tn1);
    tn2 = (double*) mem_free((char* ) tn2);
    spart = (double*) mem_free((char* ) spart);
    vsort = (double*) mem_free((char* ) vsort);
    eigval = (double*) mem_free((char* ) eigval);
    eigvec = (double*) mem_free((char* ) eigvec);
    isort = (int*) mem_free((char* ) isort);
  }

  /* Count the number of active samples */

  nutil = 0;
  for (i = 0; i < ndat; i++)
  {
    sumval = 0.;
    for (j = 0; j < ntot; j++)
      sumval += ABS(UTAB(i,j));
    if (sumval > 0.) nutil++;
  }

  /* Create the output arrays */

  rutil = (int*) mem_alloc(sizeof(int) * nutil, 0);
  if (rutil == nullptr) goto label_end;
  tutil = (double*) mem_alloc(sizeof(double) * ntot * nutil, 0);
  if (tutil == nullptr) goto label_end;
  invsig = (double*) mem_alloc(sizeof(double) * ntot * ntot, 0);
  if (invsig == nullptr) goto label_end;

  for (i = ecr = 0; i < ndat; i++)
  {
    sumval = 0.;
    for (j = 0; j < ntot; j++)
      sumval += ABS(UTAB(i,j));
    if (sumval <= 0.) continue;
    rutil[ecr] = ralls[i];
    for (j = 0; j < ntot; j++)
      TUTIL(ecr,j) = UTAB(i, j);
    ecr++;
  }
  s = model_covmat_by_ranks(model, db, nutil, rutil, db, nutil, rutil, -1, -1,
                            0, 1);
  if (s == nullptr) goto label_end;
  if (matrix_prod_norme(-1, nutil, ntot, tutil, s, invsig)) goto label_end;
  if (matrix_invert(invsig, ntot, 0)) goto label_end;
  s = (double*) mem_free((char* ) s);
  utab = (double*) mem_free((char* ) utab);
  ralls = (int*) mem_free((char* ) ralls);

  /* Returning arguments */

  *ntot_arg = ntot;
  *nutil_arg = nutil;
  *rutil_arg = rutil;
  *tutil_arg = tutil;
  *invsig_arg = invsig;

  /* Error return code */

  error = 0;

  label_end: utab = (double*) mem_free((char* ) utab);
  tl = (double*) mem_free((char* ) tl);
  xl = (double*) mem_free((char* ) xl);
  v = (double*) mem_free((char* ) v);
  s = (double*) mem_free((char* ) s);
  c = (double*) mem_free((char* ) c);
  sq = (double*) mem_free((char* ) sq);
  tn1 = (double*) mem_free((char* ) tn1);
  tn2 = (double*) mem_free((char* ) tn2);
  spart = (double*) mem_free((char* ) spart);
  vsort = (double*) mem_free((char* ) vsort);
  eigval = (double*) mem_free((char* ) eigval);
  eigvec = (double*) mem_free((char* ) eigvec);
  isort = (int*) mem_free((char* ) isort);
  ralls = (int*) mem_free((char* ) ralls);
  return (error);
}

/****************************************************************************/
/*!
 **  Perform the estimation at the data points
 **
 ** \return  Error retun code
 **
 ** \param[in]  db         Db structure
 ** \param[in]  model      Model structure
 ** \param[in]  beta       Thresholding value
 ** \param[in]  nsize1     Number of exact pivots currently selected
 ** \param[in]  ranks1     Ranks of exact pivots
 ** \param[in]  nsize2     Number of ACP pivots currently selected
 ** \param[in]  ranks2     Ranks of ACP pivots
 ** \param[in]  rother     Ranks of the idle samples
 ** \param[in]  flag_abs   1 Modify 'daata_est' to store the estimation error
 **
 ** \param[out] data_est  Array of estimation at samples
 ** \param[out] data_var  Array of estimation variance at samples
 **
 *****************************************************************************/
int st_krige_data(Db *db,
                  Model *model,
                  double beta,
                  int nsize1,
                  int *ranks1,
                  int nsize2,
                  int *ranks2,
                  int *rother,
                  int flag_abs,
                  double *data_est,
                  double *data_var)
{
  int *rutil, error, ntot, nutil, i, iech, nech;
  double *data, *tutil, *invsig, *s, *datm, *aux1, *aux2, *aux3, *aux4, *c00;
  double estim, variance, true_value;

  /* Initializations */

  error = 1;
  rutil = nullptr;
  tutil = invsig = data = datm = s = c00 = nullptr;
  aux1 = aux2 = aux3 = aux4 = nullptr;

  /* Core allocation */

  nutil = ntot = 0;
  nech = db->getSampleNumber();
  data = db_vector_alloc(db);
  if (data == nullptr) goto label_end;

  /* Perform local sampling */

  if (st_sampling_krige_data(db, model, beta, nsize1, ranks1, nsize2, ranks2,
                             rother, &ntot, &nutil, &rutil, &tutil, &invsig))
    goto label_end;

  /* Second core allocation */

  datm = (double*) mem_alloc(sizeof(double) * nutil, 0);
  if (datm == nullptr) goto label_end;
  aux1 = (double*) mem_alloc(sizeof(double) * ntot, 0);
  if (aux1 == nullptr) goto label_end;
  aux2 = (double*) mem_alloc(sizeof(double) * ntot, 0);
  if (aux2 == nullptr) goto label_end;
  aux3 = (double*) mem_alloc(sizeof(double) * ntot, 0);
  if (aux3 == nullptr) goto label_end;
  aux4 = (double*) mem_alloc(sizeof(double) * ntot, 0);
  if (aux4 == nullptr) goto label_end;

  /* Get the vector of active data and subtract the mean */

  if (db_vector_get(db, ELoc::Z, 0, data)) goto label_end;
  for (i = 0; i < nutil; i++)
    datm[i] = data[rutil[i]] - model->getMean(0);
  matrix_product(1, nutil, ntot, datm, tutil, aux1);
  matrix_product(1, ntot, ntot, aux1, invsig, aux2);

  /* Perform the estimation at all non pivot samples */

  for (iech = 0; iech < nech; iech++)
  {
    data_est[iech] = data_var[iech] = TEST;
    if (!db->isActive(iech)) continue;
    if (rother[iech] < 0) continue;
    c00 = model_covmat_by_ranks(model, db, 1, &iech, db, 1, &iech, -1, -1, 0,
                                1);
    if (c00 == nullptr) goto label_end;
    s = model_covmat_by_ranks(model, db, nutil, rutil, db, 1, &iech, -1, -1, 0,
                              1);
    if (s == nullptr) goto label_end;

    matrix_product(1, nutil, ntot, s, tutil, aux3);
    matrix_product(1, ntot, 1, aux2, aux3, &estim);
    data_est[iech] = estim + model->getMean(0);

    if (flag_abs)
    {
      true_value = db->getVariable(iech, 0);
      if (FFFF(true_value))
        data_est[iech] = TEST;
      else
        data_est[iech] = ABS(data_est[iech] - true_value);
    }

    matrix_product(1, ntot, ntot, aux3, invsig, aux4);
    matrix_product(1, ntot, 1, aux3, aux4, &variance);
    data_var[iech] = c00[0] - variance;

    s = (double*) mem_free((char* ) s);
    c00 = (double*) mem_free((char* ) c00);
  }

  /* Error return code */

  error = 0;

  label_end: data = db_vector_free(data);
  rutil = (int*) mem_free((char* ) rutil);
  tutil = (double*) mem_free((char* ) tutil);
  invsig = (double*) mem_free((char* ) invsig);
  datm = (double*) mem_free((char* ) datm);
  s = (double*) mem_free((char* ) s);
  c00 = (double*) mem_free((char* ) c00);
  aux1 = (double*) mem_free((char* ) aux1);
  aux2 = (double*) mem_free((char* ) aux2);
  aux3 = (double*) mem_free((char* ) aux3);
  aux4 = (double*) mem_free((char* ) aux4);
  return (error);
}

/****************************************************************************/
/*!
 **  Evaluate the improvement in adding a new pivot on the global score
 **
 ** \return  Error retun code
 **
 ** \param[in]  db         Db structure
 ** \param[in]  model      Model structure
 ** \param[in]  nsize1     Number of exact pivots currently selected
 ** \param[in]  ranks1     Ranks of exact pivots
 ** \param[in]  rother     Ranks of the idle samples
 **
 ** \param[out] crit       Array of criterion
 **
 *****************************************************************************/
int st_crit_global(Db *db,
                   Model *model,
                   int nsize1,
                   int *ranks1,
                   int *rother,
                   double *crit)
{
  int error, ndat, i, iech, nutil, ecr;
  double *c00, *invc, *data, *datm, *cs, *temp, *olderr, *olddiv, *aux1, *cs1;
  double *temp_loc, estim, sigma, value;

  /* Initializations */

  error = 1;
  ndat = db->getSampleNumber(true);
  nutil = ndat - nsize1;
  c00 = invc = data = datm = cs = temp = olderr = olddiv = nullptr;
  aux1 = cs1 = nullptr;

  /* Preliminary checks */

  if (nsize1 <= 0) goto label_end;

  /* Core allocation */

  data = db_vector_alloc(db);
  if (data == nullptr) goto label_end;
  datm = (double*) mem_alloc(sizeof(double) * ndat, 0);
  if (datm == nullptr) goto label_end;
  olderr = (double*) mem_alloc(sizeof(double) * nutil, 0);
  if (olderr == nullptr) goto label_end;
  olddiv = (double*) mem_alloc(sizeof(double) * nutil, 0);
  if (olddiv == nullptr) goto label_end;
  temp = (double*) mem_alloc(sizeof(double) * nsize1 * nutil, 0);
  if (temp == nullptr) goto label_end;
  aux1 = (double*) mem_alloc(sizeof(double) * nutil, 0);
  if (aux1 == nullptr) goto label_end;

  /* Establish the Kriging matrix on the pivot samples */

  invc = model_covmat_by_ranks(model, db, nsize1, ranks1, db, nsize1, ranks1,
                               -1, -1, 0, 1);
  if (invc == nullptr) goto label_end;
  if (matrix_invert(invc, nsize1, 0)) goto label_end;

  /* Set the data vector (corrected by the mean */

  if (db_vector_get(db, ELoc::Z, 0, data)) goto label_end;
  for (i = 0; i < nsize1; i++)
    datm[i] = data[ranks1[i]] - model->getMean(0);

  /* Loop on the non-pivots */

  for (iech = ecr = 0; iech < ndat; iech++)
  {
    temp_loc = &temp[ecr * nsize1];
    if (!db->isActive(iech)) continue;
    if (rother[iech] < 0) continue;

    c00 = model_covmat_by_ranks(model, db, 1, &iech, db, 1, &iech, -1, -1, 0,
                                1);
    if (c00 == nullptr) goto label_end;

    cs = model_covmat_by_ranks(model, db, nsize1, ranks1, db, 1, &iech, -1, -1,
                               0, 1);
    if (cs == nullptr) goto label_end;

    matrix_product(nsize1, nsize1, 1, invc, cs, temp_loc);
    matrix_product(1, nsize1, 1, datm, temp_loc, &estim);
    olderr[ecr] = estim + model->getMean(0) - db->getVariable(iech, 0);

    matrix_product(1, nsize1, 1, cs, temp_loc, &sigma);
    olddiv[ecr] = olderr[ecr] / (c00[0] - sigma);

    c00 = (double*) mem_free((char* ) c00);
    cs = (double*) mem_free((char* ) cs);
    ecr++;
  }

  /* Loop on the candidates */

  for (iech = ecr = 0; iech < ndat; iech++)
  {
    crit[iech] = TEST;
    if (!db->isActive(iech)) continue;
    if (rother[iech] < 0) continue;

    cs = model_covmat_by_ranks(model, db, 1, &iech, db, nsize1, ranks1, -1, -1,
                               0, 1);
    if (cs == nullptr) goto label_end;

    cs1 = model_covmat_by_ranks(model, db, 1, &iech, db, ndat, rother, -1, -1,
                                0, 1);
    if (cs1 == nullptr) goto label_end;

    matrix_product(1, nsize1, nutil, cs, temp, aux1);
    matrix_combine(nutil, 1, cs1, -1, aux1, cs1);
    matrix_combine(nutil, 1, olderr, -olddiv[ecr], cs1, cs1);

    value = 0.;
    for (i = 0; i < nutil; i++)
      value += cs1[i] * cs1[i];
    crit[iech] = value / nutil;

    cs = (double*) mem_free((char* ) cs);
    cs1 = (double*) mem_free((char* ) cs1);
    ecr++;
  }

  /* Set the error return code */

  error = 0;

  label_end: data = db_vector_free(data);
  c00 = (double*) mem_free((char* ) c00);
  invc = (double*) mem_free((char* ) invc);
  datm = (double*) mem_free((char* ) datm);
  cs = (double*) mem_free((char* ) cs);
  cs1 = (double*) mem_free((char* ) cs1);
  temp = (double*) mem_free((char* ) temp);
  aux1 = (double*) mem_free((char* ) aux1);
  olderr = (double*) mem_free((char* ) olderr);
  olddiv = (double*) mem_free((char* ) olddiv);
  return (error);
}

/****************************************************************************/
/*!
 **  Optimize the sampling design
 **
 ** \return  Error retun code
 **
 ** \param[in]  db         Db structure
 ** \param[in]  model      Model structure
 ** \param[in]  beta       Thresholding value
 ** \param[in]  method1    Criterion for choosing exact pivots
 **                        1 : Local evaluation
 **                        2 : Global evaluation
 ** \param[in]  nsize1_max Maximum number of exact pivots
 ** \param[in]  nsize1     Number of exact pivots currently selected
 ** \param[in]  ranks1     Ranks of exact pivots
 ** \param[in]  method2    Criterion for choosing ACP pivots
 **                        1 : Local evaluation
 **                        2 : Global evaluation
 ** \param[in]  nsize2_max Maximum number of ACP pivots
 ** \param[in]  nsize2     Number of ACP pivots currently selected
 ** \param[in]  ranks2     Ranks of ACP pivots
 ** \param[in]  verbose    1 for a verbose output
 **
 *****************************************************************************/
int sampling_f(Db *db,
               Model *model,
               double beta,
               int method1,
               int nsize1_max,
               int nsize1,
               int *ranks1,
               int method2,
               int nsize2_max,
               int nsize2,
               int *ranks2,
               int verbose)
{
  int *rother, error, best_rank, nech, nval;
  double *data_est, *data_var;
  double best_ecart, minimum, maximum, mean, stdv, delta;

  /* Initializations */

  error = 1;
  data_est = data_var = nullptr;
  rother = nullptr;
  nech = db->getSampleNumber();

  /* Preliminary checks */

  if (method2 != 1)
  {
    messerr("The Global Evaluation method for choosing ACP pivots");
    messerr("has not been programmed yet");
    goto label_end;
  }
  if (nsize1_max > 0 && nsize1 == 0)
  {
    messerr("The sampling requires a first sample to be defined 'ranks1'");
    goto label_end;
  }

  /* Core allocation */

  data_est = db_vector_alloc(db);
  if (data_est == nullptr) goto label_end;
  data_var = db_vector_alloc(db);
  if (data_var == nullptr) goto label_end;
  rother = st_ranks_other(nech, nsize1, ranks1, nsize2, ranks2);
  if (rother == nullptr) goto label_end;

  /* Sample the exact pivots */

  while (nsize1 < nsize1_max)
  {
    if (method1 == 1)
    {
      if (st_krige_data(db, model, beta, nsize1, ranks1, nsize2, ranks2, rother,
                        1, data_est, data_var)) goto label_end;
      best_rank = matrix_get_extreme(2, nech, data_est);
      best_ecart = data_est[best_rank];
    }
    else
    {
      if (st_crit_global(db, model, nsize1, ranks1, rother, data_est))
        goto label_end;
      best_rank = matrix_get_extreme(1, nech, data_est);
      best_ecart = data_est[best_rank];
    }
    if (verbose)
      message("Exact Pivots (%3d/%3d): Rank = %3d - value = %lf\n", nsize1 + 1,
              nsize1_max, best_rank + 1, best_ecart);
    ranks1[nsize1] = best_rank;
    rother[best_rank] = -1;
    nsize1++;
  }

  /* Sample the ACP pivots */

  while (nsize2 < nsize2_max)
  {
    if (st_krige_data(db, model, beta, nsize1, ranks1, nsize2, ranks2, rother,
                      1, data_est, data_var)) goto label_end;
    best_rank = matrix_get_extreme(2, nech, data_est);
    best_ecart = data_est[best_rank];
    if (verbose)
      message("ACP   Pivots (%3d/%3d): Rank = %3d - value = %lf\n", nsize2 + 1,
              nsize2_max, best_rank + 1, best_ecart);
    ranks2[nsize2] = best_rank;
    rother[best_rank] = -1;
    nsize2++;
  }

  /* Calculation of statistics on reproduction errors */

  if (verbose)
  {
    if (st_krige_data(db, model, beta, nsize1, ranks1, nsize2, ranks2, rother,
                      1, data_est, data_var)) goto label_end;
    ut_statistics(nech, data_est, NULL, NULL, &nval, &minimum, &maximum, &delta,
                  &mean, &stdv);
    mestitle(1, "Statistics on estimation errors");
    message("Count   = %d \n", nval);
    message("Minimum = %lf\n", minimum);
    message("Mean    = %lf\n", mean);
    message("St. Dev.= %lf\n", stdv);
    message("Maximum = %lf\n", maximum);
  }

  /* Error return code */

  error = 0;

  label_end: data_est = db_vector_free(data_est);
  data_var = db_vector_free(data_var);
  rother = (int*) mem_free((char* ) rother);
  return (error);
}

/****************************************************************************/
/*!
 **  Perform the Kriging procedure using the parcimonious search
 **  within the whole input data set
 **
 ** \return  Error retun code
 **
 ** \param[in]  dbin       Input Db structure
 ** \param[in]  dbout      Output Db structure
 ** \param[in]  model      Model structure
 ** \param[in]  beta       Thresholding value
 ** \param[in]  nsize1     Number of exact pivots currently selected
 ** \param[in]  ranks1     Ranks of exact pivots
 ** \param[in]  nsize2     Number of ACP pivots currently selected
 ** \param[in]  ranks2     Ranks of ACP pivots
 ** \param[in]  flag_std   Option for storing the standard deviation
 ** \param[in]  verbose    Verbose flag
 **
 *****************************************************************************/
int krigsampling_f(Db *dbin,
                   Db *dbout,
                   Model *model,
                   double beta,
                   int nsize1,
                   int *ranks1,
                   int nsize2,
                   int *ranks2,
                   int flag_std,
                   int verbose)
{
  int *rutil, *rother, error, nvar, ntot, nutil, i, nech;
  double *tutil, *data, *invsig, *datm, *aux1, *aux2, *aux3, *aux4, *s, *c00;
  double estim, sigma;

  /* Preliminary checks */

  error = 1;
  sigma = 0.;
  rutil = rother = nullptr;
  tutil = invsig = data = datm = s = c00 = nullptr;
  aux1 = aux2 = aux3 = aux4 = nullptr;
  st_global_init(dbin, dbout);
  FLAG_EST = 1;
  FLAG_STD = flag_std;
  if (st_check_environment(1, 1, model, NULL)) goto label_end;
  nvar = model->getVariableNumber();
  nech = dbin->getSampleNumber();

  /* Preliminary checks */

  if (nvar != 1)
  {
    messerr("This method is only programmed for monovariate case");
    goto label_end;
  }
  if (nsize1 + nsize2 <= 0)
  {
    messerr("You must specify some pivots in 'ranks1' or 'ranks2'");
    goto label_end;
  }

  /* Add the attributes for storing the results */

  if (FLAG_EST != 0)
  {
    IPTR_EST = dbout->addColumnsByConstant(nvar, 0.);
    if (IPTR_EST < 0) goto label_end;
  }
  if (FLAG_STD != 0)
  {
    IPTR_STD = dbout->addColumnsByConstant(nvar, 0.);
    if (IPTR_STD < 0) goto label_end;
  }

  /* Core allocation */

  rother = st_ranks_other(nech, nsize1, ranks1, nsize2, ranks2);
  if (rother == nullptr) goto label_end;

  /* Perform local sampling */

  if (st_sampling_krige_data(dbin, model, beta, nsize1, ranks1, nsize2, ranks2,
                             rother, &ntot, &nutil, &rutil, &tutil, &invsig))
    goto label_end;

  /* Optional printout */

  if (verbose)
  {
    message("Printout of intermediate arrays\n");
    print_imatrix("Pivot ranks", 0, 1, 1, ntot, NULL, rutil);
    print_matrix("Inv-Sigma", 0, 1, ntot, ntot, NULL, invsig);
    print_matrix("U", 0, 1, ntot, nutil, NULL, tutil);
  }

  /* Second core allocation */

  data = db_vector_alloc(dbin);
  if (data == nullptr) goto label_end;
  datm = (double*) mem_alloc(sizeof(double) * nutil, 0);
  if (datm == nullptr) goto label_end;
  aux1 = (double*) mem_alloc(sizeof(double) * ntot, 0);
  if (aux1 == nullptr) goto label_end;
  aux2 = (double*) mem_alloc(sizeof(double) * ntot, 0);
  if (aux2 == nullptr) goto label_end;
  aux3 = (double*) mem_alloc(sizeof(double) * ntot, 0);
  if (aux3 == nullptr) goto label_end;
  if (FLAG_STD != 0)
  {
    aux4 = (double*) mem_alloc(sizeof(double) * ntot, 0);
    if (aux4 == nullptr) goto label_end;
  }

  /* Get the vector of active data and substract the mean */

  if (db_vector_get(dbin, ELoc::Z, 0, data)) goto label_end;
  for (i = 0; i < nutil; i++)
    datm[i] = data[rutil[i]] - model->getMean(0);
  matrix_product(1, nutil, ntot, datm, tutil, aux1);
  matrix_product(1, ntot, ntot, aux1, invsig, aux2);

  /* Loop on the target samples */

  for (IECH_OUT = 0; IECH_OUT < DBOUT->getSampleNumber(); IECH_OUT++)
  {
    mes_process("Kriging sample", DBOUT->getSampleNumber(), IECH_OUT);
    OptDbg::setIndex(IECH_OUT + 1);
    if (!dbout->isActive(IECH_OUT)) continue;
    if (OptDbg::query(EDbg::KRIGING) || OptDbg::query(EDbg::NBGH) || OptDbg::query(EDbg::RESULTS))
    {
      mestitle(1, "Target location");
      db_sample_print(dbout, IECH_OUT, 1, 0, 0);
    }

    s = model_covmat_by_ranks(model, dbin, nutil, rutil, dbout, 1, &IECH_OUT,
                              -1, -1, 0, 1);
    if (s == nullptr) goto label_end;
    if (FLAG_STD != 0)
    {
      c00 = model_covmat_by_ranks(model, dbout, 1, &IECH_OUT, dbout, 1,
                                  &IECH_OUT, -1, -1, 0, 1);
      if (c00 == nullptr) goto label_end;
    }

    matrix_product(1, nutil, ntot, s, tutil, aux3);
    matrix_product(1, ntot, 1, aux2, aux3, &estim);
    estim += model->getMean(0);
    DBOUT->setArray(IECH_OUT, IPTR_EST, estim);

    if (FLAG_STD != 0)
    {
      matrix_product(1, ntot, ntot, aux3, invsig, aux4);
      matrix_product(1, ntot, 1, aux3, aux4, &sigma);
      sigma = c00[0] - sigma;
      sigma = (sigma > 0) ? sqrt(sigma) :
                            0.;
      DBOUT->setArray(IECH_OUT, IPTR_STD, sigma);
    }

    /* Optional printout */

    if (OptDbg::query(EDbg::RESULTS))
    {
      tab_printg(" - Estimate  = ", estim);
      message("\n");
      if (FLAG_STD != 0)
      {
        tab_printg(" - Std. Dev. = ", sigma);
        message("\n");
      }
    }

    /* Core deallocation */

    s = (double*) mem_free((char* ) s);
    c00 = (double*) mem_free((char* ) c00);
  }

  /* Error return code */

  error = 0;

  label_end: rother = (int*) mem_free((char* ) rother);
  rutil = (int*) mem_free((char* ) rutil);
  tutil = (double*) mem_free((char* ) tutil);
  invsig = (double*) mem_free((char* ) invsig);
  data = (double*) mem_free((char* ) data);
  datm = (double*) mem_free((char* ) datm);
  s = (double*) mem_free((char* ) s);
  c00 = (double*) mem_free((char* ) c00);
  aux1 = (double*) mem_free((char* ) aux1);
  aux2 = (double*) mem_free((char* ) aux2);
  aux3 = (double*) mem_free((char* ) aux3);
  aux4 = (double*) mem_free((char* ) aux4);
  return (error);
}

/****************************************************************************/
/*!
 **  Display the statistics of the target variable
 **  before or after the Declustering
 **
 ** \param[in]  mode      0 before the declustering; 1 after
 ** \param[in]  method    Method for declustering
 ** \param[in]  db        input Db structure
 ** \param[in]  iptr      Rank of the weighting variable
 **
 *****************************************************************************/
static void st_declustering_stats(int mode, int method, Db *db, int iptr)
{
  double mean, var, zval, coeff, mini, maxi, sumwgt;

  mean = var = sumwgt = 0.;
  mini = 1.e30;
  maxi = -1.e30;
  for (int iech = 0; iech < db->getSampleNumber(); iech++)
  {
    if (!db->isActive(iech)) continue;
    zval = db->getVariable(iech, 0);
    if (FFFF(zval)) continue;
    coeff = (mode == 0) ? 1. :
                          db->getArray(iech, iptr);
    coeff = ABS(coeff);
    sumwgt += coeff;
    mean += coeff * zval;
    var += coeff * zval * zval;
    if (coeff < mini) mini = coeff;
    if (coeff > maxi) maxi = coeff;
  }

  /* Scaling */

  if (sumwgt > 0)
  {
    mean /= sumwgt;
    var = var / sumwgt - mean * mean;
  }

  if (mode == 0)
    mestitle(1, "Statistics before Declustering");
  else
    mestitle(1, "Statistics after Declustering");
  if (method == 1)
    message("- Using the Number of Samples per Neighborhood\n");
  else if (method == 2)
    message("- Using the weights for Kriging the Global Mean\n");
  else
    message("- Using the average weight for Kriging cells of a Grid\n");

  message("- Sum of weights    = %lf\n", sumwgt);
  message("- Mean              = %lf\n", mean);
  message("- Variance          = %lf\n", var);
  if (mode == 1)
  {
    message("- Minimum Weight    = %lf\n", mini);
    message("- Maximum Weight    = %lf\n", maxi);
  }
}

/****************************************************************************/
/*!
 **  Truncate the negative weights and scale the remaining ones
 **
 ** \param[in]  db        Db structure
 ** \param[in]  iptr      Rank of the Weight variable
 ** \param[in]  verbose   Verbose option
 **
 *****************************************************************************/
static void st_declustering_truncate(Db *db, int iptr, int verbose)
{
  double total, coeff;

  /* Truncate the negative weights */

  total = 0;
  for (int iech = 0; iech < db->getSampleNumber(); iech++)
  {
    if (!db->isActive(iech)) continue;
    if (FFFF(db->getVariable(iech, 0))) continue;
    coeff = db->getArray(iech, iptr);
    if (coeff < 0)
    {
      if (verbose)
        messerr("Weight #%d is negative (%lf). It has been set to 0", iech + 1,
                coeff);
      db->setArray(iech, iptr, 0.);
    }
    else
      total += coeff;
  }

  /* Rescale */

  for (int iech = 0; iech < db->getSampleNumber(); iech++)
  {
    if (!db->isActive(iech)) continue;
    if (FFFF(db->getVariable(iech, 0))) continue;
    db->updArray(iech, iptr, 3, total);
  }
}

/****************************************************************************/
/*!
 **  Perform the Declustering task (Number of samples within an ellipse)
 **
 ** \return  Error return code
 **
 ** \param[in]  db        input Db structure
 ** \param[in]  iptr      Rank of the declustering weight
 ** \param[in]  radius    Array of neighborhood radius
 **
 *****************************************************************************/
static int st_declustering_1(Db *db, int iptr, double *radius)
{
  int error;
  double *vect, dist, value, total;

  /* Initializations */

  error = 1;
  vect = nullptr;

  /* Core allocation */

  vect = db_sample_alloc(db, ELoc::X);
  if (vect == nullptr) goto label_end;

  /* Loop on the target sample */

  for (int iech = 0; iech < db->getSampleNumber(); iech++)
  {
    if (!db->isActive(iech)) continue;
    if (FFFF(db->getVariable(iech, 0))) continue;

    /* Loop on the second sample */

    for (int jech = 0; jech < db->getSampleNumber(); jech++)
    {
      if (!db->isActive(jech)) continue;
      value = db->getVariable(iech, 0);
      if (FFFF(value)) continue;
      (void) distance_intra(db, iech, jech, vect);

      /* Normalize the distance */

      dist = 0.;
      for (int idim = 0; idim < db->getNDim(); idim++)
      {
        vect[idim] /= radius[idim];
        dist += vect[idim] * vect[idim];
      }
      if (dist > 1) continue;
      db->updArray(iech, iptr, 0, 1);
    }
  }

  /* Normalization step */

  total = 0.;
  for (int iech = 0; iech < db->getSampleNumber(); iech++)
  {
    if (!db->isActive(iech)) continue;
    if (FFFF(db->getVariable(iech, 0))) continue;
    value = 1. / db->getArray(iech, iptr);
    total += value;
  }
  for (int iech = 0; iech < db->getSampleNumber(); iech++)
  {
    if (!db->isActive(iech)) continue;
    if (FFFF(db->getVariable(iech, 0))) continue;
    value = 1. / db->getArray(iech, iptr) / total;
    db->setArray(iech, iptr, value);
  }

  /* Set the error return code */

  error = 0;

  label_end: vect = db_sample_free(vect);
  return (error);
}

/****************************************************************************/
/*!
 **  Perform the Declustering task as weight of the Mean Kriging (Unique Neigh)
 **
 ** \return  Error return code
 **
 ** \param[in]  db        input Db structure
 ** \param[in]  iptr      Rank of the declustering weight
 ** \param[in]  model     Model structure
 ** \param[in]  verbose   Verbose option
 **
 *****************************************************************************/
static int st_declustering_2(Db *db, int iptr, Model *model, int verbose)
{
  int error, ndim, status, nech, nred, neq, ecr, nvar;
  NeighWork nbghw;
  VectorInt nbgh_ranks;

  /* Initializations */

  error = 1;
  ndim = db->getNDim();
  NeighUnique* neighU = NeighUnique::create(ndim, false);
  nvar = model->getVariableNumber();
  st_global_init(db, db);
  FLAG_EST = 0;
  FLAG_STD = 0;
  FLAG_VARZ = 0;
  FLAG_WGT = 1;
  if (st_check_environment(1, 1, model, neighU)) goto label_end;

  /* Pre-calculations */

  nbghw.initialize(DBIN, neighU);
  nbghw.setFlagSimu(FLAG_SIMU);
  if (st_model_manage(1, model)) goto label_end;
  if (st_krige_manage(1, nvar, model, neighU)) goto label_end;
  if (krige_koption_manage(1, 1, EKrigOpt::DRIFT, 1, VectorInt()))
    goto label_end;

  /* Prepare the Neighborhood */

  nbgh_ranks = nbghw.select(DBOUT,  IECH_OUT);
  nech = (int) nbgh_ranks.size();
  status = nbgh_ranks.empty();
  if (status) goto label_end;

  /* Establish the L.H.S. */

  st_prepar(model, neighU, nbgh_ranks, &status, &nred, &neq);
  if (status) goto label_end;

  /* Loop on the targets to be processed */

  status = 0;
  IECH_OUT = 0;
  st_rhs(model, nbgh_ranks, neq, nvar, VectorVectorDouble(), &status);
  if (status) goto label_end;
  st_rhs_iso2hetero(neq, 1);
  if (OptDbg::query(EDbg::KRIGING)) krige_rhs_print(1, nech, neq, nred, flag, rhs);

  /* Derive the kriging weights */
  matrix_product(nred, nred, 1, lhs, rhs, wgt);
  if (OptDbg::query(EDbg::KRIGING))
    krige_wgt_print(status, 1, 1, model->getDriftEquationNumber(), nbgh_ranks, nred,
                    -1, flag, wgt);

  /* Store the weights */

  ecr = 0;
  for (int iech = 0; iech < db->getSampleNumber(); iech++)
  {
    if (!db->isActive(iech)) continue;
    if (FFFF(db->getVariable(iech, 0))) continue;
    db->setArray(nbgh_ranks[ecr], iptr, wgt[ecr]);
    ecr++;
  }

  /* Truncate the negative weights */

  st_declustering_truncate(db, iptr, verbose);

  /* Set the error return flag */

  error = 0;

  label_end: OptDbg::setIndex(0);
  (void) st_model_manage(-1, model);
  (void) st_krige_manage(-1, nvar, model, neighU);
  (void) krige_koption_manage(-1, 1, EKrigOpt::DRIFT, 1, VectorInt());
  delete neighU;
  return (error);
}

/****************************************************************************/
/*!
 **  Perform the Declustering task as the sum of the weight
 **  for Kriging the Cells of a grid
 **
 ** \return  Error return code
 **
 ** \param[in]  db         input Db structure
 ** \param[in]  dbgrid     output Db structure
 ** \param[in]  iptr       Rank of the declustering weight
 ** \param[in]  model      Model structure
 ** \param[in]  neighparam ANeighParam structure
 ** \param[in]  ndisc      Array of discretization counts
 ** \param[in]  verbose    Verbose option
 **
 *****************************************************************************/
static int st_declustering_3(Db *db,
                             Db *dbgrid,
                             int iptr,
                             Model *model,
                             ANeighParam *neighparam,
                             VectorInt ndisc,
                             int verbose)
{
  int error, status, nech, nred, neq, nvar, ecr;
  double ldum;
  NeighWork nbghw;
  VectorInt nbgh_ranks;

  /* Initializations */

  error = 1;
  nvar = nred = neq = 0;
  st_global_init(db, dbgrid);
  FLAG_EST = 0;
  FLAG_STD = 0;
  FLAG_VARZ = 0;
  FLAG_WGT = 1;
  if (st_check_environment(1, 1, model, neighparam)) goto label_end;
  nvar = model->getVariableNumber();

  /* Pre-calculations */

  nbghw.initialize(DBIN, neighparam);
  nbghw.setFlagSimu(FLAG_SIMU);
  if (st_model_manage(1, model)) goto label_end;
  if (st_krige_manage(1, nvar, model, neighparam)) goto label_end;
  if (krige_koption_manage(1, 1, EKrigOpt::BLOCK, 1, ndisc)) goto label_end;

  /* Loop on the grid cells */

  status = 0;
  for (IECH_OUT = 0; IECH_OUT < dbgrid->getSampleNumber(); IECH_OUT++)
  {
    if (!DBOUT->isActive(IECH_OUT)) continue;

    /* Select the Neighborhood */

    nbgh_ranks = nbghw.select(DBOUT,  IECH_OUT);
    nech = (int) nbgh_ranks.size();
    status = nbgh_ranks.empty();
    if (status) continue;

    /* Establish the kriging L.H.S. */

    if (!nbghw.isUnchanged() || neighparam->getFlagContinuous() || OptDbg::force())
    {
      st_prepar(model, neighparam, nbgh_ranks, &status, &nred, &neq);
      if (status) continue;
      st_data_dual(model, NULL, nbgh_ranks, nred, &ldum);
    }

    /* Establish the Kriging R.H.S. */

    st_rhs(model, nbgh_ranks, neq, nvar, VectorVectorDouble(), &status);
    if (status) continue;
    st_rhs_iso2hetero(neq, 1);

    /* Derive the Kriging weights */

    matrix_product(nred, nred, 1, lhs, rhs, wgt);

    /* Cumulate the weights */

    ecr = 0;
    for (int iech = 0; iech < db->getSampleNumber(); iech++)
    {
      if (!db->isActive(iech)) continue;
      if (FFFF(db->getVariable(iech, 0))) continue;
      db->updArray(nbgh_ranks[ecr], iptr, 0, wgt[ecr]);
      ecr++;
    }
  }

  /* Truncate the negative weights */

  st_declustering_truncate(db, iptr, verbose);

  /* Set the error return flag */

  error = 0;

  label_end: OptDbg::setIndex(0);
  (void) st_model_manage(-1, model);
  (void) st_krige_manage(-1, nvar, model, neighparam);
  (void) krige_koption_manage(-1, 1, EKrigOpt::BLOCK, 1, ndisc);
  return (error);
}

/****************************************************************************/
/*!
 **  Perform the Declustering task
 **
 ** \return  Error return code
 **
 ** \param[in]  dbin       input Db structure
 ** \param[in]  model      Model structure
 ** \param[in]  neighparam ANeighParam structure
 ** \param[in]  dbgrid     Grid auxiliary Db structure
 ** \param[in]  method     Method for declustering
 ** \param[in]  radius     Array of neighborhood radius
 ** \param[in]  ndisc      Array of discretization
 ** \param[in]  flag_sel   1 to mask off samples with zero weight
 ** \param[in]  verbose    Verbose option
 **
 *****************************************************************************/
int declustering_f(Db *dbin,
                   Model *model,
                   ANeighParam *neighparam,
                   DbGrid *dbgrid,
                   int method,
                   double *radius,
                   VectorInt ndisc,
                   int flag_sel,
                   int verbose)
{
  int error, iptr, iptr_sel;
  double indic;

  /* Initializations */

  error = 1;

  /* Preliminary checks */

  if (dbin->isVariableNumberComparedTo(0)) goto label_end;

  /* Add the kriging weight as a new variable */

  iptr = dbin->addColumnsByConstant(1, 0.);
  if (iptr < 0) goto label_end;

  /* Produce statistics on the target variable before declustering */

  if (verbose) st_declustering_stats(0, method, dbin, iptr);

  /* Dispatch */

  switch (method)
  {
    case 1: /* Weight proportional to nb samples */
      if (st_declustering_1(dbin, iptr, radius)) goto label_end;
      break;

    case 2: /* Weight of the Mean */
      if (st_declustering_2(dbin, iptr, model, verbose)) goto label_end;
      break;

    case 3: /* Average weight of the Block Kriging */
      if (st_declustering_3(dbin, dbgrid, iptr, model, neighparam, ndisc, verbose))
        goto label_end;
      break;

    default:
      messerr("Not yet implemented");
      break;
  }

  /* Store the selection (optional) */

  if (flag_sel)
  {
    iptr_sel = dbin->addColumnsByConstant(1, 0.);
    if (iptr_sel < 0) goto label_end;
    for (int iech = 0; iech < dbin->getSampleNumber(); iech++)
    {
      dbin->setArray(iech, iptr_sel, 0.);
      if (!dbin->isActive(iech)) continue;
      indic = (dbin->getArray(iech, iptr) > 0.);
      dbin->setArray(iech, iptr_sel, indic);
    }
  }

  /* Produce statistics on the target variable after declustering */

  if (verbose) st_declustering_stats(1, method, dbin, iptr);

  /* Set the error return code */

  error = 0;

  label_end: return (error);
}

/****************************************************************************/
/*!
 **  Establish the covariance matrix between two Dbs
 **
 ** \return  Covariance matrix (Dim: n1 * n2)
 **
 ** \param[in]  title       Title of the optional printout
 ** \param[in]  db1         First Db structure
 ** \param[in]  test_def1   1 if the first variable (ELoc::Z) must be checked
 ** \param[in]  db2         Second Db structure
 ** \param[in]  test_def2   1 if the second variable (ELoc::Z) must be checked
 ** \param[in]  model       Model structure
 **
 ** \remarks The returned argument must be freed by the calling function
 **
 *****************************************************************************/
static double* st_calcul_covmat(const char *title,
                                Db *db1,
                                int test_def1,
                                Db *db2,
                                int test_def2,
                                Model *model)
{
  int n1, n2, i1, i2;
  double *covgen;
  CovCalcMode mode;

  /* Initializations */

  n1 = (test_def1) ? db1->getActiveAndDefinedNumber(0) :
                     db1->getSampleNumber(true);
  n2 = (test_def2) ? db2->getActiveAndDefinedNumber(0) :
                     db2->getSampleNumber(true);

  /* Core allocation */

  covgen = (double*) mem_alloc(sizeof(double) * n1 * n2, 0);
  if (covgen == nullptr) return (covgen);

  for (int ii1 = i1 = 0; ii1 < db1->getSampleNumber(); ii1++)
  {
    if (test_def1)
    {
      if (!db1->isActiveAndDefined(ii1, 0)) continue;
    }
    else
    {
      if (!db1->isActive(ii1)) continue;
    }

    for (int ii2 = i2 = 0; ii2 < db2->getSampleNumber(); ii2++)
    {
      if (test_def2)
      {
        if (!db2->isActiveAndDefined(ii2, 0)) continue;
      }
      else
      {
        if (!db2->isActive(ii2)) continue;
      }

      for (int idim = 0; idim < db1->getNDim(); idim++)
        d1[idim] = db1->getDistance1D(ii1, ii2, idim);

      model_calcul_cov(NULL,model, mode, 1, 1., d1, &COVGEN(i1, i2));
      i2++;
    }
    i1++;
  }

  /* Optional printout */

  if (INH_FLAG_VERBOSE)
    print_matrix(title, INH_FLAG_LIMIT, 1, n2, n1, NULL, covgen);

  return (covgen);
}

/****************************************************************************/
/*!
 **  Establish the drift matrix for a given Db
 **
 ** \return  Drift matrix (Dim: n1 * nbfl)
 **
 ** \param[in]  title       Title of the optionla printout
 ** \param[in]  db1         First Db structure
 ** \param[in]  test_def1   1 if the first variable (ELoc::Z) must be checked
 ** \param[in]  model       Model structure
 **
 ** \remarks The returned argument must be freed by the calling function
 **
 *****************************************************************************/
static double* st_calcul_drfmat(const char *title,
                                Db *db1,
                                int test_def1,
                                Model *model)
{
  int i1, n1, nbfl;
  double *drftab;

  /* Initializations */

  n1 = (test_def1) ? db1->getActiveAndDefinedNumber(0) :
                     db1->getSampleNumber(true);
  nbfl = model->getDriftNumber();

  /* Core allocation */

  drftab = (double*) mem_alloc(sizeof(double) * n1 * nbfl, 0);
  if (drftab == nullptr) return (drftab);

  /* Loop on the samples */

  i1 = 0;
  for (int ii1 = 0; ii1 < db1->getSampleNumber(); ii1++)
  {
    if (test_def1)
    {
      if (!db1->isActiveAndDefined(ii1, 0)) continue;
    }
    else
    {
      if (!db1->isActive(ii1)) continue;
    }

    model_calcul_drift(model, ECalcMember::LHS, db1, ii1, &drftab[i1 * nbfl]);
    i1++;
  }

  /* Optional printout */

  if (INH_FLAG_VERBOSE)
    print_matrix(title, INH_FLAG_LIMIT, 1, nbfl, n1, NULL, drftab);

  return (drftab);
}

/****************************************************************************/
/*!
 **  Establish the distance matrix between two Dbs
 **
 ** \return  Covariance matrix
 **
 ** \param[in]  title       Title of the optional printout
 ** \param[in]  db1         First Db structure
 ** \param[in]  test_def1   1 if the first variable (ELoc::Z) must be checked
 ** \param[in]  db2         Second Db structure (sources)
 ** \param[in]  test_def2   1 if the second variable (ELoc::Z) must be checked
 ** \param[in]  power       Power of the Distance decay
 **
 ** \remarks The returned argument must be freed by the calling function
 **
 *****************************************************************************/
static double* st_calcul_distmat(const char *title,
                                 Db *db1,
                                 int test_def1,
                                 Db *db2,
                                 int test_def2,
                                 double power)
{
  int n1, ns, i1, is, ndim;
  double *distgen, dist;

  /* Initializations */

  n1 = (test_def1) ? db1->getActiveAndDefinedNumber(0) :
                     db1->getSampleNumber(true);
  ns = (test_def2) ? db2->getActiveAndDefinedNumber(0) :
                     db2->getSampleNumber(true);
  ndim = db1->getNDim();

  /* Core allocation */

  distgen = (double*) mem_alloc(sizeof(double) * n1 * ns, 0);
  if (distgen == nullptr) return (distgen);

  for (int ii1 = i1 = 0; ii1 < db1->getSampleNumber(); ii1++)
  {
    if (test_def1)
    {
      if (!db1->isActiveAndDefined(ii1, 0)) continue;
    }
    else
    {
      if (!db1->isActive(ii1)) continue;
    }

    for (int iis = is = 0; iis < db2->getSampleNumber(); iis++)
    {
      if (test_def2)
      {
        if (!db2->isActiveAndDefined(iis, 0)) continue;
      }
      else
      {
        if (!db2->isActive(iis)) continue;
      }

      dist = 0.;
      for (int idim = 0; idim < ndim; idim++)
      {
        d1[idim] = db1->getDistance1D(ii1, iis, idim);
        dist += d1[idim] * d1[idim];
      }

      DISTGEN(i1,is) = 1. / pow(dist, power / 2.);
      is++;
    }
    i1++;
  }

  /* Optional printout */

  if (INH_FLAG_VERBOSE)
    print_matrix(title, INH_FLAG_LIMIT, 1, ns, n1, NULL, distgen);

  return (distgen);
}

/****************************************************************************/
/*!
 **  Operate the product of the Distance by the Source covariance matrix
 **
 ** \return  Product matrix
 **
 ** \param[in]  title       Title of the optionla printout
 ** \param[in]  n1          Number of Data
 ** \param[in]  ns          Number of Sources
 ** \param[in]  covss       Covariance Matrix between Sources (Dim: ns*ns)
 ** \param[in]  distgen     Distance matrix
 **
 ** \remarks The returned argument must be freed by the calling function
 **
 *****************************************************************************/
static double* st_calcul_product(const char *title,
                                 int n1,
                                 int ns,
                                 double *covss,
                                 double *distgen)
{
  double *prodgen;

  prodgen = (double*) mem_alloc(sizeof(double) * n1 * ns, 0);
  if (prodgen == nullptr) return (prodgen);

  for (int i1 = 0; i1 < n1; i1++)
    for (int is = 0; is < ns; is++)
    {
      PRODGEN(i1,is) = 0.;
      for (int js = 0; js < ns; js++)
        PRODGEN(i1,is) += COVSS(is,js) * DISTGEN(i1, js);
    }

  /* Optional printout */

  if (INH_FLAG_VERBOSE)
    print_matrix(title, INH_FLAG_LIMIT, 1, ns, n1, NULL, prodgen);

  return (prodgen);
}

/****************************************************************************/
/*!
 **  Establish the L.H.S. for Inhomogeonous Kriging
 **
 ** \return  Covariance matrix
 **
 ** \param[in]  dbdat       Db structure containing Data
 ** \param[in]  dbsrc       Db structure containing Sources
 ** \param[in]  model_dat   Model structure for the data
 ** \param[in]  distps      Distance matrix between Data and Sources
 ** \param[in]  prodps      Product of DistPS by CovSS
 **
 *****************************************************************************/
static double* st_inhomogeneous_covpp(Db *dbdat,
                                      Db *dbsrc,
                                      Model *model_dat,
                                      double *distps,
                                      double *prodps)
{
  double *covpp;
  int np, ns, error;

  /* Initializations */

  error = 1;
  covpp = nullptr;

  np = dbdat->getActiveAndDefinedNumber(0);
  ns = dbsrc->getSampleNumber(true);

  /* Covariance matrix between Mesures */

  covpp = st_calcul_covmat("Covariance P-P", dbdat, 1, dbdat, 1, model_dat);
  if (covpp == nullptr) goto label_end;

  /* Calculate the LHS matrix */

  for (int ip = 0; ip < np; ip++)
    for (int jp = ip; jp < np; jp++)
      for (int is = 0; is < ns; is++)
      {
        COVPP(ip,jp) += DISTPS(ip,is) * PRODPS(jp, is);
        if (jp > ip) COVPP(jp,ip) = COVPP(ip, jp);
      }

  /* Set the error return code */

  error = 0;

  label_end: if (error) covpp = (double*) mem_free((char* ) covpp);
  return (covpp);
}

/****************************************************************************/
/*!
 **  Establish the R.H.S. for Inhomogeneous Kriging
 **
 ** \return  Covariance matrix
 **
 ** \param[in]  dbdat       Db structure containing Data
 ** \param[in]  dbsrc       Db structure containing Sources
 ** \param[in]  dbout       Db structure containing Targets
 ** \param[in]  flag_source If the result is the source, rather than diffusion
 ** \param[in]  model_dat   Model structure for the data
 ** \param[in]  distps      Distance matrix between Data and Sources
 ** \param[in]  prodps      Product of DistPS by CovSS
 ** \param[in]  prodgs      Product of DistGS by CovSS
 **
 ** \remarks: When used with flag_source=TRUE, 'dbsrc' and 'dbout' coincide
 **
 *****************************************************************************/
static double* st_inhomogeneous_covgp(Db *dbdat,
                                      Db *dbsrc,
                                      Db *dbout,
                                      int flag_source,
                                      Model *model_dat,
                                      double *distps,
                                      double *prodps,
                                      double *prodgs)
{
  double *covgp;
  int np, ns, ng, error;

  /* Initializations */

  error = 1;
  covgp = nullptr;

  np = dbdat->getActiveAndDefinedNumber(0);
  ns = dbsrc->getSampleNumber(true);
  ng = dbout->getSampleNumber(true);

  /* Covariance matrix between Mesures and Target */

  covgp = st_calcul_covmat("Covariance G-P", dbout, 0, dbdat, 1, model_dat);
  if (covgp == nullptr) goto label_end;

  /* Add the contribution of the source */

  if (!flag_source)
  {
    for (int ig = 0; ig < ng; ig++)
      for (int ip = 0; ip < np; ip++)
        for (int is = 0; is < ns; is++)
          COVGP(ig,ip) += DISTPS(ip,is) * PRODGS(ig, is);
  }
  else
  {
    for (int ig = 0; ig < ng; ig++)
      for (int ip = 0; ip < np; ip++)
        COVGP(ig,ip) = PRODPS(ip, ig);
  }

  /* Set the error return code */

  error = 0;

  label_end: if (error) covgp = (double*) mem_free((char* ) covgp);
  return (covgp);
}

/****************************************************************************/
/*!
 **  Establish the covariance vector at target
 **
 ** \return  Covariance vector
 **
 ** \param[in]  dbsrc       Db structure containing Sources
 ** \param[in]  dbout       Db structure for target
 ** \param[in]  flag_source If the result is the source, rather than diffusion
 ** \param[in]  model_dat   Model structure for the data
 ** \param[in]  distgs      Distance matrix between Target and Sources
 ** \param[in]  prodgs      Product of DistGS by CovSS
 **
 ** \remarks The returned argument must be freed by the calling function
 **
 *****************************************************************************/
static double* st_inhomogeneous_covgg(Db *dbsrc,
                                      Db *dbout,
                                      int flag_source,
                                      Model *model_dat,
                                      double *distgs,
                                      double *prodgs)
{
  int ng, ns, error;
  double *covgg, c00;
  CovCalcMode mode;

  /* Initializations */

  error = 1;
  covgg = nullptr;

  ns = dbsrc->getSampleNumber(true);
  ng = dbout->getSampleNumber(true);

  /* Core allocation */

  covgg = (double*) mem_alloc(sizeof(double) * ng, 0);
  if (covgg == nullptr) goto label_end;

  /* Calculate the variance term (for a zero-distance) */

  model_calcul_cov(NULL,model_dat, mode, 1, 1., VectorDouble(), &c00);

  /* Calculate the variance vector */

  if (!flag_source)
  {
    for (int ig = 0; ig < ng; ig++)
    {
      covgg[ig] = c00;
      for (int is = 0; is < ns; is++)
        covgg[ig] += DISTGS(ig,is) * PRODGS(ig, is);
    }
  }
  else
  {
    for (int ig = 0; ig < ng; ig++)
      covgg[ig] = c00;
  }

  /* Set the error return code */

  error = 0;

  label_end: if (error) covgg = (double*) mem_free((char* ) covgg);
  return (covgg);
}

/****************************************************************************/
/*!
 **  Calculate auxiliary arrays when drift is preset
 **
 ** \return  Error return code
 **
 ** \param[in]  np          Number of data
 ** \param[in]  nbfl        Number of drift functions
 ** \param[in]  covpp       Inverse Covariance between Data-Data
 ** \param[in]  drftab      Drift matrix at Data
 **
 ** \param[out] yloc        Array: t(F) %*% C^-1
 ** \param[out] zloc        Array: (t(F) %*% C^-1 %*% F)^-1
 **
 ** \remarks The returned arrays 'yloc' and 'zloc' must be freed by the
 ** \remarks calling function
 **
 *****************************************************************************/
static int st_drift_prepar(int np,
                           int nbfl,
                           double *covpp,
                           double *drftab,
                           double **yloc,
                           double **zloc)
{
  double *ymat, *zmat, value;
  int error, ecr;

  /* Initialization */

  error = 1;
  ymat = zmat = nullptr;

  /* First returned array */

  ymat = (double*) mem_alloc(sizeof(double) * nbfl * np, 0);
  if (ymat == nullptr) goto label_end;

  ecr = 0;
  for (int il = 0; il < nbfl; il++)
    for (int ip = 0; ip < np; ip++)
    {
      value = 0.;
      for (int jp = 0; jp < np; jp++)
        value += COVPP(ip,jp) * DRFTAB(jp, il);
      ymat[ecr++] = value;
    }

  /* Second retrned array */

  zmat = (double*) mem_alloc(sizeof(double) * nbfl * nbfl, 0);
  if (zmat == nullptr) goto label_end;

  ecr = 0;
  for (int il = 0; il < nbfl; il++)
    for (int jl = 0; jl < nbfl; jl++)
    {
      value = 0.;
      for (int ip = 0; ip < np; ip++)
        value += YMAT(ip,il) * DRFTAB(ip, jl);
      zmat[ecr++] = value;
    }

  /* Invert 'zmat' */

  if (matrix_invert(zmat, nbfl, -1)) goto label_end;

  /* Set the error return code */

  error = 0;

  label_end: if (error)
  {
    ymat = (double*) mem_free((char* ) ymat);
    zmat = (double*) mem_free((char* ) zmat);
  }
  *yloc = ymat;
  *zloc = zmat;
  return (error);
}

/****************************************************************************/
/*!
 **  Update the weight vector
 **
 ** \param[in]  np          Number of data
 ** \param[in]  nbfl        Number of drift functions
 ** \param[in]  covgp       Covariance matrix between Data and Target
 ** \param[in]  driftg      Drift matrix at Target
 ** \param[in]  ymat        Auxiliary array
 ** \param[in]  zmat        Auxiliary array
 ** \param[in]  maux        Auxiliary array (Dimension: nbfl)
 **
 ** \param[out] lambda      Vector of weights
 ** \param[out] mu          Vector of Lagrange parameters
 **
 *****************************************************************************/
static void st_drift_update(int np,
                            int nbfl,
                            double *covgp,
                            double *driftg,
                            double *ymat,
                            double *zmat,
                            double *maux,
                            double *lambda,
                            double *mu)
{
  double value;

  /* Calculate the Lagrange vector */

  for (int il = 0; il < nbfl; il++)
  {
    value = 0.;
    for (int ip = 0; ip < np; ip++)
      value = YMAT(ip,il) * covgp[ip] - driftg[il];
    maux[il] = value;
  }
  matrix_product(nbfl, nbfl, 1, zmat, maux, mu);

  /* Update the vector of kriging weights */

  for (int ip = 0; ip < np; ip++)
    for (int il = 0; il < nbfl; il++)
      lambda[ip] -= YMAT(ip,il) * mu[il];

  return;
}

/****************************************************************************/
/*!
 **  Inhomogeneous Kriging with Sources
 **
 ** \return  Error return code
 **
 ** \param[in]  dbdat       Db structure containing Data
 ** \param[in]  dbsrc       Db structure containing Sources
 ** \param[in]  dbout       Output Db structure
 ** \param[in]  power       Power of the Distance decay
 ** \param[in]  flag_source If the result is the source, rather than diffusion
 ** \param[in]  model_dat   Model structure for the data
 ** \param[in]  model_src   Model structure for the sources
 **
 *****************************************************************************/
int inhomogeneous_kriging(Db *dbdat,
                          Db *dbsrc,
                          Db *dbout,
                          double power,
                          int flag_source,
                          Model *model_dat,
                          Model *model_src)
{
  int error, np, ip, ns, ng, nvar, neq, nred, nfeq, nbfl, status, nech;
  double *covss, *distps, *distgs, *covpp, *covgp, *covgg, *prodps, *prodgs;
  double *data, *lambda, *driftp, *driftg, *ymat, *zmat, *mu, *maux, *rhs;
  double estim, stdev, auxval;
  NeighWork nbghw;
  VectorInt nbgh_ranks;

  /* Preliminary checks */

  error = nvar = 1;
  NeighUnique* neighU = NeighUnique::create(dbdat->getNDim(),false);
  st_global_init(dbdat, dbout);
  FLAG_EST = 1;
  FLAG_STD = 1;
  distps = distgs = prodgs = prodps = nullptr;
  covss = covpp = covgp = covgg = nullptr;
  lambda = data = driftp = driftg = nullptr;
  ymat = zmat = mu = maux = nullptr;
  if (st_check_environment(1, 1, model_dat, NULL)) goto label_end;

  /* Preliminary checks */

  if (model_dat->getVariableNumber() != nvar)
  {
    messerr("The Model for the Data must be Monovariate");
    goto label_end;
  }
  if (model_src->getVariableNumber() != nvar)
  {
    messerr("The Model for the Sources must be Monovariate");
    goto label_end;
  }

  /* Add the attributes for storing the results */

  if (FLAG_EST != 0)
  {
    IPTR_EST = dbout->addColumnsByConstant(nvar, 0.);
    if (IPTR_EST < 0) goto label_end;
  }
  if (FLAG_STD != 0)
  {
    IPTR_STD = dbout->addColumnsByConstant(nvar, 0.);
    if (IPTR_STD < 0) goto label_end;
  }
  nred = neq = np = dbdat->getActiveAndDefinedNumber(0);
  nfeq = 0;
  ns = dbsrc->getSampleNumber(true);
  ng = dbout->getSampleNumber(true);
  nbfl = model_dat->getDriftNumber();

  /* Core allocation */

  lambda = (double*) mem_alloc(sizeof(double) * np, 0);
  if (lambda == nullptr) goto label_end;
  data = (double*) mem_alloc(sizeof(double) * np, 0);
  if (data == nullptr) goto label_end;

  /* Pre-calculations */

  if (st_model_manage(1, model_dat)) goto label_end;
  if (st_krige_manage(1, nvar, model_dat, neighU)) goto label_end;
  if (krige_koption_manage(1, 1, EKrigOpt::PONCTUAL, 1, VectorInt()))
    goto label_end;
  nbghw.initialize(dbdat, neighU);

  /* Constitute the Data vector */

  for (int iip = ip = 0; iip < dbdat->getSampleNumber(); iip++)
  {
    if (!dbdat->isActiveAndDefined(iip, 0)) continue;
    data[ip] = dbdat->getVariable(iip, 0);
    ip++;
  }

  /* Establish the covariance matrix between Sources */

  covss = st_calcul_covmat("Covarance S_S", dbsrc, 0, dbsrc, 0, model_src);
  if (covss == nullptr) goto label_end;

  /* Establish the distance matrix between Data and Sources */

  distps = st_calcul_distmat("Distance P-S", dbdat, 1, dbsrc, 0, power);
  if (distps == nullptr) goto label_end;

  /* Establish the distance matrix between Target and Sources */

  if (!flag_source)
  {
    distgs = st_calcul_distmat("Distance G-S", dbout, 0, dbsrc, 0, power);
    if (distgs == nullptr) goto label_end;
  }

  /* Establish the Data-Source Product matrix */

  prodps = st_calcul_product("Convolve P-S", np, ns, covss, distps);
  if (prodps == nullptr) goto label_end;

  /* Establish the complete kriging matrix */

  covpp = st_inhomogeneous_covpp(dbdat, dbsrc, model_dat, distps, prodps);
  if (covpp == nullptr) goto label_end;
  if (OptDbg::query(EDbg::KRIGING) || OptDbg::isReferenceDefined())
    krige_lhs_print(np, neq, nred, NULL, covpp);

  /* Invert the Kriging Matrix */

  if (matrix_invert(covpp, np, -1)) goto label_end;

  /* Establish the drift at Data */

  if (nbfl > 0)
  {
    mu = (double*) mem_alloc(sizeof(double) * nbfl, 0);
    if (mu == nullptr) goto label_end;
    maux = (double*) mem_alloc(sizeof(double) * nbfl, 0);
    if (maux == nullptr) goto label_end;

    driftp = st_calcul_drfmat("Drift P", dbdat, 1, model_dat);
    if (driftp == nullptr) goto label_end;

    /* Prepare auxiliary arrays */

    if (st_drift_prepar(np, nbfl, covpp, driftp, &ymat, &zmat)) goto label_end;
  }

  /* Establish the Target-Source Product matrix */

  if (!flag_source)
  {
    prodgs = st_calcul_product("Convolve G-S", ng, ns, covss, distgs);
    if (prodgs == nullptr) goto label_end;
  }

  /* Establish the COVGP */

  covgp = st_inhomogeneous_covgp(dbdat, dbsrc, dbout, flag_source, model_dat,
                                 distps, prodps, prodgs);
  if (covgp == nullptr) goto label_end;

  /* Establish the drift at Target */

  if (nbfl > 0)
  {
    driftg = (double*) mem_alloc(sizeof(double) * nbfl, 0);
    if (driftg == nullptr) goto label_end;
  }

  /* Establish the variance at targets */

  covgg = st_inhomogeneous_covgg(dbsrc, dbout, flag_source, model_dat, distgs,
                                 prodgs);
  if (covgg == nullptr) goto label_end;

  /* Loop on the targets to be processed */

  for (IECH_OUT = 0; IECH_OUT < DBOUT->getSampleNumber(); IECH_OUT++)
  {
    mes_process("Kriging sample", DBOUT->getSampleNumber(), IECH_OUT);
    OptDbg::setIndex(IECH_OUT + 1);
    if (!dbout->isActive(IECH_OUT)) continue;
    if (OptDbg::query(EDbg::KRIGING) || OptDbg::query(EDbg::NBGH) || OptDbg::query(EDbg::RESULTS))
    {
      mestitle(1, "Target location");
      db_sample_print(dbout, IECH_OUT, 1, 0, 0);
    }

    // Neighborhood search

    nbgh_ranks = nbghw.select(DBOUT, IECH_OUT);
    nech = (int) nbgh_ranks.size();
    status = nbgh_ranks.empty();
    rhs = &COVGP(IECH_OUT, 0);

    /* Optional printout of the R.H.S */

    if (OptDbg::force()) krige_rhs_print(nvar, np, neq, nred, NULL, rhs);

    /* Fill the drift at Target point (optional) */

    if (driftp != nullptr)
      model_calcul_drift(model_dat, ECalcMember::LHS, dbout, IECH_OUT, driftg);

    /* Calculate the Kriging weights */

    matrix_product(np, np, 1, covpp, rhs, lambda);
    if (OptDbg::force())
      krige_wgt_print(0, nvar, nvar, nfeq, nbgh_ranks, nred, -1, NULL, lambda);

    /* Update vector of weights in presence of drift */

    if (nbfl > 0)
    {

      /* Evaluate the drift at Target */

      model_calcul_drift(model_dat, ECalcMember::LHS, dbout, IECH_OUT, driftg);

      /* Update the kriging weights */

      st_drift_update(np, nbfl, rhs, driftg, ymat, zmat, maux, lambda, mu);
    }

    /* Perform the estimation */

    matrix_product(1, np, 1, data, lambda, &estim);
    matrix_product(1, np, 1, rhs, lambda, &stdev);

    /* Update the variance in presence of drift */

    if (nbfl > 0)
    {
      matrix_product(1, nbfl, 1, mu, maux, &auxval);
      stdev += auxval;
    }

    /* Update the variance calculation */

    VAR0(0,0) = covgg[IECH_OUT];
    stdev = covgg[IECH_OUT] - stdev;
    stdev = (stdev > 0) ? sqrt(stdev) : 0.;

    /* Store the result */

    dbout->setArray(IECH_OUT, IPTR_EST, estim);
    dbout->setArray(IECH_OUT, IPTR_STD, stdev);

    /* Optional printout */

    if (OptDbg::query(EDbg::KRIGING) || OptDbg::force())
      st_result_kriging_print(0, nvar, 0);
  }

  /* Set the error return flag */

  error = 0;

  label_end: OptDbg::setIndex(0);
  covss = (double*) mem_free((char* ) covss);
  distps = (double*) mem_free((char* ) distps);
  distgs = (double*) mem_free((char* ) distgs);
  prodps = (double*) mem_free((char* ) prodps);
  prodgs = (double*) mem_free((char* ) prodgs);
  driftp = (double*) mem_free((char* ) driftp);
  driftg = (double*) mem_free((char* ) driftg);
  covpp = (double*) mem_free((char* ) covpp);
  covgp = (double*) mem_free((char* ) covgp);
  covgg = (double*) mem_free((char* ) covgg);
  driftp = (double*) mem_free((char* ) driftp);
  driftg = (double*) mem_free((char* ) driftg);
  ymat = (double*) mem_free((char* ) ymat);
  zmat = (double*) mem_free((char* ) zmat);
  maux = (double*) mem_free((char* ) maux);
  mu = (double*) mem_free((char* ) mu);
  data = (double*) mem_free((char* ) data);
  lambda = (double*) mem_free((char* ) lambda);
  (void) st_model_manage(-1, model_dat);
  (void) st_krige_manage(-1, 1, model_dat, neighU);
  (void) krige_koption_manage(-1, 1, EKrigOpt::PONCTUAL, 1, VectorInt());
  delete neighU;
  return (error);
}

/****************************************************************************/
/*!
 **  Kriging (Factorial) a regular grid
 **
 ** \return  Error return code
 **
 ** \param[in]  dbgrid     input and output Db grid structure
 ** \param[in]  model      Model structure
 ** \param[in]  neighparam ANeighParam structure
 ** \param[in]  namconv    Naming Convention
 **
 *****************************************************************************/
int krimage(DbGrid *dbgrid,
            Model *model,
            NeighImage *neighparam,
            const NamingConvention& namconv)
{
  int iptr_est  = 1;
  int nvar = model->getVariableNumber();

  /* Add the attributes for storing the results */

  iptr_est = dbgrid->addColumnsByConstant(nvar, 0.);
  if (iptr_est < 0) return 1;

  /* Setting options */
  // Here ALL options are set, even if most of them could keep their default values

  KrigingSystem ksys(dbgrid, dbgrid, model, neighparam);
  if (ksys.setKrigOptEstim(iptr_est, -1, -1)) return 1;
  if (! ksys.isReady()) return 1;

  /* Loop on the targets to be processed */

  for (int iech_out = 0; iech_out < dbgrid->getSampleNumber(); iech_out++)
  {
    mes_process("Image filtering", dbgrid->getSampleNumber(), iech_out);
     if (ksys.estimate(iech_out)) return 1;
  }

  /* Set the error return flag */

  namconv.setNamesAndLocators(dbgrid, ELoc::Z, nvar, dbgrid, iptr_est, "estim");

  return 0;
}
