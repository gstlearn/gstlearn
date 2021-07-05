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
#include "geoslib_e.h"
#include "Polynomials/Hermite.hpp"
#include "Basic/NamingConvention.hpp"
#include "Basic/Utilities.hpp"

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
#define MATCL(i,j)        (matCL[(i) * nvar + (j)])
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

static double  *covtab,*covaux,*drftab,*fs,*fsf,*d1_1,*d1_2,*var0;
static VectorDouble d1,d1_t;
static double  *lhs,*lhs_b,*rhs,*wgt,*zext,*zam1,*ff0,*varb;
static int     *rank,*flag;
static int      IECH_NBGH  = -1;
static int      KRIGE_INIT =  0;
static int      MODEL_INIT =  0;
static int      IECH_OUT   = -1;
static int      RAND_INDEX = -1;
static int      FLAG_EST,FLAG_STD,FLAG_WGT,FLAG_COLK,FLAG_SIMU,FLAG_LTERM;
static int      FLAG_BAYES,FLAG_PROF,FLAG_VARZ,FLAG_DGM;
static int      IPTR_EST,IPTR_STD,IPTR_VARZ,IPTR_NBGH;
static int     *RANK_COLCOK;
static Db      *DBIN,*DBOUT;
static Koption *KOPTION;
static double   R_COEFF;
static int      INH_FLAG_VERBOSE = 0;
static int      INH_FLAG_LIMIT   = 1;
static double   EPS = 1.e-5;
static char     string[100];
static int      NDIM_EXTERNAL      =  0;
static int      COV_EXTERNAL_DB_1  =  0;
static int      COV_EXTERNAL_DB_2  =  0;
static int      COV_EXTERNAL_ECH_1 = -1;
static int      COV_EXTERNAL_ECH_2 = -1;

typedef struct {
  int     ndtot;
  int     rank1;
  int     rank2;
  Model  *model;
  int     nugget_opt;
  int     nostd;
  int     member;
  int     icov_r;
  double  weight;
} Disc_Structure;

/****************************************************************************/
/*!
**  Identify the indices of the Db and the sample ranks
**  (used by the external covariance function in case of Kriging)
**
** \param[out] E_Cov    The External_Cov structure
**
** \remarks The arguments 'rank_db1', 'rank_ech1', 'rank_db2' and 'rank_ech2
** \remarks are assigned
**
*****************************************************************************/
GEOSLIB_API void fill_external_cov_kriging(External_Cov& E_Cov)
{
  E_Cov.ndim      = NDIM_EXTERNAL;
  E_Cov.rank_db1  = COV_EXTERNAL_DB_1;
  E_Cov.rank_db2  = COV_EXTERNAL_DB_2;
  E_Cov.rank_ech1 = COV_EXTERNAL_ECH_1;
  E_Cov.rank_ech2 = COV_EXTERNAL_ECH_2;
}

/****************************************************************************/
/*!
**  Local function checking if the continuous Kriging as been required.
**  Then the LHS must be updated for each target
**
** \return  1 if the Continuous Kriging is required; 0 otherwise
**
** \param[in] neigh   Neigh structure
**
*****************************************************************************/
static int flag_continuous_kriging(Neigh *neigh)
{
  return(neigh->getFlagContinuous());
}

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
static double *st_core(int nli,
                       int nco)
{
  double *tab,rsize;
  int size,i;

  /* Initialization */

  tab  = (double *) NULL;
  rsize = (double) nli * (double) nco;
  if (rsize < 0 || rsize > INT_MAX) 
  {
    messerr("Core allocation problem: Size (%d x %d) too big",nli,nco);
    return(tab);
  }
  size = nli * nco;

  /* Allocation */

  tab = (double *) mem_alloc(sizeof(double) * size,0);
  if (tab == (double *) NULL) 
  {
    messerr("Core allocation problem: Size (%d) too big",size);
    return(tab);
  }

  for (i=0; i<size; i++) tab[i] = 0.;
  return(tab);
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
static int *st_icore(int nli,
                     int nco)
{
  int *tab,size,i;

  /* Initialization */

  tab  = (int *) NULL;
  size = nli * nco;

  /* Allocation */

  tab = (int *) mem_alloc(sizeof(int) * size,0);
  if (tab != (int *) NULL) for (i=0; i<size; i++) tab[i] = 0;

  return(tab);
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
static int *st_relative_position_array(int  mode,
                                       int  neq,
                                       int *rel_arg)
{
  int *rel,i,j;

  /* Dispatch */

  if (mode > 0)
  {
    
    /* Creation */

    rel = (int *) st_icore(neq,1);
    if (rel == (int *) NULL) return(rel);
    for (i=j=0; i < neq; i++)
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
    rel = (int *) mem_free((char *) rel);
  }
  return(rel);
}

/****************************************************************************/
/*!
**  Initialize the static global variables
**
** \param[in]  dbin   input Db structure
** \param[in]  dbout  output Db structure
**
*****************************************************************************/
static void st_global_init(Db *dbin,
                           Db *dbout)
{
  FLAG_STD  = FLAG_EST   = FLAG_WGT   = FLAG_LTERM = FLAG_VARZ = 0;
  FLAG_COLK = FLAG_BAYES = FLAG_PROF  = FLAG_SIMU  = FLAG_DGM = 0;
  IPTR_EST  = IPTR_STD   = IPTR_VARZ  = IPTR_NBGH  = 0;
  IECH_OUT  = 0;
  IECH_NBGH = -1;

  /* Set the global variables */

  DBIN  = dbin;
  DBOUT = dbout;

  /* Change of support coefficient for DGM */

  R_COEFF = 1.;

  return;
}

/****************************************************************************/
/*!
**  Returns the coordinate of the data (at rank if rank >= 0)
**  or of the target (at IECH_OUT if rank < 0)
**
** \param[in]  rank   Rank of the sample
** \param[in]  idim   Rank of the coordinate
**
*****************************************************************************/
static double st_get_idim(int rank,
                          int idim)
{
  double value;

  if (rank >= 0)
  {
    value = get_IDIM(DBIN,rank,idim);
  }
  else
  {
    value = get_IDIM(DBOUT,IECH_OUT,idim);
  }
  return(value);
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
** \param[in]  member       Member of the Kriging System (::ENUM_MEMBERS)
** \param[in]  icov_r       rank of the target covariance or -1 for all
** \param[in]  weight       Weight attached to this calculation
** \param[in]  rank1        Rank of the first sample
** \param[in]  rank2        Rank of the second sample
**
** \param[out] d1           Working array
** \param[out] covtab       Output covariance array
**
*****************************************************************************/
static void st_cov(Model  *model,
                   int     flag_init,
                   int     nugget_opt,
                   int     nostd,
                   int     member,
                   int     icov_r,
                   double  weight,
                   int     rank1,
                   int     rank2,
                   VectorDouble d1,
                   double *covtab)
{
  Db *db1,*db2;
  int iech1,iech2;

  /* Initializations */
  
  if (rank1 >= 0)
  {
    db1   = DBIN;
    iech1 = rank1;
    COV_EXTERNAL_DB_1  = 1;
    COV_EXTERNAL_ECH_1 = iech1;
  }
  else
  {
    db1   = DBOUT;
    iech1 = IECH_OUT;
    COV_EXTERNAL_DB_1  = 2;
    COV_EXTERNAL_ECH_1 = iech1;
  }

  if (rank2 >= 0)
  {
    db2   = DBIN;
    iech2 = rank2;
    COV_EXTERNAL_DB_2  = 1;
    COV_EXTERNAL_ECH_2 = iech2;
  }
  else
  {
    db2   = DBOUT;
    iech2 = IECH_OUT;
    COV_EXTERNAL_DB_2  = 2;
    COV_EXTERNAL_ECH_2 = iech2;
  }

  CovCalcMode mode;
  mode.update(nugget_opt,nostd,member,icov_r,0,1);
  model_calcul_cov_nostat(model,mode,flag_init,weight,
                          db1,iech1,db2,iech2,d1,covtab);
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
static void st_data_discretize_dd(int idim,
                                  int jdim,
                                  Disc_Structure *it)
{
  double exts2,dsize,decal;

  // Initialization

  if (idim < it->model->getDimensionNumber()-1)
  {
    idim = idim + 1;

    // Loop in the current dimension

    exts2 = DBIN->getBlockExtension(it->rank1,idim) / 2.;
    dsize = KOPTION->dsize[idim];

    if (exts2 <= 0. || dsize <= 0.)
    {

      /* Ponctual support */

      d1_1[idim] = 0.;
      st_data_discretize_dd(idim,jdim,it);
    }
    else
    {
    
      /* Implicit loop until reaching the edge of the data */

      decal = -exts2 + dsize / 2.;
      do {
        d1_1[idim] = decal;
        st_data_discretize_dd(idim,jdim,it);
        decal = decal + dsize;
      } while (decal < exts2);
    }
  }
  else if (jdim < it->model->getDimensionNumber()-1)
  {
    jdim = jdim + 1;

    // Loop in the current dimension

    exts2 = DBIN->getBlockExtension(it->rank2,jdim) / 2.;
    dsize = KOPTION->dsize[jdim];
    
    if (exts2 <= 0 || dsize <= 0.)
    {
      /* Punctual support */

      d1_2[jdim] = 0.;
      st_data_discretize_dd(idim,jdim,it);
    }
    else
    {

      /* Implicit loop until reaching the edge of the data */

      decal = -exts2 + dsize / 2.;
      do {
        d1_2[jdim] = decal;
        st_data_discretize_dd(idim,jdim,it);
        decal = decal + dsize;
      } while (decal < exts2);
    }
  }
  else
  {
    
    // End of implicit loop on dimensions

    it->ndtot++;
    for (int i=0; i<it->model->getDimensionNumber(); i++)
      d1_t[i] = d1[i] + d1_1[i] + d1_2[i];
    st_cov(it->model,0,it->nugget_opt,it->nostd,it->member,
           it->icov_r,it->weight,it->rank1,it->rank2,
           d1_t,covaux);
  }
}

GEOSLIB_API int is_flag_data_disc_defined(void)
{
  return KOPTION->flag_data_disc;
}

/****************************************************************************/
/*!
**  Calculate the covariance between two Data samples
**  This covariance takes possible Data discretization into account
**
** \param[in]  model  Model structure
** \param[in]  flag_init    Initialize the array beforehand
** \param[in]  nugget_opt   Option for the nugget effect basic structure
** \li                       0 : no particular option
** \li                       1 : discard the nugget effect
** \li                      -1 : only consider the nugget effect
** \param[in]  nostd        0 standard; +-1 special; ITEST normalized
** \param[in]  member       Member of the Kriging System (::ENUM_MEMBERS)
** \param[in]  icov_r       rank of the target covariance or -1 for all
** \param[in]  weight       Weight attached to this calculation
** \param[in]  rank1        Rank of the first sample
** \param[in]  rank2        Rank of the second sample
**
** \param[out] d1           Working array
** \param[out] covtab       Output covariance array
**
*****************************************************************************/
static void st_cov_dd(Model  *model,
                      int     flag_init,
                      int     nugget_opt,
                      int     nostd,
                      int     member,
                      int     icov_r,
                      double  weight,
                      int     rank1,
                      int     rank2,
                      VectorDouble d1,
                      double *covtab)
{
  int    nvar_m;
  double scale;
  Disc_Structure int_disc;

  if (! KOPTION->flag_data_disc)
  {

    // Data is considered as ponctual

    st_cov(model,0,nugget_opt,nostd,MEMBER_LHS,
           icov_r,weight,rank1,rank2,d1,covtab);
  }
  else
  {
    // Data have a support

    nvar_m = model->getVariableNumber();
    model_covtab_init(1,model,covaux);
    
    // Implicit discretization loop for first datum

    int_disc.rank1      = rank1;
    int_disc.rank2      = rank2;
    int_disc.model      = model;
    int_disc.nugget_opt = nugget_opt;
    int_disc.nostd      = nostd;
    int_disc.member     = MEMBER_LHS;
    int_disc.icov_r     = icov_r;
    int_disc.weight     = weight;
    int_disc.ndtot      = 0;

    st_data_discretize_dd(-1,-1,&int_disc);

    // Normalization

    scale = (double) int_disc.ndtot;
    for (int ivar_m=0; ivar_m<nvar_m; ivar_m++)
      for (int jvar_m=0; jvar_m<nvar_m; jvar_m++)
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
static void st_data_discretize_dg(int idim,
                                  Disc_Structure *it)
{
  double exts2,dsize,decal;

  // Initialization

  if (idim < it->model->getDimensionNumber()-1)
  {
    idim = idim + 1;

    // Loop in the current dimension

    exts2 = DBIN->getBlockExtension(it->rank1,idim) / 2.;
    dsize = KOPTION->dsize[idim];

    if (exts2 <= 0. || dsize <= 0.) 
    {

      /* Ponctual support */

      d1_1[idim] = 0.;
      st_data_discretize_dg(idim,it);
    }
    else
    {

      /* Implicit loop until reaching the edge of the data */
      
      decal = -exts2 + dsize / 2.;
      do {
        
        d1_1[idim] = decal;
        st_data_discretize_dg(idim,it);
        decal = decal + dsize;
      } while (decal < exts2);
    }
  }
  else
  {
    
    // End of implicit loop on dimensions

    it->ndtot++;
    for (int i=0; i<it->model->getDimensionNumber(); i++) d1_t[i] = d1[i] + d1_1[i];
    st_cov(it->model,0,it->nugget_opt,it->nostd,it->member,
           it->icov_r,it->weight,it->rank1,it->rank2,
           d1_t,covaux);
  }
}

/****************************************************************************/
/*!
**  Calculate the covariance between a Data sample and a Target sample
**  This covariance takes possible Data discretization into account
**
** \param[in]  model        Model structure
** \param[in]  flag_init    Initialize the array beforehand
** \param[in]  nugget_opt   Option for the nugget effect basic structure
** \li                       0 : no particular option
** \li                       1 : discard the nugget effect
** \li                      -1 : only consider the nugget effect
** \param[in]  nostd        0 standard; +-1 special; ITEST normalized
** \param[in]  member       Member of the Kriging System (::ENUM_MEMBERS)
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
static void st_cov_dg(Model  *model,
                      int     flag_init,
                      int     nugget_opt,
                      int     nostd,
                      int     member,
                      int     icov_r,
                      double  weight,
                      int     rank1,
                      int     rank2,
                      VectorDouble d1,
                      double *covtab)
{
  int    nvar_m;
  double scale;
  Disc_Structure int_disc;

  if (! KOPTION->flag_data_disc)
  {

    // Data is considered as ponctual

    st_cov(model,0,nugget_opt,nostd,member,
           icov_r,weight,rank1,rank2,d1,covtab);
  }
  else
  {
    // Data have a support

    nvar_m = model->getVariableNumber();
    model_covtab_init(1,model,covaux);
    
    // Implicit discretization loop on the data

    int_disc.rank1      = rank1;
    int_disc.rank2      = rank2;
    int_disc.model      = model;
    int_disc.nugget_opt = nugget_opt;
    int_disc.nostd      = nostd;
    int_disc.member     = member;
    int_disc.icov_r     = icov_r;
    int_disc.weight     = weight;
    int_disc.ndtot      = 0;

    for (int i=0; i<model->getDimensionNumber(); i++) d1_1[i] = d1[i];
    st_data_discretize_dg(-1,&int_disc);

    // Normalization

    scale = (double) int_disc.ndtot;
    for (int ivar_m=0; ivar_m<nvar_m; ivar_m++)
      for (int jvar_m=0; jvar_m<nvar_m; jvar_m++)
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
static double st_get_ivar(int rank,
                          int ivar)
{
  double value;
  int    jvar;

  if (rank >= 0)
  {

    // Variable in the Input file

    if (! FLAG_SIMU)

      // Particular case of simulations

      value = DBIN->getVariable(rank,ivar);
    else

      // Case of the traditional kriging based on Z-variables

      value = DBIN->getSimvar(LOC_SIMU,rank,0,ivar,0,1,0);
  }
  else
  {

    // Variable in the Output file: colocated case

    jvar = RANK_COLCOK[ivar];
    if (jvar < 0)
      value = TEST;
    else
      value = get_ARRAY(DBOUT,IECH_OUT,jvar);
  }

  return(value);
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
static double st_get_verr(int rank,
                          int ivar)
{
  double value;

  if (rank >= 0)
  {
    value = DBIN->getVarianceError(rank,ivar);
  }
  else
  {
    value = DBOUT->getVarianceError(IECH_OUT,ivar);
  }
  return(value);
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
static double st_get_fext(int rank,
                          int ibfl)
{
  double value;

  if (rank >= 0)
  {
    value = DBIN->getExternalDrift(rank,ibfl);
  }
  else
  {
    value = DBOUT->getExternalDrift(IECH_OUT,ibfl);
  }
  return(value);
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
    value = DBIN->getSimvar(LOC_SIMU, rank, isimu, ivar, icase, nbsimu, nvar);
  else
    value = DBOUT->getSimvar(LOC_SIMU, IECH_OUT, isimu, ivar, icase, nbsimu,
                             nvar);
  return(value);
}

/****************************************************************************/
/*!
**  Checks the kriging environment
**
** \return  Error return code
**
** \param[in]  flag_in   1 if the Input Db is used
** \param[in]  flag_out  1 if the Output Db is used
** \param[in]  model     Model structure (optional)
** \param[in]  neigh     Neigh structure (optional)
**
** \remarks The address of the argument 'neigh' is memorized in a local
** \remarks static variable
**
*****************************************************************************/
static int st_check_environment(int    flag_in,
                                int    flag_out,
                                Model *model,
                                Neigh *neigh)
{
  double *dbin_mini,*dbin_maxi,*dbout_mini,*dbout_maxi;
  int     error,ndim,nvar,nfex,nparam;

  /* Initializations */

  error = 1;
  nvar = ndim = nfex = 0;
  dbin_mini = dbin_maxi = dbout_mini = dbout_maxi = (double *) NULL;

  /*********************************/
  /* Compatibility between two Dbs */
  /*********************************/

  ndim = 0;
  if (flag_in  && ndim == 0) ndim = DBIN->getNDim();
  if (flag_out && ndim == 0) ndim = DBOUT->getNDim();
  if (flag_in && flag_out && ! DBIN->hasSameDimension(DBOUT)) goto label_end;

  /**********************/
  /* Checking the model */
  /**********************/

  if (model != (Model *) NULL)
  {
    nvar = model->getVariableNumber();
    if (nvar <= 0)
    {
      messerr("The number of variables must be positive = %d",model->getVariableNumber());
      goto label_end;
    }
    // The following test is avoided in the case of simulations
    // as there may be no Z-variable defined as this stage (Gibbs)
    if (flag_in && ! FLAG_SIMU && DBIN->getVariableNumber() != nvar)
    {
      messerr("The number of variables of the Data (%d)",DBIN->getVariableNumber());
      messerr("does not match the number of variables of the Model (%d)",
              nvar);
      goto label_end;
    }
    if (model->getCovaNumber() <= 0)
    {
      messerr("The number of covariance must be positive");
      goto label_end;
    }
    if (model->getDimensionNumber() <= 0)
    {
      messerr("The Space Dimension must be positive = %d",model->getDimensionNumber());
      goto label_end;
    }
    if (model->getDimensionNumber() != ndim)
    {
      messerr("The Space Dimension of the Db structure (%d)",ndim);
      messerr("Does not correspond to the Space Dimension of the model (%d)",
              model->getDimensionNumber());
      goto label_end;
    }
    if (model->isNoStat())
    {
      nparam = model->getNoStat().getNoStatElemNumber();
      
      if (flag_out && nparam != get_LOCATOR_NITEM(DBOUT,LOC_NOSTAT))
      {
        messerr("The number of non-stationary parameters in the Model (%d)",
                nparam);
        messerr("is not equal to the number of Non-stationary locators in the Output Db (%d)",
                get_LOCATOR_NITEM(DBOUT,LOC_NOSTAT));
        goto label_end;
      }
      if (flag_in && nparam != get_LOCATOR_NITEM(DBIN,LOC_NOSTAT))
      {
        if (! (flag_out && is_grid(DBOUT)))
        {
          messerr("The number of non-stationary parameters in the Model (%d)",
                  nparam);
          messerr("is not equal to the number of Non-stationary locators in the Input Db (%d)",
                  get_LOCATOR_NITEM(DBIN,LOC_NOSTAT));
          goto label_end;
        }
      }
    }

    // External drifts
    nfex = model_nfex(model);
    if (nfex > 0)
    {
      if (flag_out && DBOUT->getExternalDriftNumber() != nfex)
      {
        messerr("The Model requires %d external drift(s)",
                model_nfex(model));
        messerr("but the output Db refers to %d external drift variables",
                DBOUT->getExternalDriftNumber());
        goto label_end;
      }
        
      if (flag_in && DBIN->getExternalDriftNumber() != nfex)
      {
        if (! (flag_out && is_grid(DBOUT)))
        {
          messerr("The Model requires %d external drift(s)",
                  model_nfex(model));
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

  if (model != (Model *) NULL)
  {

    /* Input Db structure */

    if (flag_in)
    {
      dbin_mini = db_sample_alloc(DBIN,LOC_X);
      if (dbin_mini == (double *) NULL) goto label_end;
      dbin_maxi = db_sample_alloc(DBIN,LOC_X);
      if (dbin_maxi == (double *) NULL) goto label_end;
      if (db_extension(DBIN,
                       dbin_mini,dbin_maxi,(double *) NULL)) goto label_end;
    }

    /* Output Db structure */

    if (flag_out)
    {
      dbout_mini = db_sample_alloc(DBOUT,LOC_X);
      if (dbout_mini == (double *) NULL) goto label_end;
      dbout_maxi = db_sample_alloc(DBOUT,LOC_X);
      if (dbout_maxi == (double *) NULL) goto label_end;
      if (db_extension(DBOUT,
                       dbout_mini,dbout_maxi,(double *) NULL)) goto label_end;
    }

    if (flag_in && flag_out)
      model->setField(ut_merge_extension(ndim,
                                        dbin_mini,dbin_maxi,
                                        dbout_mini,dbout_maxi));

    dbin_mini  = db_sample_free(dbin_mini);
    dbin_maxi  = db_sample_free(dbin_maxi);
    dbout_mini = db_sample_free(dbout_mini);
    dbout_maxi = db_sample_free(dbout_maxi);
  }

  /*****************************/
  /* Checking the Neighborhood */
  /*****************************/

  if (neigh != (Neigh *) NULL)
  {
    if (neigh->getNDim() != ndim)
    {
      messerr("The Space Dimension of the Neighborhood (%d)",neigh->getNDim());
      messerr("does not correspond to the Space Dimension of the first Db (%d)",
              ndim);
      goto label_end;
    }
    if (neigh->getType() == NEIGH_IMAGE && (! flag_out || ! is_grid(DBOUT)))
    {
      messerr("The Image neighborhood can only be used when the output Db is a grid");
      goto label_end;
    }
  }

  /* Set the error return code */

  error = 0;

label_end:
  return(error);
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
static int st_model_manage(int    mode,
                           Model *model)

{
  int nvar,nbfl;

  /* Initializations */

  nvar  = model->getVariableNumber();
  nbfl  = model->getDriftNumber();

  /* Dispatch */

  if (mode == 1)
  {

    /* Allocation */

    if (MODEL_INIT) return(1);
    d1.resize(DBIN->getNDim());
    d1_1 = db_sample_alloc(DBIN,LOC_X);
    if (d1_1 == (double *) NULL) return(1);
    d1_2 = db_sample_alloc(DBIN,LOC_X);
    if (d1_2 == (double *) NULL) return(1);
    d1_t.resize(DBIN->getNDim());
    covtab = st_core(nvar,nvar);
    if (covtab == (double *) NULL) return(1);
    covaux = st_core(nvar,nvar);
    if (covaux == (double *) NULL) return(1);
    if (nbfl > 0)
    {
      drftab = st_core(nbfl,1);
      if (drftab == (double *) NULL) return(1);
    }
    MODEL_INIT = 1;
  }
  else
  {
    if (! MODEL_INIT) return(1);
    d1_1   = db_sample_free(d1_1);
    d1_2   = db_sample_free(d1_2);
    covtab = (double *) mem_free((char *) covtab);
    covaux = (double *) mem_free((char *) covaux);
    drftab = (double *) mem_free((char *) drftab);
    MODEL_INIT = 0;
  }
  return(0);
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
  int neqmax,ncmax;

  /* Initializations */

  ncmax   = nmax * nvar;
  neqmax  = ncmax + nfeq;
  if (FLAG_COLK) neqmax += nvar;
  if (FLAG_COLK) nech   += 1;
  IECH_NBGH = -1;

  /* Dispatch */

  if (mode == 1)
  {

    /* Allocation */

    if (KRIGE_INIT) return(1);
    rank = st_icore(nech,1);
    if (rank == (int *) NULL) return(1);
    flag = st_icore(neqmax,1);
    if (flag == (int *) NULL) return(1);
    lhs = st_core(neqmax,neqmax);
    if (lhs == (double *) NULL) return(1);
    lhs_b = st_core(neqmax,neqmax);
    if (lhs_b == (double *) NULL) return(1);
    rhs = st_core(neqmax,nvar);
    if (rhs == (double *) NULL) return(1);
    zext = st_core(neqmax,1);
    if (zext == (double *) NULL) return(1);
    zam1 = st_core(neqmax,1);
    if (zam1 == (double *) NULL) return(1);
    wgt = st_core(neqmax,nvar);
    if (wgt == (double *) NULL) return(1);
    var0 = st_core(nvar,nvar);
    if (var0 == (double *) NULL) return(1);
    if (FLAG_BAYES)
    {
      fsf = st_core(ncmax,ncmax);
      if (fsf == (double *) NULL) return(1);
      varb   = st_core(nvar,nvar);
      if (varb   == (double *) NULL) return(1);
      if (nfeq > 0)
      {
        fs = st_core(ncmax,nfeq);
        if (fs == (double *) NULL) return(1);
        ff0   = st_core(nfeq,nvar);
        if (ff0   == (double *) NULL) return(1);
      }
    }
    KRIGE_INIT = 1;
  }
  else
  {

    /* Deallocation */

    if (! KRIGE_INIT) return(1);
    rank   = (int    *) mem_free((char *) rank);
    flag   = (int    *) mem_free((char *) flag);
    lhs    = (double *) mem_free((char *) lhs);
    lhs    = (double *) mem_free((char *) lhs);
    lhs_b  = (double *) mem_free((char *) lhs_b);
    rhs    = (double *) mem_free((char *) rhs);
    zext   = (double *) mem_free((char *) zext);
    zam1   = (double *) mem_free((char *) zam1);
    wgt    = (double *) mem_free((char *) wgt);
    var0   = (double *) mem_free((char *) var0);
    if (FLAG_BAYES)
    {
      fs     = (double *) mem_free((char *) fs);
      fsf    = (double *) mem_free((char *) fsf);
      ff0    = (double *) mem_free((char *) ff0);
    }
    KRIGE_INIT = 0;
  }

  return(0);
}

/****************************************************************************/
/*!
**  Calculate the maximum number of samples for Bench Neighborhood
**
** \return  Maximum number of samples per neighborhood
**
** \param[in]  neigh  Neigh structure
**
*****************************************************************************/
static int st_bench_nmax(Neigh *neigh)

{
  int     nech,nmax,iech,jech,nloc;
  double *tab;

  /* Initializations */

  tab   = (double *) NULL;
  nech  = get_NECH(DBIN);
  nmax  = nech;
  if (DBIN->getNDim() <= 2) return(nmax);

  /* Core allocation */
  tab = db_vector_alloc(DBIN);
  if (tab == (double *) NULL) goto label_end;

  /* Read the vector of the third coordinates */
  if (db_coorvec_get(DBIN,DBIN->getNDim()-1,tab)) goto label_end;

  /* Sort the third coordinate vector */
  ut_sort_double(0,nech,NULL,tab);

  /* Loop on the first point */
  nmax = 0;
  for (iech=0; iech<nech-1; iech++)
  {

    /* Loop on the second point */
    nloc = 1;
    for (jech=iech+1; jech<nech; jech++)
    {
      if (ABS(tab[jech] - tab[iech]) > 2. * neigh->getWidth()) break;
      nloc++;
    }

    /* Store the maximum number of samples */
    if (nloc > nmax) nmax = nloc;
  }

  if (debug_query("db"))
  {
    message("Statistics on Bench neighborhood search:\n");
    message("- Vertical Tolerance = %lf\n",neigh->getWidth());
    message("- Maximum number of samples = %d\n",nmax);
  }

label_end:
  tab = db_vector_free(tab);
  return(nmax);
}

/****************************************************************************/
/*!
**  Returns the maximum number of points per neighborhood
**
** \return  Maximum number of points per neighborhood
**
** \param[in]  neigh  Neigh structure
**
*****************************************************************************/
static int st_get_nmax(Neigh *neigh)

{
  int nmax; 

  nmax    = get_NECH(DBIN);
  if (neigh->getType() == NEIGH_MOVING)
    nmax = (neigh->getFlagSector()) ?
        neigh->getNSect() * neigh->getNSMax() : neigh->getNMaxi();
  else if (neigh->getType() == NEIGH_BENCH)
    nmax = st_bench_nmax(neigh);
  else if (neigh->getType() == NEIGH_UNIQUE)
    nmax = DBIN->getActiveSampleNumber();
  else if (neigh->getType() == NEIGH_IMAGE)
  {
    nmax = 1;
    for (int idim=0; idim<neigh->getNDim(); idim++)
      nmax *= (2. * neigh->getImageRadius(idim) + 1);
  }
  return nmax;
}

/****************************************************************************/
/*!
**  Management of internal arrays used by kriging procedure
**
** \return  Error return code
**
** \param[in]  mode   1 for allocation; -1 for deallocation
** \param[in]  nvar   Number of variables to be calculated
** \param[in]  model  Model structure
** \param[in]  neigh  Neigh structure
**
** \remarks  The number of variables corresponds to the number of variables
** \remarks  to be calculated. It is not necessarily equal to the number of
** \remarks  variables contained in Model (when kriging a linear combination
** \remarks  of variables for example): hence the use of the 'nvar' passed
** \remarks  as an argument
**
*****************************************************************************/
static int st_krige_manage(int    mode,
                           int    nvar,
                           Model *model,
                           Neigh *neigh)
{
  int nech,nfeq,nmax;

  /* Initializations */

  nvar    = model->getVariableNumber();
  nfeq    = model->getDriftEquationNumber();
  nech    = get_NECH(DBIN);
  nmax    = st_get_nmax(neigh);
  IECH_NBGH = -1;

  return(st_krige_manage_basic(mode,nech,nmax,nvar,nfeq));
}

/****************************************************************************/
/*!
**  Allocate the Target discretization 
**
** \param[in]  ndim       Space dimension
** \param[in]  ndisc      Discretization parameters (or NULL)
**
*****************************************************************************/
static int st_block_discretize_alloc(int  ndim,
                                     VectorInt ndisc)
{
  int ntot;

  ntot = 1;
  for (int idim=0; idim<ndim; idim++) ntot *= ndisc[idim];
  if (ntot <= 0) return(1);
  KOPTION->ntot = ntot;

  KOPTION->ndisc = st_icore(ndim,1);
  if (KOPTION->ndisc == (int *) NULL) return(1);
  KOPTION->disc1 = st_core(ndim,ntot);
  if (KOPTION->disc1 == (double *) NULL) return(1);
  KOPTION->disc2 = st_core(ndim,ntot);
  if (KOPTION->disc2 == (double *) NULL) return(1);
  for (int idim=0; idim<ndim; idim++) KOPTION->ndisc[idim] = ndisc[idim];
  return(0);
}

/****************************************************************************/
/*!
**  Allocate the Data discretization 
** 
** \returns Error return code
**
** \param[in] ndim    Space dimension
**
*****************************************************************************/
static void st_data_discretize_alloc(int ndim)

{
  int nrow,ncol;

  KOPTION->flag_data_disc = 0;
  if (DBIN->getBlockExtensionNumber() > 0)
  {
    if (! get_keypair("Data_Discretization",&nrow,&ncol,&KOPTION->dsize))
    {
      if (nrow * ncol != ndim)
      {
        messerr("Data discretization is defined using set_keypair mechanism");
        messerr("with keyword 'Data_Discretization'");
        messerr("But its dimension should be %d (instead of %d x %d)",
                ndim,nrow,ncol);
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
static void st_block_discretize(int mode,
                                int flag_rand,
                                int iech)
{
  int    i,j,jech,ntot,nd,nval,ndim,idim,memo;
  double taille;

  /* Initializations */

  memo = law_get_random_seed();
  ntot = KOPTION->ntot;
  ndim = KOPTION->ndim;
  law_set_random_seed(1234546);

  /* Loop on the discretization points */

  for (i=0; i<ntot; i++)
  {
    jech = i;
    nval = ntot;
    for (idim=ndim-1; idim>=0; idim--)
    {
      taille = (mode == 0) ? DBOUT->getDX(idim) : DBOUT->getBlockExtension(iech,idim);
      nd     = KOPTION->ndisc[idim];
      nval  /= nd;
      j      = jech / nval;
      jech  -= j * nval;
      DISC1(i,idim) = taille * ((j + 0.5) / nd - 0.5);
      DISC2(i,idim) = DISC1(i,idim);
      if (flag_rand) DISC2(i,idim) += 
        taille * law_uniform(-0.5,0.5) / (double) nd;
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
** \param[in]  calcul      Type of calculation (::ENUM_KOPTIONS)
** \param[in]  flag_rand   0 if the second discretization is regular
**                         1 if the second point must be randomized
** \param[in]  ndisc       Discretization parameters (or NULL)
**
** \remark  This function manages the global structure KOPTION
**
*****************************************************************************/
GEOSLIB_API int krige_koption_manage(int       mode,
                                     int       flag_check,
                                     int       calcul,
                                     int       flag_rand,
                                     VectorInt ndisc)
{
  int ndim,error;

  /* Initializations */

  error   = 1;
  ndim    = DBOUT->getNDim();

  /* Dispatch */

  if (mode == 1)
  {

    /* Allocation of the structure */

    KOPTION = (Koption *) mem_alloc(sizeof(Koption),0);
    if (KOPTION == (Koption *) NULL) return(1);
    KOPTION->calcul = calcul;

    // Target discretization
    KOPTION->ndim   = ndim;
    KOPTION->ntot   = 0;
    KOPTION->disc1  = (double *) NULL;
    KOPTION->disc2  = (double *) NULL;
    KOPTION->ndisc  = (int    *) NULL;

    // Data discretization
    KOPTION->flag_data_disc = 0;
    KOPTION->dsize = (double *) NULL;

    /* Data discretization case (optional) */

    st_data_discretize_alloc(ndim);

    /* Block discretization case */

    switch (KOPTION->calcul)
    {
      case KOPTION_PONCTUAL:
      case KOPTION_DRIFT:
        break;

      case KOPTION_BLOCK:

        /* Preliminary checks */

        if (flag_check && ! is_grid(DBOUT))
        {
          messerr("Discretization is not allowed if the Target is not a Grid");
          goto label_dealloc;
        }
        if (ndisc.empty())
        {
          messerr("For block estimation, Discretization must be provided");
          goto label_dealloc;
        }

        if (st_block_discretize_alloc(ndim,ndisc)) goto label_dealloc;

        st_block_discretize(0,flag_rand,0);

        break;
    }
    error = 0;
  }
  else
  {
    error    = 0;

    /* Deallocation procedure */

  label_dealloc:
    if (KOPTION != (Koption *) NULL)
    {
      KOPTION->ndisc = (int     *) mem_free((char *) KOPTION->ndisc);
      KOPTION->disc1 = (double  *) mem_free((char *) KOPTION->disc1);
      KOPTION->disc2 = (double  *) mem_free((char *) KOPTION->disc2);
      KOPTION->dsize = (double  *) mem_free((char *) KOPTION->dsize);
      KOPTION        = (Koption *) mem_free((char *) KOPTION);
    }
  }

  return(error);
}

/****************************************************************************/
/*!
**  Define the array flag to convert from isotropic to heterotopic case
**
** \param[in]  model  Model structure
** \param[in]  nech   Number of active samples
** \param[in]  neq    Number of equations
**
** \param[out]  nred  Reduced number of equations
**
*****************************************************************************/
static void st_flag_define(Model *model,
                           int    nech,
                           int    neq,
                           int   *nred)
{
  int i,iech,ibfl,ib,il,ivar,idim,valid,count,nvar;

  /* Initializations */

  *nred = ivar = 0;
  nvar = model->getVariableNumber();
  for (i=0; i<neq; i++) flag[i] = 1;

  /* Check on the coordinates */

  for (iech=0; iech<nech; iech++)
  {
    valid = 1;
    for (idim=0; idim<DBIN->getNDim(); idim++)
      if (FFFF(st_get_idim(rank[iech],idim))) valid = 0;
    if (! valid)
      for (ivar=0; ivar<DBIN->getVariableNumber(); ivar++) FLAG(iech,ivar) = 0;
  }

  /* Check on the data values */

  for (iech=0; iech<nech; iech++)
    for (ivar=0; ivar<nvar; ivar++)
      if (FFFF(st_get_ivar(rank[iech],ivar))) FLAG(iech,ivar) = 0;

  /* Check on the external drifts */

  for (iech=0; iech<nech; iech++)
    for (ibfl=0; ibfl<model_nfex(model); ibfl++)
      if (FFFF(st_get_fext(rank[iech],ibfl)))
        for (ivar=0; ivar<DBIN->getVariableNumber(); ivar++) FLAG(iech,ivar) = 0;

  /* Check on the drift */

  for (ib=0; ib<model->getDriftEquationNumber(); ib++)
  {
    valid = 0;
    for (il=0; il<model->getDriftNumber(); il++)
      for (ivar=0; ivar<nvar; ivar++)
      {
        if (model->getCoefDrift(ivar,il,ib) == 0.) continue;
        for (iech=0; iech<nech; iech++)
          if (! FFFF(st_get_ivar(rank[iech],ivar))) valid++;
      }
    FLAG(nech+ib,DBIN->getVariableNumber()-1) = (valid > 0);
  }
  
  /* Calculate the new number of equations */

  for (i=count=0; i<neq; i++)
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
** \param[in]  nech   Number of active samples
**
*****************************************************************************/
static int st_authorize(Model *model,
                        int    nech)
{
  int i,nvar,nfeq,n_cov,n_drf,error;

  /* Initializations */

  error = 1;
  nvar  = model->getVariableNumber();
  nfeq  = model->getDriftEquationNumber();

  /* Preliminary check */

  if (nech * nvar < nfeq) goto label_end;

  /* Check that enough information is present */

  n_cov = n_drf = 0;
  for (i=0; i<nvar * nech; i++) n_cov += flag[i];
  for (i=0; i<nfeq; i++) n_drf += flag[i + nvar * nech];
  if (n_cov <= 0 || n_cov < n_drf) goto label_end;

  /* Set the error return code */

  error = 0;

label_end:
  return(error);
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
static void st_variance0(Model  *model,
                         int     nvar,
                         VectorDouble matCL)
{
  double value;
  int    nvar_m;

  /* Initializations */

  nvar_m = model->getVariableNumber();

  /* In the non-stationary case, the calculation is postponed. */
  /* The value is temporarily set to TEST                      */
  
  if (model->isNoStat())
  {
    for (int ivar=0; ivar<nvar; ivar++)
      for (int jvar=0; jvar<nvar; jvar++)
        VAR0(ivar,jvar) = TEST;
  }
  else
  {
    COV_EXTERNAL_DB_1  = 2;
    COV_EXTERNAL_DB_2  = 2;
    COV_EXTERNAL_ECH_1 = 0;
    COV_EXTERNAL_ECH_2 = 0;
    NDIM_EXTERNAL      = model->getDimensionNumber();
    model_variance0(model,KOPTION,covtab,var0);
    
    // If 'mat' is provided, some extra-calculations are needed
    
    if (! matCL.empty())
    {
      for (int ivar=0; ivar<nvar; ivar++)
        for (int jvar=0; jvar<nvar; jvar++)
        {
          value = 0.;
          for (int ivar_m=0; ivar_m<model->getVariableNumber(); ivar_m++)
            for (int jvar_m=0; jvar_m<model->getVariableNumber(); jvar_m++)
              value += (MATCL(ivar,ivar_m) * COVTAB(ivar_m,jvar_m) * 
                        MATCL(jvar,jvar_m));
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
static double st_variance(Model *model,
                          int    ivar,
                          int    jvar,
                          int    nred)
{
  double var;
  int i,nvar;

  nvar = model->getVariableNumber();

  // In the stationary case, var0 has already been calculated
  // Otherwise if must be derived here

  if (model->isNoStat() || model->hasExternalCov())
  {
    COV_EXTERNAL_DB_1  = 2;
    COV_EXTERNAL_DB_2  = 2;
    COV_EXTERNAL_ECH_1 = IECH_OUT;
    COV_EXTERNAL_ECH_2 = IECH_OUT;
    NDIM_EXTERNAL      = model->getDimensionNumber();
    model_variance0_nostat(model,KOPTION,DBOUT,IECH_OUT,covtab,var0);
  }

  var = VAR0(ivar,jvar);
  if (FLAG_BAYES) var += VARB(ivar,ivar);
  for (i=0; i<nred; i++) 
    var -= RHS_C(i,jvar) * WGT(i,ivar);

  return(var);
}

/****************************************************************************/
/*!
**  Establish the variance of the estimator
**
** \param[in]  model   Model structure
** \param[in]  ivar    Rank of the target variable
** \param[in]  jvar    Rank of the auxiliary variable
** \param[in]  nfeq    Number of drift equations
** \param[in]  nred    Reduced number of equations
**
*****************************************************************************/
static double st_varestimate(Model *model,
                             int    ivar,
                             int    jvar,
                             int    nfeq,
                             int    nred)
{
  double var,signe;
  int i,cumflag;

  cumflag = nred - nfeq;

  var = 0.;
  for (i=0; i<nred; i++)
  {
    signe = (i < cumflag) ? 1. : -1.;
    var += signe * RHS_C(i,jvar) * WGT(i,ivar);
  }

  return(var);
}

/****************************************************************************/
/*!
**  Establish the kriging L.H.S.
**
** \param[in]  model  Model structure
** \param[in]  neigh  Neigh structure
** \param[in]  neq    Number of equations
** \param[in]  nech   Number of active data points
**
*****************************************************************************/
static void st_lhs(Model  *model,
                   Neigh  *neigh,
                   int     nech,
                   int     neq)
{
  int    i,iech,jech,idim,ivar,jvar,ib,il,nvar_m,nfeq,nbfl,code1,code2;
  double verr,cref,value;

  /* Initializations */

  nvar_m  = model->getVariableNumber();
  nfeq    = model->getDriftEquationNumber();
  nbfl    = model->getDriftNumber();
  for (i=0; i<neq*neq; i++) lhs[i] = 0.;

  /* Establish the covariance part */

  for (iech=0; iech<nech; iech++)
    for (jech=0; jech<nech; jech++)
    {
      model_covtab_init(1,model,covtab);
      for (idim=0; idim<DBIN->getNDim(); idim++)
        d1[idim] = (st_get_idim(rank[jech],idim) -
                    st_get_idim(rank[iech],idim));
      st_cov_dd(model,0,0,0,MEMBER_LHS,-1,1.,rank[iech],rank[jech],d1,covtab);

      for (ivar=0; ivar<nvar_m; ivar++)
        for (jvar=0; jvar<nvar_m; jvar++)
        {
          LHS(iech,ivar,jech,jvar) = COVTAB(ivar,jvar);

          /* Correction due to measurement errors */

          verr = 0.;
          if (FLAG_PROF)
          {
            code1 = (int) DBIN->getCode(rank[iech]);
            code2 = (int) DBIN->getCode(rank[jech]);
            if (code1 != 0 && code2 != 0 && code1 == code2)
              verr = DBIN->getVarianceError(rank[iech],0);
          }
          else
          {
            if (iech == jech)
            {
              verr = DBIN->getVarianceError(rank[iech],ivar);

              if (flag_continuous_kriging(neigh))
              {
                // In the case of continuous Kriging, we must update the LHS
                // by considering the distance between data and target
                
                cref = LHS(iech,ivar,jech,jvar);
                verr = cref * neigh_continuous_variance(neigh,
                                                        DBIN,rank[iech],
                                                        DBOUT,IECH_OUT);
              }
            }
          }
          if (! FFFF(verr) && verr > 0) LHS(iech,ivar,jech,jvar) += verr;

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
  for (iech=0; iech<nech; iech++)
  {
    if (rank[iech] >= 0)
      model_calcul_drift(model,MEMBER_LHS,DBIN,rank[iech],drftab);
    else
      model_calcul_drift(model,MEMBER_LHS,DBOUT,IECH_OUT,drftab);

    for (ivar=0; ivar<nvar_m; ivar++)
      for (ib=0; ib<nfeq; ib++)
      {
        value = 0.;
        for (il=0; il<nbfl; il++)
          value += drftab[il] * model->getCoefDrift(ivar,il,ib);
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
  int lec_lhs,ecr_lhs,i,j;

  /* Loop on the elements */

  lec_lhs = ecr_lhs = 0;
  for (i=0; i<neq; i++)
    for (j=0; j<neq; j++)
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
** \param[in]  nech  Number of active points (optional)
** \param[in]  neq   Number of equations
** \param[in]  nred  Reduced number of equations
** \param[in]  flag  Flag array (optional)
** \param[in]  lhs   Kriging L.H.S
**
*****************************************************************************/
GEOSLIB_API void krige_lhs_print(int     nech,
                                 int     neq,
                                 int     nred,
                                 int    *flag,
                                 double *lhs)
{
  int *rel,i,j,ipass,npass,ideb,ifin;

  /* Initializations */

  rel   = (int *) NULL;
  rel   = st_relative_position_array(1,neq,rel);
  npass = (nred - 1) / NBYPAS + 1;

  /* General Header */

  mestitle(0,"LHS of Kriging matrix (compressed)");
  if (nech > 0)
    message ("Number of active samples    = %d\n", nech);
  message ("Total number of equations   = %d\n", neq);
  message ("Reduced number of equations = %d\n", nred);

  /* Loop on the passes */

  for (ipass = 0; ipass < npass; ipass++)
  {
    ideb = ipass * NBYPAS;
    ifin = MIN (nred, ideb + NBYPAS);
    message("\n");

    /* Header line */

    tab_prints(NULL,1,GD_J_RIGHT,"Rank");
    tab_prints(NULL,1,GD_J_RIGHT,"    ");
    for (j = ideb; j < ifin; j++) tab_printi(NULL,1,GD_J_RIGHT,j+1);
    message("\n");

    /* Flag line */

    if (flag != NULL)
    {
      tab_prints(NULL,1,GD_J_RIGHT,"    ");
      tab_prints(NULL,1,GD_J_RIGHT,"Flag");
      for (j = ideb; j < ifin; j++) tab_printi(NULL,1,GD_J_RIGHT,rel[j]);
      message("\n");
    }

    /* Matrix lines */

    for (i = 0; i < nred; i++)
    {
      tab_printi(NULL,1,GD_J_RIGHT,i+1);
      tab_printi(NULL,1,GD_J_RIGHT,rel[i]);
      for (j = ideb; j < ifin; j++)
        tab_printg(NULL,1,GD_J_RIGHT,LHS_C(i,j));
      message("\n");
    }
  }

  rel = st_relative_position_array(-1,neq,rel);
  return;
}

/****************************************************************************/
/*!
**  Define the neighborhood
**
** \return  1 if a new Neighborhood has been found; 0 otherwise
**
** \param[in]  neigh     Neigh structure
**
** \param[out]  status   Neighborhood error status
** \param[out]  nech     Number of active data points
**
*****************************************************************************/
static int st_neigh(Neigh   *neigh,
                    int     *status,
                    int     *nech)
{
  int iech,idim,ndim,ivar,jvar,nvarin,nval,found,flag_new;
  static int nech_mem;

  /* Initializations */

  *status  = 0;
  flag_new = 0;
  nvarin   = DBIN->getVariableNumber();

  /* Should the neighborhood search be performed again */

  switch (neigh->getType())
  {
    case NEIGH_UNIQUE:
      if (IECH_NBGH >= 0) goto label_suite;
      break;

    case NEIGH_IMAGE:
      if (IECH_NBGH >= 0) goto label_suite;
      break;

    case NEIGH_BENCH:
      ndim = DBOUT->getNDim();
      if (IECH_NBGH < 0 || IECH_NBGH > get_NECH(DBOUT)) break;
      if (IECH_OUT  < 0 || IECH_OUT  > get_NECH(DBOUT)) break;
      if (is_grid(DBOUT))
      {
        nval = 1;
        for (idim=0; idim<ndim-1; idim++) nval *= DBOUT->getNX(idim);
        if ((IECH_OUT / nval) == (IECH_NBGH / nval)) goto label_suite;
      }
      else
      {
        if (get_IDIM(DBOUT,IECH_NBGH,ndim-1) ==
            get_IDIM(DBOUT,IECH_OUT ,ndim-1)) goto label_suite;
      }
      break;

    case NEIGH_MOVING:
      if (IECH_NBGH == IECH_OUT) goto label_suite;
      break;
  }

  /* Perform the neighborhood search */

  IECH_NBGH = IECH_OUT;
  *status   = neigh_select(DBIN,DBOUT,IECH_OUT,neigh,FLAG_SIMU,&nech_mem,rank);
  *nech     = nech_mem;
  flag_new  = 1;

label_suite:
  if (FLAG_COLK)
  {
    *nech = nech_mem;

    /* Do not add the target if no variable is defined */
    for (ivar=found=0; ivar<nvarin && found<=0; ivar++)
    {
      jvar = RANK_COLCOK[ivar];
      if (jvar < 0) continue;
      if (! FFFF(get_ARRAY(DBOUT,IECH_OUT,jvar))) found = 1;
    }
    if (! found) return(flag_new);

    /* Do not add the target if it coincides with a datum */
    for (iech=0; iech<nech_mem; iech++)
    {
      if (distance_inter(DBIN,DBOUT,rank[iech],IECH_OUT,NULL) <= 0.)
        return(flag_new);
    }

    /* Add the target */

    rank[*nech] = -1;
    *nech = nech_mem + 1;
    flag_new = 1;
  }

  return(flag_new);
}

/****************************************************************************/
/*!
**  Extract the valid data
**  Operate the product by the inverse covariance matrix
**
** \param[in]  model  Model structure
** \param[in]  rmean  Array giving the posterior means for the drift terms
** \param[in]  nech   Number of active samples
** \param[in]  nred   Reduced number of equations
**
** \param[out] lterm  Product Z*C-1*Z 
**                    (only produced if FLAG_LTERM is true)
**
*****************************************************************************/
static void st_data_dual(Model  *model,
                         double *rmean,
                         int     nech,
                         int     nred,
                         double *lterm)
{
  int i,iech,ivar,ecr,nvar,nfeq;
  double mean;

  /* Initializations */

  nvar = model->getVariableNumber();
  nfeq = model->getDriftEquationNumber();

  /* Set the whole array to 0 */

  for (i=0; i<nred; i++) zext[i] = 0.;

  /* Extract the data */

  for (ivar=ecr=0; ivar<nvar; ivar++)
  {
    for (iech=0; iech<nech; iech++)
    {
      if (! FLAG(iech,ivar)) continue;
      mean = 0.;
      if (nfeq <= 0) mean = model->getContext().getMean(ivar);
      if (FLAG_BAYES)
        mean = model_drift_evaluate(1,model,DBIN,rank[iech],ivar,
                                    rmean,drftab);
      zext[ecr++] = st_get_ivar(rank[iech],ivar) - mean;
    }
  }

  /* Operate the product : Z * A-1 */

  matrix_product(nred,nred,1,lhs,zext,zam1);

  /* Operate the product : Z * A-1 * Z */

  if (FLAG_LTERM) matrix_product(1,nred,1,zam1,zext,lterm);

  return;
}

/****************************************************************************/
/*!
**  Define the array flag[] and the kriging L.H.S.
**
** \param[in]  model     Model structure
** \param[in]  neigh     Neigh structure
** \param[in]  nech      Number of selected data
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
static void st_prepar(Model   *model,
                      Neigh   *neigh,
                      int      nech,
                      int     *status,
                      int     *nred_r,
                      int     *neq_r)
{
  int nred,neq;

  /* Initializations */

  *nred_r = 0;
  *neq_r  = 0;
  *status = 1;
  neq = nech * model->getVariableNumber() + model->getDriftEquationNumber();

  /* Define the array flag */

  st_flag_define(model,nech,neq,&nred);

  /* Check if the number of points is compatible with the model */

  if (st_authorize(model,nech)) return;

  /* Establish the Kriging L.H.S. */

  st_lhs(model,neigh,nech,neq);
  st_lhs_iso2hetero(neq);

  if (debug_query("kriging")) krige_lhs_print(nech,neq,nred,flag,lhs);

  /* Backup the Kriging matrix before inversion */

  (void) memcpy((char *) lhs_b,(char *) lhs,nred * nred * sizeof(double));

  /* Invert the L.H.S. matrix */

  if (matrix_invert(lhs,nred,IECH_OUT)) 
  {
    messerr("The Kriging Matrix (%d,%d) is singular",nred,nred);
    messerr("One of the usual reason is the presence of duplicates");
    messerr("To display the Kriging Matrix, run the Kriging procedure again");
    messerr("typing the following command first:");
    messerr("  set.keypair('Dump_Kriging_Matrix',1)");
    if (get_keypone("Dump_Kriging_Matrix",0))
      print_matrix("Singular matrix",0,1,nred,nred,NULL,lhs);
    return;
  }

  /* Returning arguments */

  *status = 0;
  *nred_r = nred;
  *neq_r  = neq;

  return;
}

/****************************************************************************/
/*!
**  Establish the kriging R.H.S
**
** \param[in]  model    Model structure
** \param[in]  nech     Number of active data points
** \param[in]  neq      Number of equations
** \param[in]  nvar     Number of output variables
** \param[in]  matCL    Matrix of linear combinaison (or NULL)
**                      (Dimension: nvar * model->getNVar())
**
** \param[out]  status  Kriging error status
**
** \remarks When 'matCL' is provided, 'nvar' stands for the first dimension of
** \remarks the matrix 'matCL' (its second dimension is equal to model->getNVar()).
** \remarks Otherwise nvar designates model->getNVar()
**
*****************************************************************************/
static void st_rhs(Model   *model,
                   int      nech,
                   int      neq,
                   int      nvar,
                   double  *matCL,
                   int     *status)
{
  int    i,iech,il,ib,nvar_m,nbfl,nfeq,idim,nscale;
  double value,ratio;

  /* Initializations */

  nscale = 1;
  nvar_m = model->getVariableNumber();
  nbfl   = model->getDriftNumber();
  nfeq   = model->getDriftEquationNumber();

  /* Establish the covariance part */

  for (iech=0; iech<nech; iech++)
  {
    switch (KOPTION->calcul)
    {
      case KOPTION_PONCTUAL:
        nscale = 1;
        model_covtab_init(1,model,covtab);
        for (idim=0; idim<DBIN->getNDim(); idim++)
        {
          d1[idim] = (get_IDIM(DBOUT,IECH_OUT,idim) -
                      st_get_idim(rank[iech],idim));
          // The next option is plugged for the case of target randomization
          // for the case of Point-Block Model
          if (RAND_INDEX >= 0 && KOPTION->disc1 != (double *) NULL)
            d1[idim] += DISC1(RAND_INDEX,idim);
        }
        st_cov_dg(model,0,0,0,MEMBER_RHS,-1,1.,rank[iech],-1,d1,covtab);
        break;
        
      case KOPTION_BLOCK:
        nscale = KOPTION->ntot;
        model_covtab_init(1,model,covtab);
        for (i=0; i<nscale; i++)
        {
          for (idim=0; idim<DBIN->getNDim(); idim++)
            d1[idim] = (get_IDIM(DBOUT,IECH_OUT,idim) -
                        st_get_idim(rank[iech],idim) + DISC1(i,idim));
          st_cov_dg(model,0,0,0,MEMBER_RHS,-1,1.,rank[iech],-1,d1,covtab);
        }
        break;
          
      case KOPTION_DRIFT:
        nscale = 1;
        for (int ivar_m=0; ivar_m<nvar_m; ivar_m++)
          for (int jvar_m=0; jvar_m<nvar_m; jvar_m++)
            COVTAB(ivar_m,jvar_m) = 0;
        break;
    }
    
    /* Normation */
    
    ratio = 1. / (double) nscale;
    if (FLAG_DGM) ratio *= R_COEFF;
    for (int ivar_m=0; ivar_m<nvar_m; ivar_m++)
      for (int jvar_m=0; jvar_m<nvar_m; jvar_m++)
        COVTAB(ivar_m,jvar_m) *= ratio;
    
    if (matCL == (double *) NULL)
    {
      for (int jvar=0; jvar<nvar; jvar++)
        for (int ivar=0; ivar<nvar; ivar++)
          RHS(iech,ivar,jvar) = COVTAB(ivar,jvar);
    }
    else
    {
      for (int jvar=0; jvar<nvar; jvar++)
        for (int ivar_m=0; ivar_m<nvar_m; ivar_m++)
        {
          value = 0.;
          for (int jvar_m=0; jvar_m<nvar_m; jvar_m++)
            value += MATCL(jvar,jvar_m) * COVTAB(ivar_m,jvar_m);
          RHS(iech,ivar_m,jvar) = value;
        }
    }
  }
  
  /* Establish the drift part */
  
  if (nfeq <= 0) return;
  
  model_calcul_drift(model,MEMBER_RHS,DBOUT,IECH_OUT,drftab);
  for (il=0; il<nbfl; il++)
    if (FFFF(drftab[il]))
    {
      *status = 1;
      return;
    }
  
  if (matCL != (double *) NULL)
  {
    if (model->isFlagLinked())
      messageAbort("Kriging of a Linear combination is incompatible with linked drifts");
    for (int ivar=0; ivar<nvar; ivar++)
      for (int jvar_m=ib=0; jvar_m<nvar_m; jvar_m++)
        for (int jl=0; jl<nbfl; jl++, ib++)
        {
          value = 0.;
          for (int il=0; il<nbfl; il++)
            value += drftab[il] * model->getCoefDrift(jvar_m,il,ib);
          value *= MATCL(ivar,jvar_m);
          RHS(ib,nvar_m,ivar) = value; /* nvar_m is used to calculate shift */
        }
  }
  else
  {
    
    for (int ivar=0; ivar<nvar; ivar++)
      for (ib=0; ib<nfeq; ib++)
      {
        value = 0.;
        for (int il=0; il<nbfl; il++)
          value += drftab[il] * model->getCoefDrift(ivar,il,ib);
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
static void st_ff0(Model   *model,
                   int     *status)
{
  int    ivar,ib,il,nvar,nbfl,nfeq;
  double value;

  /* Initializations */

  nvar   = model->getVariableNumber();
  nbfl   = model->getDriftNumber();
  nfeq   = model->getDriftEquationNumber();

  /* Establish the drift part */

  if (nbfl <= 0 || nfeq <= 0) return;

  model_calcul_drift(model,MEMBER_RHS,DBOUT,IECH_OUT,drftab);
  for (il=0; il<nbfl; il++)
    if (FFFF(drftab[il]))
    {
      *status = 1;
      return;
    }

  for (ivar=0; ivar<nvar; ivar++)
    for (ib=0; ib<nfeq; ib++)
    {
      value = 0.;
      for (il=0; il<nbfl; il++)
        value += drftab[il] * model->getCoefDrift(ivar,il,ib);
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
static void st_rhs_iso2hetero(int neq,
                              int nvar)
{
  int i,ivar,lec_rhs,ecr_rhs;

  /* Loop on the elements */

  lec_rhs = ecr_rhs = 0;
  for (ivar=0; ivar<nvar; ivar++)
    for (i=0; i<neq; i++, lec_rhs++)
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
GEOSLIB_API void krige_rhs_print(int      nvar,
                                 int      nech,
                                 int      neq,
                                 int      nred,
                                 int     *flag,
                                 double  *rhs)
{
  int *rel,i,ivar,idim;

  /* Initializations */

  rel = (int *) NULL;
  rel = st_relative_position_array(1,neq,NULL);

  /* General Header */

  mestitle(0,"RHS of Kriging matrix (compressed)");
  if (nech > 0)
    message ("Number of active samples    = %d\n", nech);
  message ("Total number of equations   = %d\n", neq);
  message ("Reduced number of equations = %d\n", nred);
  message ("Number of right-hand sides  = %d\n", nvar);

  /* Kriging option */

  if (KOPTION != (Koption *) NULL)
  {
    switch (KOPTION->calcul)
    {
      case KOPTION_PONCTUAL:
        message ("Ponctual Estimation\n");
        break;

      case KOPTION_BLOCK:
        message ("Block Estimation : Discretization = ");
        for (idim=0; idim<KOPTION->ndim; idim++) {
          if (idim != 0) message(" x ");
          message("%d",KOPTION->ndisc[idim]);
        }
        message("\n");
        break;

      case KOPTION_DRIFT:
        message ("Drift Estimation\n");
        break;
    }
  }
  message ("\n");

  /* Header line */

  tab_prints(NULL,1,GD_J_RIGHT,"Rank");
  if (flag != (int *) NULL) tab_prints(NULL,1,GD_J_RIGHT,"Flag");
  for (ivar = 0; ivar < nvar; ivar++) tab_printi(NULL,1,GD_J_RIGHT,ivar+1);
  message("\n");

  /* Matrix lines */

  for (i = 0; i < nred; i++)
  {
    tab_printi(NULL,1,GD_J_RIGHT,i+1);
    if (flag != (int *) NULL) tab_printi(NULL,1,GD_J_RIGHT,rel[i]);
    for (ivar = 0; ivar < nvar; ivar++)
      tab_printg(NULL,1,GD_J_RIGHT,RHS_C(i,ivar));
    message("\n");
  }

  rel = st_relative_position_array(-1,neq,rel);
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
GEOSLIB_API void krige_dual_print(int      nech,
                                  int      neq,
                                  int      nred,
                                  int     *flag,
                                  double  *dual)
{
  int *rel,i;

  /* Initializations */

  rel = (int *) NULL;
  rel = st_relative_position_array(1,neq,NULL);

  /* General Header */

  mestitle(0,"Dual Vector (completed with zeroes and compressed)");
  if (nech > 0)
    message ("Number of active samples    = %d\n", nech);
  message ("Total number of equations   = %d\n", neq);
  message ("Reduced number of equations = %d\n", nred);

  /* Header line */

  tab_prints(NULL,1,GD_J_RIGHT,"Rank");
  if (flag != (int *) NULL) tab_prints(NULL,1,GD_J_RIGHT,"Flag");
  message("\n");

  /* Matrix lines */

  for (i = 0; i < nred; i++)
  {
    tab_printi(NULL,1,GD_J_RIGHT,i+1);
    if (flag != (int *) NULL) tab_printi(NULL,1,GD_J_RIGHT,rel[i]);
    tab_printg(NULL,1,GD_J_RIGHT,dual[i]);
    message("\n");
  }

  rel = st_relative_position_array(-1,neq,rel);
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
** \param[in]  nech         Number of active samples
** \param[in]  nvar         Number of output variables
** \param[in]  nred         Reduced number of equations
**
** \remarks When 'matCL' is provided, 'nvar' stands for the first dimension of
** \remarks the matrix 'matCL' (its second dimension is equal to model->getNVar()).
** \remarks Otherwise nvar designates model->getNVar()
**
*****************************************************************************/
static void st_estimate(Model  *model,
                        double *rmean,
                        int     status,
                        int     flag_xvalid,
                        int     nech,
                        int     nvar,
                        int     nred)
{
  int    i,ivar,nfeq;
  double estim,stdv,var,valdat;

  /* Initializations */

  nfeq = model->getDriftEquationNumber();

  /* Estimation */

  if (FLAG_EST)
  {
    for (ivar=0; ivar<nvar; ivar++)
    {
      estim = 0.;
      if (nfeq <= 0) estim = model->getContext().getMean(ivar);

      if (status == 0 && (nred > 0 || nfeq <= 0 || FLAG_BAYES))
      {
        if (FLAG_BAYES)
          estim = model_drift_evaluate(0,model,DBOUT,IECH_OUT,ivar,
                                       rmean,drftab);
        for (i=0; i<nred; i++) estim += RHS_C(i,ivar) * ZAM1(i);
        
        if (ABS(flag_xvalid) == 1)
        {
          valdat = DBIN->getVariable(IECH_OUT,ivar);
          estim  = (FFFF(valdat)) ? TEST : estim - valdat;
        }
      }
      else
      {
        // In case of failure with KS, set the result to mean
        if (nfeq > 0) estim = TEST;
      }
      DBOUT->setArray(IECH_OUT,IPTR_EST+ivar,estim);
    }
  }

  /* Variance of the estimation error */

  if (FLAG_STD)
  {
    for (ivar=0; ivar<nvar; ivar++)
    {
      if (status == 0 && (nred > 0 || nfeq <= 0 || FLAG_BAYES))
      {
        stdv = st_variance(model,ivar,ivar,nred);
        if (stdv < 0) stdv = 0.;
        
        stdv = sqrt(stdv);

        if (ABS(flag_xvalid) == 1)
        {
          estim = get_ARRAY(DBOUT,IECH_OUT,IPTR_EST+ivar);
          stdv  = (FFFF(estim) || stdv <= 0.) ? TEST : estim / stdv;
        }
      }
      else
      {
        stdv = TEST;
      }
      DBOUT->setArray(IECH_OUT,IPTR_STD+ivar,stdv);
    }
  }

  /* Variance of the estimator */

  if (FLAG_VARZ)
  {
    for (ivar=0; ivar<nvar; ivar++)
    {
      if (status == 0 && (nred > 0 || nfeq <= 0))
        var = st_varestimate(model,ivar,ivar,nfeq,nred);
      else
        var = TEST;
      DBOUT->setArray(IECH_OUT,IPTR_VARZ+ivar,var);
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
** \param[in]  nech      Number of active samples
** \param[in]  nred      Reduced number of equations
**
** \remark  KS and BAYES are incompatible: we can use mean in both cases
**
*****************************************************************************/
static void st_simulate(Model  *model,
                        double *smean,
                        int     status,
                        int     icase,
                        int     nbsimu,
                        int     nech,
                        int     nred)
{
  int    isimu,iech,jech,ivar,jvar,lec,ecr,nvar,nfeq;
  double simu,mean,data,value;

  /* Initializations */

  nvar = model->getVariableNumber();
  nfeq = model->getDriftEquationNumber();

  /* Simulation */

  for (isimu=ecr=0; isimu<nbsimu; isimu++)
    for (ivar=0; ivar<nvar; ivar++,ecr++)
    {
      simu = 0.;
      if (nfeq <= 0) simu = model->getContext().getMean(ivar);

      if (status == 0)
      {
        if (FLAG_BAYES)
          simu = model_drift_evaluate(0,model,DBOUT,IECH_OUT,ivar,
                                      &SMEAN(0,isimu),drftab);
        
        lec  = ivar * nred;
        for (jvar=0; jvar<nvar; jvar++)
          for (iech=0; iech<nech; iech++)
          {
            if (! FLAG(iech,jvar)) continue;
            jech = rank[iech];
            
            mean = 0.;
            if (nfeq <= 0) mean = model->getMean(jvar);
            if (FLAG_BAYES)
              mean = model_drift_evaluate(1,model,DBIN,jech,jvar,
                                          &SMEAN(0,isimu),drftab);
            data  = st_get_array(jech,isimu,jvar,icase,nbsimu,nvar);
            simu -= wgt[lec++] * (data + mean);
          }
        
        if (debug_query("kriging"))
        {
          value = get_ARRAY(DBOUT,IECH_OUT,IPTR_EST+ecr);
          message("Non-conditional simulation #%d = %lf\n",
                  isimu+1,value);
          message("Kriged difference = %lf\n",-simu);
          message("Conditional simulation #%d = %lf\n",
                  isimu+1,value + simu);
        }
      }
      else
      {
        // In case of failure with KS, set the contidioning to mean
        if (nfeq > 0) simu = TEST;
      }
      
      /* Add the conditioning kriging to the NC simulation at target */
      DBOUT->updSimvar(LOC_SIMU, IECH_OUT, isimu, ivar, icase, nbsimu, nvar, 0,
                       simu);
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
** \param[in]  nech    Number of active points
** \param[in]  nred    Reduced number of equations
** \param[in]  icase   Rank of the PGS or GRF
** \param[in]  flag    Flag array
** \param[in]  wgt     Array of Kriging weights
**
** \remark In the case of simulations (icase>=0), the data vector is not 
** \remark printed as it changes for every sample, per simulation
**
*****************************************************************************/
static void krige_wgt_print(int     status,
                            int     nvar,
                            int     nvar_m,
                            int     nfeq,
                            int     nech,
                            int     nred,
                            int     icase,
                            int    *flag,
                            double *wgt)
{
  double *sum,value;
  int iwgt,ivar,jvar_m,ivar_m,iech,lec,cumflag,idim,ndim,ib,number,flag_value;

  /* Initializations */

  ndim  = DBIN->getNDim();
  sum   = (double *) st_core(nvar_m,1);
  if (sum == (double *) NULL) return;

  /* Header */

  mestitle (0,"(Co-) Kriging weights");

  /* First line */

  tab_prints(NULL,1,GD_J_RIGHT,"Rank");
  for (idim=0; idim<ndim; idim++)
  {
    String strloc = getLocatorName(LOC_X,idim+1);
    tab_prints(NULL,1,GD_J_RIGHT,strloc.c_str());
  }
  if (DBIN->hasCode())
    tab_prints(NULL,1,GD_J_RIGHT,"Code");
  if (DBIN->getVarianceErrorNumber() > 0)
    tab_prints(NULL,1,GD_J_RIGHT,"Err.");
  if (KOPTION->flag_data_disc)
    for (idim=0; idim<ndim; idim++)
    {
      (void) sprintf(string,"Size%d",idim+1);
      tab_prints(NULL,1,GD_J_RIGHT,string);
    }
  tab_prints(NULL,1,GD_J_RIGHT,"Data");
  for (ivar = 0; ivar < nvar; ivar++)
  {
    (void) sprintf(string,"Z%d*",ivar+1);
    tab_prints(NULL,1,GD_J_RIGHT,string);
  }
  message("\n");

  /* Display the information and the weights */

  for (jvar_m = lec = cumflag = 0; jvar_m < nvar_m; jvar_m++)
  {
    if (nvar > 1) message("Using variable Z%-2d\n", jvar_m+1);

    /* Loop on the samples */

    for (ivar_m = 0; ivar_m < nvar_m; ivar_m++) sum[ivar_m] = 0.;
    for (iech = 0; iech < nech; iech++, lec++)
    {
      flag_value = (flag != (int *) NULL) ? flag[lec] : 1;
      tab_printi(NULL,1,GD_J_RIGHT,iech+1);
      for (idim=0; idim<ndim; idim++)
        tab_printg(NULL,1,GD_J_RIGHT,st_get_idim(rank[iech],idim));
      if (DBIN->hasCode())
        tab_printg(NULL,1,GD_J_RIGHT,DBIN->getCode(rank[iech]));
      if (DBIN->getVarianceErrorNumber() > 0)
        tab_printg(NULL,1,GD_J_RIGHT,st_get_verr(rank[iech],
                                                 (FLAG_PROF) ? 0 : jvar_m));
      if (KOPTION->flag_data_disc)
      {
        for (idim=0; idim<ndim; idim++)
          tab_printg(NULL,1,GD_J_RIGHT,DBIN->getBlockExtension(rank[iech],idim));
      }
      if (icase < 0)
        tab_printg(NULL,1,GD_J_RIGHT,st_get_ivar(rank[iech],jvar_m));
      else
        tab_prints(NULL,1,GD_J_RIGHT,"   ");

      for (ivar = 0; ivar < nvar; ivar++)
      {
        iwgt  = nred * ivar + cumflag;
        value = (wgt != (double *) NULL && status == 0 && flag_value) ?
          wgt[iwgt] : TEST;
        if (! FFFF(value)) sum[ivar] += value;
        tab_printg(NULL,1,GD_J_RIGHT,value);
      }
      if (flag_value) cumflag++;
      message("\n");
    }

    number = 1 + ndim + 1;
    if (DBIN->getVarianceErrorNumber() > 0) number++;
    if (KOPTION->flag_data_disc) number += ndim;
    tab_prints(NULL,number,GD_J_LEFT,"Sum of weights");
    for (ivar = 0; ivar < nvar; ivar++)
    {
      value = (status == 0) ? sum[ivar] : TEST;
      tab_printg(NULL,1,GD_J_RIGHT,value);
    }
    message("\n");
  }

  sum = (double *) mem_free((char *) sum);
  if (nfeq <= 0 || wgt == (double *) NULL) return;

  /* Header */

  mestitle (0,"Drift coefficients");

  /* First line */

  tab_prints(NULL,1,GD_J_RIGHT,"Rank");
  tab_prints(NULL,1,GD_J_RIGHT,"Coeff");
  message("\n");

  /* Loop on the drift coefficients */

  cumflag = nred - nfeq;
  for (ib=0; ib<nfeq; ib++)
  {
    iwgt = ib + cumflag;
    tab_printi(NULL,1,GD_J_RIGHT,ib+1);
    value = (status == 0) ? zam1[iwgt] : TEST;
    tab_printg(NULL,1,GD_J_RIGHT,value);
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
static void st_store_nbgh(int     status,
                          int     ntab,
                          double *tab)
{
  double value;
  int i;

  /* Loop on the parameters */

  for (i=0; i<ntab; i++)
  {

    /* Store the parameter */

    value = (status == 0) ? tab[i] : TEST;
    DBOUT->setArray(IECH_OUT,IPTR_NBGH+i,value);
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
static void st_result_kriging_print(int flag_xvalid,
                                    int nvar,
                                    int status)
{
  int    ivar;
  double value,estim,estval,esterr,sigma,trueval,sterr,stdev;

  /* Header */

  if (flag_xvalid != 0)
    mestitle(0,"Cross-validation results");
  else
    mestitle(0,"(Co-) Kriging results");
  message("Target Sample = %d\n",IECH_OUT+1);

  /* Loop on the results */

  for (ivar = 0; ivar <nvar; ivar++)
  {
    if (flag_xvalid != 0)
    {
      message("Variable Z%-2d\n",ivar+1);
      if (FLAG_EST)
      {
        trueval = (status == 0) ? DBIN->getVariable(IECH_OUT,ivar) : TEST;
        estim   = (status == 0) ? get_ARRAY(DBOUT,IECH_OUT,IPTR_EST+ivar) :TEST;

        if (ABS(flag_xvalid) == 1)
        {
          estval  = (status == 0) ? estim + trueval : TEST;
          esterr  = (status == 0) ? estim : TEST;
        }
        else
        {
          estval  = (status == 0) ? estim : TEST;
          esterr  = (status == 0) ? estim - trueval : TEST;
        }

        tab_printg(" - True value        = ",1,GD_J_RIGHT,trueval);
        message("\n");
        tab_printg(" - Estimated value   = ",1,GD_J_RIGHT,estval);
        message("\n");
        tab_printg(" - Estimation Error  = ",1,GD_J_RIGHT,esterr);
        message("\n");

        if (FLAG_STD)
        {
          stdev = (status == 0) ? get_ARRAY(DBOUT,IECH_OUT,IPTR_STD+ivar):TEST;

          if (ABS(flag_xvalid) == 1)
          {
            sterr = stdev;
            sigma = (status == 0) ? esterr / stdev : TEST;
          }
          else
          {
            sigma = stdev;
            sterr = (status == 0) ? esterr / stdev : TEST;
          }

          tab_printg(" - Std. deviation    = ",1,GD_J_RIGHT,sigma);
          message("\n");
          tab_printg(" - Normalized Error  = ",1,GD_J_RIGHT,sterr);
          message("\n");
        }
      }
    }
    else
    {
      message("Variable Z%-2d\n",ivar+1);
      if (FLAG_EST)
      {
        value = (status == 0) ? get_ARRAY(DBOUT,IECH_OUT,IPTR_EST+ivar) : TEST;
        tab_printg(" - Estimate  = ",1,GD_J_RIGHT,value);
        message("\n");
      }
      if (FLAG_STD)
      {
        value = (status == 0) ? get_ARRAY(DBOUT,IECH_OUT,IPTR_STD+ivar) : TEST;
        tab_printg(" - Std. Dev. = ",1,GD_J_RIGHT,value);
        value = (status == 0) ? VAR0(ivar,ivar) : TEST;
        message("\n");
        tab_printg(" - Cov(h=0)  = ",1,GD_J_RIGHT,value);
        message("\n");
      }
      if (FLAG_VARZ)
      {
        value = (status == 0) ? get_ARRAY(DBOUT,IECH_OUT,IPTR_VARZ+ivar) : TEST;
        tab_printg(" - Var(Z*)   = ",1,GD_J_RIGHT,value);
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
static void st_result_simulate_print(int nbsimu,
                                     int nvar,
                                     int status)
{
  int    ivar,isimu,ecr;
  double value;

  /* Header */

  mestitle(0,"Simulation results");

  /* Loop on the results */

  for (isimu=ecr=0; isimu<nbsimu; isimu++)
    for (ivar = 0; ivar <nvar; ivar++,ecr++)
    {
      message("Simulation #%d of Z%-2d : ",isimu+1,ivar+1);
      value = (status == 0) ? get_ARRAY(DBOUT,IECH_OUT,IPTR_EST+ecr) : TEST;
      tab_printg(" = ",1,GD_J_RIGHT,value);
      message("\n");
    }
  return;
}

/****************************************************************************/
/*!
**  Print the neighborhood parameters
**
** \param[in]  status  Kriging error status
** \param[in]  ntab    Number of neighborhood parameters
** \param[in]  tab     Array of neighborhood parameters
**
*****************************************************************************/
static void st_res_nbgh_print(int     status,
                              int     ntab,
                              double *tab)
{
  if (status != 0) return;

  /* Header */

  mestitle(0,"Neighborhood Parameters");

  message ("Number of selected samples          = %d\n",(int) tab[0]);
  message ("Maximum neighborhood distance       = %lf\n",tab[1]);
  message ("Minimum neighborhood distance       = %lf\n",tab[2]);
  message ("Number of non-empty sectors         = %d\n",(int) tab[3]);
  message ("Number of consecutive empty sectors = %d\n",(int) tab[4]);

  return;
}

/****************************************************************************/
/*!
**  Create a grid which matches the grid containing the image to be
**  processed and is dimensionned to the neighborhood
**
** \return  The newly created db structure
**
** \param[in]  neigh  Neigh structure
** \param[in]  nvar   Number of variables 
**
*****************************************************************************/
static Db *st_image_build(Neigh *neigh,
                          int    nvar)

{
  int    *indg,error,ndim,i,nech,natt;
  double *coor,seuil,value;
  Db     *dbaux;
  VectorInt nx;
  VectorDouble tab;

  /* Initializations */

  dbaux = (Db  *) NULL;
  indg  = (int    *) NULL;
  coor  = (double *) NULL;
  error = 1;

  /* Preliminary checks */

  if (! is_grid(DBOUT)) goto label_end;
  ndim  = DBOUT->getNDim();
  natt  = ndim + nvar;
  seuil = 1. / neigh->getSkip();

  /* Core allocation */

  nx.resize(ndim);
  nech = 1;
  for (i=0; i<ndim; i++)
  {
    nx[i] = 2 * (int) neigh->getImageRadius(i) + 1;
    nech *= nx[i];
  }

  tab.resize(nech * natt);
  for (int i=0; i<nech * natt; i++) tab[i] = 0.;
  for (int iech=0; iech<nech; iech++) 
  {
    value = (law_uniform(0.,1.) < seuil) ? 0. : TEST;
    for (int ivar=0; ivar<nvar; ivar++)
      tab[ivar * nech + iech] = value;
  }

  /* Create the grid */

  dbaux = db_create_grid_generic(DBOUT->isGridRotated(),ndim,natt,
                                 LOAD_BY_COLUMN,1,nx,tab);

  /* Copy the grid characteristics */

  if (db_grid_copy_params(DBOUT,3,dbaux)) goto label_end;
  if (db_grid_copy_params(DBOUT,4,dbaux)) goto label_end;

  /* Set the locators */

  dbaux->setLocatorsByAttribute(nvar,0,LOC_Z);
  dbaux->setLocatorsByAttribute(ndim,nvar,LOC_X);
  if (db_grid_define_coordinates(dbaux)) goto label_end;

  /* Shift the origin */

  for (i=0; i<ndim; i++) dbaux->setX0(i,0.);
  indg = db_indg_alloc(dbaux);
  if (indg == (int *) NULL) goto label_end;
  coor = db_sample_alloc(dbaux,LOC_X);
  if (coor == (double *) NULL) goto label_end;
  db_index_sample_to_grid(dbaux,nech/2,indg);
  grid_to_point(dbaux,indg,(double *) NULL,coor);
  for (i=0; i<ndim; i++) dbaux->setX0(i, dbaux->getX0(i) -coor[i]);
  indg = db_indg_free(indg);
  coor = db_sample_free(coor);
  if (db_grid_define_coordinates(dbaux)) goto label_end;

  /* Set the error return status */

  error = 0;

label_end:
  if (error) dbaux = db_delete(dbaux);
  return(dbaux);
}

/****************************************************************************/
/*!
**  Calculate the kriging weights
**
** \return  Error return code
**
** \param[in]  model Model structure
** \param[in]  neigh Neigh structure
**
** \param[out] nred_out Number of covariance equations
** \param[out] neq_out  Number of equations per variable
**
*****************************************************************************/
static int st_image_kriging(Model *model,
                            Neigh *neigh,
                            int   *nred_out,
                            int   *neq_out)
{
  int status,nred,neq,nvar,nech,nfeq;
  double stdv;

  /* Initialization */

  nvar = model->getVariableNumber();
  nfeq = model->getDriftEquationNumber();

  /* Prepare the Koption structure */

  if (krige_koption_manage(1,1,KOPTION_PONCTUAL,1,VectorInt())) return(1);

  /* Prepare the neighborhood */

  neigh->setType(NEIGH_UNIQUE);
  (void) st_neigh(neigh,&status,&nech);
  neigh->setType(NEIGH_IMAGE);
  IECH_OUT = nech / 2;

  /* Establish the L.H.S. */

  st_prepar(model,neigh,nech,&status,&nred,&neq);
  if (status) return(1);

  /* Establish the kriging R.H.S. */

  st_rhs(model,nech,neq,nvar,NULL,&status);
  if (status) return(1);
  st_rhs_iso2hetero(neq,nvar);
  if (debug_query("kriging"))
    krige_rhs_print(nvar,nech,neq,nred,flag,rhs);

  /* Derive the kriging weights */

  matrix_product(nred,nred,nvar,lhs,rhs,wgt);
  if (debug_query("kriging"))
    krige_wgt_print(status,nvar,nvar,nfeq,nech,nred,-1,flag,wgt);

  /* Calculate the kriging variance */

  if (FLAG_STD)
  {
    st_variance0(model,nvar,VectorDouble());
    stdv = st_variance(model,0,0,nred);
    if (stdv < 0) stdv = 0.;
    stdv = sqrt(stdv);
  }

  (void) krige_koption_manage(-1,1,KOPTION_PONCTUAL,1,VectorInt());
  *nred_out = nred;
  *neq_out  = neq;
  return(0);
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
static int st_check_colcok(Db  *dbin,
                           Db  *dbout,
                           int *rank_colcok)
{
  int ivar,jvar;

  FLAG_COLK = 0;
  if (rank_colcok == (int *) NULL) return(0);

  /* Loop on the ranks of the colocated variables */

  for (ivar=0; ivar<dbin->getVariableNumber(); ivar++)
  {
    jvar = rank_colcok[ivar];
    if (IFFFF(jvar)) jvar = 0;
    if (jvar > dbout->getFieldNumber())
    {
      messerr("Error in the Colocation array:");
      messerr("Input variable (#%d): rank of the colocated variable is %d",
              ivar+1,jvar);
      messerr("But the Output file only contains %d attributes(s)",
              dbout->getFieldNumber());
      return(1);
    }
    rank_colcok[ivar] = jvar-1;
  }

  // Assign the array of ranks as a global variable
  RANK_COLCOK = rank_colcok;
  FLAG_COLK = 1;
  return(0);
}

/****************************************************************************/
/*!
**  Save the (Co-) Kriging weights using the keypair mechanism
**
** \param[in]  status   Kriging error status
** \param[in]  iech_out Rank of the output sample
** \param[in]  nvar     Number of variables
** \param[in]  nfeq     Number of drift equations
** \param[in]  nech     Number of active points
** \param[in]  nred     Reduced number of equations
** \param[in]  flag     Flag array
** \param[in]  wgt      Array of Kriging weights
**
*****************************************************************************/
static void st_save_keypair_weights(int     status,
                                    int     iech_out,
                                    int     nvar,
                                    int     nfeq,
                                    int     nech,
                                    int     nred,
                                    int    *flag,
                                    double *wgt)
{
  double wgtloc,values[5];
  int    lec,flag_value,iwgt,cumflag;

  /* Initializations */

  if (status != 0) return;
  values[0] = iech_out;

  /* Loop on the output variables */

  for (int jvar = lec = cumflag = 0; jvar < nvar; jvar++)
  {
    values[1] = jvar;
    
    /* Loop on the input samples */
    
    for (int iech = 0; iech < nech; iech++, lec++)
    {
      flag_value = (flag != (int *) NULL) ? flag[lec] : 1;
      if (flag_value)
      {
        values[2] = rank[iech];
        
        /* Loop on the input variables */
        
        for (int ivar = 0; ivar < nvar; ivar++)
        {
          iwgt   = nred * ivar + cumflag;
          wgtloc = (wgt != (double *) NULL && flag_value) ? wgt[iwgt] : TEST;
          if (! FFFF(wgtloc))
          {
            values[3] = ivar;
            values[4] = wgtloc;
            app_keypair("KrigingWeights",1,1,5,values);
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
** \param[in]  neigh       Neigh structure
** \param[in]  calcul      Kriging calculation option (::ENUM_KOPTIONS)
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
GEOSLIB_API int kriging(Db *dbin,
                        Db *dbout,
                        Model *model,
                        Neigh *neigh,
                        int calcul,
                        int flag_est,
                        int flag_std,
                        int flag_varz,
                        VectorInt ndisc,
                        VectorInt rank_colcok,
                        VectorDouble matCL,
                        NamingConvention namconv)
{
  int iext,inostat,error,status,nech,neq,nred,nvar,flag_new_nbgh,nfeq;
  int save_keypair;
  double ldum;

  /* Preliminary checks */

  error =  1;
  iext  = -1;
  nvar  =  0;
  st_global_init(dbin,dbout);
  save_keypair = (int) get_keypone("SaveKrigingWeights",0.);
  FLAG_EST  = flag_est;
  FLAG_STD  = flag_std;
  FLAG_VARZ = flag_varz;
  FLAG_WGT  = flag_std || flag_varz || save_keypair;
  if (st_check_colcok(dbin,dbout,rank_colcok.data())) goto label_end;
  if (st_check_environment(1,1,model,neigh)) goto label_end;
  if (manage_external_info(1,LOC_F,DBIN,DBOUT,&iext)) goto label_end;
  if (manage_external_info(1,LOC_NOSTAT,DBIN,DBOUT,
                           &inostat)) goto label_end;
  if (matCL.empty())
    nvar = model->getVariableNumber();
  else
    nvar = matCL.size() / model->getVariableNumber();
  nfeq = model->getDriftEquationNumber();
  nred = neq = 0;
  if (neigh->getType() == NEIGH_IMAGE)
  {
    messerr("This tool cannot function with an IMAGE neighborhood");
    goto label_end;
  }

  /* Add the attributes for storing the results */

  if (FLAG_EST)
  {
    IPTR_EST = dbout->addFields(nvar,0.);
    if (IPTR_EST < 0) goto label_end;
  }
  if (FLAG_STD)
  {
    IPTR_STD = dbout->addFields(nvar,0.);
    if (IPTR_STD < 0) goto label_end;
  }
  if (FLAG_VARZ)
  {
    IPTR_VARZ = dbout->addFields(nvar,0.);
    if (IPTR_VARZ < 0) goto label_end;
  }

  /* Pre-calculations */

  if (neigh_start(DBIN,neigh)) goto label_end;
  if (st_model_manage(1,model)) goto label_end;
  if (st_krige_manage(1,nvar,model,neigh)) goto label_end;
  if (krige_koption_manage(1,1,calcul,1,ndisc)) goto label_end;
  if (FLAG_STD) st_variance0(model,nvar,matCL);

  /* Loop on the targets to be processed */

  status = 0;
  for (IECH_OUT=0; IECH_OUT<get_NECH(DBOUT); IECH_OUT++)
  {
    mes_process("Kriging sample",get_NECH(DBOUT),IECH_OUT);
    debug_index(IECH_OUT+1);
    if (! dbout->isActive(IECH_OUT)) continue;
    if (debug_query("kriging") ||
        debug_query("nbgh")    ||
        debug_query("results"))
    {
      mestitle(1,"Target location");
      db_sample_print(dbout,IECH_OUT,1,0,0);
    }

    /* Select the Neighborhood */

    flag_new_nbgh = st_neigh(neigh,&status,&nech);
    if (status) goto label_store;

    /* Establish the kriging L.H.S. */

    if (flag_new_nbgh || flag_continuous_kriging(neigh) || debug_force())
    {
      st_prepar(model,neigh,nech,&status,&nred,&neq);
      if (status) goto label_store;
      st_data_dual(model,NULL,nech,nred,&ldum);
    }

    /* Establish the kriging R.H.S. */

    st_rhs(model,nech,neq,nvar,matCL.data(),&status);
    if (status) goto label_store;
    st_rhs_iso2hetero(neq,nvar);
    if (debug_query("kriging"))
      krige_rhs_print(nvar,nech,neq,nred,flag,rhs);

    /* Derive the kriging weights */

    if (FLAG_WGT)
    {
      matrix_product(nred,nred,nvar,lhs,rhs,wgt);
      if (debug_query("kriging"))
        krige_wgt_print(status,nvar,model->getVariableNumber(),nfeq,nech,nred,-1,flag,wgt);
    }

    /* Optional save of the Kriging/Cokriging weights */

    if (save_keypair)
      st_save_keypair_weights(status,IECH_OUT,model->getVariableNumber(),nfeq,nech,nred,
                              flag,wgt);

    /* Perform the estimation */

  label_store:
    st_estimate(model,NULL,status,neigh->getFlagXvalid(),nech,nvar,nred);
    if (debug_query("results"))
      st_result_kriging_print(neigh->getFlagXvalid(),nvar,status);
  }

  /* Set the error return flag */

  error = 0;
  namconv.setNamesAndLocators(dbin,LOC_Z,nvar,dbout,IPTR_VARZ,"varz");
  if (ABS(neigh->getFlagXvalid()) == 1)
    namconv.setNamesAndLocators(dbin,LOC_Z,-1,dbout,IPTR_STD,"stderr");
  else
    namconv.setNamesAndLocators(dbin,LOC_Z,-1,dbout,IPTR_STD,"stdev");
  if (ABS(neigh->getFlagXvalid()) == 1)
    namconv.setNamesAndLocators(dbin,LOC_Z,-1,dbout,IPTR_EST,"esterr",1,true);
  else
    namconv.setNamesAndLocators(dbin,LOC_Z,-1,dbout,IPTR_EST,"estim",1,true);

label_end:
  debug_index(0);
  (void) st_model_manage(-1,model);
  (void) st_krige_manage(-1,nvar,model,neigh);
  (void) krige_koption_manage(-1,1,calcul,1,ndisc);
  (void) manage_external_info(-1,LOC_F,DBIN,DBOUT,&iext);
  (void) manage_external_info(-1,LOC_NOSTAT,DBIN,DBOUT,&inostat);
  neigh_stop();
  return(error);
}

/****************************************************************************/
/*!
**  Perform the Cross-validation in Unique Neighborhood
**  This is a specific efficient algorithm
**
** \return  Error return code
**
** \param[in]  dbin        Input Db structure
** \param[in]  model       Model structure
** \param[in]  neigh       Neigh structrue
** \param[in]  flag_est    Option for the storing the estimation
** \param[in]  flag_std    Option for the storing the standard deviation
** \param[in]  rank_colcok Option for running Collocated Cokriging
** \param[in]  namconv     Naming Convention
**
*****************************************************************************/
static int st_xvalid_unique(Db *dbin,
                            Model *model,
                            Neigh *neigh,
                            int flag_est,
                            int flag_std,
                            VectorInt rank_colcok,
                            NamingConvention namconv)
{
  int iext,error,status,nech,neq,nred,nvar,flag_xvalid_memo,flag_new_nbgh;
  int iech,iiech,jech,jjech,inostat;
  double variance,value,stdv,valref;

  /* Preliminary checks */

  error  =  1;
  iext   = -1;
  nvar   =  0;
  st_global_init(dbin,dbin);
  FLAG_EST  = flag_est;
  FLAG_STD  = flag_std;
  FLAG_WGT  = flag_std;
  if (st_check_colcok(dbin,dbin,rank_colcok.data())) goto label_end;
  if (st_check_environment(1,1,model,neigh)) goto label_end;
  if (manage_external_info(1,LOC_F,DBIN,DBOUT,&iext)) goto label_end;
  if (manage_external_info(1,LOC_NOSTAT,DBIN,DBOUT,
                           &inostat)) goto label_end;
  nvar = model->getVariableNumber();

  /* Additional checks */

  if (neigh->getType() != NEIGH_UNIQUE)
  {
    messerr("This algorithm for Cross-Validation is dedicated to Unique Neighborhood");
    goto label_end;
  }
  if (FLAG_COLK)
  {
    messerr("The algorithm for Cross-Validation in Unique Neighborhood");
    messerr("does not allow the Colocated CoKriging Option");
    goto label_end;
  }
  if (nvar > 1)
  {
    messerr("The algorithm for Cross-Validation in Unique Neighborhood");
    messerr("is restricted to a single variable");
    goto label_end;
  }

  /* Add the attributes for storing the results */

  if (FLAG_EST)
  {
    IPTR_EST = dbin->addFields(nvar,0.);
    if (IPTR_EST < 0) goto label_end;
  }
  if (FLAG_STD)
  {
    IPTR_STD = dbin->addFields(nvar,0.);
    if (IPTR_STD < 0) goto label_end;
  }

  /* Pre-calculations */

  if (neigh_start(DBIN,neigh)) goto label_end;
  if (st_model_manage(1,model)) goto label_end;
  if (st_krige_manage(1,nvar,model,neigh)) goto label_end;
  if (krige_koption_manage(1,1,KOPTION_PONCTUAL,1,VectorInt())) goto label_end;

  /* Loop on the targets to be processed */

  status = 0;
  for (iech=iiech=0; iech<get_NECH(dbin); iech++)
  {
    IECH_OUT = iech;
    mes_process("Cross-Validation sample",get_NECH(dbin),iech);
    debug_index(iech+1);
    if (FLAG_EST) dbin->setArray(iech,IPTR_EST,TEST);
    if (FLAG_STD) dbin->setArray(iech,IPTR_STD,TEST);
    if (! dbin->isActive(iech)) continue;
    valref = dbin->getVariable(iech,0);
    if (FFFF(valref)) continue;
    if (debug_query("kriging") ||
        debug_query("nbgh")    ||
        debug_query("results"))
    {
      mestitle(1,"Target location");
      db_sample_print(dbin,iech,1,0,0);
    }

    /* Establish the Kriging matrix (without excluding any sample) */
    /* and invert it once for all */

    flag_xvalid_memo = neigh->getFlagXvalid();
    neigh->setFlagXvalid(0);
    flag_new_nbgh = st_neigh(neigh,&status,&nech);
    neigh->setFlagXvalid(flag_xvalid_memo);

    /* Establish the kriging L.H.S. */

    if (flag_new_nbgh || flag_continuous_kriging(neigh) || debug_force())
    {
      st_prepar(model,neigh,nech,&status,&nred,&neq);
      if (status) goto label_end;
    }

    /* Get the estimation variance */

    variance = 1. / LHS_C(iiech,iiech);
    stdv     = sqrt(variance);

    /* Perform the estimation */

    value = (ABS(neigh->getFlagXvalid()) == 1) ? -valref : 0.;
    for (jech=jjech=0; jech<get_NECH(dbin); jech++)
    {
      if (! dbin->isActive(jech)) continue;
      if (FFFF(dbin->getVariable(jech,0))) continue;
      if (iiech != jjech)
        value -= LHS_C(iiech,jjech) * variance * dbin->getVariable(jech,0);
      jjech++;
    }

    if (FLAG_EST) dbin->setArray(iech,IPTR_EST,value);

    if (ABS(neigh->getFlagXvalid()) == 1) stdv = value / stdv;
    if (FLAG_STD) dbin->setArray(iech,IPTR_STD,stdv);
    if (debug_query("results"))
      st_result_kriging_print(neigh->getFlagXvalid(),nvar,status);
    iiech++;
  }

  /* Set the error return flag */

  error = 0;
  if (ABS(neigh->getFlagXvalid()) == 1)
    namconv.setNamesAndLocators(dbin,LOC_Z,-1,dbin,IPTR_STD,"stderr");
  else
    namconv.setNamesAndLocators(dbin,LOC_Z,-1,dbin,IPTR_STD,"stdev");
  if (ABS(neigh->getFlagXvalid()) == 1)
    namconv.setNamesAndLocators(dbin,LOC_Z,-1,dbin,IPTR_EST,"esterr",1,true);
  else
    namconv.setNamesAndLocators(dbin,LOC_Z,-1,dbin,IPTR_EST,"estim",1,true);

label_end:
  debug_index(0);
  (void) st_model_manage(-1,model);
  (void) st_krige_manage(-1,nvar,model,neigh);
  (void) krige_koption_manage(-1,1,KOPTION_PONCTUAL,1,VectorInt());
  (void) manage_external_info(-1,LOC_F,DBIN,DBOUT,&iext);
  (void) manage_external_info(-1,LOC_NOSTAT,DBIN,DBOUT,&inostat);
  neigh_stop();
  return(error);
}

/****************************************************************************/
/*!
**  Standard Cross-Validation
**
** \return  Error return code
**
** \param[in]  db          Db structure
** \param[in]  model       Model structure
** \param[in]  neigh       Neigh structure
** \param[in]  flag_xvalid Type of output (see details)
** \param[in]  flag_code   1 if a code (K-FOLD) is used
** \param[in]  flag_est    Option for storing the estimation
** \param[in]  flag_std    Option for storing the standard deviation
** \param[in]  rank_colcok Option for running Collocated Cokriging
** \param[in]  namconv     Naming Convention
**
** \details When flag_xvalid = 1, the outputs are:
** \details - estimation: Z*-Z, st. dev: (Z*-Z)/S
** \details When flag_xvalid = 2, the outputs are:
** \details - estimation: Z*; st. dev: S
**
*****************************************************************************/
GEOSLIB_API int xvalid(Db    *db,
                       Model *model,
                       Neigh *neigh,
                       int flag_xvalid,
                       int flag_code,
                       int flag_est,
                       int flag_std,
                       VectorInt rank_colcok,
                       NamingConvention namconv)
{
  int ret_code;
  if (flag_code)
  {
    if (!db->hasCode())
      messerr("The K-FOLD option is ignored (no Code defined)");
    else if (neigh->getType() == NEIGH_UNIQUE)
      messerr("K-FOLD is not available in Unique Neighborhood");
    else
      neigh->setFlagXvalid(-flag_xvalid);
  }
  else
    neigh->setFlagXvalid(flag_xvalid);

  if (neigh->getType() == NEIGH_UNIQUE)
    ret_code = st_xvalid_unique(db, model, neigh, flag_est, flag_std,
                                rank_colcok, namconv);
  else
    ret_code = kriging(db, db, model, neigh, KOPTION_PONCTUAL, flag_est,
                       flag_std, 0, VectorInt(), rank_colcok, VectorDouble(),
                       namconv);
  neigh->setFlagXvalid(0);
  return ret_code;
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
** \param[in]  neigh       Neigh structure
** \param[in]  flag_est    Option for storing the estimation
** \param[in]  flag_std    Option for storing the standard deviation
** \param[in]  flag_varz   Option for storing the variance of the estimator
** \param[in]  rval        Change of support coefficient
**
*****************************************************************************/
GEOSLIB_API int krigdgm_f(Db     *dbin,
                        Db     *dbout,
                        Model  *model,
                        Neigh  *neigh,
                        int     flag_est,
                        int     flag_std,
                        int     flag_varz,
                        double  rval)
{
  int iext,inostat,error,status,nech,neq,nred,nvar,flag_new_nbgh,nfeq;
  int save_keypair;
  double ldum;

  /* Preliminary checks */

  error =  1;
  iext  = -1;
  nvar  =  0;
  st_global_init(dbin,dbout);
  save_keypair = (int) get_keypone("SaveKrigingWeights",0.);
  FLAG_EST  = flag_est;
  FLAG_STD  = flag_std;
  FLAG_VARZ = flag_varz;
  FLAG_WGT  = flag_std || flag_varz || save_keypair;
  FLAG_DGM  = 1;
  R_COEFF   = rval;
  if (st_check_environment(1,1,model,neigh)) goto label_end;
  if (manage_external_info(1,LOC_F,DBIN,DBOUT,&iext)) goto label_end;
  if (manage_external_info(1,LOC_NOSTAT,DBIN,DBOUT,
                           &inostat)) goto label_end;
  nvar = model->getVariableNumber();
  nfeq = model->getDriftEquationNumber();
  nred = neq = 0;
  if (neigh->getType() == NEIGH_IMAGE)
  {
    messerr("This tool cannot function with an IMAGE neighborhood");
    goto label_end;
  }

  /* Add the attributes for storing the results */

  if (FLAG_EST)
  {
    IPTR_EST = dbout->addFields(nvar,0.);
    if (IPTR_EST < 0) goto label_end;
  }
  if (FLAG_STD)
  {
    IPTR_STD = dbout->addFields(nvar,0.);
    if (IPTR_STD < 0) goto label_end;
  }
  if (FLAG_VARZ)
  {
    IPTR_VARZ = dbout->addFields(nvar,0.);
    if (IPTR_VARZ < 0) goto label_end;
  }

  /* Pre-calculations */

  if (neigh_start(DBIN,neigh)) goto label_end;
  if (st_model_manage(1,model)) goto label_end;
  if (st_krige_manage(1,nvar,model,neigh)) goto label_end;
  if (krige_koption_manage(1,1,KOPTION_PONCTUAL,1,VectorInt())) goto label_end;
  if (FLAG_STD) st_variance0(model,nvar,VectorDouble());

  /* Loop on the targets to be processed */

  status = 0;
  for (IECH_OUT=0; IECH_OUT<get_NECH(DBOUT); IECH_OUT++)
  {
    mes_process("Kriging sample",get_NECH(DBOUT),IECH_OUT);
    debug_index(IECH_OUT+1);
    if (! dbout->isActive(IECH_OUT)) continue;
    if (debug_query("kriging") ||
        debug_query("nbgh")    ||
        debug_query("results"))
    {
      mestitle(1,"Target location");
      db_sample_print(dbout,IECH_OUT,1,0,0);
    }

    /* Select the Neighborhood */

    flag_new_nbgh = st_neigh(neigh,&status,&nech);
    if (status) goto label_store;

    /* Establish the kriging L.H.S. */

    if (flag_new_nbgh || flag_continuous_kriging(neigh) || debug_force())
    {
      st_prepar(model,neigh,nech,&status,&nred,&neq);
      if (status) goto label_store;
      st_data_dual(model,NULL,nech,nred,&ldum);
    }

    /* Establish the kriging R.H.S. */

    st_rhs(model,nech,neq,nvar,NULL,&status);
    if (status) goto label_store;
    st_rhs_iso2hetero(neq,nvar);
    if (debug_query("kriging"))
      krige_rhs_print(nvar,nech,neq,nred,flag,rhs);

    /* Derive the kriging weights */

    if (FLAG_WGT)
    {
      matrix_product(nred,nred,nvar,lhs,rhs,wgt);
      if (debug_query("kriging"))
        krige_wgt_print(status,nvar,model->getVariableNumber(),nfeq,nech,nred,-1,flag,wgt);
    }

    /* Optional save of the Kriging/Cokriging weights */

    if (save_keypair)
      st_save_keypair_weights(status,IECH_OUT,model->getVariableNumber(),nfeq,nech,nred,
                              flag,wgt);

    /* Perform the estimation */

  label_store:
    st_estimate(model,NULL,status,neigh->getFlagXvalid(),nech,nvar,nred);
    if (debug_query("results"))
      st_result_kriging_print(neigh->getFlagXvalid(),nvar,status);
  }

  /* Set the error return flag */

  error = 0;

label_end:
  debug_index(0);
  (void) st_model_manage(-1,model);
  (void) st_krige_manage(-1,nvar,model,neigh);
  (void) krige_koption_manage(-1,1,KOPTION_PONCTUAL,1,VectorInt());
  (void) manage_external_info(-1,LOC_F,DBIN,DBOUT,&iext);
  (void) manage_external_info(-1,LOC_NOSTAT,DBIN,DBOUT,&inostat);
  neigh_stop();
  return(error);
}

/****************************************************************************/
/*!
**  Ponctual Kriging based on profiles
**
** \return  Error return code
**
** \param[in]  dbin      input Db structure
** \param[in]  dbout     output Db structure
** \param[in]  model     Model structure
** \param[in]  neigh     Neigh structrue
** \param[in]  ncode     Number (maximum) of codes expected
** \param[in]  flag_est  Option for the storing the estimation
** \param[in]  flag_std  Option for the storing the standard deviation
**
*****************************************************************************/
GEOSLIB_API int krigprof_f(Db    *dbin,
                         Db    *dbout,
                         Model *model,
                         Neigh *neigh,
                         int    ncode,
                         int    flag_est,
                         int    flag_std)
{
  int iext,inostat,status,nech,neq,nred,nvar,flag_new_nbgh,nfeq,iptr_dat,icode;
  int error;
  double ldum;

  /* Preliminary checks */

  if (! dbin->isVariableNumberComparedTo(1))
  {
    messerr("This method is restricted to the monovariate case");
    return(1);
  }
  if (! dbin->hasCode() || dbin->getVarianceErrorNumber() != 1)
  {
    messerr("This method requires variables CODE and V to be defined");
    return(1);
  }

  /* Preliminary checks */

  error =  1;
  iext  = iptr_dat = -1;
  st_global_init(dbin,dbout);
  FLAG_EST  = flag_est;
  FLAG_STD  = flag_std;
  FLAG_WGT  = flag_std;
  FLAG_PROF = 1;
  nvar  = dbin->getVariableNumber();
  if (st_check_environment(1,1,model,neigh)) goto label_end;
  if (manage_external_info(1,LOC_F,DBIN,DBOUT,&iext)) goto label_end;
  if (manage_external_info(1,LOC_NOSTAT,DBIN,DBOUT,
                           &inostat)) goto label_end;

  /* Add the attributes for storing the results */

  if (FLAG_EST)
  {
    IPTR_EST = dbout->addFields(ncode,0.);
    if (IPTR_EST < 0) goto label_end;
  }
  if (FLAG_STD)
  {
    IPTR_STD = dbout->addFields(ncode,0.);
    if (IPTR_STD < 0) goto label_end;
  }

  /* Preliminary transforms */

  iptr_dat = DBIN->addFields(ncode,TEST);
  if (iptr_dat < 0) goto label_end;
  nfeq  = model->getDriftEquationNumber();

  /* Pre-calculations */

  if (neigh_start(DBIN,neigh)) goto label_end;
  if (st_model_manage(1,model)) goto label_end;
  if (st_krige_manage(1,nvar,model,neigh)) goto label_end;
  if (krige_koption_manage(1,1,KOPTION_PONCTUAL,1,VectorInt())) goto label_end;
  if (FLAG_STD) st_variance0(model,nvar,VectorDouble());

  /* Loop on the targets to be processed */

  status = 0;
  for (IECH_OUT=0; IECH_OUT<get_NECH(DBOUT); IECH_OUT++)
  {
    mes_process("Kriging sample",get_NECH(DBOUT),IECH_OUT);
    debug_index(IECH_OUT+1);
    if (! dbout->isActive(IECH_OUT)) continue;
    if (debug_query("kriging") ||
        debug_query("nbgh")    ||
        debug_query("results"))
    {
      mestitle(1,"Target location");
      db_sample_print(dbout,IECH_OUT,1,0,0);
    }

    /* Select the Neighborhood */

    flag_new_nbgh = st_neigh(neigh,&status,&nech);
    if (status) goto label_store;

    /* Establish the kriging L.H.S. */

    if (flag_new_nbgh || flag_continuous_kriging(neigh) || debug_force())
    {
      st_prepar(model,neigh,nech,&status,&nred,&neq);
      if (status) goto label_store;
      st_data_dual(model,NULL,nech,nred,&ldum);
    }

    /* Establish the kriging R.H.S. */

    st_rhs(model,nech,neq,nvar,NULL,&status);
    if (status) goto label_store;
    st_rhs_iso2hetero(neq,nvar);
    if (debug_query("kriging"))
      krige_rhs_print(nvar,nech,neq,nred,flag,rhs);

    /* Derive the kriging weights */

    if (FLAG_WGT)
    {
      matrix_product(nred,nred,nvar,lhs,rhs,wgt);
      if (debug_query("kriging"))
        krige_wgt_print(status,nvar,nvar,nfeq,nech,nred,-1,flag,wgt);
    }

    /* Perform the estimation */

  label_store:
    st_estimate(model,NULL,status,neigh->getFlagXvalid(),nech,nvar,nred);
    if (debug_query("results"))
      st_result_kriging_print(neigh->getFlagXvalid(),nvar,status);
  }

  /* Set the error return flag */

  error = 0;

label_end:
  debug_index(0);
  (void) st_model_manage(-1,model);
  (void) st_krige_manage(-1,nvar,model,neigh);
  (void) krige_koption_manage(-1,1,KOPTION_PONCTUAL,1,VectorInt());
  (void) manage_external_info(-1,LOC_F,DBIN,DBOUT,&iext);
  (void) manage_external_info(-1,LOC_NOSTAT,DBIN,DBOUT,&inostat);
  if (iptr_dat >= 0)
    for (icode=0; icode<ncode; icode++)
      dbin->deleteFieldByAttribute(iptr_dat+icode);
  neigh_stop();
  return(error);
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
static int bayes_manage(int      mode,
                        int      nbsimu,
                        Model   *model,
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

    *rmean = st_core(nfeq,1);
    if (*rmean == (double *) NULL) return(1);
    *rcov  = st_core(nfeq,nfeq);
    if (*rcov  == (double *) NULL) return(1);
    if (nbsimu > 0)
    {
      *smean = st_core(nfeq,nbsimu);
      if (*smean == (double *) NULL) return(1);
    }
  }
  else
  {

    /* Deallocation */

    *rmean = (double *) mem_free((char *) (*rmean));
    *rcov  = (double *) mem_free((char *) (*rcov));
    if (nbsimu > 0) *smean = (double *) mem_free((char *) (*smean));
  }
  return(0);
}

/****************************************************************************/
/*!
**  Perform the Bayesian estimation of the Drift Coefficients
**
** \return  Error return code
**
** \param[in]  model     Model structure
** \param[in]  neigh     Neigh structrue
** \param[in]  dmean     Array giving the prior means for the drift terms
** \param[in]  dcov      Array containing the prior covariance matrix
**                       for the drift terms
**
** \param[out] rmean     Array giving the posterior means for the drift terms
** \param[out] rcov      Array containing the  posterior covariance matrix
**                       for the drift terms
**
*****************************************************************************/
static int bayes_precalc(Model  *model,
                         Neigh  *neigh,
                         double *dmean,
                         double *dcov,
                         double *rmean,
                         double *rcov)
{
  int nfeq,error,status,nech,nred,neq,shift,ib,jb,il,jl,flag_fix;
  double *ff,*smu,*sigma,*vars,value;

  /* Initializations */

  error  = 1;
  nfeq   = model->getDriftEquationNumber();
  ff     = smu = sigma = vars = (double *) NULL;

  /* Preliminary Checks */

  if (neigh->getType() != NEIGH_UNIQUE)
  {
    messerr("The Bayesian Estimation of the Drift Coefficients");
    messerr("is only available in Unique Neighborhood");
    goto label_end;
  }

  /* Core allocation */

  IECH_OUT = get_NECH(DBIN) / 2;

  /* Check that the variance-covariance matrix is symmetric */

  flag_fix = is_matrix_null(nfeq,nfeq,dcov,0);
  if (! is_matrix_symmetric(nfeq,dcov,1)) goto label_end;

  /* Prepare the Kriging matrix (without correction) */

  (void) st_neigh(neigh,&status,&nech);
  FLAG_BAYES = 0;
  st_prepar(model,neigh,nech,&status,&nred,&neq);
  FLAG_BAYES = 1;
  if (status) goto label_end;
  shift = nred - nfeq;

  /* Complementary core allocation */

  ff      = st_core(shift,nfeq);
  if (ff      == (double *) NULL) goto label_end;
  smu     = st_core(shift,1);
  if (smu     == (double *) NULL) goto label_end;
  sigma   = st_core(shift,shift);
  if (sigma   == (double *) NULL) goto label_end;
  vars    = st_core(shift,1);
  if (vars    == (double *) NULL) goto label_end;

  // Create the array of variables

  ib = 0;
  for (int iech = 0; iech < get_NECH(DBIN); iech++)
  {
    if (! DBIN->isActive(iech)) continue;
    for (int ivar = 0; ivar < DBIN->getVariableNumber(); ivar++)
    {
      double value = DBIN->getVariable(rank[iech],ivar);
      if (FFFF(value)) continue;
      vars[ib++] = value;
    }
  }

  /* Copy DCOV into S and DMEAN into RMEAN */

  (void) memcpy((char *) rcov ,(char *) dcov ,sizeof(double) * nfeq * nfeq);
  (void) memcpy((char *) rmean,(char *) dmean,sizeof(double) * nfeq);
  if (flag_fix) goto label_print;
  
  /* Establish the drift array FF */

  for (il=0; il<nfeq; il++)
    for (ib=0; ib<shift; ib++)
      FF(ib,il) = LHS_B(ib,shift+il);

  /* Calculate S-1 */

  if (matrix_invert(rcov,nfeq,-1)) goto label_end;

  /* Calculate: SMU = S-1 * MEAN */

  matrix_product(nfeq,nfeq,1,rcov,dmean,smu);
    
  /* Covariance matrix SIGMA */

  for (ib=0; ib<shift; ib++)
    for (jb=0; jb<shift; jb++)
      SIGMA(ib,jb) = LHS_B(ib,jb);
  
  /* Calculate SIGMA-1 */

  if (matrix_invert(sigma,shift,-1)) goto label_end;

  /* Inverse of posterior covariance matrix: SC-1 = FFt * SIGMA-1 * FF + S-1 */

  for (il=0; il<nfeq; il++)
    for (jl=0; jl<nfeq; jl++)
    {
      value = 0.;
      for (ib=0; ib<shift; ib++)
        for (jb=0; jb<shift; jb++)
          value += FF(ib,il) * SIGMA(ib,jb) * FF(jb,jl);
      RCOV(il,jl) += value;
    }

  /* Calculating: SMU = FFt * SIGMA-1 * Z + S-1 * MU */

  for (il=0; il<nfeq; il++)
  {
    value = 0.;
    for (ib=0; ib<shift; ib++)
      for (jb=0; jb<shift; jb++)
        value += FF(ib,il) * SIGMA(ib,jb) * vars[jb];
    SMU(il) += value;
  }

  /* Posterior mean: RMEAN = SC * SMU */

  if (matrix_invert(rcov,nfeq,-1)) goto label_end;
  matrix_product(nfeq,nfeq,1,rcov,smu,rmean);

label_print:
  if (debug_query("bayes"))
  {
    mestitle (0,"Bayesian Drift coefficients");
    print_matrix("Prior Mean",0,1,nfeq,1,NULL,dmean);
    print_matrix("Prior Variance-Covariance",0,1,nfeq,nfeq,NULL,dcov);

    print_matrix("Posterior Mean",0,1,nfeq,1,NULL,rmean);
    print_matrix("Posterior Variance-Covariance",0,1,nfeq,nfeq,NULL,rcov);
    message("\n");
  }

  /* Set the error return code */

  error = 0;
  IECH_NBGH = -1;

label_end:

  /* Core deallocation */

  ff      = (double *) mem_free((char *) ff);
  smu     = (double *) mem_free((char *) smu);
  sigma   = (double *) mem_free((char *) sigma);
  vars    = (double *) mem_free((char *) vars);

  return(error);
}

/****************************************************************************/
/*!
**  Correct the arrays RHS and VARB in Bayesian case
**
** \return  Error return code
**
** \param[in]  model  Model structure
** \param[in]  rcov   Array containing the posterior covariance matrix
**                    for the drift terms
** \param[in]  neq    Number of equations
** \param[in]  nred   Reduced number of equations
**
** \param[out] status Returned status
**
*****************************************************************************/
static void st_bayes_correct(Model  *model,
                             double *rcov,
                             int     neq,
                             int     nred,
                             int    *status)
{
  int ivar,jvar,il,jl,nvar,nfeq;

  /* Initializations */

  nvar = model->getVariableNumber();
  nfeq = model->getDriftEquationNumber();

  /* Establish the Drift matrix */

  st_ff0(model,status);
  if (*status) return;

  /* Correct the arrays */

  for (ivar=0; ivar<nvar; ivar++)
    for (jvar=0; jvar<nvar; jvar++)
    {
      VARB(ivar,jvar) = 0.;
      for (il=0; il<nfeq; il++)
        for (jl=0; jl<nfeq; jl++)
          VARB(ivar,jvar) += FF0(il,ivar) * RCOV(il,jl) * FF0(jl,jvar);
    }
  return;
}

/****************************************************************************/
/*!
**  Estimation with Bayesian Drift
**
** \return  Error return code
**
** \param[in]  dbin      input Db structure
** \param[in]  dbout     output Db structure
** \param[in]  model     Model structure
** \param[in]  neigh     Neigh structrue
** \param[in]  dmean     Array giving the prior means for the drift terms
** \param[in]  dcov      Array containing the prior covariance matrix
**                       for the drift terms
** \param[in]  flag_est  Pointer for the storing the estimation
** \param[in]  flag_std  Pointer for the storing the standard deviation
**
*****************************************************************************/
GEOSLIB_API int kribayes_f(Db *dbin,
                           Db *dbout,
                           Model *model,
                           Neigh *neigh,
                           double *dmean,
                           double *dcov,
                           int flag_est,
                           int flag_std)
{
  int      iext,error,status,nech,neq,nred,nvar;
  int      flag_new_nbgh,inostat;
  double  *rmean,*rcov,*smean,ldum;
  Model   *model_sk;

  /* Preliminary checks */

  error      =  1;
  iext       = -1;
  nvar       =  0;
  model_sk   = (Model   *) NULL;
  rmean      = smean = rcov = (double *) NULL;
  st_global_init(dbin,dbout);
  FLAG_BAYES = 1;
  FLAG_EST   = flag_est;
  FLAG_STD   = flag_std;
  FLAG_WGT   = flag_std;
  if (st_check_environment(1,1,model,neigh)) goto label_end;
  if (manage_external_info(1,LOC_F,DBIN,DBOUT,&iext)) goto label_end;
  if (manage_external_info(1,LOC_NOSTAT,DBIN,DBOUT,
                           &inostat)) goto label_end;
  nvar = model->getVariableNumber();

  /* Add the attributes for storing the results */

  if (FLAG_EST)
  {
    IPTR_EST = dbout->addFields(nvar,0.);
    if (IPTR_EST < 0) goto label_end;
  }
  if (FLAG_STD)
  {
    IPTR_STD = dbout->addFields(nvar,0.);
    if (IPTR_STD < 0) goto label_end;
  }

  /* Complementary core allocation */

  if (bayes_manage(1,0,model,&rmean,&rcov,&smean)) goto label_end;

  /* Pre-calculations */

  if (neigh_start(DBIN,neigh)) goto label_end;
  if (st_model_manage(1,model)) goto label_end;
  if (st_krige_manage(1,nvar,model,neigh)) goto label_end;
  if (krige_koption_manage(1,1,KOPTION_PONCTUAL,1,VectorInt())) goto label_end;
  if (FLAG_STD) st_variance0(model,nvar,VectorDouble());

  /* Solve the Bayesian estimation of the Drift coefficients */

  if (bayes_precalc(model,neigh,dmean,dcov,rmean,rcov)) goto label_end;

  /* Duplicate the model, suppressing the Drift terms */

  model_sk = model_duplicate(model,0.,-1);
  if (model_sk == (Model *) NULL) goto label_end;

  /* Loop on the targets to be processed */

  status = 0;
  for (IECH_OUT=0; IECH_OUT<get_NECH(DBOUT); IECH_OUT++)
  {
    mes_process("Bayesian estimation",get_NECH(DBOUT),IECH_OUT);
    debug_index(IECH_OUT+1);
    if (! dbout->isActive(IECH_OUT)) continue;
    if (debug_query("kriging") ||
        debug_query("nbgh")    ||
        debug_query("results"))
    {
      mestitle(1,"Target location");
      db_sample_print(dbout,IECH_OUT,1,0,0);
    }

    /* Select the Neighborhood */

    flag_new_nbgh = st_neigh(neigh,&status,&nech);
    if (status) goto label_store;

    /* Establish the kriging L.H.S. */
 
    if (flag_new_nbgh || flag_continuous_kriging(neigh) || debug_force())
    {
      st_prepar(model_sk,neigh,nech,&status,&nred,&neq);
      if (status) goto label_store;
      // We must use the drift initial assumption, hence model (not model_sk)
      st_data_dual(model,rmean,nech,nred,&ldum);
    }

    /* Establish the kriging R.H.S. */

    st_rhs(model_sk,nech,neq,nvar,NULL,&status);
    if (status) goto label_store;
    st_rhs_iso2hetero(neq,nvar);

    /* Modify the arrays in the Bayesian case */

    st_bayes_correct(model,rcov,nred,neq,&status);
    if (status) goto label_store;
    if (debug_query("kriging"))
      krige_rhs_print(nvar,nech,neq,nred,flag,rhs);

    /* Derive the kriging weights */

    if (FLAG_WGT)
    {
      matrix_product(nred,nred,nvar,lhs,rhs,wgt);
      if (debug_query("kriging"))
        krige_wgt_print(status,nvar,nvar,0,nech,nred,-1,flag,wgt);
    }

    /* Perform the estimation */

  label_store:
    // We must use the drift initial assumption, hence model (not model_sk)
    st_estimate(model,rmean,status,neigh->getFlagXvalid(),nech,nvar,nred);
    if (debug_query("results"))
      st_result_kriging_print(neigh->getFlagXvalid(),nvar,status);
  }

  /* Set the error return flag */

  error = 0;

label_end:
  debug_index(0);
  model_sk = model_free(model_sk);
  (void) bayes_manage(-1,0,model,&rmean,&rcov,&smean);
  (void) st_model_manage(-1,model);
  (void) st_krige_manage(-1,nvar,model,neigh);
  (void) krige_koption_manage(-1,1,KOPTION_PONCTUAL,1,VectorInt());
  (void) manage_external_info(-1,LOC_F,DBIN,DBOUT,&iext);
  (void) manage_external_info(-1,LOC_NOSTAT,DBIN,DBOUT,&inostat);
  neigh_stop();
  return(error);
}

/****************************************************************************/
/*!
**  Check the Neighborhood
**
** \return  Error return code
**
** \param[in]  dbin      input Db structure
** \param[in]  dbout     output Db structure
** \param[in]  model     Model structure (optional)
** \param[in]  neigh     Neigh structure
** \param[in]  namconv   Naming Convention
**
** \remark This procedure creates the following arrays:
** \remark 1 - The number of selected samples
** \remark 2 - The maximum neighborhood distance
** \remark 3 - The minimum neighborhood distance
** \remark 4 - The number of non-empty sectors
** \remark 5 - The number of consecutive empty sectors
**
*****************************************************************************/
GEOSLIB_API int test_neigh(Db    *dbin,
                           Db    *dbout,
                           Model *model,
                           Neigh *neigh,
                           NamingConvention namconv)
{
  int error,status,nech,ntab,iext,inostat;
  double tab[5];

  /* Preliminary checks */

  error = 1;
  iext  = -1;
  ntab  = 5;
  st_global_init(dbin,dbout);
  if (st_check_environment(1,1,model,neigh)) goto label_end;
  if (manage_external_info(1,LOC_F,DBIN,DBOUT,&iext)) goto label_end;
  if (manage_external_info(1,LOC_NOSTAT,DBIN,DBOUT,
                           &inostat)) goto label_end;

  /* Add the attributes for storing the results */

  IPTR_NBGH = dbout->addFields(ntab,0.);
  if (IPTR_NBGH < 0) goto label_end;

  /* Pre-calculations */

  if (neigh_start(DBIN,neigh)) goto label_end;
  if (st_model_manage(1,model)) goto label_end;
  if (st_krige_manage(1,model->getVariableNumber(),model,neigh)) goto label_end;

  /* Loop on the targets to be processed */

  status = 0;
  for (IECH_OUT=0; IECH_OUT<get_NECH(DBOUT); IECH_OUT++)
  {
    mes_process("Neighborhood Test",get_NECH(DBOUT),IECH_OUT);
    debug_index(IECH_OUT+1);
    if (! dbout->isActive(IECH_OUT)) continue;
    if (debug_query("kriging") ||
        debug_query("nbgh")    ||
        debug_query("results"))
    {
      mestitle(1,"Target location");
      db_sample_print(dbout,IECH_OUT,1,0,0);
    }

    /* Select the Neighborhood */

    (void) st_neigh(neigh,&status,&nech);

    /* Retrieve the neighborhood parameters */

    neigh_echo(dbin,neigh,rank,nech,tab);

    /* Store the neighborhood parameters */

    st_store_nbgh(status,ntab,tab);
    if (debug_query("nbgh")) st_res_nbgh_print(status,ntab,tab);
  }

  /* Set the error return flag */

  error = 0;
  namconv.setNamesAndLocators(NULL,LOC_UNKNOWN,1,dbout,IPTR_NBGH  ,"Number");
  namconv.setNamesAndLocators(NULL,LOC_UNKNOWN,1,dbout,IPTR_NBGH+1,"MaxDist");
  namconv.setNamesAndLocators(NULL,LOC_UNKNOWN,1,dbout,IPTR_NBGH+2,"MinDist");
  namconv.setNamesAndLocators(NULL,LOC_UNKNOWN,1,dbout,IPTR_NBGH+3,"NbNESect");
  namconv.setNamesAndLocators(NULL,LOC_UNKNOWN,1,dbout,IPTR_NBGH+4,"NbCESect");
  namconv.setLocators(dbout,IPTR_NBGH,ntab);

label_end:
  debug_index(0);
  (void) st_model_manage(-1,model);
  (void) st_krige_manage(-1,model->getVariableNumber(),model,neigh);
  (void) manage_external_info(-1,LOC_F,DBIN,DBOUT,&iext);
  (void) manage_external_info(-1,LOC_NOSTAT,DBIN,DBOUT,&inostat);
  neigh_stop();
  return(error);
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
** \param[in]  neigh      Neigh structrue
** \param[in]  dmean      Array giving the prior means for the drift terms
** \param[in]  dcov       Array containing the prior covariance matrix
**                        for the drift terms
** \param[in]  icase      Case for PGS and GRF (or -1)
** \param[in]  nbsimu     Number of simulations
** \param[in]  flag_dgm   1 if the DGM version of kriging should be used
** \param[in]  rval       Change of support coefficient
**
*****************************************************************************/
GEOSLIB_API int krigsim(const char *strloc,
                        Db     *dbin,
                        Db     *dbout,
                        Model  *model,
                        Neigh  *neigh,
                        double *dmean,
                        double *dcov,
                        int     icase,
                        int     nbsimu,
                        int     flag_dgm,
                        double  rval)
{
  int      error,status,nech,neq,nred,nvar,flag_new_nbgh,iext,nfeq,inostat;
  double  *rmean,*rcov,*smean;
  Model   *model_sk;

  /* Preliminary checks */

  error      =  1;
  iext       = -1;
  nvar       =  0;
  model_sk   = (Model   *) NULL;
  rmean      = smean = rcov = (double *) NULL;
  st_global_init(dbin,dbout);
  FLAG_BAYES = (dmean != (double *) NULL && dcov  != (double *) NULL);
  FLAG_EST   = 1;
  FLAG_STD   = 0;
  FLAG_WGT   = 1;
  FLAG_SIMU  = 1;
  FLAG_DGM   = flag_dgm;
  R_COEFF    = rval;
  IPTR_EST   = dbout->getColumnByLocator(LOC_SIMU,0);
  if (st_check_environment(1,1,model,neigh)) goto label_end;
  if (manage_external_info(1,LOC_F,DBIN,DBOUT,&iext)) goto label_end;
  if (manage_external_info(1,LOC_NOSTAT,DBIN,DBOUT,
                           &inostat)) goto label_end;
  nvar = model->getVariableNumber();
  nfeq = model->getDriftEquationNumber();

  /* Core allocation */

  if (FLAG_BAYES)
  {
    if (bayes_manage(1,nbsimu,model,&rmean,&rcov,&smean)) goto label_end;
  }

  /* Pre-calculations */

  if (neigh_start(DBIN,neigh)) goto label_end;
  if (st_model_manage(1,model)) goto label_end;
  if (st_krige_manage(1,nvar,model,neigh)) goto label_end;
  if (krige_koption_manage(1,1,KOPTION_PONCTUAL,1,VectorInt())) goto label_end;
  if (FLAG_STD) st_variance0(model,nvar,VectorDouble());

  /* Solve the Bayesian estimation of the Drift coefficients */

  if (FLAG_BAYES)
  {
    if (bayes_precalc(model,neigh,dmean,dcov,rmean,rcov)) goto label_end;
  }

  /* Simulate the drift coefficients from the posterior distributions */

  if (FLAG_BAYES)
  {
    if (bayes_simulate(model,nbsimu,rmean,rcov,smean)) goto label_end;
  }

  /* Duplicate the model, suppressing the Drift terms */

  if (FLAG_BAYES)
  {
    model_sk = model_duplicate(model,0.,-1);
    if (model_sk == (Model *) NULL) goto label_end;
  }
  else
  {
    model_sk = model;
  }

  /* Loop on the targets to be processed */

  status = 0;
  for (IECH_OUT=0; IECH_OUT<get_NECH(DBOUT); IECH_OUT++)
  {
    mes_process(strloc,get_NECH(DBOUT),IECH_OUT);
    debug_index(IECH_OUT+1);
    if (! dbout->isActive(IECH_OUT)) continue;
    if (debug_query("kriging") ||
        debug_query("nbgh")    ||
        debug_query("results"))
    {
      mestitle(1,"Target location");
      db_sample_print(dbout,IECH_OUT,1,0,0);
    }

    /* Select the Neighborhood */

    flag_new_nbgh = st_neigh(neigh,&status,&nech);
    if (status) goto label_store;

    /* Establish the kriging L.H.S. */

    if (flag_new_nbgh || flag_continuous_kriging(neigh) || debug_force())
    {
      st_prepar(model_sk,neigh,nech,&status,&nred,&neq);
      if (status) goto label_store;
    }

    /* Establish the kriging R.H.S. */

    st_rhs(model_sk,nech,neq,nvar,NULL,&status);
    if (status) goto label_store;
    st_rhs_iso2hetero(neq,nvar);

    /* Modify the arrays in the Bayesian case */

    if (FLAG_BAYES)
    {
      st_bayes_correct(model,dcov,nred,neq,&status);
      if (status) goto label_store;
    }
    if (debug_query("kriging"))
      krige_rhs_print(nvar,nech,neq,nred,flag,rhs);

    /* Derive the kriging weights */

    if (FLAG_WGT)
    {
      matrix_product(nred,nred,nvar,lhs,rhs,wgt);
      if (debug_query("kriging"))
        krige_wgt_print(status,nvar,nvar,nfeq,nech,nred,icase,flag,wgt);
    }

    /* Perform the simulation */

  label_store:
    st_simulate(model,smean,status,icase,nbsimu,nech,nred);
    if (debug_query("results"))
      st_result_simulate_print(nbsimu,nvar,status);
  }

  /* Set the error return flag */

  error = 0;

label_end:
  debug_index(0);
  if (FLAG_BAYES)
  {
    model_sk = model_free(model_sk);
    (void) bayes_manage(-1,nbsimu,model,&rmean,&rcov,&smean);
  }
  (void) st_model_manage(-1,model);
  (void) st_krige_manage(-1,nvar,model,neigh);
  (void) krige_koption_manage(-1,1,KOPTION_PONCTUAL,1,VectorInt());
  (void) manage_external_info(-1,LOC_F,DBIN,DBOUT,&iext);
  (void) manage_external_info(-1,LOC_NOSTAT,DBIN,DBOUT,&inostat);
  neigh_stop();
  return(error);
}

/****************************************************************************/
/*!
**  Kriging (Factorial) a regular grid
**
** \return  Error return code
**
** \param[in]  dbgrid    input and output Db grid structure
** \param[in]  model     Model structure
** \param[in]  neigh     Neigh structure
**
*****************************************************************************/
GEOSLIB_API int krimage_func(Db *dbgrid, Model *model, Neigh *neigh)
{
  int    i,iech,jech,error,nvar,nfeq,nb_neigh,ecr,ndim,nred,neq;
  int   *indn0,*indnl,*indg0,*indgl;
  double data,estim;
  Db    *dbaux;

  /* Preliminary checks */

  error = 1;
  dbaux = (Db *) NULL;
  indn0 = indnl = indg0 = indgl = (int *) NULL;
  st_global_init(dbgrid,dbgrid);
  nvar  = model->getVariableNumber();
  nfeq  = model->getDriftEquationNumber();
  ndim  = dbgrid->getNDim();
  FLAG_EST  = 1;
  FLAG_STD  = 1;
  FLAG_WGT  = 1;
  if (st_check_environment(1,1,model,neigh)) goto label_end;

  /* Add the attributes for storing the results */

  IPTR_EST = dbgrid->addFields(nvar,0.);
  if (IPTR_EST < 0) goto label_end;

  /* Core allocation */

  indg0 = db_indg_alloc(dbgrid);
  indgl = db_indg_alloc(dbgrid);

  /* Create the secondary grid for image processing */

  dbaux = st_image_build(neigh,nvar);
  if (dbaux == (Db *) NULL) goto label_end;

  nb_neigh = get_NECH(dbaux);
  indn0 = db_indg_alloc(dbaux);
  indnl = db_indg_alloc(dbaux);
  db_index_sample_to_grid(dbaux,nb_neigh/2,indn0);

  /* Pre-calculations */

  if (neigh_start(dbaux,neigh)) goto label_end;
  if (st_model_manage(1,model)) goto label_end;
  if (st_krige_manage(1,nvar,model,neigh)) goto label_end;

  /* Establish the kriging weights */

  DBOUT = DBIN = dbaux;
  if (st_image_kriging(model,neigh,&nred,&neq)) goto label_end;

  /* Loop on the targets to be processed */

  DBIN   = dbaux;
  DBOUT  = dbgrid;
  for (IECH_OUT=0; IECH_OUT<get_NECH(dbgrid); IECH_OUT++)
  {
    mes_process("Image filtering",get_NECH(DBOUT),IECH_OUT);
    debug_index(IECH_OUT+1);
    if (! DBOUT->isActive(IECH_OUT)) continue;
    db_index_sample_to_grid(DBOUT,IECH_OUT,indg0);

    /* Loop on the target variables */

    for (int ivar=0; ivar<nvar; ivar++)
    {

      /* Loop on the neighboring points */

      estim = 0.;
      if (nfeq <= 0) estim = model->getMean(ivar);
      
      /* Loop on the explanatory variables */

      ecr = ivar * neq;
      for (int jvar=0; jvar<nvar; jvar++)
      {
        for (iech=0; iech<nb_neigh; iech++)
        {
          if (FFFF(DBIN->getVariable(iech,0))) continue;
          db_index_sample_to_grid(DBIN,iech,indnl);
          for (i=0; i<ndim; i++)
          {
            indgl[i] = indg0[i] - indn0[i] + indnl[i];
            indgl[i] = get_mirror_sample(DBOUT->getNX(i),indgl[i]);
          }
          jech = db_index_grid_to_sample(DBOUT,indgl);
          data = DBOUT->getVariable(jech,jvar);
          if (FFFF(data))
          {
            estim = TEST;
            goto label_store;
          }
          if (nfeq <= 0) data -= model->getMean(jvar);
          estim += data * wgt[ecr];
          ecr++;
        }
      }

    label_store:
      DBOUT->setArray(IECH_OUT,IPTR_EST+ivar,estim);
    }
  }

  /* Set the error return flag */

  error = 0;

label_end:
  debug_index(0);
  (void) st_model_manage(-1,model);
  (void) st_krige_manage(-1,nvar,model,neigh);
  dbaux = db_delete(dbaux);
  indn0 = db_indg_free(indn0);
  indnl = db_indg_free(indnl);
  indg0 = db_indg_free(indg0);
  indgl = db_indg_free(indgl);
  neigh_stop();
  return(error);
}

/****************************************************************************/
/*!
**  Global estimation over a territory using arithmetic average
**
** \param[in]  dbin          Db structure
** \param[in]  dbgrid        Db grid structure
** \param[in]  model         Model structure
** \param[in]  ivar          Rank of the target variable
** \param[in]  flag_verbose  1 for a verbose output
** \param[in]  seed          Seed for the random number generator
** \param[in]  surface       surface of the integration polygon
**                           (all the grid if not defined)
**
** \param[out]  zest   global arithmetic average
** \param[out]  sse    global estimation standard devation
** \param[out]  cvgeo  CV geo
**
** \remark  This function erases any Weight variable already defined
**
*****************************************************************************/
GEOSLIB_API int global_arithmetic(Db     *dbin,
                                  Db     *dbgrid,
                                  Model  *model,
                                  int     ivar,
                                  int     flag_verbose,
                                  int     seed,
                                  double  surface,
                                  double *zest,
                                  double *sse,
                                  double *cvgeo)
{
  double cvv,cxv,cxx,wtot,ave,var,cvsam,cviid,mini,maxi;
  int    error,np,ng,iatt;

  /* Initializations */

  error = 1;
  st_global_init(dbin,dbgrid);

  /* Preliminary checks */

  if (st_check_environment(1,1,model,NULL)) goto label_end;
  if (ivar < 0 || ivar >= dbin->getVariableNumber())
  {
    messerr("The target variable (%d) must lie between 1 and the number of variables (%d)",
            ivar+1,dbin->getVariableNumber());
    goto label_end;
  }
  
  /* Preliminary assignments */

  np = dbin->getActiveSampleNumber();
  ng = dbgrid->getActiveSampleNumber();
  if (FFFF(surface)) surface = ng * db_grid_maille(dbgrid);
  if (seed != 0) law_set_random_seed(seed);

  /* Average covariance over the data */

  cxx = model_cxx(model,dbin,dbin,ivar,ivar,0,0.);

  /* Average covariance between the data and the territory */

  cxv = model_cxx(model,dbin,dbgrid,ivar,ivar,0,0.);

  /* Average covariance over the territory */

  cvv = model_cxx(model,dbgrid,dbgrid,ivar,ivar,0,db_epsilon_distance(dbgrid));

  /* Auxiliary calculations */

  iatt = db_attribute_identify(dbin,LOC_Z,ivar);
  db_monostat(dbin,iatt,&wtot,&ave,&var,&mini,&maxi);

  /* Auxiliary calculations */

  *zest  = ave;
  *sse   = cvv - 2. * cxv + cxx;
  *sse   = (*sse > 0) ? sqrt(*sse) : 0.;
  cvsam  = (ave != 0.) ? sqrt(var)  / ave : TEST;
  cviid  = cvsam      / sqrt(np);
  *cvgeo = (ave != 0.) ? (*sse) / ave : TEST;

  if (flag_verbose)
  {
    message("Global estimation by arithmetic average\n");
    message("=======================================\n");
    message("Total number of data             = %d\n" ,get_NECH(dbin));
    message("Number of active data            = %d\n" ,np);
    message("Sample variance                  = %lf\n",var);
    message("CVsample                         = %lf\n",cvsam);
    message("CViid                            = %lf\n",cviid);
    message("Cxx                              = %lf\n",cxx);
    message("Cxv                              = %lf\n",cxv);
    message("Cvv                              = %lf\n",cvv);
    if (FFFF(ave))
      message("Estimation by arithmetic average = NA\n");
    else
      message("Estimation by arithmetic average = %lf\n",ave);
    message("Estimation St. dev. of the mean  = %lf\n",*sse);
    if (FFFF(*cvgeo))
      message("CVgeo                            = NA\n");
    else
      message("CVgeo                            = %lf\n",*cvgeo);
    message("\n");
    message("Surface                          = %lf\n",surface);
    if (FFFF(ave))
      message("Q (Estimation * Surface)         = NA\n");
    else
      message("Q (Estimation * Surface)         = %lf\n",ave * surface);
    message("\n");
  }

  /* Set the error return code */

  error = 0;

label_end:
  return(error);
}

/****************************************************************************/
/*!
**  Global estimation over a territory using kriging
**
** \param[in]  dbin          Db structure
** \param[in]  dbout         Db grid structure
** \param[in]  model         Model structure
** \param[in]  ivar          Rank of the target variable
** \param[in]  flag_verbose  1 for a verbose output
** \param[in]  calcul        KOPTION_PONCTUAL or KOPTION_DRIFT
** \param[in]  seed          Seed for the random number generator
** \param[in]  surface       surface of the integration polygon
**                           (all the grid if not defined)
**
** \param[out]  zest     global estimation
** \param[out]  sse      global standard deviation
** \param[out]  cvgeo    CV geo
** \param[out]  weights  Array of weights attached to data
**                       (Dimension: nvar * nech)
**
*****************************************************************************/
GEOSLIB_API int global_kriging(Db     *dbin,
                               Db     *dbout,
                               Model  *model,
                               int     ivar,
                               int     flag_verbose,
                               int     calcul,
                               int     seed,
                               double  surface,
                               double *zest,
                               double *sse,
                               double *cvgeo,
                               double *weights)
{
  double  *rhs_tot,estim,stdv,cvv,ldum;
  int      error,i,np,ng,size,nbfl,status,nech,nred,neq,nfeq,nvar;
  int      flag_new_nbgh,ntot,lec,jvar;
  Neigh   *neigh;

  /* Initializations */

  error   = 1;
  nvar    = 0;
  rhs_tot = (double  *) NULL;
  st_global_init(dbin,dbout);
  neigh = neigh_init_unique(dbin->getNDim());
  if (st_check_environment(1,1,model,neigh)) goto label_end;
  if (ivar < 0 || ivar >= dbin->getVariableNumber())
  {
    messerr("The target variable (%d) must lie between 1 and the number of variables (%d)",
            ivar+1,dbin->getVariableNumber());
    goto label_end;
  }
  FLAG_EST  = 1;
  FLAG_STD  = 1;
  FLAG_WGT  = 1;

  /* Preliminary checks */

  np   = dbin->getActiveSampleNumber();
  ng   = dbout->getActiveSampleNumber();
  nbfl = model->getDriftNumber();
  nfeq = model->getDriftEquationNumber();
  nvar = model->getVariableNumber();
  ntot = get_NECH(dbin);
  if (FFFF(surface)) surface = ng * db_grid_maille(dbout);
  if (seed != 0) law_set_random_seed(seed);

  /* Pre-calculations */

  if (neigh_start(dbin,neigh)) goto label_end;
  if (st_model_manage(1,model)) goto label_end;
  if (st_krige_manage(1,nvar,model,neigh)) goto label_end;
  if (krige_koption_manage(1,1,calcul,1,VectorInt())) goto label_end;

  /* Core allocation */

  size = (np + nbfl) * nvar * nvar;
  rhs_tot = (double *) mem_alloc(sizeof(double) * size,0);
  if (rhs_tot == (double *) NULL) goto label_end;
  for (i=0; i<size; i++) rhs_tot[i] = 0.;

  /* Average covariance over the territory */

  cvv = model_cxx(model,dbout,dbout,ivar,ivar,0,db_epsilon_distance(dbin));

  /* Loop on the targets to be processed */

  status = 0;
  for (IECH_OUT=0; IECH_OUT<get_NECH(DBOUT); IECH_OUT++)
  {
    mes_process("Global estimation",get_NECH(DBOUT),IECH_OUT);
    debug_index(IECH_OUT+1);
    if (! dbout->isActive(IECH_OUT)) continue;
    if (debug_query("nbgh")    ||
        debug_query("results"))
    {
      mestitle(1,"Target location");
      db_sample_print(dbout,IECH_OUT,1,0,0);
    }

    /* Select the Neighborhood */

    flag_new_nbgh = st_neigh(neigh,&status,&nech);
    if (status) goto label_store;

    /* Establish the kriging L.H.S. */

    if (flag_new_nbgh || flag_continuous_kriging(neigh) || debug_force())
    {
      st_prepar(model,neigh,nech,&status,&nred,&neq);
      if (status) goto label_store;
      st_data_dual(model,NULL,nech,nred,&ldum);
    }

    /* Establish the kriging R.H.S. */

    st_rhs(model,nech,neq,nvar,NULL,&status);
    if (status) goto label_store;

    /* Cumulate the R.H.S */

    for (i=0; i<size; i++) rhs_tot[i] += rhs[i];
  }

label_store:
  if (status || ng <= 0)
  {
    estim = stdv = TEST;
  }
  else
  {

    /* Load the scaled cumulated R.H.S. in the array rhs */

    for (i=0; i<size; i++) rhs[i] = rhs_tot[i] / (double) ng;
    st_rhs_iso2hetero(neq,nvar);
    if (debug_query("kriging"))
      krige_rhs_print(nvar,nech,neq,nred,flag,rhs);

    /* Derive the kriging weights */

    matrix_product(nred,nred,nvar,lhs,rhs,wgt);
    if (debug_query("kriging"))
      krige_wgt_print(status,nvar,nvar,nfeq,nech,nred,-1,flag,wgt);

    /* Perform the estimation */

    estim = 0.;
    for (i=0; i<nred; i++) estim += RHS_C(i,ivar) * ZAM1(i);
    stdv  = cvv;
    for (i=0; i<nred; i++) stdv  -= RHS_C(i,ivar) * WGT(i,ivar);
  }

  /* Auxiliary calculations */

  *zest  = estim;
  *sse   = (stdv > 0) ? sqrt(stdv) : 0.;
  *cvgeo = (estim == 0. || FFFF(estim)) ? TEST : (*sse) / estim;

  /* Store the weights */

  for (i=0; i<nvar * ntot; i++) weights[i] = TEST;
  for (jvar=lec=0; jvar<nvar; jvar++)
    for (i=0; i<np; i++) 
      if (flag[jvar*ntot+i] && rank[i]>=0)
        weights[ntot*jvar + rank[i]] = wgt[lec++];

  /* Printout */

  if (flag_verbose)
  {
    message("Global estimation kriging\n");
    message("=========================\n");
    message("Total number of data             = %d\n" ,ntot);
    message("Number of active data            = %d\n" ,np);
    message("Number of variables              = %d\n" ,nvar);
    message("Cvv                              = %lf\n",cvv);
    if (FFFF(estim))
      message("Estimation by kriging            = NA\n");
    else
    {
      message("Estimation by kriging            = %lf\n",estim);
      if (nbfl > 0)
        for (i=0; i<nbfl; i++)
          message("Lagrange Parameter #%d            = %lf\n",i+1,wgt[np+i]);
    }
    message("Estimation St. Dev. of the mean  = %lf\n",*sse);
    if (FFFF(*cvgeo))
      message("CVgeo                            = NA\n");
    else
      message("CVgeo                            = %lf\n",*cvgeo);
    message("\n");
    message("Surface                          = %lf\n",surface);
    if (FFFF(estim))
      message("Q (Estimation * Surface)         = NA\n");
    else
      message("Q (Estimation * Surface)         = %lf\n",estim * surface);
    message("\n");
  }

  /* Set the error return code */

  error = 0;

label_end:
  debug_index(0);
  (void) st_model_manage(-1,model);
  (void) st_krige_manage(-1,nvar,model,neigh);
  (void) krige_koption_manage(-1,1,calcul,1,VectorInt());
  neigh_stop();
  rhs_tot = (double *) mem_free((char *) rhs_tot);
  neigh = neigh_free(neigh);
  return(error);
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
GEOSLIB_API int global_transitive(Db     *dbgrid,
                                  Model  *model,
                                  int     flag_verbose,
                                  int     flag_regular,
                                  int     ndisc,
                                  double *abundance,
                                  double *sse,
                                  double *cvtrans)
{
  int     idim,ndim,i,ix,iy,ix1,ix2,iy1,iy2,nx,ny,error,flag_value;
  double  c00,cvv,dx,dy,dsum,gint,dsse,wtot;
  double  value;
  VectorDouble d1;
  CovCalcMode mode;

  /* Initializations */

  error =  1;
  cvv = wtot = dsse = gint = dsum = 0.;
  flag_value = 0;
  st_global_init(dbgrid,dbgrid);
  if (st_check_environment(0,1,model,NULL)) goto label_end;
  ndim  = dbgrid->getNDim();
  d1.resize(2);

  if (ndim < 1 || ndim > 2)
  {
    messerr("The transitive global estimation is implemented for 1 and 2 space only");
    goto label_end;
  }
  if (model->getVariableNumber() != 1)
  {
    messerr("The transitive global estimation is implemented for 1 variable only");
    goto label_end;
  }
  if (! is_grid(dbgrid))
  {
    messerr("The transitive global estimation requires a Db organized as a grid");
    goto label_end;
  }

  /* Core allocation */

  for (idim=0; idim<ndim; idim++) d1[idim] = 0.;
  model_calcul_cov(model,mode,1,1.,d1,&c00);

  /* Abundance estimation */

  flag_value = 0;
  if (dbgrid->getVariableNumber() == 1)
  {
    for (i=0; i<get_NECH(dbgrid); i++)
    {
      value = dbgrid->getVariable(i,0);
      if (! FFFF(value)) dsum += value;
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

      for (ix= -nx+1; ix<=nx; ix++)
        for (iy= -ny+1; iy<=ny; iy++)
        {
          d1[0] = dx * ix;
          d1[1] = dy * iy;
          model_calcul_cov(model,mode,0,1.,d1,&dsse);
        }
      dsse *= dx * dy;
      // TODO : appeler model_integral
      // if (model_integral(model,ndisc,&gint)) goto label_end;
      *sse = dsse - gint;
    }
    else
    {

      /* Stratified case */

      for (ix1=0; ix1<ndisc; ix1++)
        for (iy1=0; iy1<ndisc; iy1++)
          for (ix2=0; ix2<ndisc; ix2++)
            for (iy2=0; iy2<ndisc; iy2++)
            {
              d1[0] = dx * (ix2 - ix1) / ndisc;
              d1[1] = dy * (iy2 - iy1) / ndisc;
              model_calcul_cov(model,mode,0,1.,d1,&cvv);
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

      for (ix= -nx+1; ix<=nx; ix++)
      {
        d1[0] = dx * ix;
        model_calcul_cov(model,mode,0,1.,d1,&dsse);
      }
      dsse *= dx;
      // TODO: appeler model_integral
      // if (model_integral(model,ndisc,&gint)) goto label_end;
      *sse = dsse - gint;
    }
    else
    {

      /* Stratified case */

      for (ix1=0; ix1<ndisc; ix1++)
        for (ix2=0; ix2<ndisc; ix2++)
        {
          d1[0] = dx * (ix2 - ix1) / ndisc;
          model_calcul_cov(model,mode,0,1.,d1,&cvv);
          wtot += 1.;
        }
      cvv /= wtot;
      *sse = dx * (c00 - cvv);
    }
  }

  if (flag_value)
  {
    *abundance = dsum;
    *cvtrans   = ((*sse) <= 0.) ? TEST : dsum / (*sse);
  }
  else
  {
    *abundance = *cvtrans = TEST;
  }
  (*sse) = (*sse > 0) ? sqrt(*sse) : 0.;

  /* Optional printout */

  if (flag_verbose)
  {
    if (flag_regular)
    {
      message("Transitive estimation (Regular case)\n");
      message("====================================\n");
      message("Space dimension           = %d \n",ndim);
      message("s * Sum[G(ks)]            = %lf\n",dsse);
      message("Integral[G(h)]            = %lf\n",gint);
    }
    else
    {
      message("Transitive estimation (Stratified case)\n");
      message("=======================================\n");
      message("Space dimension           = %d \n",ndim);
      message("G(0)                      = %lf\n",c00);
      message("G(s,s)                    = %lf\n",cvv);
    }
    message("Estimation St. Dev.       = %lf\n",(*sse));
    if (flag_value)
    {
      message("Global abundance          = %lf\n",(*abundance));
      if (FFFF(*cvtrans))
        message("Coefficient of Variation  = NA\n");
      else
        message("Coefficient of Variation  = %lf\n",(*cvtrans));
    }
  }

  /* Set the error return code */

  error = 0;

label_end:
  return(error);
}

/****************************************************************************/
/*!
**  Inverse distance estimation when Input DB is a grid
**
** \param[in]  exponent    exponent of the inverse distance
** \param[in]  flag_expand 1 for expansion option
**
** \param[out] indg        Working array
** \param[out] indref      Working array
** \param[out] coor        Working array
** \param[out] cooref      Working array
**
*****************************************************************************/
static void st_grid_invdist(int     exponent,
                            int     flag_expand,
                            int    *indg,
                            int    *indref,
                            double *coor,
                            double *cooref)
{
  int    idim,ndim,maxneigh,incorrect,ind,rank,iech_neigh;
  double result,total,val_neigh,dist,wgt,dmin;

  /* Initializations */

  ndim  = DBIN->getNDim();
  maxneigh = (int) pow(2., (double) ndim);
  (void) db_extension_diag(DBOUT,&dmin);
  dmin /= 1.e5;

  /* Loop on the targets to be processed */

  for (IECH_OUT=0; IECH_OUT<get_NECH(DBOUT); IECH_OUT++)
  {
    mes_process("Estimation by Inverse distance",get_NECH(DBOUT),IECH_OUT);
    if (! DBOUT->isActive(IECH_OUT)) continue;
    if (debug_query("kriging") ||
        debug_query("nbgh")    ||
        debug_query("results"))
    {
      mestitle(1,"Target location");
      db_sample_print(DBOUT,IECH_OUT,1,0,0);
    }

    /* Find the grid index corresponding to the target */

    for (idim=0; idim<ndim; idim++)
      cooref[idim] = get_IDIM(DBOUT,IECH_OUT,idim);
    point_to_grid(DBIN,cooref,flag_expand,indref);

    /* Loop on the neighbors */

    result = total = 0.;
    for (rank=0; rank<maxneigh; rank++)
    {
      for (idim=0; idim<ndim; idim++) indg[idim] = indref[idim];

      /* Decompose the neighborhood rank */

      idim = 0;
      ind  = rank;
      while (ind > 0)
      {
        if (ind % 2 == 1) indg[idim] += 1;
        ind /= 2;
        idim++;
      }

      /* Check that the neighboring point lies within the grid */

      for (idim=incorrect=0; idim<ndim && incorrect==0; idim++)
      {
        if (indg[idim] >= DBIN->getNX(idim))
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

        iech_neigh = db_index_grid_to_sample(DBIN,indg);
        val_neigh  = DBIN->getVariable(iech_neigh,0);
        if (FFFF(val_neigh))
        {
          result = TEST;
          break;
        }
        else
        {

          /* Calculate the distance from neighborhood to target */

          dist = 0.;
          grid_to_point(DBIN,indg,NULL,coor);
          dist = ut_distance(ndim,cooref,coor);
          if (dist < dmin)
          {
            result = val_neigh;
            total  = 1.;
            break;
          }
          wgt     = 1. / pow(dist,exponent);
          result += wgt * val_neigh;
          total  += wgt;
        }
      }
    }
    if (! FFFF(result)) result /= total;
    DBOUT->setArray(IECH_OUT,IPTR_EST,result);
  }

  return;
}

/****************************************************************************/
/*!
**  Inverse distance estimation when Input DB is a point file
**
** \param[in]  exponent    exponent of the inverse distance
** \param[in]  flag_expand 1 for expansion option
** \param[in]  dmax        Maximum search radius (only used for Point Db)
**
** \param[out] coor        Working array
** \param[out] cooref      Working array
**
*****************************************************************************/
static void st_point_invdist(int     exponent,
                             int     flag_expand,
                             double  dmax,
                             double *coor,
                             double *cooref)
{
  int    idim,iech_in,ndim;
  double result,total,val_neigh,dist,wgt,dmin;

  /* Initializations */

  ndim  = DBIN->getNDim();
  (void) db_extension_diag(DBOUT,&dmin);
  dmin /= 1.e5;

  /* Loop on the targets to be processed */

  for (IECH_OUT=0; IECH_OUT<get_NECH(DBOUT); IECH_OUT++)
  {
    mes_process("Estimation by Inverse distance",get_NECH(DBOUT),IECH_OUT);
    if (! DBOUT->isActive(IECH_OUT)) continue;
    if (debug_query("kriging") ||
        debug_query("nbgh")    ||
        debug_query("results"))
    {
      mestitle(1,"Target location");
      db_sample_print(DBOUT,IECH_OUT,1,0,0);
    }
    for (idim=0; idim<ndim; idim++)
      cooref[idim] = get_IDIM(DBOUT,IECH_OUT,idim);

    /* Loop on the data points */

    result = total = 0.;
    for (iech_in=0; iech_in<get_NECH(DBIN); iech_in++)
    {
      if (! DBIN->isActive(iech_in)) continue;
      for (idim=0; idim<ndim; idim++)
        coor[idim] = get_IDIM(DBIN,iech_in,idim);
      val_neigh = DBIN->getVariable(iech_in,0);
      if (FFFF(val_neigh)) continue;

      /* Check that the data point is a valid neighbor */

      dist = ut_distance(ndim,coor,cooref);
      if (! FFFF(dmax) && dist > dmax) continue;

      /* Process the new neighboring point */

      if (dist < dmin)
      {
        result = val_neigh;
        total  = 1.;
        break;
      }
      wgt     = 1. / pow(dist,exponent);
      result += wgt * val_neigh;
      total  += wgt;
    }
    if (! FFFF(result)) result /= total;
    DBOUT->setArray(IECH_OUT,IPTR_EST,result);
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
GEOSLIB_API int invdist_f(Db    *dbin,
                        Db    *dbout,
                        int    exponent,
                        int    flag_expand,
                        double dmax)
{
  int    *indg,*indref,error;
  double *coor,*cooref;

  /* Initializations */

  error = 1;
  indg = indref = (int *) NULL;
  coor = cooref = (double *) NULL;
  st_global_init(dbin,dbout);
  if (st_check_environment(1,1,NULL,NULL)) goto label_end;

  /* Add the attribute for storing the result */

  IPTR_EST = dbout->addFields(1,0.);
  if (IPTR_EST < 0) goto label_end;
  coor   = db_sample_alloc(dbout,LOC_X);
  if (coor   == (double *) NULL) goto label_end;
  cooref = db_sample_alloc(dbout,LOC_X);
  if (cooref == (double *) NULL) goto label_end;

  if (! is_grid(DBIN))
  {
    st_point_invdist(exponent,flag_expand,dmax,coor,cooref);
  }
  else
  {
    indg   = db_indg_alloc(dbin);
    if (indg   == (int *) NULL) goto label_end;
    indref = db_indg_alloc(dbin);
    if (indref == (int *) NULL) goto label_end;
    st_grid_invdist(exponent,flag_expand,indg,indref,coor,cooref);
  }

  /* Set the error return code */

  error = 0;

label_end:
  coor   = db_sample_free(coor);
  cooref = db_sample_free(cooref);
  if (is_grid(DBIN))
  {
    indg   = db_indg_free(indg);
    indref = db_indg_free(indref);
  }
  return(error);
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
static int st_get_limits(Db     *db,
                         double  top,
                         double  bot,
                         int    *ideb,
                         int    *ifin)
{
  int    ndim,nz,iad;
  double z0,dz;

  /* Initializations */

  ndim = db->getNDim();
  z0   = db->getX0(ndim-1);
  nz   = db->getNX(ndim-1);
  dz   = db->getDX(ndim-1);

  /* Preliminary checks */

  if (! FFFF(bot) && ! FFFF(top) && top < bot)
  {
    messerr("Error: Top(%lf) must be larger than Bottom (%lf)",top,bot);
    return(1);
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
  return(0);
}

/****************************************************************************/
/*!
**  Definition of the neighborhood
**
** \return  Error return code: 1 if the target does not belong to the
** \return  area of interest
**
** \param[in]  nz            Number of cells along Z
** \param[in]  ideb          Index of the starting sample
** \param[in]  ifin          Index of the ending sample
** \param[in]  neigh_radius  Radius of the Neighborhood
**
** \param[out] status        Neighborhood error status
** \param[out] nbefore       Number of samples in neighborhood before target
** \param[out] nafter        Number of samples in neighborhood after target
**
*****************************************************************************/
static int st_get_neigh(int  nz,
                        int  ideb,
                        int  ifin,
                        int  neigh_radius,
                        int *status,
                        int *nbefore,
                        int *nafter)
{
  int iad;

  *status = 1;
  if (IECH_OUT < ideb || IECH_OUT > ifin) return(1);

  iad = MAX(IECH_OUT - neigh_radius, ideb);
  *nbefore = IECH_OUT - iad;

  iad = MIN(IECH_OUT + neigh_radius, ifin);
  *nafter  = iad - IECH_OUT;

  *status = 0;
  return(0);
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
static double st_cov_exp(int     dist,
                         double *cov,
                         int     cov_radius,
                         int     flag_sym)
{
  double val1,val2,val;

  if (flag_sym)
  {
    val1 = cov[cov_radius - dist];
    val2 = cov[cov_radius + dist];
    val  = (val1 + val2) / 2.;
  }
  else
  {
    val  = cov[cov_radius + dist];
  }
  return(val);
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
                       int     cov_radius,
                       int     flag_sym,
                       int     nfeq,
                       int     nbefore,
                       int     nafter,
                       int     neq)
{
  int i,j;

  /* Covariance part */

  for (i= -nbefore; i<=nafter; i++)
    for (j= -nbefore; j<=nafter; j++)
    {
      LHS_EXP(i+nbefore,j+nbefore) =
        st_cov_exp(i-j,covdd,cov_radius,flag_sym);
      LHS_EXP(j+nbefore,i+nbefore) =
        st_cov_exp(j-i,covdd,cov_radius,flag_sym);
    }

  /* Drift part */

  if (nfeq == 0) return;
  for (i= -nbefore; i<=nafter; i++)
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
                       int     cov_radius,
                       int     flag_sym,
                       int     nfeq,
                       int     nbefore,
                       int     nafter,
                       int     neq)
{
  int i;

  /* Covariance part */

  for (i= -nbefore; i<=nafter; i++)
    RHS_EXP(i+nbefore) = st_cov_exp(i,covd0,cov_radius,flag_sym);

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
static double st_estim_exp(Db     *db,
                           double *wgt,
                           int     nbefore,
                           int     nafter)
{
  int i;
  double result;

  /* Perform the estimation */

  result = 0.;
  for (i= -nbefore; i<=nafter; i++)
    result += wgt[i+nbefore] * db->getVariable(IECH_OUT+i,0);

  return(result);
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
GEOSLIB_API int anakexp_f(Db     *db,
                        double *covdd,
                        double *covd0,
                        double  top,
                        double  bot,
                        int     cov_radius,
                        int     neigh_radius,
                        int     flag_sym,
                        int     nfeq)
{
  int      i,ndim,nvarin,nech,size,error,ideb,ifin,neq,status;
  int      nbefore,nafter,nbefore_mem,nafter_mem;
  double   result;

  /* Initializations */

  error = 1;
  st_global_init(db,db);
  FLAG_EST = 1;
  lhs    = rhs = wgt = (double *) NULL;
  ndim   = db->getNDim();
  nvarin = db->getVariableNumber();
  nbefore_mem = nafter_mem = -1;
  size = nech = 0;

  /* Prepare the Koption structure */

  if (krige_koption_manage(1,1,KOPTION_PONCTUAL,1,VectorInt())) return(1);

  /* Preliminary checks */

  if (ndim != 1 || ! is_grid(db))
  {
    messerr("This procedure is limited to 1-D grid");
    goto label_end;
  }
  if (nvarin != 1)
  {
    messerr("This procedure is limited to the monovariate case");
    goto label_end;
  }
  nech = db->getNX(ndim-1);
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
    messerr("to the radius of the covariance (%d)",cov_radius);
    goto label_end;
  }

  /* Add the attribute for storing the result */

  IPTR_EST = db->addFields(nvarin,0.);
  if (IPTR_EST < 0) goto label_end;
  DBOUT = db;

  /* Core allocation */

  size = 2 * neigh_radius + 1;
  st_krige_manage_basic(1,size,size,1,nfeq);
  for (i=0; i<nech; i++)
  {
    rank[i] = i;
    flag[i] = 1;
  }

  /* Get the limits of the area to be processed */

  if (st_get_limits(db,top,bot,&ideb,&ifin)) goto label_end;

  /* Loop on the grid nodes */

  status = 0;
  for (IECH_OUT=0; IECH_OUT<nech; IECH_OUT++)
  {
    mes_process("Factorial Kriging Analysis",nech,IECH_OUT);
    debug_index(IECH_OUT+1);
    if (! db->isActive(IECH_OUT)) continue;
    if (debug_query("kriging") ||
        debug_query("nbgh")    ||
        debug_query("results"))
    {
      mestitle(1,"Target location");
      db_sample_print(db,IECH_OUT,1,0,0);
    }

    /* Discard the grid nodes which doe not belong to the processed area */

    DBOUT->setArray(IECH_OUT,IPTR_EST,TEST);

    /* Look for the neighborhood */

    if (st_get_neigh(nech,ideb,ifin,neigh_radius,
                     &status,&nbefore,&nafter)) continue;

    /* If the neighborhood has changed, establish the kriging system */

    neq = nafter + nbefore + 1;
    if (nfeq == 1) neq++;
    if (nbefore_mem != nbefore || nafter_mem  != nafter || debug_force())
    {
      nbefore_mem = nbefore;
      nafter_mem  = nafter;

      /* Establish the L.H.S. of the kriging system */

      st_lhs_exp(covdd,cov_radius,flag_sym,nfeq,nbefore,nafter,neq);
      if (debug_query("kriging")) krige_lhs_print(nech,neq,neq,flag,lhs);

      /* Invert the kriging system */

      if (matrix_invert(lhs,neq,IECH_OUT))
      {
        status = 1;
        continue;
      }

      /* Establish the R.H.S. of the kriging system */

      st_rhs_exp(covd0,cov_radius,flag_sym,nfeq,nbefore,nafter,neq);
      if (debug_query("kriging"))
        krige_rhs_print(nvarin,nech,neq,neq,flag,rhs);

      /* Derive the kriging weights */

      matrix_product(neq,neq,1,lhs,rhs,wgt);
    }

    /* Calculate the estimation */

    result = st_estim_exp(db,wgt,nbefore,nafter);
    DBOUT->setArray(IECH_OUT,IPTR_EST,result);
    if (debug_query("results")) st_result_kriging_print(0,nvarin,status);
  }

  /* Set the error return flag */

  error = 0;

label_end:
  debug_index(0);
  (void) krige_koption_manage(-1,1,KOPTION_PONCTUAL,1,VectorInt());
  st_krige_manage_basic(-1,size,size,1,nfeq);
  return(error);
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
static void st_calculate_covres(Db     *db,
                                Model  *model,
                                double *cov_ref,
                                int     cov_radius,
                                int     flag_sym,
                                int     cov_ss[3],
                                int     cov_nn[3],
                                double *cov_res)
{
  double dx,dy,c00,covtot,covtab,covver;
  int    i,ix,iy,iz;
  VectorDouble d1;
  CovCalcMode mode;

  /* Initializations */

  d1.resize(3);
  dx = db->getDX(0);
  dy = db->getDX(1);
  covtot = COV_REF(0);
  for (i=0; i<3; i++) d1[i] = 0.;
  model_calcul_cov(model,mode,1,1.,d1,&c00);

  /* Evaluate the array of experimental covariance of the residual variable */

  for (ix=-cov_nn[0]; ix<=cov_nn[0]; ix++)
    for (iy=-cov_nn[1]; iy<=cov_nn[1]; iy++)
      for (iz=-cov_nn[2]; iz<=cov_nn[2]; iz++)
      {
        if (! flag_sym)
          covver = COV_REF(iz);
        else
          covver = (COV_REF(iz) + COV_REF(-iz)) / 2.;
        d1[0] = dx * ix;
        d1[1] = dy * iy;
        model_calcul_cov(model,mode,1,1.,d1,&covtab);
        COV_RES(ix,iy,iz) = covver * (covtab + covtot - c00) / covtot;
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
static void st_calculate_covtot(Db     *db,
                                int     ix0,
                                int     iy0,
                                int     flag_sym,
                                int     cov_ss[3],
                                int     cov_nn[3],
                                int    *num_tot,
                                double *cov_tot)
{
  int    ix,iy,iz,ix1,iy1,iz1,jx1,jy1,jz1,jx2,jy2,jz2,indg[3];
  int    idx,idy,idz,jdx,jdy,iad,jad;
  double val1,val2,val,ratio;

  /* Initialization */

  for (ix=-cov_nn[0]; ix<=cov_nn[0]; ix++)
    for (iy=-cov_nn[1]; iy<=cov_nn[1]; iy++)
      for (iz=-cov_nn[2]; iz<=cov_nn[2]; iz++)
      {
        COV_TOT(ix,iy,iz) = 0.;
        NUM_TOT(ix,iy,iz) = 0;
      }

  /* Loop on the first point */

  for (iz1=0; iz1<db->getNX(2); iz1++)
    for (iy1=-cov_nn[1]; iy1<=cov_nn[1]; iy1++)
      for (ix1=-cov_nn[0]; ix1<=cov_nn[0]; ix1++)
      {
        jx1 = ix0 + ix1;
        if (jx1 < 0 || jx1 >= db->getNX(0)) continue;
        jy1 = iy0 + iy1;
        if (jy1 < 0 || jy1 >= db->getNX(1)) continue;
        jz1 = iz1;

        indg[0] = jx1;
        indg[1] = jy1;
        indg[2] = jz1;
        iad  = db_index_grid_to_sample(db,indg);
        if (! db->isActive(iad)) continue;
        val1 = db->getVariable(iad,0);
        if (FFFF(val1)) continue;

        /* Loop on the second point within the covariance array */

        for (idz=-cov_nn[2]; idz<=cov_nn[2]; idz++)
          for (idy=-cov_nn[1]; idy<=cov_nn[1]; idy++)
            for (idx=-cov_nn[0]; idx<=cov_nn[0]; idx++)
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
              jad  = db_index_grid_to_sample(db,indg);
              if (! db->isActive(jad)) continue;
              val2 = db->getVariable(jad,0);
              if (FFFF(val2)) continue;

              /* Update the Covariance */

              COV_TOT(idx,idy,idz) += val1 * val2;
              NUM_TOT(idx,idy,idz) += 1;
            }
      }

  /* Scaling */

  ratio = NUM_TOT(0,0,0);
  for (ix=-cov_nn[0]; ix<=cov_nn[0]; ix++)
    for (iy=-cov_nn[1]; iy<=cov_nn[1]; iy++)
      for (iz=-cov_nn[2]; iz<=cov_nn[2]; iz++)
      {
        if (NUM_TOT(ix,iy,iz) <= 0.)
        {
          COV_TOT(ix,iy,iz) = TEST;
        }
        else
        {
          COV_TOT(ix,iy,iz) /= ratio;
        }
      }

  /* Symmetry */

  for (ix=-cov_nn[0]; ix<0; ix++)
    for (iy=-cov_nn[1]; iy<=cov_nn[1]; iy++)
      for (iz=-cov_nn[2]; iz<=cov_nn[2]; iz++)
      {
        val1 = COV_TOT( ix,iy,iz);
        val2 = COV_TOT(-ix,iy,iz);
        val  = (FFFF(val1) || FFFF(val2)) ? TEST : (val1 + val2) / 2.;
        COV_TOT( ix,iy,iz) = COV_TOT(-ix,iy,iz) = val;
      }

  for (ix=-cov_nn[0]; ix<=cov_nn[0]; ix++)
    for (iy=-cov_nn[1]; iy<0; iy++)
      for (iz=-cov_nn[2]; iz<=cov_nn[2]; iz++)
      {
        val1 = COV_TOT(ix,-iy,iz);
        val2 = COV_TOT(ix, iy,iz);
        val  = (FFFF(val1) || FFFF(val2)) ? TEST : (val1 + val2) / 2.;
        COV_TOT(ix, iy,iz) = COV_TOT(ix,-iy,iz) = val;
      }

  if (flag_sym)
    for (ix=-cov_nn[0]; ix<=cov_nn[0]; ix++)
      for (iy=-cov_nn[1]; iy<=cov_nn[1]; iy++)
        for (iz=-cov_nn[2]; iz<0; iz++)
        {
          val1 = COV_TOT(ix,iy,-iz);
          val2 = COV_TOT(ix,iy, iz);
          val  = (FFFF(val1) || FFFF(val2)) ? TEST : (val1 + val2) / 2.;
          COV_TOT(ix,iy,-iz) = COV_TOT(ix,iy, iz) = val;
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
** \param[out] nech          Number of samples in the neighborhood
**
*****************************************************************************/
static void st_neigh_find(Db     *db,
                          int     ix0,
                          int     iy0,
                          int     iz0,
                          int     nei_ss[3],
                          int     nei_nn[3],
                          int    *nech,
                          int    *nei_cur)
{
  int ix,iy,iz,jx,jy,jz,indg[3],number,locrank;

  /* Loop on the pixels of the neighborhood */

  number = 0;
  for (ix=-nei_nn[0]; ix<=nei_nn[0]; ix++)
    for (iy=-nei_nn[1]; iy<=nei_nn[1]; iy++)
      for (iz=-nei_nn[2]; iz<=nei_nn[2]; iz++)
      {
        NEI_CUR(ix,iy,iz) = -1;
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
        NEI_CUR(ix,iy,iz) = rank[number] = locrank;
        flag[number] = 1;
        number++;
      }

  /* Define the returned argument */

  *nech = number;

  return;
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
static int st_neigh_diff(int  nei_ss[3],
                         int  nei_nn[3],
                         int *nei_ref,
                         int *nei_cur)
{
  int ix,iy,iz,flag1,flag2,flag_diff;

  /* Loop on the pixels of the neighborhood */

  flag_diff = 1;
  for (ix=-nei_nn[0]; ix<=nei_nn[0]; ix++)
    for (iy=-nei_nn[1]; iy<=nei_nn[1]; iy++)
      for (iz=-nei_nn[2]; iz<=nei_nn[2]; iz++)
      {
        flag1 = NEI_REF(ix,iy,iz) < 0;
        flag2 = NEI_CUR(ix,iy,iz) < 0;
        if (flag1 != flag2) goto label_end;
      }
  flag_diff = 0;

label_end:

  /* Copy the current neighborhood into the reference neighborhood */

  for (ix=-nei_nn[0]; ix<=nei_nn[0]; ix++)
    for (iy=-nei_nn[1]; iy<=nei_nn[1]; iy++)
      for (iz=-nei_nn[2]; iz<=nei_nn[2]; iz++)
        NEI_REF(ix,iy,iz) = NEI_CUR(ix,iy,iz);

  return(flag_diff);
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
static void st_lhs_exp_3D(int     nech,
                          int     nfeq,
                          int     nei_ss[3],
                          int     nei_nn[3],
                          int     cov_ss[3],
                          int     cov_nn[3],
                          int    *nei_cur,
                          double *cov_tot,
                          double  nugget)
{
  int    ix,iy,iz,jx,jy,jz,i,j,neq;
  double value;

  /* Initializations */

  neq = nech + nfeq;

  /* Covariance part of the L.H.S. */

  i = 0;
  for (ix=-nei_nn[0]; ix<=nei_nn[0]; ix++)
    for (iy=-nei_nn[1]; iy<=nei_nn[1]; iy++)
      for (iz=-nei_nn[2]; iz<=nei_nn[2]; iz++)
      {
        if (NEI_CUR(ix,iy,iz) < 0) continue;

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
  for (i=0; i<nech; i++)
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
static void st_rhs_exp_3D(int     nech,
                          int     nfeq,
                          int     nei_ss[3],
                          int     nei_nn[3],
                          int     cov_ss[3],
                          int     cov_nn[3],
                          int    *nei_cur,
                          double *cov_res)
{
  int ix,iy,iz,neq,i;

  /* Initializations */

  neq = nech + nfeq;

  /* Covariance part of the R.H.S. */

  i = 0;
  for (ix=-nei_nn[0]; ix<=nei_nn[0]; ix++)
    for (iy=-nei_nn[1]; iy<=nei_nn[1]; iy++)
      for (iz=-nei_nn[2]; iz<=nei_nn[2]; iz++)
      {
        if (NEI_CUR(ix,iy,iz) < 0) continue;
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
static double st_estim_exp_3D(Db     *db,
                              int     nei_ss[3],
                              int     nei_nn[3],
                              int    *nei_cur,
                              double *weight)
{
  int i,ix,iy,iz;
  double result;

  /* Initializations */

  result = 0.;
  i = 0;
  for (ix=-nei_nn[0]; ix<=nei_nn[0]; ix++)
    for (iy=-nei_nn[1]; iy<=nei_nn[1]; iy++)
      for (iz=-nei_nn[2]; iz<=nei_nn[2]; iz++)
      {
        if (NEI_CUR(ix,iy,iz) < 0) continue;
        result += weight[i] * db->getVariable(NEI_CUR(ix,iy,iz),0);
        i++;
      }

  return(result);
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
static void st_vario_dump(FILE   *file,
                          int     ix0,
                          int     iy0,
                          int     cov_ss[3],
                          int     cov_nn[3],
                          int    *num_tot,
                          double *cov_tot)
{
  int    ix,iy,iz,num;
  double cov;

  fprintf(file,"*%3d %3d\n",ix0,iy0);

  for (ix=-cov_nn[0]; ix<=cov_nn[0]; ix++)
    for (iy=-cov_nn[1]; iy<=cov_nn[1]; iy++)
      for (iz=-cov_nn[2]; iz<=cov_nn[2]; iz++)
      {
        num = (num_tot == (int *) NULL) ? 0 : NUM_TOT(ix,iy,iz);
        cov = COV_TOT(ix,iy,iz);
        fprintf(file,"%3d %3d %3d %3d %lf\n",ix,iy,iz,num,cov);
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
GEOSLIB_API int anakexp_3D(Db     *db,
                           double *cov_ref,
                           int     cov_radius,
                           int     neigh_ver,
                           int     neigh_hor,
                           int     flag_sym,
                           Model  *model,
                           double  nugget,
                           int     nfeq,
                           int     dbg_ix,
                           int     dbg_iy)
{
  int      i,ix,iy,iz,ndim,nvarin,nech,error,neq,status,ecr;
  int      size_cov,size_nei,flag_new,flag_col;
  int      cov_ss[3],cov_nn[3],nei_ss[3],nei_nn[3],indg[3];
  int     *num_tot,*nei_cur,*nei_ref;
  double  *cov_tot,*cov_res,result;
  FILE    *fildmp;

  /* Initializations */

  error = 1;
  st_global_init(db,db);
  FLAG_EST = 1;
  fildmp   = (FILE *) NULL;
  cov_tot  = cov_res = (double *) NULL;
  num_tot  = nei_cur = nei_ref = (int *) NULL;
  lhs      = rhs = wgt = (double *) NULL;
  ndim     = db->getNDim();
  nvarin   = db->getVariableNumber();
  size_nei = 0;

  /* Prepare the Koption structure */

  if (krige_koption_manage(1,1,KOPTION_PONCTUAL,1,VectorInt())) return(1);

  /* Preliminary checks */

  if (ndim != 3 || ! is_grid(db))
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
    messerr("to the radius of the covariance (%d)",cov_radius);
    goto label_end;
  }

  /* Open the Variogram debugging file */

  if (dbg_ix >= -1 && dbg_ix < db->getNX(0) &&
      dbg_iy >= -1 && dbg_iy < db->getNX(1))
  {
    fildmp = fopen("Vario.dat","w");
    if (fildmp == (FILE *) NULL) goto label_end;
  }

  /* Add the attribute for storing the result */

  IPTR_EST = db->addFields(nvarin,0.);
  if (IPTR_EST < 0) goto label_end;
  DBOUT = db;

  /* Define essential variables */

  nei_nn[0] = MIN(db->getNX(0) - 1,neigh_hor);
  nei_nn[1] = MIN(db->getNX(1) - 1,neigh_hor);
  nei_nn[2] = MIN(db->getNX(2) - 1,neigh_ver);
  size_nei = size_cov = 1;
  for (i=0; i<db->getNDim(); i++)
  {
    nei_ss[i] = 2 * nei_nn[i] + 1;
    cov_nn[i] = 2 * nei_nn[i];
    cov_ss[i] = 2 * cov_nn[i] + 1;
    size_nei *= nei_ss[i];
    size_cov *= cov_ss[i];
  }
  size_nei += nfeq;

  /* Core allocation */

  num_tot = st_icore(size_cov,1);
  if (num_tot == (int    *) NULL) goto label_end;
  nei_cur = st_icore(size_nei,1);
  if (nei_cur == (int    *) NULL) goto label_end;
  nei_ref = st_icore(size_nei,1);
  if (nei_ref == (int    *) NULL) goto label_end;
  cov_tot = st_core(size_cov,1);
  if (cov_tot == (double *) NULL) goto label_end;
  cov_res = st_core(size_cov,1);
  if (cov_res == (double *) NULL) goto label_end;
  st_krige_manage_basic(1,size_nei,size_nei,1,nfeq);
  for (i=0; i<size_nei; i++) nei_ref[i] = -1;

  /* Calculate the discretized covariance of residual variable */

  st_calculate_covres(db,model,cov_ref,cov_radius,flag_sym,
                      cov_ss,cov_nn,cov_res);
  if (dbg_ix == -1 && dbg_iy == -1)
    st_vario_dump(fildmp,-1,-1,cov_ss,cov_nn,(int *) NULL,cov_res);

  /* Loop on the grid nodes */

  status = nech = neq = 0;
  IECH_OUT = ecr = 0;
  for (ix=0; ix<db->getNX(0); ix++)
    for (iy=0; iy<db->getNX(1); iy++)
    {
      flag_col = 1;

      /* Calculate the experimental covariance of total variable */

      st_calculate_covtot(db,ix,iy,flag_sym,cov_ss,cov_nn,num_tot,cov_tot);
      if (dbg_ix == ix && dbg_iy == iy)
        st_vario_dump(fildmp,ix,iy,cov_ss,cov_nn,num_tot,cov_tot);

      for (iz=0; iz<db->getNX(2); iz++,ecr++)
      {
        mes_process("3-D Factorial Kriging Analysis",get_NECH(DBOUT),ecr);
        indg[0] = ix;
        indg[1] = iy;
        indg[2] = iz;
        IECH_OUT = db_index_grid_to_sample(DBOUT,indg);
        debug_index(IECH_OUT+1);

        /* Initialize the result to TEST */

        DBOUT->setArray(IECH_OUT,IPTR_EST,TEST);

        if (FFFF(db->getVariable(IECH_OUT,0)) ||
            ! db->isActive(IECH_OUT)) continue;
        if (debug_query("kriging") ||
            debug_query("nbgh")    ||
            debug_query("results"))
        {
          mestitle(1,"Target location");
          db_sample_print(db,IECH_OUT,1,0,0);
        }

        /* Look for the neighborhood */

        st_neigh_find(db,ix,iy,iz,nei_ss,nei_nn,&nech,nei_cur);
        if (nech <= 0) continue;
        neq = (nfeq == 0) ? nech : nech + 1;

        /* Check if the neighborhood has changed */

        flag_new = flag_col || st_neigh_diff(nei_ss,nei_nn,nei_ref,nei_cur);

        /* If the neighborhood has changed, establish the kriging system */

        flag_col = 0;
        if (flag_new || debug_force())
        {

          /* Establish the L.H.S. of the kriging system */

          st_lhs_exp_3D(nech,nfeq,nei_ss,nei_nn,cov_ss,cov_nn,
                        nei_cur,cov_tot,nugget);
          if (debug_query("kriging")) krige_lhs_print(nech,neq,neq,flag,lhs);

          /* Invert the kriging system */

          if (matrix_invert(lhs,neq,IECH_OUT))
          {
            status = 1;
            continue;
          }

          /* Establish the R.H.S. of the kriging system */

          st_rhs_exp_3D(nech,nfeq,nei_ss,nei_nn,cov_ss,cov_nn,
                        nei_cur,cov_res);
          if (debug_query("kriging"))
            krige_rhs_print(nvarin,nech,neq,neq,flag,rhs);

          /* Derive the kriging weights */

          matrix_product(neq,neq,1,lhs,rhs,wgt);
          if (debug_query("kriging"))
            krige_wgt_print(status,1,1,nfeq,nech,nech,-1,flag,wgt);
        }

        /* Calculate the estimation */

        result = st_estim_exp_3D(db,nei_ss,nei_nn,nei_cur,wgt);
        DBOUT->setArray(IECH_OUT,IPTR_EST,result);
        if (debug_query("results")) st_result_kriging_print(0,nvarin,status);
      }
    }

  /* Set the error return flag */

  error = 0;
  if (fildmp != (FILE *) NULL) fclose(fildmp);

label_end:
  debug_index(0);
  (void) krige_koption_manage(-1,1,KOPTION_PONCTUAL,1,VectorInt());
  st_krige_manage_basic(-1,size_nei,size_nei,1,nfeq);
  num_tot = (int    *) mem_free((char *) num_tot);
  nei_cur = (int    *) mem_free((char *) nei_cur);
  nei_ref = (int    *) mem_free((char *) nei_ref);
  cov_tot = (double *) mem_free((char *) cov_tot);
  cov_res = (double *) mem_free((char *) cov_res);
  return(error);
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
GEOSLIB_API int bayes_simulate(Model  *model,
                               int     nbsimu,
                               double *rmean,
                               double *rcov,
                               double *smean)
{
  double *trimat,*rndmat;
  int     nfeq,il,isimu,nftri,error,rank,memo;

  /* Initializations */

  error  = 1;
  nfeq   = model->getDriftEquationNumber();
  nftri  = nfeq * (nfeq + 1) / 2;
  trimat = rndmat = (double *) NULL;
  memo   = law_get_random_seed();

  /* Core allocation */

  trimat = (double *) mem_alloc(sizeof(double) * nftri,0);
  if (trimat == (double *) NULL) goto label_end;
  rndmat = (double *) mem_alloc(sizeof(double) * nfeq,0);
  if (rndmat == (double *) NULL) goto label_end;

  /* Cholesky decomposition */

  rank = matrix_cholesky_decompose(rcov,trimat,nfeq);
  if (rank > 0)
  {
    messerr("Error in the Cholesky Decomposition of the covariance matrix");
    messerr("Rank of the Matrix = %d",rank);
    messerr("The Drift coefficients have been set to their posterior mean");
    for (isimu=0; isimu<nbsimu; isimu++)
      for (il=0; il<nfeq; il++)
        SMEAN(il,isimu) = rmean[il];
    goto label_suite;
  }

  /* Loop on the simulations */

  for (isimu=0; isimu<nbsimu; isimu++)
  {

    /* Draw a vector of gaussian independent values */

    for (il=0; il<nfeq; il++) rndmat[il] = law_gaussian();

    /* Product of the Lower triangular matrix by the random vector */

    matrix_cholesky_product(1,nfeq,1,trimat,rndmat,&SMEAN(0,isimu));

    /* Add the mean */

    for (il=0; il<nfeq; il++) SMEAN(il,isimu) += rmean[il];
  }

  /* If DEBUG option is switched ON, the values are printed out */

label_suite:
  if (debug_query("bayes"))
  {
    mestitle(1,"Simulation of Drift Coefficients (for Bayesian Simulation)");
    message("Rank     Drift Coefficients\n");
    for (isimu=0; isimu<nbsimu; isimu++)
    {
      message(" %3d ",isimu+1);
      for (il=0; il<nfeq; il++)
        message(" %lf",SMEAN(il,isimu));
      message("\n");
    }
  }

  /* Set the returned error code */

  error = 0;

label_end:
  trimat = (double *) mem_free((char *) trimat);
  rndmat = (double *) mem_free((char *) rndmat);
  law_set_random_seed(memo);
  return(error);
}

/****************************************************************************/
/*!
**  Smooth a regular grid
**
** \return  Error return code
**
** \param[in]  dbgrid    input and output Db grid structure
** \param[in]  neigh     Neigh structure
** \param[in]  type      1 for Uniform; 2 for Gaussian
** \param[in]  range     Range (used for Gaussian only)
**
*****************************************************************************/
GEOSLIB_API int image_smoother(Db    *dbgrid,
                               Neigh *neigh,
                               int    type,
                               double range)
{
  int    i,iech,jech,error,nvarin,nb_neigh,ndim,idelta;
  int   *indn0,*indnl,*indg0,*indgl;
  double data,estim,total,delta,weight,d2,r2;
  Db    *dbaux;

  /* Preliminary checks */

  error = 1;
  dbaux = (Db *) NULL;
  indn0 = indnl = indg0 = indgl = (int *) NULL;
  st_global_init(dbgrid,dbgrid);
  nvarin = dbgrid->getVariableNumber();
  ndim   = dbgrid->getNDim();
  r2     = (type == 1) ? 1. : range * range;
  if (nvarin != 1)
  {
    messerr("The Image Smoother is only programmed for a single variable");
    return(1);
  }
  if (neigh->getType() != NEIGH_IMAGE)
  {
    messerr("This tool requires an IMAGE neighborhood");
    return(1);
  }

  /* Add the attributes for storing the results */

  IPTR_EST = dbgrid->addFields(nvarin,0.);
  if (IPTR_EST < 0) goto label_end;

  /* Core allocation */

  indg0 = db_indg_alloc(dbgrid);
  indgl = db_indg_alloc(dbgrid);

  /* Create the secondary grid for image processing */

  dbaux = st_image_build(neigh,1);
  if (dbaux == (Db *) NULL) goto label_end;
  nb_neigh = get_NECH(dbaux);
  indn0 = db_indg_alloc(dbaux);
  indnl = db_indg_alloc(dbaux);
  db_index_sample_to_grid(dbaux,nb_neigh/2,indn0);

  /* Pre-calculations */

  if (neigh_start(dbaux,neigh)) goto label_end;

  /* Loop on the targets to be processed */

  DBIN   = dbaux;
  DBOUT  = dbgrid;
  for (IECH_OUT=0; IECH_OUT<get_NECH(dbgrid); IECH_OUT++)
  {
    mes_process("Image smoother",get_NECH(DBOUT),IECH_OUT);
    debug_index(IECH_OUT+1);
    if (! DBOUT->isActive(IECH_OUT)) continue;
    db_index_sample_to_grid(DBOUT,IECH_OUT,indg0);

    /* Loop on the neighboring points */

    estim = total = 0.;
    for (iech=0; iech<nb_neigh; iech++)
    {
      if (FFFF(DBIN->getVariable(iech,0))) continue;
      db_index_sample_to_grid(DBIN,iech,indnl);
      d2 = 0.;
      for (i=0; i<ndim; i++)
      {
        idelta   = (indnl[i] - indn0[i]);
        delta    = idelta * DBOUT->getDX(i);
        d2      += delta * delta;
        indgl[i] = indg0[i] +idelta;
        indgl[i] = get_mirror_sample(DBOUT->getNX(i),indgl[i]);
      }
      jech = db_index_grid_to_sample(DBOUT,indgl);
      data = DBOUT->getVariable(jech,0);
      if (! FFFF(data))
      {
        weight = (type == 1) ? 1. : exp(-d2 / r2);
        estim += data * weight;
        total += weight;
      }
    }
    estim = (total <= 0.) ? TEST : estim / total;
    DBOUT->setArray(IECH_OUT,IPTR_EST,estim);
  }

  /* Set the error return flag */

  error = 0;

label_end:
  debug_index(0);
  dbaux = db_delete(dbaux);
  indn0 = db_indg_free(indn0);
  indnl = db_indg_free(indnl);
  indg0 = db_indg_free(indg0);
  indgl = db_indg_free(indgl);
  neigh_stop();
  return(error);
}

/****************************************************************************/
/*!
**  Ponctual Multivariate Kriging under a constraint
**
** \return  Error return code
**
** \param[in]  dbin      input Db structure
** \param[in]  dbout     output Db structure
** \param[in]  model     Model structure (univariate)
** \param[in]  neigh     Neigh structrue
** \param[in]  flag_positive  1 for a positive constraints
**
** \remark  All the variables are estimated using the same model
** \remark  In this procedure, we assume that:
** \remark  - the problem is multivariate ("z" variables)
** \remark  - the constraints is stored in "f" (only used in dbout)
**
*****************************************************************************/
GEOSLIB_API int krigsum_f(Db    *dbin,
                        Db    *dbout,
                        Model *model,
                        Neigh *neigh,
                        int    flag_positive)
{
  double  *lterm,seisloc,seistot,estim;
  int     *icols,*active,error,iptr_mem,correct;
  int      nvarmod,nvarin,status,flag_new_nbgh,nech,ivar,nred,neq;

  /* Preliminary checks */

  error      = 1;
  st_global_init(dbin,dbout);
  icols      = active = (int     *) NULL;
  lterm      = (double  *) NULL;
  nvarin     = dbin->getVariableNumber();
  nvarmod    = model->getVariableNumber();
  FLAG_EST   = 1;
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
  if (neigh->getType() != NEIGH_UNIQUE)
  {
    messerr("This procedure is currently limited to the Unique Neighborhood");
    goto label_end;
  }

  /* Core allocation */

  icols  = (int    *) mem_alloc(sizeof(int)    * nvarin,0);
  if (icols  == (int    *) NULL) goto label_end;
  active = (int    *) mem_alloc(sizeof(int)    * nvarin,0);
  if (active == (int    *) NULL) goto label_end;
  lterm  = (double *) mem_alloc(sizeof(double) * nvarin,0);
  if (lterm  == (double *) NULL) goto label_end;

  /* Save the columns for variable definitions */

  for (ivar=0; ivar<nvarin; ivar++)
    icols[ivar] = db_attribute_identify(dbin,LOC_Z,ivar);
  dbin->clearLocators(LOC_Z);
  dbin->setLocatorByAttribute(icols[0],LOC_Z);
  if (st_check_environment(1,1,model,neigh)) goto label_end;

  /* Add the attributes for storing the results */

  iptr_mem = dbout->addFields(nvarin,0);
  if (iptr_mem < 0) goto label_end;

  /* Pre-calculations */

  if (neigh_start(DBIN,neigh)) goto label_end;
  if (st_model_manage(1,model)) goto label_end;
  if (st_krige_manage(1,nvarin,model,neigh)) goto label_end;
  if (krige_koption_manage(1,1,KOPTION_PONCTUAL,1,VectorInt())) goto label_end;

  /* Loop on the variables */

  status = 0;
  for (ivar=0; ivar<nvarin; ivar++)
  {
    dbin->clearLocators(LOC_Z);
    dbin->setLocatorByAttribute(icols[ivar],LOC_Z);
    IPTR_EST  = iptr_mem + ivar;
    IECH_NBGH = -1;
    (void) sprintf(string,"Kriging of variable #%d at sample",ivar+1);

    /* Loop on the targets to be processed */

    for (IECH_OUT=0; IECH_OUT<get_NECH(DBOUT); IECH_OUT++)
    {
      mes_process(string,get_NECH(DBOUT),IECH_OUT);
      debug_index(IECH_OUT+1);
      if (! dbout->isActive(IECH_OUT)) continue;
      if (debug_query("kriging") ||
          debug_query("nbgh")    ||
          debug_query("results"))
      {
        mestitle(1,"Target location");
        db_sample_print(dbout,IECH_OUT,1,0,0);
      }
      
      /* Select the Neighborhood */
      
      flag_new_nbgh = st_neigh(neigh,&status,&nech);
      if (status) goto label_store;

      /* Establish the kriging L.H.S. */
      
      if (flag_new_nbgh || flag_continuous_kriging(neigh) || debug_force())
      {
        st_prepar(model,neigh,nech,&status,&nred,&neq);
        if (status) goto label_store;
        st_data_dual(model,NULL,nech,nred,&lterm[ivar]);
      }
      
      /* Establish the kriging R.H.S. */
      
      st_rhs(model,nech,neq,nvarmod,NULL,&status);
      if (status) goto label_store;
      st_rhs_iso2hetero(neq,nvarmod);
      if (debug_query("kriging"))
        krige_rhs_print(nvarmod,nech,neq,nred,flag,rhs);
      
      /* Perform the estimation */
      
    label_store:
      st_estimate(model,NULL,status,0,nech,nvarmod,nred);
      if (debug_query("results"))
        st_result_kriging_print(neigh->getFlagXvalid(),nvarmod,status);
    }
  }

  /* Posterior scaling */
  
  for (IECH_OUT=0; IECH_OUT<get_NECH(DBOUT); IECH_OUT++)
  {
    correct = 0;
    for (ivar=0; ivar<nvarin; ivar++) active[ivar] = 0;
    
    /* Implicit loop until the solution is acceptable */
    
    while (! correct)
    {
      seistot = 0.;
      seisloc = dbout->getExternalDrift(IECH_OUT,0);
      for (ivar=0; ivar<nvarin; ivar++) 
      {
        if (active[ivar]) continue;
        seistot += lterm[ivar];
        seisloc -= get_ARRAY(dbout,IECH_OUT,iptr_mem + ivar);
      }
      if (seistot == 0.) 
      {
        messerr("The sum of scaling terms is zero. No correction is possible");
        goto label_end;
      }
      
      for (ivar=0; ivar<nvarin; ivar++)
      {
        if (active[ivar])
          estim = 0;
        else
          estim = (get_ARRAY(dbout,IECH_OUT,iptr_mem+ivar) + 
                   lterm[ivar] * seisloc / seistot);
        dbout->setArray(IECH_OUT,iptr_mem+ivar,estim);
      }
      
      correct = 1;
      for (ivar=0; ivar<nvarin; ivar++)
      {
        active[ivar] = (get_ARRAY(dbout,IECH_OUT,iptr_mem+ivar) < 0);
        if (active[ivar]) correct = 0;
      }
      if (! flag_positive) correct = 1;
    }
  }
  
  /* Set the error return flag */

  error = 0;

label_end:
  debug_index(0);
  (void) st_model_manage(-1,model);
  (void) st_krige_manage(-1,nvarin,model,neigh);
  (void) krige_koption_manage(-1,1,KOPTION_PONCTUAL,1,VectorInt());
  neigh_stop();
  icols  = (int    *) mem_free((char *) icols);
  active = (int    *) mem_free((char *) active);
  lterm  = (double *) mem_free((char *) lterm);
  return(error);
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
static int st_check_positivity(Db      *db3grid,
                               int      ix,
                               int      iy,
                               int      nvarin,
                               int      iptr_prop,
                               double  *proptab,
                               double  *lterm)
{
  int indg[3],n_wrong,iech,iz,ivar;
  double prop;

  indg[0] = ix;
  indg[1] = iy;
  n_wrong = 0;
  for (iz=0; iz<db3grid->getNX(2); iz++)
  {
    indg[2] = iz;
    iech = db_index_grid_to_sample(db3grid,indg);
    if (! db3grid->isActive(iech)) continue;

    /* No estimated proportion at this level */

    for (ivar=0; ivar<nvarin; ivar++)
    {
      prop = PROPTAB(ivar,iz);
      if (prop >= 0) continue;
      db3grid->setArray(iech,iptr_prop+ivar,0.);
      LTERM(ivar,iz) = 0.;
      n_wrong++;
    }
  }
  return(n_wrong);
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
                                       double  seisval,
                                       double *proptab)
{
  double prop;
  int    iz,error;

  error = 0;
  if (FFFF(seisval)) return(0);

  /* Calculate the average proportion corresponding to the seismic */

  prop = 0.;
  for (iz=0; iz<nz; iz++) prop += PROPTAB(fsum,iz);
  prop /= nz;

  /* Compare seismic and resulting average propotion */

  if (ABS(prop - seisval) > EPS) 
  {
    messerr("Block (%d,%d,%d) - Mismatch between proportion (%lf) and seismic (%lf)",
            ix+1,iy+1,iz+1,prop,seisval);
    error++;
  }
  return(error);
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
** \param[in]  dbin      input Db structure
** \param[in]  db3grid   output 3-D Grid Db structure
** \param[in]  db2grid   output 2-D Grid Db structure
** \param[in]  fsum      Rank of the proportion facies which average up
**                       to the seismic (no constraint if negative)
** \param[in]  model     Model structure (monovariate)
** \param[in]  neigh     Neigh structure (Unique or Bench)
**
*****************************************************************************/
GEOSLIB_API int krigmvp_f(Db    *dbin,
                        Db    *db3grid,
                        Db    *db2grid,
                        int    fsum,
                        Model *model,
                        Neigh *neigh)
{
  int     *icols,indg[3],nloc,error;
  int      status,nech,neq,nred,nvarin,nvarmod,flag_new_nbgh,nfeq,nz,pivot;
  int      ix,iy,iz,ivar,i,iech,jech,iptr_prop,n_wrong,iter,nsize,flag_correc;
  double  *lterm,*lback,*proptab,*cc,*xx,*bb;
  double   seisval,correc,proploc,lsum;

  /* Preliminary checks */

  error      =  1;
  st_global_init(dbin,db3grid);
  icols      = (int     *) NULL;
  lterm      = lback = proptab = cc = xx = bb = (double  *) NULL;
  nvarin     = dbin->getVariableNumber();
  nvarmod    = model->getVariableNumber();
  nfeq       = model->getDriftEquationNumber();
  iptr_prop  = iter = nsize = 0;
  seisval    = TEST;
  FLAG_EST   = 1;
  FLAG_WGT   = 1;
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
  if (neigh->getType() != NEIGH_UNIQUE && neigh->getType() != NEIGH_BENCH)
  {
    messerr("This procedure is currently limited to the Unique or Bench Neighborhood");
    goto label_end;
  }
  if (! is_grid(db3grid))
  {
    messerr("This procedure needs a Grid 3-D Db structure");
    goto label_end;
  }
  if (fsum >= 0)
  { 
    if (! is_grid(db2grid))
    {
      messerr("This procedure needs a Grid 2-D Db structure");
      goto label_end;
    }
    if (! db_grid_match(db2grid,db3grid)) goto label_end;
    if (db2grid->getExternalDriftNumber() != 1)
    {
      messerr("This procedure requires one External Drift in the 2-D Grid Db");
      messerr("The number of External Drift is currently equal to %d",
              db2grid->getExternalDriftNumber());
      goto label_end;
    }
    if (fsum >= nvarin)
    {
      messerr("Error in the rank of the variable whose average vertical proportion");
      messerr("Should match the seismic of the 2-D Db grid (%d)",fsum+1);
      messerr("It should lie between 1 and %d",nvarin);
      goto label_end;
    }
  }
  nz = db3grid->getNX(2);

  /* Core allocation */

  lback = (double *) mem_alloc(sizeof(double) * nvarin * nz,0);
  if (lback == (double *) NULL) goto label_end;
  lterm = (double *) mem_alloc(sizeof(double) * nvarin * nz,0);
  if (lterm == (double *) NULL) goto label_end;
  for (i=0; i<nvarin*nz; i++) lterm[i] = lback[i] = TEST;
  icols = (int    *) mem_alloc(sizeof(int)    * nvarin,0);
  if (icols == (int    *) NULL) goto label_end;
  if (fsum >= 0)
  {
    proptab = (double *) mem_alloc(sizeof(double) * nz * nvarin,0);
    if (proptab == (double *) NULL) goto label_end;
  }

  /* Save the columns for variable definitions */

  for (ivar=0; ivar<nvarin; ivar++)
    icols[ivar] = db_attribute_identify(dbin,LOC_Z,ivar);
  dbin->clearLocators(LOC_Z);
  dbin->setLocatorByAttribute(icols[0],LOC_Z);
  if (st_check_environment(1,1,model,neigh)) goto label_end;

  /* Add the attributes for storing the results */

  iptr_prop = db3grid->addFields(nvarin,0);
  if (iptr_prop < 0) goto label_end;

  /* Pre-calculations */

  if (neigh_start(DBIN,neigh)) goto label_end;
  if (st_model_manage(1,model)) goto label_end;
  if (st_krige_manage(1,nvarmod,model,neigh)) goto label_end;
  if (krige_koption_manage(1,1,KOPTION_PONCTUAL,1,VectorInt())) goto label_end;

  /* Loop on the target grid nodes to be processed */

  status = 0;
  for (ivar=0; ivar<nvarin; ivar++)
  {
    dbin->clearLocators(LOC_Z);
    dbin->setLocatorByAttribute(icols[ivar],LOC_Z);
    IPTR_EST  = iptr_prop + ivar;
    IECH_NBGH = -1;
    (void) sprintf(string,"Kriging of proportion #%d at sample",ivar+1);

    /* Loop on the target grid nodes */

    for (iz=IECH_OUT=0; iz<db3grid->getNX(2); iz++)
    {
      for (iy=0; iy<db3grid->getNX(1); iy++)
        for (ix=0; ix<db3grid->getNX(0); ix++,IECH_OUT++)
        {
          mes_process(string,get_NECH(DBOUT),IECH_OUT);
          debug_index(IECH_OUT+1);
          if (! db3grid->isActive(IECH_OUT)) continue;
          if (debug_query("kriging") ||
              debug_query("nbgh")    ||
              debug_query("results"))
          {
            mestitle(1,"Target location");
            db_sample_print(db3grid,IECH_OUT,1,0,0);
          }
          
          /* Select the Neighborhood */
          
          flag_new_nbgh = st_neigh(neigh,&status,&nech);
          if (status) goto label_store;
          
          /* Establish the kriging L.H.S. */
          
          if (flag_new_nbgh || flag_continuous_kriging(neigh) || debug_force())
          {
            st_prepar(model,neigh,nech,&status,&nred,&neq);
            if (status) goto label_store;
            st_data_dual(model,NULL,nech,nred,&LBACK(ivar,iz));
          }
          
          /* Establish the kriging R.H.S. */

          st_rhs(model,nech,neq,nvarmod,NULL,&status);
          if (status) goto label_store;
          st_rhs_iso2hetero(neq,nvarmod);
          if (debug_query("kriging"))
            krige_rhs_print(nvarmod,nech,neq,nred,flag,rhs);
          
          /* Derive the kriging weights */
          
          if (FLAG_WGT)
          {
            matrix_product(nred,nred,nvarmod,lhs,rhs,wgt);
            if (debug_query("kriging"))
              krige_wgt_print(status,nvarmod,nvarmod,nfeq,nech,nred,-1,
                              flag,wgt);
          }

          /* Perform the estimation */
          
        label_store:
          st_estimate(model,NULL,status,0,nech,nvarmod,nred);
          if (debug_query("results"))
            st_result_kriging_print(neigh->getFlagXvalid(),nvarmod,status);
        }
    }
  }

  /* Discard impossible constraints */

  lsum = 0.;
  for (ivar=0; ivar<nvarin; ivar++)
    for (iz=0; iz<nz; iz++)
      lsum += LBACK(ivar,iz);
  if (lsum <= EPS) 
  {
    error = 0;
    goto label_end;
  }

  /* Core allocation */

  nsize = (fsum >= 0) ? nz + 1 : nz;
  cc = (double *) mem_alloc(sizeof(double) * nsize * (nsize+1)/2,0);
  if (cc == (double *) NULL) goto label_end;
  xx = (double *) mem_alloc(sizeof(double) * nsize,0);
  if (xx == (double *) NULL) goto label_end;
  bb = (double *) mem_alloc(sizeof(double) * nsize,0);
  if (bb == (double *) NULL) goto label_end;

  /* Posterior corrections */

  for (iy=0; iy<db3grid->getNX(1); iy++)
    for (ix=0; ix<db3grid->getNX(0); ix++)
    {
      iter = 0;
      for (i=0; i<nvarin*nz; i++) lterm[i] = lback[i];

      /* Loop on the iterations for the same block */

    label_suite:
      iter++;
  
      indg[0] = ix;
      indg[1] = iy;
      nloc = nz;
      if (fsum >= 0)
      {
        indg[2] = 0;
        jech    = db_index_grid_to_sample(db2grid,indg);
        seisval = db2grid->getExternalDrift(jech,0);
        if (! FFFF(seisval)) nloc++;
      }

      /* Establish the system */
      
      for (i=0; i<nloc * nloc; i++) cc[i] = 0.;
      for (iz=0; iz<db3grid->getNX(2); iz++)
      {
        bb[iz]  = -1;
        indg[2] = iz;
        iech = db_index_grid_to_sample(db3grid,indg);
        for (ivar=0; ivar<nvarin; ivar++)
        {
          bb[iz] += get_ARRAY(db3grid,iech,iptr_prop+ivar);
          CC(iz,iz) += LTERM(ivar,iz);
        }
      }
      
      if (nloc > nz)
      {
        bb[nz] = -2 * seisval;
        for (iz=0; iz<db3grid->getNX(2); iz++)
        {
          indg[2] = iz;
          iech = db_index_grid_to_sample(db3grid,indg);
          bb[nz]    += get_ARRAY(db3grid,iech,iptr_prop+fsum);
          CC(iz,nz) += LTERM(fsum,iz);
          CC(nz,nz) += LTERM(fsum,iz);
        }
      }

      /* Check if correction is needed */

      flag_correc = 0;
      for (i=0; i<nloc && flag_correc==0; i++)
        if (ABS(bb[i]) > EPS) flag_correc = 1;
      if (! flag_correc) continue;
        
      /* Solve the system */

      if (matrix_solve(0,cc,bb,xx,nloc,1,&pivot))
        messageAbort("Core problem in matrix_solve");
      if (pivot > 0)
      {
        messerr("Error during the inversion of the constraint matrix");
        messerr("The constraints may be redundant");
        goto label_end;
      }

      /* Perform the correction */

      for (iz=0; iz<db3grid->getNX(2); iz++)
      {
        indg[2] = iz;
        iech = db_index_grid_to_sample(db3grid,indg);
        for (ivar=0; ivar<nvarin; ivar++)
        {
          proploc = get_ARRAY(db3grid,iech,iptr_prop+ivar);
          correc  = xx[iz];
          if (nloc > nz && ivar == fsum) correc += xx[nz];
          PROPTAB(ivar,iz) = proploc - correc * LTERM(ivar,iz);
        }
      }

      /* Check the seismic constraints */
      
      if (fsum >= 0)
        (void) st_check_constraint_seismic(ix,iy,nz,nvarin,
                                           fsum,seisval,proptab);

      /* Check positivity of the proportions */

      n_wrong = st_check_positivity(db3grid,ix,iy,nvarin,iptr_prop,
                                    proptab,lterm);
      if (n_wrong > 0) goto label_suite;

      /* Final update the proportions */

      for (iz=0; iz<nz; iz++)
        for (ivar=0; ivar<nvarin; ivar++)
        {
          indg[2] = iz;
          iech = db_index_grid_to_sample(db3grid,indg);
          db3grid->setArray(iech,iptr_prop+ivar,PROPTAB(ivar,iz));
        }
    }

  /* Set the error return flag */

  error = 0;

label_end:
  debug_index(0);
  (void) st_model_manage(-1,model);
  (void) st_krige_manage(-1,nvarmod,model,neigh);
  (void) krige_koption_manage(-1,1,KOPTION_PONCTUAL,1,VectorInt());
  neigh_stop();
  icols   = (int    *) mem_free((char *) icols);
  lback   = (double *) mem_free((char *) lback);
  lterm   = (double *) mem_free((char *) lterm);
  cc      = (double *) mem_free((char *) cc);
  xx      = (double *) mem_free((char *) xx);
  bb      = (double *) mem_free((char *) bb);
  proptab = (double *) mem_free((char *) proptab);
  return(error);
}

/****************************************************************************/
/*!
**  Perform kriging and return the dimensions
**
** \return  Error return code
**
** \param[in]  dbin      input Db structure
** \param[in]  dbout     output Db structure
** \param[in]  model     Model structure
** \param[in]  neigh     Neigh structrue
** \param[in]  iech0     Rank of the target sample
** \param[in]  calcul    Kriging calculation option (::ENUM_KOPTIONS)
** \param[in]  ndisc     Array giving the discretization counts
**
** \param[out] ndim_ret  Output space dimension
** \param[out] nech_ret  Output number of active samples
** \param[out] nred_ret  Output number of equations
** \param[out] nrhs_ret  Output number of RHS
**
*****************************************************************************/
GEOSLIB_API int krigtest_dimension(Db    *dbin,
                                   Db    *dbout,
                                   Model *model,
                                   Neigh *neigh,
                                   int    iech0,
                                   int    calcul,
                                   VectorInt ndisc,
                                   int   *ndim_ret,
                                   int   *nech_ret,
                                   int   *nred_ret,
                                   int   *nrhs_ret)
{
  int iext,error,status,nech,neq,nred,nvar,inostat;

  /* Preliminary checks */

  error =  1;
  iext  = -1;
  nvar  =  0;
  st_global_init(dbin,dbout);
  FLAG_EST  = 1;
  FLAG_STD  = 1;
  FLAG_WGT  = 1;
  if (st_check_environment(1,1,model,neigh)) goto label_end;
  if (manage_external_info(1,LOC_F,DBIN,DBOUT,&iext)) goto label_end;
  if (manage_external_info(1,LOC_NOSTAT,DBIN,DBOUT,
                           &inostat)) goto label_end;
  nvar  = model->getVariableNumber();

  /* Pre-calculations */

  if (neigh_start(DBIN,neigh)) goto label_end;
  if (st_model_manage(1,model)) goto label_end;
  if (st_krige_manage(1,nvar,model,neigh)) goto label_end;
  if (krige_koption_manage(1,1,calcul,1,ndisc)) goto label_end;

  /* Initialize the target */

  status = 0;
  IECH_OUT = iech0;
  if (iech0 < 0 || iech0 >= get_NECH(dbout)) goto label_end;

  /* Select the Neighborhood */

  (void) st_neigh(neigh,&status,&nech);
  if (status)
  {
    messerr("No valid neighborhood can be found for this target");
    goto label_end;
  }

  /* Establish the kriging L.H.S. */

  st_prepar(model,neigh,nech,&status,&nred,&neq);
  if (status) goto label_end;
  
  /* Set the error return flag */

  error     = 0;
  *ndim_ret = dbin->getNDim();
  *nech_ret = nech;
  *nred_ret = nred;
  *nrhs_ret = nvar;

label_end:
  (void) st_model_manage(-1,model);
  (void) st_krige_manage(-1,nvar,model,neigh);
  (void) manage_external_info(-1,LOC_F,DBIN,DBOUT,&iext);
  (void) manage_external_info(-1,LOC_NOSTAT,DBIN,DBOUT,&inostat);
  neigh_stop();
  return(error);
}

/****************************************************************************/
/*!
**  Perform kriging and return working arrays
**
** \return  Error return code
**
** \param[in]  dbin      input Db structure
** \param[in]  dbout     output Db structure
** \param[in]  model     Model structure
** \param[in]  neigh     Neigh structrue
** \param[in]  iech0     Rank of the target sample
** \param[in]  calcul    Kriging calculation option (::ENUM_KOPTIONS)
** \param[in]  ndisc     Array giving the discretization counts
**
** \param[out] nred_out  Output number of equations
** \param[out] nrhs_out  Output number of RHS
** \param[out] xyz_out   Output array containing the coordinates
** \param[out] data_out  Output array of (extended) data (Dimension: nred_out)
** \param[out] lhs_out   Output LHS matrix (Dimension: nred_out * nred_out)
** \param[out] rhs_out   Output RHS matrix (Dimension: nred_out * nrhs_out)
** \param[out] wgt_out   Output weight matrix (Dimension: nred_out * nrhs_out)
** \param[out] zam_out   Output ZAM matrix (Dimension: nred_out)
** \param[out] var_out   Output variance matrix (Dimension: nrhs_out * nrhs_out)
**
*****************************************************************************/
GEOSLIB_API int krigtest_f(Db     *dbin,
                         Db     *dbout,
                         Model  *model,
                         Neigh  *neigh,
                         int     iech0,
                         int     calcul,
                         VectorInt ndisc,
                         int     nred_out,
                         int     nrhs_out,
                         double *xyz_out,
                         double *data_out,
                         double *lhs_out,
                         double *rhs_out,
                         double *wgt_out,
                         double *zam_out,
                         double *var_out)
{
  int iext,ivar,jvar,ecr,inostat,status,nech,neq,nred,nvar,nfeq,iech,idim,ndim;
  int error;
  double ldum;

  /* Preliminary checks */

  error =  1;
  iext  = -1;
  nvar  =  nred = 0;
  st_global_init(dbin,dbout);
  FLAG_EST  = 1;
  FLAG_STD  = 1;
  FLAG_WGT  = 1;
  if (st_check_environment(1,1,model,neigh)) goto label_end;
  if (manage_external_info(1,LOC_F,DBIN,DBOUT,&iext)) goto label_end;
  if (manage_external_info(1,LOC_NOSTAT,DBIN,DBOUT,
                           &inostat)) goto label_end;
  nvar  = model->getVariableNumber();
  nfeq  = model->getDriftEquationNumber();
  ndim  = dbin->getNDim();

  /* Add the attributes for storing the results */

  if (FLAG_EST)
  {
    IPTR_EST = dbout->addFields(nvar,0.);
    if (IPTR_EST < 0) goto label_end;
  }
  if (FLAG_STD)
  {
    IPTR_STD = dbout->addFields(nvar,0.);
    if (IPTR_STD < 0) goto label_end;
  }

  /* Pre-calculations */

  if (neigh_start(DBIN,neigh)) goto label_end;
  if (st_model_manage(1,model)) goto label_end;
  if (st_krige_manage(1,nvar,model,neigh)) goto label_end;
  if (krige_koption_manage(1,1,calcul,1,ndisc)) goto label_end;
  if (FLAG_STD) st_variance0(model,nvar,VectorDouble());

  /* Initialize the target */

  status = 0;
  IECH_OUT = iech0;
  if (iech0 < 0 || iech0 >= get_NECH(dbout)) goto label_end;

  /* Select the Neighborhood */

  (void) st_neigh(neigh,&status,&nech);
  if (status) goto label_store;

  /* Establish the kriging L.H.S. */

  st_prepar(model,neigh,nech,&status,&nred,&neq);
  if (status) goto label_store;
  st_data_dual(model,NULL,nech,nred,&ldum);

  /* Establish the kriging R.H.S. */

  st_rhs(model,nech,neq,nvar,NULL,&status);
  if (status) goto label_store;
  st_rhs_iso2hetero(neq,nvar);

  /* Derive the kriging weights */

  if (FLAG_WGT)
  {
    matrix_product(nred,nred,nvar,lhs,rhs,wgt);
    if (debug_query("kriging"))
      krige_wgt_print(status,nvar,nvar,nfeq,nech,nred,-1,flag,wgt);
  }

  /* Perform the estimation */

label_store:
  st_estimate(model,NULL,status,neigh->getFlagXvalid(),nech,nvar,nred);
  if (debug_query("results"))
    st_result_kriging_print(neigh->getFlagXvalid(),nvar,status);

  /* Store the output arrays */

  if (nred != nred_out || nrhs_out != nvar)
  {
    messerr("The dimension of the output arrays is incorrect:");
    messerr("- Number of equations: forecast (%d) - actual (%d)",nred_out,neq);
    messerr("- Number of variables: forecast (%d) - actual (%d)",nrhs_out,nvar);
    goto label_end;
  }

  /* Store the coordinates */

  for (idim=ecr=0; idim<ndim; idim++)
    for (iech=0; iech<nech; iech++,ecr++)
      xyz_out[ecr] = st_get_idim(rank[iech],idim);
  (void) memcpy(data_out,zext ,sizeof(double) * nred);
  (void) memcpy(zam_out ,zam1 ,sizeof(double) * nred);
  (void) memcpy(lhs_out ,lhs_b,sizeof(double) * nred * nred);
  (void) memcpy(rhs_out ,rhs  ,sizeof(double) * nred * nvar);
  (void) memcpy(wgt_out ,wgt  ,sizeof(double) * nred * nvar);
  for (ivar=ecr=0; ivar<nvar; ivar++)
    for (jvar=0; jvar<nvar; jvar++,ecr++)
      var_out[ecr] = st_variance(model,ivar,jvar,nred);

  /* Set the error return flag */

  error = 0;

label_end:
  (void) st_model_manage(-1,model);
  (void) st_krige_manage(-1,nvar,model,neigh);
  (void) krige_koption_manage(-1,1,calcul,1,ndisc);
  (void) manage_external_info(-1,LOC_F,DBIN,DBOUT,&iext);
  (void) manage_external_info(-1,LOC_NOSTAT,DBIN,DBOUT,&inostat);
  neigh_stop();
  return(error);
}

/****************************************************************************/
/*!
**  Transform the Kriging results from gaussian to raw
**
** \param[in]  anam      Anam structure
**
** \remark  This procedure is designed for the monovariate case
** \remark  It assumes that the kriging estimate and variance are already
** \remark  calculated
**
*****************************************************************************/
static void st_transform_gaussian_to_raw(Anam *anam)
{
  if (anam == (Anam *) NULL) return;
  AnamHermite* anam_hermite = dynamic_cast<AnamHermite*>(anam);

  /* Get the estimation */

  double est = get_ARRAY(DBOUT,IECH_OUT,IPTR_EST);

  /* Get the variance of the kriging error */

  double std = sqrt(get_ARRAY(DBOUT,IECH_OUT,IPTR_STD));

  /* Calculate the conditional expectation */

  double condexp = hermiteCondExpElement(est,std,anam_hermite->getPsiHn());
  DBOUT->setArray(IECH_OUT,IPTR_EST,condexp);

  /* Calculate the conditional variance */

  double condvar = hermiteEvaluateZ2(est,std,anam_hermite->getPsiHn());
  condvar -= condexp * condexp;
  DBOUT->setArray(IECH_OUT,IPTR_STD,condvar);
}

/****************************************************************************/
/*!
**  Ponctual Kriging in the Anamorphosed Gaussian Model
**
** \return  Error return code
**
** \param[in]  dbin      input Db structure
** \param[in]  dbout     output Db structure
** \param[in]  anam      Anam structure
** \param[in]  model     Model structure
** \param[in]  neigh     Neigh structrue
**
*****************************************************************************/
GEOSLIB_API int kriggam_f(Db    *dbin,
                        Db    *dbout,
                        Anam  *anam,
                        Model *model,
                        Neigh *neigh)
{
  int error,status,nech,neq,nred,nvar,flag_new_nbgh,nfeq;
  double ldum;

  /* Preliminary checks */

  error = 1;
  nvar  = 0;
  st_global_init(dbin,dbout);
  FLAG_EST  = 1;
  FLAG_STD  = 1;
  if (st_check_environment(1,1,model,neigh)) goto label_end;
  nvar  = model->getVariableNumber();
  nfeq  = model->getDriftEquationNumber();

  /* Preliminary check */

  if (nvar != 1)
  {
    messerr("This procedure is limited to the monovariate case");
    goto label_end;
  }

  /* Add the attributes for storing the results */

  if (FLAG_EST)
  {
    IPTR_EST = dbout->addFields(nvar,0.);
    if (IPTR_EST < 0) goto label_end;
  }
  if (FLAG_STD)
  {
    IPTR_STD = dbout->addFields(nvar,0.);
    if (IPTR_STD < 0) goto label_end;
  }

  /* Pre-calculations */

  if (neigh_start(DBIN,neigh)) goto label_end;
  if (st_model_manage(1,model)) goto label_end;
  if (st_krige_manage(1,nvar,model,neigh)) goto label_end;
  if (krige_koption_manage(1,1,KOPTION_PONCTUAL,1,VectorInt())) goto label_end;
  if (FLAG_STD) st_variance0(model,nvar,VectorDouble());

  /* Loop on the targets to be processed */

  status = 0;
  for (IECH_OUT=0; IECH_OUT<get_NECH(DBOUT); IECH_OUT++)
  {
    mes_process("Kriging sample",get_NECH(DBOUT),IECH_OUT);
    debug_index(IECH_OUT+1);
    if (! dbout->isActive(IECH_OUT)) continue;
    if (debug_query("kriging") ||
        debug_query("nbgh")    ||
        debug_query("results"))
    {
      mestitle(1,"Target location");
      db_sample_print(dbout,IECH_OUT,1,0,0);
    }

    /* Select the Neighborhood */

    flag_new_nbgh = st_neigh(neigh,&status,&nech);
    if (status) goto label_store;

    /* Establish the kriging L.H.S. */

    if (flag_new_nbgh || flag_continuous_kriging(neigh) || debug_force())
    {
      st_prepar(model,neigh,nech,&status,&nred,&neq);
      if (status) goto label_store;
      st_data_dual(model,NULL,nech,nred,&ldum);
    }

    /* Establish the kriging R.H.S. */

    st_rhs(model,nech,neq,nvar,NULL,&status);
    if (status) goto label_store;
    st_rhs_iso2hetero(neq,nvar);
    if (debug_query("kriging"))
      krige_rhs_print(nvar,nech,neq,nred,flag,rhs);

    /* Derive the kriging weights */

    if (FLAG_WGT)
    {
      matrix_product(nred,nred,nvar,lhs,rhs,wgt);
      if (debug_query("kriging"))
        krige_wgt_print(status,nvar,nvar,nfeq,nech,nred,-1,flag,wgt);
    }

    /* Perform the estimation */

  label_store:
    st_estimate(model,NULL,status,neigh->getFlagXvalid(),nech,nvar,nred);

    /* Transform the gaussian estimates into raw estimates */

    st_transform_gaussian_to_raw(anam);

    if (debug_query("results"))
      st_result_kriging_print(neigh->getFlagXvalid(),nvar,status);
  }

  /* Set the error return flag */

  error = 0;

label_end:
  debug_index(0);
  (void) st_model_manage(-1,model);
  (void) st_krige_manage(-1,nvar,model,neigh);
  (void) krige_koption_manage(-1,1,KOPTION_PONCTUAL,1,VectorInt());
  neigh_stop();
  return(error);
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
** \param[in]  neigh       Neigh structrue
** \param[in]  ndisc       Array giving the discretization counts
** \param[in]  flag_est    Option for the storing the estimation
** \param[in]  flag_std    Option for the storing the standard deviation
** \param[in]  rank_colcok Option for running Collocated Cokriging
**
*****************************************************************************/
GEOSLIB_API int krigcell_f(Db    *dbin,
                         Db    *dbout,
                         Model *model,
                         Neigh *neigh,
                         VectorInt ndisc,
                         int    flag_est,
                         int    flag_std,
                         VectorInt rank_colcok)
{
  int iext,error,status,nech,neq,nred,nvar,flag_new_nbgh,nfeq,ndim,inostat;
  double ldum;

  /* Preliminary checks */

  error =  1;
  iext  = -1;
  nvar  =  0;
  st_global_init(dbin,dbout);
  FLAG_EST  = flag_est;
  FLAG_STD  = flag_std;
  FLAG_WGT  = flag_std;
  if (st_check_colcok(dbin,dbout,rank_colcok.data())) goto label_end;
  if (st_check_environment(1,1,model,neigh)) goto label_end;
  if (manage_external_info(1,LOC_F,DBIN,DBOUT,&iext)) goto label_end;
  if (manage_external_info(1,LOC_NOSTAT,DBIN,DBOUT,
                           &inostat)) goto label_end;
  nvar  = model->getVariableNumber();
  nfeq  = model->getDriftEquationNumber();
  ndim  = model->getDimensionNumber();

  /* Check that the variables containing the cell dimensions exist */
  if (dbout->getBlockExtensionNumber() != ndim)
  {
    messerr("The number of Extension Variables in the Output Db is %d",
            dbout->getBlockExtensionNumber());
    messerr("The Space Dimension in the Model is %d",ndim);
    messerr("These two values should be equal");
    goto label_end;
  }

  /* Add the attributes for storing the results */

  if (FLAG_EST)
  {
    IPTR_EST = dbout->addFields(nvar,0.);
    if (IPTR_EST < 0) goto label_end;
  }
  if (FLAG_STD)
  {
    IPTR_STD = dbout->addFields(nvar,0.);
    if (IPTR_STD < 0) goto label_end;
  }

  /* Pre-calculations */

  if (neigh_start(DBIN,neigh)) goto label_end;
  if (st_model_manage(1,model)) goto label_end;
  if (st_krige_manage(1,nvar,model,neigh)) goto label_end;
  if (krige_koption_manage(1,0,KOPTION_BLOCK,1,ndisc)) goto label_end;

  /* Loop on the targets to be processed */

  status = 0;
  for (IECH_OUT=0; IECH_OUT<get_NECH(DBOUT); IECH_OUT++)
  {
    mes_process("Kriging sample",get_NECH(DBOUT),IECH_OUT);
    debug_index(IECH_OUT+1);
    if (! dbout->isActive(IECH_OUT)) continue;
    if (debug_query("kriging") ||
        debug_query("nbgh")    ||
        debug_query("results"))
    {
      mestitle(1,"Target location");
      db_sample_print(dbout,IECH_OUT,1,0,0);
    }

    /* Update the discretization characteristics (per sample) */
    st_block_discretize(1,1,IECH_OUT);
    if (FLAG_STD) st_variance0(model,nvar,VectorDouble());

    /* Select the Neighborhood */

    flag_new_nbgh = st_neigh(neigh,&status,&nech);
    if (status) goto label_store;

    /* Establish the kriging L.H.S. */

    if (flag_new_nbgh || flag_continuous_kriging(neigh) || debug_force())
    {
      st_prepar(model,neigh,nech,&status,&nred,&neq);
      if (status) goto label_store;
      st_data_dual(model,NULL,nech,nred,&ldum);
    }

    /* Establish the kriging R.H.S. */

    st_rhs(model,nech,neq,nvar,NULL,&status);
    if (status) goto label_store;
    st_rhs_iso2hetero(neq,nvar);
    if (debug_query("kriging"))
      krige_rhs_print(nvar,nech,neq,nred,flag,rhs);

    /* Derive the kriging weights */

    if (FLAG_WGT)
    {
      matrix_product(nred,nred,nvar,lhs,rhs,wgt);
      if (debug_query("kriging"))
        krige_wgt_print(status,nvar,nvar,nfeq,nech,nred,-1,flag,wgt);
    }

    /* Perform the estimation */

  label_store:
    st_estimate(model,NULL,status,neigh->getFlagXvalid(),nech,nvar,nred);
    if (debug_query("results"))
      st_result_kriging_print(neigh->getFlagXvalid(),nvar,status);
  }

  /* Set the error return flag */

  error = 0;

label_end:
  debug_index(0);
  (void) st_model_manage(-1,model);
  (void) st_krige_manage(-1,nvar,model,neigh);
  (void) krige_koption_manage(-1,0,KOPTION_BLOCK,1,ndisc);
  (void) manage_external_info(-1,LOC_F,DBIN,DBOUT,&iext);
  (void) manage_external_info(-1,LOC_NOSTAT,DBIN,DBOUT,&inostat);
  neigh_stop();
  return(error);
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
static int st_calculate_hermite_factors(Db  *db,
                                        int  nfactor)
{
  VectorDouble hn;
  int iptr,error;

  /* Initializations */

  error = 1;

  /* Create the new variables */

  iptr = db->addFields(nfactor,0.);
  if (iptr < 0) goto label_end;

  /* Loop on the samples */

  for (int iech=0; iech<get_NECH(db); iech++)
  {
    if (! db->isActive(iech)) continue;

    /* Calculate the factors */

    hn = hermitePolynomials(db->getVariable(iech,0),1.,nfactor+1);

    /* Store the factors */

    for (int ih=0; ih<nfactor; ih++)
      db->setArray(iech,iptr+ih,hn[ih+1]);
  }

  /* Set the newly created variables to Z locator */

  db->setLocatorsByAttribute(nfactor,iptr,LOC_Z);

  /* Set the error return code */

  error = 0;

label_end:
  return(error);
}

/****************************************************************************/
/*!
**  Disjunctive Kriging
**
** \return  Error return code
**
** \param[in]  dbin      input Db structure (containing the factors)
** \param[in]  dbgrid    output Grid Db structure
** \param[in]  model     Model structure
** \param[in]  neigh     Neigh structrue
** \param[in]  nfactor   Number of factors to be estimated (0: all)
** \param[in]  nmult     Array of Multiplicity for Partition 
** \param[in]  ndisc     Discretization parameters (or NULL)
** \param[in]  flag_est  Option for the storing the estimation
** \param[in]  flag_std  Option for the storing the standard deviation
**
** \remark In case the Model handles a Ponctual Anamophossis, the
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
** \remark The value 'flag_calcul' must be set to:
** \remark - CALCUL_POINT for point-block estimation (flag_block = TRUE)
** \remark - CALCUL_BLOCK for point estimation
** \remark Nevertheless, for Point-Block model, if dbgrid is a grid of Panels
** \remark 'flag_calcul' is set to CALCUL_BLOCK to provoke the discretization
** \remark of Panel into SMUs.
**  
*****************************************************************************/
GEOSLIB_API int dk_f(Db *dbin,
                     Db *dbgrid,
                     Model *model,
                     Neigh *neigh,
                     int nfactor,
                     VectorInt nmult,
                     VectorInt ndisc,
                     int flag_est,
                     int flag_std)
{
  Db     *dbsmu;
  int     error,status,nech,neq,nred,nvar,nvarz,nfeq,iclass,imult,nb_mult;
  int     iptr_est_bck,iptr_std_bck,ivar,i;
  int    *varloc,flag_block,flag_panel,flag_continue,neqmax;
  double *rhs_cum,ldum;
  static double perturb = 1.e-8;

  /* Preliminary checks */

  error   = 1;
  iptr_est_bck = iptr_std_bck = -1;
  rhs_cum = (double *) NULL;
  varloc  = (int    *) NULL;
  dbsmu   = (Db     *) NULL;
  st_global_init(dbin,dbgrid);
  FLAG_EST = flag_est;
  FLAG_STD = flag_std;
  FLAG_WGT = flag_std;
  ModTrans& modtrs = model->getModTrans();

  // The model is not checked against the Data, as the number of variables
  // is not consistent: Model (1) whereas Data (nfactor-1)
  if (st_check_environment(1,1,NULL,neigh)) goto label_end;

  if (model->getVariableNumber() > 1)
  {
    messerr("This application is limited to the monovariate Model case");
    goto label_end;
  }
  if (model->getModTransMode() != MODEL_PROPERTY_ANAM)
  {
    messerr("When using Disjunctive Kriging, the Model be incremented");
    messerr("with Properties beforehad");
    messerr("For that sake, use 'model.properties'");
    goto label_end;
  }

  if (IFFFF(nfactor)) nfactor = modtrs.getAnamNClass();
  if (modtrs.getAnam()->getType() == ANAM_HERMITIAN)
  {
    /* In the gaussian case, calculate the 'nfactor-1' factors */

    if (! DBIN->isVariableNumberComparedTo(1))
    {
      messerr("In Gaussian case, Input File must contain a single variable");
      goto label_end;
    }
    if (st_calculate_hermite_factors(DBIN,nfactor-1)) goto label_end;
  }
  nvarz  = DBIN->getVariableNumber();

  if (nfactor-1 != nvarz)
  {
    messerr("The number of variables in Input Db (%d) does not match",nvarz);
    messerr("the number of factors (%d)-1",nfactor);
    goto label_end;
  }
  flag_block = 0;
  if (modtrs.getAnam()->getType() == ANAM_HERMITIAN)
  {
    AnamHermite* anam_hermite = dynamic_cast<AnamHermite*>(modtrs.getAnam());
    if (anam_hermite->getRCoef() < 1.) flag_block = 1;
  }
  else if (modtrs.getAnam()->getType() == ANAM_DISCRETE_DD)
  {
    AnamDiscreteDD* anam_discrete_DD = dynamic_cast<AnamDiscreteDD*>(modtrs.getAnam());
    if (anam_discrete_DD->getSCoef() > 0.) flag_block = 1;
  }
  else
  {
    AnamDiscreteIR* anam_discrete_IR = dynamic_cast<AnamDiscreteIR*>(modtrs.getAnam());
    if (anam_discrete_IR->getRCoef() < 1.) flag_block = 1;
  }
  flag_panel = (flag_block && ! nmult.empty());

  /* Add the attributes for storing the results */

  if (FLAG_EST)
  {
    iptr_est_bck = DBOUT->addFields(nfactor-1,TEST);
    if (iptr_est_bck < 0) goto label_end;
  }
  if (FLAG_STD)
  {
    iptr_std_bck = DBOUT->addFields(nfactor-1,TEST);
    if (iptr_std_bck < 0) goto label_end;
  }

  /* Pre-calculations */

  if (neigh_start(DBIN,neigh)) goto label_end;
  if (st_model_manage(1,model)) goto label_end;
  if (st_krige_manage(1,model->getVariableNumber(),model,neigh)) goto label_end;
  nb_mult = 1;
  if (flag_block)
  {
    if (flag_panel)
    {
      if (krige_koption_manage(1,1,KOPTION_BLOCK,1,nmult)) goto label_end;
      KOPTION->calcul = KOPTION_PONCTUAL;
      nb_mult = KOPTION->ntot;
    }
    else
    {
      if (krige_koption_manage(1,1,KOPTION_PONCTUAL,1,VectorInt())) goto label_end;
    }
  }    
  else
  {
    if (krige_koption_manage(1,1,KOPTION_BLOCK,0,ndisc)) goto label_end;
  }

  /* Centering the data */

  if (flag_block)
  {
    if (flag_panel)
    {
      dbsmu = db_create_grid_divider(dbgrid,nmult.data(),1);
      if (db_center_point_to_grid(DBIN,dbsmu,perturb)) goto label_end;
      dbsmu = db_delete(dbsmu);
    }
    else
    {
      if (db_center_point_to_grid(DBIN,dbgrid,perturb)) goto label_end;
    }
  }

  /* Core allocation */

  nvar   = model->getVariableNumber();
  nfeq   = model->getDriftEquationNumber();
  neqmax = 0;
  if (flag_panel)
  {
    neqmax  = nvar * st_get_nmax(neigh) + model->getDriftEquationNumber();
    rhs_cum = st_core(neqmax,1);
  }
  varloc = (int *) mem_alloc(sizeof(int) * nvarz,1);
  for (ivar=0; ivar<nvarz; ivar++) 
    varloc[ivar] = DBIN->getColumnByLocator(LOC_Z,ivar);

  /* Loop on the targets to be processed */
  
  status = 0;
  for (IECH_OUT=0; IECH_OUT<get_NECH(DBOUT); IECH_OUT++)
  {
    mes_process("Disjunctive Kriging for cell",get_NECH(DBOUT),IECH_OUT);
    debug_index(IECH_OUT+1);
    if (! DBOUT->isActive(IECH_OUT)) continue;
    if (debug_query("kriging") ||
        debug_query("nbgh")    ||
        debug_query("results"))
    {
      mestitle(1,"Target location");
      db_sample_print(DBOUT,IECH_OUT,1,0,0);
    }
    
    /* Select the Neighborhood */
    
    DBIN->clearLocators(LOC_Z);
    DBIN->setLocatorByAttribute(varloc[0],LOC_Z);
    (void) st_neigh(neigh,&status,&nech);
    if (status) continue;

    /* Loop on the effective factors */
    
    flag_continue = 1;
    for (iclass=1; iclass<nfactor && flag_continue; iclass++)
    {
      if (FLAG_EST) IPTR_EST = iptr_est_bck + iclass - 1;
      if (FLAG_STD) IPTR_STD = iptr_std_bck + iclass - 1;
      DBIN->clearLocators(LOC_Z);
      DBIN->setLocatorByAttribute(varloc[iclass-1],LOC_Z);
      
      /* Set the rank of the current factor in the model */
      
      if (model_anamorphosis_set_factor(model,iclass)) goto label_end;
      
      /* Constant Calculation for the variance */
      
      if (FLAG_STD) st_variance0(model,nvar,VectorDouble());
      
      /* Establish the kriging L.H.S. (always performed) */

      st_prepar(model,neigh,nech,&status,&nred,&neq);
      if (status) 
      {
        flag_continue = 0;
        continue;
      }
      st_data_dual(model,NULL,nech,nred,&ldum);

      /* Blank out the estimate and the R.H.S. */

      if (flag_panel)
        for (i=0; i<neqmax; i++) rhs_cum[i] = 0.;

      /* Loop on Panels (*1) or SMU (*nb_mult) */

      for (imult=0; imult<nb_mult; imult++)
      {
	
        /* Establish the kriging R.H.S. */

        if (flag_panel) RAND_INDEX = imult;
        st_rhs(model,nech,neq,nvar,NULL,&status);
        if (status) goto label_store;
        st_rhs_iso2hetero(neq,nvar);
	
        /* Cumulate the R.H.S. */

        if (flag_panel)
          for (i=0; i<neqmax; i++) rhs_cum[i] += rhs[i];
      }

      /* Scale the R.H.S. */

      if (flag_panel)
        for (i=0; i<neqmax; i++) rhs[i] = rhs_cum[i] / (double) nb_mult;
      if (debug_query("kriging"))
        krige_rhs_print(nvar,nech,neq,nred,flag,rhs);

      /* Derive the kriging weights */
	
      if (FLAG_WGT)
      {
        matrix_product(nred,nred,nvar,lhs,rhs,wgt);
        if (debug_query("kriging"))
          krige_wgt_print(status,nvar,nvar,nfeq,nech,nred,-1,flag,wgt);
      }
	
      /* Perform the estimation */
	
    label_store:
      st_estimate(model,NULL,status,neigh->getFlagXvalid(),nech,nvar,nred);
      if (debug_query("results"))
        st_result_kriging_print(neigh->getFlagXvalid(),nvar,status);
    }
  }

  /* Set the error return flag */

  error = 0;

label_end:
  debug_index(0);
  rhs_cum = (double *) mem_free((char *) rhs_cum);
  varloc  = (int    *) mem_free((char *) varloc);
  (void) st_model_manage(-1,model);
  (void) st_krige_manage(-1,model->getVariableNumber(),model,neigh);
  (void) krige_koption_manage(-1,1,KOPTION_PONCTUAL,1,ndisc);
  neigh_stop();
  return(error);
}

/****************************************************************************/
/*!
**  Perform the Neighborhood search
**
** \return  Array of sample indices of the target neighbors
**
** \param[in]  dbin      input Db structure
** \param[in]  model     Model structrue
** \param[in]  neigh     Neigh structrue
** \param[in]  target    Target location
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
GEOSLIB_API int *neigh_calc(Db     *dbin,
                            Model  *model,
                            Neigh  *neigh,
                            double *target,
                            int    *nech_out)
{
  int *neigh_tab,i,error,status,nech,zloc;
  Db  *dbout;

  /* Preliminary checks */

  neigh_tab = (int *) NULL;
  dbout = (Db *) NULL;
  *nech_out = 0;
  error =  1;

  /* Create a temporary dummy Db which contains the target */

  if (model == (Model *) NULL) goto label_end;
  dbout = db_create_from_target(target,model->getDimensionNumber(),1);
  if (dbout == (Db *) NULL) goto label_end;
  st_global_init(dbin,dbout);

  /* Modification of 'dbin' */

  if (dbin != (Db *) NULL && model != (Model *) NULL &&
      dbin->getVariableNumber() != model->getVariableNumber() && model->getVariableNumber() == 1)
  {
    zloc = dbin->getColumnByLocator(LOC_Z,0);
    dbin->clearLocators(LOC_Z);
    dbin->setLocatorByAttribute(zloc,LOC_Z);
  }
  if (st_check_environment(1,1,model,neigh)) goto label_end;

  /* Pre-calculations */

  if (neigh_start(DBIN,neigh)) goto label_end;
  if (st_model_manage(1,model)) goto label_end;
  if (st_krige_manage(1,model->getVariableNumber(),model,neigh)) goto label_end;

  /* Select the Neighborhood */

  IECH_OUT = 0;
  (void) st_neigh(neigh,&status,&nech);
  if (status != 0) 
  {
    messerr("Neighborhood search failed");
    goto label_end;
  }

  /* Store the neighbor indices */

  neigh_tab = (int *) mem_alloc(sizeof(int) * nech,1);
  for (i=0; i<nech; i++) neigh_tab[i] = rank[i] + 1;
  *nech_out = nech;

  /* Set the error return flag */

  error = 0;

label_end:
  if (error) neigh_tab = (int *) mem_free((char *) neigh_tab);
  dbout = db_delete(dbout);
  (void) st_model_manage(-1,model);
  (void) st_krige_manage(-1,model->getVariableNumber(),model,neigh);
  neigh_stop();
  return(neigh_tab);
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
static int *st_ranks_other(int  nech,
                           int  nsize1,
                           int *ranks1,
                           int  nsize2,
                           int *ranks2)
{
  int *rother,i;

  rother = (int *) mem_alloc(sizeof(int) * nech,0);
  if (rother == (int *) NULL) return(rother);

  for (i=0; i<nech; i++) rother[i] = i;
  for (i=0; i<nsize1; i++) 
    rother[ranks1[i]] = -1;
  for (i=0; i<nsize2; i++) 
    rother[ranks2[i]] = -1;
  return(rother);
}

/****************************************************************************/
/*!
**  Establishing the kriging system with exact and ACP points
**
** \return  Error retun code
**
** \param[in]  db        Db structure
** \param[in]  model     Model structrue
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
static int st_sampling_krige_data(Db      *db,
                                  Model   *model,
                                  double   beta,
                                  int      nsize1,
                                  int     *ranks1,
                                  int      nsize2,
                                  int     *ranks2,
                                  int     *rother,
                                  int     *ntot_arg,
                                  int     *nutil_arg,
                                  int    **rutil_arg,
                                  double **tutil_arg,
                                  double **invsig_arg)
{
  int    *isort,*ralls,*rutil;
  int     ndat,error,i,j,ntot,ntri,nother,npart,n1,ecr,nutil,nmax;
  double *utab,*s,*tl,*xl,*c,*sq,*v,*tn1,*tn2,*eigval,*eigvec,*spart,*vsort;
  double *tutil,*invsig,sumval;

  /* Initializations */

  error  = 1;
  utab   = s = tl = xl = c = v = tn1 = tn2 = sq = (double *) NULL;
  eigval = eigvec = spart = vsort = tutil = invsig = (double *) NULL;
  isort  = ralls = rutil = (int *) NULL;
  ndat   = db->getActiveSampleNumber();
  ntot   = nsize1 + nsize2;
  ntri   = nsize2 * (nsize2 + 1) / 2;
  nother = ndat - ntot;
  npart  = ndat - nsize1;
  nutil  = 0;

  /* Core allocation */

  utab = (double *) mem_alloc(sizeof(double) * ndat * ntot,0);
  if (utab == (double *) NULL) goto label_end;
  for (i=0; i<ndat * ntot; i++) utab[i] = 0.;
  ralls = (int *) mem_alloc(sizeof(int) * ndat,0);
  if (ralls == (int *) NULL) goto label_end;

  /* Defining 'utab' for exact pivots */

  for (i=0; i<nsize1; i++) UTAB(i,i) = 1.;
  ecr = 0;
  for (i=0; i<nsize1; i++) ralls[ecr++] = ranks1[i];
  for (i=0; i<nsize2; i++) ralls[ecr++] = ranks2[i];
  for (i=0; i<ndat; i++)
    if (rother[i] >= 0) ralls[ecr++] = rother[i];
    
  /* Defining 'utab' for ACP pivots */

  if (nsize2 > 0)
  {
    tl  = (double *) mem_alloc(sizeof(double) * ntri,0);
    if (tl  == (double *) NULL) goto label_end;
    xl  = (double *) mem_alloc(sizeof(double) * ntri,0);
    if (xl  == (double *) NULL) goto label_end;
    v   = (double *) mem_alloc(sizeof(double) * nother * nsize2,0);
    if (v   == (double *) NULL) goto label_end;
    sq = (double *) mem_alloc(sizeof(double) * nsize2 * nsize2,0);
    if (sq == (double *) NULL) goto label_end;
    tn1 = (double *) mem_alloc(sizeof(double) * nsize2 * nsize2,0);
    if (tn1 == (double *) NULL) goto label_end;
    tn2 = (double *) mem_alloc(sizeof(double) * nsize2 * nsize2,0);
    if (tn2 == (double *) NULL) goto label_end;
    eigval = (double *) mem_alloc(sizeof(double) * nsize2,0);
    if (eigval == (double *) NULL) goto label_end;
    eigvec = (double *) mem_alloc(sizeof(double) * nsize2 * nsize2,0);
    if (eigvec == (double *) NULL) goto label_end;
    if (beta > 0.) 
    {
      vsort = (double *) mem_alloc(sizeof(double) * npart,0);
      if (vsort == (double *) NULL) goto label_end;
      isort = (int    *) mem_alloc(sizeof(double) * npart,0);
      if (isort == (int    *) NULL) goto label_end;
    }

    s = model_covmat_by_ranks(model,db,nsize2,ranks2,db,nsize2,ranks2,
                              -1,-1,0,1);
    if (s == (double *) NULL) goto label_end;
    if (matrix_cholesky_decompose(s,tl,nsize2)) goto label_end;
    matrix_triangle_to_square(0,nsize2,tl,sq);
    matrix_cholesky_invert(nsize2,tl,xl);
    c = model_covmat_by_ranks(model,db,nsize2,ranks2,db,ndat,rother,
                              -1,-1,0,1);
    if (c == (double *) NULL) goto label_end;
    matrix_cholesky_product(4,nsize2,nother,xl,c,v);
    matrix_cholesky_norme(1,nsize2,tl,(double *) NULL,tn1);
    if (matrix_prod_norme(-1,nother,nsize2,v,NULL,tn2)) goto label_end;
    matrix_combine(nsize2 * nsize2,1,tn1,1,tn2,tn1);
    if (matrix_eigen(tn1,nsize2,eigval,eigvec)) goto label_end;
    matrix_product_by_diag(3,nsize2,eigvec,eigval,eigvec);
    spart = matrix_bind(1,nsize2,nsize2,sq,nother,nsize2,v,&npart,&n1);
    if (spart == (double *) NULL) goto label_end;
    matrix_product(npart,nsize2,nsize2,spart,eigvec,spart);

    if (beta > 0.)
    {
      for (i=0; i<npart; i++)
      {
        sumval = 0.;
        for (j=0; j<nsize2; j++) sumval = MAX(sumval,ABS(SPART(i,j)));
        vsort[i] = sumval;
        isort[i] = i;
      }
      ut_sort_double(1,npart,isort,vsort);
      nmax = MIN(npart,(int) (beta * (double) npart));
      for (i=0; i<nmax; i++)
        for (j=0; j<nsize2; j++) SPART(isort[i],j) = 0.;
    }

    for (i=0; i<npart; i++)
      for (j=0; j<nsize2; j++)
        UTAB(i+nsize1,j+nsize1) = -SPART(i,j);
    
    /* Core deallocation */
    
    tl     = (double *) mem_free((char *) tl);
    xl     = (double *) mem_free((char *) xl);
    v      = (double *) mem_free((char *) v);
    s      = (double *) mem_free((char *) s);
    c      = (double *) mem_free((char *) c);
    sq     = (double *) mem_free((char *) sq);
    tn1    = (double *) mem_free((char *) tn1);
    tn2    = (double *) mem_free((char *) tn2);
    spart  = (double *) mem_free((char *) spart);
    vsort  = (double *) mem_free((char *) vsort);
    eigval = (double *) mem_free((char *) eigval);
    eigvec = (double *) mem_free((char *) eigvec);
    isort  = (int    *) mem_free((char *) isort);
  }

  /* Count the number of active samples */

  nutil = 0;
  for (i=0; i<ndat; i++)
  {
    sumval = 0.;
    for (j=0; j<ntot; j++) sumval += ABS(UTAB(i,j));
    if (sumval > 0.) nutil++;
  }

  /* Create the output arrays */
  
  rutil = (int    *) mem_alloc(sizeof(int) * nutil,0);
  if (rutil == (int    *) NULL) goto label_end;
  tutil = (double *) mem_alloc(sizeof(double) * ntot * nutil,0);
  if (tutil == (double *) NULL) goto label_end;
  invsig = (double *) mem_alloc(sizeof(double) * ntot * ntot,0);
  if (invsig == (double *) NULL) goto label_end;

  for (i=ecr=0; i<ndat; i++)
  {
    sumval = 0.;
    for (j=0; j<ntot; j++) sumval += ABS(UTAB(i,j));
    if (sumval <= 0.) continue;
    rutil[ecr] = ralls[i];
    for (j=0; j<ntot; j++) TUTIL(ecr,j) = UTAB(i,j);
    ecr++;
  }
  s = model_covmat_by_ranks(model,db,nutil,rutil,db,nutil,rutil,
                            -1,-1,0,1);
  if (s == (double *) NULL) goto label_end;
  if (matrix_prod_norme(-1,nutil,ntot,tutil,s,invsig)) goto label_end;
  if (matrix_invert(invsig,ntot,0)) goto label_end;
  s      = (double *) mem_free((char *) s);
  utab   = (double *) mem_free((char *) utab);
  ralls  = (int    *) mem_free((char *) ralls);

  /* Returning arguments */

  *ntot_arg   = ntot;
  *nutil_arg  = nutil;
  *rutil_arg  = rutil;
  *tutil_arg  = tutil;
  *invsig_arg = invsig;

  /* Error return code */

  error = 0;

label_end:
  utab   = (double *) mem_free((char *) utab);
  tl     = (double *) mem_free((char *) tl);
  xl     = (double *) mem_free((char *) xl);
  v      = (double *) mem_free((char *) v);
  s      = (double *) mem_free((char *) s);
  c      = (double *) mem_free((char *) c);
  sq     = (double *) mem_free((char *) sq);
  tn1    = (double *) mem_free((char *) tn1);
  tn2    = (double *) mem_free((char *) tn2);
  spart  = (double *) mem_free((char *) spart);
  vsort  = (double *) mem_free((char *) vsort);
  eigval = (double *) mem_free((char *) eigval);
  eigvec = (double *) mem_free((char *) eigvec);
  isort  = (int    *) mem_free((char *) isort);
  ralls  = (int    *) mem_free((char *) ralls);
  return(error);
}

/****************************************************************************/
/*!
**  Perform the estimation at the data points
**
** \return  Error retun code
**
** \param[in]  db         Db structure
** \param[in]  model      Model structrue
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
GEOSLIB_API int st_krige_data(Db     *db,
                              Model  *model,
                              double  beta,
                              int     nsize1,
                              int    *ranks1,
                              int     nsize2,
                              int    *ranks2,
                              int    *rother,
                              int     flag_abs,
                              double *data_est,
                              double *data_var)
{
  int    *rutil,error,ntot,nutil,i,iech,nech;
  double *data,*tutil,*invsig,*s,*datm,*aux1,*aux2,*aux3,*aux4,*c00;
  double  estim,variance,true_value;

  /* Initializations */

  error  = 1;
  rutil  = (int    *) NULL;
  tutil  = invsig = data = datm = s = c00 = (double *) NULL;
  aux1   = aux2 = aux3 = aux4 = (double *) NULL;

  /* Core allocation */

  nutil = ntot = 0;
  nech  = get_NECH(db);
  data  = db_vector_alloc(db);
  if (data == (double *) NULL) goto label_end;
  
  /* Perform local sampling */
  
  if (st_sampling_krige_data(db,model,beta,nsize1,ranks1,nsize2,ranks2,rother,
                             &ntot,&nutil,
                             &rutil,&tutil,&invsig)) goto label_end;
  
  /* Second core allocation */
  
  datm = (double *) mem_alloc(sizeof(double) * nutil,0);
  if (datm == (double *) NULL) goto label_end;
  aux1 = (double *) mem_alloc(sizeof(double) * ntot,0);
  if (aux1 == ( double *) NULL) goto label_end;
  aux2 = (double *) mem_alloc(sizeof(double) * ntot,0);
  if (aux2 == ( double *) NULL) goto label_end;
  aux3 = (double *) mem_alloc(sizeof(double) * ntot,0);
  if (aux3 == ( double *) NULL) goto label_end;
  aux4 = (double *) mem_alloc(sizeof(double) * ntot,0);
  if (aux4 == ( double *) NULL) goto label_end;
  
  /* Get the vector of active data and substract the mean */
  
  if (db_vector_get(db,LOC_Z,0,data)) goto label_end;
  for (i=0; i<nutil; i++) datm[i] = data[rutil[i]] - model->getMean(0);
  matrix_product(1,nutil,ntot,datm,tutil,aux1);
  matrix_product(1,ntot,ntot,aux1,invsig,aux2);
  
  /* Perform the estimation at all non pivot samples */

  for (iech=0; iech<nech; iech++)
  {
    data_est[iech] = data_var[iech] = TEST;
    if (! db->isActive(iech)) continue;
    if (rother[iech] < 0) continue;
    c00 = model_covmat_by_ranks(model,db,1,&iech,db,1,&iech,-1,-1,0,1);
    if (c00 == (double *) NULL) goto label_end;
    s = model_covmat_by_ranks(model,db,nutil,rutil,db,1,&iech,-1,-1,0,1);
    if (s == (double *) NULL) goto label_end;

    matrix_product(1,nutil,ntot,s,tutil,aux3);
    matrix_product(1,ntot,1,aux2,aux3,&estim);
    data_est[iech] = estim + model->getMean(0);

    if (flag_abs)
    {
      true_value = db->getVariable(iech,0);
      if (FFFF(true_value))
        data_est[iech] = TEST;
      else
        data_est[iech] = ABS(data_est[iech] - true_value);
    }

    matrix_product(1,ntot,ntot,aux3,invsig,aux4);
    matrix_product(1,ntot,1,aux3,aux4,&variance);
    data_var[iech] = c00[0] - variance;

    s      = (double *) mem_free((char *) s);
    c00    = (double *) mem_free((char *) c00);
  }

  /* Error return code */

  error = 0;

label_end:
  data = db_vector_free(data);
  rutil  = (int    *) mem_free((char *) rutil);
  tutil  = (double *) mem_free((char *) tutil);
  invsig = (double *) mem_free((char *) invsig);
  datm   = (double *) mem_free((char *) datm);
  s      = (double *) mem_free((char *) s);
  c00    = (double *) mem_free((char *) c00);
  aux1   = (double *) mem_free((char *) aux1);
  aux2   = (double *) mem_free((char *) aux2);
  aux3   = (double *) mem_free((char *) aux3);
  aux4   = (double *) mem_free((char *) aux4);
  return(error);
}

/****************************************************************************/
/*!
**  Evaluate the improvement in adding a new pivot on the global score
**
** \return  Error retun code
**
** \param[in]  db         Db structure
** \param[in]  model      Model structrue
** \param[in]  nsize1     Number of exact pivots currently selected
** \param[in]  ranks1     Ranks of exact pivots
** \param[in]  rother     Ranks of the idle samples
**
** \param[out] crit       Array of criterion
**
*****************************************************************************/
GEOSLIB_API int st_crit_global(Db     *db,
                               Model  *model,
                               int     nsize1,
                               int    *ranks1,
                               int    *rother,
                               double *crit)
{
  int     error,ndat,i,iech,nutil,ecr;
  double *c00,*invc,*data,*datm,*cs,*temp,*olderr,*olddiv,*aux1,*cs1;
  double *temp_loc,estim,sigma,value;

  /* Initializations */

  error  = 1;
  ndat   = db->getActiveSampleNumber();
  nutil  = ndat - nsize1;
  c00    = invc = data = datm = cs = temp = olderr = olddiv = (double *) NULL;
  aux1   = cs1 = (double *) NULL;

  /* Preliminary checks */

  if (nsize1 <= 0) goto label_end;

  /* Core allocation */

  data  = db_vector_alloc(db);
  if (data == (double *) NULL) goto label_end;
  datm = (double *) mem_alloc(sizeof(double) * ndat,0);
  if (datm == (double *) NULL) goto label_end;
  olderr = (double *) mem_alloc(sizeof(double) * nutil,0);
  if (olderr == (double *) NULL) goto label_end;
  olddiv = (double *) mem_alloc(sizeof(double) * nutil,0);
  if (olddiv == (double *) NULL) goto label_end;
  temp   = (double *) mem_alloc(sizeof(double) * nsize1 * nutil,0);
  if (temp   == (double *) NULL) goto label_end;
  aux1   = (double *) mem_alloc(sizeof(double) * nutil,0);
  if (aux1   == (double *) NULL) goto label_end;

  /* Establish the Kriging matrix on the pivot samples */

  invc = model_covmat_by_ranks(model,db,nsize1,ranks1,db,nsize1,ranks1,
                               -1,-1,0,1);
  if (invc == (double *) NULL) goto label_end;
  if (matrix_invert(invc,nsize1,0)) goto label_end;

  /* Set the data vector (corrected by the mean */

  if (db_vector_get(db,LOC_Z,0,data)) goto label_end;
  for (i=0; i<nsize1; i++) datm[i] = data[ranks1[i]] - model->getMean(0);

  /* Loop on the non-pivots */

  for (iech=ecr=0; iech<ndat; iech++)
  {
    temp_loc = &temp[ecr * nsize1];
    if (! db->isActive(iech)) continue;
    if (rother[iech] < 0) continue;

    c00 = model_covmat_by_ranks(model,db,1,&iech,db,1,&iech,-1,-1,0,1);
    if (c00 == (double *) NULL) goto label_end;

    cs = model_covmat_by_ranks(model,db,nsize1,ranks1,db,1,&iech,-1,-1,0,1);
    if (cs == (double *) NULL) goto label_end;
  
    matrix_product(nsize1,nsize1,1,invc,cs,temp_loc);
    matrix_product(1,nsize1,1,datm,temp_loc,&estim);
    olderr[ecr] = estim + model->getMean(0) - db->getVariable(iech,0);
    
    matrix_product(1,nsize1,1,cs,temp_loc,&sigma);
    olddiv[ecr] = olderr[ecr] / (c00[0] - sigma);

    c00 = (double *) mem_free((char *) c00);
    cs  = (double *) mem_free((char *) cs);
    ecr++;
  }

  /* Loop on the candidates */

  for (iech=ecr=0; iech<ndat; iech++)
  {
    crit[iech] = TEST;
    if (! db->isActive(iech)) continue;
    if (rother[iech] < 0) continue;

    cs = model_covmat_by_ranks(model,db,1,&iech,db,nsize1,ranks1,-1,-1,0,1);
    if (cs == (double *) NULL) goto label_end;

    cs1 = model_covmat_by_ranks(model,db,1,&iech,db,ndat,rother,-1,-1,0,1);
    if (cs1 == (double *) NULL) goto label_end;

    matrix_product(1,nsize1,nutil,cs,temp,aux1);
    matrix_combine(nutil,1,cs1,-1,aux1,cs1);
    matrix_combine(nutil,1,olderr,-olddiv[ecr],cs1,cs1);

    value = 0.;
    for (i=0; i<nutil; i++) value += cs1[i] * cs1[i];
    crit[iech] = value / nutil;

    cs   = (double *) mem_free((char *) cs);
    cs1  = (double *) mem_free((char *) cs1);
    ecr++;
  }

  /* Set the error return code */

  error = 0;

label_end:
  data   = db_vector_free(data);
  c00    = (double *) mem_free((char *) c00);
  invc   = (double *) mem_free((char *) invc);
  datm   = (double *) mem_free((char *) datm);
  cs     = (double *) mem_free((char *) cs);
  cs1    = (double *) mem_free((char *) cs1);
  temp   = (double *) mem_free((char *) temp);
  aux1   = (double *) mem_free((char *) aux1);
  olderr = (double *) mem_free((char *) olderr);
  olddiv = (double *) mem_free((char *) olddiv);
  return(error);
}

/****************************************************************************/
/*!
**  Optimize the sampling design
**
** \return  Error retun code
**
** \param[in]  db         Db structure
** \param[in]  model      Model structrue
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
GEOSLIB_API int sampling_f(Db     *db,
                         Model  *model,
                         double  beta,
                         int     method1,
                         int     nsize1_max,
                         int     nsize1,
                         int    *ranks1,
                         int     method2,
                         int     nsize2_max,
                         int     nsize2,
                         int    *ranks2,
                         int     verbose)
{
  int    *rother,error,best_rank,nech,nval;
  double *data_est,*data_var;
  double  best_ecart,minimum,maximum,mean,stdv,delta;

  /* Initializations */

  error    = 1;
  data_est = data_var = (double *) NULL;
  rother   = (int *) NULL;
  nech     = get_NECH(db);

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
  if (data_est == (double *) NULL) goto label_end;
  data_var = db_vector_alloc(db);
  if (data_var == (double *) NULL) goto label_end;
  rother = st_ranks_other(nech,nsize1,ranks1,nsize2,ranks2);
  if (rother == (int *) NULL) goto label_end;
  
  /* Sample the exact pivots */

  while (nsize1 < nsize1_max)
  {
    if (method1 == 1)
    {
      if (st_krige_data(db,model,beta,nsize1,ranks1,nsize2,ranks2,rother,
                        1,data_est,data_var)) goto label_end;
      best_rank  = matrix_get_extreme(2,nech,data_est);
      best_ecart = data_est[best_rank];
    }
    else
    {
      if (st_crit_global(db,model,nsize1,ranks1,rother,
                         data_est)) goto label_end;
      best_rank  = matrix_get_extreme(1,nech,data_est);
      best_ecart = data_est[best_rank];
    }
    if (verbose)
      message("Exact Pivots (%3d/%3d): Rank = %3d - value = %lf\n",
              nsize1+1,nsize1_max,best_rank+1,best_ecart);
    ranks1[nsize1] = best_rank;
    rother[best_rank] = -1;
    nsize1++;
  }

  /* Sample the ACP pivots */

  while (nsize2 < nsize2_max)
  {
    if (st_krige_data(db,model,beta,nsize1,ranks1,nsize2,ranks2,rother,
                      1,data_est,data_var)) goto label_end;
    best_rank  = matrix_get_extreme(2,nech,data_est);
    best_ecart = data_est[best_rank];
    if (verbose)
      message("ACP   Pivots (%3d/%3d): Rank = %3d - value = %lf\n",
              nsize2+1,nsize2_max,best_rank+1,best_ecart);
    ranks2[nsize2] = best_rank;
    rother[best_rank] = -1;
    nsize2++;
  }

  /* Calculation of statistics on reproduction errors */

  if (verbose)
  {
    if (st_krige_data(db,model,beta,nsize1,ranks1,nsize2,ranks2,rother,
                      1,data_est,data_var)) goto label_end;
    ut_statistics(nech,data_est,NULL,NULL,
                  &nval,&minimum,&maximum,&delta,&mean,&stdv);
    mestitle(1,"Statistics on estimation errors");
    message("Count   = %d \n",nval);
    message("Minimum = %lf\n",minimum);
    message("Mean    = %lf\n",mean);
    message("St. Dev.= %lf\n",stdv);
    message("Maximum = %lf\n",maximum);
  }
    
  /* Error return code */

  error = 0;

label_end:
  data_est = db_vector_free(data_est);
  data_var = db_vector_free(data_var);
  rother   = (int *) mem_free((char *) rother);
  return(error);
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
** \param[in]  model      Model structrue
** \param[in]  beta       Thresholding value
** \param[in]  nsize1     Number of exact pivots currently selected
** \param[in]  ranks1     Ranks of exact pivots
** \param[in]  nsize2     Number of ACP pivots currently selected
** \param[in]  ranks2     Ranks of ACP pivots
** \param[in]  flag_std   Option for storing the standard deviation
** \param[in]  verbose    Verbose flag
**
*****************************************************************************/
GEOSLIB_API int krigsampling_f(Db     *dbin,
                             Db     *dbout,
                             Model  *model,
                             double  beta,
                             int     nsize1,
                             int    *ranks1,
                             int     nsize2,
                             int    *ranks2,
                             int     flag_std,
                             int     verbose)
{
  int    *rutil,*rother,error,nvar,ntot,nutil,i,nech;
  double *tutil,*data,*invsig,*datm,*aux1,*aux2,*aux3,*aux4,*s,*c00;
  double  estim,sigma;

  /* Preliminary checks */

  error =  1;
  rutil = rother = (int    *) NULL;
  tutil = invsig = data = datm = s = c00 = (double *) NULL;
  aux1  = aux2 = aux3 = aux4 = (double *) NULL;
  st_global_init(dbin,dbout);
  FLAG_EST  = 1;
  FLAG_STD  = flag_std;
  if (st_check_environment(1,1,model,NULL)) goto label_end;
  nvar = model->getVariableNumber();
  nech = get_NECH(dbin);

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

  if (FLAG_EST)
  {
    IPTR_EST = dbout->addFields(nvar,0.);
    if (IPTR_EST < 0) goto label_end;
  }
  if (FLAG_STD)
  {
    IPTR_STD = dbout->addFields(nvar,0.);
    if (IPTR_STD < 0) goto label_end;
  }

  /* Core allocation */

  rother = st_ranks_other(nech,nsize1,ranks1,nsize2,ranks2);
  if (rother == (int *) NULL) goto label_end;

  /* Perform local sampling */
  
  if (st_sampling_krige_data(dbin,model,
                             beta,nsize1,ranks1,nsize2,ranks2,rother,
                             &ntot,&nutil,
                             &rutil,&tutil,&invsig)) goto label_end;
  
  /* Optional printout */

  if (verbose)
  {
    message("Printout of intermediate arrays\n");
    print_imatrix("Pivot ranks",0,1,1,ntot,NULL,rutil);
    print_matrix("Inv-Sigma",0,1,ntot,ntot,NULL,invsig);
    print_matrix("U",0,1,ntot,nutil,NULL,tutil);
  }

  /* Second core allocation */
  
  data = db_vector_alloc(dbin);
  if (data == (double *) NULL) goto label_end;
  datm = (double *) mem_alloc(sizeof(double) * nutil,0);
  if (datm == (double *) NULL) goto label_end;
  aux1 = (double *) mem_alloc(sizeof(double) * ntot,0);
  if (aux1 == ( double *) NULL) goto label_end;
  aux2 = (double *) mem_alloc(sizeof(double) * ntot,0);
  if (aux2 == ( double *) NULL) goto label_end;
  aux3 = (double *) mem_alloc(sizeof(double) * ntot,0);
  if (aux3 == ( double *) NULL) goto label_end;
  if (FLAG_STD)
  {
    aux4 = (double *) mem_alloc(sizeof(double) * ntot,0);
    if (aux4 == ( double *) NULL) goto label_end;
  }

  /* Get the vector of active data and substract the mean */
  
  if (db_vector_get(dbin,LOC_Z,0,data)) goto label_end;
  for (i=0; i<nutil; i++) datm[i] = data[rutil[i]] - model->getMean(0);
  matrix_product(1,nutil,ntot,datm,tutil,aux1);
  matrix_product(1,ntot,ntot,aux1,invsig,aux2);
  
  /* Loop on the target samples */

  for (IECH_OUT=0; IECH_OUT<get_NECH(DBOUT); IECH_OUT++)
  {
    mes_process("Kriging sample",get_NECH(DBOUT),IECH_OUT);
    debug_index(IECH_OUT+1);
    if (! dbout->isActive(IECH_OUT)) continue;
    if (debug_query("kriging") ||
        debug_query("nbgh")    ||
        debug_query("results"))
    {
      mestitle(1,"Target location");
      db_sample_print(dbout,IECH_OUT,1,0,0);
    }

    s = model_covmat_by_ranks(model,dbin,nutil,rutil,dbout,1,&IECH_OUT,
                              -1,-1,0,1);
    if (s == (double *) NULL) goto label_end;
    if (FLAG_STD)
    {
      c00 = model_covmat_by_ranks(model,dbout,1,&IECH_OUT,dbout,1,&IECH_OUT,
                                  -1,-1,0,1);
      if (c00 == (double *) NULL) goto label_end;
    }

    matrix_product(1,nutil,ntot,s,tutil,aux3);
    matrix_product(1,ntot,1,aux2,aux3,&estim);
    estim += model->getMean(0);
    DBOUT->setArray(IECH_OUT,IPTR_EST,estim);

    if (FLAG_STD)
    {
      matrix_product(1,ntot,ntot,aux3,invsig,aux4);
      matrix_product(1,ntot,1,aux3,aux4,&sigma);
      sigma = c00[0] - sigma;
      sigma = (sigma > 0) ? sqrt(sigma) : 0.;
      DBOUT->setArray(IECH_OUT,IPTR_STD,sigma);
    }

    /* Optional printout */

    if (debug_query("results"))
    {
      tab_printg(" - Estimate  = ",1,GD_J_RIGHT,estim);
      message("\n");
      if (FLAG_STD) 
      {
        tab_printg(" - Std. Dev. = ",1,GD_J_RIGHT,sigma);
        message("\n");
      }
    }

    /* Core deallocation */

    s      = (double *) mem_free((char *) s);
    c00    = (double *) mem_free((char *) c00);
  }

  /* Error return code */

  error = 0;

label_end:
  rother = (int    *) mem_free((char *) rother);
  rutil  = (int    *) mem_free((char *) rutil);
  tutil  = (double *) mem_free((char *) tutil);
  invsig = (double *) mem_free((char *) invsig);
  data   = (double *) mem_free((char *) data);
  datm   = (double *) mem_free((char *) datm);
  s      = (double *) mem_free((char *) s);
  c00    = (double *) mem_free((char *) c00);
  aux1   = (double *) mem_free((char *) aux1);
  aux2   = (double *) mem_free((char *) aux2);
  aux3   = (double *) mem_free((char *) aux3);
  aux4   = (double *) mem_free((char *) aux4);
  return(error);
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
static void st_declustering_stats(int mode,
                                  int method,
                                  Db *db,
                                  int iptr)
{
  double mean,var,zval,coeff,mini,maxi,sumwgt;

  mean = var = sumwgt = 0.;
  mini =  1.e30;
  maxi = -1.e30;
  for (int iech=0; iech<get_NECH(db); iech++)
  {
    if (! db->isActive(iech)) continue;
    zval = db->getVariable(iech,0);
    if (FFFF(zval)) continue;
    coeff   = (mode == 0) ? 1. : get_ARRAY(db,iech,iptr);
    coeff   = ABS(coeff);
    sumwgt += coeff;
    mean   += coeff * zval;
    var    += coeff * zval * zval;
    if (coeff < mini) mini = coeff;
    if (coeff > maxi) maxi = coeff;
  }

  /* Scaling */

  if (sumwgt > 0)
  {
    mean /= sumwgt;
    var   = var / sumwgt - mean * mean;
  }

  if (mode == 0)
    mestitle(1,"Statistics before Declustering");
  else
    mestitle(1,"Statistics after Declustering");
  if (method == 1)
    message("- Using the Number of Samples per Neighborhood\n");
  else if (method == 2)
    message("- Using the weights for Kriging the Global Mean\n");
  else
    message("- Using the average weight for Kriging cells of a Grid\n");
  
  message("- Sum of weights    = %lf\n",sumwgt);
  message("- Mean              = %lf\n",mean);
  message("- Variance          = %lf\n",var);
  if (mode == 1)
  {
    message("- Minimum Weight    = %lf\n",mini);
    message("- Maximum Weight    = %lf\n",maxi);
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
static void st_declustering_truncate(Db     *db,
                                     int     iptr,
                                     int     verbose)
{
  double total,coeff;

  /* Truncate the negative weights */

  total = 0;
  for (int iech=0; iech<get_NECH(db); iech++)
  {
    if (! db->isActive(iech)) continue;
    if (FFFF(db->getVariable(iech,0))) continue;
    coeff = get_ARRAY(db,iech,iptr);
    if (coeff < 0) 
    {
      if (verbose) 
        messerr("Weight #%d is negative (%lf). It has been set to 0",
                iech+1,coeff);
      db->setArray(iech,iptr,0.);
    }
    else
      total += coeff;
  }

  /* Rescale */

  for (int iech=0; iech<get_NECH(db); iech++)
  {
    if (! db->isActive(iech)) continue;
    if (FFFF(db->getVariable(iech,0))) continue;
    db->updArray(iech,iptr,3,total);
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
static int st_declustering_1(Db     *db,
                             int     iptr,
                             double *radius)
{
  int error;
  double *vect,dist,value,total;

  /* Initializations */

  error = 1;
  vect  = (double *) NULL;

  /* Core allocation */

  vect = db_sample_alloc(db,LOC_X);
  if (vect == (double *) NULL) goto label_end;

  /* Loop on the target sample */
  
  for (int iech=0; iech<get_NECH(db); iech++)
  {
    if (! db->isActive(iech)) continue;
    if (FFFF(db->getVariable(iech,0))) continue;
    
    /* Loop on the second sample */
    
    for (int jech=0; jech<get_NECH(db); jech++)
    {
      if (! db->isActive(jech)) continue;
      value = db->getVariable(iech,0);
      if (FFFF(value)) continue;
      (void) distance_intra(db,iech,jech,vect);
      
      /* Normalize the distance */
      
      dist = 0.;
      for (int idim=0; idim<db->getNDim(); idim++)
      {
        vect[idim] /= radius[idim];
        dist += vect[idim] * vect[idim];
      }
      if (dist > 1) continue;
      db->updArray(iech,iptr,0,1);
    }
  }

  /* Normalization step */

  total = 0.;
  for (int iech=0; iech<get_NECH(db); iech++)
  {
    if (! db->isActive(iech)) continue;
    if (FFFF(db->getVariable(iech,0))) continue;
    value = 1. / get_ARRAY(db,iech,iptr);
    total += value;
  }
  for (int iech=0; iech<get_NECH(db); iech++)
  {
    if (! db->isActive(iech)) continue;
    if (FFFF(db->getVariable(iech,0))) continue;
    value = 1. / get_ARRAY(db,iech,iptr) / total;
    db->setArray(iech,iptr,value);
  }

  /* Set the error return code */

  error = 0;

label_end:
  vect = db_sample_free(vect);
  return(error);
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
static int st_declustering_2(Db     *db,
                             int     iptr,
                             Model  *model,
                             int     verbose)
{
  int    error,ndim,status,nech,nred,neq,ecr,nvar;
  Neigh *neigh;

  /* Initializations */

  error = 1;
  ndim  = db->getNDim();
  neigh = neigh_init_unique(ndim);
  nvar  = model->getVariableNumber();
  st_global_init(db,db);
  FLAG_EST  = 0;
  FLAG_STD  = 0;
  FLAG_VARZ = 0;
  FLAG_WGT  = 1;
  if (st_check_environment(1,1,model,neigh)) goto label_end;

  /* Pre-calculations */

  if (neigh_start(DBIN,neigh)) goto label_end;
  if (st_model_manage(1,model)) goto label_end;
  if (st_krige_manage(1,nvar,model,neigh)) goto label_end;
  if (krige_koption_manage(1,1,KOPTION_DRIFT,1,VectorInt())) goto label_end;

  /* Prepare the Neighborhood */

  (void) st_neigh(neigh,&status,&nech);
  if (status) goto label_end;

  /* Establish the L.H.S. */

  st_prepar(model,neigh,nech,&status,&nred,&neq);
  if (status) goto label_end;

  /* Loop on the targets to be processed */

  status = 0;
  IECH_OUT = 0;
  st_rhs(model,nech,neq,nvar,NULL,&status);
  if (status) goto label_end;
  st_rhs_iso2hetero(neq,1);
  if (debug_query("kriging"))
    krige_rhs_print(1,nech,neq,nred,flag,rhs);

  /* Derive the kriging weights */
  matrix_product(nred,nred,1,lhs,rhs,wgt);
  if (debug_query("kriging"))
    krige_wgt_print(status,1,1,model->getDriftEquationNumber(),nech,nred,-1,flag,wgt);

  /* Store the weights */

  ecr = 0;
  for (int iech=0; iech<get_NECH(db); iech++)
  {
    if (! db->isActive(iech)) continue;
    if (FFFF(db->getVariable(iech,0))) continue;
    db->setArray(rank[ecr],iptr,wgt[ecr]);
    ecr++;
  }
    
  /* Truncate the negative weights */

  st_declustering_truncate(db,iptr,verbose);

  /* Set the error return flag */

  error = 0;

label_end:
  debug_index(0);
  (void) st_model_manage(-1,model);
  (void) st_krige_manage(-1,nvar,model,neigh);
  (void) krige_koption_manage(-1,1,KOPTION_DRIFT,1,VectorInt());
  neigh_stop();
  neigh = neigh_free(neigh);
  return(error);
}

/****************************************************************************/
/*!
**  Perform the Declustering task as the sum of the weight 
**  for Kriging the Cells of a grid
**
** \return  Error return code
**
** \param[in]  db        input Db structure
** \param[in]  dbgrid    output Db structure
** \param[in]  iptr      Rank of the declustering weight
** \param[in]  model     Model structure
** \param[in]  neigh     Neigh structure
** \param[in]  ndisc     Array of discretization counts
** \param[in]  verbose   Verbose option
**
*****************************************************************************/
static int st_declustering_3(Db     *db,
                             Db     *dbgrid,
                             int     iptr,
                             Model  *model,
                             Neigh  *neigh,
                             VectorInt ndisc,
                             int     verbose)
{
  int    error,status,nech,nred,neq,nvar,flag_new_nbgh,ecr;
  double ldum;

  /* Initializations */

  error = 1;
  nvar  = 0;
  st_global_init(db,dbgrid);
  FLAG_EST  = 0;
  FLAG_STD  = 0;
  FLAG_VARZ = 0;
  FLAG_WGT  = 1;
  if (st_check_environment(1,1,model,neigh)) goto label_end;
  nvar = model->getVariableNumber();

  /* Pre-calculations */

  if (neigh_start(DBIN,neigh)) goto label_end;
  if (st_model_manage(1,model)) goto label_end;
  if (st_krige_manage(1,nvar,model,neigh)) goto label_end;
  if (krige_koption_manage(1,1,KOPTION_BLOCK,1,ndisc)) goto label_end;

  /* Loop on the grid cells */

  status = 0;
  for (IECH_OUT=0; IECH_OUT<get_NECH(dbgrid); IECH_OUT++)
  {
    if (! DBOUT->isActive(IECH_OUT)) continue;

    /* Select the Neighborhood */
    
    flag_new_nbgh = st_neigh(neigh,&status,&nech);
    if (status) continue;

    /* Establish the kriging L.H.S. */

    if (flag_new_nbgh || flag_continuous_kriging(neigh) || debug_force())
    {
      st_prepar(model,neigh,nech,&status,&nred,&neq);
      if (status) continue;
      st_data_dual(model,NULL,nech,nred,&ldum);
    }

    /* Establish the kriging R.H.S. */

    st_rhs(model,nech,neq,nvar,NULL,&status);
    if (status) continue;
    st_rhs_iso2hetero(neq,1);

    /* Derive the kriging weights */

    matrix_product(nred,nred,1,lhs,rhs,wgt);

    /* Cumulate the weights */

    ecr = 0;
    for (int iech=0; iech<get_NECH(db); iech++)
    {
      if (! db->isActive(iech)) continue;
      if (FFFF(db->getVariable(iech,0))) continue;
      db->updArray(rank[ecr],iptr,0,wgt[ecr]);
      ecr++;
    }
  }
    
  /* Truncate the negative weights */

  st_declustering_truncate(db,iptr,verbose);

  /* Set the error return flag */

  error = 0;

label_end:
  debug_index(0);
  (void) st_model_manage(-1,model);
  (void) st_krige_manage(-1,nvar,model,neigh);
  (void) krige_koption_manage(-1,1,KOPTION_BLOCK,1,ndisc);
  neigh_stop();
  return(error);
}

/****************************************************************************/
/*!
**  Perform the Declustering task
**
** \return  Error return code
**
** \param[in]  dbin      input Db structure
** \param[in]  model     Model structure
** \param[in]  neigh     Neigh structrue
** \param[in]  dbgrid    Grid auxiliary Db structure
** \param[in]  method    Method for declustering
** \param[in]  radius    Array of neighborhood radius
** \param[in]  ndisc     Array of discretization
** \param[in]  flag_sel  1 to mask off samples with zero weight
** \param[in]  verbose   Verbose option
**
*****************************************************************************/
GEOSLIB_API int declustering_f(Db     *dbin,
                             Model  *model,
                             Neigh  *neigh,
                             Db     *dbgrid,
                             int     method,
                             double *radius,
                             VectorInt ndisc,
                             int     flag_sel,
                             int     verbose)
{
  int error,iptr,iptr_sel;
  double indic;

  /* Initializations */

  error = 1;

  /* Preliminary checks */

  if (dbin->isVariableNumberComparedTo(0)) goto label_end;

  /* Add the kriging weight as a new variable */

  iptr = dbin->addFields(1,0.);
  if (iptr < 0) goto label_end;
  
  /* Produce statistics on the target variable before declustering */

  if (verbose) st_declustering_stats(0,method,dbin,iptr);

  /* Dispatch */

  switch (method)
  {
    case 1:			/* Weight proportional to nb samples */
      if (st_declustering_1(dbin,iptr,radius)) goto label_end;
      break;
      
    case 2:			/* Weight of the Mean */
      if (st_declustering_2(dbin,iptr,model,verbose)) goto label_end;
      break;
      
    case 3:			/* Average weight of the Block Kriging */
      if (st_declustering_3(dbin,dbgrid,iptr,model,neigh,
                            ndisc,verbose)) goto label_end;
      break;
      
    default:
      messerr("Not yet implemented");
      break;
  }
  
  /* Store the selection (optional) */

  if (flag_sel)
  {
    iptr_sel = dbin->addFields(1,0.);
    if (iptr_sel < 0) goto label_end;
    for (int iech=0; iech<get_NECH(dbin); iech++)
    {
      dbin->setArray(iech,iptr_sel,0.);
      if (! dbin->isActive(iech)) continue;
      indic = (get_ARRAY(dbin,iech,iptr) > 0.);
      dbin->setArray(iech,iptr_sel,indic);
    }
  }

  /* Produce statistics on the target variable after declustering */

  if (verbose) st_declustering_stats(1,method,dbin,iptr);

  /* Set the error return code */

  error = 0;

label_end:
  return(error);
}

/****************************************************************************/
/*!
**  Establish the covariance matrix between two Dbs
**
** \return  Covariance matrix (Dim: n1 * n2)
**
** \param[in]  title       Title of the optional printout
** \param[in]  db1         First Db structure
** \param[in]  test_def1   1 if the first variable (LOC_Z) must be checked
** \param[in]  db2         Second Db structure
** \param[in]  test_def2   1 if the second variable (LOC_Z) must be checked
** \param[in]  model       Model structure
**
** \remarks The returned argument must be freed by the calling function
**
*****************************************************************************/
static double *st_calcul_covmat(const char *title,
                                Db *db1,
                                int test_def1,
                                Db *db2,
                                int test_def2,
                                Model *model)
{
  int     n1,n2,i1,i2;
  double *covgen;
  CovCalcMode mode;

  /* Initializations */

  n1   = (test_def1) ? db1->getActiveAndDefinedNumber(0) : db1->getActiveSampleNumber();
  n2   = (test_def2) ? db2->getActiveAndDefinedNumber(0) : db2->getActiveSampleNumber();

  /* Core allocation */

  covgen = (double *) mem_alloc(sizeof(double) * n1 * n2,0);
  if (covgen == (double *) NULL) return(covgen);
  
  for (int ii1=i1=0; ii1<get_NECH(db1); ii1++)
  {
    if (test_def1)
    {
      if (! db1->isActiveAndDefined(ii1,0)) continue;
    }
    else
    {
      if (! db1->isActive(ii1)) continue;
    }
    
    for (int ii2=i2=0; ii2<get_NECH(db2); ii2++)
    {
      if (test_def2)
      {
        if (! db2->isActiveAndDefined(ii2,0)) continue;
      }
      else
      {
        if (! db2->isActive(ii2)) continue;
      }

      for (int idim=0; idim<db1->getNDim(); idim++)
        d1[idim] = (get_IDIM(db1,ii1,idim) - get_IDIM(db2,ii2,idim));

      model_calcul_cov(model,mode,1,1.,d1,&COVGEN(i1,i2));
      i2++;
    }
    i1++;
  }

  /* Optional printout */

  if (INH_FLAG_VERBOSE)
    print_matrix(title,INH_FLAG_LIMIT,1,n2,n1,NULL,covgen);

  return(covgen);
}

/****************************************************************************/
/*!
**  Establish the drift matrix for a given Db
**
** \return  Drift matrix (Dim: n1 * nbfl)
**
** \param[in]  title       Title of the optionla printout
** \param[in]  db1         First Db structure
** \param[in]  test_def1   1 if the first variable (LOC_Z) must be checked
** \param[in]  model       Model structure
**
** \remarks The returned argument must be freed by the calling function
**
*****************************************************************************/
static double *st_calcul_drfmat(const char *title,
                                Db *db1,
                                int test_def1,
                                Model *model)
{
  int     i1,n1,nbfl;
  double *drftab;

  /* Initializations */

  n1   = (test_def1) ? db1->getActiveAndDefinedNumber(0) : db1->getActiveSampleNumber();
  nbfl = model->getDriftNumber();

  /* Core allocation */

  drftab = (double *) mem_alloc(sizeof(double) * n1 * nbfl,0);
  if (drftab == (double *) NULL) return(drftab);

  /* Loop on the samples */

  i1 = 0;
  for (int ii1=0; ii1<get_NECH(db1); ii1++)
  {
    if (test_def1)
    {
      if (! db1->isActiveAndDefined(ii1,0)) continue;
    }
    else
    {
      if (! db1->isActive(ii1)) continue;
    }
    
    model_calcul_drift(model,MEMBER_LHS,db1,ii1,&drftab[i1 * nbfl]);
    i1++;
  }

  /* Optional printout */

  if (INH_FLAG_VERBOSE)
    print_matrix(title,INH_FLAG_LIMIT,1,nbfl,n1,NULL,drftab);

  return(drftab);
}


/****************************************************************************/
/*!
**  Establish the distance matrix between two Dbs
**
** \return  Covariance matrix
**
** \param[in]  title       Title of the optional printout
** \param[in]  db1         First Db structure
** \param[in]  test_def1   1 if the first variable (LOC_Z) must be checked
** \param[in]  db2         Second Db structure (sources)
** \param[in]  test_def2   1 if the second variable (LOC_Z) must be checked
** \param[in]  power       Power of the Distance decay
**
** \remarks The returned argument must be freed by the calling function
**
*****************************************************************************/
static double *st_calcul_distmat(const char *title,
                                 Db *db1,
                                 int test_def1,
                                 Db *db2,
                                 int test_def2,
                                 double power)
{
  int     n1,ns,i1,is,ndim;
  double *distgen,dist;

  /* Initializations */

  n1   = (test_def1) ? db1->getActiveAndDefinedNumber(0) : db1->getActiveSampleNumber();
  ns   = (test_def2) ? db2->getActiveAndDefinedNumber(0) : db2->getActiveSampleNumber();
  ndim = db1->getNDim();

  /* Core allocation */

  distgen = (double *) mem_alloc(sizeof(double) * n1 * ns,0);
  if (distgen == (double *) NULL) return(distgen);
  
  for (int ii1=i1=0; ii1<get_NECH(db1); ii1++)
  {
    if (test_def1)
    {
      if (! db1->isActiveAndDefined(ii1,0)) continue;
    }
    else
    {
      if (! db1->isActive(ii1)) continue;
    }
    
    for (int iis=is=0; iis<get_NECH(db2); iis++)
    {
      if (test_def2)
      {
        if (! db2->isActiveAndDefined(iis,0)) continue;
      }
      else
      {
        if (! db2->isActive(iis)) continue;
      }

      dist = 0.;
      for (int idim=0; idim<ndim; idim++)
      {
        d1[idim] = (get_IDIM(db1,ii1,idim) - get_IDIM(db2,iis,idim));
        dist += d1[idim] * d1[idim];
      }

      DISTGEN(i1,is) = 1. / pow(dist,power/2.);
      is++;
    }
    i1++;
  }

  /* Optional printout */

  if (INH_FLAG_VERBOSE)
    print_matrix(title,INH_FLAG_LIMIT,1,ns,n1,NULL,distgen);

  return(distgen);
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
static double *st_calcul_product(const char *title,
                                 int n1,
                                 int ns,
                                 double *covss,
                                 double *distgen)
{
  double *prodgen;
  
  prodgen  = (double *) mem_alloc(sizeof(double) * n1 * ns,0);
  if (prodgen  == (double *) NULL) return(prodgen);

  for (int i1=0; i1<n1; i1++)
    for (int is=0; is<ns; is++)
    {
      PRODGEN(i1,is) = 0.;
      for (int js=0; js<ns; js++) 
        PRODGEN(i1,is) += COVSS(is,js) * DISTGEN(i1,js);
    }

  /* Optional printout */

  if (INH_FLAG_VERBOSE)
    print_matrix(title,INH_FLAG_LIMIT,1,ns,n1,NULL,prodgen);

  return(prodgen);
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
** \param[in]  covss       Covariance matrix between Sources
** \param[in]  distps      Distance matrix between Data and Sources
** \param[in]  prodps      Product of DistPS by CovSS
**
*****************************************************************************/
static double *st_inhomogeneous_covpp(Db     *dbdat,
                                      Db     *dbsrc,
                                      Model  *model_dat,
                                      double *covss,
                                      double *distps,
                                      double *prodps)
{
  double *covpp;
  int     np,ns,error;

  /* Initializations */

  error = 1;
  covpp = (double *) NULL;
  
  np = dbdat->getActiveAndDefinedNumber(0);
  ns = dbsrc->getActiveSampleNumber();

  /* Covariance matrix between Mesures */

  covpp = st_calcul_covmat("Covariance P-P",dbdat,1,dbdat,1,model_dat);
  if (covpp == (double *) NULL) goto label_end;

  /* Calculate the LHS matrix */

  for (int ip=0; ip<np; ip++)
    for (int jp=ip; jp<np; jp++)
      for (int is=0; is<ns; is++)
      {
        COVPP(ip,jp) += DISTPS(ip,is) * PRODPS(jp,is);
        if (jp > ip) COVPP(jp,ip) = COVPP(ip,jp);
      }
  
  /* Set the error return code */

  error = 0;

label_end:
  if (error) covpp = (double *) mem_free((char *) covpp);
  return(covpp);
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
static double *st_inhomogeneous_covgp(Db     *dbdat,
                                      Db     *dbsrc,
                                      Db     *dbout,
                                      int     flag_source,
                                      Model  *model_dat,
                                      double *distps,
                                      double *prodps,
                                      double *prodgs)
{
  double *covgp;
  int     np,ns,ng,error;

  /* Initializations */

  error = 1;
  covgp = (double *) NULL;
  
  np = dbdat->getActiveAndDefinedNumber(0);
  ns = dbsrc->getActiveSampleNumber();
  ng = dbout->getActiveSampleNumber();

  /* Covariance matrix between Mesures and Target */

  covgp = st_calcul_covmat("Covariance G-P",dbout,0,dbdat,1,model_dat);
  if (covgp == (double *) NULL) goto label_end;

  /* Add the contribution of the source */

  if (! flag_source)
  {
    for (int ig=0; ig<ng; ig++)
      for (int ip=0; ip<np; ip++)
        for (int is=0; is<ns; is++)
          COVGP(ig,ip) += DISTPS(ip,is) * PRODGS(ig,is);
  }
  else
  {
    for (int ig=0; ig<ng; ig++)
      for (int ip=0; ip<np; ip++)
        COVGP(ig,ip) = PRODPS(ip,ig);
  }
  
  /* Set the error return code */

  error = 0;

label_end:
  if (error) covgp = (double *) mem_free((char *) covgp);
  return(covgp);
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
static double *st_inhomogeneous_covgg(Db     *dbsrc,
                                      Db     *dbout,
                                      int     flag_source,
                                      Model  *model_dat,
                                      double *distgs,
                                      double *prodgs)
{
  int     ng,ns,error;
  double *covgg,c00;
  CovCalcMode mode;

  /* Initializations */

  error = 1;
  covgg = (double *) NULL;
  
  ns = dbsrc->getActiveSampleNumber();
  ng = dbout->getActiveSampleNumber();

  /* Core allocation */

  covgg = (double *) mem_alloc(sizeof(double) * ng,0);
  if (covgg == (double *) NULL) goto label_end;
  
  /* Calculate the variance term (for a zero-distance) */
  
  model_calcul_cov(model_dat,mode,1,1.,VectorDouble(),&c00);

  /* Calculate the variance vector */
  
  if (! flag_source)
  {
    for (int ig=0; ig<ng; ig++)
    {
      covgg[ig] = c00;
      for (int is=0; is<ns; is++)
        covgg[ig] += DISTGS(ig,is) * PRODGS(ig,is);
    }
  }
  else
  {
    for (int ig=0; ig<ng; ig++)
      covgg[ig] = c00;
  }
  
  /* Set the error return code */

  error = 0;

label_end:
  if (error) covgg = (double *) mem_free((char *) covgg);
  return(covgg);
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
static int st_drift_prepar(int     np,
                           int     nbfl,
                           double *covpp,
                           double *drftab,
                           double **yloc,
                           double **zloc)
{
  double *ymat,*zmat,value;
  int     error,ecr;
  
  /* Initialization */

  error = 1;
  ymat = zmat = (double *) NULL;

  /* First returned array */
  
  ymat = (double *) mem_alloc(sizeof(double) * nbfl * np,0);
  if (ymat == (double *) NULL) goto label_end;
  
  ecr = 0;
  for (int il=0; il<nbfl; il++)
    for (int ip=0; ip<np; ip++)
    {
      value = 0.;
      for (int jp=0; jp<np; jp++)
        value += COVPP(ip,jp) * DRFTAB(jp,il);
      ymat[ecr++] = value;
    }

  /* Second retrned array */

  zmat = (double *) mem_alloc(sizeof(double) * nbfl * nbfl,0);
  if (zmat == (double *) NULL) goto label_end;
  
  ecr = 0;
  for (int il=0; il<nbfl; il++)
    for (int jl=0; jl<nbfl; jl++)
    {
      value = 0.;
      for (int ip=0; ip<np; ip++)
        value += YMAT(ip,il) * DRFTAB(ip,jl);
      zmat[ecr++] = value;
    }      
  
  /* Invert 'zmat' */

  if (matrix_invert(zmat,nbfl,-1)) goto label_end;
  
  /* Set the error return code */

  error = 0;

label_end:
  if (error)
  {
    ymat = (double *) mem_free((char *) ymat);
    zmat = (double *) mem_free((char *) zmat);
  }
  *yloc = ymat;
  *zloc = zmat;
  return(error);
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
** \param[in]  lambda      Vector of weights
**
** \param[out] lambda      Vector of weights
** \param[out] mu          Vector of Lagrange parameters
**
*****************************************************************************/
static void st_drift_update(int     np,
                            int     nbfl,
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

  for (int il=0; il<nbfl; il++)
  {
    value = 0.;
    for (int ip=0; ip<np; ip++)
      value = YMAT(ip,il) * covgp[ip] - driftg[il];
    maux[il] = value;
  }
  matrix_product(nbfl,nbfl,1,zmat,maux,mu);
  
  /* Update the vector of kriging weights */

  for (int ip=0; ip<np; ip++)
    for (int il=0; il<nbfl; il++)
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
GEOSLIB_API int inhomogeneous_kriging(Db     *dbdat,
                                      Db     *dbsrc,
                                      Db     *dbout,
                                      double  power,
                                      int     flag_source,
                                      Model  *model_dat,
                                      Model  *model_src)
{
  int     error,np,ip,ns,ng,nvar,neq,nred,nfeq,nbfl;
  double *covss,*distps,*distgs,*covpp,*covgp,*covgg,*prodps,*prodgs;
  double *data,*lambda,*driftp,*driftg,*ymat,*zmat,*mu,*maux,*rhs;
  double  estim,stdev,auxval;
  Neigh  *neigh;

  /* Preliminary checks */

  error = nvar = 1;
  neigh = neigh_init_unique(dbdat->getNDim());
  st_global_init(dbdat,dbout);
  FLAG_EST  = 1;
  FLAG_STD  = 1;
  distps = distgs = prodgs = prodps = (double *) NULL;
  covss  = covpp  = covgp  = covgg  = (double *) NULL;
  lambda = data   = driftp = driftg = (double *) NULL;
  ymat   = zmat   = mu     = maux   = (double *) NULL;
  if (st_check_environment(1,1,model_dat,NULL)) goto label_end;

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

  if (FLAG_EST)
  {
    IPTR_EST = dbout->addFields(nvar,0.);
    if (IPTR_EST < 0) goto label_end;
  }
  if (FLAG_STD)
  {
    IPTR_STD = dbout->addFields(nvar,0.);
    if (IPTR_STD < 0) goto label_end;
  }
  nred = neq = np = dbdat->getActiveAndDefinedNumber(0);
  nfeq = 0;
  ns   = dbsrc->getActiveSampleNumber();
  ng   = dbout->getActiveSampleNumber();
  nbfl = model_dat->getDriftNumber();

  /* Core allocation */

  lambda = (double *) mem_alloc(sizeof(double) * np,0);
  if (lambda == (double *) NULL) goto label_end;
  data   = (double *) mem_alloc(sizeof(double) * np,0);
  if (data   == (double *) NULL) goto label_end;
  
  /* Pre-calculations */

  if (st_model_manage(1,model_dat)) goto label_end;
  if (st_krige_manage(1,nvar,model_dat,neigh)) goto label_end;
  if (krige_koption_manage(1,1,KOPTION_PONCTUAL,1,VectorInt())) goto label_end;

  /* Constitute the Data vector */

  for (int iip=ip=0; iip<get_NECH(dbdat); iip++)
  {
    if (! dbdat->isActiveAndDefined(iip,0)) continue;
    data[ip] = dbdat->getVariable(iip,0);
    rank[ip] = iip;
    ip++;
  }
  
  /* Establish the covariance matrix between Sources */

  covss = st_calcul_covmat("Covarance S_S",dbsrc,0,dbsrc,0,model_src);
  if (covss == (double *) NULL) goto label_end;

  /* Establish the distance matrix between Data and Sources */
  
  distps = st_calcul_distmat("Distance P-S",dbdat,1,dbsrc,0,power);
  if (distps == (double *) NULL) goto label_end;

  /* Establish the distance matrix between Target and Sources */

  if (! flag_source)
  {
    distgs = st_calcul_distmat("Distance G-S",dbout,0,dbsrc,0,power);
    if (distgs == (double *) NULL) goto label_end;
  }
  
  /* Establish the Data-Source Product matrix */

  prodps = st_calcul_product("Convolve P-S",np,ns,covss,distps);
  if (prodps == (double *) NULL) goto label_end;
  
  /* Establish the complete kriging matrix */

  covpp = st_inhomogeneous_covpp(dbdat,dbsrc,model_dat,covss,distps,prodps);
  if (covpp == (double *) NULL) goto label_end;
  if (debug_query("kriging") || is_debug_reference_defined())
    krige_lhs_print(np,neq,nred,NULL,covpp);

  /* Invert the Kriging Matrix */
  
  if (matrix_invert(covpp,np,-1)) goto label_end;
  
  /* Establish the drift at Data */
     
  if (nbfl > 0)
  {
    mu   = (double *) mem_alloc(sizeof(double) * nbfl,0);
    if (mu   == (double *) NULL) goto label_end;
    maux = (double *) mem_alloc(sizeof(double) * nbfl,0);
    if (maux == (double *) NULL) goto label_end;
    
    driftp = st_calcul_drfmat("Drift P",dbdat,1,model_dat);
    if (driftp == (double *) NULL) goto label_end;

    /* Prepare auxiliary arrays */

    if (st_drift_prepar(np,nbfl,covpp,driftp,&ymat,&zmat)) goto label_end;
  }

  /* Establish the Target-Source Product matrix */

  if (! flag_source)
  {
    prodgs = st_calcul_product("Convolve G-S",ng,ns,covss,distgs);
    if (prodgs == (double *) NULL) goto label_end;
  }
  
  /* Establish the COVGP */
  
  covgp = st_inhomogeneous_covgp(dbdat,dbsrc,dbout,flag_source,model_dat,
                                 distps,prodps,prodgs);
  if (covgp == (double *) NULL) goto label_end;
  
  /* Establish the drift at Target */
     
  if (nbfl > 0)
  {
    driftg = (double *) mem_alloc(sizeof(double) * nbfl,0);
    if (driftg == (double *) NULL) goto label_end;
  }

  /* Establish the variance at targets */

  covgg = st_inhomogeneous_covgg(dbsrc,dbout,flag_source,model_dat,
                                 distgs,prodgs);
  if (covgg == (double *) NULL) goto label_end;
  
  /* Loop on the targets to be processed */

  for (IECH_OUT=0; IECH_OUT<get_NECH(DBOUT); IECH_OUT++)
  {
    mes_process("Kriging sample",get_NECH(DBOUT),IECH_OUT);
    debug_index(IECH_OUT+1);
    if (! dbout->isActive(IECH_OUT)) continue;
    if (debug_query("kriging") ||
        debug_query("nbgh")    ||
        debug_query("results"))
    {
      mestitle(1,"Target location");
      db_sample_print(dbout,IECH_OUT,1,0,0);
    }
    rhs = &COVGP(IECH_OUT,0);

    /* Optional printout of the R.H.S */
    
    if (debug_force()) krige_rhs_print(nvar,np,neq,nred,NULL,rhs);

    /* Fill the drift at Target point (optional) */
    
    if (driftp != (double *) NULL)
      model_calcul_drift(model_dat,MEMBER_LHS,dbout,IECH_OUT,driftg);
    
    /* Calculate the Kriging weights */

    matrix_product(np,np,1,covpp,rhs,lambda);
    if (debug_force())
      krige_wgt_print(0,nvar,nvar,nfeq,np,nred,-1,NULL,lambda);
      
    /* Update vector of weights in presence of drift */

    if (nbfl > 0)
    {

      /* Evaluate the drift at Target */
      
      model_calcul_drift(model_dat,MEMBER_LHS,dbout,IECH_OUT,driftg);
      
      /* Update the kriging weights */
      
      st_drift_update(np,nbfl,rhs,driftg,ymat,zmat,maux,lambda,mu);
    }
    
    /* Perform the estimation */

    matrix_product(1,np,1,data,lambda,&estim);
    matrix_product(1,np,1,rhs,lambda,&stdev);

    /* Update the variance in presence of drift */

    if (nbfl > 0)
    {
      matrix_product(1,nbfl,1,mu,maux,&auxval);
      stdev += auxval;
    }
    
    /* Update the variance calculation */
    
    VAR0(0,0) = covgg[IECH_OUT];
    stdev = covgg[IECH_OUT] - stdev;
    stdev = (stdev > 0) ? sqrt(stdev) : 0.;

    /* Store the result */

    dbout->setArray(IECH_OUT,IPTR_EST,estim);
    dbout->setArray(IECH_OUT,IPTR_STD,stdev);

    /* Optional printout */

    if (debug_query("kriging") || debug_force())
      st_result_kriging_print(0,nvar,0);
  }

  /* Set the error return flag */

  error = 0;

label_end:
  debug_index(0);
  covss  = (double *) mem_free((char *) covss);
  distps = (double *) mem_free((char *) distps);
  distgs = (double *) mem_free((char *) distgs);
  prodps = (double *) mem_free((char *) prodps);
  prodgs = (double *) mem_free((char *) prodgs);
  driftp = (double *) mem_free((char *) driftp);
  driftg = (double *) mem_free((char *) driftg);
  covpp  = (double *) mem_free((char *) covpp);
  covgp  = (double *) mem_free((char *) covgp);
  covgg  = (double *) mem_free((char *) covgg);
  driftp = (double *) mem_free((char *) driftp);
  driftg = (double *) mem_free((char *) driftg);
  ymat   = (double *) mem_free((char *) ymat);
  zmat   = (double *) mem_free((char *) zmat);
  maux   = (double *) mem_free((char *) maux);
  mu     = (double *) mem_free((char *) mu);
  data   = (double *) mem_free((char *) data);
  lambda = (double *) mem_free((char *) lambda);
  (void) st_model_manage(-1,model_dat);
  (void) st_krige_manage(-1,1,model_dat,neigh);
  (void) krige_koption_manage(-1,1,KOPTION_PONCTUAL,1,VectorInt());
  neigh_stop();
  neigh = neigh_free(neigh);
  return(error);
}
