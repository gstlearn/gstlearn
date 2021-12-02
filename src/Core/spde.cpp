/******************************************************************************/
/* COPYRIGHT ARMINES, ALL RIGHTS RESERVED                                     */
/*                                                                            */
/* THE CONTENT OF THIS WORK CONTAINS CONFIDENTIAL AND PROPRIETARY             */
/* INFORMATION OF ARMINES. ANY DUPLICATION, MODIFICATION,                     */
/* DISTRIBUTION, OR DISCLOSURE IN ANY FORM, IN WHOLE, OR IN PART, IS STRICTLY */
/* PROHIBITED WITHOUT THE PRIOR EXPRESS WRITT%EN PERMISSION OF ARMINES        */
/*                                                                            */
/* TAG_SOURCE_CG                                                              */
/******************************************************************************/
//#include "geoslib_e.h"
#include "geoslib_enum.h"
#include "geoslib_old_f.h"
#include "geoslib_f.h"
#include "Matrix/MatrixFactory.hpp"
#include "Matrix/MatrixSquareGeneral.hpp"
#include "Model/NoStatArray.hpp"
#include "Mesh/MeshEStandard.hpp"
#include "Covariances/CovAniso.hpp"
#include "Covariances/ACovAnisoList.hpp"
#include "Basic/AException.hpp"
#include "Basic/Utilities.hpp"
#include "Basic/Law.hpp"
#include "Basic/MathFunc.hpp"
#include "Basic/File.hpp"
#include "Basic/String.hpp"
#include "Db/Db.hpp"
#include "Model/Model.hpp"
#include "Mesh/tetgen.h"
#include "csparse_f.h"
#include "csparse_d.h"

#include <math.h>
#include <string.h>

/* Global symbols for SPDE */

#define NBLIN_TERMS 10
#define SPDE_MAX_NGRF 2

/*! \cond */
#define VT_NONE      -1
#define VT_IDLE       0

#define VT_FREE       1
#define VT_GIBBS      2
#define VT_HARD       4

#define VT_INPUT      8
#define VT_OUTPUT    16
#define VT_OTHER     32

#define CASE_MATRICES 0 
#define CASE_KRIGING  1
#define CASE_SIMULATE 2

#define OLD_MESHES(imesh,icorn) (s_mesh->meshes[(imesh)*s_mesh->ncorner+(icorn)]-1)
#define POINTS(ip,idim)      (s_mesh->points[(ip) * s_mesh->ndim + (idim)])
#define IADH(ndim,i,j)       (ndim * (i) + (j))
#define TEMP(ndim,i,j)       (temp[IADH(ndim,i,j)])
#define Z(ivar,nech,iech)    (z[(ivar) * nech + (iech)])
#define M(j,i)               (m[(i) * ndimp + (j)])
#define TP(j,i)              (tp[(i) * ndimp + (j)]) 
#define VEC1(ip,idim)        (vec1[(ip) * ndim + (idim)])
#define MAT(i,j)             (mat[(i) * ncorner + (j)])
#define MATU(i,j)            (matu[(i) * ncorner + (j)])
#define COTES(ip,i)          (cotes[(ip) * ncorner + (i)])
#define CORVEC(idim,ip)      (coor[(idim) * nvertex + (ip)])
#define TBLIN(ib,ip)         (tblin[nvertex * (ib) + (ip)])
#define CONTAIN(imesh,idim,i)(contain[(i) + 2 * ((idim) + ndim * (imesh))])
#define MATGRF(igrf)         (&S_ENV.SS_ENV[igrf])
#define LOCAL(ivr,jvr)       (local[(ivr) * nvr + (jvr)])
#define LOCAL0(ivar,jvar)    (local0[(ivar) * nvar + (jvar)])
#define ADM(icov,ivar,icur)  ((icur) + ncur * ((ivar) + nvar * (icov)))
#define RHS(icov,ivar,icur)  (rhs[ADM(icov,ivar,icur)])
#define XCUR(icov,ivar,icur) (xcur[ADM(icov,ivar,icur)])
#define TAB(icov,ivar,icur)  (tab[ADM(icov,ivar,icur)])
#define DATA(ivar,idata)     (data[(ivar) * ndata + (idata)])
#define GWORK(ilayer,iech)   (gwork[(ilayer) * ngrid   + (iech)])
#define YVERT(ilayer,iech)   (yvert[(ilayer) * nvertex + (iech)])
#define YDAT(ilayer,iech)    (ydat [(ilayer) * nech    + (iech)])
#define YMEAN(ilayer,iech)   (ymean[(ilayer) * nech    + (iech)])
#define DCOEF(ilayer)        (m2denv->dcoef[ilayer])

typedef struct
{
  std::vector<SPDE_Matelem> Matelems;
  int ndata; /* Number of active data */
  int *ndata1; /* Number of data per variable (icov=0) */
  int *ntarget1; /* Number of target per variable (icov=0) */
  Model *model; /* Pointer to the Model */
  double *Csill; /* Array of LU of sill matrices */
  cs **Bnugget; /* Sparse matrices for nugget effect (nvs2) */
  cs **BheteroD; /* Sparse matrices for heterotopy (nvar)*/
  cs **BheteroT; /* Sparse matrices for heterotopy (nvar)*/
} SPDE_SS_Environ;

typedef struct
{
  int ndim; /* Space Dimension */
  int nvar; /* Number of variables */
  int ngrfs; /* Number of GRFs */
  SPDE_SS_Environ SS_ENV[SPDE_MAX_NGRF];
} SPDE_Environ;

typedef struct
{
  int flag_dbin; /* Presence of an input Db */
  int flag_dbout; /* Presence of an output Db */
  int flag_mesh_dbin; /* Input points participate to meshing */
  int flag_mesh_dbout; /* Output points participate to meshing */
  int flag_est; /* Perform Estimation */
  int flag_std; /* Perform Standard deviation */
  int flag_case; /* Perform: matrices(0), est(1) or simu(2) */
  int flag_gibbs; /* Perform Gibbs sampling */
  int flag_modif; /* Post-processing simulations */
  int flag_onechol; /* Perform Simu & Kriging with same Chol */
  int flag_filnug; /* Filtering the Nugget Effect */
  int flag_mgrid; /* Use the Multigrid option */
  int flag_several; /* Perform Kriging in iterative mode */
  int simu_chol; /* Use Cholesky simulation */
  int simu_cheb; /* Use Chebychev simulation */
  int flag_Q; /* Build Q */
  int flag_Qchol; /* Perform Cholesky on global Q */
} SPDE_Decision;

typedef struct
{
  int flag_ed;
  int iatt_fd;
  int iatt_fg;
  double zmean;
  double zstdv;
  double zeps;
  double zmini;
  double zmaxi;
  double dmini;
  double dmaxi;
  double ystdv;
  double *dcoef;
} M2D_Environ;

typedef struct
{
  int flag_sphere;
  double sqdeth;
  double correc;
  double R;
  VectorDouble blin;
  VectorDouble hh;
  VectorDouble vv;
  VectorDouble srot;
} SPDE_Calcul;

static void (*SIMU_FUNC_TRANSF)(Db*, int, int, int) = NULL;
static void (*SIMU_FUNC_UPDATE)(Db*, int, int, int) = NULL;
static void (*SIMU_FUNC_SCALE)(Db*, int, int) = NULL;

/*! \endcond */
static int DEBUG = 0;
static int VERBOSE = 0;
static int FLAG_KEYPAIR = 0;
static double FACDIM[] = { 0., 1., 2., 6. };
static int SPDE_CURRENT_IGRF = 0;
static int SPDE_CURRENT_ICOV = 0;
static Db *MEM_DBIN, *MEM_DBOUT;
static SPDE_Mesh *S_EXTERNAL_MESH[3] = { NULL, NULL, NULL };
static cs *S_EXTERNAL_Q[3] = { NULL, NULL, NULL };
static cs *S_EXTERNAL_A[3] = { NULL, NULL, NULL };
static SPDE_Environ S_ENV;
static SPDE_Decision S_DECIDE;
static char NAME[100];
static char string_encode[100];
static SPDE_Calcul Calcul;

/****************************************************************************/
/*!
 **  Returns the index of a pair of variable ranks within the triangle
 **
 ** \return Absolute index
 **
 ** \param[in] ivar    Rank of the first variable
 ** \param[in] jvar    Rank of the second variable
 **
 ** \remarks The calling function does not have to bother of the relative
 ** \remarks order between 'ivar' and 'jvar'
 **
 *****************************************************************************/
static int st_get_rank(int ivar, int jvar)
{
  if (jvar > ivar)
    return (jvar * (jvar + 1) / 2 + ivar);
  else
    return (ivar * (ivar + 1) / 2 + jvar);
}

/****************************************************************************/
/*!
 **  Print the contents of one SP_Mat structure
 **
 ** \param[in] icov    Rank of the covariance
 **
 *****************************************************************************/
static void st_matelem_print(int icov)

{
  static const char *NOK[] = { "OFF", "ON" };

  const SPDE_Matelem &Matelem = spde_get_current_matelem(icov);

  mestitle(1, "Contents of Matelem structure #%d", icov + 1);
  message("S is defined:      %s\n", NOK[Matelem.S != NULL]);
  message("Aproj is defined:  %s\n", NOK[Matelem.Aproj != NULL]);
  message("QC is defined:     %s\n", NOK[Matelem.QC != NULL]);
  message("QCov are defined:  %s\n", NOK[Matelem.QCov != NULL]);
  message("Lambda is defined: %s\n", NOK[!Matelem.Lambda.empty()]);
  message("qsimu is defined:  %s\n", NOK[Matelem.qsimu != NULL]);
  message("mgs is defined:    %s\n", NOK[Matelem.mgs != NULL]);
  message("s_cheb is defined: %s\n", NOK[Matelem.s_cheb != NULL]);
  message("s_mesh is defined: %s\n", NOK[Matelem.s_mesh != NULL]);
}

/****************************************************************************/
/*!
 **  Add a new option item for an additional covariance structure
 **
 ** \param[in]  s_option    Pointer to SPDE_Option to be freed (if mode<0)
 ** \param[in]  triswitch   String defining the meshing characteristics
 **
 *****************************************************************************/
void spde_option_update(SPDE_Option &s_option, const String &triswitch)
{
  /* Add a new SPDE_SS_Option structure */

  int noption = static_cast<int>(s_option.options.size());

  /* Resize the array 'options' */

  s_option.options.resize(noption + 1);

  /* Initialize this new 'SPDE_SS_Option' structure */

  SPDE_SS_Option ss_option = s_option.options[noption];
  s_option.options[noption].mesh_dbin = 0;
  s_option.options[noption].mesh_dbout = 0;
  s_option.options[noption].triswitch = triswitch;
}

/****************************************************************************/
/*!
 **  Manage the SPDE_Option structure
 **
 *****************************************************************************/
SPDE_Option spde_option_alloc(void)
{
  SPDE_Option s_option;

  s_option.options = std::vector<SPDE_SS_Option>();
  return s_option;
}

/****************************************************************************/
/*!
 **  Set the rank of the current GRF of the environment
 **
 ** \param[in]  igrf      Rank of the current GRF
 **
 *****************************************************************************/
static void st_set_current_igrf(int igrf)
{
  SPDE_CURRENT_IGRF = igrf;
}

/****************************************************************************/
/*!
 **  Set the rank of the current Covariance of the environment
 **
 ** \param[in]  icov      Rank of the current Covariance
 **
 *****************************************************************************/
static void st_set_current_icov(int icov)
{
  SPDE_CURRENT_ICOV = icov;
}

/****************************************************************************/
/*!
 **  Get the rank of the current GRF of the environment
 **
 *****************************************************************************/
static int st_get_current_igrf(void)
{
  return (SPDE_CURRENT_IGRF);
}

/****************************************************************************/
/*!
 **  Get the pointer to the current SPDE_SS_Environ structure
 **
 *****************************************************************************/
static SPDE_SS_Environ* st_get_current_ssenv(void)
{
  return (MATGRF(st_get_current_igrf()));
}

/****************************************************************************/
/*!
 **  Get the rank of the current COV of the environment
 **
 *****************************************************************************/
static int st_get_current_icov(void)
{
  return (SPDE_CURRENT_ICOV);
}

/****************************************************************************/
/*!
 **  Get the current value for triswitch parameter
 **
 ** \param[in]  s_option   SPDE_Option structure
 **
 *****************************************************************************/
static String st_get_current_triswitch(SPDE_Option &s_option)
{
  int rank_cov = st_get_current_icov();
  int noption = static_cast<int>(s_option.options.size());
  int rank = MIN(rank_cov, noption - 1);
  String triswitch = s_option.options[rank].triswitch;
  return (triswitch);
}

/****************************************************************************/
/*!
 **  Get the pointer to the current SPDE_Matelem structure
 **
 ** \param[in] icov    Rank of the target Covariance (or -1)
 **
 *****************************************************************************/
SPDE_Matelem& spde_get_current_matelem(int icov)
{
  if (icov < 0)
    return (st_get_current_ssenv()->Matelems[st_get_current_icov()]);
  else
    return (st_get_current_ssenv()->Matelems[icov]);
}

/****************************************************************************/
/*!
 **  Update a string to include the rank of the current GRF and Covariance
 **
 ** \param[in]  flag_igrf  To add current GRF
 ** \param[in]  flag_icov  To add current COV
 ** \param[in]  rank       Rank of the highlight (see mestitle or -1)
 ** \param[in]  title      Input title
 **
 *****************************************************************************/
static void st_title(int flag_igrf, int flag_icov, int rank, const char *title)
{
  int flag_decor;

  (void) gslStrcpy(string_encode, " ");

  flag_decor = (flag_igrf || flag_icov);

  if (flag_decor)
  {
    (void) gslStrcpy(string_encode, "(");
    if (flag_igrf)
      (void) gslSPrintf(string_encode, "%s GRF:%d", string_encode,
                        st_get_current_igrf() + 1);
    if (flag_icov)
      (void) gslSPrintf(string_encode, "%s - COV:%d", string_encode,
                        st_get_current_icov() + 1);
    (void) gslSPrintf(string_encode, "%s ) %s", string_encode, title);
  }
  else
  {
    (void) gslSPrintf(string_encode, "%s", title);
  }

  if (rank >= 0)
    mestitle(rank, string_encode);
  else
  {
    (void) gslSPrintf(string_encode, "%s\n", string_encode);
    message(string_encode);
  }
}

/****************************************************************************/
/*!
 **  Return the global non-stationary characteristics

 **
 *****************************************************************************/
static Model* st_get_model(void)
{
  return (st_get_current_ssenv()->model);
}

/****************************************************************************/
/*!
 **  Returns the number of GRFs of the environment
 **
 *****************************************************************************/
static int st_get_number_grf(void)
{
  int ngrfs;
  ngrfs = MAX(1, S_ENV.ngrfs);
  return (ngrfs);
}

/****************************************************************************/
/*!
 **  Initialize the S_ENV Environment structure
 **
 *****************************************************************************/
static void st_environ_init(void)
{
  SPDE_SS_Environ *SS;

  S_ENV.ndim = 0;
  S_ENV.nvar = 0;
  S_ENV.ngrfs = 0;
  for (int igrf = 0; igrf < SPDE_MAX_NGRF; igrf++)
  {
    SS = &S_ENV.SS_ENV[igrf];
    SS->ndata = 0;
    SS->ndata1 = nullptr;
    SS->ntarget1 = nullptr;
    SS->model = nullptr;
    SS->Bnugget = nullptr;
    SS->BheteroD = nullptr;
    SS->BheteroT = nullptr;
  }
}

/****************************************************************************/
/*!
 **  Manage the Multigrid solving operations
 **
 ** \return  Error return code
 **
 ** \param[in]  mode       1 for allocation; -1 for deallocation
 ** \param[in]  mgs        cs_MGS to be freed (only for mode=-1)
 **
 *****************************************************************************/
static cs_MGS* st_mgs_manage(int mode, cs_MGS *mgs)
{
  int nlevels, path_type, flag_cg, ngc, nmg, ngs, type_coarse;
  double tolcg, tolnmg;

  /* Dispatch */

  if (mode > 0)
  {

    /* Initialize the cs_MGS structure */

    nlevels = (int) get_keypone("Multigrid_Number_Levels", 0);
    path_type = (int) get_keypone("Multigrid_Path_Type", 1);
    mgs = cs_multigrid_manage(NULL, 1, nlevels, path_type);
    if (mgs == (cs_MGS*) NULL) return (mgs);
    flag_cg = (int) get_keypone("Flag_CG", 1);
    type_coarse = (int) get_keypone("Multigrid_Coarse", 0);
    ngc = (int) get_keypone("Multigrid_ngc", 100);
    ngs = (int) get_keypone("Multigrid_ngs", 2);
    nmg = (int) get_keypone("Multigrid_nmg", 4);
    tolcg = get_keypone("Multigrid_tolcg", 1.e-7);
    tolnmg = get_keypone("Multigrid_tolnmg", 1.e-7);
    cs_multigrid_params(mgs, flag_cg, type_coarse, ngc, nmg, ngs, tolcg,
                        tolnmg);
  }
  else
  {

    /* Free the cs_MGS structure */

    mgs = cs_multigrid_manage(mgs, -1, 0, 0);
  }
  return (mgs);
}

/****************************************************************************/
/*!
 **  Manage the SPDE_Mesh structure
 **
 ** \return  The newly allocated SPDE_Mesh
 **
 ** \param[in]  mode       1 for allocation; -1 for deallocation
 ** \param[in]  s_mesh_old Pointer to SPDE_Mesh to be deallocated
 **
 *****************************************************************************/
SPDE_Mesh* spde_mesh_manage(int mode, SPDE_Mesh *s_mesh_old)
{
  SPDE_Mesh *s_mesh;

  /* Initializations */

  s_mesh = nullptr;

  /* Dispatch */

  if (mode > 0)
  {

    /* Allocation */

    s_mesh = (SPDE_Mesh*) mem_alloc(sizeof(SPDE_Mesh), 0);
    if (s_mesh == nullptr) return (s_mesh);
    s_mesh->ndim = 0;
    s_mesh->ncorner = 0;
    s_mesh->nmesh = 0;
    s_mesh->nvertex = 0;
    s_mesh->meshes = nullptr;
    s_mesh->points = nullptr;
    s_mesh->vercoloc = nullptr;
    s_mesh->vertype = nullptr;
  }
  else
  {

    /* Deallocation */

    if (s_mesh_old == nullptr) return (NULL);
    s_mesh = s_mesh_old;
    s_mesh->meshes = (int*) mem_free((char* ) s_mesh->meshes);
    s_mesh->points = (double*) mem_free((char* ) s_mesh->points);
    s_mesh->vertype = vertype_manage(-1, s_mesh->vertype, NULL, 0);
    s_mesh->vercoloc = vercoloc_manage(0, -1, NULL, NULL, 0, s_mesh->vercoloc);
    s_mesh = (SPDE_Mesh*) mem_free((char* ) s_mesh);
  }
  return (s_mesh);
}

/****************************************************************************/
/*!
 **  Define the function to transform a simulation
 **
 ** \param[in]  st_simu_transf  Pointer to the transformation function
 **
 *****************************************************************************/
void simu_define_func_transf(void (*st_simu_transf)(Db*, int, int, int))
{
  SIMU_FUNC_TRANSF = st_simu_transf;
}

/****************************************************************************/
/*!
 **  Define the function to account for the current simulation outcome
 **  in the calculation of the Modification arrays
 **
 ** \param[in]  st_simu_update  Pointer to the update function
 **
 *****************************************************************************/
void simu_define_func_update(void (*st_simu_update)(Db*, int, int, int))
{
  SIMU_FUNC_UPDATE = st_simu_update;
}

/****************************************************************************/
/*!
 **  Define the function to scale the Modification arrays
 **
 ** \param[in]  st_simu_scale  Pointer to the scaling function
 **
 *****************************************************************************/
void simu_define_func_scale(void (*st_simu_scale)(Db*, int, int))
{
  SIMU_FUNC_SCALE = st_simu_scale;
}

/****************************************************************************/
/*!
 **  Checks if there is a nugget component in the Model
 **
 ** \return 1 if a Nugget component is present; 0 otherwise
 **
 *****************************************************************************/
static int st_is_model_nugget(void)
{
  Model *model;

  model = st_get_model();

  for (int is = 0; is < model->getCovaNumber(); is++)
  {
    if (model->getCovaType(is) == ECov::NUGGET) return (1);
  }
  return (0);
}

/****************************************************************************/
/*!
 **  Returns the pointer to structure containing the Nugget Effect (or NULL)
 **
 *****************************************************************************/
static CovAniso* st_get_nugget(void)
{
  Model *model;
  CovAniso *cova;

  model = st_get_model();

  for (int is = 0; is < model->getCovaNumber(); is++)
  {
    cova = model->getCova(is);
    if (cova->getType() == ECov::NUGGET) return (cova);
  }
  return (nullptr);
}

/****************************************************************************/
/*!
 **  Returns the pointer to the structure
 **
 *****************************************************************************/
static CovAniso* st_get_cova(void)

{
  Model *model;
  CovAniso *cova;
  int is0, jcov;

  model = st_get_model();
  is0 = st_get_current_icov();

  for (int icov = jcov = 0; icov < model->getCovaNumber(); icov++)
  {
    cova = model->getCova(icov);
    if (cova->getType() == ECov::NUGGET) continue;
    if (is0 == jcov) return (cova);
    jcov++;
  }
  return (nullptr);
}

/****************************************************************************/
/*!
 **  Set the space dimension for the whole study
 **
 ** \param[in]  ndim   Space dimension
 **
 *****************************************************************************/
static void st_set_ndim(int ndim)
{
  S_ENV.ndim = ndim;
}

/****************************************************************************/
/*!
 **  Set the number of variables for the whole study
 **
 ** \param[in]  nvar   Number of variables
 **
 *****************************************************************************/
static void st_set_nvar(int nvar)
{
  S_ENV.nvar = nvar;
}

/****************************************************************************/
/*!
 **  Set the pointer to the model of the environment
 **
 ** \param[in]  model  Pointer to the Model structure
 **
 ** \remark  The pointer is copied. The contents is NOT duplicated
 **
 *****************************************************************************/
static void st_set_model(Model *model)
{
  st_get_current_ssenv()->model = model;
}

/****************************************************************************/
/*!
 **  Return the space dimension of the environment
 **
 *****************************************************************************/
static int st_get_ndim(void)
{
  return (S_ENV.ndim);
}

/****************************************************************************/
/*!
 **  Return the number of variables of the environment
 **
 *****************************************************************************/
static int st_get_nvar(void)
{
  return (S_ENV.nvar);
}

/****************************************************************************/
/*!
 **  Return the number of variables of the environment
 **
 *****************************************************************************/
static int st_get_nvs2(void)
{
  return (S_ENV.nvar * (S_ENV.nvar + 1) / 2);
}

/****************************************************************************/
/*!
 **  Return if a nugget effect component must be filtered
 **
 *****************************************************************************/
static int st_get_filnug(void)
{
  return (S_DECIDE.flag_filnug && S_DECIDE.flag_case == CASE_KRIGING);
}

/****************************************************************************/
/*!
 **  Defines if a nugget effect component must be filtered
 **
 ** \param[in]  flag_filnug  Flag to define if a nugget effect must be filtered
 **
 *****************************************************************************/
static void st_set_filnug(int flag_filnug)
{
  if (DEBUG) st_title(0, 0, -1, "(DEBUG) Set 'filnug'");
  S_DECIDE.flag_filnug = flag_filnug;
}

/****************************************************************************/
/*!
 **  Get the value of the Inverse of the sill for a given covariance and
 **  a pair of variables
 **
 *****************************************************************************/
static double st_get_isill(int icov, int ivar, int jvar)
{
  int nvar = st_get_nvar();
  const SPDE_Matelem &Maticov = spde_get_current_matelem(icov);
  double value = Maticov.Isill[(jvar) + nvar * (ivar)];
  return (value);
}

/****************************************************************************/
/*!
 **  Clean the Bhetero sparse matrices
 **
 *****************************************************************************/
static void st_clean_Bhetero(void)
{
  SPDE_SS_Environ *SS;

  SS = st_get_current_ssenv();

  /* Clean the vector of number of data / target per variable */

  SS->ndata1 = (int*) mem_free((char* ) SS->ndata1);
  SS->ntarget1 = (int*) mem_free((char* ) SS->ntarget1);

  /* Clean the sparse matrices for heterotopy at data points */

  if (SS->BheteroD != nullptr)
  {
    for (int ivar = 0; ivar < st_get_nvar(); ivar++)
      SS->BheteroD[ivar] = cs_spfree(SS->BheteroD[ivar]);
    SS->BheteroD = (cs**) mem_free((char* ) SS->BheteroD);
  }

  /* Clean the sparse matrices for heterotopy at targets */

  if (SS->BheteroT != nullptr)
  {
    for (int ivar = 0; ivar < st_get_nvar(); ivar++)
      SS->BheteroT[ivar] = cs_spfree(SS->BheteroT[ivar]);
    SS->BheteroT = (cs**) mem_free((char* ) SS->BheteroT);
  }
}

/****************************************************************************/
/*!
 **  Clean the Bnugget sparse matrices
 **
 *****************************************************************************/
static void st_clean_Bnugget(void)
{
  SPDE_SS_Environ *SS;

  SS = st_get_current_ssenv();
  if (SS->Bnugget != nullptr)
  {
    for (int i = 0; i < st_get_nvs2(); i++)
      SS->Bnugget[i] = cs_spfree(SS->Bnugget[i]);
    SS->Bnugget = (cs**) mem_free((char* ) SS->Bnugget);
  }
}

/****************************************************************************/
/*!
 **  Does the current sample correspond to a duplication
 **
 ** \return The rank of the joint sample (or -1 if no duplication)
 **
 ** \param[in]  vercoloc  Vercoloc structure
 ** \param[in]  mode      1 if the sample belongs to the Input Db
 **                       2 if the sample belongs to the Output Db
 ** \param[in]  iech      Rank of the sample
 **
 *****************************************************************************/
static int st_is_duplicated(Vercoloc *vercoloc, int mode, int iech)
{
  if (vercoloc == nullptr || vercoloc->ndupl <= 0) return (-1);

  if (mode == 1)
  {

    /* Sample belongs to the Input Db */

    if (S_DECIDE.flag_mesh_dbin) for (int i = 0; i < vercoloc->ndupl; i++)
      if (vercoloc->dupl_dabs[i] == iech) return (vercoloc->dupl_grid[i]);
  }
  else
  {

    /* Sample belongs to the Output Db */

    if (S_DECIDE.flag_mesh_dbout) for (int i = 0; i < vercoloc->ndupl; i++)
      if (vercoloc->dupl_grid[i] == iech) return (vercoloc->dupl_dabs[i]);
  }
  return (-1);
}

/****************************************************************************/
/*!
 **  Returns the array of indices in the meshing corresponding to
 **  the input Db
 **
 ** \return The array of indices
 **
 ** \param[in]  vertype      Vertype structure
 ** \param[in]  vercoloc     Vercoloc structure
 **
 ** \param[out] nbnodup      Number of not duplicated samples
 **
 ** \remarks The Array allocated must be freed by the calling function
 ** \remarks The returned array is numbered starting from 1
 **
 *****************************************************************************/
int* vercoloc_get_dbin_indices(Vertype *vertype,
                               Vercoloc *vercoloc,
                               int *nbnodup)
{
  int *indice, nech, pos, ndupl;

  /* Initializations */

  indice = nullptr;
  ndupl = vercoloc->ndupl;
  nech = (vertype->order == 1) ? vertype->nb1 :
                                 vertype->nb2 + ndupl;
  *nbnodup = nech - ndupl;

  /* Core allocation */

  indice = (int*) mem_alloc(sizeof(int) * nech, 0);
  if (indice == nullptr) return (indice);
  for (int i = 0; i < nech; i++)
    indice[i] = 0;

  /* Dispatch */

  if (vertype->order == 1)
  {

    // If order==1, samples of Dbin come first in the meshing

    for (int i = 0; i < nech; i++)
      indice[i] = i + 1;
  }
  else
  {

    // If order==2, samples of Dbout come first in the meshing
    // Therefore the samples (collocated to grid nodes) may come together
    // with the Dbout samples; followed by the remaining samples

    for (int i = 0; i < ndupl; i++)
      indice[vercoloc->dupl_data[i]] = vercoloc->dupl_grid[i] + 1;

    pos = vertype->nb1;
    for (int i = 0; i < nech; i++)
    {
      if (indice[i] != 0) continue;
      indice[i] = pos + 1;
      pos++;
    }
  }
  return (indice);
}

/****************************************************************************/
/*!
 **  Print the Mesh information
 **
 ** \param[in]  s_mesh    SPDE_Mesh structure
 **
 *****************************************************************************/
static void st_print_mesh(SPDE_Mesh *s_mesh)
{
  int ndim, nmesh, ncorner, nvertex;

  /* Initializations */

  ndim = s_mesh->ndim;
  nmesh = s_mesh->nmesh;
  ncorner = s_mesh->ncorner;
  nvertex = s_mesh->nvertex;

  // Title

  mestitle(0, "Mesh Information");
  message("Number of meshes   = %d\n", nmesh);
  message("Number of vertices = %d\n", nvertex);
  message("Space dimension    = %d\n", ndim);
  message("Number of vertices per Mesh = %d\n", ncorner);

  // List of Meshes
  for (int imesh = 0; imesh < nmesh; imesh++)
  {
    message("Mesh #%5d :", imesh + 1);
    for (int icorn = 0; icorn < ncorner; icorn++)
      message(" %5d", OLD_MESHES(imesh,icorn) + 1);
    message("\n");
  }

  message("\n");

  // List of Vertices
  for (int ip = 0; ip < nvertex; ip++)
  {
    message("Vertex #%5d : ", ip + 1);
    for (int idim = 0; idim < ndim; idim++)
      message(" %8.2lf", POINTS(ip, idim));
    message("\n");
  }
}

/****************************************************************************/
/*!
 **  Encode the status of the variable
 **
 ** \param[in]  auth    Status option
 **
 *****************************************************************************/
static void st_print_status(int auth)

{
  if (auth & VT_OTHER) message("OTHER ");
  if (auth & VT_FREE) message("FREE ");
  if (auth & VT_GIBBS) message("GIBBS ");
  if (auth & VT_HARD) message("DATA ");
  if (auth & VT_INPUT) message("INPUT ");
  if (auth & VT_OUTPUT) message("OUTPUT ");
}

/****************************************************************************/
/*!
 **  Print information on the Filter
 **
 ** \param[in]  title   Optional title
 ** \param[in]  auth    Filter option
 **
 *****************************************************************************/
static void st_qchol_filter(const char *title, int auth)
{
  message("%s = ", title);
  st_print_status(auth);
  message("\n");
}

/****************************************************************************/
/*!
 **  Print information on the Sparse matrix and its cholesky decomposition
 **
 ** \param[in]  title   Optional title
 ** \param[in]  QC      Pointer to the QChol structure
 **
 *****************************************************************************/
static void st_qchol_print(const char *title, QChol *QC)
{
  int nrows, ncols, count;
  double percent;

  if (QC == nullptr) return;

  if (title != NULL) message("%s\n", title);
  cs_rowcol(QC->Q, &nrows, &ncols, &count, &percent);
  message("- Nrows(%d) x Ncols(%d) - Non-zeros(%d) [%6.2lf (percent)]", nrows,
          ncols, count, percent);
  if (QC->S != NULL || QC->N != NULL) message(" (Cholesky)");
  message("\n");
}

/****************************************************************************/
/*!
 **  Check if the selection criterion is matched
 **
 ** \return  1 if the criterion is matched; 0 otherwise
 **
 ** \param[in]  vertype      Vertype value of the current sample
 ** \param[in]  vertype_auth Authorized vertype
 **
 *****************************************************************************/
static int st_ok(int vertype, int vertype_auth)
{
  if (vertype & vertype_auth) return (1);
  return (0);
}

/****************************************************************************/
/*!
 **  Return the count of samples belonging to a given type
 **
 ** \return  Number of samples of the target type of a given type
 **
 ** \param[in]  vertype      Vertype structure
 ** \param[in]  vertype_auth Authorized vertype
 **
 *****************************************************************************/
static int st_count_vertype(Vertype *vertype, int vertype_auth)
{
  int number;

  number = 0;
  for (int i = 0; i < vertype->nvertex; i++)
  {
    if (st_ok(vertype->vt[i], vertype_auth)) number++;
  }
  return (number);
}

/****************************************************************************/
/*!
 **  Construct the sparse matrix Q from another sparse matrix
 **
 ** \return The Q structure or NULL
 **
 ** \param[in]  Q_in      Input sparse matrix
 ** \param[in]  vertype   Vertype structure
 ** \param[in]  row_auth  Specification for rows extraction
 ** \param[in]  col_auth  Specification for columns extraction
 **
 ** \remarks The Cholesky decomposition is performed (if possible)
 **
 *****************************************************************************/
static cs* st_extract_Q_from_Q(cs *Q_in,
                               Vertype *vertype,
                               int row_auth,
                               int col_auth)
{
  int *rank_rows, *rank_cols, error, ecr_row, ecr_col;
  cs *Q = nullptr;

  /* Initializations */

  error = 1;
  rank_rows = rank_cols = nullptr;

  /* Core allocation */

  rank_rows = (int*) mem_alloc(sizeof(int) * vertype->nvertex, 0);
  if (rank_rows == nullptr) goto label_end;
  rank_cols = (int*) mem_alloc(sizeof(int) * vertype->nvertex, 0);
  if (rank_cols == nullptr) goto label_end;

  /* Fill the index vectors */

  ecr_row = ecr_col = 0;
  for (int i = 0; i < vertype->nvertex; i++)
  {
    rank_rows[i] = (st_ok(vertype->vt[i], row_auth)) ? ecr_row++ :
                                                       -1;
    rank_cols[i] = (st_ok(vertype->vt[i], col_auth)) ? ecr_col++ :
                                                       -1;
  }

  /* Extract the submatrix */

  Q = cs_extract_submatrix_by_ranks(Q_in, rank_rows, rank_cols);
  if (Q == nullptr) goto label_end;

  /* Set the error return code */

  error = 0;

  label_end: if (error) Q = cs_spfree(Q);
  rank_rows = (int*) mem_free((char* ) rank_rows);
  rank_cols = (int*) mem_free((char* ) rank_cols);
  return (Q);
}

/****************************************************************************/
/*!
 **  Construct the QChol sub-structure from another QChol structure
 **
 ** \return The QChol structure or NULL
 **
 ** \param[in]  title     Name of the QChol item
 ** \param[in]  QC_in     Input QChol structure
 ** \param[in]  vertype   Vertype structure
 ** \param[in]  row_auth  Specification for rows extraction
 ** \param[in]  col_auth  Specification for columns extraction
 **
 *****************************************************************************/
static QChol* st_extract_QC_from_Q(const char *title,
                                   QChol *QC_in,
                                   Vertype *vertype,
                                   int row_auth,
                                   int col_auth)
{
  int error;
  QChol *QC;

  /* Initializations */

  error = 1;
  QC = qchol_manage(1, nullptr);

  /* Extract the submatrix */

  QC->Q = st_extract_Q_from_Q(QC_in->Q, vertype, row_auth, col_auth);
  if (QC->Q == nullptr) goto label_end;

  /* Optional printout */

  if (VERBOSE)
  {
    message("Extracting a part of Q for '%s'\n", title);
    st_qchol_filter("- Row authorization code   ", row_auth);
    st_qchol_filter("- Column authorization code", col_auth);
    st_qchol_print(NULL, QC);
  }

  /* Set the error return code */

  error = 0;

  label_end: if (error) QC = qchol_manage(-1, QC);
  return (QC);
}

/****************************************************************************/
/*!
 **  Print the contents of the Vertype structure
 **
 ** \param[in]  vertype    Vertype structure
 **
 *****************************************************************************/
static void st_vertype_print(Vertype *vertype)

{
  if (vertype == nullptr) return;
  st_title(0, 0, 1, "Vertices Designation");
  message("- Total number      = %d\n", vertype->nvertex);
  message("- Free              = %d\n", st_count_vertype(vertype, VT_FREE));
  message("- Gibbs             = %d\n", st_count_vertype(vertype, VT_GIBBS));
  message("- Hard              = %d\n", st_count_vertype(vertype, VT_HARD));
  message("- From Input File   = %d\n", st_count_vertype(vertype, VT_INPUT));
  message("- From Output File  = %d\n", st_count_vertype(vertype, VT_OUTPUT));
  message("- Steiner Points    = %d\n", st_count_vertype(vertype, VT_OTHER));
  message("- Gibbs Points      = %d\n", vertype->ngibbs);
  if (vertype->order == 1)
    message(
        "Meshing order: 1) Input file, 2) Output file, 3) Steiner points\n");
  else
    message(
        "Meshing order: 1) Output file, 2) Input file, 3) Steiner points\n");
}

/****************************************************************************/
/*!
 **  Manage the QSimu structure
 **
 ** \return  Pointer to the QSimu structure
 **
 ** \param[in]  mode        Type of operation
 **                          1 : Allocation
 **                         -1 : Deallocation
 ** \param[in]  s_mesh      SPDE_Mesh structure
 ** \param[in]  qsimu       QSimu structure
 **
 *****************************************************************************/
static QSimu* st_qsimu_manage(int mode, SPDE_Mesh *s_mesh, QSimu *qsimu)
{
  int error;

  /* Initializations */

  error = 1;

  /* Dispatch */

  switch (mode)
  {
    case 1: /* Allocation */
      if (VERBOSE) st_title(0, 0, 1, "Building Environment");
      qsimu = (QSimu*) mem_alloc(sizeof(QSimu), 0);
      if (qsimu == nullptr) return (qsimu);

      /* Extract sub-matrices */

      if (S_DECIDE.flag_dbin)
      {
        qsimu->QCtt = st_extract_QC_from_Q("f_f",
                                           spde_get_current_matelem(-1).QC,
                                           s_mesh->vertype, VT_FREE, VT_FREE);
        if (qsimu->QCtt == nullptr) goto label_end;
        if (S_DECIDE.flag_mesh_dbin)
        {
          qsimu->QCtd = st_extract_QC_from_Q("f_gd",
                                             spde_get_current_matelem(-1).QC,
                                             s_mesh->vertype, VT_FREE,
                                             VT_GIBBS | VT_HARD);
          if (qsimu->QCtd == nullptr) goto label_end;
        }
      }
      break;

    case -1: /* Deallocation */
      if (qsimu == nullptr) return (qsimu);
      qsimu->QCtt = qchol_manage(-1, qsimu->QCtt);
      qsimu->QCtd = qchol_manage(-1, qsimu->QCtd);
      qsimu = (QSimu*) mem_free((char* ) qsimu);
      break;
  }

  /* Set the error return code */

  error = 0;

  label_end: if (error) qsimu = st_qsimu_manage(-1, NULL, qsimu);
  return (qsimu);
}

/****************************************************************************/
/*!
 **  Manage the Sparse matrix and its cholesky decomposition
 **
 ** \return  Pointer to the QChol structure
 **
 ** \param[in]  mode       Management mode
 **                         1 : Allocation
 **                        -1 : Total deallocation (NULL is returned)
 ** \param[in]  QC         Pointer to the QChol structure (when mode < 0)
 **
 *****************************************************************************/
QChol* qchol_manage(int mode, QChol *QC)
{

  /* Dispatch */

  switch (mode)
  {
    case 1: /* Allocation */
      QC = (QChol*) mem_alloc(sizeof(QChol), 1);
      QC->Q = nullptr;
      QC->S = nullptr;
      QC->N = nullptr;
      break;

    case -1: /* Total deallocation */
      if (QC == nullptr) return (QC);
      QC->Q = cs_spfree(QC->Q);
      QC->S = cs_sfree(QC->S);
      QC->N = cs_nfree(QC->N);
      QC = (QChol*) mem_free((char* ) QC);
      break;
  }
  return (QC);
}

/****************************************************************************/
/*!
 **  Return the sill of the Nugget Effect (or TEST)
 **
 ** \param[in] ivar    Rank of the first variable
 ** \param[in] jvar    Rank of the second variable
 **
 ** \remarks To save time, no check is performed with respect to the rank
 ** \remarks of the variables
 **
 *****************************************************************************/
static double st_get_nugget_sill(int ivar, int jvar)
{
  CovAniso *cova = st_get_nugget();
  if (cova == nullptr)
    return (TEST);
  else
    return (cova->getSill(ivar, jvar));
}

/****************************************************************************/
/*!
 **  Return the sill of the model (or TEST)
 **
 ** \param[in] ivar    Rank of the first variable
 ** \param[in] jvar    Rank of the second variable
 **
 ** \remarks To save time, no check is performed with respect to the rank
 ** \remarks of the structure or of the variables
 **
 *****************************************************************************/
static double st_get_cova_sill(int ivar, int jvar)
{
  CovAniso *cova = st_get_cova();
  return (cova->getSill(ivar, jvar));
}

/****************************************************************************/
/*!
 **  Return the param of the model
 **
 *****************************************************************************/
static double st_get_cova_param(void)
{
  return (st_get_cova()->getParam());
}

/****************************************************************************/
/*!
 **  Return the maximum number of covariances in the different Models
 **  for all GRFS (nugget included)
 **
 *****************************************************************************/
static int st_get_ncova_max(void)
{
  int ncova, ncova_max, igrf_memo;

  ncova_max = 0;
  igrf_memo = st_get_current_igrf();
  for (int igrf = 0; igrf < st_get_number_grf(); igrf++)
  {
    st_set_current_igrf(igrf);
    ncova = st_get_model()->getCovaNumber();
    if (ncova > ncova_max) ncova_max = ncova;
  }
  st_set_current_igrf(igrf_memo);
  return (ncova_max);
}

/****************************************************************************/
/*!
 **  Returns the number of structures in the Model (nugget excluded)
 **
 *****************************************************************************/
static int st_get_ncova(void)

{
  Model *model;
  CovAniso *cova;
  int ncova;

  ncova = 0;
  model = st_get_model();
  if (model == nullptr) return (ncova);
  for (int is = 0; is < model->getCovaNumber(); is++)
  {
    cova = model->getCova(is);
    if (cova->getType() != ECov::NUGGET) ncova++;
  }
  return (ncova);
}

/****************************************************************************/
/*!
 **  Return the number of vertices for the current Matelem
 **
 ** \param[in] icov    Rank of the target Covariance (or -1)
 **
 *****************************************************************************/
static int st_get_nvertex(int icov)

{
  return (spde_get_current_matelem(icov).s_mesh->nvertex);
}

/****************************************************************************/
/*!
 **  Return the maximum number of vertices for all meshes for all GRFs
 **
 *****************************************************************************/
static int st_get_nvertex_max(void)

{
  int nvertex, nvertex_max, igrf_memo;

  nvertex_max = 0;
  igrf_memo = st_get_current_igrf();
  for (int igrf = 0; igrf < st_get_number_grf(); igrf++)
  {
    st_set_current_igrf(igrf);
    for (int icov = 0; icov < st_get_ncova(); icov++)
    {
      nvertex = st_get_nvertex(icov);
      if (nvertex > nvertex_max) nvertex_max = nvertex;
    }
  }
  st_set_current_igrf(igrf_memo);
  return (nvertex_max);
}

/****************************************************************************/
/*!
 **  Return the dimension for array allocation
 **  This dimension is the aximum of the maximum number of vertices
 **  and the number of data
 **
 *****************************************************************************/
static int st_get_dimension(void)
{
  SPDE_SS_Environ *SS;
  int size;

  size = st_get_nvertex_max();
  for (int igrf = 0; igrf < SPDE_MAX_NGRF; igrf++)
  {
    SS = &S_ENV.SS_ENV[igrf];
    if (size < SS->ndata) size = SS->ndata;
  }
  return (size);
}

/****************************************************************************/
/*!
 **  Return the mean of the model
 **
 ** \param[in]  ivar    Rank of the target variable (for its mean)
 **
 *****************************************************************************/
static double st_get_model_mean(int ivar)
{
  return (st_get_model()->getMean(0));
}

/****************************************************************************/
/*!
 **  Get the normalized range
 **
 *****************************************************************************/
static double st_get_cova_range(void)
{
  return (st_get_cova()->getRange());
}

/****************************************************************************/
/*!
 **  Return the total sill of the model
 **
 ** \param[in] ivar    Rank of the first variable
 ** \param[in] jvar    Rank of the second variable
 **
 ** \remarks To save time, no check is performed with respect to the rank
 ** \remarks of the variables
 **
 *****************************************************************************/
static double st_get_sill_total(int ivar, int jvar)
{
  double total = 0.;

  CovAniso *cova = st_get_nugget();
  if (cova != nullptr) total += cova->getSill(ivar, jvar);

  for (int icov = 0; icov < st_get_ncova(); icov++)
  {
    st_set_current_icov(icov);
    cova = st_get_cova();
    total += cova->getSill(ivar, jvar);
  }
  return (total);
}

/****************************************************************************/
/*!
 **  Use the keypair mechanism to output some arrays
 **
 ** \param[in]  name   Name of the output
 ** \param[in]  iter   Concatenated in name if >= 0; ignored otherwise
 ** \param[in]  tab    Tab to be output
 **
 *****************************************************************************/
static void st_keypair_array(const char *name, int iter, double *tab)
{
  int nvar, ncova, ncur;

  if (!FLAG_KEYPAIR) return;
  nvar = st_get_nvar();
  ncova = st_get_ncova();
  ncur = st_get_nvertex_max();

  for (int icov = 0; icov < ncova; icov++)
    for (int ivar = 0; ivar < nvar; ivar++)
    {
      if (iter < 0)
        (void) gslSPrintf(NAME, "%s.%d.%d", name, icov + 1, ivar + 1);
      else
        (void) gslSPrintf(NAME, "%s.%d.%d.%d", name, iter + 1, icov + 1,
                          ivar + 1);
      set_keypair(NAME, 1, ncur, 1, &TAB(icov, ivar, 0));
    }
}

/****************************************************************************/
/*!
 **  Use the keypair mechanism to output some Sparse Matrix
 **
 ** \param[in]  name   Name of the output
 ** \param[in]  cs     Pointer to the sparse matrix
 ** \param[in]  i1     First argument
 ** \param[in]  i2     Second argument
 ** \param[in]  i3     Third argument
 ** \param[in]  i4     Fourth argument
 ** \param[in]  i5     Fifth argument
 **
 *****************************************************************************/
static void st_keypair_cs(const char *name,
                          cs *cs,
                          int i1,
                          int i2,
                          int i3,
                          int i4,
                          int i5)
{
  if (!FLAG_KEYPAIR) return;
  (void) gslStrcpy(NAME, name);
  if (i1 > 0) (void) gslSPrintf(NAME, "%s.%d", NAME, i1);
  if (i2 > 0) (void) gslSPrintf(NAME, "%s.%d", NAME, i2);
  if (i3 > 0) (void) gslSPrintf(NAME, "%s.%d", NAME, i3);
  if (i4 > 0) (void) gslSPrintf(NAME, "%s.%d", NAME, i4);
  if (i5 > 0) (void) gslSPrintf(NAME, "%s.%d", NAME, i5);
  cs_keypair(NAME, cs, 1);
}

/****************************************************************************/
/*!
 **  Print the Matelem characteristics (for given GRF and COV)
 **
 ** \param[in]  title   Title to be printed
 **
 *****************************************************************************/
static void st_print_all(const char *title)
{

  /* Initializations */

  int ndim = st_get_ndim();
  CovAniso *cova = st_get_cova();

  /* Print the title */

  st_title(1, 1, 1, title);

  /* Global parameters */

  message("Rank of the GRF       = %d\n", st_get_current_igrf() + 1);
  message("Rank of the structure = %d\n", st_get_current_icov() + 1);
  message("Param                 = %lf\n", st_get_cova_param());
  message("Alpha                 = %lf\n", st_get_cova_param() + ndim / 2.);
  message("Total Sill            = %lf\n", st_get_sill_total(0, 0));
  message("Ranges                = ");
  for (int idim = 0; idim < ndim; idim++)
    message("%lf ", st_get_cova_range() * cova->getAnisoCoeffs(idim));
  message("\n");

  /* 'H' Rotation */

  print_matrix("Anisotropy H matrix", 0, 1, ndim, ndim, NULL, Calcul.hh.data());
  message("Square root of Determinant                    = %lf\n",
          Calcul.sqdeth);
  message("Correction factor                             = %lf\n",
          Calcul.correc);

  /* Linear combination */

  int nblin = static_cast<int>(Calcul.blin.size());
  message("Number of terms in Linear Combination         = %d\n", nblin);
  print_matrix("Coefficients of the Linear Combination", 0, 1, 1, nblin, NULL,
               Calcul.blin.data());
}

/****************************************************************************/
/*!
 **  Compute the variance correction term
 **  Store in the SPDE_Calcul structure
 **
 *****************************************************************************/
static void st_compute_correc(void)

{
  int ndim = st_get_ndim();
  double param = st_get_cova_param();
  double value = spde_compute_correc(ndim, param);
  Calcul.correc = value;
}

double spde_compute_correc(int ndim, double param)

{
  double g0, ndims2, gammap, gammaa, value;

  ndims2 = ((double) ndim) / 2.;
  gammap = exp(loggamma(param));
  gammaa = exp(loggamma(param + ndims2));
  g0 = pow(4. * GV_PI, ndims2);
  value = gammap / (g0 * gammaa);
  return value;
}

/****************************************************************************/
/*!
 **  Compute the coefficients of the linear combination
 **
 ** \remarks This function stores the coefficients 'blin' in SPDE_Calcul
 **
 *****************************************************************************/
static void st_compute_blin(void)
{
  double ndims2, alpha, lambda, delta, correc, *m, *tp, *v;
  int p, ndimp;

  /* Initializations */

  int ndim = st_get_ndim();
  double param = st_get_cova_param();
  ndims2 = ((double) ndim) / 2.;
  alpha = param + ndims2;
  p = (int) ceil(alpha);
  ndimp = p + 1;
  lambda = alpha - floor(alpha);
  delta = lambda - alpha;
  correc = Calcul.correc;
  m = v = tp = nullptr;

  Calcul.blin.resize(NBLIN_TERMS, 0);

  if (lambda > 0.)
  {
    /* Core allocation */

    v = (double*) mem_alloc(sizeof(double) * ndimp, 1);
    m = (double*) mem_alloc(sizeof(double) * ndimp * ndimp, 1);
    tp = ut_pascal(ndimp);

    for (int idim = 0; idim < ndimp; idim++)
    {
      v[idim] = 1. / (2. * p - idim + delta);
      for (int jdim = 0; jdim < ndimp; jdim++)
        M(idim,jdim) = 1. / (2. * p - idim - jdim + lambda);
    }
    (void) matrix_invert(m, ndimp, -1);
    matrix_product(ndimp, ndimp, 1, m, v, v);
    matrix_product(ndimp, ndimp, 1, tp, v, Calcul.blin.data());
  }
  else
  {
    for (int i = 0; i <= p; i++)
      Calcul.blin[i] = ut_cnp(p, i) * correc;
  }

  Calcul.blin.resize(ndimp);

  /* Core deallocation */

  v = (double*) mem_free((char* ) v);
  m = (double*) mem_free((char* ) m);
  tp = (double*) mem_free((char* ) tp);
}

/****************************************************************************/
/*!
 **  Compute H matrix for anisotropic case and the square root of determinant
 **  Requires the knowledge of the actual parameters of the current Covariance
 **  Fills the SPDE_Calcul structure
 **
 *****************************************************************************/
static void st_compute_hh()
{

  /* Initializations */

  int ndim = st_get_ndim();
  CovAniso *cova = st_get_cova();
  VectorDouble temp(ndim * ndim, 0.);

  /* Processing */

  for (int i = 0; i < ndim; i++)
  {
    double scale = cova->getScale(i);
    if (Calcul.flag_sphere) scale /= Calcul.R;
    TEMP(ndim,i,i) = scale * scale;
  }
  matrix_prod_norme(1, ndim, ndim, cova->getAnisoRotMatVec().data(),
                    temp.data(), Calcul.hh.data());
}

/****************************************************************************/
/*!
 **  Initialize the contents of SPDE_Calcul structure
 **
 *****************************************************************************/
static void st_calcul_init(int ndim)
{
  variety_query(&Calcul.flag_sphere);

  Calcul.sqdeth = 0.;
  Calcul.correc = 0.;
  Calcul.R = 0.;
  Calcul.hh.resize(ndim * ndim, 0.);
  if (Calcul.flag_sphere)
  {
    variety_get_characteristics(&Calcul.R);
    Calcul.srot.resize(2, 0.);
  }
  Calcul.vv.resize(ndim, 0.);
}

/****************************************************************************/
/*!
 **  Update the contents of SPDE_Calcul structure
 **
 *****************************************************************************/
static void st_calcul_update(void)
{
  int ndim = st_get_ndim();

  // Check that the structure has already been initiated

  if (Calcul.hh.size() <= 0)
  my_throw("You should run 'st_calcul_init' beforehand");

  // Calculate the 'correc' term (from 'param')
  st_compute_correc();

  // Calculate the set of 'blin' coefficients (from 'param' and 'correc')
  st_compute_blin();

  // Calculate the 'HH' matrix
  st_compute_hh();

  // Calculate the determinant of HH
  Calcul.sqdeth = sqrt(matrix_determinant(ndim, Calcul.hh.data()));
}

/****************************************************************************/
/*!
 **  Modify the Exponential into a Bessel_K
 **
 ** \param[in]  cova         Covariance sructure
 **
 *****************************************************************************/
static void st_convert_exponential2bessel(CovAniso *cova)
{
  double scale_exp, range_exp, scale_bes, range_bes;

  if (cova->getType() != ECov::EXPONENTIAL) return;

  range_exp = cova->getRange();
  scale_exp = model_range2scale(ECov::EXPONENTIAL, range_exp, 0.);

  scale_bes = scale_exp;
  range_bes = model_scale2range(ECov::BESSEL_K, scale_bes, 0.5);

  cova->setType(ECov::BESSEL_K);
  cova->setParam(0.5);
  cova->setRange(range_bes);

  /* Optional printout */

  if (VERBOSE)
  {
    message("Convert from Exponential to Bessel-K\n");
    message("- Exponential: Range=%lf Scale=%lf\n", range_exp, scale_exp);
    message("- Bessel_K   : Range=%lf Scale=%lf\n", range_bes, scale_bes);
  }
  return;
}

/****************************************************************************/
/*!
 **  Attach the model (used to perform the assignments for external calls)
 **
 ** \return Error returned code
 **
 ** \param[in]  model        Model structure
 **
 *****************************************************************************/
int spde_attach_model(Model *model)

{
  CovAniso *cova;
  int ndim, nvar;
  double silltot;

  /* Check space dimension */

  if (model == nullptr) return (1);

  ndim = model->getDimensionNumber();
  nvar = model->getVariableNumber();

  if (ndim > 3)
  {
    messerr("The SPDE Methodology is implemented up to 3-D");
    return (1);
  }
  st_set_ndim(ndim);
  st_set_nvar(nvar);
  st_set_model(model);

  /* Checking the Model contents */

  silltot = 0.;
  for (int icov = 0; icov < model->getCovaNumber(); icov++)
  {
    cova = model->getCova(icov);
    silltot += cova->getSill(0, 0);
    if (cova->getType() == ECov::BESSEL_K)
    {
      continue;
    }
    else if (cova->getType() == ECov::EXPONENTIAL)
    {
      st_convert_exponential2bessel(cova);
      continue;
    }
    else if (cova->getType() == ECov::NUGGET)
    {
      if (model->getCova(icov)->getSill(0, 0) > 0)
        st_set_filnug(model->isCovaFiltered(icov));
    }
    else
    {
      messerr("SPDE Model can only support:");
      messerr("- Bessel_K basic structures");
      messerr("- Exponential basic structures");
      messerr("- A complementary Neugget Effect");
      return (1);
    }
  }
  if (st_get_ncova() <= 0)
  {
    messerr("The SPDE procedure requires at least one Bessel structure");
    return (1);
  }
  /* Check incompatibility between non-stationary and multivariate */

  if (st_get_nvar() > 1)
  {
    const ANoStat *nostat = st_get_model()->getNoStat();
    if (nostat != nullptr && nostat->isDefinedByType(-1, EConsElem::SILL))
    {
      messerr("Non-stationary Sill parameter incompatible with multivariate");
      return (1);
    }
  }

  return (0);
}

/****************************************************************************/
/*!
 **  Check that the Model is authorized for SPDE
 **
 ** \return Error returned code
 **
 ** \param[in]  dbin         Input Db structure
 ** \param[in]  dbout        Output Db structure
 ** \param[in]  model        Model structure
 **
 *****************************************************************************/
static int st_check_model(const Db *dbin, const Db *dbout, Model *model)
{
  CovAniso *cova;
  int ndim, nvar, flag_mult_data, flag_nugget;
  double silltot, nugval;

  /* Check space dimension */

  if (model == nullptr) return (1);

  ndim = model->getDimensionNumber();
  nvar = model->getVariableNumber();
  if (dbin != nullptr)
  {
    if (dbin->getNDim() != ndim)
    {
      messerr("Model (%d) and Input Db (%d) must have the same space dimension",
              ndim, dbin->getNDim());
      return (1);
    }
    flag_mult_data = (int) get_keypone("Flag_Mult_Data", 0);
    if (flag_mult_data)
    {
      if (nvar != 1)
      {
        messerr("The multiple variable used as entry");
        messerr("is only valid in the monovariate case");
      }
    }
    else
    {
      if (dbin->getVariableNumber() != nvar && S_DECIDE.flag_case
          != CASE_MATRICES
          && !S_DECIDE.flag_gibbs)
      {
        messerr(
            "Model (%d) and Input Db (%d) must refer to the same number of variables",
            nvar, dbin->getVariableNumber());
        return (1);
      }
    }
  }
  if (dbout != nullptr)
  {
    if (dbout->getNDim() != ndim)
    {
      messerr("Model(%d) and output Db(%d) must have same space dimension",
              ndim, dbout->getNDim());
      return (1);
    }
  }
  if (ndim != 2 && ndim != 3)
  {
    messerr("The SPDE Methodology is implemented for 2-D or 3-D case only");
    return (1);
  }
  st_set_ndim(ndim);
  st_set_nvar(nvar);
  st_set_model(model);

  /* Checking the Model contents */

  silltot = 0.;
  flag_nugget = 0;
  for (int icov = 0; icov < model->getCovaNumber(); icov++)
  {
    cova = model->getCova(icov);
    silltot += cova->getSill(0, 0);
    if (cova->getType() == ECov::BESSEL_K)
    {
      continue;
    }
    else if (cova->getType() == ECov::EXPONENTIAL)
    {
      st_convert_exponential2bessel(cova);
      continue;
    }
    else if (cova->getType() == ECov::NUGGET)
    {
      flag_nugget = 1;
      if (model->getSill(icov, 0, 0) > 0)
        st_set_filnug(model->isCovaFiltered(icov));
    }
    else
    {
      messerr("SPDE Model can only support:");
      messerr("- Bessel_K basic structures");
      messerr("- Exponential basic structures");
      messerr("- A complementary Nugget Effect");
      return (1);
    }
  }
  if (st_get_ncova() <= 0)
  {
    messerr("The SPDE procedure requires at least one Bessel structure");
    return (1);
  }

  /* If 'flag_mesh_dbin' is switched ON, Model must contain nugget Effect */

  if (S_DECIDE.flag_mesh_dbin && !flag_nugget)
  {
    nugval = silltot / 1000.;
    VectorDouble sill;
    sill.resize(nvar * nvar);
    int ecr = 0;
    for (int ivar = 0; ivar < nvar; ivar++)
      for (int jvar = 0; jvar < nvar; jvar++)
        sill[ecr++] = (ivar == jvar) ? nugval :
                                       0.;
    if (model_add_cova(model, ECov::NUGGET, 0, 0, 0., 0., VectorDouble(),
                       VectorDouble(), sill)) return (1);
  }

  /* Check incompatibility between non-stationary and multivariate */

  if (st_get_nvar() > 1)
  {
    const ANoStat *nostat = model->getNoStat();
    if (nostat != nullptr && nostat->isDefinedByType(-1, EConsElem::SILL))
    {
      messerr("Non-stationary Sill parameter incompatible with multivariate");
      return (1);
    }
  }

  if (st_get_ncova() > 1 || st_get_nvar() > 1 || st_is_model_nugget())
    S_DECIDE.flag_several = 1;

  return (0);
}

/****************************************************************************/
/*!
 **  Identify a parameter among the non-stationary ones
 **
 ** \return The rank of the parameter of -1 (not found)
 **
 ** \param[in]  icov0     Rank of the target covariance
 ** \param[in]  type0     Type of parameter (EConsElem)
 ** \param[in]  ivar0     Rank of the target variable (only when type=EConsElem::SILL)
 ** \param[in]  jvar0     Rank of the target variable (only when type=EConsElem::SILL)
 **
 ** \remark The covariance are ranked from 0 for non-nugget ones
 ** \remark The nugget effect corresponds to rank (-1)
 **
 *****************************************************************************/
static int st_identify_nostat_param(int icov0,
                                    const EConsElem &type0,
                                    int ivar0,
                                    int jvar0)
{
  const ANoStat *nostat = st_get_model()->getNoStat();
  if (nostat == nullptr) return -1;
  int igrf0 = st_get_current_igrf();
  int ipar = nostat->getRank(igrf0, icov0, type0, ivar0, jvar0);
  return ipar;
}

/****************************************************************************/
/*!
 **  Return the initial value given to a data under constraints
 **
 ** \return Initial value
 **
 ** \param[in]  db      Db structure
 ** \param[in]  igrf    Rank of the GRF
 ** \param[in]  iech    Sample rank
 **
 ** \remarks The test that Interval exists has already been performed
 **
 *****************************************************************************/
static double st_get_data_constraints(Db *db, int igrf, int iech)
{
  double vmin, vmax, pmin, pmax, simu, sk;

  vmin = db->getLowerBound(iech, igrf);
  vmax = db->getUpperBound(iech, igrf);
  sk = sqrt(st_get_sill_total(0, 0));
  pmin = law_cdf_gaussian(vmin);
  pmax = law_cdf_gaussian(vmax);
  simu = sk * law_invcdf_gaussian((pmin + pmax) / 2.);
  return (simu);
}

/****************************************************************************/
/*!
 **  Return the simulated value under constraints
 **
 ** \return Simulated value
 **
 ** \param[in]  db               Db structure
 ** \param[in]  igrf             Rank of the GRF
 ** \param[in]  iech             Rank of the sample
 ** \param[in]  iter0            Gibbs iteration rank (starting from 0)
 ** \param[in]  ngibbs_burn      Number of iterations (Burning step)
 ** \param[in]  yk               Kriged value
 ** \param[in]  sk               Standard deviation of the kriged error
 **
 ** \remarks In the gradual case, the constraints is calculated as a function
 ** \remarks of the iteration rank
 ** \remarks Otherwise the constant bounds are used
 **
 *****************************************************************************/
static double st_simu_constraints(Db *db,
                                  int igrf,
                                  int iech,
                                  int iter0,
                                  int ngibbs_burn,
                                  double yk,
                                  double sk)
{
  double vmin, vmax, ys, ratio, rndval, delta;

  /* Perform the simulation: read the (final) bound values */

  vmin = db->getLowerBound(iech, igrf);
  vmax = db->getUpperBound(iech, igrf);
  ratio =
      (iter0 < ngibbs_burn) ? (double) (ngibbs_burn - iter0 - 1) / (double) (iter0
                                  + 1) :
                              0.;
  if (FFFF(vmin))
    delta = ABS(vmax);
  else if (FFFF(vmax))
    delta = ABS(vmin);
  else
    delta = ABS(vmax - vmin);

  /* Lower bound */

  if (!FFFF(vmin))
  {
    vmin = vmin - delta * ratio;
    vmin = (vmin - yk) / sk;
  }

  /* Upper bound */

  if (!FFFF(vmax))
  {
    vmax = vmax + delta * ratio;
    vmax = (vmax - yk) / sk;
  }

  /* Draw the random gaussian value under constraints */

  rndval = law_gaussian_between_bounds(vmin, vmax);

  /* Simulate the value under constraints */

  ys = yk + sk * rndval;

  return (ys);
}

/****************************************************************************/
/*!
 **  Perform one iteration of the Gibbs sampler
 **
 ** \param[in]  igrf        Rank of the GRF
 ** \param[in]  vertype     Vertype structure
 ** \param[in]  ngibbs_int  Number of internal Gibbs iterations
 ** \param[in]  iter0       Gibbs iteration rank (starting from 0)
 ** \param[in]  ngibbs_burn Number of iterations (Burning step)
 ** \param[in]  dbin        Input Db
 ** \param[in]  dbout       Output Db
 **
 ** \param[out] zcur        Vector of simulation (VT_FREE|VT_GIBBS|VT_HARD)
 **
 ** \remarks The argument 'iter0' refers to the iteration rank, including
 ** \remarks 'ngibbs_burn' and 'ngibbs_iter'
 **
 *****************************************************************************/
static void st_gibbs(int igrf,
                     Vertype *vertype,
                     int ngibbs_int,
                     int iter0,
                     int ngibbs_burn,
                     Db *dbin,
                     Db *dbout,
                     double *zcur)
{
  int iech, jech, jg, p, niter, *Ap, *Ai;
  double *Ax, coeff, yk, sk;
  QChol *QC;
  cs *A;

  /* Initializations */

  QC = spde_get_current_matelem(-1).QC;
  A = QC->Q;
  Ap = A->p;
  Ai = A->i;
  Ax = A->x;
  sk = yk = 0.;
  iech = jech = 0;
  niter = MAX(1, ngibbs_int);

  /* Loop on the Gibbs samples */

  for (int iter = 0; iter < niter; iter++)
    for (int ig = 0; ig < vertype->ngibbs; ig++)
    {
      yk = 0.;
      iech = vertype->r_g[ig];
      for (p = Ap[iech]; p < Ap[iech + 1]; p++)
      {
        coeff = Ax[p];
        if (ABS(coeff) <= 0.) continue;
        jech = Ai[p];

        if (iech == jech)
          sk = coeff;
        else
          yk -= coeff * zcur[jech];
      }
      yk /= sk;
      sk = sqrt(1. / sk);
      jg = vertype->r_abs[ig];
      if (jg > 0)
        zcur[iech] = st_simu_constraints(dbout, igrf, jg - 1, iter0,
                                         ngibbs_burn, yk, sk);
      else
        zcur[iech] = st_simu_constraints(dbin, igrf, -jg - 1, iter0,
                                         ngibbs_burn, yk, sk);
    }

  if (DEBUG)
  {
    message("(DEBUG) Gibbs\n");
    print_range("- Result", vertype->nvertex, zcur, NULL);
  }
}

/****************************************************************************/
/*!
 **  Copy an array into the Db
 **
 ** \param[in]  vertype      Vertype structure
 ** \param[in]  z            Array of values
 ** \param[in]  dbout        Output Db structure
 ** \param[in]  locatorType       Rank of the pointer
 ** \param[in]  iatt_simu    Pointer to the current simulation
 **
 *****************************************************************************/
static void st_save_result(Vertype *vertype,
                           double *z,
                           Db *dbout,
                           const ELoc &locatorType,
                           int iatt_simu)
{
  int iech, lec, ecr;

  /* Loop on all the vertices */

  lec = ecr = 0;
  for (int ivar = 0; ivar < st_get_nvar(); ivar++)
  {
    iech = 0;
    for (int i = 0; i < vertype->nvertex; i++, lec++)
    {
      if (!st_ok(vertype->vt[i], VT_OUTPUT)) continue;
      while (!dbout->isActive(iech))
        iech++;
      set_LOCATOR_ITEM(dbout, locatorType, iatt_simu + ivar, iech, z[lec]);
      iech++;
      ecr++;
    }
  }
  if (DEBUG)
  {
    message("(DEBUG) Save ");
    st_print_status(VT_OUTPUT);
    message("\n");
    message("- Writing %d values (%d variable)\n", ecr, st_get_nvar());
  }
}

/****************************************************************************/
/*!
 **  Merge an auxiliary array into an permanent array
 **
 ** \param[in]  vertype      Vertype structure
 ** \param[in]  vertype_auth Authorized vertype
 ** \param[in]  z            Array of values to be merged
 **                          (only the ones whose vertype matches vertype_auth)
 **
 ** \param[out] zperm        Fill permanent array
 **
 *****************************************************************************/
static void st_merge(Vertype *vertype,
                     int vertype_auth,
                     double *z,
                     double *zperm)
{
  int ecr, lec;

  /* Loop on all the vertices */

  lec = ecr = 0;
  for (int ivar = 0; ivar < st_get_nvar(); ivar++)
  {
    for (int i = 0; i < vertype->nvertex; i++, ecr++)
    {
      if (st_ok(vertype->vt[i], vertype_auth))
      {
        zperm[ecr] = z[lec];
        lec++;
      }
    }
  }

  if (DEBUG)
  {
    message("(DEBUG) Merge ");
    st_print_status(vertype_auth);
    message("\n");
    print_range("- From  ", lec, z, NULL);
    print_range("- To    ", ecr, zperm, NULL);
  }
}

/****************************************************************************/
/*!
 **  Copy an auxiliary array into an permanent array
 **
 ** \param[in]  vertype      Vertype structure
 ** \param[in]  vertype_auth Authorized vertype
 ** \param[in]  z            Array of values to be merged
 **
 ** \param[out] zperm        Permanent array
 **
 ** \remarks This operation takes place for all samples located within the Part
 **
 *****************************************************************************/
static void st_copy(Vertype *vertype,
                    int vertype_auth,
                    double *z,
                    double *zperm)
{
  int ecr, lec;

  /* Loop on all the vertices */

  lec = ecr = 0;
  for (int ivar = 0; ivar < st_get_nvar(); ivar++)
  {
    for (int i = 0; i < vertype->nvertex; i++, ecr++)
    {
      if (vertype_auth == VT_NONE || st_ok(vertype->vt[i], vertype_auth))
      {
        zperm[ecr] = z[lec];
        lec++;
      }
    }
  }

  if (DEBUG)
  {
    message("(DEBUG) Copy ");
    st_print_status(vertype_auth);
    message("\n");
    print_range("- From  ", lec, z, NULL);
    print_range("- To    ", ecr, zperm, NULL);
  }
}

/****************************************************************************/
/*!
 **  Restores the conditional simulation array at mesh vertices
 **
 ** \param[in] number     : Number of terms in arrays 'perm' and 'aux'
 ** \param[in] zsnc       : Array containing the non-conditional simulation
 **
 ** \param[in,out] zcur   : Array containing the simulation error in input
 **                         and the conditional simulation in output
 **
 *****************************************************************************/
static void st_simu_add_vertices(int number, double *zsnc, double *zcur)
{
  for (int i = 0; i < number; i++)
    zcur[i] += zsnc[i];

  if (DEBUG)
  {
    message("(DEBUG) Add non-conditional simulation\n");
    print_range("- Result", number, zcur, NULL);
  }
}

/****************************************************************************/
/*!
 **  Load an auxiliary array into a permanent array: perm <- aux
 **
 ** \param[in] number     : Number of terms in arrays 'perm' and 'aux'
 ** \param[in] aux        : Auxiliary array
 ** \param[in,out] perm   : Input/Output array
 **
 *****************************************************************************/
static void st_load_array(int number, double *aux, double *perm)
{
  for (int i = 0; i < number; i++)
    perm[i] = aux[i];

  if (DEBUG)
  {
    message("(DEBUG) Loading\n");
    print_range("- Result", number, perm, NULL);
  }
}

/****************************************************************************/
/*!
 **  Initialize an array with a constant mean value
 **
 ** \param[in] ncova      Number of structures
 ** \param[in] nvar       Number of variables
 ** \param[in] ncur       Number of terms in arrays 'perm' and 'aux'
 ** \param[in] flag_var   1 : to initialize a variable;
 **                       0 : to initialize variance
 **
 ** \param[out] zperm     Output array
 **
 *****************************************************************************/
static void st_init_array(int ncova,
                          int nvar,
                          int ncur,
                          int flag_var,
                          double *zperm)
{
  double mean;
  int ecr;

  ecr = 0;
  for (int icov = 0; icov < ncova; icov++)
    for (int ivar = 0; ivar < nvar; ivar++)
    {
      mean = (flag_var) ? st_get_model_mean(ivar) / ncova :
                          0.;
      for (int icur = 0; icur < ncur; icur++)
      {
        zperm[ecr] = mean;
        ecr++;
      }
    }

  if (DEBUG)
  {
    message("(DEBUG) Initialize array\n");
    print_range("- Init  ", ncova * nvar * ncur, zperm, NULL);
  }
}

/****************************************************************************/
/*!
 **  Constitutes the array containing the simulation error
 **
 ** \return Error return code
 **
 ** \param[in]  vertype      Vertype structure
 ** \param[in]  vertype_auth Authorized vertype
 ** \param[in]  zsnc         Permanent array
 ** \param[in]  data         Array of datas
 **
 ** \param[in]  zerr         Array of extracted values
 **
 *****************************************************************************/
static int st_simu_subtract_data(Vertype *vertype,
                                 int vertype_auth,
                                 double *zsnc,
                                 double *data,
                                 double *zerr)
{
  int ecr;

  /* Loop on all the vertices */

  for (int i = ecr = 0; i < vertype->nvertex; i++)
  {
    if (!st_ok(vertype->vt[i], vertype_auth)) continue;
    zerr[ecr] = data[ecr] - zsnc[i];
    ecr++;
  }

  if (DEBUG)
  {
    message("(DEBUG) Subtracting non-conditional simulation\n");
    print_range("- Simu.Error", ecr, zerr, NULL);
  }

  return (0);
}

/****************************************************************************/
/*!
 **  Load the Vertype structure
 **
 ** \param[in]  vertype    Vertype structure
 ** \param[in]  vercoloc   Vercoloc structure
 ** \param[in]  dbin       Db structure for the conditioning data
 ** \param[in]  dbout      Db structure of the grid
 ** \param[in]  s_option   SPDE_Option structure
 **
 *****************************************************************************/
static void st_vertype_load(Vertype *vertype,
                            Vercoloc *vercoloc,
                            Db *dbin,
                            Db *dbout,
                            SPDE_Option &s_option)
{
  int ecr, ijoint, vertype_loc, ngibbs;

  /* Initializations */

  ecr = ngibbs = 0;
  vertype->order = 2;
  vertype->nb1 = 0;
  vertype->nb2 = 0;

  /* From the Output Db if present */

  if (dbout != nullptr && S_DECIDE.flag_mesh_dbout)
  {
    for (int iech = 0; iech < dbout->getSampleNumber(); iech++)
    {
      if (!dbout->isActive(iech)) continue;
      vertype_loc = VT_FREE;
      ijoint = st_is_duplicated(vercoloc, 2, iech);
      if (ijoint < 0)
      {
        if (S_DECIDE.flag_gibbs && dbout->getIntervalNumber() > 0)
          vertype_loc = VT_GIBBS;
        else
          vertype_loc = VT_FREE;
      }
      else
      {
        if (S_DECIDE.flag_gibbs && dbin->getIntervalNumber() > 0)
          vertype_loc = VT_GIBBS;
        else if (S_DECIDE.flag_mesh_dbin) vertype_loc = VT_HARD | VT_INPUT;
      }
      if (vertype_loc == VT_GIBBS)
      {
        vertype->r_g[ngibbs] = ecr;
        vertype->r_abs[ngibbs] = (ijoint < 0) ? iech + 1 :
                                                -(ijoint + 1);
        ngibbs++;
      }
      vertype->vt[ecr++] = VT_OUTPUT | vertype_loc;
      vertype->nb1++;
    }
  }

  /* Input file */

  if (dbin != nullptr && S_DECIDE.flag_mesh_dbin)
  {
    for (int iech = 0; iech < dbin->getSampleNumber(); iech++)
    {
      if (!dbin->isActive(iech)) continue;
      if (st_is_duplicated(vercoloc, 1, iech) >= 0) continue;
      if (S_DECIDE.flag_gibbs && dbin->getIntervalNumber() > 0)
        vertype_loc = VT_GIBBS;
      else
        vertype_loc = VT_HARD;
      if (vertype_loc == VT_GIBBS)
      {
        vertype->r_g[ngibbs] = ecr;
        vertype->r_abs[ngibbs] = -(iech + 1);
        ngibbs++;
      }
      vertype->vt[ecr++] = VT_INPUT | vertype_loc;
      vertype->nb2++;
    }
  }
  vertype->ngibbs = ngibbs;

  /* Optional printout */

  if (VERBOSE) st_vertype_print(vertype);
}

/****************************************************************************/
/*!
 **  Manage the Vertype structure
 **
 ** \return  Pointer to the newly created array (or NULL if failure)
 **
 ** \param[in]  mode       Operation
 **                         1 : Create
 **                        -1 : Delete
 ** \param[in]  vertype    Vertype structure
 ** \param[in]  vercoloc   Vercoloc structure
 ** \param[in]  nvertex    Number of vertices
 **
 *****************************************************************************/
Vertype* vertype_manage(int mode,
                        Vertype *vertype,
                        Vercoloc *vercoloc,
                        int nvertex)
{
  int error, ngibbs;

  /* Dispatch */

  error = 1;
  if (mode > 0)
  {

    /* Allocation */

    vertype = (Vertype*) mem_alloc(sizeof(Vertype), 0);
    if (vertype == nullptr) return (vertype);
    vertype->nvertex = nvertex;
    vertype->ngibbs = ngibbs = 0;
    vertype->vt = nullptr;
    vertype->r_g = nullptr;
    vertype->r_abs = nullptr;
    vertype->vt = (int*) mem_alloc(sizeof(int) * nvertex, 0);
    if (vertype->vt == nullptr) goto label_end;
    vertype->r_g = (int*) mem_alloc(sizeof(int) * nvertex, 0);
    if (vertype->r_g == nullptr) goto label_end;
    vertype->r_abs = (int*) mem_alloc(sizeof(int) * nvertex, 0);
    if (vertype->r_abs == nullptr) goto label_end;
    for (int i = 0; i < vertype->nvertex; i++)
      vertype->vt[i] = VT_OTHER | VT_FREE;
  }
  else
  {

    /* Deallocation */

    if (vertype == nullptr) return (vertype);
    vertype->r_g = (int*) mem_free((char* ) vertype->r_g);
    vertype->r_abs = (int*) mem_free((char* ) vertype->r_abs);
    vertype->vt = (int*) mem_free((char* ) vertype->vt);
    vertype = (Vertype*) mem_free((char* ) vertype);
  }

  /* Set the error return code */

  error = 0;

  label_end: if (error) vertype = vertype_manage(-1, vertype, NULL, 0);
  return (vertype);
}

/****************************************************************************/
/*!
 **  Prepare the RHS for Kriging
 **
 ** \param[in]  QCtd        Pointer to QChol structure (target-data)
 ** \param[in]  data        Input array (Dimension: ndata)
 ** \param[in]  ntarget     Number of target values
 **
 ** \param[out] rhs         Output array  (Dimension: ntarget)
 **
 *****************************************************************************/
static void st_kriging_one_rhs(QChol *QCtd,
                               double *data,
                               int ntarget,
                               double *rhs)
{
  cs_mulvec(QCtd->Q, ntarget, data, rhs);
  for (int i = 0; i < ntarget; i++)
    rhs[i] = -rhs[i];
}

/****************************************************************************/
/*!
 **  Perform the Calculation of the Kriging estimate
 **
 ** \return  Error return code
 **
 ** \param[in]  QC          Pointer to QChol structure
 ** \param[in]  rhs         R.H.S. array (Dimension: ntarget)
 **
 ** \param[out] work        Working array (Dimension: ntarget)
 ** \param[out] z           Output array  (Dimension: ntarget)
 **
 *****************************************************************************/
static int st_kriging_cholesky(QChol *QC, double *rhs, double *work, double *z)
{
  int ntarget;

  /* Initializations */

  ntarget = QC->Q->n;
  for (int icur = 0; icur < ntarget; icur++)
    work[icur] = 0.;

  /* Prepare Cholesky decomposition (if not already performed) */

  if (QC->S == nullptr)
  {
    if (qchol_cholesky(VERBOSE, QC)) return (1);
  }

  /* Process the Cholesky inversion */

  cs_chol_invert(QC, z, rhs, work);

  /* Optional debugging information */

  if (DEBUG)
  {
    message("(DEBUG) Kriging (Cholesky)\n");
    print_range("- Result", ntarget, z, NULL);
  }
  return (0);
}

/****************************************************************************/
/*!
 **  Perform the Filtering of Nugget Effect
 **
 ** \return  Error return code
 **
 ** \param[out] work        Working array (Dimension: nvertex)
 ** \param[in]  y           Output Array (Dimension: nvertex)
 **
 *****************************************************************************/
static int st_filter(double *work, double *y)
{
  QChol *QC;
  int ntarget;

  /* Initializations */

  QC = spde_get_current_matelem(-1).QC;
  ntarget = QC->Q->n;
  for (int icur = 0; icur < ntarget; icur++)
    work[icur] = 0.;

  /* Prepare Cholesky decomposition (if not already performed) */

  if (QC->S == nullptr)
  {
    if (qchol_cholesky(VERBOSE, QC)) return (1);
  }

  /* Process the Cholesky inversion */

  cs_chol_invert(QC, y, y, work);

  /* Optional debugging information */

  if (DEBUG)
  {
    message("(DEBUG) Filtering\n");
    print_range("- Result", ntarget, y, NULL);
  }
  return (0);
}

/****************************************************************************/
/*!
 **  Perform the calculation of the Standard Deviation of Estimation Error
 **
 ** \return  Error return code
 **
 ** \param[out] vcur     Output array
 **
 *****************************************************************************/
int spde_build_stdev(double *vcur)
{
  int *wZdiagp, *wLmunch, error, nzmax, ntarget;
  double *d2, *wz, *diag, *z;
  cs *Dinv, *LDinv, *TLDinv, *Pattern;
  QChol *QCtt;
  SPDE_Mesh *s_mesh;

  /* Initializations */

  error = 1;
  SPDE_Matelem &Matelem = spde_get_current_matelem(-1);
  s_mesh = Matelem.s_mesh;
  QCtt = Matelem.qsimu->QCtt;
  wZdiagp = wLmunch = nullptr;
  d2 = wz = diag = z = nullptr;
  Dinv = LDinv = TLDinv = Pattern = nullptr;
  ntarget = QCtt->Q->n;

  // Perform the Cholesky (if not already done) */

  if (qchol_cholesky(0, QCtt)) goto label_end;

  /* Pre-processing */

  d2 = csd_extract_diag(QCtt->N->L, 2);
  if (d2 == nullptr) goto label_end;
  Dinv = cs_extract_diag(QCtt->N->L, -1);
  if (Dinv == nullptr) goto label_end;
  LDinv = cs_multiply(QCtt->N->L, Dinv);
  if (LDinv == nullptr) goto label_end;
  TLDinv = cs_transpose(LDinv, 1);
  if (TLDinv == nullptr) goto label_end;
  Pattern = cs_add(LDinv, TLDinv, 1, 1);
  if (Pattern == nullptr) goto label_end;
  if (cs_sort_i(Pattern)) goto label_end;
  if (cs_sort_i(LDinv)) goto label_end;

  /* Core allocation */

  nzmax = Pattern->nzmax;
  z = (double*) mem_alloc(sizeof(double) * ntarget, 0);
  if (z == nullptr) goto label_end;
  wz = (double*) mem_alloc(sizeof(double) * nzmax, 0);
  if (wz == nullptr) goto label_end;
  wZdiagp = (int*) mem_alloc(sizeof(int) * nzmax, 0);
  if (wZdiagp == nullptr) goto label_end;
  wLmunch = (int*) mem_alloc(sizeof(int) * nzmax, 0);
  if (wLmunch == nullptr) goto label_end;
  for (int i = 0; i < nzmax; i++)
    wz[i] = 0.;

  if (sparseinv(ntarget, LDinv->p, LDinv->i, LDinv->x, d2, LDinv->p, LDinv->i,
                LDinv->x, Pattern->p, Pattern->i, Pattern->x, wz, wZdiagp,
                wLmunch)
      == -1) goto label_end;

  /* Extracting the diagonal of wz */

  diag = csd_extract_diag(Pattern, 1);
  cs_pvec(ntarget, QCtt->S->Pinv, diag, z);
  for (int iech = 0; iech < ntarget; iech++)
    z[iech] = sqrt(z[iech]);
  st_merge(s_mesh->vertype, VT_FREE, z, vcur);

  /* Set the error return code */

  error = 0;

  label_end: wZdiagp = (int*) mem_free((char* ) wZdiagp);
  wLmunch = (int*) mem_free((char* ) wLmunch);
  wz = (double*) mem_free((char* ) wz);
  d2 = (double*) mem_free((char* ) d2);
  diag = (double*) mem_free((char* ) diag);
  z = (double*) mem_free((char* ) z);
  Dinv = cs_spfree(Dinv);
  LDinv = cs_spfree(LDinv);
  TLDinv = cs_spfree(TLDinv);
  Pattern = cs_spfree(Pattern);
  return (error);
}

/****************************************************************************/
/*!
 **  Calculate the array of dimensions of the meshes
 **
 ** \return Pointer to the newly allocated array containing mesh dimensions
 **
 ** \param[in]  amesh    MeshEStandard structure
 **
 ** \remark The array returned by this function must be deallocated
 **
 *****************************************************************************/
double* _spde_get_mesh_dimension(MeshEStandard *amesh)

{
  double *units, mat[9];
  int flag_sphere;

  /* Initializations */

  units = nullptr;
  int ndim = amesh->getNDim();
  int nmesh = amesh->getNMeshes();
  int ncorner = amesh->getNApexPerMesh();
  variety_query(&flag_sphere);

  /* Core allocation */

  units = (double*) mem_alloc(sizeof(double) * nmesh, 0);
  if (units == nullptr) return (units);

  /* Dispatch */

  if (flag_sphere)
  {
    for (int imesh = 0; imesh < nmesh; imesh++)
    {
      units[imesh] = ut_geodetic_triangle_surface(amesh->getCoor(imesh, 0, 0),
                                                  amesh->getCoor(imesh, 0, 1),
                                                  amesh->getCoor(imesh, 1, 0),
                                                  amesh->getCoor(imesh, 1, 1),
                                                  amesh->getCoor(imesh, 2, 0),
                                                  amesh->getCoor(imesh, 2, 1));
    }
  }
  else
  {
    for (int imesh = 0; imesh < nmesh; imesh++)
    {
      int ecr = 0;
      for (int icorn = 1; icorn < ncorner; icorn++)
        for (int idim = 0; idim < ndim; idim++)
          mat[ecr++] = (amesh->getCoor(imesh, icorn, idim)
              - amesh->getCoor(imesh, 0, idim));
      units[imesh] = ABS(matrix_determinant(ndim,mat)) / FACDIM[ndim];
    }
  }
  return (units);
}

/****************************************************************************/
/*!
 **  Update parameters in S_ENV structure in non-stationary case
 **
 ** \param[in]  amesh     MeshEStandard structure
 ** \param[in]  imesh0    Rank of the current mesh
 **
 *****************************************************************************/
static void st_calcul_update_nostat(MeshEStandard *amesh, int imesh0)

{
  Model *model = st_get_model();
  const ANoStat *nostat = model->getNoStat();

  /* Initializations */

  int ndim = st_get_ndim();
  int igrf0 = st_get_current_igrf();
  int icov0 = st_get_current_icov();
  int ncorner = amesh->getNApexPerMesh();

  /* Update the Tensor 'hh' */

  if (nostat->isDefinedforAnisotropy(igrf0, icov0))
  {
    VectorDouble hhtot(ndim * ndim, 0.);
    for (int ic = 0; ic < ncorner; ic++)
    {
      nostat->updateModel(model, amesh->getApex(imesh0, ic));
      st_compute_hh();
      ut_vector_cumul(hhtot, Calcul.hh, 1.);
    }
    ut_vector_divide_inplace(hhtot, (double) ncorner);
    ut_vector_copy(Calcul.hh, hhtot);
    Calcul.sqdeth = sqrt(matrix_determinant(ndim, Calcul.hh.data()));
  }

  /* Update the Spherical Rotation array */

  if (nostat->isDefined(igrf0, icov0, EConsElem::SPHEROT, -1, -1))
  {
    VectorDouble srot(2, 0.);
    for (int i = 0; i < 2; i++)
    {
      int ipar = nostat->getRank(igrf0, icov0, EConsElem::SPHEROT, i, -1);
      if (ipar < 0) continue;
      double total = 0.;
      for (int ic = 0; ic < ncorner; ic++)
        total += nostat->getValue(ipar, 0, amesh->getApex(imesh0, ic));
      Calcul.srot[i] = total / (double) ncorner;
    }
  }

  /* Update the Velocity array */

  if (nostat->isDefined(igrf0, icov0, EConsElem::VELOCITY, -1, -1))
  {
    VectorDouble vv(ndim, 0.);
    for (int idim = 0; idim < ndim; idim++)
    {
      int ipar = nostat->getRank(igrf0, icov0, EConsElem::VELOCITY, idim, -1);
      if (ipar < 0) continue;
      double total = 0.;
      for (int ic = 0; ic < ncorner; ic++)
        total += nostat->getValue(ipar, 0, amesh->getApex(imesh0, ic));
      Calcul.vv[idim] = total / (double) ncorner;
    }
  }
}

/****************************************************************************/
/*!
 **  Fill the Isill matrix linked to the covariance of the Model
 **
 ** \return Error returned code
 **
 ** \remark The matrix 'Isill' is dimensioned to nvar * nvar where
 **
 *****************************************************************************/
static int st_fill_Isill(void)
{
  double *mcova;
  int nvar, nvar2, error, icov, ecr;

  /* Initializations */

  error = 1;
  nvar = st_get_nvar();
  nvar2 = nvar * nvar;
  mcova = nullptr;
  icov = st_get_current_icov();
  SPDE_Matelem &Matelem = spde_get_current_matelem(icov);

  /* Core allocation */

  mcova = (double*) mem_alloc(sizeof(double) * nvar2, 0);
  if (mcova == nullptr) goto label_end;

  /* Load the sill of the covariance */

  ecr = 0;
  for (int ivar = 0; ivar < nvar; ivar++)
    for (int jvar = 0; jvar < nvar; jvar++)
      mcova[ecr++] = st_get_cova_sill(ivar, jvar);

  /* Loop on the structures to invert the sill matrices */

  if (matrix_invert(mcova, nvar, -1)) goto label_end;

  /* Optional printout */

  if (VERBOSE) message("Calculation of Isill\n");

  /* Set the error return code */

  error = 0;

  label_end: if (error) mcova = (double*) mem_free((char* ) mcova);
  Matelem.Isill = mcova;
  return (error);
}

/****************************************************************************/
/*!
 **  Fill the Csill matrix linked to the continuous parts of the Model
 **
 ** \return Error returned code
 **
 ** \remark The matrix 'Csill' is dimensioned to ncova * nvar * (nvar+1)/2 where
 ** \remark - ncova designates the number of continuous structures of the Model
 **
 *****************************************************************************/
static int st_fill_Csill(void)
{
  Model *model;
  double *mcova;
  int nvar, nvs2, error, icov;

  /* Initializations */

  error = 1;
  model = st_get_model();
  nvar = st_get_nvar();
  nvs2 = nvar * (nvar + 1) / 2;
  mcova = nullptr;
  icov = st_get_current_icov();
  SPDE_Matelem &Matelem = spde_get_current_matelem(icov);

  /* Core allocation */

  mcova = (double*) mem_alloc(sizeof(double) * nvs2, 0);
  if (mcova == nullptr) goto label_end;

  /* Load the sills of continuous covariance elements */

  if (matrix_cholesky_decompose(
      model->getCova(icov)->getSill().getValues().data(), mcova, nvar))
    goto label_end;

  /* Optional printout */

  if (VERBOSE) message("Calculation of Csill\n");

  /* Set the error return code */

  error = 0;

  label_end: if (error) mcova = (double*) mem_free((char* ) mcova);
  Matelem.Csill = mcova;
  return (error);
}

/****************************************************************************/
/*!
 **  Fill the Bnugget sparse matrix linked to nugget effect
 **
 ** \return Error returned code
 **
 ** \param[in]  dbin      Db structure
 **
 ** \remark This function allocates 'nvs2' sparse matrices of dimension 'ndata'.
 ** \remark where nvs2 is the product nvar * (nvar+1) / 2
 **
 *****************************************************************************/
static int st_fill_Bnugget(Db *dbin)

{
  double *mat, *local, *local0;
  int *ind, error, ndata, nvar, nvs2, nvar2, size, ecr, nvr, ivar, jvar, iad;
  int flag_nostat_sillnug;
  Model *model;
  cs **Bnugget;

  /* Initializations */

  error = 1;
  model = st_get_model();
  ndata = dbin->getActiveSampleNumber();
  nvar = model->getVariableNumber();
  nvar2 = nvar * nvar;
  nvs2 = nvar * (nvar + 1) / 2;
  mat = local = local0 = nullptr;
  ind = nullptr;
  Bnugget = nullptr;

  /* In the non-stationary case, identify the rank of the parameter */
  /* which corresponds to the sill of the nugget effect */

  flag_nostat_sillnug = st_identify_nostat_param(-1, EConsElem::SILL, -1, -1)
      >= 0;
  if (flag_nostat_sillnug)
  {
    messerr("Non-stationarity on nugget sill values not programmed yet");
    goto label_end;
  }

  /* Core allocation */

  size = ndata * nvs2;
  local = (double*) mem_alloc(sizeof(double) * nvar2, 0);
  if (local == nullptr) goto label_end;
  local0 = (double*) mem_alloc(sizeof(double) * nvar2, 0);
  if (local0 == nullptr) goto label_end;
  ind = (int*) mem_alloc(sizeof(int) * ndata, 0);
  if (ind == nullptr) goto label_end;
  mat = (double*) mem_alloc(sizeof(double) * size, 0);
  if (mat == nullptr) goto label_end;
  for (int i = 0; i < size; i++)
    mat[i] = 0.;

  /* Establish the nugget sill matrix for isotopic case (only in stationary) */

  if (!flag_nostat_sillnug)
  {
    for (ivar = 0; ivar < nvar; ivar++)
      for (jvar = 0; jvar < nvar; jvar++)
        LOCAL0(ivar,jvar) = st_get_nugget_sill(ivar, jvar);
    if (matrix_invert(local0, nvar, -1))
    {
      messerr("Problem when inverting the Global Nugget matrix of sill");
      goto label_end;
    }
  }

  /* Loop on the active samples */

  ecr = 0;
  for (int iech = 0; iech < dbin->getSampleNumber(); iech++)
  {
    if (!dbin->isActive(iech)) continue;

    /* Check the heterotopy for the nugget effect */

    nvr = 0;
    for (ivar = 0; ivar < nvar; ivar++)
    {
      if (FFFF(dbin->getVariable(iech, ivar))) continue;
      ind[nvr] = ivar;
      nvr++;
    }
    if (nvr <= 0)
    {
      messerr("For sample %#d, no variable is defined", iech + 1);
      goto label_end;
    }

    /* Dispatch */

    if (nvr == nvar && !flag_nostat_sillnug)
    {

      /* Isotopic case: Store the sill partial matrix */

      for (ivar = 0; ivar < nvar; ivar++)
        for (jvar = 0; jvar <= ivar; jvar++)
        {
          iad = st_get_rank(ivar, jvar);
          mat[iad * ndata + ecr] = LOCAL0(ivar, jvar);
        }
    }
    else
    {

      /* Constitute the sill matrix for the nugget effect */

      for (int ivr = 0; ivr < nvr; ivr++)
        for (int jvr = 0; jvr < nvr; jvr++)
          LOCAL(ivr,jvr) = st_get_nugget_sill(ind[ivr], ind[jvr]);

      /* Invert the sill partial matrix */

      if (matrix_invert(local, nvr, -1))
      {
        messerr("Problem when inverting Nugget matrix of sill at sample #%d",
                iech + 1);
        goto label_end;
      }

      /* Store the sill partial matrix */

      for (int ivr = 0; ivr < nvr; ivr++)
        for (int jvr = 0; jvr <= ivr; jvr++)
        {
          ivar = ind[ivr];
          jvar = ind[jvr];
          iad = st_get_rank(ivar, jvar);
          mat[iad * ndata + ecr] = LOCAL(ivr, jvr);
        }
    }
    ecr++;
  }

  /* Define the sparse matrices */

  Bnugget = (cs**) mem_alloc(sizeof(cs*) * nvs2, 0);
  if (Bnugget == nullptr) goto label_end;
  for (int ivs2 = 0; ivs2 < nvs2; ivs2++)
    Bnugget[ivs2] = nullptr;
  ecr = 0;
  for (ivar = 0; ivar < nvar; ivar++)
    for (jvar = 0; jvar <= ivar; jvar++, ecr++)
      Bnugget[ecr] = cs_eye_tab(ndata, &mat[ecr * ndata]);

  /* Optional printout */

  if (VERBOSE) message("Calculation of Bnugget (%d sparse matrices)\n", nvs2);

  /* Set the error return code */

  error = 0;

  label_end: if (error) st_clean_Bnugget();
  st_get_current_ssenv()->Bnugget = Bnugget;
  st_get_current_ssenv()->ndata = ndata;
  ind = (int*) mem_free((char* ) ind);
  local = (double*) mem_free((char* ) local);
  local0 = (double*) mem_free((char* ) local0);
  mat = (double*) mem_free((char* ) mat);
  return (error);
}

/****************************************************************************/
/*!
 **  Return the list of (target+data) indices for a given mesh
 **
 ** \return An array of vertex identification (Dimension: nvertex) or NULL
 **
 ** \param[in]  s_mesh      SPDE_Mesh structure
 ** \param[in]  dbin        Db input structure
 ** \param[in]  dbout       Db output structure
 **
 ** \remarks The array ranks is filled as follows:
 ** \remarks - Its contents follows the mesh numbering
 ** \remarks - If positive, its value provides the rank of the data
 ** \remarks - If negative, its absolute valu provides the rank of the target
 ** \remarks - If zero, theses are Steiner points
 ** \remarks Warning: Ranks are counted from 1
 **
 ** \remarks The returned array must be freed by the calling function
 ** \remarks Dimension: nvertex
 **
 *****************************************************************************/
static int* st_get_vertex_ranks(SPDE_Mesh *s_mesh, Db *dbin, Db *dbout)
{
  Vertype *vertype;
  int *ranks, nvertex, ndata, ngrid;

  /* Initializations */

  ranks = nullptr;
  nvertex = s_mesh->nvertex;
  vertype = s_mesh->vertype;

  /* Core allocation */

  ranks = (int*) mem_alloc(sizeof(int) * nvertex, 0);
  if (ranks == nullptr) return (ranks);

  /* Identify the vertices */

  ndata = ngrid = 0;
  for (int i = 0; i < nvertex; i++)
  {
    ranks[i] = 0;
    if (vertype->vt[i] & VT_INPUT)
    {
      ranks[i] = (ndata + 1);
      ndata++;
    }
    else if (vertype->vt[i] & VT_OUTPUT)
    {
      ranks[i] = -(ngrid + 1);
      ngrid++;
    }
  }
  return (ranks);
}

/****************************************************************************/
/*!
 **  Fill some matrices for Kriging in the case of a model without nugget effect
 **  Constitute the Bhetero sparse matrices
 **
 ** \return Error returned code
 **
 ** \param[in]  dbin      Input Db structure
 ** \param[in]  dbout     Output Db structure
 **
 *****************************************************************************/
static int st_fill_Bhetero(Db *dbin, Db *dbout)

{
  int *ranks, *ndata1, *ntarget1;
  int ndata, nvar, ecrT, nvertex, flag_add, iech, error;
  double value;
  Model *model;
  cs **BheteroD, **BheteroT, *Btriplet;
  SPDE_Mesh *s_mesh;

  /* Initializations */

  error = 1;
  model = st_get_model();
  ndata = dbin->getActiveSampleNumber();
  nvar = model->getVariableNumber();
  BheteroD = BheteroT = nullptr;
  Btriplet = nullptr;
  ranks = ndata1 = ntarget1 = nullptr;
  SPDE_Matelem &Mat1 = spde_get_current_matelem(0);
  s_mesh = Mat1.s_mesh;
  nvertex = s_mesh->nvertex;

  /* Core allocation */

  ranks = st_get_vertex_ranks(s_mesh, dbin, dbout);
  if (ranks == nullptr) goto label_end;

  /* Define the sparse matrices */

  ndata1 = (int*) mem_alloc(sizeof(int) * nvar, 0);
  if (ndata1 == nullptr) goto label_end;
  for (int ivar = 0; ivar < nvar; ivar++)
    ndata1[ivar] = 0;
  ntarget1 = (int*) mem_alloc(sizeof(int) * nvar, 0);
  if (ntarget1 == nullptr) goto label_end;
  for (int ivar = 0; ivar < nvar; ivar++)
    ntarget1[ivar] = 0;
  BheteroD = (cs**) mem_alloc(sizeof(cs*) * nvar, 0);
  if (BheteroD == nullptr) goto label_end;
  for (int ivar = 0; ivar < nvar; ivar++)
    BheteroD[ivar] = nullptr;
  BheteroT = (cs**) mem_alloc(sizeof(cs*) * nvar, 0);
  if (BheteroT == nullptr) goto label_end;
  for (int ivar = 0; ivar < nvar; ivar++)
    BheteroT[ivar] = nullptr;

  /**************************************************************/
  /* Creating the sparse matrix for handling heterotopy on data */
  /**************************************************************/
  /* This matrix dimension is [NDmax , Nvertex]                 */
  /* where NDmax is the number of active samples of Dbin        */
  /* regardless of their contents (heterotopy)                  */
  /* A line (dbin sample) contains the barycenter weights       */
  /* assigned to the vertices of the mesh to which it belongs   */
  /* (if the variable is defined for this sample); 0 otherwise  */

  for (int ivar = 0; ivar < nvar; ivar++)
  {
    Btriplet = cs_spalloc(0, 0, 1, 1, 1);
    if (Btriplet == nullptr) goto label_end;

    for (int i = 0; i < nvertex; i++)
    {
      if (ranks[i] <= 0) continue; // Target or Steiner
      ndata1[ivar]++;
      iech = ranks[i] - 1;
      value = (FFFF(dbin->getVariable(iech, ivar))) ? 0. :
                                                      1.;
      if (!cs_entry(Btriplet, iech, i, value)) goto label_end;
    }
    // Add a fictitious sample (zero value) as a dimension constraint
    if (!cs_entry(Btriplet, ndata1[ivar] - 1, nvertex - 1, 0.)) goto label_end;
    BheteroD[ivar] = cs_triplet(Btriplet);
    Btriplet = cs_spfree(Btriplet);
    st_keypair_cs("HeteroD", BheteroD[ivar], st_get_current_icov() + 1,
                  ivar + 1, 0, 0, 0);
  }

  /* Optional printout */

  if (VERBOSE)
    message("Calculation of Bhetero for Data (%d sparse matrices)\n", nvar);

  /**********************************************************************/
  /* Creating the sparse matrix for handling heterotopy on target       */
  /**********************************************************************/
  /* Per variable, this matrix should mask off all mesh vertices which  */
  /* correspond to a data sample which is defined for this variable but */
  /* not defined for all variables                                      */

  for (int ivar = 0; ivar < nvar; ivar++)
  {
    Btriplet = cs_spalloc(0, 0, 1, 1, 1);
    if (Btriplet == nullptr) goto label_end;

    ecrT = 0;
    for (int i = 0; i < nvertex; i++)
    {
      flag_add = 0;
      if (ranks[i] <= 0)
      {
        // This is a target (from dbout) or a Steiner
        flag_add = 1;
      }
      else
      {
        // This is a data
        iech = ranks[i] - 1;

        // Could the data be considered as a target (heterotopic case)
        if (FFFF(dbin->getVariable(iech, ivar)))
        {
          // The sample is not defined for the current variable: it is a target
          flag_add = 1;
        }
      }
      if (!flag_add) continue;
      ntarget1[ivar]++;
      if (!cs_entry(Btriplet, ecrT, i, 1.)) goto label_end;
      ecrT++;
    }

    // Add a fictitious sample (zero value) as a dimension constraint
    if (!cs_entry(Btriplet, ntarget1[ivar] - 1, nvertex - 1, 0.))
      goto label_end;
    BheteroT[ivar] = cs_triplet(Btriplet);
    Btriplet = cs_spfree(Btriplet);
    st_keypair_cs("HeteroT", BheteroT[ivar], st_get_current_icov() + 1,
                  ivar + 1, 0, 0, 0);
  }

  /* Optional printout */

  if (VERBOSE)
    message("Calculation of Bhetero for Target (%d sparse matrices)\n", nvar);

  /* Set the error return code */

  error = 0;

  label_end: st_get_current_ssenv()->BheteroD = BheteroD;
  st_get_current_ssenv()->BheteroT = BheteroT;
  st_get_current_ssenv()->ndata = ndata;
  st_get_current_ssenv()->ndata1 = ndata1;
  st_get_current_ssenv()->ntarget1 = ntarget1;
  ranks = (int*) mem_free((char* ) ranks);
  if (error) st_clean_Bhetero();
  return (error);
}

/****************************************************************************/
/*!
 **  Get the 3-D coordinates of the center of a triangle on the sphere
 **
 ** \param[in]  amesh    MeshEStandard structure
 ** \param[in]  ncorner  Number of vertices per element
 ** \param[in]  imesh    Rank of the current mesh
 **
 ** \param[out] center   Coordinates of the center point (Dimension: 3)
 ** \param[out] xyz      Coordinate of the point (Dimension: 3x3)
 **
 *****************************************************************************/
static void st_triangle_center(MeshEStandard *amesh,
                               int ncorner,
                               int imesh,
                               double center[3],
                               double xyz[3][3])
{
  double ratio;

  for (int i = 0; i < 3; i++)
    center[i] = 0.;
  for (int icorn = 0; icorn < ncorner; icorn++)
  {
    util_convert_sph2cart(amesh->getCoor(imesh, icorn, 0),
                          amesh->getCoor(imesh, icorn, 1), &xyz[icorn][0],
                          &xyz[icorn][1], &xyz[icorn][2]);
    for (int i = 0; i < 3; i++)
      center[i] += xyz[icorn][i];
  }

  ratio = 0.;
  for (int i = 0; i < 3; i++)
  {
    center[i] /= 3.;
    ratio += center[i] * center[i];
  }
  ratio = 1. / sqrt(ratio);
  for (int i = 0; i < 3; i++)
    center[i] *= ratio;
}

/****************************************************************************/
/*!
 **  Project a point on the tangent plane
 **
 ** \param[in]  center    Coordinates of the reference point (Dimension: 3)
 ** \param[in]  axes      Coordinates of the endpoints (Dimension: 2 * 3)
 ** \param[in]  xyz       Coordinates of the target point (Dimension: 3)
 **
 ** \param[out] coeff     Coordinate of point in the local system (Dimension: 2)
 **
 *****************************************************************************/
static void st_project_plane(double center[3],
                             double axes[2][3],
                             double xyz[3],
                             double coeff[2])
{
  double v[3];

  /* Projection */

  for (int j = 0; j < 2; j++)
  {
    coeff[j] = 0.;
    for (int i = 0; i < 3; i++)
      coeff[j] += (axes[j][i] - center[i]) * (xyz[i] - center[i]);
  }

  /* Projected vector */

  for (int i = 0; i < 3; i++)
  {
    v[i] = 0.;
    for (int j = 0; j < 2; j++)
      v[i] += coeff[j] * (axes[j][i] - center[i]);
  }

  /* Returned coordinates */

  vector_translate(3, center, v, xyz);
}

/****************************************************************************/
/*!
 **  Get the coordinates of the axis endpoints in the tangent plane
 **
 ** \param[in]  center  Coordinates of the reference point (Dimension: 3)
 ** \param[in]  srot    Rotation angles on sphere (Dimension: 2)
 **
 ** \param[out] axes    Coordinates of the endpoints (Dimension: 2 * 3)
 **
 *****************************************************************************/
static void st_tangent_calculate(double center[3],
                                 double srot[2],
                                 double axes[2][3])
{
  double sinphi, cosphi, sintet, costet, theta, phi, v[3], w[3];

  // Center gives the vector joining the origin to the center of triangle
  phi = srot[1] * GV_PI / 180.;
  theta = srot[0] * GV_PI / 180.;
  sinphi = sin(phi);
  cosphi = cos(phi);
  sintet = sin(theta);
  costet = cos(theta);
  // W is the Pole vector
  w[0] = sinphi * costet;
  w[1] = sinphi * sintet;
  w[2] = cosphi;
  // V = Center ^ w: first axis
  vector_product(center, w, v);
  ut_normalize(3, v);
  // W = Center ^ V: second axis
  vector_product(center, v, w);
  ut_normalize(3, w);
  // Get the end points from Unit vectors
  vector_translate(3, center, v, axes[0]);
  vector_translate(3, center, w, axes[1]);
}

/****************************************************************************/
/*!
 **  Fill the sparse matrix S linked to mesh vertices
 **
 ** \return G sparse matrix
 **
 ** \param[in]  amesh     MeshEStandard_Mesh structure
 ** \param[in]  model     Model structure
 ** \param[in]  units     Array containing the mesh dimensions
 **
 *****************************************************************************/
cs* _spde_fill_S(MeshEStandard *amesh, Model *model, double *units)
{
  double vald, mat[16], matu[16], matw[16], matinvw[16], mat1[16];
  double xyz[3][3], center[3], axes[2][3], matv[3], coeff[3][2];
  int ecr, errcod, error, ndim, ncorner, flag_sphere, flag_nostat;
  long ip1, ip2;
  cs *Gtriplet, *G;
  std::map<std::pair<int, int>, double> tab;
  std::pair<std::map<std::pair<int, int>, double>::iterator, bool> ret;
  std::map<std::pair<int, int>, double>::iterator it;

  /* Initializations */

  error = 1;
  Gtriplet = G = nullptr;
  ndim = amesh->getNDim();
  ncorner = amesh->getNApexPerMesh();
  Gtriplet = cs_spalloc(0, 0, 1, 1, 1);
  model = st_get_model();
  if (Gtriplet == nullptr) goto label_end;
  variety_query(&flag_sphere);
  flag_nostat = model->isNoStat();
  if (!flag_nostat) st_calcul_update();

  /* Loop on the meshes */

  for (int imesh = 0; imesh < amesh->getNMeshes(); imesh++)
  {

    /* Get parameters in the non-stationary case */

    if (flag_nostat) st_calcul_update_nostat(amesh, imesh);

    // Processing on the Sphere

    if (flag_sphere)
    {

      // Case of the calculations on the Sphere

      st_triangle_center(amesh, ncorner, imesh, center, xyz);
      if (ncorner < 0 || ncorner > 3)
      {
        messerr("Error in st_triangle_center: wrong number or corners: %d",
                ncorner);
        goto label_end;
      }

      /* Look for the tangent plane and its axes */

      st_tangent_calculate(center, Calcul.srot.data(), axes);

      /* Project corner points on the Tangent plane */

      for (int icorn = 0; icorn < ncorner; icorn++)
        st_project_plane(center, axes, xyz[icorn], coeff[icorn]);

      for (int icorn = 0; icorn < ncorner; icorn++)
      {
        for (int idim = 0; idim < ndim; idim++)
          MATU(idim,icorn) = coeff[icorn][idim];
        MATU(ncorner-1,icorn) = 1.;
      }
    }
    else
    {

      // Case of Euclidean geometry

      for (int icorn = 0; icorn < ncorner; icorn++)
      {
        for (int idim = 0; idim < ndim; idim++)
          MATU(idim,icorn) = amesh->getCoor(imesh, icorn, idim);
        MATU(ncorner-1,icorn) = 1.;
      }
    }

    /* Invert the matrix 'matu'*/

    errcod = matrix_invreal(matu, ncorner);
    if (errcod)
    {
      messerr("Error in Mesh #%3d - Its volume is zero", imesh + 1);
      for (int icorn = 0; icorn < ncorner; icorn++)
      {
        message("Sample #%4d - Coordinates (", amesh->getApex(imesh, icorn));
        for (int idim = 0; idim < ndim; idim++)
          message(" %lf", amesh->getCoor(imesh, icorn, idim));
        message(")\n");
      }
      print_matrix("MATU", 0, 1, ncorner, ncorner, NULL, matu);
      goto label_end;
    }

    ecr = 0;
    for (int icorn = 0; icorn < ncorner; icorn++)
      for (int idim = 0; idim < ndim; idim++)
        matw[ecr++] = MATU(icorn, idim);
    matrix_transpose(ndim, ncorner, matw, matinvw);

    matrix_product(ncorner, ndim, ndim, matinvw, Calcul.hh.data(), mat1);
    if (flag_nostat)
      matrix_product(ncorner, ndim, 1, matinvw, Calcul.vv.data(), matv);
    matrix_product(ncorner, ndim, ncorner, mat1, matw, mat);

    for (int j0 = 0; j0 < ncorner; j0++)
      for (int j1 = 0; j1 < ncorner; j1++)
      {
        ip1 = amesh->getApex(imesh, j0);
        ip2 = amesh->getApex(imesh, j1);
        std::pair<int, int> key(ip1, ip2);
        vald = units[imesh] * MAT(j0, j1);
        if (flag_nostat) vald += matv[j1] * units[imesh];
        ret = tab.insert(std::pair<std::pair<int, int>, double>(key, vald));
        if (!ret.second) ret.first->second += vald;
      }
  }

  it = tab.begin();
  while (it != tab.end())
  {
    ip1 = it->first.first;
    ip2 = it->first.second;
    if (!cs_entry(Gtriplet, ip1, ip2, it->second)) goto label_end;
    it++;
  }

  /* Optional printout */

  G = cs_triplet(Gtriplet);
  if (G == nullptr) goto label_end;
  if (VERBOSE) message("Filling G Sparse Matrix performed successfully\n");

  /* Set the error return code */

  error = 0;

  label_end: Gtriplet = cs_spfree(Gtriplet);
  if (error) G = cs_spfree(G);
  return (G);
}

/****************************************************************************/
/*!
 **  Fill the vector TildeC (Dimension: nvertex)
 **
 ** \return Error return code
 **
 ** \param[in]  amesh     MeshEStandard structure
 ** \param[in]  units     Array containing the element units
 **
 *****************************************************************************/
VectorDouble _spde_fill_TildeC(MeshEStandard *amesh, double *units)
{
  VectorDouble tildec, cumunit;
  int nvertex = amesh->getNApices();
  int ncorner = amesh->getNApexPerMesh();
  cumunit.resize(nvertex, 0);

  /* Loop on the meshes */

  for (int imesh = 0; imesh < amesh->getNMeshes(); imesh++)
  {

    /* Loop on the vertices */

    for (int icorn = 0; icorn < ncorner; icorn++)
    {
      int ip = amesh->getApex(imesh, icorn);
      cumunit[ip] += units[imesh];
    }
  }

  /* Scale */

  double factor = (double) ncorner;
  for (int ip = 0; ip < nvertex; ip++)
  {
    double value = cumunit[ip] / factor;
    if (ABS(value) <= 0.)
    {
      messerr("Meshing unit (%d) has a zero volume", ip + 1);
      return VectorDouble();
    }
    tildec.push_back(value);
  }
  return tildec;
}

/****************************************************************************/
/*!
 **  Fill the vector for sill correction factors
 **  Works for both stationary and non-stationary cases
 **
 ** \param[in]  model     Model structure
 ** \param[in]  amesh     MeshEStandard structure
 ** \param[in]  TildeC    Vector TildeC
 **
 *****************************************************************************/
VectorDouble _spde_fill_Lambda(Model *model,
                               MeshEStandard *amesh,
                               const VectorDouble &TildeC)
{
  const ANoStat *nostat = model->getNoStat();
  VectorDouble Lambda;
  int igrf0 = st_get_current_igrf();
  int icov0 = st_get_current_icov();
  int ndim = st_get_ndim();
  int nvertex = amesh->getNApices();
  double sill = st_get_cova_sill(0, 0);

  /* Fill the array */

  if (st_get_model()->isNoStat() && nostat->isDefinedforAnisotropy(igrf0,
                                                                   icov0))
  {
    for (int ip = 0; ip < nvertex; ip++)
    {
      nostat->updateModel(model, ip);
      st_compute_hh();
      double sqdeth = sqrt(matrix_determinant(ndim, Calcul.hh.data()));
      Lambda.push_back(sqrt((TildeC[ip]) / (sqdeth * sill)));
    }
  }
  else
  {
    double sqdeth = Calcul.sqdeth;
    for (int ip = 0; ip < nvertex; ip++)
      Lambda.push_back(sqrt((TildeC[ip]) / (sqdeth * sill)));
  }

  return (Lambda);
}

/****************************************************************************/
/*!
 **  Extract the sparse matrix from the Q matrix (case of nugget effect)
 **
 ** \return The extracted sparse matrix or NULL
 **
 ** \param[in]  row_var    Rank of the variable for the row
 ** \param[in]  col_var    Rank of the variable for the column
 **
 ** \param[out] nrows      Number of rows
 ** \param[out] ncols      Number of columns
 **
 ** \remarks Extracts a part of Bnugget matrix for:
 ** \remarks - a given pair of variables
 ** \remarks - for Data-Data operators
 **
 *****************************************************************************/
static cs* st_extract_Q1_nugget(int row_var,
                                int col_var,
                                int *nrows,
                                int *ncols)
{
  SPDE_SS_Environ *SS;
  cs *B0;

  SS = st_get_current_ssenv();
  B0 = cs_duplicate(SS->Bnugget[st_get_rank(row_var, col_var)]);
  if (B0 != nullptr) *nrows = *ncols = B0->n;

  /* Keypair storage (optional) */

  st_keypair_cs("ExtractNug", B0, st_get_current_icov() + 1, row_var + 1,
                col_var + 1, 0, 0);

  return (B0);
}

/****************************************************************************/
/*!
 **  Extract the sparse matrix from the Q matrix (coninuous structure)
 **
 ** \return The extracted sparse matrix or NULL
 **
 ** \param[in]  row_var    Rank of the variable for the row
 ** \param[in]  col_var    Rank of the variable for the column
 ** \param[in]  row_oper   Operator type for row (1:Data or 2:Target)
 ** \param[in]  col_oper   Operator type for column (1:Data or 2:Target)
 **
 ** \param[out] nrows      Number of rows
 ** \param[out] ncols      Number of columns
 **
 ** \remarks Extracts a part of Q matrix (for the first structure) for:
 ** \remarks - a given pair of variables
 ** \remarks - a given pair of operators (Data or Target)
 ** \remarks The returned matrix is multipled by the inverse of the Sill
 **
 *****************************************************************************/
static cs* st_extract_Q1_hetero(int row_var,
                                int col_var,
                                int row_oper,
                                int col_oper,
                                int *nrows,
                                int *ncols)
{
  int error;
  cs *Q, *Brow, *Bcol, *B1, *Bt, *Qn;
  SPDE_SS_Environ *SS;

  /* Initializations */

  error = 1;
  Q = Brow = Bcol = B1 = Bt = Qn = nullptr;
  SS = st_get_current_ssenv();
  SPDE_Matelem &Matelem1 = spde_get_current_matelem(0);

  /* Identify the operating matrices */

  Brow = (row_oper == 1) ? SS->BheteroD[row_var] :
                           SS->BheteroT[row_var];
  if (Brow == nullptr) goto label_end;
  Bcol = (col_oper == 1) ? SS->BheteroD[col_var] :
                           SS->BheteroT[col_var];
  if (Bcol == nullptr) goto label_end;
  Bt = cs_transpose(Bcol, 1);
  if (Bt == nullptr) goto label_end;
  B1 = cs_multiply(Brow, Matelem1.QC->Q);
  if (B1 == nullptr) goto label_end;
  Qn = cs_multiply(B1, Bt);
  if (Qn == nullptr) goto label_end;

  /* Multiply by the corresponding sill */

  Q = cs_add(Qn, Qn, st_get_isill(0, row_var, col_var), 0.);
  if (Q == nullptr) goto label_end;

  /* Set the error return code */

  error = 0;
  *nrows = (row_oper == 1) ? SS->ndata1[row_var] :
                             SS->ntarget1[row_var];
  *ncols = (col_oper == 1) ? SS->ndata1[col_var] :
                             SS->ntarget1[col_var];

  /* Keypair storage (optional) */

  st_keypair_cs("Extract", Q, st_get_current_icov() + 1, row_var + 1,
                col_var + 1, row_oper, col_oper);

  label_end: B1 = cs_spfree(B1);
  Bt = cs_spfree(Bt);
  Qn = cs_spfree(Qn);
  if (error) Q = cs_spfree(Q);
  return (Q);
}

/****************************************************************************/
/*!
 **  Construct the sparse matrix QCov (used in multistructure - multivariable)
 **
 ** \return Error return code
 **
 ** \param[in]  Matelem     SPDE_Matelem structure
 **
 ** \remarks This function requires the Q matrices to be established already,
 ** \remarks as well as the Aproj matrices.
 ** \remarks In case of presence of nugget effect, we also need 'Bnugget'
 ** \remarks Otherwise, we need 'BheteroD' and 'BheteroT'
 **
 *****************************************************************************/
static int st_build_QCov(SPDE_Matelem &Matelem)

{
  int error, nvar, icov0, nrows, ncols;
  cs *B0, *Bi;
  QChol **QCov;
  SPDE_SS_Environ *SS;

  /* Initializations */

  if (!S_DECIDE.flag_several) return (0);
  error = 1;
  nvar = st_get_nvar();
  Bi = B0 = nullptr;
  SS = st_get_current_ssenv();
  icov0 = st_get_current_icov();

  /* Core allocation */

  QCov = (QChol**) mem_alloc(sizeof(QChol*) * nvar, 1);
  for (int ivar = 0; ivar < nvar; ivar++)
    QCov[ivar] = qchol_manage(1, NULL);

  /* Dispatch */

  if (st_is_model_nugget())
  {

    /****************************************/
    /* Case when a nugget effect is present */
    /****************************************/

    if (Matelem.Aproj == NULL || SS->Bnugget == NULL) return (1);

    for (int ivar = 0; ivar < nvar; ivar++)
    {
      // Sill(icov)_ii * Q(icov) + A^t(icov) * E_ii * A(icov)
      B0 = st_extract_Q1_nugget(ivar, ivar, &nrows, &ncols);
      if (B0 == nullptr) goto label_end;
      Bi = cs_prod_norm(1, B0, Matelem.Aproj);
      if (Bi == nullptr) goto label_end;
      QCov[ivar]->Q = cs_add(Matelem.QC->Q, Bi, st_get_isill(icov0, ivar, ivar),
                             1.);
      if (QCov[ivar]->Q == nullptr) goto label_end;
      Bi = cs_spfree(Bi);
      B0 = cs_spfree(B0);
    }
  }
  else
  {

    /***************************************/
    /* Case when there is no nugget effect */
    /***************************************/

    if (Matelem.Aproj == NULL || SS->BheteroD == NULL || SS->BheteroT == NULL)
      return (1);

    for (int ivar = 0; ivar < nvar; ivar++)
    {
      if (icov0 == 0)
      {
        // Q1_tt_ii
        QCov[ivar]->Q = st_extract_Q1_hetero(ivar, ivar, 2, 2, &nrows, &ncols);
        if (QCov[ivar]->Q == nullptr) goto label_end;
      }
      else
      {
        // Sill(icov)_ii * Q(icov) + A^t(icov) * Q1_dd_ii * A(icov)
        B0 = st_extract_Q1_hetero(ivar, ivar, 1, 1, &nrows, &ncols);
        if (B0 == nullptr) goto label_end;
        Bi = cs_prod_norm(1, B0, Matelem.Aproj);
        if (Bi == nullptr) goto label_end;
        QCov[ivar]->Q = cs_add(Matelem.QC->Q, Bi,
                               st_get_isill(icov0, ivar, ivar), 1.);
        if (QCov[ivar]->Q == nullptr) goto label_end;
        Bi = cs_spfree(Bi);
        B0 = cs_spfree(B0);
      }

      st_keypair_cs("QCov", QCov[ivar]->Q, icov0 + 1, 0, 0, 0, 0);
    }
  }
  Matelem.QCov = QCov;

  /* Optional printout */

  if (VERBOSE) message("Building QCov (%d sparse matrices)\n", nvar);

  /* Set the error return code */

  error = 0;

  label_end: B0 = cs_spfree(B0);
  Bi = cs_spfree(Bi);
  if (error)
  {
    if (QCov != NULL)
    {
      for (int ivar = 0; ivar < nvar; ivar++)
        QCov[ivar] = qchol_manage(-1, QCov[ivar]);
    }
  }
  return (error);
}

/****************************************************************************/
/*!
 **  Construct the final sparse matrix Q from the Model
 **
 ** \return Error return code
 **
 ** \param[in] S        Shift operator
 ** \param[in] Lambda   Lambda vector
 ** \param[in] nblin    Number of blin coefficients
 ** \param[in] blin     Array of coefficients for Linear combinaison
 **
 *****************************************************************************/
cs* _spde_build_Q(cs *S, const VectorDouble &Lambda, int nblin, double *blin)
{
  int error, iterm, nvertex;
  double *work, *tblin;
  cs *Be, *Q, *Bi;

  /* Initializations */

  error = 1;
  iterm = 0;
  nvertex = S->n;
  work = tblin = nullptr;
  Q = Be = Bi = nullptr;

  // Preliminary checks

  if (nvertex <= 0)
  {
    messerr("You must define a valid Meshing beforehand");
    return Q;
  }
  if (nblin <= 0)
  {
    messerr("You must have a set of already available 'blin' coefficients");
    messerr("These coefficients come from the decomposition in series for Q");
    messerr("This decomposition is available only if 'alpha' is an integer");
    messerr("where: alpha = param + ndim/2");
    return Q;
  }
  /* Build the tblin array */

  tblin = (double*) mem_alloc(sizeof(double) * nvertex * nblin, 0);
  if (tblin == nullptr) goto label_end;
  for (int i = 0; i < nvertex * nblin; i++)
    tblin[i] = 0;

  /* Stationary Case */

  for (int ib = 0; ib < nblin; ib++)
    for (int ip = 0; ip < nvertex; ip++)
      tblin[nvertex * (ib) + (ip)] = sqrt(blin[ib]);

  /* First step */

  work = (double*) mem_alloc(sizeof(double) * nvertex, 0);
  if (work == nullptr) goto label_end;
  for (int i = 0; i < nvertex; i++)
    work[i] = TBLIN(0,i) * TBLIN(0, i);
  Q = cs_eye_tab(nvertex, work);
  if (Q == nullptr) goto label_end;
  work = (double*) mem_free((char* ) work);
  Bi = cs_duplicate(S);
  if (Bi == nullptr) goto label_end;

  /* Loop on the different terms */

  for (iterm = 1; iterm < nblin; iterm++)
  {
    Be = cs_matvecnorm(Bi, &TBLIN(iterm, 0), 0);
    if (Be == nullptr) goto label_end;
    Q = cs_add_and_release(Q, Be, 1., 1., 1);
    if (Q == nullptr) goto label_end;
    Be = cs_spfree(Be);

    if (iterm < nblin - 1)
    {
      Bi = cs_multiply_and_release(Bi, S, 1);
      if (Bi == nullptr) goto label_end;
    }
  }
  Bi = cs_spfree(Bi);

  /* Final scaling */

  cs_matvecnorm_inplace(Q, Lambda.data(), 0);

  /* Set the error return code */

  error = 0;

  label_end: if (error) Q = cs_spfree(Q);

  Bi = cs_spfree(Bi);
  work = (double*) mem_free((char* ) work);
  tblin = (double*) mem_free((char* ) tblin);

  return Q;
}

/****************************************************************************/
/*!
 **  Construct the final sparse matrix Q from the Model
 **
 ** \return Error return code
 **
 ** \param[in]  Matelem    SPDE_Matelem structure
 **
 *****************************************************************************/
static int st_build_Q(SPDE_Matelem &Matelem)

{
  int error;
  QChol *QC;

  /* Initializations */

  error = 1;
  QC = nullptr;

  /* Core allocation */

  Matelem.QC = qchol_manage(1, NULL);

  int nblin = static_cast<int>(Calcul.blin.size());
  Matelem.QC->Q = _spde_build_Q(Matelem.S, Matelem.Lambda, nblin,
                                Calcul.blin.data());
  if (Matelem.QC->Q == nullptr) goto label_end;

  /* Optional printout */

  if (VERBOSE) message("Building Global Q matrix\n");

  /* Set the error return code */

  error = 0;

  label_end: if (error) QC = qchol_manage(-1, QC);
  return (error);
}

/****************************************************************************/
/*!
 **  Build all matrices needed for establishing the Q sparse matrix
 **
 ** \return Error return code
 **
 ** \param[in]  model      Model structure
 ** \param[in]  verbose    Verbose option
 **
 ** \remarks Contents of SP_MAT (sparse matrices or vectors) is allocated here
 ** \remarks It must be freed by the calling functions
 **
 *****************************************************************************/
int spde_build_matrices(Model *model, int verbose)
{
  int error;
  double *units;
  VectorDouble tildec;

  /* Initializations */

  error = 1;
  units = nullptr;
  VERBOSE = verbose;
  SPDE_Matelem &Matelem = spde_get_current_matelem(-1);
  MeshEStandard amesh;
  amesh.convertFromOldMesh(Matelem.s_mesh, 0);

  /* Calculate the units of the meshes */

  units = _spde_get_mesh_dimension(&amesh);
  if (units == nullptr) goto label_end;

  /* Fill S sparse matrix */

  Matelem.S = _spde_fill_S(&amesh, model, units);
  if (Matelem.S == nullptr) goto label_end;
  if (VERBOSE) message("Filling S Sparse Matrix performed successfully\n");

  /* Fill the TildeC vector */

  tildec = _spde_fill_TildeC(&amesh, units);
  if (VERBOSE) message("Filling TildeC Sparse Matrix performed successfully\n");

  /* Construct the matrix for the sill correction array */

  Matelem.Lambda = _spde_fill_Lambda(model, &amesh, tildec);
  if (VERBOSE) message("Filling Lambda Sparse Matrix performed successfully\n");

  /* Build the sparse matrix B */

  cs_matvecnorm_inplace(Matelem.S, tildec.data(), 2);

  /* Build the sparse matrix Q */

  if (S_DECIDE.flag_Q)
  {
    if (st_build_Q(Matelem)) goto label_end;
  }

  /* Set the error return code */

  error = 0;

  label_end: units = (double*) mem_free((char* ) units);
  return (error);
}

/****************************************************************************/
/*!
 **  Load the data information within the complete output array
 **
 ** \param[in]  s_mesh      SPDE_Mesh structure
 ** \param[in]  dbin        Db input structure
 ** \param[in]  dbout       Db output structure
 ** \param[in]  s_option    SPDE_Option structure
 ** \param[in]  ivar0       Rank of the variable or -1
 **
 ** \param[out] data        Vector of active data
 ** \param[out] zcur        Output array (Dimension: number of meshes)
 **
 ** \remarks When 'ivar0' is defined, this corresponds to the flag_mult_data
 ** \remarks flag where the target is the simulation
 ** \remarks Otherwise the Z-variables must all be loaded
 **
 *****************************************************************************/
static void st_load_data(SPDE_Mesh *s_mesh,
                         Db *dbin,
                         Db *dbout,
                         SPDE_Option &s_option,
                         int ivar0,
                         double *data,
                         double *zcur)
{
  Vertype *vertype;
  double zloc;
  int ecr, ijoint, ecrd, nvar, ivar, nvertex, igrf;

  /* Initializations */

  vertype = s_mesh->vertype;
  nvertex = vertype->nvertex;
  igrf = st_get_current_igrf();
  nvar = (ivar0 >= 0) ? 1 :
                        st_get_nvar();
  MEM_DBIN = dbin;           // This horrible assignment is to allows
  MEM_DBOUT = dbout;          // passing information to Kriging sub-procedures

  /* Loop on the variables */

  ecrd = ecr = 0;
  for (int jvar = 0; jvar < nvar; jvar++)
  {
    ecr = jvar * nvertex;
    ivar = (ivar0 >= 0) ? ivar0 :
                          jvar;

    /* From the Output Db if present */

    if (dbout != nullptr)
    {
      for (int iech = 0; iech < dbout->getSampleNumber(); iech++)
      {
        if (!dbout->isActive(iech)) continue;
        ijoint = st_is_duplicated(s_mesh->vercoloc, 2, iech);
        if (ijoint < 0)
        {

          /* This target is not collocated to any Datum */

          if (S_DECIDE.flag_gibbs && dbout->getIntervalNumber() > 0)
            zloc = st_get_data_constraints(dbout, igrf, iech);
          else
            zloc = TEST;
        }
        else
        {

          /* This target is collocated to a Datum */

          if (S_DECIDE.flag_gibbs && dbin->getIntervalNumber() > 0)
            zloc = st_get_data_constraints(dbin, igrf, ijoint);
          else
          {
            zloc = dbin->getVariable(ijoint, ivar);
            if (!S_DECIDE.flag_several) data[ecrd++] = zloc;
          }
        }
        if (S_DECIDE.flag_mesh_dbout)
        {
          if (!FFFF(zloc))
          {
            zcur[ecr] = zloc;
            if (st_get_filnug()) zcur[ecr] /= st_get_nugget_sill(0, 0);
          }
          else
          {
            if (st_get_filnug()) zcur[ecr] = 0.;
          }
          ecr++;
        }
      }
    }

    /* From the Input Db (if present) */
    /* Note that 'ecr' is only incremented for non-duplicate sample */

    if (dbin != nullptr)
    {
      for (int iech = 0; iech < dbin->getSampleNumber(); iech++)
      {
        if (!dbin->isActive(iech)) continue;
        if (S_DECIDE.flag_several) data[ecrd++] = dbin->getVariable(iech, ivar);

        /* If collocated with a target, vertex has already been loaded */

        if (st_is_duplicated(s_mesh->vercoloc, 1, iech) >= 0) continue;

        /* No collocation */

        if (S_DECIDE.flag_gibbs && dbin->getIntervalNumber() > 0)
          zloc = st_get_data_constraints(dbin, igrf, iech);
        else
          zloc = dbin->getVariable(iech, ivar);

        if (!S_DECIDE.flag_several) data[ecrd++] = zloc;

        if (S_DECIDE.flag_mesh_dbin)
        {
          if (!FFFF(zloc))
          {
            zcur[ecr] = zloc;
            if (st_get_filnug()) zcur[ecr] /= st_get_nugget_sill(0, 0);
          }
          else
          {
            if (st_get_filnug()) zcur[ecr] = 0.;
          }
          ecr++;
        }
      }
    }
  }

  if (DEBUG)
  {
    message("(DEBUG) Load Data\n");
    for (ivar = 0; ivar < nvar; ivar++)
      print_range("- Load  ", nvertex, &zcur[ivar * nvertex], NULL);
  }
}

/****************************************************************************/
/*!
 **  Internal function used for the Chebychev approximation
 **
 ** \return  Returned value
 **
 ** \param[in]  x           Input value
 ** \param[in]  power       Parameter used in the Chebychev approximation
 ** \param[in]  nblin       Number of blin coefficients
 ** \param[in]  blin        Array of coefficients for Linear combination
 **
 *****************************************************************************/
static double st_chebychev_function(double x,
                                    double power,
                                    int nblin,
                                    double *blin)
{
  double value, total;

  value = 1.;
  total = blin[0];
  for (int i = 1; i < nblin; i++)
  {
    value *= x;
    total += blin[i] * value;
  }
  if (power == 0.)
    return (log(total));
  else
    return (pow(total, power));
}

/****************************************************************************/
/*!
 **  Evaluate the number of coefficients necessary to evaluate a function
 **  (at a sample location) at a given approximation
 **
 ** \return Error return code
 **
 ** \param[in]  cheb_elem  Cheb_Elem structure to be filled
 ** \param[in]  verbose    Verbose flag
 ** \param[in]  nblin      Number of blin coefficients
 ** \param[in]  blin       Array of coefficients for Linear combination
 **
 *****************************************************************************/
static int st_chebychev_calculate_coeffs(Cheb_Elem *cheb_elem,
                                         int verbose,
                                         int nblin,
                                         double *blin)

{
  double value, a, b;
  int error, number, numloc, ndisc;

  /* Initializations */

  error = 1;
  a = cheb_elem->a;
  b = cheb_elem->b;
  ndisc = cheb_elem->ndisc;

  /* Calculate the polynomials */

  cheb_elem->coeffs = (double*) mem_alloc(sizeof(double) * cheb_elem->ncmax, 0);
  if (cheb_elem->coeffs == nullptr) goto label_end;

  /* Evaluate the coefficients of the Chebychev approximation */

  if (ut_chebychev_coeffs(st_chebychev_function, cheb_elem, nblin, blin))
    goto label_end;

  /* Loop on some discretized samples of the interval */

  number = 0;
  for (int idisc = 1; idisc < ndisc; idisc++)
  {
    value = a + (b - a) * idisc / ndisc;
    numloc = ut_chebychev_count(st_chebychev_function, cheb_elem, value, nblin,
                                blin);
    if (numloc > number) number = numloc;
  }

  /* Optional printout */

  if (verbose)
  {
    message("Chebychev Polynomial Approximation:\n");
    message("- Power = %lf\n", cheb_elem->power);
    message("- Performed using %d terms\n", number);
    message("- between %lf and %lf (Nb. discretization steps=%d)\n", a, b,
            ndisc);
    message("- with a tolerance of %lg\n", cheb_elem->tol);
  }

  /* Core Reallocation */

  cheb_elem->coeffs = (double*) mem_realloc((char* ) cheb_elem->coeffs,
                                            sizeof(double) * number, 0);
  if (cheb_elem->coeffs == nullptr) goto label_end;
  cheb_elem->ncoeffs = number;

  /* Set the error return code */

  error = 0;

  label_end: return (error);
}

/****************************************************************************/
/*!
 **  Simulate the nugget effect component
 **
 ** \param[in]  ncur       Number of target to be simulated
 ** \param[out] zsnc       Output array (Dimension: nvertex)
 **
 *****************************************************************************/
static void st_simulate_nugget(int ncur, double *zsnc)

{
  double sill;

  sill = st_get_nugget_sill(0, 0);
  if (FFFF(sill) || sill <= 0.) return;
  sill = sqrt(sill);

  for (int icur = 0; icur < ncur; icur++)
    zsnc[icur] += sill * law_gaussian();
}

/****************************************************************************/
/*!
 **  Perform the basic non-conditional Simulation
 **  using the Cholesky decomposition method
 **
 ** \return Error return code
 **
 ** \param[in]  QC         Pointer to the QChol structure (finalized)
 **
 ** \param[out] work       Working array (Dimension: nvertex)
 ** \param[out] zsnc       Output array (Dimension: nvertex)
 **
 *****************************************************************************/
static int st_simulate_cholesky(QChol *QC, double *work, double *zsnc)
{
  int nvertex;

  /* Initializations */

  nvertex = QC->Q->n;
  for (int ip = 0; ip < nvertex; ip++)
    work[ip] = law_gaussian();

  /* Prepare Cholesky decomposition (if not already performed) */

  if (QC->S == nullptr)
  {
    if (qchol_cholesky(VERBOSE, QC)) return (1);
  }

  /* Perform the simulation */

  cs_chol_simulate(QC, zsnc, work);

  if (DEBUG)
  {
    message("(DEBUG) Simulate (Cholesky)\n");
    print_range("- Result", nvertex, zsnc, NULL);
  }
  return (0);
}

/****************************************************************************/
/*!
 **  Perform the Chebychev polynomial procedure on an input vector
 **
 ** \return Error return code
 **
 ** \param[in]     S          Shift operator
 ** \param[in]     cheb_elem  Cheb_Elem structure
 ** \param[in]     lambda     Scaling vector
 ** \param[in]     x          Input array (Dimension: nvertex)
 **
 ** \param[out]    y          Output array (Dimension: nvertex)
 **
 *****************************************************************************/
int spde_chebychev_operate(cs *S,
                           Cheb_Elem *cheb_elem,
                           const VectorDouble &lambda,
                           const double *x,
                           double *y)
{
  double *coeffs, *tm1, *tm2, *px, *tx, v1, v2, coeff_ib, power;
  int error, ncoeffs, nvertex;
  cs *T1;

  /* Initializations */

  error = 1;
  T1 = nullptr;
  tm1 = tm2 = px = tx = nullptr;
  nvertex = S->n;
  v1 = cheb_elem->v1;
  v2 = cheb_elem->v2;
  power = cheb_elem->power;
  ncoeffs = cheb_elem->ncoeffs;
  coeffs = cheb_elem->coeffs;

  /* Core allocation */

  tm1 = (double*) mem_alloc(sizeof(double) * nvertex, 0);
  if (tm1 == nullptr) goto label_end;
  tm2 = (double*) mem_alloc(sizeof(double) * nvertex, 0);
  if (tm2 == nullptr) goto label_end;
  px = (double*) mem_alloc(sizeof(double) * nvertex, 0);
  if (px == nullptr) goto label_end;
  tx = (double*) mem_alloc(sizeof(double) * nvertex, 0);
  if (tx == nullptr) goto label_end;

  /* Create the T1 sparse matrix */

  T1 = cs_eye(nvertex, 1.);
  if (T1 == nullptr) goto label_end;
  T1 = cs_add_and_release(T1, S, v2, v1, 1);
  if (T1 == nullptr) goto label_end;

  /* Initialize the simulation */

  for (int i = 0; i < nvertex; i++)
  {
    tm1[i] = 0.;
    y[i] = x[i];
  }
  if (!cs_gaxpy(T1, y, tm1)) goto label_end;
  for (int i = 0; i < nvertex; i++)
  {
    px[i] = coeffs[0] * y[i] + coeffs[1] * tm1[i];
    tm2[i] = y[i];
  }

  /* Loop on the Chebychev polynomials */

  for (int ib = 2; ib < ncoeffs; ib++)
  {
    cs_mulvec(T1, nvertex, tm1, tx);
    coeff_ib = coeffs[ib];
    for (int i = 0; i < nvertex; i++)
    {
      tx[i] = 2. * tx[i] - tm2[i];
      px[i] += coeff_ib * tx[i];
      tm2[i] = tm1[i];
      tm1[i] = tx[i];
    }
  }

  /* Save the results */

  for (int i = 0; i < nvertex; i++)
    y[i] = px[i] * pow(lambda[i], 2. * power);

  /* Core deallocation */
  tm1 = (double*) mem_free((char* ) tm1);
  tm2 = (double*) mem_free((char* ) tm2);
  px = (double*) mem_free((char* ) px);
  tx = (double*) mem_free((char* ) tx);

  /* Set the error return code */

  error = 0;

  label_end: return (error);
}

/****************************************************************************/
/*!
 **  Perform the basic non-conditional Simulation
 **  using the Chebychev Polynomial procedure
 **
 ** \return Error return code
 **
 ** \param[out] zsnc       Output array (Dimension: nvertex)
 **
 *****************************************************************************/
static int st_simulate_chebychev(double *zsnc)
{
  Cheb_Elem *cheb_elem;
  double *x;
  int error, nvertex;

  /* Initializations */

  error = 1;
  SPDE_Matelem &Matelem = spde_get_current_matelem(-1);
  cheb_elem = Matelem.s_cheb;
  nvertex = Matelem.s_mesh->nvertex;
  x = nullptr;

  /* Core allocation */

  x = (double*) mem_alloc(sizeof(double) * nvertex, 0);
  if (x == nullptr) goto label_end;

  /* Initialize the simulation */

  for (int i = 0; i < nvertex; i++)
    x[i] = law_gaussian();

  /* Operate the Chebychev polynomials */

  if (spde_chebychev_operate(Matelem.S, cheb_elem, Matelem.Lambda, x, zsnc))
    goto label_end;

  /* Set the error return code */

  error = 0;

  label_end: x = (double*) mem_free((char* ) x);
  return (error);
}

/****************************************************************************/
/*!
 **  Perform the Calculation of the Kriging estimate (Multigrid case)
 **
 ** \return  Error return code
 **
 ** \param[in]  QC          Pointer to QChol structure (target-target)(finalized)
 ** \param[in]  rhs         R.H.S. array (Dimension: ntarget)
 **
 ** \param[out] work        Working array (Dimension: ntarget)
 ** \param[out] z           Output array  (Dimension: ntarget)
 **
 *****************************************************************************/
static int st_kriging_multigrid(QChol *QC, double *rhs, double *work, double *z)
{
  int ntarget = QC->Q->n;
  SPDE_Matelem &Matelem = spde_get_current_matelem(-1);

  if (cs_multigrid_process(Matelem.mgs, QC, VERBOSE, z, rhs, work)) return (1);

  if (DEBUG)
  {
    message("(DEBUG) Kriging (Multigrid)\n");
    print_range("- Result", ntarget, z, NULL);
  }
  return (0);
}

/****************************************************************************/
/*!
 **  Solve the linear system (either using Cholesky or Multigrid)
 **
 ** \return  Error return code
 **
 ** \param[in]  QC          Qchol structure
 ** \param[in]  rhs         R.H.S.
 **
 ** \param[out] work        Working array (Dimension: ntarget)
 ** \param[out] z           Output Array (Dimension: ntarget)
 **
 *****************************************************************************/
static int st_solve(QChol *QC, double *rhs, double *work, double *z)
{
  if (S_DECIDE.flag_mgrid)
  {
    if (st_kriging_multigrid(QC, rhs, work, z)) return (1);
  }
  else
  {
    if (st_kriging_cholesky(QC, rhs, work, z)) return (1);
  }
  return (0);
}

/****************************************************************************/
/*!
 **  Perform the Kriging operation (monovariate - monostructure)
 **
 ** \return  Error return code
 **
 ** \param[in]  data        Vector of data
 **
 ** \param[out] rhs         R.H.S.
 ** \param[out] work        Working array
 ** \param[out] z           Output Array
 **
 *****************************************************************************/
static int st_kriging_one(double *data, double *rhs, double *work, double *z)
{
  QSimu *qsimu;
  int ntarget;

  /* Initializations */

  qsimu = spde_get_current_matelem(-1).qsimu;
  ntarget = qsimu->QCtt->Q->n;

  /* Initialize the resulting array */

  st_init_array(1, 1, ntarget, 1, z);

  /* Define the Right-hand side */

  st_kriging_one_rhs(qsimu->QCtd, data, ntarget, rhs);

  /* Solve the Kriging system */

  if (st_solve(qsimu->QCtt, rhs, work, z)) return (1);

  return (0);
}

/****************************************************************************/
/*!
 **  Prepare the R.H.S. for the Kriging operation (several)
 **
 ** \return  Error return code
 **
 ** \param[in]  data        Vector of data
 **
 ** \param[out] rhs         R.H.S.
 ** \param[out] work        Working array
 ** \param[out] ss_arg      Global returned criterion
 **
 *****************************************************************************/
static int st_kriging_several_rhs(double *data,
                                  double *rhs,
                                  double *work,
                                  double *ss_arg)
{
  SPDE_SS_Environ *SS;
  cs *B0;
  double *temp, ss;
  int ndata, ncova, nvar, ncur, error, nrows, ncols, size;

  /* Initializations */

  error = 1;
  SS = st_get_current_ssenv();
  ncova = st_get_ncova();
  nvar = st_get_nvar();
  ncur = st_get_nvertex_max();
  size = st_get_dimension();
  ndata = SS->ndata;
  temp = nullptr;
  B0 = nullptr;

  /* Core allocation */

  temp = (double*) mem_alloc(sizeof(double) * size, 0);
  if (temp == nullptr) goto label_end;

  /* Constitute the RHS */

  for (int icov = 0; icov < ncova; icov++)
  {
    SPDE_Matelem &Matelem = spde_get_current_matelem(icov);
    ncur = st_get_nvertex(icov);
    if (FLAG_KEYPAIR)
    {
      for (int ivar = 0; ivar < nvar; ivar++)
      {
        (void) gslSPrintf(NAME, "DATA.%d", ivar);
        set_keypair(NAME, 1, ndata, 1, &DATA(ivar, 0));
      }
    }

    if (st_is_model_nugget())
    {

      /* Case with nugget effect */

      for (int ivar = 0; ivar < nvar; ivar++)
      {
        for (int i = 0; i < size; i++)
          work[i] = 0.;
        for (int jvar = 0; jvar < nvar; jvar++)
        {
          // E_ij * Z_j
          B0 = st_extract_Q1_nugget(ivar, jvar, &nrows, &ncols);
          if (B0 == nullptr) goto label_end;
          cs_mulvec(B0, nrows, &DATA(jvar, 0), temp);
          for (int i = 0; i < nrows; i++)
            work[i] += temp[i];
          B0 = cs_spfree(B0);
        }
        cs_tmulvec(Matelem.Aproj, ncur, work, &RHS(icov, ivar, 0));
      }
    }
    else
    {

      /* Case without nugget effect */

      for (int ivar = 0; ivar < nvar; ivar++)
      {
        for (int i = 0; i < size; i++)
          work[i] = 0.;
        if (icov == 0)
        {
          for (int jvar = 0; jvar < nvar; jvar++)
          {
            // Q1_td_ij * Z_j
            B0 = st_extract_Q1_hetero(ivar, jvar, 2, 1, &nrows, &ncols);
            if (B0 == nullptr) goto label_end;
            cs_print_dim("Apres extraction de Qtd", B0);
            cs_mulvec(B0, nrows, &DATA(jvar, 0), temp);
            for (int i = 0; i < nrows; i++)
              work[i] += temp[i];
            B0 = cs_spfree(B0);
          }
          // -sum_j { Q1_td_ij * Z_j }
          for (int i = 0; i < ncur; i++)
            RHS(icov,ivar,i) = -work[i];
        }
        else
        {
          for (int jvar = 0; jvar < nvar; jvar++)
          {
            // Q1_dd_ij * Z_j
            B0 = st_extract_Q1_hetero(ivar, jvar, 1, 1, &nrows, &ncols);
            if (B0 == nullptr) goto label_end;
            cs_mulvec(B0, nrows, &DATA(jvar, 0), temp);
            for (int i = 0; i < nrows; i++)
              work[i] += temp[i];
            B0 = cs_spfree(B0);
          }
          // A^t(icov) * sum{ Q1_dd_ij * Z_j } 
          cs_tmulvec(Matelem.Aproj, ncur, work, &RHS(icov, ivar, 0));
        }
      }
    }
  }

  /* Evaluate the reference score */

  ss = 0.;
  for (int icov = 0; icov < ncova; icov++)
  {
    ncur = st_get_nvertex(icov);
    for (int ivar = 0; ivar < nvar; ivar++)
      for (int icur = 0; icur < ncur; icur++)
        ss += RHS(icov,ivar,icur) * RHS(icov, ivar, icur);
  }

  /* Set the error return code */

  error = 0;
  *ss_arg = ss;

  label_end: B0 = cs_spfree(B0);
  temp = (double*) mem_free((char* ) temp);
  return (error);
}

/****************************************************************************/
/*!
 **  Perform the Loop for the Kriging operation (several)
 **
 ** \return  Error return code
 **
 ** \param[in]  flag_crit   1 if used for criterion evaluation;
 **                         0 when used for calculation loop
 **
 ** \param[out] work        Working array
 ** \param[out] rhscur      Working array
 ** \param[out] rhsloc      Working array
 ** \param[in,out] xcur        Current solution vector
 ** \param[in,out] rhs         R.H.S.
 ** \param[in,out] crit        Criterion
 **
 *****************************************************************************/
static int st_kriging_several_loop(int flag_crit,
                                   double *work,
                                   double *rhscur,
                                   double *rhsloc,
                                   double *xcur,
                                   double *rhs,
                                   double *crit)
{
  cs *tAicov, *Ajcov, *Bf, *B2, *B0, *Qicov;
  int ncova, nvar, ncur, nicur, njcur, error, nrows, ncols, signe;

  /* Initializations */

  error = 1;
  ncova = st_get_ncova();
  nvar = st_get_nvar();
  ncur = st_get_nvertex_max();
  tAicov = Bf = B2 = B0 = Qicov = nullptr;

  /* Loop on the first covariance */

  (*crit) = 0.;
  for (int icov = 0; icov < ncova; icov++)
  {
    SPDE_Matelem &Maticov = spde_get_current_matelem(icov);
    Qicov = Maticov.QC->Q;
    nicur = st_get_nvertex(icov);
    tAicov = cs_transpose(Maticov.Aproj, 1);
    if (tAicov == nullptr) goto label_end;

    st_keypair_cs("Qicov", Qicov, icov + 1, 0, 0, 0, 0);
    st_keypair_cs("Aproj", Maticov.Aproj, icov + 1, 0, 0, 0, 0);

    /* Loop on the first variable */

    for (int ivar = 0; ivar < nvar; ivar++)
    {
      for (int icur = 0; icur < nicur; icur++)
        rhscur[icur] = RHS(icov, ivar, icur);

      /* Loop on the second variable */

      for (int jvar = 0; jvar < nvar; jvar++)
      {

        // Suppress the Q terms

        if (flag_crit || ivar != jvar)
        {
          for (int icur = 0; icur < nicur; icur++)
            rhsloc[icur] = 0.;

          if (st_is_model_nugget())
          {
            // Sill(icov)_ij * Q(icov) * X_j
            cs_mulvec(Qicov, nicur, &XCUR(icov, jvar, 0), rhsloc);
            for (int icur = 0; icur < nicur; icur++)
              rhsloc[icur] *= st_get_isill(icov, ivar, jvar);
          }
          else
          {
            if (icov == 0)
            {
              // Q1_tt_ij * X_j
              B0 = st_extract_Q1_hetero(ivar, jvar, 2, 2, &nrows, &ncols);
              cs_mulvec(B0, nrows, &XCUR(icov, jvar, 0), rhsloc);
              B0 = cs_spfree(B0);
            }
            else
            {
              // Sill(icov)_ij * Q(icov) * X_j
              cs_mulvec(Qicov, nicur, &XCUR(icov, jvar, 0), rhsloc);
              for (int icur = 0; icur < nicur; icur++)
                rhsloc[icur] *= st_get_isill(icov, ivar, jvar);
            }
          }
          for (int icur = 0; icur < nicur; icur++)
            rhscur[icur] -= rhsloc[icur];
        }

        // Pre-calculation of factorized term

        if (st_is_model_nugget())
        {
          // A^t(icov) * E_ij
          B0 = st_extract_Q1_nugget(ivar, jvar, &nrows, &ncols);
          if (B0 == nullptr) goto label_end;
          Bf = cs_multiply(tAicov, B0);
          if (Bf == nullptr) goto label_end;
          B0 = cs_spfree(B0);
        }
        else
        {
          if (icov == 0)
          {
            // Q1_td_ij
            Bf = st_extract_Q1_hetero(ivar, jvar, 2, 1, &nrows, &ncols);
            if (Bf == nullptr) goto label_end;
          }
          else
          {
            // A^t(icov) * Q1_dd_ij
            B0 = st_extract_Q1_hetero(ivar, jvar, 1, 1, &nrows, &ncols);
            if (B0 == nullptr) goto label_end;
            Bf = cs_multiply(tAicov, B0);
            if (Bf == nullptr) goto label_end;
            B0 = cs_spfree(B0);
          }
        }

        /* Loop on the other structures */

        for (int jcov = 0; jcov < ncova; jcov++)
        {
          SPDE_Matelem &Matjcov = spde_get_current_matelem(jcov);
          njcur = st_get_nvertex(jcov);
          Ajcov = Matjcov.Aproj;

          // Suppress the AQAt terms

          if (flag_crit || (icov != jcov || ivar != jvar))
          {
            signe = 0;
            for (int icur = 0; icur < nicur; icur++)
              rhsloc[icur] = 0.;

            if (st_is_model_nugget())
            {
              // A^t(icov) * E_ij * A(jcov) * X_j
              cs_mulvec(Ajcov, njcur, &XCUR(jcov, jvar, 0), work);
              cs_mulvec(Bf, nicur, work, rhsloc);
              signe = 1;
            }
            else
            {
              if (icov > 0 && jcov == 0)
              {
                // -A^t(icov) * Q1_dt_ij * X_j 
                B0 = st_extract_Q1_hetero(ivar, jvar, 1, 2, &nrows, &ncols);
                cs_mulvec(B0, nrows, &XCUR(jcov, jvar, 0), work);
                cs_mulvec(tAicov, nicur, work, rhsloc);
                B0 = cs_spfree(B0);
                signe = -1;
              }
              else if (icov == 0 && jcov > 0)
              {
                // -Q1_td_ij * A(jcov) * X_j
                cs_mulvec(Ajcov, njcur, &XCUR(jcov, jvar, 0), work);
                cs_mulvec(Bf, nicur, work, rhsloc);
                signe = -1;
              }
              else if (icov > 0 && jcov > 0)
              {
                // A^t(icov) * Q1_dd_ij * A(jcov) * X_j
                cs_mulvec(Ajcov, njcur, &XCUR(jcov, jvar, 0), work);
                cs_mulvec(Bf, nicur, work, rhsloc);
                signe = +1;
              }
            }
            for (int icur = 0; icur < nicur; icur++)
              rhscur[icur] -= rhsloc[icur] * signe;
          }
        }
        Bf = cs_spfree(Bf);
      }

      if (flag_crit)
      {
        for (int icur = 0; icur < nicur; icur++)
          (*crit) += rhscur[icur] * rhscur[icur];
      }
      else
      {
        if (st_solve(Maticov.QCov[ivar], rhscur, work, &XCUR(icov, ivar, 0)))
          goto label_end;
      }
    }
  }

  /* Set the error return code */

  error = 0;

  label_end: tAicov = cs_spfree(tAicov);
  B0 = cs_spfree(B0);
  B2 = cs_spfree(B2);
  Bf = cs_spfree(Bf);
  return (error);
}

/****************************************************************************/
/*!
 **  Perform the final reconstitution for multivariate or multistructure
 **  Kriging
 **
 ** \return  Error return code
 **
 ** \param[in]  xcur        Solution vector (Dimension: ncur * ncova * nvar)
 **
 ** \param[out] z           Output Array
 **
 *****************************************************************************/
static int st_kriging_several_results(double *xcur, double *z)
{
  int *ranks, ncova, nvar, ncur, error, flag_data, ecr, lec;
  double valdat, valsum;
  SPDE_Mesh *s_mesh;

  /* Initializations */

  error = 1;
  ncova = st_get_ncova();
  nvar = st_get_nvar();
  ncur = st_get_nvertex_max();
  ranks = nullptr;
  valdat = TEST;
  s_mesh = spde_get_current_matelem(0).s_mesh;
  flag_data = 0;

  /* Sample designation */

  ranks = st_get_vertex_ranks(s_mesh, MEM_DBIN, MEM_DBOUT);
  if (ranks == nullptr) goto label_end;

  /* Constitute the resulting vector */

  if (st_is_model_nugget())
  {

    /* Case with Nugget Effect */

    ecr = 0;
    for (int ivar = 0; ivar < nvar; ivar++)
    {
      for (int icur = 0; icur < ncur; icur++, ecr++)
      {

        /* Check if the vertex corresponds to a target (or a Steiner point) */
        /* If it coincides with datum, we must still check that the variable */
        /* is defined (heterotopic case) */

        if (ranks[icur] <= 0)
        {
          flag_data = 0;
        }
        else
        {
          valdat = MEM_DBIN->getVariable(ranks[icur] - 1, ivar);
          flag_data = !FFFF(valdat);
        }

        /* Store the sum of kriged components */

        if (flag_data)
          valsum = valdat;
        else
        {
          valsum = 0.;
          for (int icov = 0; icov < ncova; icov++)
            valsum += XCUR(icov, ivar, icur);
        }
        z[ecr] = valsum;
      }
    }
  }
  else
  {

    /* Case with No Nugget Effect */

    ecr = 0;
    for (int ivar = 0; ivar < nvar; ivar++)
    {
      lec = 0;
      for (int icur = 0; icur < ncur; icur++, ecr++)
      {

        /* Check if the vertex corresponds to a target (or a Steiner point) */
        /* If it coincides with datum, we must still check that the variable */
        /* is defined (heterotopic case) */

        if (ranks[icur] <= 0)
        {
          flag_data = 0;
        }
        else
        {
          valdat = MEM_DBIN->getVariable(ranks[icur] - 1, ivar);
          flag_data = !FFFF(valdat);
        }

        if (flag_data)
        {
          valsum = valdat;
        }
        else
        {
          valsum = XCUR(0, ivar, lec++);
          for (int icov = 1; icov < ncova; icov++)
            valsum += XCUR(icov, ivar, icur);
        }
        z[ecr] = valsum;
      }
    }
  }

  if (FLAG_KEYPAIR) set_keypair("Res.final", 1, ecr / nvar, nvar, z);

  /* Set the error return code */

  error = 0;

  label_end: ranks = (int*) mem_free((char* ) ranks);
  return (error);
}

/****************************************************************************/
/*!
 **  Perform the Kriging operation (multivariate and/or multistructure)
 **
 ** \return  Error return code
 **
 ** \param[in]  data        Vector of data
 **
 ** \param[out] rhs         R.H.S.
 ** \param[out] work        Working array
 ** \param[out] z           Output Array
 **
 ** \remark This code is only provided for the Cholesky solver
 **
 *****************************************************************************/
static int st_kriging_several(double *data,
                              double *rhs,
                              double *work,
                              double *z)
{
  double *rhsloc, *rhscur, *xcur, ss, crit, tolmult;
  int ncur, ncova, error, maxiter, nvar;

  /* Initializations */

  error = 1;
  maxiter = (int) get_keypone("Multi_Structures_Maxiter", 200);
  tolmult = get_keypone("Multi_Structures_Tolerance", 1.e-8);
  rhsloc = rhscur = xcur = nullptr;
  if (S_DECIDE.flag_mgrid)
  {
    messerr("The multi-structure Kriging is not programmed in multigrid");
    return (1);
  }
  ncova = st_get_ncova();
  nvar = st_get_nvar();
  ncur = st_get_nvertex_max();
  ss = 0.;

  /* Core allocation */

  xcur = (double*) mem_alloc(sizeof(double) * ncur * ncova * nvar, 0);
  if (xcur == nullptr) goto label_end;
  rhsloc = (double*) mem_alloc(sizeof(double) * ncur, 0);
  if (rhsloc == nullptr) goto label_end;
  rhscur = (double*) mem_alloc(sizeof(double) * ncur, 0);
  if (rhscur == nullptr) goto label_end;

  /* Define the Initial solution */

  st_init_array(ncova, nvar, ncur, 0, xcur);

  /* Define the Right-hand sides */

  if (st_kriging_several_rhs(data, rhs, work, &ss)) goto label_end;
  st_keypair_array("RHS_init", -1, rhs);

  /* Loop on iterations */

  for (int iter = 0; iter < maxiter; iter++)
  {

    /* Forward process */

    if (st_kriging_several_loop(0, work, rhscur, rhsloc, xcur, rhs, &crit))
      goto label_end;
    if (iter < 2)
    {
      st_keypair_array("RHS", iter, rhs);
      st_keypair_array("XCUR", iter, xcur);
    }

    /* Calculate the convergence criterion */

    if (st_kriging_several_loop(1, work, rhscur, rhsloc, xcur, rhs, &crit))
      goto label_end;

    crit /= ss;
    if (VERBOSE) message("Iteration %3d - Criterion = %lg\n", iter + 1, crit);
    if (crit < tolmult) break;
  }

  /* Constitute the resulting vector */

  st_keypair_array("XCUR_results", -1, xcur);
  if (st_kriging_several_results(xcur, z)) goto label_end;

  /* Set the error return code */

  error = 0;

  label_end: xcur = (double*) mem_free((char* ) xcur);
  rhsloc = (double*) mem_free((char* ) rhsloc);
  rhscur = (double*) mem_free((char* ) rhscur);
  return (error);
}

/****************************************************************************/
/*!
 **  Perform the Kriging operation
 **
 ** \return  Error return code
 **
 ** \param[in]  s_mesh      SPDE_Mesh structure
 ** \param[in]  data        Vector of data
 **
 ** \param[out] zkrig       Output Array
 **
 *****************************************************************************/
static int st_kriging(SPDE_Mesh *s_mesh, double *data, double *zkrig)
{
  double *work, *zkdat, *rhs;
  int error, ncur, size, ncova, nvar;

  /* Initializations */

  if (VERBOSE || DEBUG)
  {
    if (S_DECIDE.flag_several)
      message("Kriging step (multiple path algorithm)\n");
    else
      message("Kriging step (single path algorithm)\n");
  }
  error = 1;
  work = zkdat = rhs = nullptr;
  nvar = st_get_nvar();
  ncova = st_get_ncova();
  ncur = st_get_nvertex_max();
  size = st_get_dimension();

  /* Core allocation */

  work = (double*) mem_alloc(sizeof(double) * size, 0);
  if (work == nullptr) goto label_end;
  zkdat = (double*) mem_alloc(sizeof(double) * ncur * nvar, 0);
  if (zkdat == nullptr) goto label_end;
  rhs = (double*) mem_alloc(sizeof(double) * ncur * ncova * nvar, 0);
  if (rhs == nullptr) goto label_end;
  for (int i = 0; i < ncur * ncova * nvar; i++)
    rhs[i] = 0.;
  for (int i = 0; i < ncur * nvar; i++)
    zkdat[i] = 0.;

  /* Define the Initial solution and Right-hand side */

  if (S_DECIDE.flag_several)
  {
    if (st_kriging_several(data, rhs, work, zkdat)) goto label_end;
    st_copy(s_mesh->vertype, VT_NONE, zkdat, zkrig);
  }
  else
  {
    if (st_kriging_one(data, rhs, work, zkdat)) goto label_end;
    st_merge(s_mesh->vertype, VT_FREE, zkdat, zkrig);
  }

  /* Set the error return code */

  error = 0;

  label_end: work = (double*) mem_free((char* ) work);
  rhs = (double*) mem_free((char* ) rhs);
  zkdat = (double*) mem_free((char* ) zkdat);
  return (error);
}

/****************************************************************************/
/*!
 **  Manage Cheb_Elem structure
 **
 ** \return  Error return code
 **
 ** \param[in]  mode       1 for allocation; -1 for deallocation
 ** \param[in]  verbose    Verbosity flag
 ** \param[in]  power      Parameter passed to Chebychev function
 ** \param[in]  nblin      Number of blin coefficients
 ** \param[in]  blin       Array of coefficients for Linear combinaison
 ** \param[in]  S          Shift operator
 ** \param[in]  cheb_old   Cheb_Elem to be freed (only for mode=-1)
 **
 ** \remarks Arguments 'power', 'nblin', 'blin' and 'B' are used if mode=1
 ** \remarks Argument 'cheb_old' is used if mode=-1
 **
 *****************************************************************************/
Cheb_Elem* spde_cheb_manage(int mode,
                            int verbose,
                            double power,
                            int nblin,
                            double *blin,
                            cs *S,
                            Cheb_Elem *cheb_old)
{
  Cheb_Elem *cheb_elem;
  double a, b, v1, v2, tol;
  int error, ncmax, ndisc;

  /* Initializations */

  error = 1;

  /* Dispatch */

  if (mode > 0)
  {

    // Allocation

    cheb_elem = new Cheb_Elem();
    if (cheb_elem == nullptr) goto label_end;
    cheb_elem->coeffs = nullptr;

    ncmax = (int) get_keypone("Number_Polynomials_Chebychev", 10001.);
    ndisc = (int) get_keypone("Number_Discretization_Chebychev", 100.);
    tol = get_keypone("Chebychev_Tolerance", 5.e-3);

    /* Calculate key values */

    a = 0.;
    b = cs_norm(S);
    v1 = 2. / (b - a);
    v2 = -(b + a) / (b - a);

    /* Store the values */

    cheb_elem->a = a;
    cheb_elem->b = b;
    cheb_elem->v1 = v1;
    cheb_elem->v2 = v2;
    cheb_elem->power = power;
    cheb_elem->ncmax = ncmax;
    cheb_elem->ndisc = ndisc;
    cheb_elem->tol = tol;
    cheb_elem->ncoeffs = 0;
    cheb_elem->coeffs = nullptr;

    /* Get the optimal count of Chebychev coefficients */

    if (st_chebychev_calculate_coeffs(cheb_elem, verbose, nblin, blin))
      goto label_end;
  }
  else
  {

    // Deallocation

    cheb_elem = cheb_old;
    if (cheb_elem != nullptr)
      cheb_elem->coeffs = (double*) mem_free((char* ) cheb_elem->coeffs);
    delete cheb_elem;
    cheb_elem = nullptr;
  }

  // Set the error return code

  error = 0;

  label_end: if (error)
    cheb_elem = spde_cheb_manage(-1, 0, 0, 0, NULL, NULL, cheb_elem);
  return (cheb_elem);
}

/****************************************************************************/
/*!
 **  Duplicate a Cheb_Elem structure
 **
 ** \return  Pointer to the newly created Cheb_Eleme structure
 **
 ** \param[in]  cheb_in    Input Cheb_Eleme structure
 **
 *****************************************************************************/
Cheb_Elem* _spde_cheb_duplicate(Cheb_Elem *cheb_in)
{
  Cheb_Elem *cheb_out;
  int error;

  // Initializations

  error = 1;
  cheb_out = nullptr;
  if (cheb_in == nullptr) return (cheb_out);

  // Allocation

  cheb_out = new Cheb_Elem();
  if (cheb_out == nullptr) goto label_end;

  cheb_out->ncoeffs = cheb_in->ncoeffs;
  cheb_out->ncmax = cheb_in->ncmax;
  cheb_out->ndisc = cheb_in->ndisc;
  cheb_out->power = cheb_in->power;
  cheb_out->a = cheb_in->a;
  cheb_out->b = cheb_in->b;
  cheb_out->v1 = cheb_in->v1;
  cheb_out->v2 = cheb_in->v2;
  cheb_out->tol = cheb_in->tol;
  cheb_out->coeffs = nullptr;

  cheb_out->coeffs = (double*) mem_alloc(sizeof(double) * cheb_in->ncoeffs, 0);
  if (cheb_out->coeffs == nullptr) goto label_end;
  for (int i = 0; i < cheb_in->ncoeffs; i++)
    cheb_out->coeffs[i] = cheb_in->coeffs[i];

  // Set the error return code

  error = 0;

  label_end: if (error)
    cheb_out = spde_cheb_manage(-1, 0, 0, 0, NULL, NULL, cheb_out);
  return (cheb_out);
}

/****************************************************************************/
/*!
 **  Initialize one SP_Mat structure
 **
 ** \param[in] mode    Type of the action
 **                    1 for allocation;
 **                    0 for partial deallocation (of current Matelem)
 **                   -1 for deallocation
 **
 ** \remarks This function is called when the current IGRF has been chosen
 **
 *****************************************************************************/
static void st_matelem_manage(int mode)

{
  int ncova = st_get_ncova();
  SPDE_SS_Environ *SS = st_get_current_ssenv();

  /* Dispatch */

  switch (mode)
  {
    case 1:                     // Allocation
      SS->Matelems.resize(ncova);

      for (int is = 0; is < ncova; is++)
      {
        SPDE_Matelem &Matelem = SS->Matelems[is];
        Matelem.S = nullptr;
        Matelem.Aproj = nullptr;
        Matelem.QC = nullptr;
        Matelem.QCov = nullptr;
        Matelem.Isill = nullptr;
        Matelem.Csill = nullptr;
        Matelem.qsimu = nullptr;
        Matelem.mgs = (cs_MGS*) NULL;
        if (S_DECIDE.flag_mgrid) Matelem.mgs = st_mgs_manage(1, NULL);
        Matelem.s_cheb = nullptr;
        Matelem.s_mesh = nullptr;
      }
      break;

    case -1:                    // Deallocation
      for (int icov = 0; icov < ncova; icov++)
      {
        SPDE_Matelem &Matelem = spde_get_current_matelem(icov);
        Matelem.S = cs_spfree(Matelem.S);
        Matelem.Aproj = cs_spfree(Matelem.Aproj);
        Matelem.QC = qchol_manage(-1, Matelem.QC);
        if (Matelem.QCov != NULL)
        {
          for (int ivar = 0; ivar < st_get_nvar(); ivar++)
            Matelem.QCov[ivar] = qchol_manage(-1, Matelem.QCov[ivar]);
        }
        Matelem.Isill = (double*) mem_free((char* ) Matelem.Isill);
        Matelem.Csill = (double*) mem_free((char* ) Matelem.Csill);
        Matelem.qsimu = st_qsimu_manage(-1, NULL, Matelem.qsimu);
        Matelem.mgs = st_mgs_manage(-1, Matelem.mgs);
        Matelem.s_cheb = spde_cheb_manage(-1, 0, 0, 0, NULL, NULL,
                                          Matelem.s_cheb);
        Matelem.s_mesh = spde_mesh_manage(-1, Matelem.s_mesh);
      }
      break;
  }
  return;
}

/****************************************************************************/
/*!
 **  Perform the Simulations
 **
 ** \return  Error return code
 **
 ** \param[in]  QC          Pointer to the QChol structure (finalized)
 **
 ** \param[out] zsnc        Output simulation array (Dimension: nvertex * nvar)
 **
 *****************************************************************************/
static int st_simulate(QChol *QC, double *zsnc)
{
  int error, ncur, nvar, ncova, iad, nvs2;
  double *work, *zloc;

  /* Initializations */

  error = 1;
  ncur = st_get_nvertex_max();
  nvar = st_get_nvar();
  ncova = st_get_ncova();
  nvs2 = nvar * (nvar + 1) / 2;
  zloc = work = nullptr;

  /* Core allocation */

  work = (double*) mem_alloc(sizeof(double) * ncur, 0);
  if (work == nullptr) goto label_end;
  zloc = (double*) mem_alloc(sizeof(double) * ncur, 0);
  if (zloc == nullptr) goto label_end;
  for (int i = 0; i < nvar * ncur; i++)
    zsnc[i] = 0.;

  /* Loop on the covariances */

  for (int icov = 0; icov < ncova; icov++)
  {
    st_set_current_icov(icov);
    SPDE_Matelem &Matelem = spde_get_current_matelem(icov);

    /* Loop on the variables */

    for (int jvar = 0; jvar < nvar; jvar++)
    {
      /* Simulate the continuous part */

      if (S_DECIDE.simu_chol)
      {
        if (st_simulate_cholesky(QC, work, zloc)) return (1);
      }
      else
      {
        if (st_simulate_chebychev(zloc)) return (1);
      }

      /* Simulate the nugget effect */

      st_simulate_nugget(ncur, zloc);

      /* Update the array for non-conditional simulation */

      for (int ivar = 0; ivar < nvar; ivar++)
      {
        iad = st_get_rank(ivar, jvar);
        for (int icur = 0; icur < ncur; icur++)
          zsnc[icur + ivar * ncur] += zloc[icur]
              * Matelem.Csill[iad + nvs2 * icov];
      }
    }
  }

  /* Set the error return code */

  error = 0;

  label_end: zloc = (double*) mem_free((char* ) zloc);
  work = (double*) mem_free((char* ) work);
  return (error);
}

/****************************************************************************/
/*!
 **  Perform the main Simulations steps after Q has been constructed
 **
 ** \return  Error return code
 **
 ** \param[in]  dbin        Input Db grid structure (optional)
 ** \param[in]  dbout       Output Db grid structure
 ** \param[in]  s_option    SPDE_Option structure
 ** \param[in]  nbsimu      Number of simulations
 ** \param[in]  ngibbs_burn Number of iterations (Burning step)
 ** \param[in]  ngibbs_iter Number of iterations
 ** \param[in]  ngibbs_int  Number of iterations internal to Gibbs (SPDE)
 **
 ** \remarks  If the number of simulations 'nbsimu' is set to 0,
 ** \remarks  the simulation algorithm is turned into a kriging one
 **
 *****************************************************************************/
int spde_process(Db *dbin,
                 Db *dbout,
                 SPDE_Option &s_option,
                 int nbsimu,
                 int ngibbs_burn,
                 int ngibbs_iter,
                 int ngibbs_int)
{
  int ncur, ndata, nbsimuw, isimuw, error, iatt_simu, ivar0;
  int flag_mult_data, ngrf, nv_krige;
  int ngibbs_total, ngtime, nvar;
  double *data, *zcur, *zkrig, *zout, *vcur, *zsnc, *zdat;
  SPDE_Mesh *s_mesh;

  /* Initializations */

  flag_mult_data = (int) get_keypone("Flag_Mult_Data", 0);
  error = 1;
  data = zcur = zkrig = zout = vcur = zsnc = zdat = nullptr;
  ndata = 0;
  ngrf = st_get_number_grf();
  nvar = st_get_nvar();
  ncur = st_get_nvertex_max();

  /* Title (optional) */

  if (VERBOSE) st_title(0, 0, 1, "Processing");

  /* Core allocation and global initializations */

  zcur = (double*) mem_alloc(sizeof(double) * ncur * nvar, 0);
  if (zcur == nullptr) goto label_end;

  if (S_DECIDE.flag_case == CASE_SIMULATE)
  {
    zsnc = (double*) mem_alloc(sizeof(double) * ncur * nvar, 0);
    if (zsnc == nullptr) goto label_end;
    zout = (double*) mem_alloc(sizeof(double) * ncur * nvar, 0);
    if (zout == nullptr) goto label_end;
  }
  if (S_DECIDE.flag_dbin)
  {
    ndata = dbin->getActiveSampleNumber();
    zkrig = (double*) mem_alloc(sizeof(double) * ncur * nvar, 0);
    if (zkrig == nullptr) goto label_end;
    zdat = (double*) mem_alloc(sizeof(double) * ncur * nvar, 0);
    if (zdat == nullptr) goto label_end;
    data = (double*) mem_alloc(sizeof(double) * ndata * nvar, 0);
    if (data == nullptr) goto label_end;
    if (S_DECIDE.flag_std)
    {
      vcur = (double*) mem_alloc(sizeof(double) * ncur * nvar, 0);
      if (vcur == nullptr) goto label_end;
    }
  }

  /***********************/
  /* Preliminary Kriging */
  /***********************/

  if (S_DECIDE.flag_onechol && !flag_mult_data && !S_DECIDE.flag_gibbs
      && ngrf == 1)
  {
    st_set_current_igrf(0);
    s_mesh = spde_get_current_matelem(-1).s_mesh;
    st_init_array(1, nvar, ncur, 1, zkrig);
    st_load_data(s_mesh, dbin, dbout, s_option, -1, data, zkrig);
    if (st_kriging(s_mesh, data, zkrig)) goto label_end;
  }

  if (S_DECIDE.flag_case == CASE_KRIGING)
  {
    // Calculation of the Kriging Variance (optional)
    if (S_DECIDE.flag_std) spde_build_stdev(vcur);

    // Saving operation
    nv_krige = 0;
    s_mesh = spde_get_current_matelem(-1).s_mesh;
    if (S_DECIDE.flag_est)
      st_save_result(s_mesh->vertype, zkrig, dbout, ELoc::Z, nv_krige++);
    if (S_DECIDE.flag_std)
      st_save_result(s_mesh->vertype, vcur, dbout, ELoc::Z, nv_krige++);
  }
  else if (S_DECIDE.flag_case == CASE_SIMULATE)
  {
    /***************/
    /* Simulations */
    /***************/

    nbsimuw = (S_DECIDE.flag_modif) ? 1 :
                                      nbsimu;

    /***************************/
    /* Loop on the simulations */
    /***************************/

    for (int isimu = 0; isimu < nbsimu; isimu++)
    {
      if (VERBOSE || DEBUG) message("Simulation #%d/%d\n", isimu + 1, nbsimu);
      isimuw = (S_DECIDE.flag_modif) ? 0 :
                                       isimu;
      st_init_array(1, nvar, ncur, 1, zcur);
      if (S_DECIDE.flag_std) st_init_array(1, nvar, ncur, 0, vcur);

      /****************************************************/
      /* Loop on the underlying Gaussian Random Functions */
      /****************************************************/

      for (int igrf = 0; igrf < ngrf; igrf++)
      {
        if ((VERBOSE || DEBUG) && st_get_number_grf() > 1)
          message("GRF iteration #%d/%d\n", igrf + 1, ngrf);
        st_set_current_igrf(igrf);
        s_mesh = spde_get_current_matelem(-1).s_mesh;

        /* Conditional Simulation */
        /* Load the data set corresponding to the new variable */

        if (S_DECIDE.flag_dbin)
        {
          st_init_array(1, nvar, ncur, 1, zdat);
          ivar0 = (flag_mult_data) ? isimu :
                                     -1;
          st_load_data(s_mesh, dbin, dbout, s_option, ivar0, data, zdat);
        }

        /***********************************************************/
        /* Loop on the gibbs iterations (performed at end of loop) */
        /***********************************************************/

        ngtime = 1;
        ngibbs_total = ngibbs_burn + ngibbs_iter;
        if (s_mesh->vertype->ngibbs <= 0) ngibbs_total = 0;
        if (S_DECIDE.flag_gibbs) ngtime = MAX(1, ngibbs_total);

        for (int igtime = 0; igtime < ngtime; igtime++)
        {
          if ((VERBOSE || DEBUG) && S_DECIDE.flag_gibbs)
            message("Gibbs iteration #%d/%d\n", igtime + 1, ngtime);

          /* Perform the non-conditional simulation */
          st_init_array(1, nvar, ncur, 1, zsnc);
          if (!S_DECIDE.flag_onechol)
          {
            // Simulate the continuous part
            if (st_simulate(spde_get_current_matelem(-1).QC, zsnc))
              goto label_end;

            if (S_DECIDE.flag_dbin)
            {
              // Calculate the simulation error
              st_simu_subtract_data(s_mesh->vertype, VT_GIBBS | VT_HARD, zsnc,
                                    data, zdat);

              // Perform Conditional Kriging
              if (st_kriging(s_mesh, zdat, zcur)) goto label_end;

              // Add the non-conditional simulation
              st_simu_add_vertices(ncur, zsnc, zcur);
            }
            else
            {
              st_load_array(ncur, zsnc, zcur);
            }
          }
          else
          {
            // Perform simulation of the residuals
            if (st_simulate(spde_get_current_matelem(-1).qsimu->QCtt, zout))
              goto label_end;
            st_merge(s_mesh->vertype, VT_FREE, zout, zsnc);

            if (S_DECIDE.flag_gibbs)
            {
              // Perform Kriging of the Gibbsed Gaussian
              if (st_kriging(s_mesh, data, zcur)) goto label_end;

              // Add the non-conditional simulation
              st_simu_add_vertices(ncur, zsnc, zcur);
            }
            else
            {
              st_load_array(ncur, zkrig, zcur);
              st_simu_add_vertices(ncur, zsnc, zcur);
            }
          }

          /* Gibbs iteration */

          if (S_DECIDE.flag_gibbs && s_mesh->vertype->ngibbs > 0)
          {
            st_gibbs(igrf, s_mesh->vertype, ngibbs_int, igtime, ngibbs_burn,
                     dbin, dbout, zcur);
          }
        }

        /* Saving operation */

        iatt_simu = dbout->getSimvarRank(isimuw, 0, igrf, nbsimuw, 1);
        st_save_result(s_mesh->vertype, zcur, dbout, ELoc::SIMU, iatt_simu);
      }

      /* Perform the transformation */

      if (SIMU_FUNC_TRANSF != NULL)
        SIMU_FUNC_TRANSF(dbout, VERBOSE, isimuw, nbsimuw);

      /* Transforming the result (flag_modif) */

      if (S_DECIDE.flag_modif && SIMU_FUNC_UPDATE != NULL)
        SIMU_FUNC_UPDATE(dbout, VERBOSE, isimuw, nbsimuw);
    }

    /* Scale the simulations */

    if (S_DECIDE.flag_modif && SIMU_FUNC_SCALE != NULL)
      SIMU_FUNC_SCALE(dbout, VERBOSE, nbsimu);
  }

  /* Set the error return code */

  error = 0;

  label_end: data = (double*) mem_free((char* ) data);
  zdat = (double*) mem_free((char* ) zdat);
  zkrig = (double*) mem_free((char* ) zkrig);
  zout = (double*) mem_free((char* ) zout);
  zsnc = (double*) mem_free((char* ) zsnc);
  zcur = (double*) mem_free((char* ) zcur);
  vcur = (double*) mem_free((char* ) vcur);
  return (error);
}

/****************************************************************************/
/*!
 **  Load the segments and vertex coordinates
 **
 ** \return  Error return code
 **
 ** \param[in]  verbose    Verbose flag
 ** \param[in]  dbin       Db structure for the conditioning data
 ** \param[in]  dbout      Db structure of the grid
 ** \param[in]  gext       Array of domain dilation
 ** \param[in]  nb_dupl    Number of samples to be masked (optional)
 ** \param[in]  is_dupl    Array of indices of the samples to be masked
 ** \param[in]  s_option   SPDE_Option structure
 **
 ** \param[in,out]  s_mesh SPDE_Mesh structure
 **
 *****************************************************************************/
static int st_load_meshes_1D(int verbose,
                             Db *dbin,
                             Db *dbout,
                             const VectorDouble &gext,
                             int nb_dupl,
                             int *is_dupl,
                             SPDE_Option &s_option,
                             SPDE_Mesh *s_mesh)
{
  int error;
  segmentio in, out;

  /* Initializations */

  error = 1;

  /* Initialize the Meshing output structure */

  meshes_1D_init(&in);
  meshes_1D_init(&out);

  /* Set the control points for the triangulation */

  if (dbout != nullptr && S_DECIDE.flag_mesh_dbout)
  {
    if (meshes_1D_from_db(dbout, 0, NULL, &in)) goto label_end;
  }
  if (dbin != nullptr && S_DECIDE.flag_mesh_dbin)
  {
    if (meshes_1D_from_db(dbin, nb_dupl, is_dupl, &in)) goto label_end;
  }
  if (!(dbin != nullptr && S_DECIDE.flag_mesh_dbin) && !(dbout != nullptr
      && S_DECIDE.flag_mesh_dbout))
  {
    meshes_1D_default(dbin, dbout, &in);
  }

  /* Extend the domain if gext is specified */

  if (!gext.empty())
  {
    meshes_1D_extended_domain(dbout, gext.data(), &in);
  }

  /* Perform the triangulation */

  meshes_1D_create(VERBOSE, &in, &out);

  /* Coordinates of the triangle vertices */

  meshes_1D_load_vertices(&out, "Points", &s_mesh->nvertex, &s_mesh->ndim,
                          (void**) &s_mesh->points);

  meshes_1D_load_vertices(&out, "Segments", &s_mesh->nmesh, &s_mesh->ncorner,
                          (void**) &s_mesh->meshes);

  /* Set the error return code */

  error = 0;

  label_end: meshes_1D_free(&in, 1);
  meshes_1D_free(&out, 0);
  return (error);
}

/****************************************************************************/
/*!
 **  Load the triangles and vertex coordinates
 **
 ** \return  Error return code
 **
 ** \param[in]  verbose    Verbose flag
 ** \param[in]  dbin       Db structure for the conditioning data
 ** \param[in]  dbout      Db structure of the grid
 ** \param[in]  gext       Array of domain dilation
 ** \param[in]  nb_dupl    Number of samples to be masked (optional)
 ** \param[in]  is_dupl    Array of indices of the samples to be masked
 ** \param[in]  s_option   SPDE_Option structure
 **
 ** \param[in,out]  s_mesh SPDE_Mesh structure
 **
 *****************************************************************************/
static int st_load_meshes_2D(int verbose,
                             Db *dbin,
                             Db *dbout,
                             const VectorDouble &gext,
                             int nb_dupl,
                             int *is_dupl,
                             SPDE_Option &s_option,
                             SPDE_Mesh *s_mesh)
{
  int error;
  triangulateio in, out, vorout;

  /* Initializations */

  error = 1;

  /* Initialize the Meshing output structure */

  meshes_2D_init(&in);
  meshes_2D_init(&out);
  meshes_2D_init(&vorout);

  /* Set the control points for the triangulation */

  if (dbout != nullptr && S_DECIDE.flag_mesh_dbout)
  {
    if (meshes_2D_from_db(dbout, 1, 0, NULL, &in)) goto label_end;
  }
  if (dbin != nullptr && S_DECIDE.flag_mesh_dbin)
  {
    if (meshes_2D_from_db(dbin, 1, nb_dupl, is_dupl, &in)) goto label_end;
  }
  if (!(dbin != nullptr && S_DECIDE.flag_mesh_dbin) && !(dbout != nullptr
      && S_DECIDE.flag_mesh_dbout))
  {
    meshes_2D_default(dbin, dbout, &in);
  }

  /* Extend the domain if gext is specified */

  if (!gext.empty())
  {
    meshes_2D_extended_domain(dbout, gext.data(), &in);
  }

  /* Perform the triangulation */

  meshes_2D_create(VERBOSE, st_get_current_triswitch(s_option), &in, &out,
                   &vorout);

  /* Coordinates of the triangle vertices */

  meshes_2D_load_vertices(&out, "Points", &s_mesh->nvertex, &s_mesh->ndim,
                          (void**) &s_mesh->points);

  meshes_2D_load_vertices(&out, "Triangles", &s_mesh->nmesh, &s_mesh->ncorner,
                          (void**) &s_mesh->meshes);

  /* Set the error return code */

  error = 0;

  label_end: meshes_2D_free(&in, 1);
  meshes_2D_free(&out, 0);
  meshes_2D_free(&vorout, 0);
  return (error);
}

/****************************************************************************/
/*!
 **  Load the spherical triangles and vertex coordinates
 **
 ** \return  Error return code
 **
 ** \param[in]  verbose    Verbose flag
 ** \param[in]  dbin       Db structure for the conditioning data
 ** \param[in]  dbout      Db structure of the grid
 ** \param[in]  nb_dupl    Number of samples to be masked (optional)
 ** \param[in]  is_dupl    Array of indices of the samples to be masked
 ** \param[in]  s_option   SPDE_Option structure
 **
 ** \param[in,out]  s_mesh SPDE_Mesh structure
 **
 *****************************************************************************/
static int st_load_meshes_2D_sph(int verbose,
                                 Db *dbin,
                                 Db *dbout,
                                 int nb_dupl,
                                 int *is_dupl,
                                 SPDE_Option &s_option,
                                 SPDE_Mesh *s_mesh)
{
  int error;
  SphTriangle in;

  /* Initializations */

  error = 1;

  /* Initialize the Meshing output structure */

  meshes_2D_sph_init(&in);

  /* Set the control points for the triangulation */

  if (dbout != nullptr && S_DECIDE.flag_mesh_dbout)
  {
    if (meshes_2D_sph_from_db(dbout, 0, NULL, &in)) goto label_end;
  }
  if (dbin != nullptr && S_DECIDE.flag_mesh_dbin)
  {
    if (meshes_2D_sph_from_db(dbin, nb_dupl, is_dupl, &in)) goto label_end;
  }

  /* Add auxiliary random points */

  if (meshes_2D_sph_from_auxiliary(st_get_current_triswitch(s_option), &in))
    goto label_end;

  /* Perform the triangulation */

  if (meshes_2D_sph_create(VERBOSE, &in)) goto label_end;

  /* Coordinates of the triangle vertices */

  meshes_2D_sph_load_vertices(&in, "Points", &s_mesh->nvertex, &s_mesh->ndim,
                              (void**) &s_mesh->points);

  meshes_2D_sph_load_vertices(&in, "Triangles", &s_mesh->nmesh,
                              &s_mesh->ncorner, (void**) &s_mesh->meshes);

  /* Set the error return code */

  error = 0;

  label_end: meshes_2D_sph_free(&in, 0);
  return (error);
}

/****************************************************************************/
/*!
 **  Load the tetrahedra and vertex coordinates
 **
 ** \return  Error return code
 **
 ** \param[in]  verbose    Verbose flag
 ** \param[in]  dbin       Db structure for the conditioning data
 ** \param[in]  dbout      Db structure of the grid
 ** \param[in]  gext       Array of domain dilation
 ** \param[in]  nb_dupl    Number of samples to be masked (optional)
 ** \param[in]  is_dupl    Array of indices of the samples to be masked
 ** \param[in]  s_option   SPDE_Option structure
 **
 ** \param[in,out]  s_mesh SPDE_Mesh structure
 **
 *****************************************************************************/
static int st_load_meshes_3D(int verbose,
                             Db *dbin,
                             Db *dbout,
                             const VectorDouble &gext,
                             int nb_dupl,
                             int *is_dupl,
                             SPDE_Option &s_option,
                             SPDE_Mesh *s_mesh)
{
  int error;
  tetgenio in, out;

  /* Initializations */

  error = 1;

  /* Set the control points for the tetrahedralization */

  if (dbout != nullptr && S_DECIDE.flag_mesh_dbout)
  {
    if (meshes_3D_from_db(dbout, 0, NULL, &in)) goto label_end;
  }
  if (dbin != nullptr && S_DECIDE.flag_mesh_dbin)
  {
    if (meshes_3D_from_db(dbin, nb_dupl, is_dupl, &in)) goto label_end;
  }
  if (!(dbin != nullptr && S_DECIDE.flag_mesh_dbin) && !(dbout != nullptr
      && S_DECIDE.flag_mesh_dbout))
  {
    meshes_3D_default(dbin, dbout, &in);
  }

  /* Extend the domain if gext is specified */

  if (!gext.empty())
  {
    meshes_3D_extended_domain(dbout, gext.data(), &in);
  }

  /* Perform the tetrahedralization */

  meshes_3D_create(VERBOSE, st_get_current_triswitch(s_option), &in, &out);

  /* Coordinates of the vertices */

  meshes_3D_load_vertices(&out, "Points", &s_mesh->nvertex, &s_mesh->ndim,
                          (void**) &s_mesh->points);

  meshes_3D_load_vertices(&out, "Tetrahedra", &s_mesh->nmesh, &s_mesh->ncorner,
                          (void**) &s_mesh->meshes);

  /* Set the error return code */

  error = 0;

  label_end: meshes_3D_free(&in);
  meshes_3D_free(&out);
  return (error);
}

/****************************************************************************/
/*!
 **  Save the meshing in keypair
 **
 ** \param[in]  s_mesh    SPDE_Mesh structure
 ** \param[in]  icov0     Rank of the covariance corresponding to the Mesh
 **
 *****************************************************************************/
static void st_save_meshing_keypair(SPDE_Mesh *s_mesh, int icov0)
{
  int flag_save;

  flag_save = (int) get_keypone("Meshing_External_Save", 0);
  if (!flag_save) return;

  (void) gslSPrintf(NAME, "Meshing_External_Points.%d", icov0 + 1);
  set_keypair(NAME, 1, s_mesh->nvertex, s_mesh->ndim, s_mesh->points);
  (void) gslSPrintf(NAME, "Meshing_External_Meshes.%d", icov0 + 1);
  set_keypair_int(NAME, 1, s_mesh->nmesh, s_mesh->ncorner, s_mesh->meshes);
}

/****************************************************************************/
/*!
 **  Load the meshes
 **
 ** \return  Error return code
 **
 ** \param[in]  dbin       Db structure for the conditioning data
 ** \param[in]  dbout      Db structure of the grid
 ** \param[in]  gext       Array of domain dilation
 ** \param[in]  s_option   SPDE_Option structure
 **
 ** \param[in,out]  s_mesh SPDE_Mesh structure
 **
 ** \remarks The option 'flag_force' forces to use the regular meshing rather
 ** \remarks than the Turbo one
 **
 *****************************************************************************/
static int st_load_all_meshes(Db *dbin,
                              Db *dbout,
                              const VectorDouble &gext,
                              SPDE_Option &s_option,
                              SPDE_Mesh *s_mesh)
{
  int *is_dupl, ndim_loc, flag_sphere, nb_dupl, flag_force;
  Db *dbloc;

  flag_force = (int) get_keypone("Force_Regular_Meshing", 0);
  if (VERBOSE)
  {
    message("Generating the meshes\n");
    if (!S_DECIDE.flag_mesh_dbin)
      message("- Input data do not participate to the Meshing\n");
    if (!S_DECIDE.flag_mesh_dbout)
      message("- Output targets do not participate to the Meshing\n");
  }

  ndim_loc = 0;
  if (dbin != nullptr) ndim_loc = MAX(ndim_loc, dbin->getNDim());
  if (dbout != nullptr) ndim_loc = MAX(ndim_loc, dbout->getNDim());
  variety_query(&flag_sphere);

  // Manage the mask 

  if (s_mesh->vercoloc != nullptr)
  {
    nb_dupl = s_mesh->vercoloc->ndupl;
    is_dupl = s_mesh->vercoloc->dupl_dabs;
  }
  else
  {
    nb_dupl = 0;
    is_dupl = nullptr;
  }

  // Processing

  if (flag_sphere)
  {

    /* Particular case of data on the sphere */

    if (st_load_meshes_2D_sph(VERBOSE, dbin, dbout, nb_dupl, is_dupl, s_option,
                              s_mesh)) return (1);
  }
  else
  {

    /* Standard case */

    /* Check that a single file must be meshed and that it corresponds */
    /* to a grid */

    dbloc = NULL;
    if (!flag_force)
    {
      if (((!S_DECIDE.flag_dbin || !S_DECIDE.flag_mesh_dbin) && is_grid(dbout)))
        dbloc = dbout;
      if (((!S_DECIDE.flag_dbout || !S_DECIDE.flag_mesh_dbout) && is_grid(dbin)))
        dbloc = dbin;
    }

    if (dbloc != NULL)
    {
      if (VERBOSE) message("Using Turbo Meshing\n");

      /* Regular meshing */

      if (ndim_loc == 1)
      {
        if (meshes_turbo_1D_grid_build(VERBOSE, dbloc, s_mesh)) return (1);
      }
      else if (ndim_loc == 2)
      {
        if (meshes_turbo_2D_grid_build(VERBOSE, dbloc, s_mesh)) return (1);
      }
      else if (ndim_loc == 3)
      {
        if (meshes_turbo_3D_grid_build(VERBOSE, dbloc, s_mesh)) return (1);
      }

    }
    else
    {
      if (VERBOSE) message("Using Regular Meshing\n");

      if (ndim_loc == 1)
      {
        if (st_load_meshes_1D(VERBOSE, dbin, dbout, gext, nb_dupl, is_dupl,
                              s_option, s_mesh)) return (1);
      }
      else if (ndim_loc == 2)
      {
        if (st_load_meshes_2D(VERBOSE, dbin, dbout, gext, nb_dupl, is_dupl,
                              s_option, s_mesh)) return (1);
      }
      else if (ndim_loc == 3)
      {
        if (st_load_meshes_3D(VERBOSE, dbin, dbout, gext, nb_dupl, is_dupl,
                              s_option, s_mesh)) return (1);
      }
    }
  }
  if (s_mesh->nvertex <= 0 || s_mesh->nmesh <= 0) return (1);
  return (0);
}

/****************************************************************************/
/*!
 **  Check if External Meshing has been defined for the current covariance
 **
 ** \return 1 if an external mesh has been defined; 0 otherwise
 **
 ** \param[in] icov0    Rank of the current Covariance
 **
 *****************************************************************************/
static int st_is_external_mesh_defined(int icov0)
{
  return (S_EXTERNAL_MESH[icov0] != nullptr);
}

/****************************************************************************/
/*!
 **  Check if External Q has been defined
 **
 ** \return 1 if External A-Q has been defined (for 'icov'); 0 otherwise
 **
 ** \param[in]  icov0     Rank of the (non-nugget) covariance
 **
 *****************************************************************************/
static int st_is_external_AQ_defined(int icov0)
{
  return (S_EXTERNAL_MESH[icov0] != nullptr && S_EXTERNAL_Q[icov0] != nullptr
          && S_EXTERNAL_A[icov0] != nullptr);
}

/****************************************************************************/
/*!
 **  Copy the contents of the internal S_EXTERNAL_MESH into an output SPDE_Mesh
 **
 **  Error retur code
 **
 ** \param[in]  s_mesh   Output SPDE_Mesh
 ** \param[in]  icov0    Rank of the current Covariance
 **
 *****************************************************************************/
int spde_external_mesh_copy(SPDE_Mesh *s_mesh, int icov0)
{
  int size;

  if (S_EXTERNAL_MESH[icov0] == nullptr)
  {
    messerr("The Internal SPDE_Mesh must be allocated before using it");
    return (1);
  }
  if (s_mesh == nullptr)
  {
    messerr("The output SPDE_Mesh must already exist");
    return (1);
  }

  s_mesh->ndim = S_EXTERNAL_MESH[icov0]->ndim;
  s_mesh->ncorner = S_EXTERNAL_MESH[icov0]->ncorner;
  s_mesh->nvertex = S_EXTERNAL_MESH[icov0]->nvertex;
  s_mesh->nmesh = S_EXTERNAL_MESH[icov0]->nmesh;

  /* Copy the array 'points' */

  s_mesh->points = (double*) mem_free((char* ) s_mesh->points);
  size = S_EXTERNAL_MESH[icov0]->nvertex * S_EXTERNAL_MESH[icov0]->ndim;
  s_mesh->points = (double*) mem_alloc(sizeof(double) * size, 0);
  if (s_mesh->points == nullptr) return (1);
  for (int i = 0; i < size; i++)
    s_mesh->points[i] = S_EXTERNAL_MESH[icov0]->points[i];

  /* Copy the array 'meshes' */

  s_mesh->meshes = (int*) mem_free((char* ) s_mesh->meshes);
  size = S_EXTERNAL_MESH[icov0]->nmesh * S_EXTERNAL_MESH[icov0]->ncorner;
  s_mesh->meshes = (int*) mem_alloc(sizeof(double) * size, 0);
  if (s_mesh->meshes == nullptr) return (1);
  for (int i = 0; i < size; i++)
    s_mesh->meshes[i] = S_EXTERNAL_MESH[icov0]->meshes[i];

  if (S_EXTERNAL_MESH[icov0]->vertype != nullptr)
  {
    s_mesh->vertype = vertype_manage(1, NULL, NULL, s_mesh->nvertex);
    if (s_mesh->vertype == nullptr) return (1);

    for (int i = 0; i < S_EXTERNAL_MESH[icov0]->nvertex; i++)
      s_mesh->vertype->vt[i] = S_EXTERNAL_MESH[icov0]->vertype->vt[i];
  }

  return (0);
}

/****************************************************************************/
/*!
 **  Copy the contents of the internal S_EXTERNAL_AQ into an output SPDE_Mesh
 **
 **  Error retur code
 **
 ** \param[in]  matelem  Output SPDE_Matelem structure
 ** \param[in]  icov0    Rank of the current Covariance
 **
 *****************************************************************************/
int spde_external_AQ_copy(SPDE_Matelem &matelem, int icov0)
{
  SPDE_Mesh *s_mesh;

  s_mesh = matelem.s_mesh;
  if (S_EXTERNAL_MESH[icov0] == nullptr)
  {
    messerr("The Internal SPDE_Mesh must be allocated before using it");
    return (1);
  }
  if (S_EXTERNAL_A[icov0] == nullptr)
  {
    messerr("The External A must be allocated before using it");
    return (1);
  }
  if (S_EXTERNAL_Q[icov0] == nullptr)
  {
    messerr("The External Q must be allocated before using it");
    return (1);
  }
  if (s_mesh == nullptr)
  {
    messerr("The output SPDE_Mesh must already exist");
    return (1);
  }

  s_mesh->ndim = S_EXTERNAL_MESH[icov0]->ndim;
  s_mesh->nvertex = S_EXTERNAL_MESH[icov0]->nvertex;
  s_mesh->nmesh = S_EXTERNAL_MESH[icov0]->nmesh;

  /* Copy the sparse matrix 'QC' */

  matelem.QC = qchol_manage(1, NULL);
  matelem.QC->Q = cs_duplicate(S_EXTERNAL_Q[icov0]);
  matelem.Aproj = cs_duplicate(S_EXTERNAL_A[icov0]);

  if (S_EXTERNAL_MESH[icov0]->vertype != nullptr)
  {
    s_mesh->vertype = vertype_manage(1, NULL, NULL, s_mesh->nvertex);
    if (s_mesh->vertype == nullptr) return (1);

    for (int i = 0; i < S_EXTERNAL_MESH[icov0]->nvertex; i++)
      s_mesh->vertype->vt[i] = S_EXTERNAL_MESH[icov0]->vertype->vt[i];
  }

  return (0);
}

/****************************************************************************/
/*!
 **  Load the SPDE_Mesh structure
 **
 ** \return  Error return code
 **
 ** \param[in,out] s_mesh  Pointer to SPDE_Mesh to be loaded
 ** \param[in]  verbose    Verbose option
 ** \param[in]  dbin       Db structure for the conditioning data
 ** \param[in]  dbout      Db structure of the grid
 ** \param[in]  gext       Array of domain dilation
 ** \param[in]  s_option   SPDE_Option structure
 **
 *****************************************************************************/
int spde_mesh_load(SPDE_Mesh *s_mesh,
                   int verbose,
                   Db *dbin,
                   Db *dbout,
                   const VectorDouble &gext,
                   SPDE_Option &s_option)
{
  int icov0;

  if (s_mesh == nullptr) return (0);
  VERBOSE = verbose;
  icov0 = st_get_current_icov();

  /* Load the meshes */

  if (st_is_external_mesh_defined(icov0))
  {
    spde_external_mesh_copy(s_mesh, icov0);
  }
  else
  {

    /* Look for duplicates */

    s_mesh->vercoloc = vercoloc_manage(verbose, 1, dbin, dbout,
                                       S_DECIDE.flag_mesh_dbin, NULL);
    if (s_mesh->vercoloc == nullptr) return (1);

    /* Load the meshing */

    if (st_load_all_meshes(dbin, dbout, gext, s_option, s_mesh)) return (1);

    /* Manage the Vertype structure */

    s_mesh->vertype = vertype_manage(1, NULL, s_mesh->vercoloc,
                                     s_mesh->nvertex);
    if (s_mesh->vertype == nullptr) return (1);

    /* Load the vertex identification array */

    st_vertype_load(s_mesh->vertype, s_mesh->vercoloc, dbin, dbout, s_option);
  }
  st_save_meshing_keypair(s_mesh, icov0);

  return (0);
}

/****************************************************************************/
/*!
 **  Find the index of a duplicate
 **
 ** \return Rank of the matching duplicate or -1
 **
 ** \param[in]  rank      Rank of the sample
 ** \param[in]  ndupl     Number of duplicates
 ** \param[in]  dupl      Array of duplicates
 **
 *****************************************************************************/
static int st_find_in_duplicate(int rank, int ndupl, int *dupl)
{
  if (ndupl <= 0 || dupl == nullptr) return (-1);
  if (!S_DECIDE.flag_mesh_dbin || !S_DECIDE.flag_mesh_dbout) return (-1);
  for (int i = 0; i < ndupl; i++)
  {
    if (dupl[i] == rank) return (i);
  }
  return (-1);
}

/****************************************************************************/
/*!
 **  Define the Vertype structure from External information
 **
 ** \param[in]  s_mesh    Pointer on the already existing SPDE_Mesh structure
 ** \param[in]  nvertex   Number of vertices
 ** \param[in]  nbin      Number of vertices dedicated to the Input File
 ** \param[in]  nbout     Number of vertices dedicated to the Output File
 ** \param[in]  ndupl     Number of duplicates
 ** \param[in]  order     1: dbin then dbout; 2: dbout then dbin
 ** \param[in]  dupl_in   Array of duplicates in dbin (Dimension: ndupl)
 ** \param[in]  dupl_out  Array of duplicates in dbout (Dimension: ndupl)
 **
 *****************************************************************************/
static int st_external_vertype_define(SPDE_Mesh *s_mesh,
                                      int nvertex,
                                      int nbin,
                                      int nbout,
                                      int ndupl,
                                      int order,
                                      int *dupl_in,
                                      int *dupl_out)
{
  int vertype_loc, ijoint, ecr;
  Vertype *vertype;

  /* Initializations */

  s_mesh->vertype = vertype = vertype_manage(1, NULL, NULL, nvertex);
  if (s_mesh->vertype == nullptr) return (1);
  s_mesh->vercoloc = vercoloc_from_external(ndupl, dupl_in, dupl_out);
  if (s_mesh->vercoloc == nullptr) return (1);
  vertype->order = order;
  vertype->nb1 = 0;
  vertype->nb2 = 0;

  /* Dispatch according to the order */

  vertype_loc = VT_FREE;
  ecr = 0;
  if (order == 1)
  {

    /* Dbin followed by Dbout */

    if (S_DECIDE.flag_mesh_dbin)
    {
      for (int i = 0; i < nbin; i++)
      {
        ijoint = st_find_in_duplicate(i, ndupl, dupl_in);
        if (ijoint < 0)
        {
          vertype_loc = VT_HARD;
        }
        else
        {
          vertype_loc = VT_HARD | VT_OUTPUT;
        }
        vertype->vt[ecr++] = VT_INPUT | vertype_loc;
        vertype->nb1++;
      }
    }

    for (int i = 0; i < nbout; i++)
    {
      ijoint = st_find_in_duplicate(i, ndupl, dupl_out);
      if (ijoint >= 0) continue;
      vertype_loc = VT_FREE;
      vertype->vt[ecr++] = VT_OUTPUT | vertype_loc;
      vertype->nb2++;
    }
  }
  else
  {

    /* Dbout followed by Dbin */

    for (int i = 0; i < nbout; i++)
    {
      ijoint = st_find_in_duplicate(i, ndupl, dupl_out);
      if (ijoint < 0)
      {
        vertype_loc = VT_FREE;
      }
      else if (S_DECIDE.flag_mesh_dbin)
      {
        vertype_loc = VT_HARD | VT_INPUT;
      }
      vertype->vt[ecr++] = VT_OUTPUT | vertype_loc;
      vertype->nb1++;
    }

    if (S_DECIDE.flag_mesh_dbin)
    {
      for (int i = 0; i < nbin; i++)
      {
        ijoint = st_find_in_duplicate(i, ndupl, dupl_in);
        if (ijoint >= 0) continue;
        vertype_loc = VT_HARD;
        vertype->vt[ecr++] = VT_INPUT | vertype_loc;
        vertype->nb2++;
      }
    }
  }

  /* Optional printout */

  if (VERBOSE) st_vertype_print(vertype);

  return (0);
}

/****************************************************************************/
/*!
 **  Manage the contents of the External SPDE_Mesh structure used for
 **  storing external meshing information
 **
 ** \param[in]  mode      1 for storing; -1 for deallocating
 ** \param[in]  icov0     Rank of the current covariance (from 0 to 2)
 ** \param[in]  ndim      Space dimension
 ** \param[in]  ncorner   Number of vertices per element
 ** \param[in]  nvertex   Number of points
 ** \param[in]  nmesh     Number of meshes
 ** \param[in]  nbin      Number of vertices dedicated to the Input File
 ** \param[in]  nbout     Number of vertices dedicated to the Output File
 ** \param[in]  ndupl     Number of duplicates
 ** \param[in]  order     1: dbin then dbout; 2: dbout then dbin
 ** \param[in]  dupl_in   Array of duplicates in dbin (Dimension: ndupl)
 ** \param[in]  dupl_out  Array of duplicates in dbout (Dimension: ndupl)
 ** \param[in]  meshes    Array containing the meshes
 ** \param[in]  points    Array containing the vertex coordinates
 **
 *****************************************************************************/
int spde_external_mesh_define(int mode,
                              int icov0,
                              int ndim,
                              int ncorner,
                              int nvertex,
                              int nmesh,
                              int nbin,
                              int nbout,
                              int ndupl,
                              int order,
                              int *dupl_in,
                              int *dupl_out,
                              int *meshes,
                              double *points)
{
  int error, size;

  /* Initializations */

  error = 1;

  /* Dispatch */

  if (mode > 0)
  {
    S_EXTERNAL_MESH[icov0] = spde_mesh_manage(1, NULL);
    if (S_EXTERNAL_MESH[icov0] == nullptr) return (1);

    S_EXTERNAL_MESH[icov0]->ndim = ndim;
    S_EXTERNAL_MESH[icov0]->ncorner = ncorner;
    S_EXTERNAL_MESH[icov0]->nvertex = nvertex;
    S_EXTERNAL_MESH[icov0]->nmesh = nmesh;

    size = nvertex * ndim;
    S_EXTERNAL_MESH[icov0]->points = (double*) mem_alloc(sizeof(double) * size,
                                                         0);
    if (S_EXTERNAL_MESH[icov0]->points == nullptr) goto label_end;
    for (int i = 0; i < size; i++)
      S_EXTERNAL_MESH[icov0]->points[i] = points[i];

    size = nmesh * ncorner;
    S_EXTERNAL_MESH[icov0]->meshes = (int*) mem_alloc(sizeof(double) * size, 0);
    if (S_EXTERNAL_MESH[icov0]->meshes == nullptr) goto label_end;
    for (int i = 0; i < size; i++)
      S_EXTERNAL_MESH[icov0]->meshes[i] = meshes[i];

    if (st_external_vertype_define(S_EXTERNAL_MESH[icov0], nvertex, nbin, nbout,
                                   ndupl, order, dupl_in, dupl_out))
      goto label_end;
  }
  else
  {
    S_EXTERNAL_MESH[icov0] = spde_mesh_manage(-1, S_EXTERNAL_MESH[icov0]);
  }

  /* Set the error retun code */

  error = 0;

  label_end: if (error)
    spde_external_mesh_define(-1, icov0, ndim, ncorner, nvertex, nmesh, nbin,
                              nbout, ndupl, order, dupl_in, dupl_out, meshes,
                              points);
  return (error);
}

/****************************************************************************/
/*!
 **  Manage the contents of the External Q structure used for
 **  storing external information
 **
 ** \param[in]  mode      1 for storing; -1 for deallocating
 ** \param[in]  icov0     Rank of the current covariance (from 0 to 2)
 ** \param[in]  ndim      Space dimension
 ** \param[in]  nvertex   Number of points
 ** \param[in]  nmesh     Number of meshes
 ** \param[in]  nbin      Number of vertices dedicated to the Input File
 ** \param[in]  nbout     Number of vertices dedicated to the Output File
 ** \param[in]  ndupl     Number of duplicates
 ** \param[in]  order     1: dbin then dbout; 2: dbout then dbin
 ** \param[in]  dupl_in   Array of duplicates in dbin (Dimension: ndupl)
 ** \param[in]  dupl_out  Array of duplicates in dbout (Dimension: ndupl)
 ** \param[in]  A         Sparse matrix
 ** \param[in]  Q         Sparse matrix
 **
 *****************************************************************************/
int spde_external_AQ_define(int mode,
                            int icov0,
                            int ndim,
                            int nvertex,
                            int nmesh,
                            int nbin,
                            int nbout,
                            int ndupl,
                            int order,
                            int *dupl_in,
                            int *dupl_out,
                            cs *A,
                            cs *Q)
{
  int error;

  /* Initializations */

  error = 1;

  /* Dispatch */

  if (mode > 0)
  {
    S_EXTERNAL_MESH[icov0] = spde_mesh_manage(1, NULL);
    if (S_EXTERNAL_MESH[icov0] == nullptr) return (1);

    S_EXTERNAL_MESH[icov0]->ndim = ndim;
    S_EXTERNAL_MESH[icov0]->nvertex = nvertex;
    S_EXTERNAL_MESH[icov0]->nmesh = nmesh;

    S_EXTERNAL_Q[icov0] = cs_duplicate(Q);
    if (S_EXTERNAL_Q[icov0] == nullptr) return (1);

    S_EXTERNAL_A[icov0] = cs_duplicate(A);
    if (S_EXTERNAL_A[icov0] == nullptr) return (1);

    if (st_external_vertype_define(S_EXTERNAL_MESH[icov0], nvertex, nbin, nbout,
                                   ndupl, order, dupl_in, dupl_out))
      goto label_end;
  }
  else
  {
    S_EXTERNAL_Q[icov0] = cs_spfree(S_EXTERNAL_Q[icov0]);
    S_EXTERNAL_MESH[icov0] = spde_mesh_manage(-1, S_EXTERNAL_MESH[icov0]);
  }

  /* Set the error retun code */

  error = 0;

  label_end: if (error)
    spde_external_AQ_define(-1, icov0, ndim, nvertex, nmesh, nbin, nbout, ndupl,
                            order, dupl_in, dupl_out, A, Q);
  return (error);
}

/****************************************************************************/
/*!
 **  Assign fields of the SPDE_Mesh structure which would have been calculated
 **  elsewhere
 **
 ** \param[in]  ndim      Space dimension
 ** \param[in]  ncorner   Number of vertices per element
 ** \param[in]  nvertex   Number of points
 ** \param[in]  nmesh     Number of meshes
 ** \param[in]  meshes    Array containing the meshes
 ** \param[in]  points    Array containing the vertex coordinates
 ** \param[in]  verbose   Verbose flag
 **
 ** \param[in,out]  s_mesh   Pointer to SPDE_Mesh to be assigned
 **
 *****************************************************************************/
void spde_mesh_assign(SPDE_Mesh *s_mesh,
                      int ndim,
                      int ncorner,
                      int nvertex,
                      int nmesh,
                      int *meshes,
                      double *points,
                      int verbose)
{
  int number;
  static int debug = 0;

  s_mesh->ndim = ndim;
  s_mesh->ncorner = ncorner;
  s_mesh->nvertex = nvertex;
  s_mesh->nmesh = nmesh;

  // Arrays meshes and points must be copied in the OLD version
  // as the input argument may be a VectorDouble (freed automatically)

  number = nmesh * ncorner;
  s_mesh->meshes = (int*) mem_alloc(sizeof(int) * number, 1);
  for (int i = 0; i < number; i++)
    s_mesh->meshes[i] = meshes[i];
  number = ndim * nvertex;
  s_mesh->points = (double*) mem_alloc(sizeof(double) * number, 1);
  for (int i = 0; i < number; i++)
    s_mesh->points[i] = points[i];

  if (debug) st_print_mesh(s_mesh);
}

/****************************************************************************/
/*!
 **  Preparation using SPDE (for all GRF and COV)
 **
 ** \return  Error return code
 **
 ** \param[in]  dbin        Db structure for the conditioning data
 ** \param[in]  dbout       Db structure of the grid
 ** \param[in]  gext        Array of domain dilation
 ** \param[in]  s_option    SPDE_Option structure
 **
 *****************************************************************************/
int spde_prepar(Db *dbin,
                Db *dbout,
                const VectorDouble &gext,
                SPDE_Option &s_option)
{
  int error, nblin, flag_AQ_defined;

  /* Initializations */

  error = 1;
  st_calcul_init(st_get_ndim());

  /* Title (optional) */

  if (VERBOSE) st_title(0, 0, 1, "Preparing the Environment");

  /* Loop on the GRFs */

  for (int igrf = 0; igrf < st_get_number_grf(); igrf++)
  {
    st_set_current_igrf(igrf);

    /* Prepare the array of inverse of nugget sill matrices */

    if (S_DECIDE.flag_dbin && S_DECIDE.flag_several && st_is_model_nugget())
    {
      if (st_fill_Bnugget(dbin)) goto label_end;
    }

    /* Loop on the covariances */

    for (int icov = 0; icov < st_get_ncova(); icov++)
    {
      st_set_current_icov(icov);
      SPDE_Matelem &Matelem = spde_get_current_matelem(icov);
      flag_AQ_defined = st_is_external_AQ_defined(icov);

      /* Title (optional) */

      if (VERBOSE) st_title(1, 1, 1, "Preparing the Process");

      /* Initialize the structures */

      Matelem.s_mesh = spde_mesh_manage(1, NULL);
      if (Matelem.s_mesh == nullptr) goto label_end;

      /* Load the SPDE_Mesh structure */

      if (!flag_AQ_defined)
      {
        if (spde_mesh_load(Matelem.s_mesh, VERBOSE, dbin, dbout, gext,
                           s_option)) goto label_end;
      }

      // Locally convert from old to new Meshing

      MeshEStandard *amesh = new MeshEStandard();
      amesh->convertFromOldMesh(Matelem.s_mesh, 0);

      /* Load External Q (if any) */

      if (flag_AQ_defined)
      {
        if (spde_external_AQ_copy(Matelem, icov)) goto label_end;
      }

      /* Prepare the array of sparse matrices (without nugget effect) */

      if (S_DECIDE.flag_dbin && S_DECIDE.flag_several && !st_is_model_nugget())
      {
        if (st_fill_Bhetero(dbin, dbout)) goto label_end;
      }

      /* Preparation in non-stationary case */

      if (st_get_model()->isNoStat() && !flag_AQ_defined)
      {
        const ANoStat *nostat = st_get_model()->getNoStat();
        nostat->attachToMesh(amesh);
      }

      /* Prepare the projection matrix */

      if (S_DECIDE.flag_dbin)
      {
        if ((S_DECIDE.flag_several && !flag_AQ_defined) || !S_DECIDE.flag_mesh_dbin)
        {
          Matelem.Aproj = db_mesh_sparse(dbin, amesh, 0);
          if (Matelem.Aproj == nullptr) goto label_end;
          st_keypair_cs("Aproj", Matelem.Aproj, icov + 1, 0, 0, 0, 0);
        }
      }

      /* Prepare the kriging environment per structure */

      if (S_DECIDE.flag_dbin && S_DECIDE.flag_several)
      {
        if (st_fill_Isill()) goto label_end;
      }

      /* Prepare the simulation environment per structure */

      if (S_DECIDE.flag_case == CASE_SIMULATE)
      {
        if (st_fill_Csill()) goto label_end;
      }

      /* Build all relevant matrices */

      if (!flag_AQ_defined)
      {
        if (spde_build_matrices(st_get_model(), VERBOSE)) goto label_end;
      }

      /* Build additional matrices */

      if (S_DECIDE.flag_Q && S_DECIDE.flag_dbin)
      {
        if (st_build_QCov(Matelem)) goto label_end;
      }

      /* Partially free the SP_Mat structure */

      st_matelem_manage(0);

      /* Building simulation or Kriging environment */

      Matelem.qsimu = st_qsimu_manage(1, Matelem.s_mesh, NULL);
      if (Matelem.qsimu == nullptr) goto label_end;

      /* Prepare the Chebychev simulation environment */

      if (S_DECIDE.simu_cheb)
      {
        nblin = static_cast<int>(Calcul.blin.size());
        Matelem.s_cheb = spde_cheb_manage(1, VERBOSE, -0.5, nblin,
                                          Calcul.blin.data(), Matelem.S,
                                          NULL);
        if (Matelem.s_cheb == nullptr) goto label_end;
      }

      /* Verbose output (optional) */

      if (DEBUG && VERBOSE) st_matelem_print(icov);
    }
  }

  /* Set the error return code */

  error = 0;

  label_end: st_set_current_igrf(0);
  st_set_current_icov(0);
  return (error);
}

/****************************************************************************/
/*!
 **  Cleaning operation after SPDE
 **
 ** \return  Error return code
 **
 ** \param[in]  dbin        Db structure for the conditioning data
 ** \param[in]  dbout       Db structure of the grid
 ** \param[in]  gext        Array of domain dilation
 ** \param[in]  s_option    SPDE_Option structure
 **
 *****************************************************************************/
int spde_posterior(Db *dbin,
                   Db *dbout,
                   const VectorDouble &gext,
                   SPDE_Option &s_option)
{
  if (st_get_model()->isNoStat())
  {
    const ANoStat *nostat = st_get_model()->getNoStat();
    nostat->detachFromMesh();
  }
  return 0;
}

/****************************************************************************/
/*!
 **  Print the environment
 **
 ** \param[in]  dbout         Db output structure
 ** \param[in]  gext          Array of domain dilation
 **
 *****************************************************************************/
static void st_environ_print(const Db *dbout, const VectorDouble &gext)
{
  if (S_DECIDE.flag_case == CASE_KRIGING)
  {
    if (S_DECIDE.flag_est) message("Estimation\n");
    if (S_DECIDE.flag_std)
      message("Standard Deviation of the Estimation Error\n");
    if (st_get_filnug())
      message("Filtering Nugget effect (Sill=%lg)\n", st_get_nugget_sill(0, 0));
  }
  else
  {
    if (S_DECIDE.flag_dbin)
      message("Conditional Simulation\n");
    else
      message("Non-Conditional Simulation\n");
  }
  if (S_DECIDE.flag_onechol)
    message("- Single Cholesky option: ON\n");
  else
    message("- Single Cholesky option: OFF\n");

  if (S_DECIDE.flag_filnug) message("- Filter component must be filtered\n");

  if (S_DECIDE.flag_gibbs) message("- Gibbs iterations\n");

  if (is_grid(dbout) && !gext.empty())
  {
    message("- The resulting Grid is dilated: %lf", gext[0]);
    for (int idim = 1; idim < dbout->getNDim(); idim++)
      message(" * %lf", gext[idim]);
    message("\n");
  }
}

/****************************************************************************/
/*!
 **  True if a point belongs to a mesh and compute barycentric coordinates
 **
 ** \returns 1 if the point belongs to the mesh; 0 otherwise
 **
 ** \param[in]  amesh     MeshEStandard structure
 ** \param[in]  coor      Array of target coordinates
 ** \param[in]  imesh     Mesh Index
 ** \param[in]  units     Array of the meshes
 **
 ** \param[out] cotes     Array of barycentric coordinates
 **
 *****************************************************************************/
static bool is_in_mesh(MeshEStandard *amesh,
                       double *coor,
                       int imesh,
                       double *units,
                       double *cotes)
{
  double mat[9], ratio, total;
  int ecr, ncorner, ndim, flag_sphere;

  /* Initializations */

  ncorner = amesh->getNApexPerMesh();
  ndim = amesh->getNDim();
  variety_query(&flag_sphere);

  /* Calculate the barycentric coordinates */

  if (flag_sphere)
  {
    double pts[3][2];
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 2; j++)
        pts[i][j] = amesh->getCoor(imesh, i, j);
    return (is_in_spherical_triangle_optimized(coor, pts[0], pts[1], pts[2],
                                               cotes));
  }
  else
  {
    total = 0.;
    for (int icorn = 0; icorn < ncorner; icorn++)
    {

      /* Fill the matrix for determinant */

      ecr = 0;
      for (int jcorn = 0; jcorn < ncorner; jcorn++)
      {
        if (icorn == jcorn) continue;
        for (int idim = 0; idim < ndim; idim++)
          mat[ecr++] = amesh->getCoor(imesh, jcorn, idim) - coor[idim];
      }
      ratio = ABS(matrix_determinant(ndim,mat)) / units[imesh] / FACDIM[ndim];
      if (ratio < 0 || ratio > 1) return (0);
      cotes[icorn] = ratio;
      total += ratio;
    }
    if (ABS(total - 1) > 1.e-3) return (0);
  }
  return (1);
}

/****************************************************************************/
/*!
 **  True if a (dilated) point intersects a mesh
 **
 ** \returns 1 if the intersection is not empty; 0 otherwise
 **
 ** \param[in]  s_mesh    SPDE_Mesh structure
 ** \param[in]  coor      Point coordinates (Dimension: ndim)
 ** \param[in]  caux      Working array (Dimension: ndim)
 ** \param[in]  ndim      Comon space dimension between coor and s_mesh
 ** \param[in]  imesh     Mesh Index
 ** \param[in]  radius    Dilation radius
 **
 ** \remarks This function must not be called if on a sphere
 **
 *****************************************************************************/
static bool is_in_mesh_neigh(SPDE_Mesh *s_mesh,
                             double *coor,
                             double *caux,
                             int ndim,
                             int imesh,
                             double radius)
{
  double top, bot, alpha, delta, ref;
  int ncorner;

  /* Initializations */

  ncorner = s_mesh->ncorner;

  /* Check that at least one mesh vertex belongs to the point neighborhood */

  for (int icorn = 0; icorn < ncorner; icorn++)
  {
    for (int idim = 0; idim < ndim; idim++)
      caux[idim] = POINTS(OLD_MESHES(imesh,icorn), idim);
    if (ut_distance(ndim, coor, caux) <= radius) return (1);
  }

  /* Project the point on each mesh edge */

  for (int icorn = 0; icorn < ncorner - 1; icorn++)
    for (int jcorn = icorn + 1; jcorn < ncorner; jcorn++)
    {

      // Get the coordinates of the projection

      top = bot = 0.;
      for (int idim = 0; idim < ndim; idim++)
      {
        delta = (POINTS(OLD_MESHES(imesh,jcorn), idim)
            - POINTS(OLD_MESHES(imesh,icorn), idim));
        top += delta * (POINTS(OLD_MESHES(imesh,icorn),idim) - coor[idim]);
        bot += delta * delta;
      }
      alpha = top / bot;

      // Check that the projection lie between two vertices

      if (alpha < 0 || alpha > 1) continue;

      // Exhibit the coordinates of the projection

      for (int idim = 0; idim < ndim; idim++)
      {
        ref = POINTS(OLD_MESHES(imesh,icorn), idim);
        caux[idim] = ref + alpha * (POINTS(OLD_MESHES(imesh,jcorn),idim) - ref);
      }

      // Calculate the distance between the projection and the point

      if (ut_distance(ndim, coor, caux) <= radius) return (1);
    }

  return (0);
}

/****************************************************************************/
/*!
 **  Define the container for each mesh
 **
 ** \return Pointer to the array containing the containers
 **
 ** \param[in]  amesh    MeshEStandard structure
 **
 ** \remarks This function allocated memory which must be freed
 ** \remarks by calling function
 **
 *****************************************************************************/
static double* st_get_containers(MeshEStandard *amesh)
{
  double *contain, value, vmin, vmax;
  int ncorner, ndim, nmesh;

  /* Initializations */

  ndim = amesh->getNDim();
  nmesh = amesh->getNMeshes();
  ncorner = amesh->getNApexPerMesh();

  /* Allocation */

  contain = (double*) mem_alloc(sizeof(double) * 2 * ndim * nmesh, 0);
  if (contain == nullptr) return (contain);

  /* Loop on the meshes */

  for (int imesh = 0; imesh < nmesh; imesh++)
  {

    /* Loop on the dimensions */

    for (int idim = 0; idim < ndim; idim++)
    {
      vmin = 1.e30;
      vmax = -1.e30;

      /* Loop on the corners */

      for (int icorn = 0; icorn < ncorner; icorn++)
      {
        value = amesh->getCoor(imesh, icorn, idim);
        if (value < vmin) vmin = value;
        if (value > vmax) vmax = value;
      }

      CONTAIN(imesh,idim,0) = vmin;
      CONTAIN(imesh,idim,1) = vmax;
    }
  }
  return (contain);
}

/****************************************************************************/
/*!
 **  True if a point belongs to parallelepiped containing the mesh (container)
 **
 ** \return 1 if the point belongs to the container; 0 otherwise
 **
 ** \param[in]  coor      Array of target coordinates
 ** \param[in]  ndim      Space dimension
 ** \param[in]  imesh     Mesh Index
 ** \param[in]  contain   Array of containers
 **
 ** \remarks If the array 'contain' does not exist, the check is not performed
 **
 *****************************************************************************/
static bool is_in_mesh_container(double *coor,
                                 int ndim,
                                 int imesh,
                                 double *contain)
{
  double center, vmin, vmax;

  /* Initializations */

  if (contain == nullptr) return (1);

  /* Loop on the space dimension */

  for (int idim = 0; idim < ndim; idim++)
  {
    center = coor[idim];
    vmin = CONTAIN(imesh, idim, 0);
    vmax = CONTAIN(imesh, idim, 1);

    /* Check the rejection */

    if ((center - vmin) * (center - vmax) > 0) return (0);
  }
  return (1);
}

/****************************************************************************/
/*!
 **  Returns the set of barycenter weights as a sparse matrix
 **
 ** \return Pointer to the newly created sparse matrix (or NULL)
 **
 ** \param[in]  db         Db structure
 ** \param[in]  amesh      MeshEStandard structure
 ** \param[in]  verbose    Verbose flag
 **
 ** \remarks In the returned array 'Aproj', the samples are ordered as they
 ** \remarks appear within the array vertype (rather than in their natural
 ** \remarks order in Db)
 **
 *****************************************************************************/
cs* db_mesh_sparse(Db *db, MeshEStandard *amesh, int verbose)
{
  double *units, *coor, *contain, *weight;
  int error, imesh, imesh0, ip, flag_sphere, found, iech, ip_max;
  cs *A, *Atriplet;

  /* Initializations */

  error = 1;
  units = nullptr;
  coor = nullptr;
  contain = nullptr;
  weight = nullptr;
  Atriplet = A = nullptr;
  int ncorner = amesh->getNApexPerMesh();
  int nvertex = amesh->getNApices();
  int nmesh = amesh->getNMeshes();
  int ndim = amesh->getNDim();
  variety_query(&flag_sphere);

  Atriplet = cs_spalloc(0, 0, 1, 1, 1);
  if (Atriplet == nullptr) goto label_end;

  /* Core allocation */

  coor = db_sample_alloc(db, ELoc::X);
  if (coor == nullptr) goto label_end;
  weight = (double*) mem_alloc(sizeof(double) * ncorner, 0);
  if (weight == nullptr) goto label_end;
  if (!flag_sphere)
  {
    contain = st_get_containers(amesh);
    if (contain == nullptr) goto label_end;
  }

  /* Calculate the mesh units */

  units = _spde_get_mesh_dimension(amesh);
  if (units == nullptr) goto label_end;

  /* Loop on the samples */

  imesh0 = ip_max = iech = 0;
  for (int jech = 0; jech < db->getSampleNumber(); jech++)
  {
    if (!db->isActive(jech)) continue;

    /* Identification */

    for (int idim = 0; idim < ndim; idim++)
      coor[idim] = db->getCoordinate(jech, idim);

    /* Loop on the meshes */

    found = -1;
    for (int jmesh = 0; jmesh < nmesh; jmesh++)
    {
      imesh = imesh0 + jmesh;
      if (imesh >= nmesh) imesh -= nmesh;
      if (!is_in_mesh_container(coor, ndim, imesh, contain)) continue;
      if (!is_in_mesh(amesh, coor, imesh, units, weight)) continue;

      /* Store the items in the sparse matrix */

      for (int icorn = 0; icorn < ncorner; icorn++)
      {
        ip = amesh->getApex(imesh, icorn);
        if (ip > ip_max) ip_max = ip;
        if (!cs_entry(Atriplet, iech, ip, weight[icorn])) goto label_end;
      }
      imesh0 = found = imesh;
      break;
    }

    /* Printout if a point does not belong to any mesh */

    if (found < 0)
    {
      messerr("Point %d does not belong to any mesh", jech + 1);
      for (int idim = 0; idim < ndim; idim++)
        messerr(" Coordinate #%d = %lf", idim + 1, coor[idim]);
    }
    iech++;
  }

  /* Add the extreme value to force dimension */

  if (ip_max < nvertex - 1)
  {
    if (!cs_entry(Atriplet, db->getActiveSampleNumber() - 1, nvertex - 1, 0.))
      goto label_end;
  }

  /* Convert the triplet into a sparse matrix */

  A = cs_triplet(Atriplet);

  /* Set the error return code */

  error = 0;

  label_end: Atriplet = cs_spfree(Atriplet);
  contain = (double*) mem_free((char* ) contain);
  units = (double*) mem_free((char* ) units);
  weight = (double*) mem_free((char* ) weight);
  coor = db_sample_free(coor);
  if (error) A = cs_spfree(A);
  return (A);
}

/****************************************************************************/
/*!
 **  Construct the array of coordinates of triangles in 3-D space
 **  by gluing the original 2-D coordinates (turned back in the original space)
 **  with the information provided by the third dimension array
 **
 ** \return  Array of 3-D coordinates
 **
 ** \param[in]  s_mesh    SPDE_Mesh structure
 ** \param[in]  zcur      Array of results
 **
 ** \param[out] np3d      Number of 3-D points
 **
 *****************************************************************************/
static double* st_get_coords_3D(SPDE_Mesh *s_mesh, double *zcur, int *np3d)
{
  double *p3d, *points, value;
  int error, ndim, np, nvertex;

  /* Initializations */

  error = 1;
  ndim = 3;
  *np3d = 0;
  p3d = nullptr;
  nvertex = s_mesh->nvertex;
  points = s_mesh->points;

  /* Core allocation */

  p3d = (double*) mem_alloc(sizeof(double) * ndim * nvertex, 0);
  if (p3d == nullptr) goto label_end;

  /* Load the array of 3-D coordinates */

  np = 0;
  for (int i = 0; i < nvertex; i++)
  {
    value = zcur[i];
    if (FFFF(value)) continue;
    p3d[3 * np + 0] = points[2 * i + 0];
    p3d[3 * np + 1] = points[2 * i + 1];
    p3d[3 * np + 2] = value;
    np++;
  }

  /* Reallocation (if necessary) */

  if (np != nvertex)
  {
    p3d = (double*) mem_realloc((char* ) p3d, sizeof(double) * ndim * np, 0);
    if (p3d == nullptr) goto label_end;
  }
  *np3d = np;

  /* Set the error return code */

  error = 0;

  label_end: if (error) p3d = (double*) mem_free((char* ) p3d);
  return (p3d);
}

/****************************************************************************/
/*!
 **  Free all memory used in SPDE
 **
 *****************************************************************************/
void spde_free_all(void)

{
  for (int igrf = 0; igrf < SPDE_MAX_NGRF; igrf++)
  {
    st_set_current_igrf(igrf);
    st_matelem_manage(-1);
    st_clean_Bnugget();
    st_clean_Bhetero();
  }
}

/****************************************************************************/
/*!
 **  Define the main options
 **
 ** \return Error return code
 **
 ** \param[in]  dbin          Pointer to the input Db
 ** \param[in]  dbout         Pointer to the output Db
 ** \param[in]  model1        Model structure (first)
 ** \param[in]  model2        Model structure (second)
 ** \param[in]  verbose       Verbose flag
 ** \param[in]  gext          Array of domain dilation
 ** \param[in]  mesh_dbin     1 if Input points must participate to meshing
 ** \param[in]  mesh_dbout    1 if Output points must participate to meshing
 ** \param[in]  flag_advanced 1 for advanced calculus (estimation or simulation)
 **                           0 if only matrices are required
 ** \param[in]  flag_est      1 for estimation
 ** \param[in]  flag_std      1 for standard deviation
 ** \param[in]  flag_gibbs    1 for Gibbs sampler
 ** \param[in]  flag_modif    1 for post-processing simulations
 **
 *****************************************************************************/
int spde_check(const Db *dbin,
               const Db *dbout,
               Model *model1,
               Model *model2,
               int verbose,
               const VectorDouble &gext,
               int mesh_dbin,
               int mesh_dbout,
               int flag_advanced,
               int flag_est,
               int flag_std,
               int flag_gibbs,
               int flag_modif)
{
  Model *models[2];
  int nlevels, ncova;

  st_environ_init();

  VERBOSE = verbose;
  models[0] = model1;
  models[1] = model2;

  FLAG_KEYPAIR = (int) get_keypone("SPDE_FLAG_KEYPAIR", 0);
  nlevels = (int) get_keypone("Multigrid_Number_Levels", 0);
  DEBUG = (int) get_keypone("SPDE_DEBUG", DEBUG);
  S_DECIDE.simu_chol = (int) get_keypone("Flag_Simu_Chol", 0);
  S_DECIDE.simu_cheb = !S_DECIDE.simu_chol;

  S_DECIDE.flag_dbin = (dbin != nullptr);
  S_DECIDE.flag_dbout = (dbout != nullptr);
  S_DECIDE.flag_mesh_dbin = mesh_dbin;
  S_DECIDE.flag_mesh_dbout = mesh_dbout;
  S_DECIDE.flag_est = flag_est;
  S_DECIDE.flag_std = flag_std;
  S_DECIDE.flag_gibbs = flag_gibbs;
  S_DECIDE.flag_mgrid = nlevels > 0;
  S_DECIDE.flag_several = 0;

  S_DECIDE.flag_case = 0;
  if (!flag_advanced)
    S_DECIDE.flag_case = CASE_MATRICES;
  else if (S_DECIDE.flag_est == 0 && S_DECIDE.flag_std == 0)
    S_DECIDE.flag_case = CASE_SIMULATE;
  else
    S_DECIDE.flag_case = CASE_KRIGING;
  S_DECIDE.flag_Q = 1;
  if (!S_DECIDE.flag_dbin && S_DECIDE.simu_cheb) S_DECIDE.flag_Q = 0;
  if (!flag_advanced) S_DECIDE.flag_Q = 1;
  S_DECIDE.flag_Q = (int) get_keypone("Flag_Q", S_DECIDE.flag_Q);
  S_DECIDE.flag_Qchol = (S_DECIDE.flag_case == CASE_SIMULATE && S_DECIDE.flag_Q
                         && S_DECIDE.simu_chol);
  S_DECIDE.flag_modif = (S_DECIDE.flag_case == CASE_SIMULATE && flag_modif);
  S_DECIDE.flag_onechol = (!S_DECIDE.flag_mgrid && S_DECIDE.simu_chol);
  S_DECIDE.flag_onechol = (int) get_keypone("Flag_OneChol",
                                            S_DECIDE.flag_onechol);
  if (!S_DECIDE.flag_dbin) S_DECIDE.flag_onechol = 0;
  if (S_DECIDE.flag_est) S_DECIDE.flag_onechol = 1;
  if (S_DECIDE.flag_onechol) S_DECIDE.flag_Qchol = 0;

  /* Checks */

  if (S_DECIDE.flag_case != CASE_SIMULATE && S_DECIDE.flag_gibbs)
  {
    messerr(
        "'flag_gibbs' requires simulation ('flag_est' and 'flag_std' must be FALSE)");
    return (1);
  }
  if (S_DECIDE.flag_case != CASE_SIMULATE && flag_modif)
  {
    messerr(
        "'flag_modif' is limited to simulations ('flag_est' and 'flag_std' must be FALSE)");
    return (1);
  }
  if (S_DECIDE.flag_case == CASE_KRIGING && !S_DECIDE.flag_dbin)
  {
    messerr("You need to define an input Db to perform Estimation");
    return (1);
  }

  S_ENV.ngrfs = 0;
  for (int igrf = 0; igrf < SPDE_MAX_NGRF; igrf++)
  {
    if (models[igrf] != nullptr)
    {
      st_set_current_igrf(igrf);
      if (st_check_model(dbin, dbout, models[igrf])) return (1);
      st_calcul_init(st_get_ndim());
      st_matelem_manage(1);
      ncova = st_get_ncova();

      for (int icov = 0; icov < ncova; icov++)
      {
        st_set_current_icov(icov);
        st_calcul_update();
        if (VERBOSE) st_print_all("Model (Stationary) Parameters");
      }
      S_ENV.ngrfs++;
    }
  }
  S_DECIDE.flag_Qchol = S_DECIDE.flag_Qchol
      || (S_DECIDE.flag_case == CASE_MATRICES && S_DECIDE.flag_std);
  if (S_DECIDE.flag_std && st_get_nvar() > 1)
  {
    messerr(
        "Calculation of Kriging Variance is incompatible with Multivariate");
    return (1);
  }

  /* Optional printout */

  if (verbose)
  {
    st_title(0, 0, 1, "Environment for SPDE processing");
    message("Space Dimension          = %d\n", st_get_ndim());
    message("Number of variables      = %d\n", st_get_nvar());
    message("Presence of an input Db  = %d\n", S_DECIDE.flag_dbin);
    message("Presence of an output Db = %d\n", S_DECIDE.flag_dbout);
    message("Calculate estimation     = %d\n", S_DECIDE.flag_est);
    message("Calculate st. deviation  = %d\n", S_DECIDE.flag_std);
    message("Perform estimation       = %d\n",
            S_DECIDE.flag_case == CASE_KRIGING);
    message("Perform simulations      = %d\n",
            S_DECIDE.flag_case == CASE_SIMULATE);
    message("Perform gibbs sampler    = %d\n", S_DECIDE.flag_gibbs);
    message("Post-process simulation  = %d\n", S_DECIDE.flag_modif);
    if (S_DECIDE.flag_case != CASE_MATRICES)
    {
      message("Cholesky Simulations     = %d\n", S_DECIDE.simu_chol);
      message("Chebychev Simulations    = %d\n", S_DECIDE.simu_cheb);
    }
    message("Build Q                  = %d\n", S_DECIDE.flag_Q);
    message("Perform Cholesky of Q    = %d\n", S_DECIDE.flag_Qchol);

    st_environ_print(dbout, gext);
  }

  S_DECIDE.flag_Qchol = S_DECIDE.flag_Qchol
      || (S_DECIDE.flag_case == CASE_KRIGING && st_get_filnug());

  return (0);
}

/****************************************************************************/
/*!
 **  Perform Kriging using SPDE on the set of constructed 2-D triangles
 **
 ** \return  Error return code
 **
 ** \param[in]  dbin        Db input structure
 ** \param[in]  model       Model structure
 ** \param[in]  s_option    SPDE_Option structure
 ** \param[in]  verbose     Verbose option
 **
 ** \param[out] nmesh_arg    Number of triangles generated
 ** \param[out] nvertex_arg  Number of vertices
 ** \param[out] meshes_arg   Array of triangulate vertex indices
 ** \param[out] points_arg   Array containing the 2-D vertices
 **
 ** \remark It is the responsability of the calling function to free
 ** \remark the returned arrays 'points_arg' and 'meshes_arg'.
 **
 ** \remark This function assumes that there is only one Model defined
 ** \remark Hence the test on the number of models
 **
 *****************************************************************************/
int kriging2D_spde(Db *dbin,
                   Model *model,
                   SPDE_Option &s_option,
                   int verbose,
                   int *nmesh_arg,
                   int *nvertex_arg,
                   int **meshes_arg,
                   double **points_arg)
{
  int error, ncur, size, ndata, nvar, ncova;
  double *zcur, *work, *data;
  SPDE_Mesh *s_mesh;

  /* Initializations */

  error = 1;
  zcur = work = data = nullptr;
  *nmesh_arg = *nvertex_arg = 0;
  *meshes_arg = nullptr;
  *points_arg = nullptr;

  /* Preliminary checks */

  if (spde_check(dbin, NULL, model, NULL, verbose, VectorDouble(), 1, 1, 1, 1,
                 0, 0, 0)) goto label_end;
  if (st_get_number_grf() != 1)
  {
    messerr("This function should be called in the case of a single Model");
    messerr("In your case: %d\n", st_get_number_grf());
    goto label_end;
  }

  /* Preliminary checks */

  if (model->getDimensionNumber() != 2)
  {
    messerr("This application is restricted to the 2-D case (ndim=%d)",
            model->getDimensionNumber());
    goto label_end;
  }

  /* Prepare all material */

  if (spde_prepar(NULL, dbin, VectorDouble(), s_option)) goto label_end;
  st_set_current_igrf(0);
  {
    SPDE_Matelem &Matelem = spde_get_current_matelem(-1);
    s_mesh = Matelem.s_mesh;

    /* Core allocation */

    nvar = st_get_nvar();
    ncova = st_get_ncova_max();
    ndata = dbin->getActiveSampleNumber();
    ncur = s_mesh->nvertex;
    zcur = (double*) mem_alloc(sizeof(double) * ncur * nvar, 0);
    if (zcur == nullptr) goto label_end;
    work = (double*) mem_alloc(sizeof(double) * ncur, 0);
    if (work == nullptr) goto label_end;
    data = (double*) mem_alloc(sizeof(double) * ndata * nvar, 0);
    if (data == nullptr) goto label_end;

    /* Load the data */

    st_init_array(ncova, nvar, ncur, 1, zcur);
    st_load_data(s_mesh, dbin, NULL, s_option, -1, data, zcur);

    /* for estimation */

    if (st_get_filnug())
    {
      if (st_filter(work, zcur)) goto label_end;
    }
    else
    {
      if (st_kriging(s_mesh, data, zcur)) goto label_end;
    }
  }

  /* Create the returned array */

  *points_arg = st_get_coords_3D(s_mesh, zcur, nvertex_arg);
  size = s_mesh->nmesh * s_mesh->ncorner;
  *meshes_arg = (int*) mem_alloc(sizeof(int) * size, 0);
  if (*meshes_arg == nullptr) goto label_end;
  (void) memcpy((char*) *meshes_arg, (char*) s_mesh->meshes,
                sizeof(int) * size);
  *nmesh_arg = s_mesh->nmesh;

  /* Cleaning procedure */

  spde_posterior(NULL, dbin, VectorDouble(), s_option);

  /* Set the error code */

  error = 0;

  label_end: zcur = (double*) mem_free((char* ) zcur);
  work = (double*) mem_free((char* ) work);
  data = (double*) mem_free((char* ) data);
  return (error);
}

/****************************************************************************/
/*!
 **  Perform Estimation / Simulations using SPDE
 **
 ** \return  Error return code
 **
 ** \param[in]  dbin          Db input structure
 ** \param[in]  dbout         Db output structure
 ** \param[in]  model         Model structure
 ** \param[in]  gext          Array of domain dilation
 ** \param[in]  s_option      SPDE_Option structure
 ** \param[in]  mesh_dbin     1 if Data Samples belong to meshing vertices
 ** \param[in]  mesh_dbout    1 if Target Nodes belong to meshing vertices
 ** \param[in]  seed          Seed value for the random number generator
 ** \param[in]  nbsimu        Number of simulations
 ** \param[in]  ngibbs_burn   Number of iterations (Burning step)
 ** \param[in]  ngibbs_iter   Maximum number of iterations
 ** \param[in]  ngibbs_int    Number of iterations internal to Gibbs (SPDE)
 ** \param[in]  flag_est      1 for estimation
 ** \param[in]  flag_std      1 for standard deviation
 ** \param[in]  flag_gibbs    1 if the iterative Gibbs method must be used
 ** \param[in]  flag_modif    1 if the simulations must be transformed
 ** \param[in]  verbose       1 for a verbose processing
 **
 ** \remarks  If the number of simulations 'nbsimu' is set to 0,
 ** \remarks  the simulation algorithm is turned into a kriging one
 **
 *****************************************************************************/
int spde_f(Db *dbin,
           Db *dbout,
           Model *model,
           const VectorDouble &gext,
           SPDE_Option &s_option,
           int mesh_dbin,
           int mesh_dbout,
           int seed,
           int nbsimu,
           int ngibbs_burn,
           int ngibbs_iter,
           int ngibbs_int,
           int flag_est,
           int flag_std,
           int flag_gibbs,
           int flag_modif,
           int verbose)
{
  int error, iad, nvar, nv_krige;

  /* Initializations */

  error = 1;

  if (spde_check(dbin, dbout, model, NULL, verbose, gext, mesh_dbin, mesh_dbout,
                 1, flag_est, flag_std, flag_gibbs, flag_modif)) return (1);
  simu_define_func_transf(NULL);
  simu_define_func_update(simu_func_continuous_update);
  simu_define_func_scale(simu_func_continuous_scale);
  nvar = st_get_nvar();

  /* Preliminary checks */

  nv_krige = 0;
  if (S_DECIDE.flag_case == CASE_KRIGING)
  {
    if (S_DECIDE.flag_est) nv_krige += nvar;
    if (S_DECIDE.flag_std) nv_krige += nvar;
  }

  /* Initial assignments */

  law_set_random_seed(seed);
  if (S_DECIDE.flag_est || S_DECIDE.flag_std) nbsimu = 0;

  /* Add the attributes */

  if (S_DECIDE.flag_case == CASE_SIMULATE)
  {
    if (!S_DECIDE.flag_modif)
    {
      if (db_locator_attribute_add(dbout, ELoc::SIMU, MAX(1,nbsimu) * nvar, 0,
                                   0., &iad)) goto label_end;
    }
    else
    {
      if (db_locator_attribute_add(dbout, ELoc::SIMU, nvar, 0, 0., &iad))
        goto label_end;
      if (db_locator_attribute_add(dbout, ELoc::Z, 2 * nvar, 0, 0., &iad))
        goto label_end;
    }
  }
  else
  {
    if (db_locator_attribute_add(dbout, ELoc::Z, nv_krige, 0, 0., &iad))
      goto label_end;
  }

  /* Prepare all the material */

  if (spde_prepar(dbin, dbout, gext, s_option)) goto label_end;

  /* Perform the simulation */

  if (spde_process(dbin, dbout, s_option, nbsimu, ngibbs_burn, ngibbs_iter,
                   ngibbs_int)) goto label_end;

  /* Garbage collector */

  spde_posterior(dbin, dbout, gext, s_option);

  /* Set the error return code */

  error = 0;

  label_end: if (S_DECIDE.flag_modif) dbout->deleteFieldByLocator(ELoc::SIMU);
  return (error);
}

/****************************************************************************/
/*!
 **  Perform the product of x by Q using the blin decomposition
 **
 ** \param[in]  nblin     Number of blin coefficients
 ** \param[in]  blin      Array of coefficients for Linear combinaison
 ** \param[in]  S         Shift operator
 ** \param[in]  Lambda    Vector Lambda
 ** \param[in]  TildeC    Vector TildeC
 ** \param[in]  x         Input array
 **
 ** \param[out] y         Output array
 **
 *****************************************************************************/
static void st_product_Q(int nblin,
                         double *blin,
                         cs *S,
                         const VectorDouble &Lambda,
                         const VectorDouble &TildeC,
                         double *x,
                         double *y)
{
  double *x1, *x2;
  int n;

  // Initializations

  n = S->n;

  // Core allocation 

  x1 = (double*) mem_alloc(sizeof(double) * n, 1);
  x2 = (double*) mem_alloc(sizeof(double) * n, 1);

  for (int i = 0; i < n; i++)
    y[i] = 0.;
  for (int i = 0; i < n; i++)
    x1[i] = x[i] * sqrt(TildeC[i]);

  for (int ilin = 0; ilin < nblin; ilin++)
  {
    for (int i = 0; i < n; i++)
      y[i] += blin[ilin] * x1[i];
    cs_vecmult(S, n, x1, x2);
    for (int i = 0; i < n; i++)
      x1[i] = x2[i];
  }

  for (int i = 0; i < n; i++)
    y[i] *= pow(Lambda[i] / sqrt(TildeC[i]), 2.) * sqrt(TildeC[i]);

  x1 = (double*) mem_free((char* ) x1);
  x2 = (double*) mem_free((char* ) x2);
}

/****************************************************************************/
/*!
 **  Perform the product of a vector by the inverse of the power
 **  of a sparse matrix using the Chebychev Polynomial procedure
 **
 ** \return Error return code
 **
 ** \param[in]  nblin     Number of blin coefficients
 ** \param[in]  blin      Array of coefficients for Linear combinaison
 ** \param[in]  S         Shift operator
 ** \param[in]  Lambda    Vector Lambda
 ** \param[in]  TildeC    Vector TildeC
 ** \param[in]  power     Parameter used in the Chebychev approximation
 ** \param[in]  x         Input array
 **
 ** \param[out] y         Output array
 **
 *****************************************************************************/
int spde_eval(int nblin,
              double *blin,
              cs *S,
              const VectorDouble &Lambda,
              const VectorDouble &TildeC,
              double power,
              double *x,
              double *y)
{
  Cheb_Elem *cheb_elem;
  int error, n;

  /* Initializations */

  error = 1;
  cheb_elem = nullptr;
  n = S->n;
  if (power != 1.0 && power != -1.0 && power != -0.5)
  {
    messerr("Invalid value for the 'power' argument (%lf)", power);
    messerr("It should be either -1, -0.5 or 1");
    goto label_end;
  }

  // Dispatch 

  if (power == 1.)
  {
    st_product_Q(nblin, blin, S, Lambda, TildeC, x, y);
  }
  else
  {
    /* Pre-processing for the case power=-1 */

    if (power == -1) for (int i = 0; i < n; i++)
      x[i] /= sqrt(TildeC[i]);

    /* Create the Cheb_Elem structure */

    cheb_elem = spde_cheb_manage(1, VERBOSE, power, nblin, blin, S, NULL);
    if (cheb_elem == nullptr) goto label_end;

    /* Operate the Chebychev polynomials */

    if (spde_chebychev_operate(S, cheb_elem, Lambda, x, y)) goto label_end;

    /* Post-processing for the case power=-1 */

    if (power == -1) for (int i = 0; i < n; i++)
      y[i] *= sqrt(TildeC[i]);
  }

  /* Set the error return code */

  error = 0;

  label_end: cheb_elem = spde_cheb_manage(-1, 0, 0, 0, NULL, NULL, cheb_elem);
  return (error);
}

/****************************************************************************/
/*!
 **  Check the pinchout variable
 **
 ** \return Error returned code
 **
 ** \param[in]  dbgrid      Grid structure
 ** \param[in]  icol_pinch  Pointer to the pinchout variabe
 **
 *****************************************************************************/
static int st_m2d_check_pinchout(Db *dbgrid, int icol_pinch)
{
  int nech, error;
  double *tab;

  // Preliminary checks

  if (dbgrid == nullptr) return (0);
  if (icol_pinch < 0) return (0);

  // Initializations

  error = 1;
  nech = dbgrid->getSampleNumber();
  tab = db_vector_alloc(dbgrid);
  if (tab == nullptr) return (1);
  if (db_vector_get_att(dbgrid, icol_pinch, tab)) goto label_end;

  // Check that values are within [0,1] interval

  for (int iech = 0; iech < nech; iech++)
  {
    if (!dbgrid->isActive(iech)) continue;
    if (FFFF(tab[iech])) continue;
    if (tab[iech] < 0 || tab[iech] > 1)
    {
      messerr("Pinchout variable should lie in [0,1]");
      messerr("At grid node %d/%d, the value is %lf", iech + 1, nech,
              tab[iech]);
      goto label_end;
    }
  }

  // Set the error return code

  error = 0;

  label_end: tab = db_vector_free(tab);
  return (error);
}

/****************************************************************************/
/*!
 **  Get the elevation within bounds
 **
 ** \return The value assigned to this inequality
 **
 ** \param[in]  m2denv      M2D_Environ structure
 ** \param[in]  nlayer      Number of layers
 ** \param[in]  ilayer      Rank of the target layer
 ** \param[in]  lower       Lower bound
 ** \param[in]  upper       Upper bound
 **
 *****************************************************************************/
static double st_m2d_draw_elevation(M2D_Environ *m2denv,
                                    int nlayer,
                                    int ilayer,
                                    double lower,
                                    double upper)
{
  double value, lowloc, upploc, mean, stdv;

  mean = m2denv->zmean;
  stdv = m2denv->zstdv;
  lowloc = lower;
  upploc = upper;
  value = 0.;
  if (!FFFF(lower)) lowloc = (lower - mean) / stdv;
  if (!FFFF(upper)) upploc = (upper - mean) / stdv;
  if (!FFFF(lower) && !FFFF(upper))
    value = mean + stdv * law_gaussian_between_bounds(lowloc, upploc);
  else if (FFFF(lower) && FFFF(upper))
    value = mean;
  else if (FFFF(lower))
    value = mean + stdv * law_gaussian_between_bounds(TEST, upploc);
  else if (FFFF(upper))
    value = mean + stdv * law_gaussian_between_bounds(lowloc, TEST);

  return (value);
}

/****************************************************************************/
/*!
 **  Print (concatenate) the printout of an interval
 **
 ** \param[in]  title       Optional title
 ** \param[in]  lower       Lower bound or FFFF
 ** \param[in]  upper       Upper bound or FFFF
 ** \param[in]  tail        0: blank character; 1: "\n" character
 **
 ** \remarks The printed string starts with a blank character
 ** \remarks It ends with either a blank or a <SR/LF> character (see 'tail')
 **
 *****************************************************************************/
static void st_print_concatenate_interval(const char *title,
                                          double lower,
                                          double upper,
                                          int tail)
{
  if (title != NULL) message("%s", title);
  message(" [");
  if (FFFF(lower))
    message("    NA");
  else
    message("%6.2lf", lower);
  message(" ; ");
  if (FFFF(upper))
    message("    NA");
  else
    message("%6.2lf", upper);
  message("]");

  if (tail == 0)
    message(" ");
  else
    message("\n");
}

/****************************************************************************/
/*!
 **  Print the constraints information for a single point
 **
 ** \param[in]  ilayer      Rank of the layer
 ** \param[in]  iech        Rank of the sample
 ** \param[in]  value       Current value
 ** \param[in]  drift       Drift value (or TEST)
 ** \param[in]  vgaus       Current Gaussian value (or TEST)
 ** \param[in]  lower       Lower bound or FFFF
 ** \param[in]  upper       Upper bound or FFFF
 **
 *****************************************************************************/
static void st_print_constraints_per_point(int ilayer,
                                           int iech,
                                           double value,
                                           double drift,
                                           double vgaus,
                                           double lower,
                                           double upper)
{
  message("Sample (%d) - Layer (%3d) in", iech + 1, ilayer + 1);
  st_print_concatenate_interval(NULL, lower, upper, 0);
  if (!FFFF(drift)) message("- Drift=%8.3lf ", drift);
  if (!(FFFF(value) && FFFF(vgaus)))
  {
    message("->");
    if (FFFF(value))
      message("       NA");
    else
      message(" %8.4lf", value);
    if (!FFFF(vgaus)) message(" (Gaus=%8.4lf)", vgaus);
  }
  message("\n");
}

/****************************************************************************/
/*!
 **  Check the validity of the Mean and Variance values
 **
 ** \return Error return code
 **
 ** \param[in]  db            Db structure containing the constraints
 ** \param[in]  ilayer        Rank of the layer of interest
 ** \param[in]  iech          Rank of the sample of interest
 ** \param[in]  flag_positive Positivity check
 ** \param[in]  flag_verbose  Verbose output
 ** \param[in]  M             Value for the Mean
 ** \param[in]  S             Value for the Variance
 **
 *****************************************************************************/
static int st_check_validity_MS(Db *db,
                                int ilayer,
                                int iech,
                                int flag_positive,
                                int flag_verbose,
                                double M,
                                double S)
{
  int error;
  static double eps = 1.e-3;

  error = 0;
  if (FFFF(M) || FFFF(S)) error = 1;
  if (flag_positive)
  {
    if (M < eps || S < eps) error = 1;
  }
  if (error == 0) return (0);
  if (flag_verbose)
  {
    messerr("Error at Sample #%d/%d for Layer #%d", iech + 1,
            db->getSampleNumber(), ilayer + 1);
    if (FFFF(M))
      messerr("- Mean is undefined");
    else
    {
      if (flag_positive && M < eps)
        messerr("- Mean has a too small value (%lf)", M);
    }
    if (FFFF(S))
      messerr("- Variance is undefined");
    else
    {
      if (flag_positive && S < eps)
        messerr("- Variance has a too small value (%lf)", S);
    }
  }
  return (1);
}

/****************************************************************************/
/*!
 **  Returns the value of the drift increment at a sample (mean)
 **
 ** \return The mean value or TEST value
 **
 ** \param[in]  m2denv      M2D_Environ structure
 ** \param[in]  db          Db structure containing the constraints
 ** \param[in]  type        1 for the constraining Db
 **                         2 for the grid output Db
 ** \param[in]  ilayer      Rank of the layer of interest
 ** \param[in]  iech        Rank of the sample of interest
 **
 *****************************************************************************/
static double st_m2d_get_M(M2D_Environ *m2denv,
                           Db *db,
                           int type,
                           int ilayer,
                           int iech)
{
  double value;
  int iatt;

  if (type == 1)
    iatt = m2denv->iatt_fd;
  else
    iatt = m2denv->iatt_fg;
  value = db->getArray(iech, iatt + ilayer);
  return (value);
}

/****************************************************************************/
/*!
 **  Returns the value of the gaussian standard deviation
 **
 ** \return The mean value or TEST value
 **
 ** \param[in]  m2denv      M2D_Environ structure
 ** \param[in]  db          Db structure containing the constraints
 ** \param[in]  type        1 for the constraining Db
 **                         2 for the grid output Db
 ** \param[in]  ilayer      Rank of the layer of interest
 ** \param[in]  iech        Rank of the sample of interest
 **
 *****************************************************************************/
static double st_m2d_get_S(M2D_Environ *m2denv,
                           Db *db,
                           int type,
                           int ilayer,
                           int iech)
{
  double value;

  value = m2denv->ystdv;
  return (value);
}

/****************************************************************************/
/*!
 **  At a point, returns the external drift increment from previous layer
 **
 ** \return The external drift increment
 **
 ** \param[in]  m2denv      M2D_Environ structure
 ** \param[in]  db          Db structure containing the constraints
 ** \param[in]  ilayer0     Rank of the layer of interest
 ** \param[in]  iech0       Rank of the sample of interest
 **
 *****************************************************************************/
static double st_m2d_external_drift_increment(M2D_Environ *m2denv,
                                              Db *db,
                                              int ilayer0,
                                              int iech0)
{
  double value, previous;

  value = db->getExternalDrift(iech0, ilayer0);
  if (FFFF(value)) return (TEST);
  if (ilayer0 > 1)
    previous = db->getExternalDrift(iech0, ilayer0 - 1);
  else
    previous = m2denv->dmini;
  if (FFFF(previous)) return (TEST);
  value -= previous;
  return (value);
}

/****************************************************************************/
/*!
 **  Returns the value of the drift contribution at a sample
 **  This value is a weighted combinaison of constant and external drift term
 **
 ** \return The drift interval value
 **
 ** \param[in]  m2denv      M2D_Environ structure
 ** \param[in]  db          Db structure containing the constraints
 ** \param[in]  ilayer0     Rank of the layer of interest
 ** \param[in]  iech0       Rank of the sample of interest
 **
 *****************************************************************************/
static double st_m2d_get_drift(M2D_Environ *m2denv,
                               Db *db,
                               int ilayer0,
                               int iech0)
{
  double coeff, value, drift;

  coeff = DCOEF(ilayer0);
  if (m2denv->flag_ed)
    drift = st_m2d_external_drift_increment(m2denv, db, ilayer0, iech0);
  else
    drift = 1.;
  if (FFFF(drift)) return (TEST);
  value = coeff * drift;
  return (value);
}

/****************************************************************************/
/*!
 **  Calculate the drift increment in a Db
 **
 ** \param[in]  m2denv      M2D_Environ structure
 ** \param[in]  nlayer      Number of layers
 ** \param[in]  icol_pinch  Pointer to the pinchout variabe
 ** \param[in]  db          Db structure
 ** \param[in]  iatt        Pointer to the drift vector
 **
 *****************************************************************************/
static void st_m2d_set_M(M2D_Environ *m2denv,
                         int nlayer,
                         int icol_pinch,
                         Db *db,
                         int iatt)
{
  double drift;

  for (int ilayer = 0; ilayer < nlayer; ilayer++)
  {
    for (int iech = 0; iech < db->getSampleNumber(); iech++)
    {
      if (db->isActive(iech))
      {
        drift = st_m2d_get_drift(m2denv, db, ilayer, iech);
        if (!FFFF(drift) && ilayer > 0 && icol_pinch >= 0)
          drift *= db->getArray(iech, icol_pinch);
      }
      else
      {
        drift = TEST;
      }
      db->setArray(iech, iatt + ilayer, drift);
    }
  }
}

/****************************************************************************/
/*!
 **  Locally migrate the pinchout distance from grid to point
 **
 ** \return  Address of the newly added vector in 'dbc'
 **
 ** \param[in]  dbout       Db output structure
 ** \param[in]  dbc         Db constraints structure
 ** \param[in]  icol_pinch  Pointer to the pinchout variabe
 **
 *****************************************************************************/
static int st_m2d_migrate_pinch_to_point(Db *dbout, Db *dbc, int icol_pinch)
{
  double *tab;
  int iptr, error;
  VectorInt cols(1);
  cols[0] = icol_pinch;

  // Initializations

  error = 1;
  iptr = -1;
  tab = nullptr;
  if (dbout == nullptr) return (0);
  if (icol_pinch < 0) return (0);

  // Add an attribute

  iptr = dbc->addFields(1, TEST);
  if (iptr < 0) goto label_end;

  // Core allocation

  tab = db_vector_alloc(dbc);
  if (tab == nullptr) goto label_end;

  // Migrate information from grid to point

  if (migrateByAttribute(dbout, dbc, cols, 0, VectorDouble(), false, false))
    goto label_end;

  // Store the resulting array in the file

  dbc->setFieldByAttribute(tab, iptr);

  // Set the error returned code

  error = 0;

  label_end: if (error && iptr >= 0) dbc->deleteFieldByAttribute(iptr);
  tab = (double*) mem_free((char* ) tab);
  return (iptr);
}

/****************************************************************************/
/*!
 **  Calculate and store drift value per point in constraints and output Db
 **  Check the validity of the drift at points
 **
 ** \return  Error returned code
 **
 ** \param[in]  m2denv      M2D_Environ structure
 ** \param[in]  mode        1 adding; -1 deleting
 ** \param[in]  nlayer      Number of layers
 ** \param[in]  icol_pinch  Pointer to the pinchout variabe
 ** \param[in]  dbc         Db constraints structure
 ** \param[in]  dbout       Db output structure
 **
 *****************************************************************************/
static int st_m2d_drift_inc_manage(M2D_Environ *m2denv,
                                   int mode,
                                   int nlayer,
                                   int icol_pinch,
                                   Db *dbc,
                                   Db *dbout)
{
  double M, S;
  int iptr;

  /* Initializations */

  if (m2denv == (M2D_Environ*) NULL) return (1);
  iptr = -1;

  /* Dispatch */

  if (mode > 0)
  {

    /* Identify the drift at the constraining samples */

    m2denv->iatt_fd = dbc->addFields(nlayer, TEST);
    if (m2denv->iatt_fd < 0) return (1);

    /* If pinch-out is defined, interpolate it at well data */

    iptr = st_m2d_migrate_pinch_to_point(dbout, dbc, icol_pinch);
    st_m2d_set_M(m2denv, nlayer, iptr, dbc, m2denv->iatt_fd);
    if (iptr >= 0) dbc->deleteFieldByAttribute(iptr);

    /* Check validity of drift at data points */

    for (int iech = 0; iech < dbc->getSampleNumber(); iech++)
    {
      if (!dbc->isActive(iech)) continue;
      for (int ilayer = 0; ilayer < nlayer; ilayer++)
      {
        M = st_m2d_get_M(m2denv, dbc, 1, ilayer, iech);
        S = st_m2d_get_S(m2denv, dbc, 1, ilayer, iech);
        if (st_check_validity_MS(dbc, ilayer, iech, 1, 1, M, S)) return (1);
      }
    }

    /* Identify the drift at the target grid nodes */

    m2denv->iatt_fg = dbout->addFields(nlayer, TEST);
    if (m2denv->iatt_fg < 0) return (1);
    st_m2d_set_M(m2denv, nlayer, icol_pinch, dbout, m2denv->iatt_fg);
  }
  else
  {

    /* Deleting the drift at the constraining samples */

    if (m2denv->iatt_fd >= 0)
      (void) db_attribute_del_mult(dbc, m2denv->iatt_fd, nlayer);

    /* Deleting the drift at the target grid */

    if (m2denv->iatt_fg >= 0)
      (void) db_attribute_del_mult(dbout, m2denv->iatt_fg, nlayer);
  }
  return (0);
}

/****************************************************************************/
/*!
 **  Calculate global statistics on elevations
 **
 ** \param[in]  m2denv      M2D_Environ structure
 ** \param[in]  dbin        Db input structure
 ** \param[in]  nlayer      Number of layers
 ** \param[in]  verbose     Verbose flag
 **
 *****************************************************************************/
static void st_m2d_stats_init(M2D_Environ *m2denv,
                              Db *dbin,
                              int nlayer,
                              int verbose)
{
  int nech;
  double lower, upper, nb, mm, vv, mini, maxi, delta;
  static double percent = 0.05;

  /* Initializations */

  nech = dbin->getSampleNumber();
  nb = mm = vv = 0.;
  mini = 1.e30;
  maxi = -1.e30;

  /* Loop on the layers */

  for (int ilayer = 0; ilayer < nlayer; ilayer++)
  {

    /* Loop on the samples */

    for (int iech = 0; iech < nech; iech++)
    {
      if (!dbin->isActive(iech)) continue;
      lower = dbin->getLowerBound(iech, ilayer);
      upper = dbin->getUpperBound(iech, ilayer);

      // Process the minimum bound 

      if (!FFFF(lower))
      {
        nb += 1.;
        mm += lower;
        vv += lower * lower;
        if (lower < mini) mini = lower;
        if (lower > maxi) maxi = lower;
      }

      // Process the maximum bound 

      if (!FFFF(upper))
      {
        nb += 1.;
        mm += upper;
        vv += upper * upper;
        if (upper < mini) mini = upper;
        if (upper > maxi) maxi = upper;
      }
    }
  }

  /* Normation */

  if (nb > 0)
  {
    mm /= nb;
    vv = vv / nb - mm * mm;
  }
  else
  {
    mm = 0.;
    vv = 1.;
    mini = -0.5;
    maxi = 0.5;
  }

  delta = maxi - mini;
  if (delta <= 0) delta = ABS(mm) / 10.;
  if (delta <= 0) delta = 1.;
  m2denv->zmean = mm;
  m2denv->zeps = ABS(mm) / 1.e4;
  m2denv->zstdv = (vv > 0) ? sqrt(vv) :
                             1.;
  m2denv->zmini = mini - delta * percent;
  m2denv->zmaxi = maxi + delta * percent;

  if (verbose)
  {
    mestitle(2, "Global Statistics on Raw Elevations (extended by %4.2lf)",
             percent);
    message("Statistics are derived from compiling bounds (when defined)\n");
    message("Number of valid bounds = %d\n", (int) nb);
    message("Mean                   = %lf\n", m2denv->zmean);
    message("St. Deviation          = %lf\n", m2denv->zstdv);
    message("Tolerance              = %lf\n", m2denv->zeps);
    message("Minimum                = %lf\n", m2denv->zmini);
    message("Maximum                = %lf\n", m2denv->zmaxi);
    message("Range                  = %lf\n", m2denv->zmaxi - m2denv->zmini);
  }
}

/****************************************************************************/
/*!
 **  Update global statistics on the raw information
 **
 ** \param[in]  m2denv      M2D_Environ structure
 ** \param[in]  dbc         Db structure for constraints
 ** \param[in]  nlayer      Number of layers
 ** \param[in]  verbose     Verbose flag
 **
 *****************************************************************************/
static void st_m2d_stats_updt(M2D_Environ *m2denv,
                              Db *dbc,
                              int nlayer,
                              int verbose)
{
  int nech;
  double nb, mm, vv, mini, maxi, zval, delta;
  static double percent = 0.05;

  /* Initializations */

  nech = dbc->getSampleNumber();
  nb = mm = vv = 0.;
  mini = 1.e30;
  maxi = -1.e30;

  /* Loop on the layers */

  for (int ilayer = 0; ilayer < nlayer; ilayer++)
  {

    /* Loop on the samples */

    for (int iech = 0; iech < nech; iech++)
    {
      zval = dbc->getVariable(iech, ilayer);

      nb += 1.;
      mm += zval;
      vv += zval * zval;
      if (zval < mini) mini = zval;
      if (zval > maxi) maxi = zval;
    }
  }

  /* Normation */

  if (nb > 0)
  {
    mm /= nb;
    vv = vv / nb - mm * mm;
  }
  else
  {
    mm = 0.;
    vv = 1.;
    mini = -0.5;
    maxi = 0.5;
  }

  delta = maxi - mini;
  if (delta <= 0.) delta = ABS(mm) / 10.;
  if (delta <= 0) delta = 1.;
  m2denv->zmean = mm;
  m2denv->zeps = ABS(mm) / 1.e4;
  m2denv->zstdv = (vv > 0) ? sqrt(vv) :
                             1.;
  m2denv->zmini = mini - delta * percent;
  m2denv->zmaxi = maxi + delta * percent;

  if (verbose)
  {
    mestitle(2, "Global Statistics on Centered Elevations");
    message("Statistics are compiled from initial values within bounds\n");
    message("Number of values = %d\n", (int) nb);
    message("Mean             = %lf\n", m2denv->zmean);
    message("St. Deviation    = %lf\n", m2denv->zstdv);
    message("Tolerance        = %lf\n", m2denv->zeps);
    message("Minimum          = %lf\n", m2denv->zmini);
    message("Maximum          = %lf\n", m2denv->zmaxi);
    message("Range            = %lf\n", m2denv->zmaxi - m2denv->zmini);
  }
}

/****************************************************************************/
/*!
 **  Set the initial elevations at the constraining information
 **
 ** \return  Error returned code
 **
 ** \param[in]  m2denv      M2D_Environ structure
 ** \param[in]  dbc         Db structure for constraints
 ** \param[in]  nlayer      Number of layers
 **
 ** \param[out] work        Array of tentative values (Dimension: nlayer)
 **
 ** \remarks This function also add the attributes to 'dbin' per layer:
 ** \remarks - the initial value (ELoc::Z)
 **
 *****************************************************************************/
static int st_m2d_initial_elevations(M2D_Environ *m2denv,
                                     Db *dbc,
                                     int nlayer,
                                     double *work)
{
  int nech, flag_jter;
  double zmin, zmax, zval, eps;
  static int njter_max = 20;

  /* Initializations */

  nech = dbc->getSampleNumber();
  eps = m2denv->zeps;

  /* Loop on the samples */

  for (int iech = 0; iech < nech; iech++)
  {

    /* Define the values at sample as unconstrained information */

    for (int ilayer = 0; ilayer < nlayer; ilayer++)
    {
      zmin = dbc->getLowerBound(iech, ilayer);
      zmax = dbc->getUpperBound(iech, ilayer);
      work[ilayer] = st_m2d_draw_elevation(m2denv, nlayer, ilayer, zmin, zmax);
    }

    /* Loop on iterations for ordering the values */

    for (int jter = 0; jter < njter_max; jter++)
    {
      flag_jter = 0;

      /* Loop on the layers */

      for (int ilayer = 0; ilayer < nlayer; ilayer++)
      {

        /* Determine the bounds at data locations */

        zmin = dbc->getLowerBound(iech, ilayer);
        zmax = dbc->getUpperBound(iech, ilayer);

        /* Loop on the other layers */

        for (int jlayer = 0; jlayer < nlayer; jlayer++)
        {
          if (ilayer == jlayer) continue;
          zval = work[jlayer];

          if (jlayer < ilayer)
          {

            // Comparing with a layer located shallower than the current one

            if (!FFFF(zmax) && zval > zmax)
              flag_jter = 1;
            else
            {
              if (FFFF(zmin))
                zmin = zval + eps;
              else
                zmin = MAX(zmin, zval);
            }
          }
          else
          {

            // Comparing with a layer located deeper than the current one

            if (!FFFF(zmin) && zval < zmin)
              flag_jter = 1;
            else
            {
              if (FFFF(zmax))
                zmax = zval - eps;
              else
                zmax = MIN(zmax, zval);
            }
          }
        }

        // Update target value according to constraints

        work[ilayer] = st_m2d_draw_elevation(m2denv, nlayer, ilayer, zmin,
                                             zmax);
      }

      /* Interrupt iterations */

      if (!flag_jter) break;
    }

    // Run abort in case of lack of convergence

    if (flag_jter)
    {
      messerr("At constraining sample #%d/%d, correct interval ordering",
              iech + 1, nech);
      messerr("has not been reached after %d iterations. Run is aborted",
              njter_max);
      for (int ilayer = 0; ilayer < nlayer; ilayer++)
      {
        zmin = dbc->getLowerBound(iech, ilayer);
        zmax = dbc->getUpperBound(iech, ilayer);
        st_print_constraints_per_point(ilayer, iech, work[ilayer],
        TEST,
                                       TEST, zmin, zmax);
      }
      messerr("\n");
      messerr(">>> You should check the ordering of your bound variables");
      messerr(">>>in the Well File");
      messerr("\n");
      return (1);
    }

    /* Store the resulting values */

    for (int ilayer = 0; ilayer < nlayer; ilayer++)
      dbc->setVariable(iech, ilayer, work[ilayer]);
  }

  return (0);
}

/****************************************************************************/
/*!
 **  Fit the coefficients of the trend terms for each layer
 **
 ** \return  Error returned code
 **
 ** \param[in]  m2denv      M2D_Environ structure
 ** \param[in]  dbin        Db input structure
 ** \param[in]  dbout       Db output structure
 ** \param[in]  nlayer      Number of layers
 ** \param[in]  verbose     Verbose flag
 **
 ** \param[out] iatt_f      Pointer in dbin to the added variables ELoc::F
 **
 ** \remarks This function also add the attributes to 'dbin' per layer:
 ** \remarks - the external drift values (ELoc::F)
 **
 *****************************************************************************/
static int st_m2d_drift_manage(M2D_Environ *m2denv,
                               Db *dbin,
                               Db *dbout,
                               int nlayer,
                               int verbose,
                               int *iatt_f)
{
  int nechin, error, nb;
  double *dval, value, delta;
  static double percent = 0.05;
  VectorInt cols(1);

  /* Initializations */

  error = 1;
  nechin = dbin->getSampleNumber();
  dval = nullptr;
  (*iatt_f) = -1;

  /* Core allocation */

  if (m2denv->flag_ed)
  {
    dval = (double*) mem_alloc(sizeof(double) * nechin, 0);
    if (dval == nullptr) goto label_end;
  }

  /* Add attributes to 'dbin' */
  /* - the external drift value at data points (optional) */
  /* - the initial value at data points */

  if (m2denv->flag_ed)
  {
    if (db_locator_attribute_add(dbin, ELoc::F, nlayer, 0, TEST, iatt_f))
      goto label_end;
  }

  /* Loop on the layers */

  nb = 0;
  for (int ilayer = 0; ilayer < nlayer; ilayer++)
  {

    /* Export External drift from 'dbout' to 'dbin' (optional) */

    if (m2denv->flag_ed)
    {
      cols[0] = dbout->getColumnByLocator(ELoc::F, ilayer);

      // Migrate the information from Grid to Wells

      migrateByAttribute(dbout, dbin, cols, 0, VectorDouble(), false, false);

      // Calculate the statistics of the external drift on the grid

      for (int iech = 0; iech < dbout->getSampleNumber(); iech++)
      {
        if (!dbout->isActive(iech)) continue;
        value = dbout->getExternalDrift(iech, ilayer);
        if (FFFF(value)) continue;
        nb++;
        if (FFFF(m2denv->dmini) || value < m2denv->dmini) m2denv->dmini = value;
        if (FFFF(m2denv->dmaxi) || value > m2denv->dmaxi) m2denv->dmaxi = value;
      }
    }

    /* Loop on the samples */

    for (int iech = 0; iech < nechin; iech++)
    {
      if (!dbin->isActive(iech)) continue;
      if (m2denv->flag_ed)
      {
        if (FFFF(dval[iech])) continue;
        dbin->setExternalDrift(iech, ilayer, dval[iech]);
      }
    }
  }

  /* Patch the statistics on drift if no external drift */

  if (!m2denv->flag_ed)
  {
    m2denv->dmini = 0.;
    m2denv->dmaxi = 1.;
  }
  else
  {
    delta = m2denv->dmaxi - m2denv->dmini;
    m2denv->dmini -= delta * percent;
    m2denv->dmaxi += delta * percent;
  }

  if (verbose)
  {
    mestitle(2, "Global Statistics on Trends (extended by %4.2lf)", percent);
    message("Statistics are derived from compiling drift at grid nodes\n");
    message("Number of valid nodes  = %d\n", (int) nb);
    message("Minimum Drift          = %lf\n", m2denv->dmini);
    message("Maximum Drift          = %lf\n", m2denv->dmaxi);
    message("Range of Drift         = %lf\n", m2denv->dmaxi - m2denv->dmini);
  }

  /* Set the error return code */

  error = 0;

  label_end: dval = (double*) mem_free((char* ) dval);
  return (error);
}

/****************************************************************************/
/*!
 **  Print the details of the constraints
 **
 ** \param[in]  dbc         Db structure for constraints
 ** \param[in]  nech        Number of hard data
 ** \param[in]  ilayer      Rank of the target layer
 **
 *****************************************************************************/
static void st_print_details(Db *dbc, int nech, int ilayer)
{
  double value, lower, upper;
  int nvar, nbdmin, nbdmax;

  nvar = nbdmin = nbdmax = 0;
  for (int iech = 0; iech < nech; iech++)
  {
    value = dbc->getVariable(iech, ilayer);
    if (!FFFF(value)) nvar++;
    lower = dbc->getLowerBound(iech, ilayer);
    upper = dbc->getUpperBound(iech, ilayer);
    if (!FFFF(lower)) nbdmin++;
    if (!FFFF(upper)) nbdmax++;
  }

  // Printout

  message("  . Number of hard data    = %d\n", nvar);
  message("  . Number of lower limits = %d\n", nbdmin);
  message("  . Number of upper limits = %d\n", nbdmax);
}

/****************************************************************************/
/*!
 **  Fit the coefficients of the trend terms for each layer
 **
 ** \return  Error returned code
 **
 ** \param[in]  m2denv      M2D_Environ structure
 ** \param[in]  dbc         Db structure for constraints
 ** \param[in]  nlayer      Number of layers
 ** \param[in]  number_hard Number of hard data used to fit the drift
 ** \param[in]  verbose     Verbose flag
 **
 ** \remarks The drift is only established on data where lower and upper bounds
 ** \remarks are both defined. The drift coefficients are assumed to be the same
 ** \remarks for all layers
 ** \remarks The impact of areal constraints being to important, it has been
 ** \remarks chosen to base the drift fitting only on the first 'number_hard'
 ** \remarks samples (which correspond to constraints coming from 'dbin'.
 **
 *****************************************************************************/
static int st_m2d_drift_fitting(M2D_Environ *m2denv,
                                Db *dbc,
                                int nlayer,
                                int number_hard,
                                int verbose)
{
  int nech, error, numb, nbfl;
  double ff, *a, *b, mean, ffmean, stdv, epais, mini, maxi, ffmini, ffmaxi;

  /* Initializations */

  error = 1;
  nech = MIN(number_hard, dbc->getSampleNumber());
  nbfl = 1;
  a = b = nullptr;

  /* Core allocation */

  m2denv->dcoef = (double*) mem_alloc(sizeof(double) * nlayer, 0);
  if (m2denv->dcoef == nullptr) goto label_end;
  a = (double*) mem_alloc(sizeof(double) * nbfl * nbfl, 0);
  if (a == nullptr) goto label_end;
  b = (double*) mem_alloc(sizeof(double) * nbfl, 0);
  if (b == nullptr) goto label_end;

  /* Loop on the layers */

  for (int ilayer = 0; ilayer < nlayer; ilayer++)
  {

    /* Initializations */

    numb = 0;
    mean = ffmean = stdv = 0.;
    mini = ffmini = 1.e30;
    maxi = ffmaxi = -1.e30;
    for (int i = 0; i < nbfl; i++)
      b[i] = 0.;
    for (int i = 0; i < nbfl * nbfl; i++)
      a[i] = 0.;

    /* Loop on the samples */

    for (int iech = 0; iech < nech; iech++)
    {

      /* Get the values at the data point */

      epais = dbc->getVariable(iech, ilayer);
      if (ilayer > 0)
        epais -= dbc->getVariable(iech, ilayer - 1);
      else
        epais -= m2denv->zmini;

      /* Set the drift vector at data point */

      if (m2denv->flag_ed)
        ff = st_m2d_external_drift_increment(m2denv, dbc, ilayer, iech);
      else
        ff = 1.;

      /* Update statistics */

      numb += 1;
      mean += epais;
      stdv += epais * epais;
      ffmean += ff;
      if (epais < mini) mini = epais;
      if (epais > maxi) maxi = epais;
      if (ff < ffmini) ffmini = ff;
      if (ff > ffmaxi) ffmaxi = ff;

      /* Fill the linear system */

      b[0] = ff * epais;
      a[0] = ff * ff;
    }

    /* Save the results */

    DCOEF(ilayer) = b[0] / a[0];

    /* Normalize statistics */

    if (numb > 0)
    {
      mean /= numb;
      ffmean /= numb;
      stdv = stdv / numb - mean * mean;
      stdv = (stdv > 0) ? sqrt(stdv) :
                          0.;
    }

    /* Print statistics (optional) */

    if (verbose)
    {
      message("\nLayer #%d\n", ilayer + 1);
      message("- Number of Constraints = %d \n", numb);
      st_print_details(dbc, nech, ilayer);
      message("- Drift:\n");
      if (m2denv->flag_ed)
      {
        message("  . Mean          = %lf\n", ffmean);
        message("  . Minimum       = %lf\n", ffmini);
        message("  . Maximum       = %lf\n", ffmaxi);
      }
      message("  . Coefficient   = %lg\n", DCOEF(ilayer));
      message("- Residual:\n");
      message("  . Mean          = %lf\n", mean);
      message("  . St. Deviation = %lf\n", stdv);
      message("  . Minimum       = %lf\n", mini);
      message("  . Maximum       = %lf\n", maxi);
    }
  }

  /* Set the error return cde */

  error = 0;

  label_end: a = (double*) mem_free((char* ) a);
  b = (double*) mem_free((char* ) b);
  return (error);
}

/****************************************************************************/
/*!
 **  Save the drift at the grid nodes
 **
 ** \param[in]  m2denv      M2D_Environ structure
 ** \param[in]  dbout       Db otput structure
 ** \param[in]  nlayer      Number of layers
 **
 ** \param[out] gwork       Working array
 **
 ** \remarks The drift returned as a surface uses directly the coefficients
 ** \remarks of the linear combinaison.
 **
 *****************************************************************************/
static void st_m2d_drift_save(M2D_Environ *m2denv,
                              Db *dbout,
                              int nlayer,
                              double *gwork)
{
  double drift, value;
  int ngrid;

  /* Initializations */

  ngrid = dbout->getSampleNumber();

  /* Loop on the target nodes */

  for (int igrid = 0; igrid < ngrid; igrid++)
  {
    if (!dbout->isActive(igrid)) continue;
    drift = m2denv->zmini;

    /* Loop no the layers */

    for (int ilayer = 0; ilayer < nlayer; ilayer++)
    {
      value = st_m2d_get_drift(m2denv, dbout, ilayer, igrid);
      if (FFFF(value))
        drift = TEST;
      else
        drift += value;
      GWORK(ilayer,igrid) = drift;
    }
  }
}

/****************************************************************************/
/*!
 **  Check if a sample must be considered as an active constraint
 **
 ** \return  1 if the sample is active; 0 otherwise
 **
 ** \param[in]  db          Db input structure
 ** \param[in]  ndim        Space dimension
 ** \param[in]  nlayer      Number of layers
 ** \param[in]  iech        Rank of the sample
 ** \param[in]  bypass      1 to bypass check that at least one bound is defined
 **
 ** \remark A sample is an active constraint if at least one constraint
 ** \remark is defined
 **
 *****************************************************************************/
static int st_active_sample(Db *db, int ndim, int nlayer, int iech, int bypass)
{
  double vmin, vmax;

  /* Check on the coordinates */

  for (int idim = 0; idim < ndim; idim++)
    if (FFFF(db->getCoordinate(iech, idim))) return (0);

  /* Check on the inequality bounds */

  for (int ilayer = 0; ilayer < nlayer; ilayer++)
  {
    vmin = db->getLowerBound(iech, ilayer);
    vmax = db->getUpperBound(iech, ilayer);
    if (!bypass)
    {
      if (FFFF(vmin) && FFFF(vmax)) continue;
      if (!FFFF(vmin) && !FFFF(vmax) && vmin > vmax) continue;
    }
    return (1);
  }
  return (0);
}

/****************************************************************************/
/*!
 **  Record a new active point
 **
 ** \return  Error returned code
 **
 ** \param[in]  m2denv      M2D_Environ structure
 ** \param[in]  db          Db structure
 ** \param[in]  iech        Sample rank in 'db'
 ** \param[in]  ndim        Space dimension
 ** \param[in]  natt        Number of attributes
 ** \param[in]  nlayer      Number of layers
 ** \param[in]  bypass      1 to bypass check that at least one bound is defined
 **
 ** \param[in,out] number_arg Number of samples
 ** \param[in,out] tab        Array of samples
 **
 *****************************************************************************/
static int st_record_sample(M2D_Environ *m2denv,
                            Db *db,
                            int iech,
                            int ndim,
                            int natt,
                            int nlayer,
                            int bypass,
                            int *number_arg,
                            double *tab)
{
  double lower, upper;
  int ecr, number;

  // Skip the record 

  number = *number_arg;
  if (!db->isActive(iech)) return (0);
  if (!st_active_sample(db, ndim, nlayer, iech, bypass)) return (0);

  // Perform the different assignments

  ecr = number * natt;

  // Set the rank

  tab[ecr++] = (double) number + 1;

  // Set the coordinates

  for (int idim = 0; idim < ndim; idim++)
    tab[ecr++] = db->getCoordinate(iech, idim);

  // For each layer, set the bounds and the initial value

  for (int ilayer = 0; ilayer < nlayer; ilayer++)
  {
    lower = db->getLowerBound(iech, ilayer);
    upper = db->getUpperBound(iech, ilayer);

    tab[ecr++] = lower;
    tab[ecr++] = upper;
    tab[ecr++] = TEST;
  }

  // For each layer, set the External Drift value (optional) 

  if (m2denv->flag_ed) for (int ilayer = 0; ilayer < nlayer; ilayer++)
    tab[ecr++] = db->getExternalDrift(iech, ilayer);

  /* Increment the number of records by 1 */

  number++;

  /* Set the returned arguments */

  *number_arg = number;
  return (0);
}

/****************************************************************************/
/*!
 **  Define the locators on the newly created Db
 **
 ** \param[in]  m2denv      M2D_Environ structure
 ** \param[in]  db          Db constraints structure
 ** \param[in]  ndim        Number of coodinates
 ** \param[in]  nlayer      Number of layers
 ** \param[in]  nvar        Number of variables
 **
 *****************************************************************************/
static void st_define_locators(M2D_Environ *m2denv,
                               Db *db,
                               int ndim,
                               int nvar,
                               int nlayer)
{
  int ivar;

  ivar = 1;
  db->setLocatorsByAttribute(ndim, ivar, ELoc::X);
  ivar += ndim;
  for (int ilayer = 0; ilayer < nlayer; ilayer++)
  {
    db->setLocatorByAttribute(ivar++, ELoc::L, ilayer);
    db->setLocatorByAttribute(ivar++, ELoc::U, ilayer);
    if (ilayer < nvar) db->setLocatorByAttribute(ivar, ELoc::Z, ilayer);
    ivar++;
  }
  if (m2denv->flag_ed) db->setLocatorsByAttribute(nlayer, ivar, ELoc::F);
}

/****************************************************************************/
/*!
 **  Print the Environnement
 **
 *****************************************************************************/
static void st_m2d_print_environ(const char *title, M2D_Environ *m2denv)
{
  mestitle(1, title);

  if (m2denv->flag_ed)
    message("Use of External Drift\n");
  else
    message("No External Drift\n");
  message("Z Minimum               = %lf\n", m2denv->zmini);
  message("Z Maximum               = %lf\n", m2denv->zmaxi);
  message("Z Mean                  = %lf\n", m2denv->zmean);
  message("Z St. Deviation         = %lf\n", m2denv->zstdv);
  message("Z Tolerance             = %lf\n", m2denv->zeps);
  message("Drift Minimum           = %lf\n", m2denv->dmini);
  message("Drift Maximum           = %lf\n", m2denv->dmaxi);
  message("Y St. Deviation         = %lf\n", m2denv->ystdv);
}

/****************************************************************************/
/*!
 **  Create a Db containing all the constraining information
 **
 ** \return  Pointer to the newly created Db or NULL
 **
 ** \param[in]  m2denv      M2D_Environ structure
 ** \param[in]  dbin        Db input structure
 ** \param[in]  dbout       Db output structure
 ** \param[in]  ndim        Space dimension
 ** \param[in]  nlayer      Number of layers
 **
 ** \param[out] number_hard Number of hard data which will serve for
 **                         seting the optimal drift
 **
 ** \remark Note that the file constructed here contains as many samples
 ** \remark as the number of ACTIVE samples of the input Db
 **
 *****************************************************************************/
static Db* st_m2d_create_constraints(M2D_Environ *m2denv,
                                     Db *dbin,
                                     Db *dbout,
                                     int ndim,
                                     int nlayer,
                                     int *number_hard)
{
  Db *db;
  VectorDouble tab;
  int nechin, nechout, nech, natt, number, error, ecr;

  /* Initializations */

  error = 1;
  db = nullptr;
  nechin = nechout = 0;

  nechin = dbin->getActiveSampleNumber();
  nechout = dbout->getActiveSampleNumber();
  nech = nechin + nechout;
  natt = 1;                  // Rank
  natt += ndim;               // Coordinates
  natt += 3 * nlayer;         // LowBound, UppBound and Variable per layer
  if (m2denv->flag_ed) natt += nlayer;  // External Drift

  /* Core allocation */

  tab.resize(nech * natt);

  /* Load information from 'dbin' */

  number = 0;
  for (int iech = 0; iech < dbin->getSampleNumber(); iech++)
  {
    if (!dbin->isActive(iech)) continue;
    if (st_record_sample(m2denv, dbin, iech, ndim, natt, nlayer, 0, &number,
                         tab.data())) goto label_end;
  }
  *number_hard = number;

  /* Load information from 'dbout' */

  for (int iech = 0; iech < dbout->getSampleNumber(); iech++)
  {
    if (!dbout->isActive(iech)) continue;
    if (st_record_sample(m2denv, dbout, iech, ndim, natt, nlayer, 0, &number,
                         tab.data())) goto label_end;
  }

  /* When forcing, the first active sample is used */

  if (number <= 0)
  {
    for (int iech = 0; iech < dbin->getSampleNumber(); iech++)
    {
      if (!dbin->isActive(iech)) continue;
      if (st_record_sample(m2denv, dbin, iech, ndim, natt, nlayer, 1, &number,
                           tab.data())) goto label_end;
      if (number > 0) break;
    }
  }

  /* Core reallocation */

  if (number < nech) tab.resize(number * natt);

  /* Create the output Db */

  db = db_create_point(number, natt, ELoadBy::SAMPLE, 0, tab);
  if (db == nullptr) goto label_end;

  // Assigning names to the variables (not pointers yet)

  ecr = 0;
  db_name_set(db, ecr++, "rank");
  for (int idim = 0; idim < ndim; idim++)
  {
    (void) gslSPrintf(string_encode, "X%d", idim + 1);
    db_name_set(db, ecr++, string_encode);
  }
  for (int ilayer = 0; ilayer < nlayer; ilayer++)
  {
    (void) gslSPrintf(string_encode, "Lower%d", ilayer + 1);
    db_name_set(db, ecr++, string_encode);
    (void) gslSPrintf(string_encode, "Upper%d", ilayer + 1);
    db_name_set(db, ecr++, string_encode);
    (void) gslSPrintf(string_encode, "Value%d", ilayer + 1);
    db_name_set(db, ecr++, string_encode);
  }
  if (m2denv->flag_ed)
  {
    for (int ilayer = 0; ilayer < nlayer; ilayer++)
    {
      (void) gslSPrintf(string_encode, "Drift%d", ilayer + 1);
      db_name_set(db, ecr++, string_encode);
    }
  }

  /* Set the error return code */

  error = 0;

  label_end: if (error) db = db_delete(db);
  return (db);
}

/****************************************************************************/
/*!
 **  Calculate the inverse of (s2 * Q + B %*% Bt) and store it
 **  into a new QChol object
 **
 ** \return  The calculated sparse matrix
 **
 ** \param[in]  s2          Nugget effect value
 ** \param[in]  Qc          Qc structure (already existing)
 ** \param[in]  Matelem     Matelem structure
 **
 *****************************************************************************/
static QChol* st_derive_Qc(double s2, QChol *Qc, SPDE_Matelem &Matelem)
{
  cs *Bt, *B2, *Q, *B;
  int error;

  // Initializations 

  error = 1;
  Bt = B2 = nullptr;
  Q = Matelem.QC->Q;
  B = Matelem.Aproj;

  // Clean the previous Qc (if it exists)

  if (Qc != nullptr) Qc = qchol_manage(-1, Qc);

  // Calculate: Q + t(B) %*% B

  message("Building Q (Size:%d) with additional nugget effect (%lf) ... ", Q->n,
          s2);
  Bt = cs_transpose(B, 1);
  if (Bt == nullptr) goto label_end;
  B2 = cs_multiply(Bt, B);
  if (B2 == nullptr) goto label_end;

  Qc = qchol_manage(1, NULL);
  if (Qc == nullptr) goto label_end;
  Qc->Q = cs_add(Q, B2, s2, 1.);
  if (Qc->Q == nullptr) goto label_end;

  // Perform the Cholesky transform 

  error = qchol_cholesky(0, Qc);

  // Free memory

  label_end: message("Done\n");
  Bt = cs_spfree(Bt);
  B2 = cs_spfree(B2);
  if (error) Qc = qchol_manage(-1, Qc);
  return (Qc);
}

/****************************************************************************/
/*!
 **  Returns the projection matrix of a set of points (contained in a Db)
 **  onto a meshing
 **
 ** \return Pointer to the newly created sparse matrix (or NULL)
 **
 ** \param[in]  db         Db structure
 ** \param[in]  s_mesh     SPDE_Mesh structure
 ** \param[in]  flag_exact Type of test for intersection (See remarks)
 ** \param[in]  radius     Neighborhood radius
 ** \param[in]  verbose    Verbose flag
 **
 ** \param[out] nactive_arg Number of active samples from the Db
 ** \param[out] ranks_arg   Ranks of the active samples
 **
 ** \remarks The calling function must free the argument 'ranks'
 **
 ** \remarks When flag_exact is TRUE, for each active sample of Db, a vertex
 ** \remarks of the mesh is active as soon as it lies within the vicinity
 ** \remarks of the sample.
 ** \remarks If flag_exact is FALSE, all vertices of a mesh are considered as
 ** \remarks active as soon as the mesh intersects the ball around a sample.
 **
 ** \remarks The vicinity is defined as any point located at a distance
 ** \remarks from the sample smaller than 'radius'. The distance is calculated
 ** \remarks as the Euclidean distance over the space whose dimension is
 ** \remarks if the smallest value between the Db et Mesh space dimensions.
 **
 *****************************************************************************/
cs* db_mesh_neigh(const Db *db,
                  SPDE_Mesh *s_mesh,
                  double radius,
                  int flag_exact,
                  int verbose,
                  int *nactive_arg,
                  int **ranks_arg)
{
  double *coor, *caux, total;
  int *pts, *ranks, error, ncorner, ip, flag_sphere, ndimd, ndimv, ndim, jech;
  int jech_max, ip_max, nech, nactive;
  cs *A, *Atriplet;

  /* Initializations */

  error = 1;
  coor = nullptr;
  caux = nullptr;
  pts = nullptr;
  Atriplet = A = nullptr;
  ncorner = s_mesh->ncorner;
  nech = db->getSampleNumber();
  variety_query(&flag_sphere);
  if (flag_sphere)
  {
    messerr("The function 'db_mesh_neigh' is not programmed on sphere");
    goto label_end;
  }

  /* Create the Triplet container */

  Atriplet = cs_spalloc(0, 0, 1, 1, 1);
  if (Atriplet == nullptr) goto label_end;

  /* Core allocation */

  ranks = (int*) mem_alloc(sizeof(int) * nech, 0);
  if (ranks == nullptr) goto label_end;
  for (int iech = 0; iech < nech; iech++)
    ranks[iech] = -1;

  /* Core allocation */

  ndimd = db->getNDim();
  ndimv = s_mesh->ndim;
  ndim = MIN(ndimd, ndimv);
  coor = db_sample_alloc(db, ELoc::X);
  if (coor == nullptr) goto label_end;
  caux = db_sample_alloc(db, ELoc::X);
  if (caux == nullptr) goto label_end;
  pts = (int*) mem_alloc(sizeof(int) * s_mesh->nvertex, 0);
  if (pts == nullptr) goto label_end;

  /* Loop on the samples */

  ip_max = jech_max = jech = 0;
  for (int iech = 0; iech < nech; iech++)
  {
    if (!db->isActive(iech)) continue;

    /* Identification of the sample in the meshing */

    for (int idim = 0; idim < ndim; idim++)
      coor[idim] = db->getCoordinate(iech, idim);

    /* Blank out the array of hitting points */

    for (ip = 0; ip < s_mesh->nvertex; ip++)
      pts[ip] = 0;

    /* Loop on the meshes */

    for (int imesh = 0; imesh < s_mesh->nmesh; imesh++)
    {
      if (flag_exact)
      {
        for (int icorn = 0; icorn < s_mesh->ncorner; icorn++)
        {
          ip = OLD_MESHES(imesh, icorn);
          for (int idim = 0; idim < ndim; idim++)
            caux[idim] = POINTS(ip, idim);
          if (ut_distance(ndim, coor, caux) <= radius)
          {
            pts[ip] = 1;
            break;
          }
        }
      }
      else
      {
        if (!is_in_mesh_neigh(s_mesh, coor, caux, ndim, imesh, radius))
          continue;

        /* The meshing element is in the neighborhood of the sample */

        for (int icorn = 0; icorn < ncorner; icorn++)
        {
          ip = OLD_MESHES(imesh, icorn);
          pts[ip] = 1;
          if (ip > ip_max) ip_max = ip;
        }
      }
    }

    /* Count the number of vertices hit by the point neighborhood */

    total = 0.;
    for (ip = 0; ip < s_mesh->nvertex; ip++)
      total += pts[ip];
    if (total <= 0.) continue;

    /* Add the active vertices to the triplet */

    for (ip = 0; ip < s_mesh->nvertex; ip++)
    {
      if (pts[ip] <= 0) continue;
      if (ip > ip_max) ip_max = ip;
      if (!cs_entry(Atriplet, jech, ip, 1. / total)) goto label_end;
    }
    ranks[jech] = iech;
    if (jech > jech_max) jech_max = jech;
    jech++;
  }

  /* Add the extreme value to force dimension */

  if (ip_max < s_mesh->nvertex - 1)
  {
    if (!cs_entry(Atriplet, jech_max, s_mesh->nvertex - 1, 0.)) goto label_end;
  }

  /* Core reallocation */

  nactive = jech_max + 1;
  ranks = (int*) mem_realloc((char* ) ranks, sizeof(int) * nactive, 0);
  if (ranks == nullptr) goto label_end;

  /* Convert the triplet into a sparse matrix */

  A = cs_triplet(Atriplet);

  /* Set the error return code */

  error = 0;
  *nactive_arg = nactive;
  *ranks_arg = ranks;

  label_end: Atriplet = cs_spfree(Atriplet);
  pts = (int*) mem_free((char* ) pts);
  coor = db_sample_free(coor);
  caux = db_sample_free(caux);
  if (error) A = cs_spfree(A);
  return (A);
}

/****************************************************************************/
/*!
 **  Draw a Z-value within bounds
 **
 ** \return Z-Value in the working domain
 **
 ** \param[in]  m2denv      M2D_Environ structure
 ** \param[in]  dbc         Db structure
 ** \param[in]  verbose     Verbose flag
 ** \param[in]  iter        Rank of the iteration
 ** \param[in]  ilayer      Rank of the layer
 ** \param[in]  iech        Rank of the sample
 ** \param[in]  Zval        Input value
 ** \param[in]  Zcum        Cumulated Z value of layers above
 ** \param[in]  Zmin        Lower bound in Z
 ** \param[in]  Zmax        Upper bound in Z
 ** \param[in]  Ymean       Mean of the Y Law
 ** \param[in]  Ysigma      Standard deviation of the Y Law
 **
 *****************************************************************************/
static double st_m2d_draw_gaussian(M2D_Environ *m2denv,
                                   Db *dbc,
                                   int verbose,
                                   int iter,
                                   int ilayer,
                                   int iech,
                                   double Zval,
                                   double Zcum,
                                   double Zmin,
                                   double Zmax,
                                   double Ymean,
                                   double Ysigma)
{
  double M, S, Yval, Ymin, Ymax, Zminc, Zmaxc;
  static int verif = 1;

  /* Initializations */

  M = st_m2d_get_M(m2denv, dbc, 1, ilayer, iech);
  S = st_m2d_get_S(m2denv, dbc, 1, ilayer, iech);
  if (st_check_validity_MS(dbc, ilayer, iech, 1, 1, M, S))
    messageAbort("- Impossible to have M or S undefined");

  if (verbose)
  {
    message("Input Z elevation=%lf", Zval);
    st_print_concatenate_interval(NULL, Zmin, Zmax, 1);
  }

  /* Centering in Z */

  Zminc = Zmin;
  if (!FFFF(Zminc)) Zminc -= Zcum;
  if (Zminc < 0) Zminc = 0.;
  Zmaxc = Zmax;
  if (!FFFF(Zmaxc)) Zmaxc -= Zcum;
  if (Zmaxc < 0) Zmaxc = 0.;
  if (verbose) st_print_concatenate_interval("Z thickness", Zminc, Zmaxc, 1);

  /* Converting from Z to Y */

  if (FFFF(Zminc))
    Ymin = TEST;
  else
    Ymin = (Zminc == 0) ? TEST :
                          (S * S / 2. + log(Zminc / M)) / S;
  if (FFFF(Zmaxc))
    Ymax = TEST;
  else
    Ymax = (Zmaxc == 0) ? TEST :
                          (S * S / 2. + log(Zmaxc / M)) / S;
  if (verbose) st_print_concatenate_interval("Y gaussian", Ymin, Ymax, 1);

  /* Centering in Y */

  if (!FFFF(Ymin)) Ymin = (Ymin - Ymean) / Ysigma;
  if (!FFFF(Ymax)) Ymax = (Ymax - Ymean) / Ysigma;
  if (verbose) st_print_concatenate_interval("Y centered", Ymin, Ymax, 1);

  Yval = law_gaussian_between_bounds(Ymin, Ymax);
  // Two next lines are there for robustification: they should not be removed
  if (!FFFF(Ymin) && Yval < Ymin) Yval = Ymin;
  if (!FFFF(Ymax) && Yval > Ymax) Yval = Ymax;

  Yval = Ymean + Ysigma * Yval;
  Zval = Zcum + M * exp(S * Yval - S * S / 2.);

  if (verif)
  {
    if (std::isinf(Zval))
    {
      message("Iteration #%d - Layer #%d - Sample #%d\n", iter + 1, ilayer + 1,
              iech + 1);
      message("  Zval=Inf");
      st_print_concatenate_interval(NULL, Zmin, Zmax, 1);
      messageAbort("Strange output value for Zval");
    }
    if (!FFFF(Zmin) && Zval < Zmin - m2denv->zeps)
    {
      message("Iteration #%d - Layer #%d - Sample #%d\n", iter + 1, ilayer + 1,
              iech + 1);
      message(" Zval=%lf", Zval);
      st_print_concatenate_interval(NULL, Zmin, Zmax, 1);
      message(" Yval=%lf", Yval);
      st_print_concatenate_interval(NULL, Ymin, Ymax, 1);
      messageAbort("Zval should not be smaller than Zmin");
    }
    if (!FFFF(Zmax) && Zval > Zmax + m2denv->zeps)
    {
      message("Iteration #%d - Layer #%d - Sample #%d\n", iter + 1, ilayer + 1,
              iech + 1);
      message(" Zval=%lf", Zval);
      st_print_concatenate_interval(NULL, Zmin, Zmax, 1);
      message(" Yval=%lf", Yval);
      st_print_concatenate_interval(NULL, Ymin, Ymax, 1);
      messageAbort("Zval should not be larger than Zmax");
    }
  }

  if (verbose)
  {
    message("Output Z elevation=%lf in", Zval);
    st_print_concatenate_interval(NULL, Zmin, Zmax, 1);
  }

  return (Zval);
}

/****************************************************************************/
/*!
 **  Convert a layer-pile at a datum from the working to the true domain
 **
 ** \param[in]  m2denv      M2D_Environ structure
 ** \param[in]  dbc         Db structure
 ** \param[in]  nlayer      Number of layers
 ** \param[in]  type        1 for the constraining Db
 **                         2 for the grid output Db
 ** \param[in]  iech        Rank of the sample
 ** \param[in,out] tab      Input/Output array of Z-values (Dimension: nlayer)
 **
 *****************************************************************************/
static void st_convert_Z2Y(M2D_Environ *m2denv,
                           Db *dbc,
                           int nlayer,
                           int type,
                           int iech,
                           double *tab)
{
  double M, S, Yval, Zval, Zcur, Zcum;
  int flag_undef;

  flag_undef = 0;
  Zcum = m2denv->zmini;
  for (int ilayer = 0; ilayer < nlayer; ilayer++)
  {
    M = st_m2d_get_M(m2denv, dbc, type, ilayer, iech);
    S = st_m2d_get_S(m2denv, dbc, type, ilayer, iech);
    if (st_check_validity_MS(dbc, ilayer, iech, 1, 1, M, S) || flag_undef)
    {
      flag_undef = 1;
      Yval = TEST;
    }
    else
    {
      Zcur = tab[ilayer];
      Zval = Zcur - Zcum;
      if (Zval <= 0)
      {
        flag_undef = 1;
        Yval = TEST;
      }
      else
        Yval = (S * S / 2. + log(Zval / M)) / S;
      Zcum = Zcur;
    }
    tab[ilayer] = Yval;
  }
  return;
}

/****************************************************************************/
/*!
 **  Convert a layer-pile at a datum from the true to the working domain
 **
 ** \param[in]  m2denv      M2D_Environ structure
 ** \param[in]  db          Db structure
 ** \param[in]  nlayer      Number of layers
 ** \param[in]  type        1 for the constraining Db
 **                         2 for the grid output Db
 ** \param[in]  iech        Rank of the sample
 ** \param[in,out] tab      Input/Ouput array of Y-values (Dimension: nlayer)
 **
 *****************************************************************************/
static void st_convert_Y2Z(M2D_Environ *m2denv,
                           Db *db,
                           int nlayer,
                           int type,
                           int iech,
                           double *tab)
{
  double M, S, Zval, Yval, Zcur;
  int flag_undef;

  flag_undef = 0;
  Zcur = m2denv->zmini;
  for (int ilayer = 0; ilayer < nlayer; ilayer++)
  {
    M = st_m2d_get_M(m2denv, db, type, ilayer, iech);
    S = st_m2d_get_S(m2denv, db, type, ilayer, iech);
    if (st_check_validity_MS(db, ilayer, iech, 0, 0, M, S) || flag_undef)
    {
      flag_undef = 1;
      Zcur = TEST;
    }
    else
    {
      Yval = tab[ilayer];
      Zval = M * exp(S * Yval - S * S / 2.);
      Zcur += Zval;
    }
    tab[ilayer] = Zcur;
  }
  return;
}

/****************************************************************************/
/*!
 **  Print the values at a sample location
 **
 ** \param[in]  m2denv      M2D_Environ structure
 ** \param[in]  title       Title
 ** \param[in]  dbc         Db structure containing the constraints
 ** \param[in]  nlayer      Number of layers
 ** \param[in]  iech        Sample rank
 ** \param[in]  work        Array of values (defined in Z)
 **
 *****************************************************************************/
static void st_print_sample(const char *title,
                            M2D_Environ *m2denv,
                            Db *dbc,
                            int nlayer,
                            int iech,
                            double *work)
{
  int nech;
  double zmin, zmax;

  nech = dbc->getSampleNumber();
  message("%s - Sample #%d/%d\n", title, iech + 1, nech);

  for (int ilayer = 0; ilayer < nlayer; ilayer++)
  {
    zmin = dbc->getLowerBound(iech, ilayer);
    zmax = dbc->getUpperBound(iech, ilayer);
    message("Z(%d)=%lf in", ilayer, work[ilayer]);
    st_print_concatenate_interval(NULL, zmin, zmax, 1);
  }
}

/****************************************************************************/
/*!
 **  Perform the Gibbs iterations
 **
 ** \return Error returned code
 **
 ** \param[in]  m2denv      M2D_Environ structure
 ** \param[in]  dbc         Db structure containing the constraints
 ** \param[in]  verbose     Verbose flag
 ** \param[in]  iter        Rank of the iteration
 ** \param[in]  nlayer      Number of layers
 ** \param[in]  sigma       Standard deviation of the nugget value
 ** \param[in]  ymean       Array of mean values at constraints
 ** \param[in,out] ydat     Array of values at constraints samples
 **
 ** \param[out] work        Array of tentative values (Dimension: nlayer)
 **
 ** \remarks The [in,out] argument 'ydat' is expressed in working domain
 ** \remarks It needs to be locally transformed in real domain in order to
 ** \remarks be compared to the bounds.
 **
 *****************************************************************************/
static int st_global_gibbs(M2D_Environ *m2denv,
                           Db *dbc,
                           int verbose,
                           int iter,
                           int nlayer,
                           double sigma,
                           double *ymean,
                           double *ydat,
                           double *work)
{
  int nech;
  double zval, zmin, zmax, zcum;

  // Initializations

  nech = dbc->getSampleNumber();

  // Loop on the samples

  for (int iech = 0; iech < nech; iech++)
  {

    // Set the initial values

    for (int ilayer = 0; ilayer < nlayer; ilayer++)
      work[ilayer] = YDAT(ilayer, iech);
    st_convert_Y2Z(m2denv, dbc, nlayer, 1, iech, work);
    if (verbose)
      st_print_sample("Entering in Gibbs", m2denv, dbc, nlayer, iech, work);

    // Loop on the layers

    for (int ilayer = 0; ilayer < nlayer; ilayer++)
    {

      // Get the elevation of the previous layer

      zcum = m2denv->zmini;

      // Getting the elevation and the bounds for the current layer

      zmin = dbc->getLowerBound(iech, ilayer);
      zmax = dbc->getUpperBound(iech, ilayer);
      if (verbose)
      {
        message("ilayer=%d", ilayer);
        st_print_concatenate_interval(NULL, zmin, zmax, 1);
      }

      // Loop on the other layers

      for (int jlayer = 0; jlayer < nlayer; jlayer++)
      {
        if (ilayer == jlayer) continue;
        zval = work[jlayer];
        if (verbose)
          message("Constrained by jlayer=%d zval=%lf\n", jlayer, zval);

        if (jlayer < ilayer)
        {

          // Comparing with a layer located shallower than the current one

          if (FFFF(zmin))
            zmin = zval;
          else
            zmin = MAX(zmin, zval);
          zcum = zval;
        }
        else
        {

          // Comparing with a layer located deeper than the current one

          if (FFFF(zmax))
            zmax = zval;
          else
            zmax = MIN(zmax, zval);
        }
      }

      // Drawing plausible values according to constraints

      work[ilayer] = st_m2d_draw_gaussian(m2denv, dbc, verbose, iter, ilayer,
                                          iech, work[ilayer], zcum, zmin, zmax,
                                          YMEAN(ilayer, iech), sigma);
    }

    // Load the new values 

    if (verbose)
      st_print_sample("Exiting Gibbs", m2denv, dbc, nlayer, iech, work);
    st_convert_Z2Y(m2denv, dbc, nlayer, 1, iech, work);

    /* Store in the extracted vector */

    for (int ilayer = 0; ilayer < nlayer; ilayer++)
      YDAT(ilayer,iech) = work[ilayer];
  }
  return (0);
}

/****************************************************************************/
/*!
 **  Check that the Gibbs constraints are fullfilled at datum locations
 **
 ** \return  Error return code
 **
 ** \param[in]  title       Title for the printout (if error)
 ** \param[in]  m2denv      M2D_Environ structure
 ** \param[in]  dbc         Db structure containing the constraints
 ** \param[in]  nlayer      Number of layers
 ** \param[in]  verbose     Verbosity flag
 ** \param[in]  ydat        Array of simulations on the data
 **
 ** \param[out] work        Array of tentative values (Dimension: nlayer)
 **
 *****************************************************************************/
static int st_check_gibbs_data(const char *title,
                               M2D_Environ *m2denv,
                               Db *dbc,
                               int nlayer,
                               int verbose,
                               double *ydat,
                               double *work)
{
  int error, nech;
  double zmin, zmax, depth, eps;

  // Initializations

  error = 0;
  nech = dbc->getSampleNumber();
  eps = m2denv->zeps;

  // Loop on the constraints samples

  for (int iech = 0; iech < nech; iech++)
  {
    for (int ilayer = 0; ilayer < nlayer; ilayer++)
      work[ilayer] = YDAT(ilayer, iech);
    st_convert_Y2Z(m2denv, dbc, nlayer, 1, iech, work);

    // Loop on the layers 

    for (int ilayer = 0; ilayer < nlayer; ilayer++)
    {
      depth = work[ilayer];

      // Getting the elevation and the bounds for the current layer

      zmin = dbc->getLowerBound(iech, ilayer);
      zmax = dbc->getUpperBound(iech, ilayer);

      // Check consistency

      if (!FFFF(zmin))
      {
        if (depth < zmin - eps)
        {
          messerr("%s: Sample(%d/%d) of Layer(%d/%d): Depth(%lf) < Lower(%lf)",
                  title, iech + 1, nech, ilayer + 1, nlayer, depth, zmin);
          error++;
        }
      }
      if (!FFFF(zmax))
      {
        if (depth > zmax + eps)
        {
          messerr("%s: Sample(%d/%d) of Layer(%d/%d): Depth(%lf) > Upper(%lf)",
                  title, iech + 1, nech, ilayer + 1, nlayer, depth, zmax);
          error++;
        }
      }
    }
  }

  if (verbose)
  {
    if (error == 0)
      message("%s: No inconsistency\n", title);
    else
      message("%s: %d error(s) found\n", title, error);
  }
  return (error);
}

/****************************************************************************/
/*!
 **  Manage the M2D_Environ structure
 **
 ** \return  Pointer to the M2D_Environ structure
 **
 ** \param[in]  mode        1 for allocation; -1 for deallocation
 ** \param[in]  flag_ed     1 if external drift is used; 0 otherwise
 ** \param[in]  ystdv       Stamdard deviation of the Gaussian Transformed
 ** \param[in]  m2denv_old  Pointer to the already existing M2D_Environ
 **                         (only used when mode==-1)
 **
 *****************************************************************************/
static M2D_Environ* m2denv_manage(int mode,
                                  int flag_ed,
                                  double ystdv,
                                  M2D_Environ *m2denv_old)
{
  M2D_Environ *m2denv;

  /* Dispatch */

  if (mode > 0)
  {

    // Allocation

    m2denv = (M2D_Environ*) mem_alloc(sizeof(M2D_Environ), 0);
    if (m2denv == (M2D_Environ*) NULL) return (m2denv);
    m2denv->flag_ed = flag_ed;
    m2denv->iatt_fd = -1;
    m2denv->iatt_fg = -1;
    m2denv->zmean = 0.;
    m2denv->zstdv = 1.;
    m2denv->zeps = 0.;
    m2denv->zmini = TEST;
    m2denv->zmaxi = TEST;
    m2denv->dmini = TEST;
    m2denv->dmaxi = TEST;
    m2denv->ystdv = ystdv;
    m2denv->dcoef = nullptr;
  }
  else
  {
    m2denv = m2denv_old;
    if (m2denv != (M2D_Environ*) NULL)

    {
      m2denv->dcoef = (double*) mem_free((char* ) m2denv->dcoef);
      m2denv = (M2D_Environ*) mem_free((char* ) m2denv);
    }
  }
  return (m2denv);
}

/****************************************************************************/
/*!
 **  Extract a vector containing the constraints
 **
 ** \param[in]  m2denv      M2D_Environ structure
 ** \param[in]  dbc         Db constraints structure
 ** \param[in]  nlayer      Number of layers
 **
 ** \param[out] ydat        Array of values at constraints samples
 **                         (Dimension: nech * nlayer)
 ** \param[out] work        Array of tentative values (Dimension: nlayer)
 **
 *****************************************************************************/
static void st_m2d_vector_extract(M2D_Environ *m2denv,
                                  Db *dbc,
                                  int nlayer,
                                  double *ydat,
                                  double *work)
{
  int nech;

  /* Initializations */

  nech = dbc->getSampleNumber();

  /* Loop on the samples */

  for (int iech = 0; iech < nech; iech++)
  {

    /* Loop on the layers */

    for (int ilayer = 0; ilayer < nlayer; ilayer++)
      work[ilayer] = dbc->getVariable(iech, ilayer);

    /* Convert from the depth to thickness */

    st_convert_Z2Y(m2denv, dbc, nlayer, 1, iech, work);

    /* Store in the extracted vector */

    for (int ilayer = 0; ilayer < nlayer; ilayer++)
      YDAT(ilayer,iech) = work[ilayer];
  }
}

/****************************************************************************/
/*!
 **  Print the set of constraints
 **
 ** \param[in]  title       Title
 ** \param[in]  db          Db constraints structure
 ** \param[in]  ydat        Array of gaussian values at constraints (optional)
 ** \param[in]  nlayer      Number of layers
 ** \param[in]  verbose     Verbose flag
 **
 ** \remarks This function tends to produce verbose outputs.
 ** \remarks This is the reason why it has been conditioned to print only
 ** \remarks the values of the first samples. This is controled by the
 ** \remarks internal parameter 'nprint' which can be ruled by keypair:
 ** \remarks     set.keypair("Print_Data",0)
 ** \remarks The default number of samples i s 0 (no printout)
 **
 *****************************************************************************/
static void st_print_db_constraints(const char *title,
                                    Db *db,
                                    double *ydat,
                                    int nlayer,
                                    int verbose)
{
  double value, lower, drift, upper, vgaus;
  int nech, nprint;

  // Initializations

  nprint = (int) get_keypone("Print_Data", 10.);
  if (!verbose || nprint == 0) return;

  // Printout

  mestitle(1, title);
  nech = db->getSampleNumber();
  if (nprint > 0) nech = MIN(nech, nprint);
  for (int iech = 0; iech < nech; iech++)
  {
    if (!db->isActive(iech)) continue;
    for (int ilayer = 0; ilayer < nlayer; ilayer++)
    {
      lower = db->getLowerBound(iech, ilayer);
      upper = db->getUpperBound(iech, ilayer);
      value = db->getVariable(iech, ilayer);
      drift = db->getExternalDrift(iech, ilayer);
      vgaus = (ydat != nullptr) ? YDAT(ilayer, iech) :
                                  TEST;
      st_print_constraints_per_point(ilayer, iech, value, drift, vgaus, lower,
                                     upper);
    }
  }
}

/****************************************************************************/
/*!
 **  Print the statistics on the current array
 **
 ** \param[in]  title       Title attache to the printou
 ** \param[in]  nlayer      Number of layers
 ** \param[in]  nech        Number of samples per layer in the target array
 ** \param[in]  ydat        Target (generic) array
 **
 ** \remarks   This function can be used for any array
 **
 *****************************************************************************/
static void st_m2d_stats_gaus(const char *title,
                              int nlayer,
                              int nech,
                              double *ydat)
{
  if (!DEBUG) return;
  for (int ilayer = 0; ilayer < nlayer; ilayer++)
  {
    (void) gslSPrintf(string_encode, "%s (Layer #%d)", title, ilayer + 1);
    ut_stats_mima_print(string_encode, nech, &YDAT(ilayer, 0), NULL);
  }
}

/****************************************************************************/
/*!
 **  Perform Gibbs on a multilayer setup
 **
 ** \return  Error return code
 **
 ** \param[in]  dbin        Db input structure
 ** \param[in]  dbout       Db output structure
 ** \param[in]  model       Model structure
 ** \param[in]  flag_ed     1 if External Drit is used
 ** \param[in]  nlayer      Number of layers
 ** \param[in]  niter       Number of iterations
 ** \param[in]  seed        Seed for random number generator
 ** \param[in]  nbsimu      Number of simulaations
 ** \param[in]  icol_pinch  Address of the variable containing the pinchout
 ** \param[in]  flag_drift  1 to return the drift only
 **                         0 the simulations
 ** \param[in]  flag_ce     1 if the conditional expectation
 **                         should be returned instead of simulations
 ** \param[in]  flag_cstd   1 if the conditional standard deviation
 **                         should be returned instead of simulations
 ** \param[in]  verbose     Verbose option
 **
 ** \remarks In 'dbin':
 ** \remarks - the lower and upper bounds must be defined for each datum
 ** \remarks   (set to the locator ELoc::L and ELoc::U
 ** \remarks In 'dbout':
 ** \remarks - the trend (if flag_ed is 1) must be defined and set to
 ** \remarks   the locator ELoc::F
 ** \remarks When defined, the pinchout should be defined as a grid variable
 ** \remarks with values ranging between 0 and 1 (FFFF are admitted).
 ** \remarks It will serve as a multiplier to the Mean thickness maps.
 **
 *****************************************************************************/
int m2d_gibbs_spde(Db *dbin,
                   Db *dbout,
                   Model *model,
                   int flag_ed,
                   int nlayer,
                   int niter,
                   int seed,
                   int nbsimu,
                   int icol_pinch,
                   int flag_drift,
                   int flag_ce,
                   int flag_cstd,
                   int verbose)
{
  int error, iatt_f, iatt_out, nvertex, nech, ngrid, ndim, number_hard, nfois;
  int iptr_ce, iptr_cstd, ecr;
  double *ydat, *ymean, *yvert, *rhs, *zkrig, *gwork, *vwork, *lwork;
  double *ydat_loc, *ymean_loc, *yvert_loc, nugget, ysigma, vartot;
  cs *Bproj;
  Db *dbc;
  QChol *Qc;
  M2D_Environ *m2denv;
  SPDE_Option s_option;

  /* Initializations */

  error = 1;
  iatt_f = iatt_out = -1;
  ydat = ymean = yvert = rhs = zkrig = gwork = vwork = nullptr;
  ydat_loc = ymean_loc = yvert_loc = lwork = nullptr;
  dbc = nullptr;
  Qc = nullptr;
  Bproj = nullptr;
  m2denv = (M2D_Environ*) NULL;
  ysigma = 0.;
  number_hard = 0;

  /* Preliminary checks */

  if (dbin == nullptr)
  {
    messerr("The function requires an input Db argument");
    goto label_end;
  }
  if (dbout == nullptr)
  {
    messerr("The function requires an output Db argument");
    goto label_end;
  }
  ndim = model->getDimensionNumber();
  if (model->getVariableNumber() != 1)
  {
    messerr("This function should be called in the case of a single Model");
    messerr("In your case: %d\n", model->getVariableNumber());
    goto label_end;
  }
  if (nlayer <= 0)
  {
    messerr("This application requires the Number of Layers to be positive");
    goto label_end;
  }
  if (dbin->getIntervalNumber() < nlayer)
  {
    messerr("This application requires Lower and Upper variables");
    messerr("to be defined in the Input Db for each layer (nint=%d)",
            dbin->getIntervalNumber());
    goto label_end;
  }
  if (!is_grid(dbout))
  {
    messerr("This application is restricted to a Grid output Db");
    goto label_end;
  }
  if (ndim != 2)
  {
    messerr("This application is restricted to the 2-D case (ndim=%d)", ndim);
    goto label_end;
  }
  if (flag_ed && nlayer > dbout->getExternalDriftNumber())
  {
    messerr("External Drifts are used for Drift definition");
    messerr("- Count of F-variables (%d) must match Count of layers (%d)",
            dbout->getExternalDriftNumber(), nlayer);
    goto label_end;
  }
  if (nbsimu <= 0)
  {
    if (!flag_drift)
    {
      messerr("When 'nbsimu=0', the option 'flag.drift' is set to TRUE");
      messerr("Then the Optimal Drift is calculated only");
    }
    flag_drift = 1;
  }
  if (st_m2d_check_pinchout(dbout, icol_pinch)) goto label_end;

  law_set_random_seed(seed);
  ngrid = dbout->getSampleNumber();

  /* Prepare the M2D_Environ structure */

  vartot = model->getCovAnisoList()->getTotalSill(0, 0);

  m2denv = m2denv_manage(1, flag_ed, sqrt(vartot), NULL);
  if (m2denv == (M2D_Environ*) NULL) goto label_end;

  /* Preparing the variables in 'dbout' */

  nfois = (flag_drift) ? 1 :
                         nbsimu;
  iatt_out = dbout->addFields(nlayer * nfois, TEST);
  if (iatt_out < 0) goto label_end;

  /* Core allocation */

  lwork = (double*) mem_alloc(sizeof(double) * nlayer, 0);
  if (lwork == nullptr) goto label_end;

  /* Global statistics on Raw elevations */

  st_m2d_stats_init(m2denv, dbin, nlayer, verbose);

  /* Manage the Drift: define External Drift on input and output Db */

  if (verbose)
    message("\n==> Migrating Drift Information from Grid to Wells\n");
  if (st_m2d_drift_manage(m2denv, dbin, dbout, nlayer, verbose, &iatt_f))
    goto label_end;
  st_print_db_constraints("List of Initial Constraining Data", dbin, NULL,
                          nlayer, verbose);

  /* Constitute the new Db containing all the inequality constraints */
  /* whether they belong to 'dbin' or to 'dbout' */

  if (verbose)
    message("\n==> Creating a Temporary Data Base with all constraints\n");
  dbc = st_m2d_create_constraints(m2denv, dbin, dbout, ndim, nlayer,
                                  &number_hard);
  if (dbc == nullptr) goto label_end;
  nech = dbc->getActiveSampleNumber();

  /* Check SPDE environment */
  // At the first call, only one variable is Z_locatorized in order to
  // let the checks be performed on a mono-variate case (as all variables
  // will share the same Q matrix)
  // Then the environment is set to the multivariate case
  if (verbose) message("\n==> Checking SPDE Environment\n");
  st_define_locators(m2denv, dbc, ndim, 1, nlayer);
  if (spde_check(dbc, dbout, model, NULL, 0, VectorDouble(), 0, 1, 1, 0, 0, 0,
                 0)) goto label_end;
  st_define_locators(m2denv, dbc, ndim, nlayer, nlayer);

  /* Define initial values at constraints and set in Db */

  if (verbose) message("\n==> Creating Initial Value within bounds at Wells\n");
  if (st_m2d_initial_elevations(m2denv, dbc, nlayer, lwork)) goto label_end;

  /* Global statistics on Centered Elevations */

  st_m2d_stats_updt(m2denv, dbc, nlayer, verbose);

  /* Fitting the coefficients of the drift (external or not) */

  if (verbose) message("\n==> Fitting the optimal Drift(s)\n");
  if (st_m2d_drift_fitting(m2denv, dbc, nlayer, number_hard, verbose))
    goto label_end;

  /* Save the drift only (optional) */

  if (flag_drift)
  {
    gwork = (double*) mem_alloc(sizeof(double) * ngrid * nlayer, 0);
    if (gwork == nullptr) goto label_end;
    st_m2d_drift_save(m2denv, dbout, nlayer, gwork);
    for (int ilayer = 0; ilayer < nlayer; ilayer++)
    {
      dbout->setFieldByAttribute(&GWORK(ilayer, 0), iatt_out + ilayer);
      (void) gslSPrintf(string_encode, "Drift%d", ilayer + 1);
      db_name_set(dbout, iatt_out + ilayer, string_encode);
    }
    error = 0;
    goto label_end;
  }

  /**********************************************************************/
  /* From now on, the information is stored as drift increment          */
  /**********************************************************************/

  /* Manage Drift: */
  /* Drift (corrected by pinch-out) is stored in 'dbc' and 'dbout' */

  if (verbose) message("\n==> Transforming Drift information as Thickness\n");
  if (st_m2d_drift_inc_manage(m2denv, 1, nlayer, icol_pinch, dbc, dbout))
    goto label_end;

  /* Prepare all material */

  if (verbose) message("\n==> Preparing SPDE\n");
  s_option = spde_option_alloc();
  if (spde_prepar(dbc, dbout, VectorDouble(), s_option)) goto label_end;
  {
    SPDE_Matelem &Matelem = spde_get_current_matelem(0);
    MeshEStandard amesh;
    amesh.convertFromOldMesh(Matelem.s_mesh, 0);
    nvertex = st_get_nvertex(0);

    /* Core allocation */

    ydat = (double*) mem_alloc(sizeof(double) * nech * nlayer, 0);
    if (ydat == nullptr) goto label_end;
    ymean = (double*) mem_alloc(sizeof(double) * nech * nlayer, 0);
    if (ymean == nullptr) goto label_end;
    yvert = (double*) mem_alloc(sizeof(double) * nlayer * nvertex, 0);
    if (yvert == nullptr) goto label_end;
    rhs = (double*) mem_alloc(sizeof(double) * nvertex, 0);
    if (rhs == nullptr) goto label_end;
    vwork = (double*) mem_alloc(sizeof(double) * nvertex, 0);
    if (vwork == nullptr) goto label_end;
    zkrig = (double*) mem_alloc(sizeof(double) * nvertex, 0);
    if (zkrig == nullptr) goto label_end;

    for (int i = 0; i < nlayer * nvertex; i++)
      yvert[i] = 0.;
    for (int i = 0; i < nlayer * nech; i++)
      ymean[i] = 0.;

    /* Extract the vector of current data */

    if (verbose) message("\n==> Extracting the Initial Values at Wells\n");
    st_print_db_constraints("List of Initial Constraining Data", dbc, NULL,
                            nlayer, verbose);
    st_m2d_vector_extract(m2denv, dbc, nlayer, ydat, lwork);
    st_print_db_constraints("List of Constraining Data at Wells", dbc, ydat,
                            nlayer, verbose);
    st_m2d_stats_gaus("G-vect (initial)", nlayer, nech, ydat);

    /* Print environment just before entering in iterative process */

    if (verbose) st_m2d_print_environ("Environment before Simulations", m2denv);

    /* Loop on the simulations */

    for (int isimu = 0; isimu < nbsimu; isimu++)
    {
      message("Simulation #%d/%d\n", isimu + 1, nbsimu);

      /* Loop on Gibbs iterations */

      if (verbose)
        message("\n==> Launching the Simulations (%d iterations)\n", niter);
      nugget = vartot;
      for (int iter = 0; iter < niter; iter++)
      {
        if (verbose) message(">>>> Iteration #%d/%d\n", iter + 1, niter);

        // Update the Cholesky matrix

        if (iter == 0 && isimu == 0)
        {
          nugget /= 100.;
          ysigma = sqrt(nugget);
          Qc = st_derive_Qc(nugget, Qc, Matelem);
          if (Qc == nullptr) goto label_end;
        }

        // Perform the conditional simulation at meshing vartices

        for (int ilayer = 0; ilayer < nlayer; ilayer++)
        {
          ydat_loc = &YDAT(ilayer, 0);
          yvert_loc = &YVERT(ilayer, 0);
          ymean_loc = &YMEAN(ilayer, 0);

          // Non-conditional simulation

          st_simulate_cholesky(Qc, vwork, yvert_loc);
          for (int i = 0; i < nvertex; i++)
            yvert_loc[i] *= ysigma;

          // Conditional simulation

          for (int i = 0; i < nvertex; i++)
            zkrig[i] = vwork[i] = 0.;
          cs_tmulvec(Matelem.Aproj, nvertex, ydat_loc, rhs);
          st_kriging_cholesky(Qc, rhs, vwork, zkrig);
          for (int i = 0; i < nvertex; i++)
            yvert_loc[i] += zkrig[i];

          // Project the Simulation from the vertices onto the Data

          cs_mulvec(Matelem.Aproj, nech, yvert_loc, ymean_loc);
        }

        // Perform a Gibbs iteration on the constraints

        st_m2d_stats_gaus("G-Mean before Gibbs", nlayer, nech, ymean);
        if (st_global_gibbs(m2denv, dbc, 0, iter, nlayer, ysigma, ymean, ydat,
                            lwork)) goto label_end;
        st_m2d_stats_gaus("G-vect after Gibbs", nlayer, nech, ydat);
      }

      /* Check that the Constraints on the Wells are honored */

      if (verbose) message("\n==> Checking the Constraints at Wells\n");
      if (st_check_gibbs_data("Checking Constraints at Wells", m2denv, dbc,
                              nlayer, verbose, ydat, lwork)) goto label_end;
      st_m2d_stats_gaus("G-vect final", nlayer, nech, ydat);

      /* Store the conditional simulation on the grid */

      Bproj = db_mesh_sparse(dbout, &amesh, 0);
      if (Bproj == nullptr) goto label_end;
      gwork = (double*) mem_alloc(sizeof(double) * ngrid * nlayer, 0);
      if (gwork == nullptr) goto label_end;

      /* Project from vertices to grid nodes */

      for (int ilayer = 0; ilayer < nlayer; ilayer++)
        cs_mulvec(Bproj, ngrid, &YVERT(ilayer, 0), &GWORK(ilayer, 0));

      /* Convert from Gaussian to Depth */

      for (int igrid = 0; igrid < ngrid; igrid++)
      {
        if (!dbout->isActive(igrid)) continue;
        for (int ilayer = 0; ilayer < nlayer; ilayer++)
          lwork[ilayer] = GWORK(ilayer, igrid);
        st_convert_Y2Z(m2denv, dbout, nlayer, 2, igrid, lwork);
        for (int ilayer = 0; ilayer < nlayer; ilayer++)
          GWORK(ilayer,igrid) = lwork[ilayer];
      }

      st_m2d_stats_gaus("Depth on grid", nlayer, ngrid, gwork);
      for (int ilayer = 0; ilayer < nlayer; ilayer++)
      {
        dbout->setFieldByAttribute(&GWORK(ilayer, 0),
                                   iatt_out + isimu * nlayer + ilayer);
      }
    }

    // Renaming the simulation outcomes

    ecr = 0;
    for (int isimu = 0; isimu < nbsimu; isimu++)
    {
      for (int ilayer = 0; ilayer < nlayer; ilayer++)
      {
        (void) gslSPrintf(string_encode, "Layer-%d_Simu-%d", ilayer + 1,
                          isimu + 1);
        db_name_set(dbout, iatt_out + ecr, string_encode);
        ecr++;
      }
    }

    /* Convert the simulations into the mean and variance */

    if (flag_ce || flag_cstd)
    {
      // Modify the locator to ELoc::GAUSFAC before grouping to CE estimation

      dbout->setLocatorsByAttribute(nbsimu * nlayer, iatt_out, ELoc::GAUSFAC);

      if (db_simulations_to_ce(dbout, ELoc::GAUSFAC, nbsimu, nlayer, &iptr_ce,
                               &iptr_cstd)) goto label_end;

      // We release the attributes dedicated to simulations on Dbout

      if (!flag_ce)
      {
        (void) db_attribute_del_mult(dbout, iptr_ce, nlayer);
        iptr_ce = -1;
      }
      if (!flag_cstd)
      {
        (void) db_attribute_del_mult(dbout, iptr_cstd, nlayer);
        iptr_cstd = -1;
      }
      dbout->deleteFieldByLocator(ELoc::GAUSFAC);

      // Renaming the resulting variables

      if (iptr_ce >= 0) for (int ilayer = 0; ilayer < nlayer; ilayer++)
      {
        (void) gslSPrintf(string_encode, "Layer-%d_CE", ilayer + 1);
        db_name_set(dbout, iptr_ce + ilayer, string_encode);
      }
      if (iptr_cstd >= 0) for (int ilayer = 0; ilayer < nlayer; ilayer++)
      {
        (void) gslSPrintf(string_encode, "Layer-%d_CStd", ilayer + 1);
        db_name_set(dbout, iptr_cstd + ilayer, string_encode);
      }
    }
  }

  /* Set the error code */

  spde_posterior(dbc, dbout, VectorDouble(), s_option);
  error = 0;

  label_end: (void) st_m2d_drift_inc_manage(m2denv, -1, nlayer, icol_pinch, dbc,
                                            dbout);
  m2denv = m2denv_manage(-1, flag_ed, 0., m2denv);
  Qc = qchol_manage(-1, Qc);
  Bproj = cs_spfree(Bproj);
  ydat = (double*) mem_free((char* ) ydat);
  ymean = (double*) mem_free((char* ) ymean);
  yvert = (double*) mem_free((char* ) yvert);
  gwork = (double*) mem_free((char* ) gwork);
  vwork = (double*) mem_free((char* ) vwork);
  lwork = (double*) mem_free((char* ) lwork);
  rhs = (double*) mem_free((char* ) rhs);
  zkrig = (double*) mem_free((char* ) zkrig);
  if (iatt_f >= 0) (void) db_attribute_del_mult(dbin, iatt_f, nlayer);
  if (error && iatt_out >= 0)
    (void) db_attribute_del_mult(dbout, iatt_out, nlayer);
  return (error);
}
