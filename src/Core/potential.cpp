/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://github.com/gstlearn                                       */
/* License: BSD 3 clauses                                                     */
/*                                                                            */
/******************************************************************************/
#include "geoslib_f.h"
#include "geoslib_old_f.h"

#include "Enum/EDrift.hpp"

#include "Basic/Utilities.hpp"
#include "Basic/Law.hpp"
#include "Basic/File.hpp"
#include "Basic/String.hpp"
#include "Basic/OptDbg.hpp"
#include "Covariances/CovAniso.hpp"
#include "Covariances/CovLMGradient.hpp"
#include "Drifts/DriftList.hpp"
#include "Db/Db.hpp"
#include "Db/DbGrid.hpp"
#include "Model/Model.hpp"
#include "Neigh/ANeigh.hpp"
#include <Simulation/CalcSimuTurningBands.hpp>

#include <math.h>
#include <string.h>

/*! \cond */

#define DRF(il)          (TAB_DRF[il])
#define MAT(i,j)         (mat[(i) * nequa + (j)])
#define B(isol,i)        (b[(isol) * number + (i)])
#define POTVAL(isimu,il) (potval[(isimu) * pot_env->nlayers + (il)])
#define POTSIM(isimu,il) (potsim[(isimu) * nlayers + (il)])
#define ZDUALS(isimu,i)  (zduals[(isimu) * nequa + (i)])

typedef struct
{
  int ndim; /* Space dimension */
  int niso; /* Number of Iso-potential information */
  int nlayers; /* Number of Iso-potential values */
  int ngrd; /* Number of gradient information */
  int ntgt; /* Number of tangent information */
  int next; /* Number of external drifts */
  int nequa; /* Number of equations in the System */
  int order; /* Order of the drift */
  int size_iso; /* Matrix size linked to iso-potential */
  int size_grd; /* Matrix size linked to gradient */
  int size_tgt; /* Matrix size linked to tangent */
  int size_drf; /* Matrix size linked to Drift functions */
  int size_ext; /* Matrix size linked to External Drifts */
  int start_iso; /* Address of the first iso-potential */
  int start_grd; /* Address of the first gradient */
  int start_tgt; /* Address of the first tangent */
  int start_drf; /* Address of the first drift */
  int start_ext; /* Address of the first external drift */
  VectorInt nb_per_layer; /* Array of counts of samples per layer */
  VectorInt ptr_per_layer; /* Array of ptr per layer */
  VectorInt rank_iso; /* Array of ranks for Iso-potential */
  VectorInt rank_grd; /* Array of ranks for Gradients */
  VectorInt rank_tgt; /* Array of ranks for Tangents */
  int  opt_part;
  bool flag_pot;
  bool flag_grad;
  bool flag_trans;
} Pot_Env;

typedef struct
{
  int ndim;
  int nring;
  int nfull;
  double range;
  DbGrid *db;
  Model *model;
  int *indg;
  int *indg0;
  double *data;
  double *weight;
} Pot_Ext;

static int TAB_DRF[9];
static bool VERBOSE = false;
static Pot_Env* POTENV = nullptr;
static Pot_Ext* POTEXT = nullptr;
static Db* DBISO = nullptr;
static Db* DBGRD = nullptr;
static Db* DBTGT = nullptr;

/*! \endcond */

static void st_potenv_define(Pot_Env* pot_env,
                             Pot_Ext* pot_ext,
                             Db* dbiso,
                             Db* dbgrd,
                             Db* dbtgt,
                             Db* dbout)
{
  POTENV = pot_env;
  POTEXT = pot_ext;
  DBISO = dbiso;
  DBGRD = dbgrd;
  DBTGT = dbtgt;

  set_DBIN(dbiso);
  set_DBOUT(dbout);

  pot_env->ndim = dbiso->getNDim();
}

static int GRX(int i)
{
  if (POTENV->ndim < 1)
    return -1;
  else
    return i;
}
static int GRY(int i)
{
  if (POTENV->ndim < 2)
    return -1;
  else
    return i + POTENV->ngrd;
}
static int GRZ(int i)
{
  if (POTENV->ndim < 3)
    return -1;
  else
    return i + 2 * POTENV->ngrd;
}
static int TGT(int i)
{
  return POTENV->start_tgt + i;
}
static int ISC(int ic,int i)
{
  return POTENV->start_iso + POTENV->ptr_per_layer[ic] + (i) - (ic) - 1;
}
static int IAD_GRD(int ig)
{
  return POTENV->rank_grd[ig];
}
static void set_IAD_GRD(int ig, int value)
{
  POTENV->rank_grd[ig] = value;
}
static int IAD_TGT(int it)
{
  return POTENV->rank_tgt[it];
}
static void set_IAD_TGT(int it, int value)
{
  POTENV->rank_tgt[it] = value;
}
static double IAD_ISO(int ic,int i)
{
  return POTENV->rank_iso[POTENV->ptr_per_layer[ic] + (i)];
}
static double TGT_COO(int it,int i)
{
  return DBTGT->getCoordinate(IAD_TGT(it),i);
}
static double TGT_VAL(int it,int idim)
{
  if (idim >= POTENV->ndim) return TEST;
  return DBTGT->getLocVariable(ELoc::TGTE,IAD_TGT(it),idim);
}
static double GRD_COO(int ig,int idim)
{
  if (idim >= POTENV->ndim) return TEST;
  return DBGRD->getCoordinate(IAD_GRD(ig),idim);
}
static double GRD_VAL(int ig,int idim)
{
  if (idim >= POTENV->ndim) return TEST;
  return DBGRD->getLocVariable(ELoc::G,IAD_GRD(ig),idim);
}
static double ISO_COO(int ic,int j,int idim)
{
  if (idim >= POTENV->ndim) return TEST;
  return DBISO->getCoordinate(IAD_ISO(ic,j),idim);
}
static int EXT(int iext)
{
  return POTENV->start_ext + (iext);
}

/****************************************************************************/
/*!
 **  Cehck if the Model can be used for Potentials
 **
 ** \return  Error return code
 **
 ** \param[in]  model      Model structure
 **
 *****************************************************************************/
static int st_model_invalid(Model *model)

{
  for (int icov = 0; icov < model->getCovaNumber(); icov++)
  {
    ECov type = model->getCovaType(icov);
    if (type != ECov::GAUSSIAN && type != ECov::CUBIC
        && type != ECov::SPLINE2_GC && type != ECov::NUGGET)
    {
      messerr("The Model is invalid for Potential calculations");
      messerr("It may only contain:");
      messerr("- Cubic covariance");
      messerr("- Gaussian covariance");
      messerr("- Duchon Spline generalized covariance");
      messerr("An additional nugget effect can also be considered");
      return (1);
    }
    if (type == ECov::SPLINE2_GC && model->getMaximumOrder() < 2)
    {
      messerr("The Model includes Second Order Spline Generalized Covariance");
      messerr("This requires a second order drift");
      return (1);
    }
  }
  return (0);
}

/****************************************************************************/
/*!
 **  Calculate the covariance and the derivatives
 **
 ** \param[in] model     Model structure
 ** \param[in] flag_grad True if the Gradients must be calculated
 ** \param[in] dx        Increment along X
 ** \param[in] dy        Increment along Y
 ** \param[in] dz        Increment along Z
 **
 ** \param[out] covar    Covariance of potential
 ** \param[out] covGp    Covariance between potential and gradient
 ** \param[out] covGG    Covariance between gradient and gradient
 **
 *****************************************************************************/
static void st_cov(Model* model,
                   bool flag_grad,
                   double dx,
                   double dy,
                   double dz,
                   double& covar,
                   VectorDouble& covGp,
                   VectorDouble& covGG)
{
  int ndim = model->getDimensionNumber();
  VectorDouble vec(ndim);
  if (ndim >= 1) vec[0] = dx;
  if (ndim >= 2) vec[1] = dy;
  if (ndim >= 3) vec[2] = dz;

  model->evalZAndGradients(vec, covar, covGp, covGG, CovCalcMode(), flag_grad);

  return;
}

/****************************************************************************/
/*!
 **  Establish the neighborhood data Db file
 **
 ** \return  Error return code
 **
 ** \param[in]  dbout      Output Db structure
 **
 ** \param[out] pot_ext     Pot_Ext structure
 **
 *****************************************************************************/
static int st_extdrift_create_db(DbGrid *dbout, Pot_Ext *pot_ext)
{
  int error, nech;
  VectorInt nx;
  VectorDouble x0;

  /* Initializations */

  error = 1;
  nech = 0;

  /* Core allocation */

  nx.resize(pot_ext->ndim);
  x0.resize(pot_ext->ndim);

  /* Creating the attributes from the output grid */

  nech = 1;
  for (int idim = 0; idim < pot_ext->ndim; idim++)
  {
    nx[idim] = 2 * pot_ext->nring + 1;
    x0[idim] = -dbout->getDX(idim) * pot_ext->nring;
    nech *= nx[idim];
  }

  /* Creating the data grid */

  pot_ext->db = DbGrid::create(nx, dbout->getDXs(), x0, dbout->getAngles(),
                               ELoadBy::COLUMN, VectorDouble(),
                               VectorString(), VectorString(), 1);
  if (pot_ext->db == nullptr) goto label_end;
  pot_ext->nfull = nech;

  /* Add the selection */

  pot_ext->db->addColumnsByConstant(1, 0., String(), ELoc::SEL);

  /* Complementary core allocation */

  pot_ext->data = (double*) mem_alloc(sizeof(double) * nech, 0);
  if (pot_ext->data == nullptr) goto label_end;
  pot_ext->weight = (double*) mem_alloc(sizeof(double) * nech * 4, 0);
  if (pot_ext->weight == nullptr) goto label_end;
  pot_ext->indg0 = (int*) mem_alloc(sizeof(int) * 3, 0);
  if (pot_ext->indg0 == nullptr) goto label_end;
  pot_ext->indg = (int*) mem_alloc(sizeof(int) * 3, 0);
  if (pot_ext->indg == nullptr) goto label_end;

  /* Set the error return code */

  error = 0;

  /* Returning arguments */

  label_end: if (error)
  {
    if (pot_ext->db != nullptr) delete pot_ext->db;
    if (pot_ext->data != nullptr) delete pot_ext->data;
    if (pot_ext->weight != nullptr) delete pot_ext->weight;
  }
  return (error);
}

/****************************************************************************/
/*!
 **  Establish the Model
 **
 ** \return  Error return code
 **
 ** \param[out] pot_ext     Pot_Ext structure
 **
 *****************************************************************************/
static int st_extdrift_create_model(Pot_Ext *pot_ext)
{
  double sill = 1.;

  /* Creating the model */

  CovContext ctxt = CovContext(1, pot_ext->ndim, 1.);
  pot_ext->model = new Model(ctxt);
  if (pot_ext->model == nullptr) return 1;

  // Covariance part
  CovLMGradient covs(ctxt.getSpace());
  CovAniso cov(ECov::CUBIC, pot_ext->range, 0., sill, ctxt);
  covs.addCov(&cov);
  pot_ext->model->setCovList(&covs);

  // Drift part
  DriftList drifts(ctxt.getSpace());
  drifts.setFlagLinked(true);
  pot_ext->model->setDriftList(&drifts);

  return 0;
}

/****************************************************************************/
/*!
 **  Establish kriging system
 **
 ** \return  Error return code
 **
 ** \param[out] pot_ext    Pot_Ext structure
 ** \param[out] number     Matrix dimension
 ** \param[out] a          LHS matrix
 ** \param[out] b          RHS matrix
 ** \param[out] wgt        Matrix of weights
 **
 *****************************************************************************/
static int st_extdrift_solve(Pot_Ext *pot_ext,
                             int number,
                             double *a,
                             double *b,
                             double *wgt)
{
  int ecr;
  double covar = 0.;
  VectorDouble covGp(3, 0.);
  VectorDouble covGG(9, 0.);

  /* Establish the kriging matrix */

  model_covmat(pot_ext->model, pot_ext->db, pot_ext->db, -1, -1, 0, 1, a);

  /* Establish the Right-Hand side */

  ecr = 0;
  for (int iech = 0; iech < pot_ext->nfull; iech++)
  {
    if (!pot_ext->db->isActive(iech)) continue;
    st_cov(pot_ext->model, 1,
           pot_ext->db->getCoordinate(iech, 0),
           pot_ext->db->getCoordinate(iech, 1),
           pot_ext->db->getCoordinate(iech, 2),
           covar, covGp, covGG);
    B(0,ecr) = covar;
    B(1,ecr) = -covGp[0];
    B(2,ecr) = -covGp[1];
    B(3,ecr) = -covGp[2];
    ecr++;
  }

  /* Perform the inversion and store the weights */

  if (matrix_invert(a, number, 0)) return (1);
  matrix_product_safe(number, number, 4, a, b, wgt);

  return (0);
}

/****************************************************************************/
/*!
 **  Establish kriging system for calculating Gradient on External Drift
 **
 ** \return  Error return code
 **
 ** \param[in]  dbout      Output Db structure
 **
 ** \param[out] pot_ext    Pot_Ext structure
 **
 *****************************************************************************/
static int st_extdrift_calc_init(DbGrid *dbout, Pot_Ext *pot_ext)
{
  int number, error;
  double *a, *b;

  /* Initializations */

  error = 1;
  a = b = nullptr;

  /* Creating the Db for neighborhood */

  if (st_extdrift_create_db(dbout, pot_ext)) goto label_end;
  number = pot_ext->nfull;

  /* Creating the model */

  if (st_extdrift_create_model(pot_ext)) goto label_end;

  /* Core allocation */

  a = (double*) mem_alloc(sizeof(double) * number * number, 0);
  if (a == nullptr) goto label_end;
  b = (double*) mem_alloc(sizeof(double) * number * 4, 0);
  if (b == nullptr) goto label_end;

  /* Solve the kriging system */

  if (st_extdrift_solve(pot_ext, number, a, b, pot_ext->weight)) goto label_end;

  /* Set the error return code */

  error = 0;

  label_end:

  /* Core deallocation */

  a = (double*) mem_free((char* ) a);
  b = (double*) mem_free((char* ) b);
  return (error);
}

/****************************************************************************/
/*!
 **  Manage the Pot_Ext structure
 **
 ** \param[in]  mode      1 for allocation; 0 for save; -1 for deallocation
 ** \param[in]  pot_ext   Pot_Ext structure to be managed
 ** \param[in]  nring     Number of rings used for Ext. Drift calculations
 ** \param[in]  range     Range of the structure
 ** \param[in]  dbout     Output Db structure
 **
 *****************************************************************************/
static int st_potext_manage(int mode,
                            Pot_Ext *pot_ext,
                            int nring,
                            double range,
                            DbGrid *dbout)
{
  /* Dispatch */

  switch (mode)
  {
    case 0: /* Initialization */
      pot_ext->ndim = 0;
      pot_ext->nring = 0;
      pot_ext->nfull = 0;
      pot_ext->range = 0.;
      pot_ext->db = nullptr;
      pot_ext->model = nullptr;
      pot_ext->indg = nullptr;
      pot_ext->indg0 = nullptr;
      pot_ext->data = nullptr;
      pot_ext->weight = nullptr;
      return (0);

    case 1: /* Allocation */
      pot_ext->ndim = dbout->getNDim();
      pot_ext->nring = nring;
      pot_ext->range = range;
      if (st_extdrift_calc_init(dbout, pot_ext)) return (1);
      return (0);

    case -1: /* Deletion */
      delete pot_ext->db;
      pot_ext->db = nullptr;
      delete pot_ext->model;
      pot_ext->model = nullptr;
      pot_ext->indg = db_indg_free(pot_ext->indg);
      pot_ext->indg0 = db_indg_free(pot_ext->indg0);
      pot_ext->data = (double*) mem_free((char* ) pot_ext->data);
      pot_ext->weight = (double*) mem_free((char* ) pot_ext->weight);
      return (0);
  }
  return (1);
}

bool st_potenv_valid(Pot_Env* pot_env,
                     Pot_Ext* pot_ext,
                     Db *dbiso,
                     Db *dbgrd,
                     Db *dbtgt,
                     DbGrid* dbout,
                     Model *model,
                     ANeigh *neigh)
{
  static int nring = 1;

  if (pot_env->ndim > 3)
  {
    messerr("The input Db must be defined in Space with dimension < 3");
    return false;
  }
  if (dbgrd != nullptr && dbgrd->getNDim() != pot_env->ndim)
  {
    messerr("The Gradient and Data Db must share the same space dimension");
    return false;
  }
  if (dbtgt != nullptr && dbtgt->getNDim() != pot_env->ndim)
  {
    messerr("The Tangent and Data Db must share the same space dimension");
    return false;
  }
  if (model->getDimensionNumber() != pot_env->ndim)
  {
    messerr("The Model and Data Db must have the same space dimension");
    return false;
  }
  if (dbout != NULL && dbout->getNDim() != pot_env->ndim)
  {
    messerr("The Db files 'dbin' and 'dbout' should have the same dimension");
    return false;
  }
  if (st_model_invalid(model)) return false;
  if (!exist_LOCATOR(dbiso, ELoc::LAYER))
  {
    messerr("The input Db must contain a LAYER locator");
    return false;
  }
  if (model->getVariableNumber() != 1)
  {
    messerr("The Model must be monovariate");
    return false;
  }
  if (neigh->getNeighParam()->getType() != ENeigh::UNIQUE)
  {
    messerr("This procedure is only available in Unique Neighborhood");
    return false;
  }

  int next = model_nfex(model);
  if (dbout != NULL && next != dbout->getLocNumber(ELoc::F))
  {
    messerr("Inconsistency for External Drift between Model and Dbout");
    return false;
  }
  if (dbout == NULL && next > 0)
  {
    messerr("Usage of External drift is forbidden without Output Grid");
    return false;
  }
  if (next > 0)
  {
    if (next > 1)
    {
      messerr("This application cannot deal with more than 1 External Drift");
      messerr("Check your output file");
      return false;
    }
    if (! dbout->isGrid())
    {
      messerr("The External Drift requires an Output Grid File");
      return false;
    }
    double range = 3. * MAX(dbout->getDX(0), dbout->getDX(1));
    if (st_potext_manage(1, pot_ext, nring, range, dbout)) return false;
  }

  return true;
}

/****************************************************************************/
/*!
 **  Manage the Pot_Env structure
 **
 ** \param[in,out]  pot_env    Pot_env structure
 ** \param[in]      flag_pot   True if the potential must be calculated
 ** \param[in]      flag_grad  True if the gradients must be calculated
 ** \param[in]      flag_trans True if the estimation result must be translated
 **                            into layer number
 ** \param[in]      opt_part   Option to exhibit only a part of estimation:
 ** \li                        0 : the whole estimation
 ** \li                        1 : the gradient contribution only
 ** \li                        2 : the tangent contribution only
 ** \li                        3 : the isovalues contribution only
 ** \li                        4 : the drift contribution only
 ** \li                        5 : the external drift contribution only
 ** \param[in]      verbose    Verbose flag
 **
 *****************************************************************************/
static void st_potenv_manage(Pot_Env *pot_env,
                             bool flag_pot,
                             bool flag_grad,
                             bool flag_trans,
                             int opt_part,
                             bool verbose)
{
  VERBOSE = verbose;
  if (opt_part) flag_trans = false;

  pot_env->ndim = 0;
  pot_env->niso = 0;
  pot_env->nlayers = 0;
  pot_env->ngrd = 0;
  pot_env->ntgt = 0;
  pot_env->next = 0;
  pot_env->nequa = 0;
  pot_env->order = 0;
  pot_env->size_iso = 0;
  pot_env->size_grd = 0;
  pot_env->size_tgt = 0;
  pot_env->size_drf = 0;
  pot_env->size_ext = 0;
  pot_env->start_iso = 0;
  pot_env->start_grd = 0;
  pot_env->start_tgt = 0;
  pot_env->start_drf = 0;
  pot_env->start_ext = 0;
  pot_env->nb_per_layer = VectorInt();
  pot_env->ptr_per_layer = VectorInt();
  pot_env->rank_iso = VectorInt();
  pot_env->rank_grd = VectorInt();
  pot_env->rank_tgt = VectorInt();
  pot_env->flag_pot = flag_pot;
  pot_env->flag_grad = flag_grad;
  pot_env->flag_trans = flag_trans;
  pot_env->opt_part = opt_part;
}

/****************************************************************************/
/*!
 **  Load the information linked to the Iso_potentials
 **
 ** \return  Error return code
 **
 ** \param[in]  dbiso      Input Db structure for Iso-Potential
 ** \param[in,out] pot_env The Pot_Env structure
 **
 *****************************************************************************/
static int st_update_isopot(Db *dbiso, Pot_Env *pot_env)

{
  int *layval, *laycnt, error;
  int i, ival, nech, nlayers, niso, ipos, j, ecr, found;
  double value;

  // Initializations

  if (dbiso == nullptr) return (0);
  error = 1;
  nech = dbiso->getSampleNumber();
  nlayers = niso = 0;
  layval = laycnt = nullptr;

  // Count the number of different iso-potential values 

  for (int iech = 0; iech < nech; iech++)
  {
    if (!dbiso->isActive(iech)) continue;
    value = get_LOCATOR_ITEM(dbiso, ELoc::LAYER, 0, iech);
    if (FFFF(value)) continue;
    ival = (int) value;

    // Look for an already registered layer value 

    for (i = 0, found = -1; i < nlayers && found < 0; i++)
      if (ival == layval[i]) found = i;

    if (found < 0)
    {
      layval = (int*) mem_realloc((char* ) layval, sizeof(int) * (nlayers + 1),
                                  1);
      laycnt = (int*) mem_realloc((char* ) laycnt, sizeof(int) * (nlayers + 1),
                                  1);
      layval[nlayers] = ival;
      laycnt[nlayers] = 1;
      nlayers++;
    }
    else
    {
      laycnt[found]++;
    }
  }

  // Eliminate layers with not enough samples 

  niso = 0;
  for (i = j = 0; i < nlayers; i++)
  {
    if (laycnt[i] < 2) continue;
    layval[j] = layval[i];
    laycnt[j] = laycnt[i];
    niso += laycnt[i];
    j++;
  }
  pot_env->nlayers = nlayers = j;
  pot_env->niso = niso;
  pot_env->size_iso = niso - nlayers;

  // Core allocation 

  pot_env->nb_per_layer.resize(nlayers);
  pot_env->ptr_per_layer.resize(nlayers);
  pot_env->rank_iso.resize(niso);
  if (pot_env->rank_iso.empty()) goto label_end;

  // Set the final length and pointers 

  ipos = 0;
  for (i = 0; i < nlayers; i++)
  {
    pot_env->nb_per_layer[i] = laycnt[i];
    pot_env->ptr_per_layer[i] = ipos;
    ipos += pot_env->nb_per_layer[i];
  }

  // Sort the samples per iso-potential value

  for (i = ecr = 0; i < nlayers; i++)
  {
    for (int iech = 0; iech < nech; iech++)
    {
      if (!dbiso->isActive(iech)) continue;
      value = get_LOCATOR_ITEM(dbiso, ELoc::LAYER, 0, iech);
      if (FFFF(value)) continue;
      ival = (int) value;
      if (ival != layval[i]) continue;
      pot_env->rank_iso[ecr++] = iech;
    }
  }

  // Reading failure

  if (niso < 1 || nlayers < 1)
  {
    messerr("The number of iso-potential informations cannot be null");
    goto label_end;
  }

  // Set the error return code

  error = 0;

  label_end: layval = (int*) mem_free((char* ) layval);
  laycnt = (int*) mem_free((char* ) laycnt);
  return (error);
}

/****************************************************************************/
/*!
 **  Load the information linked to the Gradients
 **
 ** \return  Array of indices of the active gradient information
 **
 ** \param[in]  dbgrd      Input Db structure
 ** \param[in,out] pot_env The Pot_Env structure
 **
 *****************************************************************************/
static int st_update_gradient(Db *dbgrd, Pot_Env *pot_env)
{
  int found, ngrd, nech;

  if (dbgrd == nullptr) return (0);
  nech = dbgrd->getSampleNumber();
  ngrd = 0;

  // Core allocation 

  pot_env->rank_grd.resize(nech);

  // Loop on the gradients

  for (int iech = 0; iech < nech; iech++)
  {
    if (!dbgrd->isActive(iech)) continue;
    found = 0;
    for (int idim = 0; idim < pot_env->ndim && found == 0; idim++)
      if (FFFF(dbgrd->getLocVariable(ELoc::G,iech, idim))) found = 1;
    if (found) continue;
    set_IAD_GRD(ngrd++,iech);
  }

  // Core reallocation

  pot_env->rank_grd.resize(ngrd);
  pot_env->ngrd = ngrd;
  pot_env->size_grd = ngrd * pot_env->ndim;

  if (ngrd < 1)
  {
    messerr("The number of gradient informations cannot be null");
    return (1);
  }

  return (0);
}

/****************************************************************************/
/*!
 **  Load the information linked to the Tangents
 **
 ** \return  Error return code
 **
 ** \param[in]  dbtgt      Input Db structure
 ** \param[in,out] pot_env The Pot_Env structure
 **
 *****************************************************************************/
static int st_update_tangent(Db *dbtgt, Pot_Env *pot_env)

{
  int found, ntgt, nech;

  if (dbtgt == nullptr) return (0);
  nech = dbtgt->getSampleNumber();
  ntgt = 0;

  // Core allocation 

  pot_env->rank_tgt.resize(nech);

  // Loop on the tangents

  for (int iech = 0; iech < nech; iech++)
  {
    if (!dbtgt->isActive(iech)) continue;
    found = 0;
    for (int idim = 0; idim < pot_env->ndim && found == 0; idim++)
      if (FFFF(dbtgt->getLocVariable(ELoc::TGTE,iech, idim))) found = 1;
    if (found) continue;
    set_IAD_TGT(ntgt++, iech);
  }

  // Core reallocation

  pot_env->rank_tgt.resize(ntgt);
  pot_env->ntgt = ntgt;
  pot_env->size_tgt = ntgt;

  return (0);
}

/****************************************************************************/
/*!
 **  Load the information linked to the Model
 **
 ** \return  Error return code
 **
 ** \param[in]  model      Model structure
 ** \param[in,out] pot_env The Pot_Env structure
 **
 *****************************************************************************/
static int st_update_model(Model *model, Pot_Env *pot_env)
{
  int nbfl;

  nbfl = model->getDriftNumber();
  if (model->isDriftDefined(EDrift::UC)) nbfl--;
  pot_env->order =  model->getMaximumOrder();
  pot_env->size_drf = nbfl;
  pot_env->next = pot_env->size_ext = model_nfex(model);

  return (0);
}

/****************************************************************************/
/*!
 **  Make the final checks and define the addresses
 **
 ** \return  Error return code
 **
 ** \param[in]  model      Model structure
 ** \param[in,out] pot_env The Pot_Env structure
 **
 *****************************************************************************/
static int st_update_final(Model *model, Pot_Env *pot_env)

{
  int pos;

  // Compute the starting addresses

  pos = 0;
  pot_env->start_grd = pos;
  pos += pot_env->size_grd;
  pot_env->start_tgt = pos;
  pos += pot_env->size_tgt;
  pot_env->start_iso = pos;
  pos += pot_env->size_iso;
  pot_env->start_drf = pos;
  pos += pot_env->size_drf;
  pot_env->start_ext = pos;

  // Compute the number of equations in the CoKriging System

  pot_env->nequa = (pot_env->size_iso + pot_env->size_grd + pot_env->size_tgt
                    + pot_env->size_drf + pot_env->size_ext);

  // Define the addresses for the drift functions

  pos = pot_env->start_drf;
  for (int i = 0; i < 9; i++)
    TAB_DRF[i] = -1;

  if (model->isDriftDefined(EDrift::X)) TAB_DRF[0] = pos++;
  if (model->isDriftDefined(EDrift::Y)) TAB_DRF[1] = pos++;
  if (model->isDriftDefined(EDrift::Z)) TAB_DRF[2] = pos++;
  if (model->isDriftDefined(EDrift::X2)) TAB_DRF[3] = pos++;
  if (model->isDriftDefined(EDrift::Y2)) TAB_DRF[4] = pos++;
  if (model->isDriftDefined(EDrift::Z2)) TAB_DRF[5] = pos++;
  if (model->isDriftDefined(EDrift::XY)) TAB_DRF[6] = pos++;
  if (model->isDriftDefined(EDrift::XZ)) TAB_DRF[7] = pos++;
  if (model->isDriftDefined(EDrift::YZ)) TAB_DRF[8] = pos++;

  /* Optional output */

  if (VERBOSE)
  {
    mestitle(0, "Environment summary");
    message("Space dimension         = %d\n", pot_env->ndim);
    message("Number of Iso-Potential = %d\n", pot_env->nlayers);
    message("Number of Gradients     = %d\n", pot_env->ngrd);
    message("Number of Tangents      = %d\n", pot_env->ntgt);
    message("Number of Isovalues     = %d\n", pot_env->niso);
    message("Order of the drift      = %d\n", pot_env->order);
    message("Number of Drifts        = %d\n", pot_env->size_drf);
    message("Number of Ext. Drifts   = %d\n", pot_env->size_ext);
    message("Number of Equations     = %d\n", pot_env->nequa);
  }

  return (0);
}

/****************************************************************************/
/*!
 **  Calculate the inner product of two vectors
 **
 **  \param[in]   ndim   : Space dimension
 **  \param[in]   ux     : First coordinate of the first vector
 **  \param[in]   uy     : Second coordinate of the first vector
 **  \param[in]   uz     : Third coordinate of the first vector
 **  \param[in]   vx     : First coordinate of the second vector
 **  \param[in]   vy     : Second coordinate of the second vector
 **  \param[in]   vz     : Third coordinate of the second vector
 **
 *****************************************************************************/
static double matrix_UV(int ndim,
                        double ux,
                        double uy,
                        double uz,
                        double vx,
                        double vy,
                        double vz)
{
  double prod;

  prod = 0.;
  if (ndim >= 1 && !FFFF(ux) && !FFFF(vx)) prod += ux * vx;
  if (ndim >= 2 && !FFFF(uy) && !FFFF(vy)) prod += uy * vy;
  if (ndim >= 3 && !FFFF(uz) && !FFFF(vz)) prod += uz * vz;
  return (prod);
}

/****************************************************************************/
/*!
 **  Calculate the norm product of two vectors by a matrix
 **
 **  \param[in]   ndim   : Space dimension
 **  \param[in]   a      : Matrix
 **  \param[in]   ux     : First coordinate of the first vector
 **  \param[in]   uy     : Second coordinate of the first vector
 **  \param[in]   uz     : Third coordinate of the first vector
 **  \param[in]   vx     : First coordinate of the second vector
 **  \param[in]   vy     : Second coordinate of the second vector
 **  \param[in]   vz     : Third coordinate of the second vector
 **
 *****************************************************************************/
static double matrix_UAV(int ndim,
                         double *a,
                         double ux,
                         double uy,
                         double uz,
                         double vx,
                         double vy,
                         double vz)
{
  double prod;

  prod = 0.;
  if (ndim >= 1 && !FFFF(ux) && !FFFF(vx)) prod += ux * vx * a[0];
  if (ndim >= 2 && !FFFF(ux) && !FFFF(vy)) prod += ux * vy * a[1];
  if (ndim >= 3 && !FFFF(ux) && !FFFF(vz)) prod += ux * vz * a[2];
  if (ndim >= 2 && !FFFF(uy) && !FFFF(vx)) prod += uy * vx * a[3];
  if (ndim >= 2 && !FFFF(uy) && !FFFF(vy)) prod += uy * vy * a[4];
  if (ndim >= 3 && !FFFF(uy) && !FFFF(vz)) prod += uy * vz * a[5];
  if (ndim >= 3 && !FFFF(uz) && !FFFF(vx)) prod += uz * vx * a[6];
  if (ndim >= 3 && !FFFF(uz) && !FFFF(vy)) prod += uz * vy * a[7];
  if (ndim >= 3 && !FFFF(uz) && !FFFF(vz)) prod += uz * vz * a[8];
  return (prod);
}

/****************************************************************************/
/*!
 **  Set and element in the Kriging R.H.S. vector
 **
 ** \param[in] rhs      Vector to be filled
 ** \param[in] nequa    Number of equations
 ** \param[in] i        Row number
 ** \param[in] isol     Column number
 ** \param[in] value    Value to be assigned to this cell
 **
 *****************************************************************************/
static void set_rhs(double *rhs, int nequa, int i, int isol, double value)
{
  if (i < 0 || isol < 0) return;
  if (i >= nequa)
    messageAbort("Fatal error in set_rhs(%d,%d) with nequa=%d", i, isol, nequa);
  rhs[(isol) * nequa + (i)] = value;
}

/****************************************************************************/
/*!
 **  Set and element in the Kriging L.H.S. matrix
 **
 ** \param[in] lhs      Matrix to be filled
 ** \param[in] nequa    Number of equations
 ** \param[in] i        Row number
 ** \param[in] j        Column number
 ** \param[in] value    Value to be assigned to this cell
 **
 *****************************************************************************/
static void set_lhs(double *lhs, int nequa, int i, int j, double value)
{
  if (i < 0 || j < 0) return;
  if (i >= nequa || j >= nequa)
    messageAbort("Fatal error in set_lhs(%d,%d) with nequa=%d", i, j, nequa);
  lhs[(i) * nequa + (j)] = value;
  if (i != j)
    lhs[(j) * nequa + (i)] = value;
}

/****************************************************************************/
/*!
 **  Get one element from the Kriging L.H.S. matrix
 **
 ** \return The returned value
 **
 ** \param[in] lhs      Matrix to be filled
 ** \param[in] nequa    Number of equations
 ** \param[in] i        Row number
 ** \param[in] j        Column number
 **
 *****************************************************************************/
static double get_lhs(double *lhs, int nequa, int i, int j)
{
  if (i < 0 || j < 0) return (0.);
  if (i >= nequa || j >= nequa)
    messageAbort("Fatal error in get_lhs(%d,%d) with nequa=%d", i, j, nequa);
  return (lhs[(i) * nequa + (j)]);
}

/****************************************************************************/
/*!
 **  Establish the local neighborhood
 **
 ** \return  Error return code (Neighborhood not complete)
 **
 ** \param[in]  dbgrid     Output Db Grid structure
 ** \param[in]  pot_ext    Pot_Ext structure
 **
 *****************************************************************************/
static int st_extdrift_neigh(DbGrid *dbgrid, Pot_Ext *pot_ext)
{
  int ecr, iech;
  double drift;

  /* Loop on the neighboring samples defined in the neighboring grid */

  ecr = 0;
  for (int iz = 0; iz < pot_ext->db->getNX(2); iz++)
    for (int iy = 0; iy < pot_ext->db->getNX(1); iy++)
      for (int ix = 0; ix < pot_ext->db->getNX(0); ix++)
      {

        /* Calculate the index of the sample within the Ext Drift grid */

        pot_ext->indg[0] = pot_ext->indg0[0] + ix - pot_ext->nring;
        if (pot_ext->indg[0] < 0 || pot_ext->indg[0] > dbgrid->getNX(0))
          return (1);
        pot_ext->indg[1] = pot_ext->indg0[1] + iy - pot_ext->nring;
        if (pot_ext->indg[1] < 0 || pot_ext->indg[1] > dbgrid->getNX(1))
          return (1);
        pot_ext->indg[2] = pot_ext->indg0[2] + iz - pot_ext->nring;
        if (pot_ext->indg[2] < 0 || pot_ext->indg[2] > dbgrid->getNX(2))
          return (1);
        iech = db_index_grid_to_sample(dbgrid, pot_ext->indg);

        /* Check that the external drift value is defined */

        drift = dbgrid->getLocVariable(ELoc::F,iech, 0);
        if (FFFF(drift)) return (1);
        pot_ext->data[ecr] = drift;
        ecr++;
      }
  return (0);
}

/****************************************************************************/
/*!
 **  Calculate the external drift contribution
 **
 ** \return  Error return code (target not within the grid or target on the
 ** \return  edge of the grid of external drift definition)
 **
 ** \param[in] target   Type of the target
 ** \param[in] x0       Coordinate along X
 ** \param[in] y0       Coordinate along Y
 ** \param[in] z0       Coordinate along Z
 ** \param[in] dbgrid   Output Db Grid structure
 ** \param[in] pot_ext  Pot_Ext structure
 **
 ** \param[out] extval  Value of the external drift
 ** \param[out] extgrd  Gradient components of the external drift
 **
 *****************************************************************************/
static int st_extdrift_eval(const char *target,
                            double x0,
                            double y0,
                            double z0,
                            DbGrid *dbgrid,
                            Pot_Ext *pot_ext,
                            double *extval,
                            VectorDouble& extgrd)
{
  VectorDouble coor(3);
  VectorDouble result(4);
  int error;

  error = 1;
  if (dbgrid == nullptr) return 1;
  coor[0] = x0;
  coor[1] = y0;
  coor[2] = z0;

  /* Find the location of the target within the external drift grid */

  if (point_to_grid(dbgrid, coor.data(), 0, pot_ext->indg0) < 0) goto label_end;

  /* Find the neighborhood around the target grid node */

  if (st_extdrift_neigh(dbgrid, pot_ext)) goto label_end;

  /* Perform the estimation */

  matrix_product_safe(1, pot_ext->nfull, 4, pot_ext->data, pot_ext->weight,
                      result.data());

  /* Retrieve the results */

  *extval = result[0];
  for (int idim = 0; idim < pot_ext->ndim; idim++)
    extgrd[idim] = result[1 + idim];

  error = 0;

  label_end:
  if (error && VERBOSE)
  {
    messerr("The External Drift cannot be estimated at %s point (%lf %lf %lf)",
            target, x0, y0, z0);
  }
  return (error);
}

/****************************************************************************/
/*!
 **  Establish the cokriging system
 **
 ** \return  Error return code
 **
 ** \param[in]  pot_env       Pot_Env structure
 ** \param[in]  pot_ext       Pot_Ext structure
 ** \param[in]  dbout         Target Db structure (only used for external drift)
 ** \param[in]  model         Model structure
 ** \param[in]  nugget_grd    Nugget effect for Gradients
 ** \param[in]  nugget_tgt    Nugget effect for Tangents
 **
 ** \param[out] lhs           Cokriging matrix
 **
 ** \remark   Organization of the cokriging system
 ** \remark
 ** \remark   |   A11 =        |  A12  =        | A13  =              | F1 =  |
 ** \remark   | <Gu(i),Gv(i')> | <Gu(i) ,(T,G)> | <Gu(i),Pl(j)-Pl(0)> | G(Fl) |
 ** \remark   -----------------------------------------------------------------
 ** \remark   |   A21  =       |  A22  =        | A23  =              | F2 =  |
 ** \remark   | <(T,G),Gv(i')  | <(T,G) ,(T,G)> | <(T,G),Pl(j)-Pl(0)> | (T,Fl)|
 ** \remark   -----------------------------------------------------------------
 ** \remark   |   A31  =       |  A32  =        | A33 =               | F3 =  |
 ** \remark   | <Pl(i)-Pl(0),G>| <Pli-Pl0,(T,g)>| <Pli-Pl0,Pl'i-Pl'0> |Fli-Fl0|
 ** \remark   -----------------------------------------------------------------
 ** \remark   | F1t            | F2t            | F3t                 |  0    |
 ** \remark
 ** \remark   The matrix A11 is subdivided as follows:
 ** \remark
 ** \remark         |  <Gx,Gx> |  <Gx,Gy>  |  <Gx,Gz> |
 ** \remark         |          |           |          |
 ** \remark   A11 = |  <Gy,Gx> |  <Gy,Gy>  |  ....    |
 ** \remark         |          |           |          |
 ** \remark         |  <Gz,Gx> |  ...      |  ..      |
 ** \remark
 ** \remark   each one of the 9 blocks has dimension = the number of gradients
 **
 *****************************************************************************/
static int st_build_lhs(Pot_Env *pot_env,
                        Pot_Ext *pot_ext,
                        DbGrid *dbout,
                        Model *model,
                        double nugget_grd,
                        double nugget_tgt,
                        double *lhs)
{
  int iext, nequa, ndim;
  double extval, extval1, extval2;

  ndim = pot_env->ndim;
  VectorDouble covGp(3, 0.);
  VectorDouble covGG(9, 0.);
  VectorDouble cov2Gp(3, 0.);
  VectorDouble cov2GG(9, 0.);
  VectorDouble center(3, 0.);
  VectorDouble extgrd(3, 0.);
  double covar = 0.;
  double covar1 = 0.;
  double covar2 = 0.;
  double covar3 = 0.;
  double covar4 = 0.;

  // Blank out the cokriging matrix

  nequa = pot_env->nequa;
  for (int i = 0; i < nequa * nequa; i++)
    lhs[i] = 0.;

  /******************************/
  /* PART RELATIVE TO GRADIENTS */
  /******************************/

  for (int ig = 0; ig < pot_env->ngrd; ig++)
  {
    for (int jg = 0; jg < ig; jg++)
    {
      st_cov(model, 1,
             GRD_COO(ig,0) - GRD_COO(jg, 0),
             GRD_COO(ig,1) - GRD_COO(jg, 1),
             GRD_COO(ig,2) - GRD_COO(jg, 2),
             covar, covGp, covGG);

      set_lhs(lhs, nequa, GRX(ig), GRX(jg), covGG[0]);
      set_lhs(lhs, nequa, GRX(ig), GRY(jg), covGG[1]);
      set_lhs(lhs, nequa, GRX(ig), GRZ(jg), covGG[2]);
      set_lhs(lhs, nequa, GRY(ig), GRX(jg), covGG[3]);
      set_lhs(lhs, nequa, GRY(ig), GRY(jg), covGG[4]);
      set_lhs(lhs, nequa, GRY(ig), GRZ(jg), covGG[5]);
      set_lhs(lhs, nequa, GRZ(ig), GRX(jg), covGG[6]);
      set_lhs(lhs, nequa, GRZ(ig), GRY(jg), covGG[7]);
      set_lhs(lhs, nequa, GRZ(ig), GRZ(jg), covGG[8]);
    }
    st_cov(model, 1, 0., 0., 0., covar, covGp, covGG);
    set_lhs(lhs, nequa, GRX(ig), GRX(ig), covGG[0] + nugget_grd);
    set_lhs(lhs, nequa, GRX(ig), GRY(ig), covGG[1]);
    set_lhs(lhs, nequa, GRX(ig), GRZ(ig), covGG[2]);
    set_lhs(lhs, nequa, GRY(ig), GRX(ig), covGG[3]);
    set_lhs(lhs, nequa, GRY(ig), GRY(ig), covGG[4] + nugget_grd);
    set_lhs(lhs, nequa, GRY(ig), GRZ(ig), covGG[5]);
    set_lhs(lhs, nequa, GRZ(ig), GRX(ig), covGG[6]);
    set_lhs(lhs, nequa, GRZ(ig), GRY(ig), covGG[7]);
    set_lhs(lhs, nequa, GRZ(ig), GRZ(ig), covGG[8] + nugget_grd);
  }

  /*****************************/
  /* PART RELATIVE TO TANGENTS */
  /*****************************/

  for (int it = 0; it < pot_env->ntgt; it++)
  {

    /* block tangents-gradients */

    for (int ig = 0; ig < pot_env->ngrd; ig++)
    {
      st_cov(model, 1,
             TGT_COO(it, 0) - GRD_COO(ig, 0),
             TGT_COO(it, 1) - GRD_COO(ig, 1),
             TGT_COO(it, 2) - GRD_COO(ig, 2),
             covar, covGp, covGG);

      set_lhs(lhs, nequa, TGT(it), GRX(ig),
              matrix_UV(ndim, TGT_VAL(it, 0), TGT_VAL(it, 1), TGT_VAL(it, 2),
                        covGG[0], covGG[1], covGG[2]));
      set_lhs(lhs, nequa, TGT(it), GRY(ig),
              matrix_UV(ndim, TGT_VAL(it, 0), TGT_VAL(it, 1), TGT_VAL(it, 2),
                        covGG[3], covGG[4], covGG[5]));
      set_lhs(lhs, nequa, TGT(it), GRZ(ig),
              matrix_UV(ndim, TGT_VAL(it, 0), TGT_VAL(it, 1), TGT_VAL(it, 2),
                        covGG[6], covGG[7], covGG[8]));
    }

    /* block diagonal tangents */

    for (int jt = 0; jt < it; jt++)
    {
      st_cov(model, 1,
             TGT_COO(it,0) - TGT_COO(jt, 0),
             TGT_COO(it,1) - TGT_COO(jt, 1),
             TGT_COO(it,2) - TGT_COO(jt, 2),
             covar, covGp, covGG);

      set_lhs(lhs, nequa, TGT(it), TGT(jt),
          matrix_UAV(ndim, covGG.data(),
                     TGT_VAL(it, 0), TGT_VAL(it, 1), TGT_VAL(it, 2),
                     TGT_VAL(jt, 0), TGT_VAL(jt, 1), TGT_VAL(jt, 2)));
    }
    st_cov(model, 1, 0., 0., 0., covar, covGp, covGG);
    set_lhs(lhs, nequa, TGT(it), TGT(it),
        matrix_UAV(ndim, covGG.data(),
                   TGT_VAL(it, 0), TGT_VAL(it, 1), TGT_VAL(it, 2),
                   TGT_VAL(it, 0), TGT_VAL(it, 1), TGT_VAL(it, 2))
        + nugget_tgt);
  }

  /***********************************/
  /* PART RELATIVE TO ISO-POTENTIALS */
  /***********************************/

  for (int ic1 = 0; ic1 < pot_env->nlayers; ic1++)
  {
    for (int j = 1; j < pot_env->nb_per_layer[ic1]; j++)
    {

      /* Interactions isopotentials - gradients */

      for (int ig = 0; ig < pot_env->ngrd; ig++)
      {
        st_cov(model, 1,
               GRD_COO(ig,0) - ISO_COO(ic1, 0, 0),
               GRD_COO(ig,1) - ISO_COO(ic1, 0, 1),
               GRD_COO(ig,2) - ISO_COO(ic1, 0, 2),
               covar, covGp, covGG);
        st_cov(model, 1,
               GRD_COO(ig,0) - ISO_COO(ic1, j, 0),
               GRD_COO(ig,1) - ISO_COO(ic1, j, 1),
               GRD_COO(ig,2) - ISO_COO(ic1, j, 2),
               covar, cov2Gp, cov2GG);
        set_lhs(lhs, nequa, ISC(ic1, j), GRX(ig), cov2Gp[0] - covGp[0]);
        set_lhs(lhs, nequa, ISC(ic1, j), GRY(ig), cov2Gp[1] - covGp[1]);
        set_lhs(lhs, nequa, ISC(ic1, j), GRZ(ig), cov2Gp[2] - covGp[2]);
      }

      /* Interactions increments-tangentes */

      for (int it = 0; it < pot_env->ntgt; it++)
      {
        st_cov(model, 1,
               TGT_COO(it,0) - ISO_COO(ic1, 0, 0),
               TGT_COO(it,1) - ISO_COO(ic1, 0, 1),
               TGT_COO(it,2) - ISO_COO(ic1, 0, 2),
               covar, covGp, covGG);
        st_cov(model, 1,
               TGT_COO(it,0) - ISO_COO(ic1, j, 0),
               TGT_COO(it,1) - ISO_COO(ic1, j, 1),
               TGT_COO(it,2) - ISO_COO(ic1, j, 2),
               covar, cov2Gp, cov2GG);
        set_lhs(lhs, nequa, ISC(ic1, j), TGT(it),
            matrix_UV(ndim, TGT_VAL(it, 0), TGT_VAL(it, 1), TGT_VAL(it, 2),
                      cov2Gp[0] - covGp[0], cov2Gp[1] - covGp[1],
                      cov2Gp[2] - covGp[2]));
      }

      /* Block diagonal for iso-potentials */

      for (int ic2 = 0; ic2 <= ic1; ic2++)
      {
        for (int j2 = 1; j2 < pot_env->nb_per_layer[ic2]; j2++)
        {
          st_cov(model, 0,
                 ISO_COO(ic2,j2,0) - ISO_COO(ic1, j, 0),
                 ISO_COO(ic2,j2,1) - ISO_COO(ic1, j, 1),
                 ISO_COO(ic2,j2,2) - ISO_COO(ic1, j, 2),
                 covar1, covGp, covGG);
          st_cov(model, 0,
                 ISO_COO(ic2,j2,0) - ISO_COO(ic1, 0, 0),
                 ISO_COO(ic2,j2,1) - ISO_COO(ic1, 0, 1),
                 ISO_COO(ic2,j2,2) - ISO_COO(ic1, 0, 2),
                 covar2, covGp, covGG);
          st_cov(model, 0,
                 ISO_COO(ic2,0,0) - ISO_COO(ic1, j, 0),
                 ISO_COO(ic2,0,1) - ISO_COO(ic1, j, 1),
                 ISO_COO(ic2,0,2) - ISO_COO(ic1, j, 2),
                 covar3, covGp, covGG);
          st_cov(model, 0,
                 ISO_COO(ic2,0,0) - ISO_COO(ic1, 0, 0),
                 ISO_COO(ic2,0,1) - ISO_COO(ic1, 0, 1),
                 ISO_COO(ic2,0,2) - ISO_COO(ic1, 0, 2),
                 covar4, covGp, covGG);
          set_lhs(lhs, nequa, ISC(ic1, j), ISC(ic2, j2),
                  covar1 - covar2 - covar3 + covar4);
        }
      }
    }
  }

  /****************************/
  /* Part linked to the drift */
  /****************************/

  /* Part relative to gradients */

  for (int ig = 0; ig < pot_env->ngrd; ig++)
  {
    if (pot_env->order >= 1)
    {

      /* Function x, y and z */

      set_lhs(lhs, nequa, DRF(0), GRX(ig), 1.);
      set_lhs(lhs, nequa, DRF(0), GRY(ig), 0.);
      set_lhs(lhs, nequa, DRF(0), GRZ(ig), 0.);

      set_lhs(lhs, nequa, DRF(1), GRX(ig), 0.);
      set_lhs(lhs, nequa, DRF(1), GRY(ig), 1.);
      set_lhs(lhs, nequa, DRF(1), GRZ(ig), 0.);

      set_lhs(lhs, nequa, DRF(2), GRX(ig), 0.);
      set_lhs(lhs, nequa, DRF(2), GRY(ig), 0.);
      set_lhs(lhs, nequa, DRF(2), GRZ(ig), 1.);
    }

    if (pot_env->order >= 2)
    {

      /* Functions x^2, y^2 et z^2 */

      set_lhs(lhs, nequa, DRF(3), GRX(ig), 2. * GRD_COO(ig, 0));
      set_lhs(lhs, nequa, DRF(3), GRY(ig), 0.);
      set_lhs(lhs, nequa, DRF(3), GRZ(ig), 0.);

      set_lhs(lhs, nequa, DRF(4), GRX(ig), 0.);
      set_lhs(lhs, nequa, DRF(4), GRY(ig), 2. * GRD_COO(ig, 1));
      set_lhs(lhs, nequa, DRF(4), GRZ(ig), 0.);

      set_lhs(lhs, nequa, DRF(5), GRX(ig), 0.);
      set_lhs(lhs, nequa, DRF(5), GRY(ig), 0.);
      set_lhs(lhs, nequa, DRF(5), GRZ(ig), 2. * GRD_COO(ig, 2));

      /* Functions xy, xz, et yz */

      set_lhs(lhs, nequa, DRF(6), GRX(ig), GRD_COO(ig, 1));
      set_lhs(lhs, nequa, DRF(6), GRY(ig), GRD_COO(ig, 0));
      set_lhs(lhs, nequa, DRF(6), GRZ(ig), 0.);

      set_lhs(lhs, nequa, DRF(7), GRX(ig), GRD_COO(ig, 2));
      set_lhs(lhs, nequa, DRF(7), GRY(ig), 0.);
      set_lhs(lhs, nequa, DRF(7), GRZ(ig), GRD_COO(ig, 0));

      set_lhs(lhs, nequa, DRF(8), GRX(ig), 0.);
      set_lhs(lhs, nequa, DRF(8), GRY(ig), GRD_COO(ig, 2));
      set_lhs(lhs, nequa, DRF(8), GRZ(ig), GRD_COO(ig, 1));
    }

    /* External drift(s) */

    for (iext = 0; iext < pot_env->next; iext++)
    {
      if (st_extdrift_eval("Gradient",
                           GRD_COO(ig, 0), GRD_COO(ig, 1), GRD_COO(ig, 2),
                           dbout, pot_ext, &extval, extgrd))
        return (1);
      set_lhs(lhs, nequa, EXT(iext), GRX(ig), extgrd[0]);
      set_lhs(lhs, nequa, EXT(iext), GRY(ig), extgrd[1]);
      set_lhs(lhs, nequa, EXT(iext), GRZ(ig), extgrd[2]);
    }
  }

  /* Part relative to tangents : Tx*f'x +Ty*f'y +Tz*f'z  */

  for (int it = 0; it < pot_env->ntgt; it++)
  {
    if (pot_env->order >= 1)
    {

      /* Derivates f = x, y, et z */

      set_lhs(lhs, nequa, DRF(0), TGT(it), TGT_VAL(it, 0));
      set_lhs(lhs, nequa, DRF(1), TGT(it), TGT_VAL(it, 1));
      set_lhs(lhs, nequa, DRF(2), TGT(it), TGT_VAL(it, 2));
    }

    if (pot_env->order >= 2)
    {

      /* Derivates f = x^2, y^2, et z^2 */

      set_lhs(lhs, nequa, DRF(3), TGT(it),
              2. * TGT_COO(it, 0) * TGT_VAL(it, 0));
      set_lhs(lhs, nequa, DRF(4), TGT(it),
              2. * TGT_COO(it, 1) * TGT_VAL(it, 1));
      set_lhs(lhs, nequa, DRF(5), TGT(it),
              2. * TGT_COO(it, 2) * TGT_VAL(it, 2));

      /* Derivates f = xy, xz, et yz */

      set_lhs(lhs, nequa, DRF(6), TGT(it), (TGT_COO(it,1) * TGT_VAL(it, 0) +
      TGT_COO(it,0) * TGT_VAL(it, 1)));
      set_lhs(lhs, nequa, DRF(7), TGT(it), (TGT_COO(it,2) * TGT_VAL(it, 0) +
      TGT_COO(it,0) * TGT_VAL(it, 2)));
      set_lhs(lhs, nequa, DRF(8), TGT(it), (TGT_COO(it,2) * TGT_VAL(it, 1) +
      TGT_COO(it,1) * TGT_VAL(it, 2)));
    }

    /* External drift(s) */

    for (iext = 0; iext < pot_env->next; iext++)
    {
      if (st_extdrift_eval("Tangent", TGT_COO(it, 0), TGT_COO(it, 1),
                           TGT_COO(it, 2), dbout, pot_ext, &extval, extgrd))
        return (1);
      set_lhs(lhs, nequa, EXT(iext), TGT(it),
          matrix_UV(ndim, TGT_VAL(it, 0), TGT_VAL(it, 1), TGT_VAL(it, 2),
                    extgrd[0], extgrd[1], extgrd[2]));
    }
  }

  /* Part relative to the iso-potentials */

  for (int ic1 = 0; ic1 < pot_env->nlayers; ic1++)
  {
    for (int j = 1; j < pot_env->nb_per_layer[ic1]; j++)
    {
      if (pot_env->order >= 1)
      {

        /* Functions x, y, z */

        set_lhs(lhs, nequa, DRF(0), ISC(ic1, j),
        ISO_COO(ic1,j,0) - ISO_COO(ic1, 0, 0));
        set_lhs(lhs, nequa, DRF(1), ISC(ic1, j),
        ISO_COO(ic1,j,1) - ISO_COO(ic1, 0, 1));
        set_lhs(lhs, nequa, DRF(2), ISC(ic1, j),
        ISO_COO(ic1,j,2) - ISO_COO(ic1, 0, 2));
      }

      /* Functions x^2, y^2, z^2 */

      if (pot_env->order >= 2)
      {
        set_lhs(lhs, nequa, DRF(3), ISC(ic1, j),
        ISO_COO(ic1,j,0) * ISO_COO(ic1, j, 0) -
        ISO_COO(ic1,0,0) * ISO_COO(ic1, 0, 0));
        set_lhs(lhs, nequa, DRF(4), ISC(ic1, j),
        ISO_COO(ic1,j,1) * ISO_COO(ic1, j, 1) -
        ISO_COO(ic1,0,1) * ISO_COO(ic1, 0, 1));
        set_lhs(lhs, nequa, DRF(5), ISC(ic1, j),
        ISO_COO(ic1,j,2) * ISO_COO(ic1, j, 2) -
        ISO_COO(ic1,0,2) * ISO_COO(ic1, 0, 2));

        /* Functions xy,xz, yz */

        set_lhs(lhs, nequa, DRF(6), ISC(ic1, j),
        ISO_COO(ic1,j,0) * ISO_COO(ic1, j, 1) -
        ISO_COO(ic1,0,0) * ISO_COO(ic1, 0, 1));
        set_lhs(lhs, nequa, DRF(7), ISC(ic1, j),
        ISO_COO(ic1,j,0) * ISO_COO(ic1, j, 2) -
        ISO_COO(ic1,0,0) * ISO_COO(ic1, 0, 2));
        set_lhs(lhs, nequa, DRF(8), ISC(ic1, j),
        ISO_COO(ic1,j,1) * ISO_COO(ic1, j, 2) -
        ISO_COO(ic1,0,1) * ISO_COO(ic1, 0, 2));
      }

      /* External drift(s) */

      for (iext = 0; iext < pot_env->next; iext++)
      {
        if (st_extdrift_eval("Iso-potential", ISO_COO(ic1, j, 0),
                             ISO_COO(ic1, j, 1), ISO_COO(ic1, j, 2), dbout,
                             pot_ext, &extval2, extgrd)) return (1);
        if (st_extdrift_eval("Iso-potential", ISO_COO(ic1, 0, 0),
                             ISO_COO(ic1, 0, 1), ISO_COO(ic1, 0, 2), dbout,
                             pot_ext, &extval1, extgrd)) return (1);
        set_lhs(lhs, nequa, EXT(iext), ISC(ic1, j), extval2 - extval1);
      }
    }
  }

  return (0);
}

/****************************************************************************/
/*!
 **  Establish the data vector
 **
 ** \param[in]  pot_env       Pot_Env structure
 **
 ** \param[out] zval          Data vector
 **
 *****************************************************************************/
static void st_fill_dual(Pot_Env *pot_env, double *zval)
{
  int nequa;

  // Initializations 

  nequa = pot_env->nequa;

  // Blank out the vector

  for (int i = 0; i < nequa; i++) zval[i] = 0.;

  for (int ig = 0; ig < pot_env->ngrd; ig++)
  {
    if (GRX(ig) >= 0) zval[GRX(ig)] = GRD_VAL(ig, 0);
    if (GRY(ig) >= 0) zval[GRY(ig)] = GRD_VAL(ig, 1);
    if (GRZ(ig) >= 0) zval[GRZ(ig)] = GRD_VAL(ig, 2);
  }
}

/****************************************************************************/
/*!
 **  Establish the simulation errors
 **
 ** \param[in]  pot_env       Pot_Env structure
 ** \param[in]  dbiso         Db containing the iso-values
 ** \param[in]  dbgrd         Gradient Db structure
 ** \param[in]  dbtgt         Tangent Db structure (optional)
 ** \param[in]  nbsimu        Number of simulations
 **
 ** \param[out] zval          Simulated errors
 **
 *****************************************************************************/
static void st_fill_dual_simulation(Pot_Env *pot_env,
                                    Db *dbiso,
                                    Db *dbgrd,
                                    Db *dbtgt,
                                    int nbsimu,
                                    double *zval)
{
  int nequa, shift, ndim;

  // Initializations 

  nequa = pot_env->nequa;
  ndim = dbgrd->getNDim();

  // Blank out the vector

  for (int i = 0; i < nequa * nbsimu; i++)
    zval[i] = 0.;

  // Loop on the simulations */

  shift = 0;
  for (int isimu = 0; isimu < nbsimu; isimu++)
  {

    // Load the gradient simulation errors

    for (int ig = 0; ig < pot_env->ngrd; ig++)
    {
      if (ndim >= 1)
        zval[shift + GRX(ig)] = dbgrd->getSimvar(ELoc::SIMU, IAD_GRD(ig),
                                                 isimu + 0 * nbsimu, 0, 0,
                                                 ndim * nbsimu, 1)
                                - GRD_VAL(ig, 0);
      if (ndim >= 2)
        zval[shift + GRY(ig)] = dbgrd->getSimvar(ELoc::SIMU, IAD_GRD(ig),
                                                 isimu + 1 * nbsimu, 0, 0,
                                                 ndim * nbsimu, 1)
                                - GRD_VAL(ig, 1);
      if (ndim >= 3)
        zval[shift + GRZ(ig)] = dbgrd->getSimvar(ELoc::SIMU, IAD_GRD(ig),
                                                 isimu + 2 * nbsimu, 0, 0,
                                                 ndim * nbsimu, 1)
                                - GRD_VAL(ig, 2);
    }

    // Load the tangent simulation errors 

    for (int it = 0; it < pot_env->ntgt; it++)
    {
      zval[shift + TGT(it)] = dbtgt->getSimvar(ELoc::SIMU, IAD_TGT(it), isimu,
                                               0, 0, nbsimu, 1);
    }

    // Load the iso-potential simulation errors

    for (int ic = 0; ic < pot_env->nlayers; ic++)
      for (int j = 1; j < pot_env->nb_per_layer[ic]; j++)
      {
        zval[shift + ISC(ic, j)] = dbiso->getSimvar(ELoc::SIMU, IAD_ISO(ic, j),
                                                    isimu, 0, 0, nbsimu, 1)
                                   - dbiso->getSimvar(ELoc::SIMU,
                                                      IAD_ISO(ic, 0), isimu, 0,
                                                      0, nbsimu, 1);
      }

    shift += nequa;
  }
}

/****************************************************************************/
/*!
 **  Blank out part the R.H.S. according to 'flag.part'
 **
 ** \param[in]  pot_env       Pot_Env structure
 **
 ** \param[in,out] rhs        Array for the R.H.S.
 **
 *****************************************************************************/
static void st_rhs_part(Pot_Env *pot_env,
                        double *rhs)
{
  int nequa = pot_env->nequa;
  int opt_part = pot_env->opt_part;
  int ideb = 0;
  int ifin = nequa;
  if (opt_part == 0) return;

  /* Dispatch */

  switch (opt_part)
  {
    case 1: /* Reveal Gradient */
      ideb = pot_env->start_grd;
      ifin = ideb + pot_env->size_grd;
      break;

    case 2: /* Reveal Tangent */
      ideb = pot_env->start_tgt;
      ifin = ideb + pot_env->size_tgt;
      break;

    case 3: /* Reveal Isovalues */
      ideb = pot_env->start_iso;
      ifin = ideb + pot_env->size_iso;
      break;

    case 4: /* Reveal internal drift */
      ideb = pot_env->start_drf;
      ifin = ideb + pot_env->size_drf;
      break;

    case 5: /* Reveal external drift */
      ideb = pot_env->start_ext;
      ifin = ideb + pot_env->size_ext;
      break;
  }

  /* Blank out the R.H.S. */

  for (int i = 0; i < nequa; i++)
  {
    if (i >= ideb && i < ifin) continue;
    set_rhs(rhs, nequa, i, 0, 0.);
    if (pot_env->flag_grad)
      for (int igrad = 1; igrad < 4; igrad++)
        set_rhs(rhs, nequa, i, igrad, 0.);
  }
  return;
}

/****************************************************************************/
/*!
 **  Calculate the estimation and gradient components at one target location
 **
 ** \param[in]  pot_env       Pot_Env structure
 ** \param[in]  pot_ext       Pot_Ext structure
 ** \param[in]  flag_grad     True if the gradients must also be calculated
 ** \param[in]  dbgrid        Output Grid Db structure (for External Drift)
 ** \param[in]  model         Model structure
 ** \param[in]  coor          Coordinates of the target
 **
 ** \param[out] rhs           R.H.S. vector (Dimension: nequa * 4)
 **
 *****************************************************************************/
static void st_build_rhs(Pot_Env *pot_env,
                         Pot_Ext *pot_ext,
                         bool flag_grad,
                         DbGrid *dbgrid,
                         Model *model,
                         VectorDouble& coor,
                         double *rhs)
{
  int nequa, nsol, ndim;
  double extval;

  /* Initializations */

  double covar = 0.;
  double covar1 = 0.;
  ndim = pot_env->ndim;
  nequa = pot_env->nequa;
  nsol = (flag_grad) ? 1 + pot_env->ndim : 1;

  VectorDouble covGp(3, 0);
  VectorDouble covGG(9, 0.);
  VectorDouble cov1Gp(3, 0.);
  VectorDouble cov1GG(9, 0.);
  VectorDouble center(3, 0.);
  VectorDouble ccor(3, 0.);
  VectorDouble extgrd(3, 0.);
  for (int i = 0; i < nequa * nsol; i++)
    rhs[i] = 0.;

  /*******************/
  /* Covariance part */
  /*******************/

  /* Part relative to gradients */

  for (int ig = 0; ig < pot_env->ngrd; ig++)
  {
    st_cov(model, flag_grad,
           GRD_COO(ig,0) - coor[0],
           GRD_COO(ig,1) - coor[1],
           GRD_COO(ig,2) - coor[2],
           covar, covGp, covGG);
    set_rhs(rhs, nequa, GRX(ig), 0, covGp[0]);
    set_rhs(rhs, nequa, GRY(ig), 0, covGp[1]);
    set_rhs(rhs, nequa, GRZ(ig), 0, covGp[2]);
    if (flag_grad)
    {
      set_rhs(rhs, nequa, GRX(ig), 1, covGG[0]);
      set_rhs(rhs, nequa, GRY(ig), 1, covGG[1]);
      set_rhs(rhs, nequa, GRZ(ig), 1, covGG[2]);
      set_rhs(rhs, nequa, GRX(ig), 2, covGG[3]);
      set_rhs(rhs, nequa, GRY(ig), 2, covGG[4]);
      set_rhs(rhs, nequa, GRZ(ig), 2, covGG[5]);
      set_rhs(rhs, nequa, GRX(ig), 3, covGG[6]);
      set_rhs(rhs, nequa, GRY(ig), 3, covGG[7]);
      set_rhs(rhs, nequa, GRZ(ig), 3, covGG[8]);
    }
  }

  /* Part relative to tangents */

  for (int it = 0; it < pot_env->ntgt; it++)
  {
    st_cov(model, flag_grad,
           TGT_COO(it,0) - coor[0],
           TGT_COO(it,1) - coor[1],
           TGT_COO(it,2) - coor[2],
           covar, covGp, covGG);
    set_rhs(rhs, nequa, TGT(it), 0,
        matrix_UV(ndim, TGT_VAL(it, 0), TGT_VAL(it, 1), TGT_VAL(it, 2),
                  covGp[0], covGp[1], covGp[2]));
    if (flag_grad)
    {
      set_rhs(rhs, nequa, TGT(it), 1,
          matrix_UV(ndim, TGT_VAL(it, 0), TGT_VAL(it, 1), TGT_VAL(it, 2),
                    covGG[0], covGG[1], covGG[2]));
      set_rhs(rhs, nequa, TGT(it), 2,
          matrix_UV(ndim, TGT_VAL(it, 0), TGT_VAL(it, 1), TGT_VAL(it, 2),
                    covGG[3], covGG[4], covGG[5]));
      set_rhs(rhs, nequa, TGT(it), 3,
          matrix_UV(ndim, TGT_VAL(it, 0), TGT_VAL(it, 1), TGT_VAL(it, 2),
                    covGG[6], covGG[7], covGG[8]));
    }
  }

  /* Part relative to iso-potentials */

  for (int ic = 0; ic < pot_env->nlayers; ic++)
  {
    for (int j = 1; j < pot_env->nb_per_layer[ic]; j++)
    {
      st_cov(model, flag_grad,
             ISO_COO(ic,j,0) - coor[0],
             ISO_COO(ic,j,1) - coor[1],
             ISO_COO(ic,j,2) - coor[2],
             covar1, cov1Gp, cov1GG);
      st_cov(model, flag_grad,
             ISO_COO(ic,0,0) - coor[0],
             ISO_COO(ic,0,1) - coor[1],
             ISO_COO(ic,0,2) - coor[2],
             covar, covGp, covGG);
      set_rhs(rhs, nequa, ISC(ic, j), 0, covar1 - covar);
      if (flag_grad)
      {
        set_rhs(rhs, nequa, ISC(ic, j), 1, -(cov1Gp[0] - covGp[0]));
        set_rhs(rhs, nequa, ISC(ic, j), 2, -(cov1Gp[1] - covGp[1]));
        set_rhs(rhs, nequa, ISC(ic, j), 3, -(cov1Gp[2] - covGp[2]));
      }
    }
  }

  /****************************/
  /* Part linked to the drift */
  /****************************/

  for (int i = 0; i < 3; i++)
    ccor[i] = coor[i] - center[i];

  if (pot_env->order >= 1)
  {
    set_rhs(rhs, nequa, DRF(0), 0, ccor[0]);
    set_rhs(rhs, nequa, DRF(1), 0, ccor[1]);
    set_rhs(rhs, nequa, DRF(2), 0, ccor[2]);
    if (flag_grad)
    {
      set_rhs(rhs, nequa, DRF(0), 1, 1.);
      set_rhs(rhs, nequa, DRF(1), 2, 1.);
      set_rhs(rhs, nequa, DRF(2), 3, 1.);
    }
  }

  if (pot_env->order >= 2)
  {
    set_rhs(rhs, nequa, DRF(3), 0, ccor[0] * ccor[0]);
    set_rhs(rhs, nequa, DRF(4), 0, ccor[1] * ccor[1]);
    set_rhs(rhs, nequa, DRF(5), 0, ccor[2] * ccor[2]);
    set_rhs(rhs, nequa, DRF(6), 0, ccor[0] * ccor[1]);
    set_rhs(rhs, nequa, DRF(7), 0, ccor[0] * ccor[2]);
    set_rhs(rhs, nequa, DRF(8), 0, ccor[1] * ccor[2]);
    if (flag_grad)
    {
      set_rhs(rhs, nequa, DRF(3), 1, ccor[0] * 2.);
      set_rhs(rhs, nequa, DRF(4), 2, ccor[1] * 2.);
      set_rhs(rhs, nequa, DRF(5), 3, ccor[2] * 2.);
      set_rhs(rhs, nequa, DRF(6), 1, ccor[1]);
      set_rhs(rhs, nequa, DRF(6), 2, ccor[0]);
      set_rhs(rhs, nequa, DRF(7), 1, ccor[2]);
      set_rhs(rhs, nequa, DRF(7), 3, ccor[0]);
      set_rhs(rhs, nequa, DRF(8), 2, ccor[2]);
      set_rhs(rhs, nequa, DRF(8), 3, ccor[1]);
    }
  }

  for (int iext = 0; iext < pot_env->next; iext++)
  {
    if (st_extdrift_eval("Target", coor[0], coor[1], coor[2], dbgrid, pot_ext,
                         &extval, extgrd)) return;
    set_rhs(rhs, nequa, EXT(iext), 0, extval);
    if (flag_grad)
    {
      set_rhs(rhs, nequa, EXT(iext), 1, extgrd[0]);
      set_rhs(rhs, nequa, EXT(iext), 2, extgrd[1]);
      set_rhs(rhs, nequa, EXT(iext), 3, extgrd[2]);
    }
  }

  // Blank out the R.H.S. according to masking option 

  st_rhs_part(pot_env, rhs);

  // Printout (optional) 

  if (OptDbg::query(EDbg::KRIGING))
    krige_rhs_print(nsol, 0, nequa, nequa, NULL, rhs);

  return;
}

/****************************************************************************/
/*!
 **  Calculate the estimation and gradient components at one target
 **
 ** \param[in]  pot_env       Pot_Env structure
 ** \param[in]  pot_ext       Pot_Ext structure
 ** \param[in]  flag_grad     True if the gradients must also be calculated
 ** \param[in]  dbiso         Iso-potential Db structure (not used)
 ** \param[in]  dbgrd         Gradient Db structure (not used)
 ** \param[in]  dbtgt         Tangent Db structure (not used)
 ** \param[in]  dbgrid        Output Db structure (for Ext Drift)
 ** \param[in]  model         Model structure
 ** \param[in]  zdual         Dual vector (Dimension: nequa)
 ** \param[in]  rhs           R.H.S. vector (Dimension: nequa * 4)
 ** \param[in]  db_target     Db corresponding to the target
 ** \param[in]  iech0         Rank of the target sample
 **
 ** \param[out] result        Array of results (Dimension: nsol)
 **
 *****************************************************************************/
static void st_calc_point(Pot_Env *pot_env,
                          Pot_Ext *pot_ext,
                          bool flag_grad,
                          Db *dbiso,
                          Db *dbgrd,
                          Db *dbtgt,
                          DbGrid *dbgrid,
                          Model *model,
                          double *zdual,
                          double *rhs,
                          Db *db_target,
                          int iech0,
                          VectorDouble& result)
{
  DECLARE_UNUSED(dbiso, dbgrd, dbtgt);
  int nsol;
  VectorDouble coor(3,0.);

  /* Initializations */

  nsol = (flag_grad) ? 1 + pot_env->ndim : 1;

  /* Load the coordinates */

  for (int idim = 0; idim < pot_env->ndim; idim++)
    coor[idim] = db_target->getCoordinate(iech0, idim);

  /* Optional printout */

  if (OptDbg::query(EDbg::KRIGING) || OptDbg::query(EDbg::NBGH))
  {
    mestitle(1, "Target location");
    db_sample_print(db_target, iech0, 1, 0, 0);
  }

  /* Establish the R.H.S */

  st_build_rhs(pot_env, pot_ext, flag_grad, dbgrid, model, coor, rhs);

  /* Perform the estimation */

  for (int i = 0; i < nsol; i++) result[i] = TEST;
  matrix_product_safe(1, pot_env->nequa, nsol, zdual, rhs, result.data());

  // Printout (optional) 

  if (OptDbg::query(EDbg::KRIGING))
  {
    print_matrix("Results", 0, 1, 1, nsol, NULL, result.data());
    message("\n");
  }

  return;
}

/****************************************************************************/
/*!
 **  Translate potential into layer value or center it
 **
 ** \param[in]  pot_env       Pot_Env structure
 ** \param[in]  isimu         Rank of the simulation
 ** \param[in]  potval        Array of potential values at different layers
 ** \param[in]  result        Resulting value (in potential scale)
 **                           On output, Resulting value in layer scale
 **
 ** \remarks The potential values at iso-potential samples are assumed
 ** \remarks to be ordered
 ** \remarks It is assumed that the potential has already been centered
 ** \remarks Therefore the 'potval' values must also be centered (locally)
 **
 *****************************************************************************/
static void st_potential_to_layer(Pot_Env *pot_env,
                                  int isimu,
                                  double *potval,
                                  VectorDouble& result)
{
  double minval = -1.e30;
  double potref = POTVAL(isimu, 0);

  int ilayer = -1;
  for (int i = 0; i < pot_env->nlayers && ilayer < 0; i++)
  {
    if (result[0] > minval && result[0] <= (POTVAL(isimu,i) - potref))
      ilayer = i;
    minval = (POTVAL(isimu,i) - potref);
  }
  result[0] = ilayer + 1;
}

/****************************************************************************/
/*!
 **  Calculate the estimation and/or gradient components at target samples
 **
 ** \param[in]  pot_env       Pot_Env structure
 ** \param[in]  pot_ext       Pot_Ext structure
 ** \param[in]  flag_grad     True if the gradients must also be calculated
 ** \param[in]  dbiso         Iso-potential Db structure
 ** \param[in]  dbgrd         Gradient Db structure
 ** \param[in]  dbtgt         Tangent Db structure (optional)
 ** \param[in]  dbout         Output Db structure
 ** \param[in]  model         Model structure
 ** \param[in]  refpot        Potential at Reference point
 ** \param[in]  zdual         Dual vector (Dimension: nequa)
 ** \param[in]  rhs           R.H.S. vector (Dimension: nequa * 4)
 ** \param[in]  potval        Potential values at iso-potential samples
 **
 *****************************************************************************/
static void st_estimate_result(Pot_Env *pot_env,
                               Pot_Ext *pot_ext,
                               bool flag_grad,
                               Db *dbiso,
                               Db *dbgrd,
                               Db *dbtgt,
                               DbGrid *dbout,
                               Model *model,
                               double refpot,
                               double *zdual,
                               double *rhs,
                               double *potval)
{
  VectorDouble result(4);

  for (int iech = 0; iech < dbout->getSampleNumber(); iech++)
  {
    mes_process("Potential Estimation on Grid", dbout->getSampleNumber(),iech);
    OptDbg::setCurrentIndex(iech);
    if (!dbout->isActive(iech)) continue;

    // Perform the estimation

    st_calc_point(pot_env, pot_ext, flag_grad, dbiso, dbgrd, dbtgt,
                  dbout, model, zdual, rhs, dbout, iech, result);

    // Center to the reference potential

    result[0] -= refpot;

    // Printout (optional) 

    if (OptDbg::query(EDbg::KRIGING))
      message("Centered estimation = %lf\n", result[0]);

    // Translate from potential into layer

    if (pot_env->flag_trans)
      st_potential_to_layer(pot_env, 0, potval, result);

    // Store the results

    dbout->setLocVariable(ELoc::Z,iech, 0, result[0]);
    if (flag_grad)
      for (int idim = 0; idim < pot_env->ndim; idim++)
        dbout->setLocVariable(ELoc::G,iech, idim, result[idim + 1]);
  }
  OptDbg::setCurrentIndex(-1);
  return;
}

static void st_estimate_data(Pot_Env *pot_env,
                             Pot_Ext *pot_ext,
                             Db *dbiso,
                             Db *dbgrd,
                             Db *dbtgt,
                             DbGrid *dbout,
                             Model *model,
                             double refpot,
                             double *zdual,
                             double *rhs,
                             Db *db_target,
                             VectorInt& uid_pot,
                             VectorInt& uid_grad)
{
  VectorDouble result(4);
  if (db_target == nullptr) return;

  for (int iech = 0; iech < db_target->getSampleNumber(); iech++)
  {
    if (! db_target->isActive(iech)) continue;

    // Perform the estimation

    st_calc_point(pot_env, pot_ext, true,
                  dbiso, dbgrd, dbtgt, dbout, model,
                  zdual, rhs, db_target, iech, result);

    // Center to the reference potential

    result[0] -= refpot;

    // Store the results

    if (! uid_pot.empty())
    {
      db_target->setArray(iech, uid_pot[0], result[0]);
      db_target->setLocatorsByUID(uid_pot, ELoc::Z);
    }
    if (! uid_grad.empty())
    {
      for (int idim = 0; idim < pot_env->ndim; idim++)
        db_target->setArray(iech, uid_grad[idim], result[idim + 1]);
      db_target->setLocatorsByUID(uid_grad, ELoc::G);
    }
  }
  return;
}

/****************************************************************************/
/*!
 **  Calculate the cross-validation at the iso-potential samples
 **
 ** \param[in]  pot_env       Pot_Env structure
 ** \param[in]  pot_ext       Pot_Ext structure
 ** \param[in]  dbiso         Iso-potential Db structure (not used)
 ** \param[in]  dbgrd         Gradient Db structure  (not used)
 ** \param[in]  dbtgt         Tangent Db structure  (not used)
 ** \param[in]  model         Model structure
 ** \param[in]  ic0           Rank of the isoline
 ** \param[in]  j0            Rank of the sample within this isoline
 ** \param[in]  zval          Data vector
 ** \param[in]  lhs_orig      Copy of the initial LHS (non inverted)
 ** \param[in]  lhs_aux       Working array for storing the new LHS
 ** \param[in]  rhs           Right-hand side
 ** \param[in]  zdual         Dual vector (Dimension: nequa)
 **
 ** \param[out] dist_euc      Error converted into Euclidean distance
 ** \param[out] dist_geo      Error converted into along surface distance
 **
 ** \remark We assume that the new data set (with one sample OFF) is still
 ** \remark contained in 'zval' as the Gradient information (coming first in
 ** \remark this vector) is never excluded.
 **
 *****************************************************************************/
static void st_dist_convert(Pot_Env *pot_env,
                            Pot_Ext *pot_ext,
                            Db *dbiso,
                            Db *dbgrd,
                            Db *dbtgt,
                            Model *model,
                            int ic0,
                            int j0,
                            double *zval,
                            double *lhs_orig,
                            double *lhs_aux,
                            double *rhs,
                            double *zdual,
                            double *dist_euc,
                            double *dist_geo)
{
  DECLARE_UNUSED(dbiso, dbgrd, dbtgt);
  double potval, delta;
  int nsol, nequa, neqm1, icol0;
  VectorDouble result(4);
  static int niter_max = 50;
  static double eps = 1.e-3;

  VectorDouble coor(3, 0.);
  VectorDouble coor0(3, 0.);
  int ndim = pot_env->ndim;
  nequa = pot_env->nequa;
  neqm1 = nequa - 1;
  nsol = 1 + ndim;
  icol0 = ISC(ic0, j0);
  VectorDouble deuc(ndim, 0.);
  VectorDouble dgeo(ndim, 0.);

  /* Update the L.H.S. by dropping the current data point */

  matrix_manage(nequa, nequa, -1, -1, &icol0, &icol0, lhs_orig, lhs_aux);

  /* Invert the new LHS */

  if (matrix_invert(lhs_aux, neqm1, 0)) return;

  /* Calculate the dual system */

  matrix_product_safe(neqm1, neqm1, 1, lhs_aux, zval, zdual);

  /* Evaluate the reference point */

  for (int idim = 0; idim < pot_env->ndim; idim++)
    coor0[idim] = ISO_COO(ic0, 0, idim);
  st_build_rhs(pot_env, pot_ext, 0, nullptr, model, coor0, rhs);
  matrix_manage(nequa, 1, -1, 0, &icol0, NULL, rhs, rhs);
  matrix_product_safe(1, neqm1, 1, zdual, rhs, &potval);

  /* Evaluate the target point */

  for (int idim = 0; idim < pot_env->ndim; idim++)
    coor0[idim] = coor[idim] = ISO_COO(ic0, j0, idim);
  st_build_rhs(pot_env, pot_ext, 1, nullptr, model, coor0, rhs);
  matrix_manage(nequa, nsol, -1, 0, &icol0, NULL, rhs, rhs);
  matrix_product_safe(1, neqm1, nsol, zdual, rhs, result.data());
  if (OptDbg::query(EDbg::CONVERGE))
  {
    message("Sample:%2d/%2d Iter:%2d Potential:%lf", j0 + 1, ic0 + 1, 0,
            result[0]);
    for (int idim = 0; idim < pot_env->ndim; idim++)
      message(" %lf", coor0[idim]);
    message("\n");
  }

  /* Move the target and estimate again */

  for (int iter = 0; iter < niter_max; iter++)
  {
    if (ABS(result[0]) < eps) break;
    for (int idim = 0; idim < pot_env->ndim; idim++)
    {
      if (ABS(result[1+idim]) < eps) continue;
      delta = 0.1 * result[0] / result[1 + idim];
      coor[idim] -= delta;
      dgeo[idim] += delta * delta;
    }
    st_build_rhs(pot_env, pot_ext, 1, nullptr, model, coor, rhs);
    matrix_manage(nequa, nsol, -1, 0, &icol0, NULL, rhs, rhs);
    matrix_product_safe(1, neqm1, nsol, zdual, rhs, result.data());
    if (OptDbg::query(EDbg::CONVERGE))
    {
      message("Sample:%2d/%2d Iter:%2d Potential:%lf", j0 + 1, ic0 + 1, iter,
              result[0]);
      for (int idim = 0; idim < pot_env->ndim; idim++)
        message(" %lf", coor[idim]);
      message("\n");
    }
  }

  /* Determine the euclidean distance */

  for (int idim = 0; idim < pot_env->ndim; idim++)
  {
    delta = coor[idim] - coor0[idim];
    deuc[idim] = delta * delta;
  }

  /* Find both distances */

  (*dist_euc) = (*dist_geo) = 0.;
  for (int idim = 0; idim < pot_env->ndim; idim++)
  {
    (*dist_euc) += deuc[idim];
    (*dist_geo) += dgeo[idim];
  }
  (*dist_euc) = sqrt(*dist_euc);
  (*dist_geo) = sqrt(*dist_geo);
  return;
}

/****************************************************************************/
/*!
 **  Calculate the cross-validation at the iso-potential samples
 **
 ** \param[in]  pot_env        Pot_Env structure
 ** \param[in]  pot_ext        Pot_Ext structure
 ** \param[in]  dbiso          Iso-potential Db structure
 ** \param[in]  dbgrd          Gradient Db structure
 ** \param[in]  dbtgt          Tangent Db structure (optional)
 ** \param[in]  model          Model structure
 ** \param[in]  lhs            Inverted Kriging system
 ** \param[in]  flag_dist_conv Flag for converting into distance
 ** \param[in]  zval           Data vector
 ** \param[in]  lhs_orig       Copy of the Initial LHS
 ** \param[in]  lhs_aux        Working array for LHS
 ** \param[in]  rhs            Right-hand side
 ** \param[in]  zdual          Dual vector (Dimension: nequa)
 **
 ** \remarks Arguments from 'zval' are only used to convert into distance
 **
 *****************************************************************************/
static void st_xvalid_potential(Pot_Env *pot_env,
                                Pot_Ext *pot_ext,
                                Db *dbiso,
                                Db *dbgrd,
                                Db *dbtgt,
                                Model *model,
                                double *lhs,
                                int flag_dist_conv,
                                double *zval,
                                double *lhs_orig,
                                double *lhs_aux,
                                double *rhs,
                                double *zdual)
{
  double variance, value, stats[4][2], stdev, dist_euc, dist_geo;
  int nequa, icol0, iech0, number, nitem;
  VectorDouble result(4);

  /* Initializations */

  nequa = pot_env->nequa;
  nitem = (flag_dist_conv) ? 4 :
                             2;
  for (int i = 0; i < nitem; i++)
    for (int j = 0; j < 2; j++)
      stats[i][j] = 0.;

  /* Loop on the Iso-potential samples */

  for (int ic = 0; ic < pot_env->nlayers; ic++)
  {
    number = 0;
    for (int j = 1; j < pot_env->nb_per_layer[ic]; j++)
    {
      iech0 = IAD_ISO(ic, j);
      mes_process("Potential Estimation on Iso-Potential %d of %d", j+1,ic+1);
      OptDbg::setCurrentIndex(iech0);

      // Get the variance and the weights from the inverted L.H.S.

      icol0 = ISC(ic, j);
      variance = 1. / get_lhs(lhs, nequa, icol0, icol0);
      stdev = sqrt(variance);
      dist_geo = dist_euc = 0.;

      // Perform the estimation

      value = 0.;
      for (int ig = 0; ig < pot_env->ngrd; ig++)
      {
        value += get_lhs(lhs, nequa, icol0, GRX(ig)) * GRD_VAL(ig, 0);
        value += get_lhs(lhs, nequa, icol0, GRY(ig)) * GRD_VAL(ig, 1);
        value += get_lhs(lhs, nequa, icol0, GRZ(ig)) * GRD_VAL(ig, 2);
      }
      result[0] = -value * variance;

      // Finding the closest distance to the Isoline

      if (flag_dist_conv)
        st_dist_convert(pot_env, pot_ext, dbiso, dbgrd, dbtgt, model, ic, j,
                        zval, lhs_orig, lhs_aux, rhs, zdual, &dist_geo,
                        &dist_euc);

      // Debugging option

      if (OptDbg::query(EDbg::RESULTS))
      {
        message("Sample %d/%d (%d): Error=%lf - Variance=%lf", j + 1, ic + 1,
                iech0 + 1, result[0], variance);
        if (flag_dist_conv)
          message(" - D-Geo=%lf - D-Surf=%lf", dist_euc, dist_geo);
        message("\n");
      }

      // Storing the results 

      dbiso->setLocVariable(ELoc::Z,iech0, 0, result[0]);
      dbiso->setLocVariable(ELoc::Z,iech0, 1, variance);
      if (flag_dist_conv)
      {
        dbiso->setLocVariable(ELoc::Z,iech0, 2, dist_euc);
        dbiso->setLocVariable(ELoc::Z,iech0, 3, dist_geo);
      }

      // Update statistics */

      value = result[0];
      stats[0][0] += value;
      stats[0][1] += value * value;
      value = result[0] / stdev;
      stats[1][0] += value;
      stats[1][1] += value * value;
      if (flag_dist_conv)
      {
        value = dist_geo;
        stats[2][0] += value;
        stats[2][1] += value * value;
        value = dist_euc;
        stats[3][0] += value;
        stats[3][1] += value * value;
      }
      number++;
    }

    // Print the global statistics (optinal) 

    if (VERBOSE && number > 0)
    {
      for (int i = 0; i < nitem; i++)
      {
        for (int j = 0; j < 2; j++)
          stats[i][j] /= number;
        stats[i][1] -= stats[i][0] * stats[i][0];
        stats[i][1] = (stats[i][1] > 0) ? sqrt(stats[i][1]) :
                                          0.;
      }
      message("\nIso-Potential #%d\n", ic + 1);
      message("Cross-validation Error: Mean=%lf St. Dev.=%lf\n", stats[0][0],
              stats[0][1]);
      message("Standardized Error    : Mean=%lf St. Dev.=%lf\n", stats[1][0],
              stats[1][1]);
      if (flag_dist_conv)
      {
        message("Euclidean Distance    : Mean=%lf St. Dev.=%lf\n", stats[2][0],
                stats[2][1]);
        message("Geodetic Distance     : Mean=%lf St. Dev.=%lf\n", stats[3][0],
                stats[3][1]);
      }
    }
  }
  OptDbg::setCurrentIndex(-1);
  return;
}

/****************************************************************************/
/*!
 **  Amortize the conditional simulations
 **
 ** \param[in]  dbout         Output Db structure
 ** \param[in]  iech          Rank of the sample
 ** \param[in]  dist_tempere  Distance for tempering simulations (or TEST)
 ** \param[in]  reskrige      Kriging result
 ** \param[in]  result        Conditional Simulation result
 **                           On output, Conditional Simulation tempered result
 **
 ** \remarks This function does nothing if 'dist_tempere' is undefined
 **
 *****************************************************************************/
static void st_tempere(DbGrid *dbout,
                       int iech,
                       double dist_tempere,
                       double reskrige,
                       VectorDouble& result)
{
  double simerr, amortval, kdist;
  static int test = 0;

  simerr = result[0] - reskrige;
  kdist = dbout->getLocVariable(ELoc::Z,iech, 0);

  switch (test)
  {
    case 0: /* Simulation amortie */
      amortval = MIN(1., exp(-kdist / dist_tempere));
      result[0] = reskrige + simerr * amortval;
      break;

    case 1: /* Distance normee */
      result[0] = kdist / dist_tempere;
      break;

    case 2: /* Simulation Conditionnelle */
      break;

    case 3: /* Simulation non-conditionnelle */
      result[0] = simerr;
      break;

    case 4: /* Krigeage */
      result[0] = reskrige;
      break;
  }
}

/****************************************************************************/
/*!
 **  Calculate the conditional simulation at target samples
 **
 ** \param[in]  pot_env       Pot_Env structure
 ** \param[in]  pot_ext       Pot_Ext structure
 ** \param[in]  dist_tempere  Distance for tempering simulations (or TEST)
 ** \param[in]  flag_trans    True if the estimation result must be translated
 **                           into layer number
 ** \param[in]  nbsimu        Number of simulation
 ** \param[in]  dbiso         Iso-potential Db structure
 ** \param[in]  dbgrd         Gradient Db structure
 ** \param[in]  dbtgt         Tangent Db structure (optional)
 ** \param[in]  dbout         Output Db structure
 ** \param[in]  model         Model structure
 ** \param[in]  refpot        Potential at Reference point
 ** \param[in]  potsim        Potential simulated values at iso-potential samples
 ** \param[in]  zdual         Dual estimated vector (Dimension: nequa)
 ** \param[in]  zduals        Dual simulated vector (Dimension: nequa * nbsimu)
 ** \param[in]  rhs           R.H.S. vector (Dimension: nequa * 4)
 **
 *****************************************************************************/
static void st_simcond(Pot_Env *pot_env,
                       Pot_Ext *pot_ext,
                       double dist_tempere,
                       bool flag_trans,
                       int nbsimu,
                       Db *dbiso,
                       Db *dbgrd,
                       Db *dbtgt,
                       DbGrid *dbout,
                       Model *model,
                       double refpot,
                       double *potsim,
                       double *zdual,
                       double *zduals,
                       double *rhs)
{
  VectorDouble resest(4), result(4);
  int nequa, ndim;

  nequa = pot_env->nequa;
  ndim = dbgrd->getNDim();
  for (int iech = 0; iech < dbout->getSampleNumber(); iech++)
  {
    mes_process("Potential Simulation on Grid", dbout->getSampleNumber(),iech);
    OptDbg::setCurrentIndex(iech);
    if (!dbout->isActive(iech)) continue;

    if (!FFFF(dist_tempere))
    {

      // Perform the estimation

      st_calc_point(pot_env, pot_ext, 1, dbiso, dbgrd, dbtgt, dbout, model,
                    zdual, rhs, dbout, iech, resest);

      // Center to the reference potential

      resest[0] -= refpot;
    }

    for (int isimu = 0; isimu < nbsimu; isimu++)
    {

      // Perform the estimation of the simulated error

      st_calc_point(pot_env, pot_ext, 0, dbiso, dbgrd, dbtgt, dbout, model,
                    &ZDUALS(isimu, 0), rhs, dbout, iech, result);

      // Convert into simulation error

      result[0] = (dbout->getSimvar(ELoc::SIMU, iech, isimu, 0, 0, nbsimu, 1)
          - result[0]);
      for (int idim = 0; idim < ndim; idim++)
        result[1 + idim] = (dbgrd->getSimvar(ELoc::SIMU, iech,
                                             isimu + idim * nbsimu, 0, 0,
                                             ndim * nbsimu, 1)
                            - result[1 + idim]);

      // Center to the reference potential

      result[0] -= refpot;

      // Amortize the variance for conditional simulation

      if (!FFFF(dist_tempere))
        st_tempere(dbout, iech, dist_tempere, resest[0], result);

      // Translate from potential into layer

      if (flag_trans)
        st_potential_to_layer(pot_env, isimu, potsim, result);

      // Store the results

      dbout->setSimvar(ELoc::SIMU, iech, isimu, 0, 0, nbsimu, 1, result[0]);
    }
  }
  OptDbg::setCurrentIndex(-1);
  return;
}

/****************************************************************************/
/*!
 **  Print the estimation at a target sample
 **
 ** \param[in]  pot_env    Pot_env structure
 ** \param[in]  isimu      Rank of the simulation (or -1)
 ** \param[in]  result     Array of results
 ** \param[in]  tgtval     Value of the tangent (or TEST)
 **
 *****************************************************************************/
static void st_print_result(Pot_Env *pot_env,
                            int isimu,
                            double *result,
                            double tgtval)
{
  if (isimu >= 0) message("Simulation %2d - ", isimu + 1);

  message(" - Pot* =%10.5lf", result[0]);

  message(" - Grad* =");
  for (int idim = 0; idim < pot_env->ndim; idim++)
    message(" %10.5lf", result[1 + idim]);

  if (!FFFF(tgtval)) message(" - Tangent= %10.5lf", tgtval);

  message("\n");
}

/****************************************************************************/
/*!
 **  Calculate the estimation and gradient components at data information
 **
 ** \param[in]  pot_env       Pot_Env structure
 ** \param[in]  pot_ext       Pot_Ext structure
 ** \param[in]  dbiso         Iso-potential Db structure
 ** \param[in]  dbgrd         Gradient Db structure
 ** \param[in]  dbtgt         Tangent Db structure (optional)
 ** \param[in]  dbgrid        Output Db structure (for External drift)
 ** \param[in]  model         Model structure
 ** \param[in]  isimu         Rank of the simulation (or -1)
 ** \param[in]  nbsimu        Number of simulations
 ** \param[in]  refpot        Potential at Reference point
 ** \param[in]  zdual         Dual vector (Dimension: nequa)
 ** \param[in]  rhs           R.H.S. vector (Dimension: nequa * 4)
 **
 *****************************************************************************/
static void st_check_data(Pot_Env *pot_env,
                          Pot_Ext *pot_ext,
                          Db *dbiso,
                          Db *dbgrd,
                          Db *dbtgt,
                          DbGrid *dbgrid,
                          Model *model,
                          int isimu,
                          int nbsimu,
                          double refpot,
                          double *zdual,
                          double *rhs)
{
  VectorDouble result(4);

  /* Preliminary check */

  if (VERBOSE) mestitle(0, "Information completed at Data Points");

  /* For the Iso-Potential file */

  if (dbiso != nullptr)
  {
    if (VERBOSE) mestitle(1, "Iso-Potential Information");

    int rank = 0;
    for (int ic = 0; ic < pot_env->nlayers; ic++)
    {
      for (int j = 0; j < pot_env->nb_per_layer[ic]; j++, rank++)
      {
        OptDbg::setCurrentIndex(rank);
        int iech = dbiso->getRankRelativeToAbsolute(rank);
        st_calc_point(pot_env, pot_ext, 1, dbiso, dbgrd, dbtgt, dbgrid, model,
                      zdual, rhs, dbiso, iech, result);
        result[0] -= refpot;

        // Printout (conditional)

        if (VERBOSE)
        {
          // Convert into simulation error

          if (nbsimu > 0)
            result[0] = dbiso->getSimvar(ELoc::SIMU,iech,isimu,0,0,nbsimu,1)
                        - result[0];

          // Print the results

          message(" %d - %d - Coor =",ic+1, j+1);
          for (int idim = 0; idim < pot_env->ndim; idim++)
            message(" %lf", ISO_COO(ic, j, idim));
          st_print_result(pot_env, isimu, result.data(), TEST);
        }
      }
      OptDbg::setCurrentIndex(-1);
    }
  }

  /* For the Gradient file */

  if (dbgrd != nullptr)
  {
    if (VERBOSE) mestitle(1, "Gradient Information");

    for (int ig = 0; ig < pot_env->ngrd; ig++)
    {
      OptDbg::setCurrentIndex(ig);
      int iech = dbgrd->getRankRelativeToAbsolute(ig);
      st_calc_point(pot_env, pot_ext, 1, dbiso, dbgrd, dbtgt, dbgrid, model,
                    zdual, rhs, dbgrd, iech, result);
      result[0] -= refpot;

      // Printout (optional)

      if (VERBOSE)
      {

        // Convert into simulation error

        if (nbsimu > 0)
        {
          for (int idim = 0; idim < pot_env->ndim; idim++)
            result[1 + idim] = (dbgrd->getSimvar(ELoc::SIMU, iech,
                                                 isimu + idim * nbsimu, 0, 0,
                                                 pot_env->ndim * nbsimu, 1)
                                - result[1 + idim]);
        }

        // Print the results

        message(" %2d - Coor =", ig+1);
        for (int idim = 0; idim < pot_env->ndim; idim++)
          message(" %lf", GRD_COO(ig, idim));
        st_print_result(pot_env, isimu, result.data(), TEST);
      }
    }
    OptDbg::setCurrentIndex(-1);
  }

  /* For the Tangent file */

  if (dbtgt != nullptr)
  {
    if (VERBOSE) mestitle(1, "Tangent Information");

    for (int it = 0; it < pot_env->ntgt; it++)
    {
      OptDbg::setCurrentIndex(it);
      int iech = dbtgt->getRankRelativeToAbsolute(it);
      if (!dbtgt->isActive(iech)) continue;
      st_calc_point(pot_env, pot_ext, 1, dbiso, dbgrd, dbtgt, dbgrid, model,
                    zdual, rhs, dbtgt, iech, result);
      result[0] -= refpot;

      // Printout (conditional) 

      if (VERBOSE)
      {
        double tgte = 0.;
        for (int idim = 0; idim < pot_env->ndim; idim++)
          tgte += result[1 + idim] * dbtgt->getLocVariable(ELoc::TGTE,iech, idim);

        // Print the results

        message(" %2d - Coor =", it+1);
        for (int idim = 0; idim < pot_env->ndim; idim++)
          message(" %lf", TGT_COO(it, idim));
        st_print_result(pot_env, isimu, result.data(), tgte);
      }
    }
    OptDbg::setCurrentIndex(-1);
  }

  return;
}

/****************************************************************************/
/*!
 **  Calculate the estimation at the potential at first point of first potential
 **
 ** \return The Potential value at first point of first iso-potential
 **
 ** \param[in]  pot_env       Pot_Env structure
 ** \param[in]  pot_ext       Pot_Ext structure
 ** \param[in]  dbiso         Iso-potential Db structure
 ** \param[in]  dbgrd         Gradient Db structure
 ** \param[in]  dbtgt         Tangent Db structure (optional)
 ** \param[in]  dbgrid        Ouput Db structure (for external drift)
 ** \param[in]  model         Model structure
 ** \param[in]  zdual         Dual vector (Dimension: nequa)
 ** \param[in]  rhs           R.H.S. vector (Dimension: nequa * 4)
 **
 *****************************************************************************/
static double st_evaluate_refpot(Pot_Env *pot_env,
                                 Pot_Ext *pot_ext,
                                 Db *dbiso,
                                 Db *dbgrd,
                                 Db *dbtgt,
                                 DbGrid *dbgrid,
                                 Model *model,
                                 double *zdual,
                                 double *rhs)
{
  VectorDouble result(4);
  int ip1, ic;

  if (dbiso == nullptr) return (TEST);

  // Calculate the reference values for iso-potentials

  ic = 0;
  ip1 = IAD_ISO(ic, 0);
  st_calc_point(pot_env, pot_ext, 0, dbiso, dbgrd, dbtgt, dbgrid, model,
                zdual, rhs, dbiso, ip1, result);
  return (result[0]);
}

/****************************************************************************/
/*!
 **  Calculate the estimation at the iso-potential samples
 **
 ** \param[in]  pot_env       Pot_Env structure
 ** \param[in]  pot_ext       Pot_Ext structure
 ** \param[in]  dbiso         Iso-potential Db structure
 ** \param[in]  dbgrd         Gradient Db structure
 ** \param[in]  dbtgt         Tangent Db structure (optional)
 ** \param[in]  dbgrid        Ouput Db structure (for external drift)
 ** \param[in]  model         Model structure
 ** \param[in]  refpot        Potential at Reference point
 ** \param[in]  isimu         Rank of the simulation (or -1)
 ** \param[in]  nbsimu        Number of simulations (or 0)
 ** \param[in]  zdual         Dual vector (Dimension: nequa)
 ** \param[in]  rhs           R.H.S. vector (Dimension: nequa * 4)
 **
 ** \param[out] potval        Array of Potential values
 **
 *****************************************************************************/
static void st_evaluate_potval(Pot_Env *pot_env,
                               Pot_Ext *pot_ext,
                               Db *dbiso,
                               Db *dbgrd,
                               Db *dbtgt,
                               DbGrid *dbgrid,
                               Model *model,
                               double refpot,
                               int isimu,
                               int nbsimu,
                               double *zdual,
                               double *rhs,
                               double *potval)
{
  VectorDouble result(4);
  int ip1;

  if (dbiso == nullptr) return;

  // Calculate the reference values for isopotentials 

  for (int ic = 0; ic < pot_env->nlayers; ic++)
  {
    ip1 = IAD_ISO(ic, 0);
    st_calc_point(pot_env, pot_ext, 0, dbiso, dbgrd, dbtgt, dbgrid, model,
                  zdual, rhs, dbiso, ip1, result);

    // Convert into simulation error

    if (nbsimu > 0)
      result[0] = (dbiso->getSimvar(ELoc::SIMU, ip1, isimu, 0, 0, nbsimu, 1)
          - result[0]);

    // Center to the reference potential

    result[0] -= refpot;

    // Store in the 'potval' array

    potval[ic] = result[0];
  }

  // Sort them by ascending order 

  ut_sort_double(0, pot_env->nlayers, NULL, potval);

  return;
}

static void st_save_result_on_data(Pot_Env* pot_env,
                                   Db* db,
                                   int nvar,
                                   double value,
                                   const ELoc& loctype_pot,
                                   const ELoc& loctype_grad,
                                   VectorInt& uid_pot,
                                   VectorInt& uid_grad)
{
  int uid;
  uid_pot.clear();
  uid_grad.clear();
  if (db == nullptr) return;

  if (pot_env->flag_pot)
  {
    uid = db->addColumnsByConstant(nvar, value, "Potential",
                                   loctype_pot);
    uid_pot.push_back(uid);
  }
  if (pot_env->flag_grad)
  {
    uid = db->addColumnsByConstant(nvar * pot_env->ndim, value, "Gradients",
                                   loctype_grad);
    for (int idim = 0; idim < pot_env->ndim; idim++)
      uid_grad.push_back(uid + idim);
  }
}

/****************************************************************************/
/*!
 **  Potential estimation
 **
 ** \return  Error return code
 **
 ** \param[in]  dbiso      Iso-potential Db structure
 ** \param[in]  dbgrd      Gradient Db structure
 ** \param[in]  dbtgt      Tangent Db structure (optional)
 ** \param[in]  dbout      Output Db structure
 ** \param[in]  model      Model structure
 ** \param[in]  neigh      ANeigh structure
 ** \param[in]  nugget_grd Nugget effect for Gradients
 ** \param[in]  nugget_tgt Nugget effect for Tangents
 ** \param[in]  flag_pot   True if the Potential must be estimated
 ** \param[in]  flag_grad  True if the gradient must also be estimated
 ** \param[in]  flag_trans True if the estimation result must be translated
 **                        into layer number
 ** \param[in]  flag_save_data True if the Potential / Gradient must be
 **                        saved on any Information file
 ** \param[in]  opt_part   Option to exhibit only a part of estimation:
 ** \li                       0 : the whole estimation
 ** \li                       1 : the gradient contribution only
 ** \li                       2 : the tangent contribution only
 ** \li                       3 : the isovalues contribution only
 ** \li                       4 : the drift contribution only
 ** \li                       5 : the external drift contribution only
 ** \param[in]  verbose    Verbose option
 **
 ** \remark The results will be stored in the dbout file
 ** \remark - the estimation in the variable ELoc::Z
 ** \remark - the gradient components in the variables ELoc::GRD
 **
 *****************************************************************************/
int potential_kriging(Db *dbiso,
                      Db *dbgrd,
                      Db *dbtgt,
                      DbGrid *dbout,
                      Model *model,
                      ANeigh *neigh,
                      double nugget_grd,
                      double nugget_tgt,
                      bool flag_pot,
                      bool flag_grad,
                      bool flag_trans,
                      bool flag_save_data,
                      int opt_part,
                      bool verbose)
{
  int error, nequa;
  VectorInt uid_iso_pot, uid_iso_grad;
  VectorInt uid_grd_pot, uid_grd_grad;
  VectorInt uid_tgt_pot, uid_tgt_grad;
  VectorInt uid_out_pot, uid_out_grad;
  double *lhs, *zval, *zdual, *rhs, *potval, refpot;
  Pot_Env pot_env;
  Pot_Ext pot_ext;

  // Initialization

  error = 1;
  lhs = zval = zdual = rhs = potval = nullptr;
  st_potenv_manage(&pot_env, flag_pot, flag_grad, flag_trans, opt_part, verbose);
  st_potext_manage(0, &pot_ext, 0, 0., NULL);
  st_potenv_define(&pot_env, &pot_ext, dbiso, dbgrd, dbtgt, dbout);

  // Preliminary checks

  if (krige_koption_manage(1, 1, EKrigOpt::POINT, 1, VectorInt()))
    goto label_end;
  if (!st_potenv_valid(&pot_env, &pot_ext, dbiso, dbgrd, dbtgt, dbout, model,
                       neigh)) goto label_end;

  // Count the gradients and the tangents

  if (st_update_isopot(dbiso, &pot_env)) goto label_end;
  if (st_update_gradient(dbgrd, &pot_env)) goto label_end;
  if (st_update_tangent(dbtgt, &pot_env)) goto label_end;
  if (st_update_model(model, &pot_env)) goto label_end;
  if (st_update_final(model, &pot_env)) goto label_end;

  // Allocating the output variables

  st_save_result_on_data(&pot_env, dbout, 1, TEST, ELoc::Z, ELoc::G,
                         uid_out_pot, uid_out_grad);
  if (flag_save_data)
  {
    st_save_result_on_data(&pot_env, dbiso, 1, TEST, ELoc::UNKNOWN,
                           ELoc::UNKNOWN, uid_iso_pot, uid_iso_grad);
    st_save_result_on_data(&pot_env, dbgrd, 1, TEST, ELoc::UNKNOWN,
                           ELoc::UNKNOWN, uid_grd_pot, uid_grd_grad);
    st_save_result_on_data(&pot_env, dbtgt, 1, TEST, ELoc::UNKNOWN,
                           ELoc::UNKNOWN, uid_tgt_pot, uid_tgt_grad);
  }

  // Core allocation

  nequa = pot_env.nequa;
  lhs = (double*) mem_alloc(sizeof(double) * nequa * nequa, 0);
  if (lhs == nullptr) goto label_end;
  zval = (double*) mem_alloc(sizeof(double) * nequa, 0);
  if (zval == nullptr) goto label_end;
  zdual = (double*) mem_alloc(sizeof(double) * nequa, 0);
  if (zdual == nullptr) goto label_end;
  rhs = (double*) mem_alloc(sizeof(double) * nequa * 4, 0);
  if (rhs == nullptr) goto label_end;
  potval = (double*) mem_alloc(sizeof(double) * pot_env.nlayers, 0);
  if (potval == nullptr) goto label_end;

  // Establish the cokriging system

  if (st_build_lhs(&pot_env, &pot_ext, dbout, model,
                   nugget_grd, nugget_tgt, lhs)) goto label_end;
  if (OptDbg::isReferenceDefined() || OptDbg::query(EDbg::KRIGING))
    krige_lhs_print(0, nequa, nequa, NULL, lhs);

  // Invert the matrix

  if (matrix_invert(lhs, nequa, -1)) goto label_end;
  if (OptDbg::isReferenceDefined() || OptDbg::query(EDbg::KRIGING))
    print_matrix("[LHS]-1", 0, 1, nequa, nequa, NULL, lhs);

  // Establish the data vector and get the dual form

  st_fill_dual(&pot_env, zval);
  if (OptDbg::isReferenceDefined() || OptDbg::query(EDbg::KRIGING))
    print_matrix("\n[Z]", 0, 1, 1, nequa, NULL, zval);
  matrix_product_safe(nequa, nequa, 1, lhs, zval, zdual);
  if (OptDbg::isReferenceDefined() || OptDbg::query(EDbg::KRIGING))
    print_matrix("\n[Z] *%* [LHS]-1", 0, 1, 1, nequa, NULL, zdual);

  // Evaluate Potential at Reference point

  refpot = st_evaluate_refpot(&pot_env, &pot_ext, dbiso, dbgrd, dbtgt, dbout,
                              model, zdual, rhs);

  // Check that the information is fulfilled correctly

  if (VERBOSE)
    st_check_data(&pot_env, &pot_ext, dbiso, dbgrd, dbtgt, dbout, model, -1, 0,
                  refpot, zdual, rhs);

  // Get the Potential value at the iso-potential samples

  st_evaluate_potval(&pot_env, &pot_ext, dbiso, dbgrd, dbtgt, dbout,
                     model, refpot, -1, 0, zdual, rhs, potval);

  // Perform the estimation

  st_estimate_result(&pot_env, &pot_ext, flag_grad, dbiso,
                     dbgrd, dbtgt, dbout, model, refpot, zdual, rhs, potval);
  if (flag_save_data)
  {
    st_estimate_data(&pot_env, &pot_ext, dbiso, dbgrd, dbtgt, dbout,
                     model, refpot, zdual, rhs,
                     dbiso, uid_iso_pot, uid_iso_grad);
    st_estimate_data(&pot_env, &pot_ext, dbiso, dbgrd, dbtgt, dbout,
                     model, refpot, zdual, rhs,
                     dbgrd, uid_grd_pot, uid_grd_grad);
    st_estimate_data(&pot_env, &pot_ext, dbiso, dbgrd, dbtgt, dbout,
                     model, refpot, zdual, rhs,
                     dbtgt, uid_tgt_pot, uid_tgt_grad);
  }

  // Set the error return code

  error = 0;

  label_end:
  st_potext_manage(-1, &pot_ext, 0, 0., NULL);
  (void) krige_koption_manage(-1, 1, EKrigOpt::POINT, 1, VectorInt());
  lhs = (double*) mem_free((char* ) lhs);
  zval = (double*) mem_free((char* ) zval);
  zdual = (double*) mem_free((char* ) zdual);
  rhs = (double*) mem_free((char* ) rhs);
  potval = (double*) mem_free((char* ) potval);
  return (error);
}

/****************************************************************************/
/*!
 **  Transform the Estimation variable into a distance to the isoline
 **
 ** \return  Error return code
 **
 ** \param[in]  dbout        Output Db structure
 **
 *****************************************************************************/
static int st_distance_to_isoline(DbGrid *dbout)

{
  int radius, seed, memo;
  double value, eps;

  // Initializations 
  radius = 1;
  eps = 1.e-3;
  seed = 3432521;
  memo = law_get_random_seed();

  // Highlight the isoline of interest
  for (int iech = 0; iech < dbout->getSampleNumber(); iech++)
  {
    value = dbout->getLocVariable(ELoc::Z,iech, 0);
    if (!FFFF(value) && ABS(value) > eps) dbout->setLocVariable(ELoc::Z,iech, 0, TEST);
  }

  // Calculate the distance
  if (db_grid_fill(dbout, 3, seed, radius)) return (1);

  law_set_random_seed(memo);
  return (0);
}

/****************************************************************************/
/*!
 **  Potential simulations
 **
 ** \return  Error return code
 **
 ** \param[in]  dbiso        Iso-potential Db structure
 ** \param[in]  dbgrd        Gradient Db structure
 ** \param[in]  dbtgt        Tangent Db structure (optional)
 ** \param[in]  dbout        Output Db structure
 ** \param[in]  model        Model structure
 ** \param[in]  neigh        ANeigh structure
 ** \param[in]  nugget_grd   Nugget effect for Gradients
 ** \param[in]  nugget_tgt   Nugget effect for Tangents
 ** \param[in]  dist_tempere Distance for tempering simulations (or TEST)
 ** \param[in]  flag_trans   True if the estimation result must be translated
 **                          into layer number
 ** \param[in]  seed         Seed for the random number generator
 ** \param[in]  nbsimu       Number of simulations
 ** \param[in]  nbtuba       Number of turning bands
 ** \param[in]  verbose      Verbose option
 **
 ** \remark The simulations will be stored in the dbout file (ELoc::SIMU)
 **
 *****************************************************************************/
int potential_simulate(Db *dbiso,
                       Db *dbgrd,
                       Db *dbtgt,
                       DbGrid *dbout,
                       Model *model,
                       ANeigh *neigh,
                       double nugget_grd,
                       double nugget_tgt,
                       double dist_tempere,
                       bool flag_trans,
                       int seed,
                       int nbsimu,
                       int nbtuba,
                       bool verbose)
{
  int error, nequa, nlayers, flag_tempere;
  double *lhs, *zval, *zduals, *zdual, *rhs, *potsim, *potval, delta, refpot;
  Pot_Env pot_env;
  Pot_Ext pot_ext;
  VectorInt uid_out_pot, uid_out_grad;
  VectorInt uid_iso_pot, uid_iso_grad;
  VectorInt uid_grd_pot, uid_grd_grad;
  VectorInt uid_tgt_pot, uid_tgt_grad;

  // Initialization

  error = 1;
  lhs = zval = zduals = zdual = rhs = potsim = potval = nullptr;
  st_potenv_manage(&pot_env, true, false, flag_trans, 0, verbose);
  st_potext_manage(0, &pot_ext, 0, 0., NULL);
  st_potenv_define(&pot_env, &pot_ext, dbiso, dbgrd, dbtgt, dbout);
  law_set_random_seed(seed);
  flag_tempere = !FFFF(dist_tempere);
  refpot = 0.;

  // Preliminary checks

  if (krige_koption_manage(1, 1, EKrigOpt::POINT, 1, VectorInt()))
    goto label_end;
  if (db_extension_diag(dbiso, &delta)) goto label_end;
  delta /= 1000.;

  if (!st_potenv_valid(&pot_env, &pot_ext, dbiso, dbgrd, dbtgt, dbout, model,
                       neigh)) goto label_end;

  // Count the gradients and the tangents

  if (st_update_isopot(dbiso, &pot_env)) goto label_end;
  if (st_update_gradient(dbgrd, &pot_env)) goto label_end;
  if (st_update_tangent(dbtgt, &pot_env)) goto label_end;
  if (st_update_model(model, &pot_env)) goto label_end;
  if (st_update_final(model, &pot_env)) goto label_end;
  nlayers = pot_env.nlayers;

  /* Add the attributes for storing the results */

  st_save_result_on_data(&pot_env, dbout, nbsimu, 0., ELoc::SIMU, ELoc::UNKNOWN,
                         uid_out_pot, uid_out_grad);
  st_save_result_on_data(&pot_env, dbiso, 2*nbsimu, 0., ELoc::SIMU, ELoc::UNKNOWN,
                         uid_iso_pot, uid_iso_grad);
  st_save_result_on_data(&pot_env, dbgrd, 2*nbsimu, 0., ELoc::SIMU, ELoc::UNKNOWN,
                         uid_grd_pot, uid_grd_grad);
  st_save_result_on_data(&pot_env, dbtgt, 2*nbsimu, 0., ELoc::SIMU, ELoc::UNKNOWN,
                         uid_tgt_pot, uid_tgt_grad);
  if (flag_tempere)
    (void) dbout->addColumnsByConstant(1, TEST, String(), ELoc::Z);

  /* Processing the non-conditional simulation over the iso-values */

  {
    CalcSimuTurningBands situba_new(nbsimu, nbtuba, seed);
    if (situba_new.simulatePotential(dbiso, dbgrd, dbtgt, dbout, model, delta))
      goto label_end;
  }

  // Core allocation

  nequa = pot_env.nequa;
  lhs = (double*) mem_alloc(sizeof(double) * nequa * nequa, 0);
  if (lhs == nullptr) goto label_end;
  zval = (double*) mem_alloc(sizeof(double) * nequa * nbsimu, 0);
  if (zval == nullptr) goto label_end;
  zduals = (double*) mem_alloc(sizeof(double) * nequa * nbsimu, 0);
  if (zduals == nullptr) goto label_end;
  rhs = (double*) mem_alloc(sizeof(double) * nequa * 4, 0);
  if (rhs == nullptr) goto label_end;
  potsim = (double*) mem_alloc(sizeof(double) * nlayers * nbsimu, 0);
  if (potsim == nullptr) goto label_end;
  potval = (double*) mem_alloc(sizeof(double) * nlayers, 0);
  if (potval == nullptr) goto label_end;
  if (flag_tempere)
  {
    zdual = (double*) mem_alloc(sizeof(double) * nequa, 0);
    if (zdual == nullptr) goto label_end;
  }

  // Establish the cokriging system

  if (st_build_lhs(&pot_env, &pot_ext, dbout, model,
                   nugget_grd, nugget_tgt, lhs)) goto label_end;

  // Invert the matrix

  if (matrix_invert(lhs, nequa, -1)) goto label_end;
  if (OptDbg::isReferenceDefined() || OptDbg::query(EDbg::KRIGING))
    print_matrix("Inverted LHS", 0, 1, nequa, nequa, NULL, lhs);

  if (flag_tempere)
  {

    // Establish the data vector and get the dual form

    st_fill_dual(&pot_env, zval);
    if (OptDbg::isReferenceDefined() || OptDbg::query(EDbg::KRIGING))
      print_matrix("\n[Z]", 0, 1, 1, nequa, NULL, zval);
    matrix_product_safe(nequa, nequa, 1, lhs, zval, zdual);
    if (OptDbg::isReferenceDefined() || OptDbg::query(EDbg::KRIGING))
      print_matrix("\n[Z] *%* [A]-1", 0, 1, 1, nequa, NULL, zdual);

    // Evaluate Potential at Reference point

    refpot = st_evaluate_refpot(&pot_env, &pot_ext, dbiso, dbgrd, dbtgt, dbout,
                                model, zdual, rhs);

    // Get the Estimated Potential value at the iso-potential samples

    st_evaluate_potval(&pot_env, &pot_ext, dbiso, dbgrd, dbtgt, dbout, model,
                       refpot, -1, 0, zdual, rhs, potval);

    // Perform the estimation 

    st_estimate_result(&pot_env, &pot_ext, 0, dbiso, dbgrd, dbtgt, dbout, model,
                       refpot, zdual, rhs, potval);

    // Transform the Estimation variable into a distance to the isoline
    if (st_distance_to_isoline(dbout)) goto label_end;
  }

  // Establish the simulated error vector and get the dual form

  st_fill_dual_simulation(&pot_env, dbiso, dbgrd, dbtgt, nbsimu, zval);
  if (OptDbg::isReferenceDefined() || OptDbg::query(EDbg::KRIGING))
    print_matrix("\n[Simu-Err]", 0, 1, nbsimu, nequa, NULL, zval);
  matrix_product_safe(nequa, nequa, nbsimu, lhs, zval, zduals);
  if (OptDbg::isReferenceDefined() || OptDbg::query(EDbg::KRIGING))
    print_matrix("\n[Simu-Err] *%* [A]-1", 0, 1, nbsimu, nequa, NULL, zduals);

  // Get the Simulated Potential value at the iso-potential samples

  for (int isimu = 0; isimu < nbsimu; isimu++)
  {

    // Calculate the simulated value at the reference point

    refpot = st_evaluate_refpot(&pot_env, &pot_ext, dbiso, dbgrd, dbtgt, dbout,
                                model, &ZDUALS(isimu, 0), rhs);

    // Check that the information is fulfilled correctly

    st_check_data(&pot_env, &pot_ext, dbiso, dbgrd, dbtgt, dbout, model, isimu,
                  nbsimu, refpot, &ZDUALS(isimu, 0), rhs);

    // Calculate the simulated iso-value

    st_evaluate_potval(&pot_env, &pot_ext, dbiso, dbgrd, dbtgt, dbout, model,
                       refpot, isimu, nbsimu, &ZDUALS(isimu, 0), rhs,
                       &POTSIM(isimu, 0));
  }

  // Perform the conditional simulations on the grid

  st_simcond(&pot_env, &pot_ext, dist_tempere, flag_trans, nbsimu, dbiso, dbgrd,
             dbtgt, dbout, model, refpot, potsim, zdual, zduals, rhs);

  // Set the error return code

  error = 0;

  label_end: if (flag_tempere) dbout->deleteColumnsByLocator(ELoc::Z);
  st_potext_manage(-1, &pot_ext, 0, 0., NULL);
  (void) krige_koption_manage(-1, 1, EKrigOpt::POINT, 1, VectorInt());
  lhs = (double*) mem_free((char* ) lhs);
  zval = (double*) mem_free((char* ) zval);
  zdual = (double*) mem_free((char* ) zdual);
  zduals = (double*) mem_free((char* ) zduals);
  rhs = (double*) mem_free((char* ) rhs);
  potval = (double*) mem_free((char* ) potval);
  potsim = (double*) mem_free((char* ) potsim);
  return (error);
}

/****************************************************************************/
/*!
 **  Potential cross-validation
 **
 ** \return  Error return code
 **
 ** \param[in]  dbiso          Iso-potential Db structure
 ** \param[in]  dbgrd          Gradient Db structure
 ** \param[in]  dbtgt          Tangent Db structure (optional)
 ** \param[in]  model          Model structure
 ** \param[in]  neigh          ANeigh structure
 ** \param[in]  nugget_grd     Nugget effect for Gradients
 ** \param[in]  nugget_tgt     Nugget effect for Tangents
 ** \param[in]  flag_dist_conv Flag for converting into distance
 ** \param[in]  verbose        Verbose option
 **
 *****************************************************************************/
int potential_xvalid(Db *dbiso,
                     Db *dbgrd,
                     Db *dbtgt,
                     Model *model,
                     ANeigh *neigh,
                     double nugget_grd,
                     double nugget_tgt,
                     int flag_dist_conv,
                     bool verbose)
{
  int error, nequa, nvar;
  double *lhs, *zval, *zdual, *rhs, *lhs_orig, *lhs_aux;
  Pot_Env pot_env;
  Pot_Ext pot_ext;

  // Initialization

  error = 1;
  lhs = zval = zdual = rhs = lhs_orig = lhs_aux = nullptr;
  st_potenv_manage(&pot_env, true, false, false, 0, verbose);
  st_potext_manage(0, &pot_ext, 0, 0., NULL);
  st_potenv_define(&pot_env, &pot_ext, dbiso, dbgrd, dbtgt, NULL);

  // Preliminary checks

  if (krige_koption_manage(1, 1, EKrigOpt::POINT, 1, VectorInt()))
    goto label_end;
  if (!st_potenv_valid(&pot_env, &pot_ext, dbiso, dbgrd, dbtgt, NULL, model,
                       neigh)) goto label_end;

  // Count the gradients and the tangents

  if (st_update_isopot(dbiso, &pot_env)) goto label_end;
  if (st_update_gradient(dbgrd, &pot_env)) goto label_end;
  if (st_update_tangent(dbtgt, &pot_env)) goto label_end;
  if (st_update_model(model, &pot_env)) goto label_end;
  if (st_update_final(model, &pot_env)) goto label_end;

  // Allocating the output variables

  nvar = 2;
  if (flag_dist_conv) nvar = 4;
  (void) dbiso->addColumnsByConstant(nvar, TEST, String(), ELoc::Z);

  // Core allocation

  nequa = pot_env.nequa;
  lhs = (double*) mem_alloc(sizeof(double) * nequa * nequa, 0);
  if (lhs == nullptr) goto label_end;
  zval = (double*) mem_alloc(sizeof(double) * nequa, 0);
  if (zval == nullptr) goto label_end;
  zdual = (double*) mem_alloc(sizeof(double) * nequa, 0);
  if (zdual == nullptr) goto label_end;
  rhs = (double*) mem_alloc(sizeof(double) * nequa * 4, 0);
  if (rhs == nullptr) goto label_end;
  if (flag_dist_conv)
  {
    lhs_orig = (double*) mem_alloc(sizeof(double) * nequa * nequa, 0);
    if (lhs_orig == nullptr) goto label_end;
    lhs_aux = (double*) mem_alloc(sizeof(double) * nequa * nequa, 0);
    if (lhs_aux == nullptr) goto label_end;
  }

  // Establish the cokriging system

  if (st_build_lhs(&pot_env, &pot_ext, nullptr, model,
                   nugget_grd, nugget_tgt, lhs)) goto label_end;

  // Save the matrix (used for converting into distance)

  if (flag_dist_conv)
    (void) memcpy(lhs_orig, lhs, sizeof(double) * nequa * nequa);

  // Invert the matrix

  if (matrix_invert(lhs, nequa, -1)) goto label_end;
  if (OptDbg::isReferenceDefined() || OptDbg::query(EDbg::KRIGING))
    print_matrix("Inverted LHS", 0, 1, nequa, nequa, NULL, lhs);

  // Establish the data vector and get the dual form

  st_fill_dual(&pot_env, zval);
  if (OptDbg::isReferenceDefined() || OptDbg::query(EDbg::KRIGING))
    print_matrix("\n[Z]", 0, 1, 1, nequa, NULL, zval);
  matrix_product_safe(nequa, nequa, 1, lhs, zval, zdual);
  if (OptDbg::isReferenceDefined() || OptDbg::query(EDbg::KRIGING))
    print_matrix("\n[Z] *%* [A]-1", 0, 1, 1, nequa, NULL, zdual);

  /* Process the estimate at masked-off isovalues */

  st_xvalid_potential(&pot_env, &pot_ext, dbiso, dbgrd, dbtgt, model, lhs,
                      flag_dist_conv, zval, lhs_orig, lhs_aux, rhs, zdual);

  // Set the error return code

  error = 0;

  label_end:
  st_potext_manage(-1, &pot_ext, 0, 0., NULL);
  (void) krige_koption_manage(-1, 1, EKrigOpt::POINT, 1, VectorInt());
  lhs = (double*) mem_free((char* ) lhs);
  zval = (double*) mem_free((char* ) zval);
  zdual = (double*) mem_free((char* ) zdual);
  rhs = (double*) mem_free((char* ) rhs);
  lhs_orig = (double*) mem_free((char* ) lhs_orig);
  lhs_aux = (double*) mem_free((char* ) lhs_aux);
  return (error);
}

/****************************************************************************/
/*!
 **  Print the type of information for the Potential covariance
 **
 ** \param[in]  rank      Rank of the point
 ** \param[in]  type      Type of the first point
 **                       1 for gradient; 2 for tangent; 3 for isopotential
 **
 *****************************************************************************/
static void st_print_type(int rank, int type)
{
  message("Data Set #%d: ", rank);
  switch (type)
  {
    case -1:
      message("Target Gradient\n");
      break;
    case 1:
      message("Gradient\n");
      break;
    case 2:
      message("Tangent\n");
      break;
    case -3:
      message("Target IsoPotential\n");
      break;
    case 3:
      message("IsoPotential\n");
      break;
  }
}

/****************************************************************************/
/*!
 **  Potential covariance
 **
 ** \return Error return code
 **
 ** \param[in]  model      Model structure
 ** \param[in]  verbose    Verbose flag
 ** \param[in]  type1      Type of the first point
 **                        1 for gradient; 2 for tangent; 3 for isopotential
 ** \param[in]  x10        Coordinates of the centering for first point
 ** \param[in]  x1p        Coordinates of the first point
 ** \param[in]  tx1        Tangent values at the first point
 ** \param[in]  type2      Type of the second point
 **                        1 for gradient; 2 for tangent; 3 for isopotential
 **                        (Sign is negative for target point)
 ** \param[in]  x20        Coordinates of the centering for second point
 ** \param[in]  x2p        Coordinates of the second point
 ** \param[in]  tx2        Tangent values at the second point
 **
 ** \param[out] covtab     Array of returned values (dimensionned to ndim*ndim)
 **
 *****************************************************************************/
int potential_cov(Model *model,
                  bool verbose,
                  int type1,
                  const VectorDouble& x10,
                  const VectorDouble& x1p,
                  const VectorDouble& tx1,
                  int type2,
                  const VectorDouble& x20,
                  const VectorDouble& x2p,
                  const VectorDouble& tx2,
                  VectorDouble& covtab)

{
  int idim, jdim, ecr, lec, i;
  double dd[3] = { 0., 0., 0. };
  VectorDouble covGp(3, 0.);
  VectorDouble cov2Gp(3, 0.);
  VectorDouble covGG(9, 0.);
  double covar = 0;
  double covar1 = 0;
  double covar2 = 0;
  double covar3 = 0;
  double covar4 = 0;

  // Preliminary checks

  VERBOSE = verbose;
  int ndim = model->getDimensionNumber();
  covtab.resize(ndim * ndim, TEST);

  /* Preliminary checks */

  if (type1 < 1 || type1 > 3)
  {
    messerr("Argument 'type1'(%d) must be equal to 1, 2 or 3", type1);
    return (1);
  }
  if (type2 < -3 || type2 > 3 || type2 == -2 || type2 == 0)
  {
    messerr("Argument 'type2'(%d) must be equal to -3,-1,1,2 or 3", type2);
    return (1);
  }

  /* Optional printout */

  if (VERBOSE)
  {
    st_print_type(1, type1);
    if (! x10.empty()) print_matrix("x10", 0, 1, 1, ndim, NULL, x10.data());
    if (! x1p.empty()) print_matrix("x1p", 0, 1, 1, ndim, NULL, x1p.data());
    if (! tx1.empty()) print_matrix("tx1", 0, 1, 1, ndim, NULL, tx1.data());
    st_print_type(2, type2);
    if (! x20.empty()) print_matrix("x20", 0, 1, 1, ndim, NULL, x20.data());
    if (! x2p.empty()) print_matrix("x2p", 0, 1, 1, ndim, NULL, x2p.data());
    if (! tx2.empty()) print_matrix("tx2", 0, 1, 1, ndim, NULL, tx2.data());
  }

  /* Dispatch */

  int n1 = 1;
  int n2 = 1;
  switch (type1)
  {
    case 1:                     // 1-Gradient
      n1 = ndim;
      switch (type2)
      {
        case 1:                 // 2-Gradient
        case -1:                // 2-Gradient-Target
          for (idim = 0; idim < ndim; idim++)
            dd[idim] = x1p[idim] - x2p[idim];
          st_cov(model, 1, dd[0], dd[1], dd[2], covar, covGp, covGG);
          ecr = lec = 0;
          for (idim = 0; idim < 3; idim++)
            for (jdim = 0; jdim < 3; jdim++, lec++)
            {
              if (idim < ndim && jdim < ndim) covtab[ecr++] = covGG[lec];
            }
          n2 = ndim;
          break;

        case 2:                 // 2-Tangent
          for (idim = 0; idim < ndim; idim++)
            dd[idim] = x1p[idim] - x2p[idim];
          st_cov(model, 1, dd[0], dd[1], dd[2], covar, covGp, covGG);
          for (idim = 0; idim < ndim; idim++)
          {
            i = 3 * idim;
            covtab[idim] = matrix_UV(ndim, tx2[0], tx2[1], tx2[2], covGG[i + 0],
                                     covGG[i + 1], covGG[i + 2]);
          }
          n2 = 1;
          break;

        case 3:                 // 2-IsoPotential
          for (idim = 0; idim < ndim; idim++)
            dd[idim] = x1p[idim] - x2p[idim];
          st_cov(model, 1, dd[0], dd[1], dd[2], covar, cov2Gp, covGG);
          for (idim = 0; idim < ndim; idim++)
            dd[idim] = x1p[idim] - x20[idim];
          st_cov(model, 1, dd[0], dd[1], dd[2], covar, covGp, covGG);
          for (idim = 0; idim < ndim; idim++)
            covtab[idim] = cov2Gp[idim] - covGp[idim];
          n2 = 1;
          break;

        case -3:                 // 2-IsoPotential-Target
          for (idim = 0; idim < ndim; idim++)
            dd[idim] = x1p[idim] - x2p[idim];
          st_cov(model, 1, dd[0], dd[1], dd[2], covar, cov2Gp, covGG);
          for (idim = 0; idim < ndim; idim++)
            covtab[idim] = cov2Gp[idim];
          n2 = 1;
          break;
      }
      break;

    case 2:                     // 1-Tangent
      n1 = 1;
      switch (type2)
      {
        case 1:                 // 2-Gradient
        case -1:                // 2-Gradient-Target
          for (idim = 0; idim < ndim; idim++)
            dd[idim] = x1p[idim] - x2p[idim];
          st_cov(model, 1, dd[0], dd[1], dd[2], covar, covGp, covGG);
          for (idim = 0; idim < ndim; idim++)
          {
            i = 3 * idim;
            covtab[idim] = matrix_UV(ndim,
                                     tx1[0], tx1[1], tx1[2], covGG[i + 0],
                                     covGG[i + 1], covGG[i + 2]);
          }
          n2 = ndim;
          break;

        case 2:                 // 2-Tangent
          for (idim = 0; idim < ndim; idim++)
            dd[idim] = x2p[idim] - x1p[idim];
          st_cov(model, 1, dd[0], dd[1], dd[2], covar, covGp, covGG);
          covtab[0] = matrix_UAV(ndim, covGG.data(),
                                 tx1[0], tx1[1], tx1[2],
                                 tx2[0], tx2[1], tx2[2]);
          n2 = 1;
          break;

        case 3:                 // 2-IsoPotential
          for (idim = 0; idim < ndim; idim++)
            dd[idim] = x1p[idim] - x2p[idim];
          st_cov(model, 1, dd[0], dd[1], dd[2], covar, cov2Gp, covGG);
          for (idim = 0; idim < ndim; idim++)
            dd[idim] = x1p[idim] - x20[idim];
          st_cov(model, 1, dd[0], dd[1], dd[2], covar, covGp, covGG);
          covtab[0] = matrix_UV(ndim, tx1[0], tx1[1], tx1[2],
                                cov2Gp[0] - covGp[0], cov2Gp[1] - covGp[1],
                                cov2Gp[2] - covGp[2]);
          n2 = 1;
          break;

        case -3:                 // 2-IsoPotential-Target
          for (idim = 0; idim < ndim; idim++)
            dd[idim] = x1p[idim] - x2p[idim];
          st_cov(model, 1, dd[0], dd[1], dd[2], covar, cov2Gp, covGG);
          covtab[0] = matrix_UV(ndim, tx1[0], tx1[1], tx1[2], cov2Gp[0],
                                cov2Gp[1], cov2Gp[2]);
          n2 = 1;
          break;
      }
      break;

    case 3:                     // 1-IsoPotential
      n1 = 1;
      switch (type2)
      {
        case 1:                 // 2-Gradient
        case -1:                 // 2-Gradient-Target
          for (idim = 0; idim < ndim; idim++)
            dd[idim] = x10[idim] - x2p[idim];
          st_cov(model, 1, dd[0], dd[1], dd[2], covar, covGp, covGG);
          for (idim = 0; idim < ndim; idim++)
            dd[idim] = x1p[idim] - x2p[idim];
          st_cov(model, 1, dd[0], dd[1], dd[2], covar, cov2Gp, covGG);
          for (idim = 0; idim < ndim; idim++)
            covtab[idim] = covGp[idim] - cov2Gp[idim];
          n2 = ndim;
          break;

        case 2:                 // 2-Tangent
          for (idim = 0; idim < ndim; idim++)
            dd[idim] = x1p[idim] - x2p[idim];
          st_cov(model, 1, dd[0], dd[1], dd[2], covar, covGp, covGG);
          for (idim = 0; idim < ndim; idim++)
            dd[idim] = x10[idim] - x2p[idim];
          st_cov(model, 1, dd[0], dd[1], dd[2], covar, cov2Gp, covGG);
          covtab[0] = matrix_UV(ndim, tx2[0], tx2[1], tx2[2],
                                cov2Gp[0] - covGp[0], cov2Gp[1] - covGp[1],
                                cov2Gp[2] - covGp[2]);
          n2 = 1;
          break;

        case 3:                 // 2-IsoPotential
          for (idim = 0; idim < ndim; idim++)
            dd[idim] = x2p[idim] - x1p[idim];
          st_cov(model, 0, dd[0], dd[1], dd[2], covar1, covGp, covGG);
          for (idim = 0; idim < ndim; idim++)
            dd[idim] = x2p[idim] - x10[idim];
          st_cov(model, 0, dd[0], dd[1], dd[2], covar2, covGp, covGG);
          for (idim = 0; idim < ndim; idim++)
            dd[idim] = x20[idim] - x1p[idim];
          st_cov(model, 0, dd[0], dd[1], dd[2], covar3, covGp, covGG);
          for (idim = 0; idim < ndim; idim++)
            dd[idim] = x20[idim] - x10[idim];
          st_cov(model, 0, dd[0], dd[1], dd[2], covar4, covGp, covGG);
          covtab[0] = covar1 - covar2 - covar3 + covar4;
          n2 = 1;
          break;

        case -3:                 // 2-IsoPotential-Target
          for (idim = 0; idim < ndim; idim++)
            dd[idim] = x2p[idim] - x1p[idim];
          st_cov(model, 0, dd[0], dd[1], dd[2], covar1, covGp, covGG);
          for (idim = 0; idim < ndim; idim++)
            dd[idim] = x2p[idim] - x10[idim];
          st_cov(model, 0, dd[0], dd[1], dd[2], covar2, covGp, covGG);
          covtab[0] = covar1 - covar2;
          n2 = 1;
          break;
      }
      break;
  }

  /* Printout (verbose option) */

  if (VERBOSE) print_matrix("Covariance", 0, 1, n2, n1, NULL, covtab.data());

  return (0);
}
