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
#include "Anamorphosis/AnamDiscreteIR.hpp"
#include "Anamorphosis/AnamHermite.hpp"
#include "Basic/File.hpp"
#include "Basic/Law.hpp"
#include "Basic/NamingConvention.hpp"
#include "Basic/OptDbg.hpp"
#include "Basic/String.hpp"
#include "Basic/Utilities.hpp"
#include "Basic/VectorNumT.hpp"
#include "Core/Keypair.hpp"
#include "Covariances/CovContext.hpp"
#include "Db/Db.hpp"
#include "Db/DbGrid.hpp"
#include "Drifts/DriftList.hpp"
#include "Enum/ECalcMember.hpp"
#include "Estimation/KrigingSystem.hpp"
#include "Matrix/AMatrix.hpp"
#include "Matrix/MatrixFactory.hpp"
#include "Matrix/MatrixSquare.hpp"
#include "Matrix/MatrixSymmetric.hpp"
#include "Model/CovInternal.hpp"
#include "Model/Model.hpp"
#include "Neigh/ANeigh.hpp"
#include "Neigh/NeighImage.hpp"
#include "Neigh/NeighMoving.hpp"
#include "Neigh/NeighUnique.hpp"
#include "geoslib_define.h"
#include "geoslib_f.h"
#include "geoslib_f_private.h"
#include "geoslib_old_f.h"
#include <climits>
#include <cmath>
#include <cstdio>
#include <cstring>

/*! \cond */
#define NBYPAS                  5
#define VAR0(iv, jv)            (var0_global[(jv) + nvar * (iv)])
#define LHS_C(i, j)             (lhs[(i) + nred * (j)])
#define RHS_C(i, iv)            (rhs[(i) + nred * (iv)])
#define DISC1(i, idim)          (KOPTION.disc1[(idim) * KOPTION.ntot + (i)])
#define DISC2(i, idim)          (KOPTION.disc2[(idim) * KOPTION.ntot + (i)])
#define LHS_EXP(i, j)           (lhs_global[(i) * neq + (j)])
#define RHS_EXP(i)              (rhs_global[(i)])
#define COV_REF(iz)             (cov_ref[cov_radius + (iz)])
#define SMEAN(i, isimu)         (smean[(isimu) * nfeq + (i)])
#define IAD(ix, iy, iz, nn, ss) (((iz) + nn[2]) + ss[2] * (((iy) + nn[1]) + ss[1] * ((ix) + nn[0])))
#define COV_RES(ix, iy, iz)     cov_res[IAD(ix, iy, iz, cov_nn, cov_ss)]
#define COV_TOT(ix, iy, iz)     cov_tot[IAD(ix, iy, iz, cov_nn, cov_ss)]
#define NUM_TOT(ix, iy, iz)     num_tot[IAD(ix, iy, iz, cov_nn, cov_ss)]
#define NEI_CUR(ix, iy, iz)     nei_cur[IAD(ix, iy, iz, nei_nn, nei_ss)]
#define NEI_REF(ix, iy, iz)     nei_ref[IAD(ix, iy, iz, nei_nn, nei_ss)]
#define CC(iz, jz)              (cc[(jz) * ((jz) + 1) / 2 + (iz)])
#define UTAB(i, j)              (utab[(i) + ndat * (j)])
#define SPART(i, j)             (spart[(i) + npart * (j)])
#define COVSS(is, js)           (covss[(js) + ns * (is)])
#define COVGEN(i1, i2)          (covgen[(i2) + n2 * (i1)])
#define COVPP(ip, jp)           (covpp[(jp) + np * (ip)])
#define COVGP(ig, ip)           (covgp[(ip) + np * (ig)])
#define DISTGEN(i, is)          (distgen[(is) + ns * (i)])
#define DISTPS(ip, is)          (distps[(is) + ns * (ip)])
#define DISTGS(ig, is)          (distgs[(is) + ns * (ig)])
#define PRODGEN(i, is)          (prodgen[(is) + ns * (i)])
#define PRODPS(ip, is)          (prodps[(is) + ns * (ip)])
#define PRODGS(ig, is)          (prodgs[(is) + ns * (ig)])
#define DRFTAB(ip, il)          (drftab[(il) + nbfl * (ip)])
#define YMAT(ip, il)            (ymat[(il) + nbfl * (ip)])
/*! \endcond */

namespace gstlrn
{
// TODO : remove all these static stuffs !
static VectorDouble d1_1_global;
static VectorDouble d1_2_global;
static VectorDouble covaux_global, var0_global;
static VectorDouble d1_global, d1_t_global;
static VectorDouble lhs_global, rhs_global, wgt_global, zam1_global;
static VectorInt flag_global;
static Id KRIGE_INIT = 0;
static Id MODEL_INIT = 0;
static Id IECH_OUT   = -1;
static Id FLAG_COLK, FLAG_SIMU, FLAG_EST, FLAG_STD, FLAG_VARZ, FLAG_PROF;
static Id IPTR_EST, IPTR_STD, IPTR_VARZ, IPTR_NBGH;
static Id* RANK_COLCOK;
static Db *DBIN, *DBOUT;
static Koption KOPTION;
static Id INH_FLAG_VERBOSE = 0;
static Id INH_FLAG_LIMIT   = 1;
static String string;

static CovInternal COVINT;

typedef struct
{
  Id ndtot;
  Id rank1;
  Id rank2;
  Model* model;
  Id nugget_opt;
  Id nostd;
  ECalcMember member;
  Id icov_r;
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
static VectorDouble st_core(Id nli, Id nco)
{
  VectorDouble tab;
  double rsize;
  Id size;

  /* Initialization */

  rsize = static_cast<double>(nli) * static_cast<double>(nco);
  if (rsize < 0 || rsize > INT_MAX)
  {
    messerr("Core allocation problem: Size (%d x %d) too big", nli, nco);
    return (tab);
  }
  size = nli * nco;

  /* Allocation */

  tab.resize(size);
  return (tab);
}

/****************************************************************************/
/*!
 **  Manage the relative position array
 **
 ** \param[in]  neq     Number of kriging equations
 **
 *****************************************************************************/
static VectorInt st_relative_position_array(Id neq)
{
  /* Creation */

  VectorInt rel(neq);
  Id j = 0;
  for (Id i = 0; i < neq; i++)
  {
    // Comment: the next code has been commented out as it seems to be unused
    // The whole function should probably disappear soon.
    // if (flag_global != NULL && flag_global[i])
    //   rel[j++] = i + 1;
    // else
    rel[j++] = i + 1;
  }
  return rel;
}

/****************************************************************************/
/*!
 **  Initialize the static global variables
 **
 ** \param[in]  dbin   input Db structure
 ** \param[in]  dbout  output Db structure
 **
 *****************************************************************************/
static void st_global_init(Db* dbin, Db* dbout)
{
  FLAG_COLK = FLAG_PROF = FLAG_SIMU = 0;
  IPTR_EST = IPTR_STD = IPTR_VARZ = IPTR_NBGH = 0;
  IECH_OUT                                    = 0;
  FLAG_EST = FLAG_STD = FLAG_VARZ = false;

  /* Set the global variables */

  DBIN  = dbin;
  DBOUT = dbout;

  /* Change of support coefficient for DGM */

  COVINT = CovInternal();
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
static double st_get_idim(Id loc_rank, Id idim)
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

///****************************************************************************/
///*!
// **  Calculate the covariance between two samples from two Db
// **
// ** \param[in]  model        Model structure
// ** \param[in]  flag_init    Initialize the array beforehand
// ** \param[in]  nugget_opt   Option for the nugget effect basic structure
// ** \li                       0 : no particular option
// ** \li                       1 : discard the nugget effect
// ** \li                      -1 : only consider the nugget effect
// ** \param[in]  nostd        0 standard; +-1 special; ITEST normalized
// ** \param[in]  member       Member of the Kriging System (ECalcMember)
// ** \param[in]  icov_r       rank of the target covariance or -1 for all
// ** \param[in]  weight       Weight attached to this calculation
// ** \param[in]  rank1        Rank of the first sample
// ** \param[in]  rank2        Rank of the second sample
// **
// ** \param[out] d1loc        Working array
// ** \param[out] covtab_loc   Output covariance array
// **
// *****************************************************************************/
// static void st_cov(Model *model,
//                   Id flag_init,
//                   Id nugget_opt,
//                   Id nostd,
//                   const ECalcMember &member,
//                   Id icov_r,
//                   double weight,
//                   Id rank1,
//                   Id rank2,
//                   VectorDouble d1loc,
//                   double *covtab_loc)
//{
//  DECLARE_UNUSED(nostd);
//  DECLARE_UNUSED(nugget_opt);
//
//  /* Initializations */
//
//  if (rank1 >= 0)
//  {
//    COVINT.setDb1(DBIN);
//    COVINT.setIcas1(1);
//    COVINT.setIech1(rank1);
//  }
//  else
//  {
//    COVINT.setDb1(DBOUT);
//    COVINT.setIcas1(2);
//    COVINT.setIech1(IECH_OUT);
//  }
//
//  if (rank2 >= 0)
//  {
//    COVINT.setDb2(DBIN);
//    COVINT.setIcas2(1);
//    COVINT.setIech2(rank2);
//  }
//  else
//  {
//    COVINT.setDb2(DBOUT);
//    COVINT.setIcas2(2);
//    COVINT.setIech2(IECH_OUT);
//  }
//
//  CovCalcMode mode(member);
//  model->setActiveCovListFromOne(icov_r);
//  model->evaluateMatInPlace(&COVINT, d1loc, covtab_loc, flag_init, weight, &mode);
//}
//
///****************************************************************************/
///*!
// **  Internal recursive function for calculating covariance between data
// **  and data, when data are discretized
// **
// ** \param[in]  idim   Space dimension for current iteration (first point)
// ** \param[in]  jdim   Space dimension for current iteration (second point)
// ** \param[in]  it     Pointer to the Internal Disc_Structure
// **
// *****************************************************************************/
// static void st_data_discretize_dd(Id idim, Id jdim, Disc_Structure *it)
//{
//  double exts2, dsize, decal;
//
//  // Initialization
//
//  if (idim < it->model->getNDim() - 1)
//  {
//    idim = idim + 1;
//
//    // Loop in the current dimension
//
//    exts2 = DBIN->getLocVariable(ELoc::BLEX,it->rank1, idim) / 2.;
//    dsize = KOPTION->dsize[idim];
//
//    if (exts2 <= 0. || dsize <= 0.)
//    {
//
//      /* Punctual support */
//
//      d1_1_global[idim] = 0.;
//      st_data_discretize_dd(idim, jdim, it);
//    }
//    else
//    {
//
//      /* Implicit loop until reaching the edge of the data */
//
//      decal = -exts2 + dsize / 2.;
//      do
//      {
//        d1_1_global[idim] = decal;
//        st_data_discretize_dd(idim, jdim, it);
//        decal = decal + dsize;
//      }
//      while (decal < exts2);
//    }
//  }
//  else if (jdim < it->model->getNDim() - 1)
//  {
//    jdim = jdim + 1;
//
//    // Loop in the current dimension
//
//    exts2 = DBIN->getLocVariable(ELoc::BLEX,it->rank2, jdim) / 2.;
//    dsize = KOPTION->dsize[jdim];
//
//    if (exts2 <= 0 || dsize <= 0.)
//    {
//      /* Punctual support */
//
//      d1_2_global[jdim] = 0.;
//      st_data_discretize_dd(idim, jdim, it);
//    }
//    else
//    {
//
//      /* Implicit loop until reaching the edge of the data */
//
//      decal = -exts2 + dsize / 2.;
//      do
//      {
//        d1_2_global[jdim] = decal;
//        st_data_discretize_dd(idim, jdim, it);
//        decal = decal + dsize;
//      }
//      while (decal < exts2);
//    }
//  }
//  else
//  {
//
//    // End of implicit loop on dimensions
//
//    it->ndtot++;
//    for (Id i = 0; i < it->model->getNDim(); i++)
//      d1_t_global[i] = d1_global[i] + d1_1_global[i] + d1_2_global[i];
//    st_cov(it->model, 0, it->nugget_opt, it->nostd, it->member, it->icov_r,
//           it->weight, it->rank1, it->rank2, d1_t_global, covaux_global);
//  }
//}

void set_DBIN(Db* dbin)
{
  DBIN = dbin;
}

void set_DBOUT(Db* dbout)
{
  DBOUT = dbout;
}

///****************************************************************************/
///*!
// **  Internal recursive function for calculating covariance between data
// **  and target, when data is discretized
// **
// ** \param[in]  idim   Space dimension for current iteration
// ** \param[in]  it     Pointer to the Internal Disc_Structure
// **
// *****************************************************************************/
// static void st_data_discretize_dg(Id idim, Disc_Structure *it)
//{
//  double exts2, dsize, decal;
//
//  // Initialization
//
//  if (idim < it->model->getNDim() - 1)
//  {
//    idim = idim + 1;
//
//    // Loop in the current dimension
//
//    exts2 = DBIN->getLocVariable(ELoc::BLEX,it->rank1, idim) / 2.;
//    dsize = KOPTION->dsize[idim];
//
//    if (exts2 <= 0. || dsize <= 0.)
//    {
//
//      /* Punctual support */
//
//      d1_1_global[idim] = 0.;
//      st_data_discretize_dg(idim, it);
//    }
//    else
//    {
//
//      /* Implicit loop until reaching the edge of the data */
//
//      decal = -exts2 + dsize / 2.;
//      do
//      {
//
//        d1_1_global[idim] = decal;
//        st_data_discretize_dg(idim, it);
//        decal = decal + dsize;
//      }
//      while (decal < exts2);
//    }
//  }
//  else
//  {
//
//    // End of implicit loop on dimensions
//
//    it->ndtot++;
//    for (Id i = 0; i < it->model->getNDim(); i++)
//      d1_t_global[i] = d1_global[i] + d1_1_global[i];
//    st_cov(it->model, 0, it->nugget_opt, it->nostd, it->member, it->icov_r,
//           it->weight, it->rank1, it->rank2, d1_t_global, covaux_global);
//  }
//}

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
static double st_get_ivar(Id rank, Id ivar)
{
  double value;
  Id jvar;

  if (rank >= 0)
  {

    // Variable in the Input file

    if (!FLAG_SIMU)

      // Particular case of simulations

      value = DBIN->getZVariable(rank, ivar);
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
static double st_get_verr(Id rank, Id ivar)
{
  double value;

  if (rank >= 0)
  {
    value = DBIN->getLocVariable(ELoc::V, rank, ivar);
  }
  else
  {
    value = DBOUT->getLocVariable(ELoc::V, IECH_OUT, ivar);
  }
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
 **
 ** \remarks The address of the argument 'neigh' is memorized in a local
 ** \remarks static variable
 **
 *****************************************************************************/
static Id st_check_environment(Id flag_in,
                               Id flag_out,
                               Model* model)
{
  Id error, ndim, nvar, nfex;

  /* Initializations */

  error = 1;
  ndim = nfex = 0;

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
    nvar = model->getNVar();
    if (nvar <= 0)
    {
      messerr("The number of variables must be positive = %d",
              model->getNVar());
      goto label_end;
    }
    // The following test is avoided in the case of simulations
    // as there may be no Z-variable defined as this stage (Gibbs)
    if (flag_in && !FLAG_SIMU && DBIN->getNLoc(ELoc::Z) != nvar)
    {
      messerr("The number of variables of the Data (%d)",
              DBIN->getNLoc(ELoc::Z));
      messerr("does not match the number of variables of the Model (%d)", nvar);
      goto label_end;
    }
    if (model->getNCov() <= 0)
    {
      messerr("The number of covariance must be positive");
      goto label_end;
    }
    if (model->getNDim() <= 0)
    {
      messerr("The Space Dimension must be positive = %d",
              model->getNDim());
      goto label_end;
    }
    if (static_cast<Id>(model->getNDim()) != ndim)
    {
      messerr("The Space Dimension of the Db structure (%d)", ndim);
      messerr("Does not correspond to the Space Dimension of the model (%d)",
              model->getNDim());
      goto label_end;
    }

    // External drifts
    nfex = model->getNExtDrift();
    if (nfex > 0)
    {
      if (flag_out && DBOUT->getNLoc(ELoc::F) != nfex)
      {
        messerr("The Model requires %d external drift(s)", model->getNExtDrift());
        messerr("but the output Db refers to %d external drift variables",
                DBOUT->getNLoc(ELoc::F));
        goto label_end;
      }

      if (flag_in && DBIN->getNLoc(ELoc::F) != nfex)
      {
        if (!(flag_out && DBOUT->isGrid()))
        {
          messerr("The Model requires %d external drift(s)", model->getNExtDrift());
          messerr("but the input Db refers to %d external drift variables",
                  DBIN->getNLoc(ELoc::F));
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
    db_mini.resize(ndim, TEST);
    db_maxi.resize(ndim, TEST);

    /* Input Db structure */

    if (flag_in)
      DBIN->getExtensionInPlace(db_mini, db_maxi, true);

    /* Output Db structure */

    if (flag_out)
      DBOUT->getExtensionInPlace(db_mini, db_maxi, true);

    model->setField(VH::extensionDiagonal(db_mini, db_maxi));
  }

  /* Set the error return code */

  error = 0;
label_end:
  return (error);
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
static Id st_model_manage(Id mode, Model* model)

{
  Id nvar;

  /* Initializations */

  nvar = model->getNVar();

  /* Dispatch */

  if (mode == 1)
  {

    /* Allocation */

    if (MODEL_INIT) return (1);
    Id ndim = DBIN->getNDim();
    d1_global.resize(ndim);
    d1_1_global.resize(ndim);
    d1_2_global.resize(ndim);
    d1_t_global.resize(ndim);
    covaux_global = st_core(nvar, nvar);
    if (covaux_global.empty()) return (1);
    MODEL_INIT = 1;
  }
  else
  {
    if (!MODEL_INIT) return (1);
    covaux_global.clear();
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
static Id st_krige_manage_basic(Id mode,
                                Id nech,
                                Id nmax,
                                Id nvar,
                                Id nfeq)
{
  DECLARE_UNUSED(nech);
  Id neqmax, ncmax;

  /* Initializations */

  ncmax  = nmax * nvar;
  neqmax = ncmax + nfeq;
  if (FLAG_COLK) neqmax += nvar;
  //  if (FLAG_COLK) nech += 1;

  /* Dispatch */

  if (mode == 1)
  {

    /* Allocation */

    if (KRIGE_INIT) return (1);
    flag_global.resize(neqmax);
    if (flag_global.empty()) return (1);
    lhs_global = st_core(neqmax, neqmax);
    if (lhs_global.empty()) return (1);
    rhs_global = st_core(neqmax, nvar);
    if (rhs_global.empty()) return (1);
    zam1_global = st_core(neqmax, 1);
    if (zam1_global.empty()) return (1);
    wgt_global = st_core(neqmax, nvar);
    if (wgt_global.empty()) return (1);
    var0_global = st_core(nvar, nvar);
    if (var0_global.empty()) return (1);
    KRIGE_INIT = 1;
  }
  else
  {

    /* Deallocation */

    if (!KRIGE_INIT) return (1);
    flag_global.clear();
    lhs_global.clear();
    rhs_global.clear();
    zam1_global.clear();
    wgt_global.clear();
    var0_global.clear();
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
 ** \param[in]  neigh ANeigh structure
 **
 *****************************************************************************/
static Id st_get_nmax(ANeigh* neigh)
{
  return neigh->getNSampleMax(DBIN);
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
 ** \param[in]  neigh      ANeigh structure
 **
 ** \remarks  The number of variables corresponds to the number of variables
 ** \remarks  to be calculated. It is not necessarily equal to the number of
 ** \remarks  variables contained in Model (when kriging a linear combination
 ** \remarks  of variables for example): hence the use of the 'nvar' passed
 ** \remarks  as an argument
 **
 *****************************************************************************/
static Id st_krige_manage(Id mode,
                          Id nvar,
                          Model* model,
                          ANeigh* neigh)
{
  Id nech, nfeq, nmax;

  /* Initializations */

  nvar = model->getNVar();
  nfeq = model->getNDriftEquation();
  nech = DBIN->getNSample();
  nmax = st_get_nmax(neigh);

  return (st_krige_manage_basic(mode, nech, nmax, nvar, nfeq));
}

/****************************************************************************/
/*!
 **  Allocate the Target discretization
 **
 ** \param[in]  ndim       Space dimension
 ** \param[in]  ndiscs     Discretization parameters (or NULL)
 **
 *****************************************************************************/
static Id st_block_discretize_alloc(Id ndim, const VectorInt& ndiscs)
{
  Id ntot;

  ntot = 1;
  for (Id idim = 0; idim < ndim; idim++)
    ntot *= ndiscs[idim];
  if (ntot <= 0) return (1);
  KOPTION.ntot = ntot;

  KOPTION.ndisc = ndiscs;
  KOPTION.disc1.resize(ndim * ntot);
  KOPTION.disc2.resize(ndim * ntot);
  return (0);
}

/****************************************************************************/
/*!
 **  Allocate the Data discretization
 **
 ** \param[in] ndim    Space dimension
 **
 *****************************************************************************/
static void st_data_discretize_alloc(Id ndim)

{
  Id nrow, ncol;

  KOPTION.flag_data_disc = 0;
  if (DBIN->getNLoc(ELoc::BLEX) > 0)
  {
    if (!get_keypair("Data_Discretization", &nrow, &ncol, KOPTION.dsize))
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
        KOPTION.flag_data_disc = 1;
      }
    }
    else
    {
      if (DBIN->getNLoc(ELoc::BLEX) > 0)
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
static void st_block_discretize(Id mode, Id flag_rand, Id iech)
{
  Id i, j, jech, ntot, nd, nval, ndim, idim, memo;
  double taille;

  /* Initializations */

  memo = law_get_random_seed();
  ntot = KOPTION.ntot;
  ndim = KOPTION.ndim;
  law_set_random_seed(1234546);
  DbGrid* dbgrid = dynamic_cast<DbGrid*>(DBOUT);

  /* Loop on the discretization points */

  for (i = 0; i < ntot; i++)
  {
    jech = i;
    nval = ntot;
    for (idim = ndim - 1; idim >= 0; idim--)
    {
      taille = (mode == 0) ? dbgrid->getDX(idim) : DBOUT->getLocVariable(ELoc::BLEX, iech, idim);
      nd     = KOPTION.ndisc[idim];
      nval /= nd;
      j = jech / nval;
      jech -= j * nval;
      DISC1(i, idim) = taille * ((j + 0.5) / nd - 0.5);
      DISC2(i, idim) = DISC1(i, idim);
      if (flag_rand)
        DISC2(i, idim) += taille * law_uniform(-0.5, 0.5) / static_cast<double>(nd);
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
 ** \param[in]  ndiscs      Discretization parameters (or NULL)
 **
 ** \remark  This function manages the global structure KOPTION
 **
 *****************************************************************************/
Id krige_koption_manage(Id mode,
                        Id flag_check,
                        const EKrigOpt& calcul,
                        Id flag_rand,
                        const VectorInt& ndiscs)
{
  Id ndim, error;

  /* Initializations */

  error = 1;
  ndim  = DBOUT->getNDim();

  /* Dispatch */

  if (mode == 1)
  {

    /* Allocation of the structure */
    KOPTION.calcul = calcul;

    // Target discretization
    KOPTION.ndim = ndim;
    KOPTION.ntot = 0;
    KOPTION.disc1.clear();
    KOPTION.disc2.clear();
    KOPTION.ndisc.clear();

    // Data discretization
    KOPTION.flag_data_disc = 0;
    KOPTION.dsize.clear();

    /* Data discretization case (optional) */

    st_data_discretize_alloc(ndim);

    /* Block discretization case */

    switch (KOPTION.calcul.toEnum())
    {
      case EKrigOpt::E_POINT:
      case EKrigOpt::E_DRIFT:
      case EKrigOpt::E_DGM:
        break;

      case EKrigOpt::E_BLOCK:

        /* Preliminary checks */

        if (flag_check && !DBOUT->isGrid())
        {
          messerr("Discretization is not allowed if the Target is not a Grid");
          return error;
        }
        if (ndiscs.empty())
        {
          messerr("For block estimation, Discretization must be provided");
          return error;
        }

        if (st_block_discretize_alloc(ndim, ndiscs)) return error;

        st_block_discretize(0, flag_rand, 0);

        break;
    }
    error = 0;
  }
  else
  {
    error = 0;
  }
  return (error);
}

/****************************************************************************/
/*!
 **  Print the L.H.S. matrix
 **
 ** \param[in]  nech    Number of active points (optional)
 ** \param[in]  neq     Number of equations
 ** \param[in]  nred    Reduced number of equations
 ** \param[in]  flag    Flag array (optional)
 ** \param[in]  lhs     Kriging L.H.S
 **
 *****************************************************************************/
void krige_lhs_print(Id nech,
                     Id neq,
                     Id nred,
                     const Id* flag,
                     const double* lhs)
{
  Id i, j, ipass, npass, ideb, ifin;

  /* Initializations */

  VectorInt rel = st_relative_position_array(neq);
  npass         = (nred - 1) / NBYPAS + 1;

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

    if (flag != nullptr)
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
void krige_rhs_print(Id nvar,
                     Id nech,
                     Id neq,
                     Id nred,
                     const Id* flag,
                     double* rhs)
{
  VectorInt rel = st_relative_position_array(neq);

  /* General Header */

  mestitle(0, "RHS of Kriging matrix (compressed)");
  if (nech > 0) message("Number of active samples    = %d\n", nech);
  message("Total number of equations   = %d\n", neq);
  message("Reduced number of equations = %d\n", nred);
  message("Number of right-hand sides  = %d\n", nvar);

  /* Kriging option */

  {
    switch (KOPTION.calcul.toEnum())
    {
      case EKrigOpt::E_POINT:
        message("Punctual Estimation\n");
        break;

      case EKrigOpt::E_BLOCK:
        message("Block Estimation : Discretization = ");
        for (Id idim = 0; idim < KOPTION.ndim; idim++)
        {
          if (idim != 0) message(" x ");
          message("%ld", KOPTION.ndisc[idim]);
        }
        message("\n");
        break;

      case EKrigOpt::E_DRIFT:
        message("Drift Estimation\n");
        break;

      case EKrigOpt::E_DGM:
        message("DGM Estimation\n");
        break;
    }
  }
  message("\n");

  /* Header line */

  tab_prints(NULL, "Rank");
  if (flag != nullptr) tab_prints(NULL, "Flag");
  for (Id ivar = 0; ivar < nvar; ivar++)
    tab_printi(NULL, ivar + 1);
  message("\n");

  /* Matrix lines */

  for (Id i = 0; i < nred; i++)
  {
    tab_printi(NULL, i + 1);
    if (flag != nullptr) tab_printi(NULL, rel[i]);
    for (Id ivar = 0; ivar < nvar; ivar++)
      tab_printg(NULL, RHS_C(i, ivar));
    message("\n");
  }
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
static void st_krige_wgt_print(Id status,
                               Id nvar,
                               Id nvar_m,
                               Id nfeq,
                               const VectorInt& nbgh_ranks,
                               Id nred,
                               Id icase,
                               const Id* flag,
                               const double* wgt)
{
  double value;
  Id iwgt, ivar, jvar_m, ivar_m, iech, lec, cumflag, idim, ndim, ib, number,
    flag_value, nech;

  /* Initializations */

  nech     = static_cast<Id>(nbgh_ranks.size());
  ndim     = DBIN->getNDim();
  auto sum = st_core(nvar_m, 1);
  if (sum.empty()) return;

  /* Header */

  mestitle(0, "(Co-) Kriging weights");

  /* First line */

  tab_prints(NULL, "Rank");
  for (idim = 0; idim < ndim; idim++)
  {
    String strloc = getLocatorName(ELoc::X, idim);
    tab_prints(NULL, strloc.c_str());
  }
  if (DBIN->hasLocVariable(ELoc::C)) tab_prints(NULL, "Code");
  if (DBIN->getNLoc(ELoc::V) > 0)
    tab_prints(NULL, "Err.");
  if (KOPTION.flag_data_disc)
    for (idim = 0; idim < ndim; idim++)
    {
      (void)gslSPrintf(string, "Size%d", idim + 1);
      tab_prints(NULL, string.data());
    }
  tab_prints(NULL, "Data");
  for (ivar = 0; ivar < nvar; ivar++)
  {
    (void)gslSPrintf(string, "Z%d*", ivar + 1);
    tab_prints(NULL, string.data());
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
      flag_value = (flag != nullptr) ? flag[lec] : 1;
      tab_printi(NULL, iech + 1);
      for (idim = 0; idim < ndim; idim++)
        tab_printg(NULL, st_get_idim(nbgh_ranks[iech], idim));
      if (DBIN->hasLocVariable(ELoc::C))
        tab_printg(NULL, DBIN->getLocVariable(ELoc::C, nbgh_ranks[iech], 0));
      if (DBIN->getNLoc(ELoc::V) > 0)
        tab_printg(NULL, st_get_verr(nbgh_ranks[iech], (FLAG_PROF) ? 0 : jvar_m));
      if (KOPTION.flag_data_disc)
      {
        for (idim = 0; idim < ndim; idim++)
          tab_printg(NULL, DBIN->getLocVariable(ELoc::BLEX, nbgh_ranks[iech], idim));
      }
      if (icase < 0)
        tab_printg(NULL, st_get_ivar(nbgh_ranks[iech], jvar_m));
      else
        tab_prints(NULL, "   ");

      for (ivar = 0; ivar < nvar; ivar++)
      {
        iwgt  = nred * ivar + cumflag;
        value = (wgt != nullptr && status == 0 && flag_value) ? wgt[iwgt] : TEST;
        if (!FFFF(value)) sum[ivar] += value;
        tab_printg(NULL, value);
      }
      if (flag_value) cumflag++;
      message("\n");
    }

    number = 1 + ndim + 1;
    if (DBIN->getNLoc(ELoc::V) > 0) number++;
    if (KOPTION.flag_data_disc) number += ndim + 1;
    tab_prints(NULL, "Sum of weights", number, EJustify::LEFT);
    for (ivar = 0; ivar < nvar; ivar++)
    {
      value = (status == 0) ? sum[ivar] : TEST;
      tab_printg(NULL, value);
    }
    message("\n");
  }

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
    value = (status == 0) ? zam1_global[iwgt] : TEST;
    tab_printg(NULL, value);

    message("\n");
  }
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
static void st_result_kriging_print(Id flag_xvalid, Id nvar, Id status)
{
  Id ivar;
  double value;

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
      message("Printout for Cross-validation should not be performed anymore\n");
    }
    else
    {
      message("Variable Z%-2d\n", ivar + 1);
      if (FLAG_EST)
      {
        value = (status == 0) ? DBOUT->getArray(IECH_OUT, IPTR_EST + ivar) : TEST;
        tab_printg(" - Estimate  = ", value);
        message("\n");
      }
      if (FLAG_STD)
      {
        value = (status == 0) ? DBOUT->getArray(IECH_OUT, IPTR_STD + ivar) : TEST;
        tab_printg(" - Std. Dev. = ", value);
        value = (status == 0) ? VAR0(ivar, ivar) : TEST;
        message("\n");
        tab_printg(" - Cov(h=0)  = ", value);
        message("\n");
      }
      if (FLAG_VARZ)
      {
        value = (status == 0) ? DBOUT->getArray(IECH_OUT, IPTR_VARZ + ivar) : TEST;
        tab_printg(" - Var(Z*)   = ", value);
        message("\n");
      }
    }
  }
}

/****************************************************************************/
/*!
 **  Conditioning Kriging
 **
 ** \return  Error return code
 **
 ** \param[in]  dbin       input Db structure
 ** \param[in]  dbout      output Db structure
 ** \param[in]  model      Model structure
 ** \param[in]  neigh      ANeigh structure
 ** \param[in]  flag_bayes 1 if Bayes option is switched ON
 ** \param[in]  dmean      Array giving the prior means for the drift terms
 ** \param[in]  dcov       Array containing the prior covariance matrix
 **                        for the drift terms
 ** \param[in]  icase      Case for PGS and GRF (or -1)
 ** \param[in]  nbsimu     Number of simulations
 ** \param[in]  flag_dgm   1 if the DGM version of kriging should be used
 **
 ** \remark: The model contains an anamorphosis with a change of support
 ** \remark: coefficient as soon as flag_dgm is TRUE
 **
 *****************************************************************************/
Id _krigsim(Db* dbin,
            Db* dbout,
            const Model* model,
            ANeigh* neigh,
            bool flag_bayes,
            const VectorDouble& dmean,
            const MatrixSymmetric& dcov,
            Id icase,
            Id nbsimu,
            bool flag_dgm)
{
  // Preliminary checks

  if (neigh->getType() == ENeigh::IMAGE)
  {
    messerr("This tool cannot function with an IMAGE neighborhood");
    return 1;
  }

  /* Add the attributes for storing the results */

  Id iptr_est = dbout->getColIdxByLocator(ELoc::SIMU, 0);
  if (iptr_est < 0) return 1;

  /* Setting options */

  KrigOpt krigopt;
  krigopt.setOptionDGM(flag_dgm);

  KrigingSystem ksys(dbin, dbout, model, neigh, krigopt);
  if (ksys.setKrigOptFlagSimu(true, nbsimu, icase)) return 1;
  if (ksys.updKrigOptEstim(iptr_est, -1, -1, true)) return 1;
  if (ksys.setKrigOptBayes(flag_bayes, dmean, dcov)) return 1;
  if (!ksys.isReady()) return 1;

  /* Loop on the targets to be processed */

  for (Id iech_out = 0; iech_out < dbout->getNSample(); iech_out++)
  {
    mes_process("Conditional Simulation", dbout->getNSample(), iech_out);
    if (ksys.estimate(iech_out)) return 1;
  }

  ksys.conclusion();

  return 0;
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
Id global_transitive(DbGrid* dbgrid,
                     Model* model,
                     Id flag_verbose,
                     Id flag_regular,
                     Id ndisc,
                     double* abundance,
                     double* sse,
                     double* cvtrans)
{
  Id i, ix, iy, ix1, ix2, iy1, iy2, nx, ny, flag_value;
  double c00, cvv, dx, dy, dsum, gint, dsse, wtot, value;
  CovCalcMode mode;

  /* Initializations */

  cvv = wtot = dsse = gint = dsum = 0.;
  st_global_init(dbgrid, dbgrid);
  if (st_check_environment(0, 1, model)) return 1;
  ;
  Id ndim = dbgrid->getNDim();
  VectorDouble d1(ndim, 0.);

  if (ndim < 1 || ndim > 2)
  {
    messerr("The transitive global estimation is implemented for 1 and 2 space only");
    return 1;
  }
  if (model->getNVar() != 1)
  {
    messerr("The transitive global estimation is implemented for 1 variable only");
    return 1;
  }

  /* Core allocation */

  c00 = model->evaluateOneGeneric(nullptr, d1);

  /* Abundance estimation */

  flag_value = 0;
  if (dbgrid->getNLoc(ELoc::Z) == 1)
  {
    for (i = 0; i < dbgrid->getNSample(); i++)
    {
      value = dbgrid->getZVariable(i, 0);
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

      dsse = 0.;
      for (ix = -nx + 1; ix <= nx; ix++)
        for (iy = -ny + 1; iy <= ny; iy++)
        {
          d1[0] = dx * ix;
          d1[1] = dy * iy;
          dsse += model->evaluateOneGeneric(nullptr, d1);
        }
      dsse *= dx * dy;
      // TODO : appeler model_integral
      // if (model_integral(model,ndisc,&gint)) goto label_end;
      *sse = dsse - gint;
    }
    else
    {

      /* Stratified case */

      cvv = 0.;
      for (ix1 = 0; ix1 < ndisc; ix1++)
        for (iy1 = 0; iy1 < ndisc; iy1++)
          for (ix2 = 0; ix2 < ndisc; ix2++)
            for (iy2 = 0; iy2 < ndisc; iy2++)
            {
              d1[0] = dx * (ix2 - ix1) / ndisc;
              d1[1] = dy * (iy2 - iy1) / ndisc;
              cvv += model->evaluateOneGeneric(nullptr, d1);
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

      dsse = 0.;
      for (ix = -nx + 1; ix <= nx; ix++)
      {
        d1[0] = dx * ix;
        dsse += model->evaluateOneGeneric(nullptr, d1);
      }
      dsse /= dx;
      // TODO: appeler model_integral
      // if (model_integral(model,ndisc,&gint)) goto label_end;
      *sse = dsse - gint;
    }
    else
    {

      /* Stratified case */

      cvv = 0.;
      for (ix1 = 0; ix1 < ndisc; ix1++)
        for (ix2 = 0; ix2 < ndisc; ix2++)
        {
          d1[0] = dx * (ix2 - ix1) / ndisc;
          cvv += model->evaluateOneGeneric(nullptr, d1);
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
  return 0;
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
static Id st_get_limits(DbGrid* db, double top, double bot, Id* ideb, Id* ifin)
{
  Id ndim, nz, iad;
  double z0, dz;

  /* Initializations */

  ndim = db->getNDim();
  z0   = db->getX0(ndim - 1);
  nz   = db->getNX(ndim - 1);
  dz   = db->getDX(ndim - 1);

  /* Preliminary checks */

  if (!FFFF(bot) && !FFFF(top) && top < bot)
  {
    messerr("Error: Top(%lf) must be larger than Bottom (%lf)", top, bot);
    return (1);
  }

  if (FFFF(bot))
    iad = 0;
  else
    iad = static_cast<Id>((bot - z0) / dz);
  if (bot > z0 + iad * dz) iad++;
  *ideb = MAX(0, MIN(iad, nz - 1));

  if (FFFF(top))
    iad = nz - 1;
  else
    iad = static_cast<Id>((top - z0) / dz);
  *ifin = MAX(0, MIN(iad, nz - 1));
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
static Id st_get_neigh(Id ideb,
                       Id ifin,
                       Id neigh_radius,
                       Id* status,
                       Id* nbefore,
                       Id* nafter)
{
  Id iad;

  *status = 1;
  if (IECH_OUT < ideb || IECH_OUT > ifin) return (1);

  iad      = MAX(IECH_OUT - neigh_radius, ideb);
  *nbefore = IECH_OUT - iad;

  iad     = MIN(IECH_OUT + neigh_radius, ifin);
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
static double st_cov_exp(Id dist, const double* cov, Id cov_radius, Id flag_sym)
{
  double val1, val2, val;

  if (flag_sym)
  {
    val1 = cov[cov_radius - dist];
    val2 = cov[cov_radius + dist];
    val  = (val1 + val2) / 2.;
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
static void st_lhs_exp(double* covdd,
                       Id cov_radius,
                       Id flag_sym,
                       Id nfeq,
                       Id nbefore,
                       Id nafter,
                       Id neq)
{
  Id i, j;

  /* Covariance part */

  for (i = -nbefore; i <= nafter; i++)
    for (j = -nbefore; j <= nafter; j++)
    {
      LHS_EXP(i + nbefore, j + nbefore) = st_cov_exp(i - j, covdd, cov_radius,
                                                     flag_sym);
      LHS_EXP(j + nbefore, i + nbefore) = st_cov_exp(j - i, covdd, cov_radius,
                                                     flag_sym);
    }

  /* Drift part */

  if (nfeq == 0) return;
  for (i = -nbefore; i <= nafter; i++)
  {
    LHS_EXP(i + nbefore, neq - 1) = 1.;
    LHS_EXP(neq - 1, i + nbefore) = 1.;
  }
  LHS_EXP(neq - 1, neq - 1) = 0.;
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
static void st_rhs_exp(double* covd0,
                       Id cov_radius,
                       Id flag_sym,
                       Id nfeq,
                       Id nbefore,
                       Id nafter,
                       Id neq)
{
  Id i;

  /* Covariance part */

  for (i = -nbefore; i <= nafter; i++)
    RHS_EXP(i + nbefore) = st_cov_exp(i, covd0, cov_radius, flag_sym);

  /* Drift part */

  if (nfeq == 0) return;
  RHS_EXP(neq - 1) = 1.;
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
static double st_estim_exp(Db* db, const double* wgt, Id nbefore, Id nafter)
{
  Id i;
  double result;

  /* Perform the estimation */

  result = 0.;
  for (i = -nbefore; i <= nafter; i++)
    result += wgt[i + nbefore] * db->getZVariable(IECH_OUT + i, 0);

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
Id anakexp_f(DbGrid* db,
             double* covdd,
             double* covd0,
             double top,
             double bot,
             Id cov_radius,
             Id neigh_radius,
             Id flag_sym,
             Id nfeq)
{
  Id i, ndim, nvarin, nech, size, error, ideb, ifin, neq, status;
  Id nbefore, nafter, nbefore_mem, nafter_mem;
  double result;
  VectorInt ranks;

  /* Initializations */

  error = 1;
  st_global_init(db, db);
  FLAG_EST    = true;
  ndim        = db->getNDim();
  nvarin      = db->getNLoc(ELoc::Z);
  nbefore_mem = nafter_mem = -1;
  size                     = 0;

  /* Prepare the Koption structure */

  if (krige_koption_manage(1, 1, EKrigOpt::POINT, 1, VectorInt())) return (1);

  /* Preliminary checks */

  if (ndim != 1 || !db->isGrid())
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
    ranks[i]       = i;
    flag_global[i] = 1;
  }

  /* Get the limits of the area to be processed */

  if (st_get_limits(db, top, bot, &ideb, &ifin)) goto label_end;

  /* Loop on the grid nodes */

  status = 0;
  for (IECH_OUT = 0; IECH_OUT < nech; IECH_OUT++)
  {
    mes_process("Factorial Kriging Analysis", nech, IECH_OUT);
    OptDbg::setCurrentIndex(IECH_OUT + 1);
    if (!db->isActive(IECH_OUT)) continue;
    if (OptDbg::query(EDbg::KRIGING) || OptDbg::query(EDbg::NBGH) || OptDbg::query(EDbg::RESULTS))
    {
      mestitle(1, "Target location");
      db_sample_print(db, IECH_OUT, 1, 0, 0, 0);
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
      nafter_mem  = nafter;

      /* Establish the L.H.S. of the kriging system */

      st_lhs_exp(covdd, cov_radius, flag_sym, nfeq, nbefore, nafter, neq);
      if (OptDbg::query(EDbg::KRIGING))
        krige_lhs_print(nech, neq, neq, flag_global.data(), lhs_global.data());

      /* Invert the kriging system */

      if (matrix_invert(lhs_global.data(), neq, IECH_OUT))
      {
        status = 1;
        continue;
      }

      /* Establish the R.H.S. of the kriging system */

      st_rhs_exp(covd0, cov_radius, flag_sym, nfeq, nbefore, nafter, neq);
      if (OptDbg::query(EDbg::KRIGING))
        krige_rhs_print(nvarin, nech, neq, neq, flag_global.data(), rhs_global.data());

      /* Derive the kriging weights */

      matrix_product_safe(neq, neq, 1, lhs_global.data(), rhs_global.data(), wgt_global.data());
    }

    /* Calculate the estimation */

    result = st_estim_exp(db, wgt_global.data(), nbefore, nafter);
    DBOUT->setArray(IECH_OUT, IPTR_EST, result);
    if (OptDbg::query(EDbg::RESULTS)) st_result_kriging_print(0, nvarin, status);
  }

  /* Set the error return flag */

  error = 0;

label_end:
  OptDbg::setCurrentIndex(0);
  (void)krige_koption_manage(-1, 1, EKrigOpt::POINT, 1, VectorInt());
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
static void st_calculate_covres(DbGrid* db,
                                Model* model,
                                const double* cov_ref,
                                Id cov_radius,
                                Id flag_sym,
                                const Id cov_ss[3],
                                const Id cov_nn[3],
                                double* cov_res)
{
  double covtab, covver;

  /* Initializations */

  VectorDouble d1(3, 0.);
  double dx     = db->getDX(0);
  double dy     = db->getDX(1);
  double covtot = COV_REF(0);
  double c00    = model->evaluateOneGeneric(nullptr, d1);

  /* Evaluate the array of experimental covariance of the residual variable */

  for (Id ix = -cov_nn[0]; ix <= cov_nn[0]; ix++)
    for (Id iy = -cov_nn[1]; iy <= cov_nn[1]; iy++)
      for (Id iz = -cov_nn[2]; iz <= cov_nn[2]; iz++)
      {
        if (!flag_sym)
          covver = COV_REF(iz);
        else
          covver = (COV_REF(iz) + COV_REF(-iz)) / 2.;
        d1[0]               = dx * ix;
        d1[1]               = dy * iy;
        covtab              = model->evaluateOneGeneric(nullptr, d1);
        COV_RES(ix, iy, iz) = covver * (covtab + covtot - c00) / covtot;
      }
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
static void st_calculate_covtot(DbGrid* db,
                                Id ix0,
                                Id iy0,
                                Id flag_sym,
                                const Id cov_ss[3],
                                const Id cov_nn[3],
                                Id* num_tot,
                                double* cov_tot)
{
  Id ix, iy, iz, ix1, iy1, iz1, jx1, jy1, jz1, jx2, jy2, jz2;
  Id idx, idy, idz, jdx, jdy, iad, jad;
  double val1, val2, val, ratio;

  /* Initialization */

  Id ndim = db->getNDim();
  VectorInt indg(ndim, 0);
  for (ix = -cov_nn[0]; ix <= cov_nn[0]; ix++)
    for (iy = -cov_nn[1]; iy <= cov_nn[1]; iy++)
      for (iz = -cov_nn[2]; iz <= cov_nn[2]; iz++)
      {
        COV_TOT(ix, iy, iz) = 0.;
        NUM_TOT(ix, iy, iz) = 0;
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
        iad     = db->indiceToRank(indg);
        if (!db->isActive(iad)) continue;
        val1 = db->getZVariable(iad, 0);
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
              jad     = db->indiceToRank(indg);
              if (!db->isActive(jad)) continue;
              val2 = db->getZVariable(jad, 0);
              if (FFFF(val2)) continue;

              /* Update the Covariance */

              COV_TOT(idx, idy, idz) += val1 * val2;
              NUM_TOT(idx, idy, idz) += 1;
            }
      }

  /* Scaling */

  ratio = NUM_TOT(0, 0, 0);
  for (ix = -cov_nn[0]; ix <= cov_nn[0]; ix++)
    for (iy = -cov_nn[1]; iy <= cov_nn[1]; iy++)
      for (iz = -cov_nn[2]; iz <= cov_nn[2]; iz++)
      {
        if (NUM_TOT(ix, iy, iz) <= 0.)
        {
          COV_TOT(ix, iy, iz) = TEST;
        }
        else
        {
          COV_TOT(ix, iy, iz) /= ratio;
        }
      }

  /* Symmetry */

  for (ix = -cov_nn[0]; ix < 0; ix++)
    for (iy = -cov_nn[1]; iy <= cov_nn[1]; iy++)
      for (iz = -cov_nn[2]; iz <= cov_nn[2]; iz++)
      {
        val1                = COV_TOT(ix, iy, iz);
        val2                = COV_TOT(-ix, iy, iz);
        val                 = (FFFF(val1) || FFFF(val2)) ? TEST : (val1 + val2) / 2.;
        COV_TOT(ix, iy, iz) = COV_TOT(-ix, iy, iz) = val;
      }

  for (ix = -cov_nn[0]; ix <= cov_nn[0]; ix++)
    for (iy = -cov_nn[1]; iy < 0; iy++)
      for (iz = -cov_nn[2]; iz <= cov_nn[2]; iz++)
      {
        val1                = COV_TOT(ix, -iy, iz);
        val2                = COV_TOT(ix, iy, iz);
        val                 = (FFFF(val1) || FFFF(val2)) ? TEST : (val1 + val2) / 2.;
        COV_TOT(ix, iy, iz) = COV_TOT(ix, -iy, iz) = val;
      }

  if (flag_sym)
    for (ix = -cov_nn[0]; ix <= cov_nn[0]; ix++)
      for (iy = -cov_nn[1]; iy <= cov_nn[1]; iy++)
        for (iz = -cov_nn[2]; iz < 0; iz++)
        {
          val1                 = COV_TOT(ix, iy, -iz);
          val2                 = COV_TOT(ix, iy, iz);
          val                  = (FFFF(val1) || FFFF(val2)) ? TEST : (val1 + val2) / 2.;
          COV_TOT(ix, iy, -iz) = COV_TOT(ix, iy, iz) = val;
        }
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
static VectorInt st_neigh_find(DbGrid* db,
                               Id ix0,
                               Id iy0,
                               Id iz0,
                               const Id nei_ss[3],
                               const Id nei_nn[3],
                               Id* nei_cur)
{
  Id ix, iy, iz, jx, jy, jz, number, locrank;
  VectorInt nbgh_ranks;
  Id ndim = db->getNDim();
  VectorInt indg(ndim, 0);

  /* Loop on the pixels of the neighborhood */

  number = 0;
  for (ix = -nei_nn[0]; ix <= nei_nn[0]; ix++)
    for (iy = -nei_nn[1]; iy <= nei_nn[1]; iy++)
      for (iz = -nei_nn[2]; iz <= nei_nn[2]; iz++)
      {
        NEI_CUR(ix, iy, iz) = -1;
        jx                  = ix0 + ix;
        if (jx < 0 || jx >= db->getNX(0)) continue;
        jy = iy0 + iy;
        if (jy < 0 || jy >= db->getNX(1)) continue;
        jz = iz0 + iz;
        if (jz < 0 || jz >= db->getNX(2)) continue;
        indg[0] = jx;
        indg[1] = jy;
        indg[2] = jz;
        locrank = db->indiceToRank(indg);
        if (FFFF(db->getZVariable(locrank, 0))) continue;
        NEI_CUR(ix, iy, iz) = locrank;
        nbgh_ranks.push_back(locrank);
        flag_global[number] = 1;
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
static Id st_neigh_diff(const Id nei_ss[3],
                        const Id nei_nn[3],
                        Id* nei_ref,
                        const Id* nei_cur)
{
  Id ix, iy, iz, flag1, flag2, flag_diff;

  /* Loop on the pixels of the neighborhood */

  flag_diff = 1;
  for (ix = -nei_nn[0]; ix <= nei_nn[0]; ix++)
    for (iy = -nei_nn[1]; iy <= nei_nn[1]; iy++)
      for (iz = -nei_nn[2]; iz <= nei_nn[2]; iz++)
      {
        flag1 = NEI_REF(ix, iy, iz) < 0;
        flag2 = NEI_CUR(ix, iy, iz) < 0;
        if (flag1 != flag2) goto label_end;
      }
  flag_diff = 0;

label_end:

  /* Copy the current neighborhood into the reference neighborhood */

  for (ix = -nei_nn[0]; ix <= nei_nn[0]; ix++)
    for (iy = -nei_nn[1]; iy <= nei_nn[1]; iy++)
      for (iz = -nei_nn[2]; iz <= nei_nn[2]; iz++)
        NEI_REF(ix, iy, iz) = NEI_CUR(ix, iy, iz);

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
static void st_lhs_exp_3D(Id nech,
                          Id nfeq,
                          const Id nei_ss[3],
                          const Id nei_nn[3],
                          const Id cov_ss[3],
                          const Id cov_nn[3],
                          const Id* nei_cur,
                          const double* cov_tot,
                          double nugget)
{
  Id ix, iy, iz, jx, jy, jz, i, j, neq;
  double value;

  /* Initializations */

  neq = nech + nfeq;

  /* Covariance part of the L.H.S. */

  i = 0;
  for (ix = -nei_nn[0]; ix <= nei_nn[0]; ix++)
    for (iy = -nei_nn[1]; iy <= nei_nn[1]; iy++)
      for (iz = -nei_nn[2]; iz <= nei_nn[2]; iz++)
      {
        if (NEI_CUR(ix, iy, iz) < 0) continue;

        j = 0;
        for (jx = -nei_nn[0]; jx <= nei_nn[0]; jx++)
          for (jy = -nei_nn[1]; jy <= nei_nn[1]; jy++)
            for (jz = -nei_nn[2]; jz <= nei_nn[2]; jz++)
            {
              if (NEI_CUR(jx, jy, jz) < 0) continue;
              value         = COV_TOT(ix - jx, iy - jy, iz - jz);
              LHS_EXP(i, j) = LHS_EXP(j, i) = value;
              if (i == j) LHS_EXP(i, j) += nugget;
              j++;
            }
        i++;
      }

  /* Drift part */

  if (nfeq == 0) return;
  for (i = 0; i < nech; i++)
  {
    LHS_EXP(i, neq - 1) = 1.;
    LHS_EXP(neq - 1, i) = 1.;
  }
  LHS_EXP(neq - 1, neq - 1) = 0.;
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
static void st_rhs_exp_3D(Id nech,
                          Id nfeq,
                          const Id nei_ss[3],
                          const Id nei_nn[3],
                          const Id cov_ss[3],
                          const Id cov_nn[3],
                          const Id* nei_cur,
                          const double* cov_res)
{
  Id ix, iy, iz, neq, i;

  /* Initializations */

  neq = nech + nfeq;

  /* Covariance part of the R.H.S. */

  i = 0;
  for (ix = -nei_nn[0]; ix <= nei_nn[0]; ix++)
    for (iy = -nei_nn[1]; iy <= nei_nn[1]; iy++)
      for (iz = -nei_nn[2]; iz <= nei_nn[2]; iz++)
      {
        if (NEI_CUR(ix, iy, iz) < 0) continue;
        RHS_EXP(i) = COV_RES(ix, iy, iz);
        i++;
      }

  /* Drift part */

  if (nfeq == 0) return;
  RHS_EXP(neq - 1) = 1.;
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
static double st_estim_exp_3D(Db* db,
                              const Id nei_ss[3],
                              const Id nei_nn[3],
                              Id* nei_cur,
                              const double* weight)
{
  Id i, ix, iy, iz;
  double result;

  /* Initializations */

  result = 0.;
  i      = 0;
  for (ix = -nei_nn[0]; ix <= nei_nn[0]; ix++)
    for (iy = -nei_nn[1]; iy <= nei_nn[1]; iy++)
      for (iz = -nei_nn[2]; iz <= nei_nn[2]; iz++)
      {
        if (NEI_CUR(ix, iy, iz) < 0) continue;
        result += weight[i] * db->getZVariable(NEI_CUR(ix, iy, iz), 0);
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
static void st_vario_dump(FILE* file,
                          Id ix0,
                          Id iy0,
                          const Id cov_ss[3],
                          const Id cov_nn[3],
                          const Id* num_tot,
                          const double* cov_tot)
{
  Id ix, iy, iz, num;
  double cov;

  fprintf(file, "*%3ld %3ld\n", ix0, iy0);

  for (ix = -cov_nn[0]; ix <= cov_nn[0]; ix++)
    for (iy = -cov_nn[1]; iy <= cov_nn[1]; iy++)
      for (iz = -cov_nn[2]; iz <= cov_nn[2]; iz++)
      {
        num = (num_tot == nullptr) ? 0 : NUM_TOT(ix, iy, iz);
        cov = COV_TOT(ix, iy, iz);
        fprintf(file, "%3ld %3ld %3ld %3ld %lf\n", ix, iy, iz, num, cov);
      }
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
Id anakexp_3D(DbGrid* db,
              double* cov_ref,
              Id cov_radius,
              Id neigh_ver,
              Id neigh_hor,
              Id flag_sym,
              Model* model,
              double nugget,
              Id nfeq,
              Id dbg_ix,
              Id dbg_iy)
{
  Id i, ix, iy, iz, ndim, nvarin, nech, error, neq, status, ecr;
  Id size_cov, size_nei, flag_new, flag_col;
  Id cov_ss[3], cov_nn[3], nei_ss[3], nei_nn[3];
  VectorInt num_tot, nei_cur, nei_ref;
  VectorDouble cov_tot, cov_res;
  double result;
  FILE* fildmp;
  VectorInt nbgh_ranks;

  /* Initializations */

  error = 1;
  st_global_init(db, db);
  FLAG_EST = true;
  fildmp   = nullptr;
  ndim     = db->getNDim();
  nvarin   = db->getNLoc(ELoc::Z);
  size_nei = 0;
  VectorInt indg(ndim, 0);

  /* Prepare the Koption structure */

  if (krige_koption_manage(1, 1, EKrigOpt::POINT, 1, VectorInt())) return (1);

  /* Preliminary checks */

  if (ndim != 3 || !db->isGrid())
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

  if (dbg_ix >= -1 && dbg_ix < db->getNX(0) && dbg_iy >= -1 && dbg_iy < db->getNX(1))
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

  num_tot.resize(size_cov);
  nei_cur.resize(size_nei);
  nei_ref.resize(size_nei);
  cov_tot = st_core(size_cov, 1);
  if (cov_tot.empty()) goto label_end;
  cov_res = st_core(size_cov, 1);
  if (cov_res.empty()) goto label_end;
  st_krige_manage_basic(1, size_nei, size_nei, 1, nfeq);
  for (i = 0; i < size_nei; i++)
    nei_ref[i] = -1;

  /* Calculate the discretized covariance of residual variable */

  st_calculate_covres(db, model, cov_ref, cov_radius, flag_sym, cov_ss, cov_nn,
                      cov_res.data());
  if (dbg_ix == -1 && dbg_iy == -1)
    st_vario_dump(fildmp, -1, -1, cov_ss, cov_nn, nullptr, cov_res.data());

  /* Loop on the grid nodes */

  status = nech = neq = 0;
  IECH_OUT = ecr = 0;
  for (ix = 0; ix < db->getNX(0); ix++)
    for (iy = 0; iy < db->getNX(1); iy++)
    {
      flag_col = 1;

      /* Calculate the experimental covariance of total variable */

      st_calculate_covtot(db, ix, iy, flag_sym, cov_ss, cov_nn, num_tot.data(),
                          cov_tot.data());
      if (dbg_ix == ix && dbg_iy == iy)
        st_vario_dump(fildmp, ix, iy, cov_ss, cov_nn, num_tot.data(), cov_tot.data());

      for (iz = 0; iz < db->getNX(2); iz++, ecr++)
      {
        mes_process("3-D Factorial Kriging Analysis", DBOUT->getNSample(),
                    ecr);
        indg[0]  = ix;
        indg[1]  = iy;
        indg[2]  = iz;
        IECH_OUT = db->indiceToRank(indg);
        OptDbg::setCurrentIndex(IECH_OUT + 1);

        /* Initialize the result to TEST */

        DBOUT->setArray(IECH_OUT, IPTR_EST, TEST);

        if (FFFF(db->getZVariable(IECH_OUT, 0)) || !db->isActive(IECH_OUT))
          continue;
        if (OptDbg::query(EDbg::KRIGING) || OptDbg::query(EDbg::NBGH) || OptDbg::query(EDbg::RESULTS))
        {
          mestitle(1, "Target location");
          db_sample_print(db, IECH_OUT, 1, 0, 0, 0);
        }

        /* Look for the neighborhood */

        nbgh_ranks = st_neigh_find(db, ix, iy, iz, nei_ss, nei_nn, nei_cur.data());
        nech       = static_cast<Id>(nbgh_ranks.size());
        if (nech <= 0) continue;
        neq = (nfeq == 0) ? nech : nech + 1;

        /* Check if the neighborhood has changed */

        flag_new = flag_col || st_neigh_diff(nei_ss, nei_nn, nei_ref.data(), nei_cur.data());

        /* If the neighborhood has changed, establish the kriging system */

        flag_col = 0;
        if (flag_new || OptDbg::force())
        {

          /* Establish the L.H.S. of the kriging system */

          st_lhs_exp_3D(nech, nfeq, nei_ss, nei_nn, cov_ss, cov_nn, nei_cur.data(),
                        cov_tot.data(), nugget);
          if (OptDbg::query(EDbg::KRIGING))
            krige_lhs_print(nech, neq, neq, flag_global.data(), lhs_global.data());

          /* Invert the kriging system */

          if (matrix_invert(lhs_global.data(), neq, IECH_OUT))
          {
            status = 1;
            continue;
          }

          /* Establish the R.H.S. of the kriging system */

          st_rhs_exp_3D(nech, nfeq, nei_ss, nei_nn, cov_ss, cov_nn, nei_cur.data(), cov_res.data());
          if (OptDbg::query(EDbg::KRIGING))
            krige_rhs_print(nvarin, nech, neq, neq, flag_global.data(), rhs_global.data());

          /* Derive the kriging weights */

          matrix_product_safe(neq, neq, 1, lhs_global.data(), rhs_global.data(), wgt_global.data());
        }

        /* Calculate the estimation */

        result = st_estim_exp_3D(db, nei_ss, nei_nn, nei_cur.data(), wgt_global.data());
        DBOUT->setArray(IECH_OUT, IPTR_EST, result);
        if (OptDbg::query(EDbg::RESULTS)) st_result_kriging_print(0, nvarin, status);
      }
    }

  /* Set the error return flag */

  error = 0;
  if (fildmp != nullptr) fclose(fildmp);

label_end:
  OptDbg::setCurrentIndex(0);
  (void)krige_koption_manage(-1, 1, EKrigOpt::POINT, 1, VectorInt());
  st_krige_manage_basic(-1, size_nei, size_nei, 1, nfeq);
  return (error);
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
 ** \param[in]  neigh      ANeigh structure
 ** \param[in]  flag_positive  1 for a positive constraints
 ** \param[in]  namconv    Naming convention
 **
 ** \remark  All the variables are estimated using the same model
 ** \remark  In this procedure, we assume that:
 ** \remark  - the problem is multivariate ("z" variables)
 ** \remark  - the constraints is stored in "sum" (only used in dbout)
 **
 *****************************************************************************/
Id krigsum(Db* dbin,
           Db* dbout,
           Model* model,
           ANeigh* neigh,
           bool flag_positive,
           const NamingConvention& namconv)
{
  Id nvar = dbin->getNLoc(ELoc::Z);
  if (model->getNVar() != 1)
  {
    messerr("This procedure requires a monovariate model");
    return 1;
  }
  if (dbout->getNLoc(ELoc::SUM) != 1)
  {
    messerr("This procedure requires one Variable with Locator SUM in the Output Db");
    messerr("The number of such variable is currently equal to %d",
            dbout->getNLoc(ELoc::SUM));
    return 1;
  }

  /* Add the attributes for storing the results */

  Id iptr_est = dbout->addColumnsByConstant(nvar, 0.);
  if (iptr_est < 0) return 1;
  VectorInt active(nvar);
  VectorDouble lterm(nvar);
  VectorInt iuids = dbin->getUIDsByLocator(ELoc::Z);

  /* Setting options */

  // Locally turn the problem to a Monovariate case to have it accepted
  dbin->clearLocators(ELoc::Z);
  dbin->setLocatorByUID(iuids[0], ELoc::Z, 0);
  KrigingSystem ksys(dbin, dbout, model, neigh);
  if (ksys.updKrigOptEstim(iptr_est, -1, -1)) return 1;
  if (ksys.setKrigOptFlagLTerm(true)) return 1;
  if (!ksys.isReady()) return 1;

  /* Loop on the variables */

  for (Id ivar = 0; ivar < nvar; ivar++)
  {
    dbin->clearLocators(ELoc::Z);
    dbin->setLocatorByUID(iuids[ivar], ELoc::Z, 0);
    if (ksys.resetData()) return 1;
    if (ksys.updKrigOptEstim(iptr_est + ivar, -1, -1)) return 1;
    (void)gslSPrintf(string, "Kriging of variable #%d at sample", ivar + 1);

    /* Loop on the targets to be processed */

    for (Id iech_out = 0; iech_out < dbout->getNSample(); iech_out++)
    {
      mes_process(string.data(), dbout->getNSample(), iech_out);
      if (ksys.estimate(iech_out)) return 1;
    }

    // Retrieve Lterm only once per variable (Unique Neighborhood)

    lterm[ivar] = ksys.getLTerm();
  }

  ksys.conclusion();

  // Posterior scaling

  for (Id iech_out = 0; iech_out < dbout->getNSample(); iech_out++)
  {
    bool correct = false;
    for (Id ivar = 0; ivar < nvar; ivar++) active[ivar] = 0;

    /* Implicit loop until the solution is acceptable */

    while (!correct)
    {
      double seistot = 0.;
      double seisloc = dbout->getFromLocator(ELoc::SUM, iech_out);
      for (Id ivar = 0; ivar < nvar; ivar++)
      {
        if (active[ivar]) continue;
        double estim = dbout->getArray(iech_out, iptr_est + ivar);
        seistot += lterm[ivar];
        seisloc -= estim;
      }
      if (isZero(seistot))
      {
        messerr("The sum of scaling terms is zero. No correction is possible");
        return 1;
      }

      for (Id ivar = 0; ivar < nvar; ivar++)
      {
        double estim = 0.;
        if (!active[ivar])
          estim = (dbout->getArray(iech_out, iptr_est + ivar) + lterm[ivar] * seisloc / seistot);
        dbout->setArray(iech_out, iptr_est + ivar, estim);
      }
      correct = true;

      // Correct if negative values are not allowed

      if (flag_positive)
      {
        for (Id ivar = 0; ivar < nvar; ivar++)
        {
          active[ivar] = (dbout->getArray(iech_out, iptr_est + ivar) < 0);
          if (active[ivar]) correct = false;
        }
      }
    }
  }

  /* Reset the locators and rename the output variables */

  dbin->clearLocators(ELoc::Z);
  for (Id ivar = 0; ivar < nvar; ivar++)
    dbin->setLocatorByUID(iuids[ivar], ELoc::Z, ivar);
  namconv.setNamesAndLocators(dbin, VectorString(), ELoc::Z, nvar, dbout, iptr_est, "estim");

  return 0;
}

/****************************************************************************/
/*!
 **  Allocate a vector of sample ranks excluding already selected pivots
 **
 ** \return  Pointer to the newly create integer vector
 **
 ** \param[in]  nech       Number of samples
 ** \param[in]  ranks1     Ranks of exact pivots
 ** \param[in]  ranks2     Ranks of ACP pivots
 **
 ** \remarks The output array must be free by the calling function
 **
 *****************************************************************************/
static VectorInt st_ranks_other(Id nech,
                                const VectorInt& ranks1,
                                const VectorInt& ranks2)
{
  VectorInt rother(nech, 0);
  for (Id i = 0; i < nech; i++)
    rother[i] = i;
  for (Id i = 0, nsize1 = static_cast<Id>(ranks1.size()); i < nsize1; i++)
    rother[ranks1[i]] = -1;
  for (Id i = 0, nsize2 = static_cast<Id>(ranks2.size()); i < nsize2; i++)
    rother[ranks2[i]] = -1;
  return rother;
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
 ** \param[in]  ranks1    Ranks of samples (exact)
 ** \param[in]  ranks2    Ranks of samples (ACP)
 ** \param[in]  rother    Ranks of the idle samples (modified by routine)
 **
 ** \param[out] ntot_arg   Number of pivots
 ** \param[out] nutil_arg  Number of active samples
 ** \param[out] rutil      Rank of the active samples
 ** \param[out] tutil      Returned Matrix for the U array
 ** \param[out] invsig     Returned array for Inverse Sigma
 **
 *****************************************************************************/
static Id st_sampling_krige_data(Db* db,
                                 Model* model,
                                 double beta,
                                 VectorInt& ranks1,
                                 VectorInt& ranks2,
                                 VectorInt& rother,
                                 Id* ntot_arg,
                                 Id* nutil_arg,
                                 VectorInt& rutil,
                                 MatrixDense& tutil,
                                 MatrixSquare& invsig)
{
  Id i, j, ecr, nmax;
  double sumval;
  VectorInt isort;
  VectorDouble vsort;

  /* Initializations */

  Id error  = 1;
  Id ndat   = db->getNSample(true);
  Id nsize1 = static_cast<Id>(ranks1.size());
  Id nsize2 = static_cast<Id>(ranks2.size());
  Id ntot   = nsize1 + nsize2;
  Id npart  = ndat - nsize1;
  Id nutil  = 0;
  MatrixSymmetric mat_s;

  /* Core allocation */

  VectorDouble utab(ndat * ntot, 0.);
  VectorDouble ralls(ndat, 0.);

  /* Defining 'utab' for exact pivots */

  for (i = 0; i < nsize1; i++)
    UTAB(i, i) = 1.;
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
    if (beta > 0.)
    {
      vsort.resize(npart, 0);
      isort.resize(npart, 0);
    }

    mat_s = model->evalCovMatSym(db);

    CholeskyDense mat_s_Chol(mat_s);
    if (!mat_s_Chol.isReady()) goto label_end;
    VectorDouble tl = mat_s_Chol.getLowerTriangle();

    MatrixSymmetric* sq = MatrixSymmetric::createFromTriangle(0, nsize2, tl);

    VectorDouble xl = mat_s_Chol.getUpperTriangleInverse();

    MatrixDense mat_c = model->evalCovMat(db, db, -1, -1, ranks2, rother);

    MatrixDense v;
    mat_s_Chol.matProductInPlace(4, mat_c, v);
    MatrixSymmetric tn1;
    mat_s_Chol.normMatInPlace(1, nsize2, MatrixSymmetric(), tn1);

    auto* tn2 = dynamic_cast<MatrixSymmetric*>(MatrixFactory::prodMatMat(&v, &v, true, false));

    tn1.linearCombination(1, &tn1, 1, tn2);

    if (tn1.computeEigen()) goto label_end;
    VectorDouble eigval  = tn1.getEigenValues();
    MatrixSquare* eigvec = tn1.getEigenVectors()->clone();

    eigvec->prodByDiagInPlace(3, eigval);
    auto* spart = dynamic_cast<MatrixDense*>(MatrixFactory::createGlue(sq, &v, true, false));
    spart->prodMatMatInPlace(spart, eigvec);
    delete eigvec;

    if (beta > 0.)
    {
      for (i = 0; i < npart; i++)
      {
        sumval = 0.;
        for (j = 0; j < nsize2; j++)
          sumval = MAX(sumval, ABS(spart->getValue(i, j)));
        vsort[i] = sumval;
        isort[i] = i;
      }
      ut_sort_double(1, npart, isort.data(), vsort.data());
      nmax = MIN(npart, (Id)(beta * (double)npart));
      for (i = 0; i < nmax; i++)
        for (j = 0; j < nsize2; j++)
          spart->setValue(isort[i], j, 0.);
    }

    for (i = 0; i < npart; i++)
      for (j = 0; j < nsize2; j++)
        UTAB(i + nsize1, j + nsize1) = -spart->getValue(i, j);

    /* Core deallocation */

    delete sq;
    delete tn2;
    delete spart;
  }

  /* Count the number of active samples */

  nutil = 0;
  for (i = 0; i < ndat; i++)
  {
    sumval = 0.;
    for (j = 0; j < ntot; j++)
      sumval += ABS(UTAB(i, j));
    if (sumval > 0.) nutil++;
  }

  /* Create the output arrays */

  rutil.resize(nutil, 0);
  tutil.reset(nutil, ntot);
  invsig.reset(ntot, ntot);

  for (i = ecr = 0; i < ndat; i++)
  {
    sumval = 0.;
    for (j = 0; j < ntot; j++)
      sumval += ABS(UTAB(i, j));
    if (sumval <= 0.) continue;
    rutil[ecr] = ralls[i];
    for (j = 0; j < ntot; j++)
      tutil.setValue(ecr, j, UTAB(i, j));
    ecr++;
  }
  mat_s = model->evalCovMatSym(db, rutil, -1);
  if (matrix_prod_norme(-1, nutil, ntot, tutil.getValues().data(), mat_s.getValues().data(),
                        invsig.getValues().data())) goto label_end;
  if (matrix_invert(invsig.getValues().data(), ntot, 0)) goto label_end;
  utab.clear();

  /* Returning arguments */

  *ntot_arg  = ntot;
  *nutil_arg = nutil;

  /* Error return code */

  error = 0;

label_end:
  return (error);
}

/****************************************************************************/
/*!
 **  Perform the estimation at the data points
 **
 ** \return  Error return code
 **
 ** \param[in]  db         Db structure
 ** \param[in]  model      Model structure
 ** \param[in]  beta       Thresholding value
 ** \param[in]  ranks1     Ranks of exact pivots
 ** \param[in]  ranks2     Ranks of ACP pivots
 ** \param[in]  rother     Ranks of the idle samples
 ** \param[in]  flag_abs   1 Modify 'daata_est' to store the estimation error
 **
 ** \param[out] data_est  Array of estimation at samples
 ** \param[out] data_var  Array of estimation variance at samples
 **
 *****************************************************************************/
Id st_krige_data(Db* db,
                 Model* model,
                 double beta,
                 VectorInt& ranks1,
                 VectorInt& ranks2,
                 VectorInt& rother,
                 Id flag_abs,
                 double* data_est,
                 double* data_var)
{
  VectorInt rutil;
  MatrixDense tutil;
  MatrixSquare invsig;

  /* Core allocation */

  Id nutil = 0;
  Id ntot  = 0;
  Id nech  = db->getNSample();

  /* Perform local sampling */

  if (st_sampling_krige_data(db, model, beta, ranks1, ranks2, rother,
                             &ntot, &nutil, rutil, tutil, invsig))
    return 1;

  /* Second core allocation */

  VectorDouble datm(nutil);
  VectorDouble aux1(ntot);
  VectorDouble aux2(ntot);
  VectorDouble aux3(ntot);
  VectorDouble aux4(ntot);
  VectorDouble s;
  VectorDouble c00;

  /* Get the vector of active data and subtract the mean */

  VectorDouble data = db->getColumnByLocator(ELoc::Z);
  for (Id i = 0; i < nutil; i++)
    datm[i] = data[rutil[i]] - model->getMean(0);
  tutil.prodVecMatInPlace(datm, aux1);
  invsig.prodVecMatInPlace(aux1, aux2);

  /* Perform the estimation at all non pivot samples */

  for (Id iech = 0; iech < nech; iech++)
  {
    data_est[iech] = data_var[iech] = TEST;
    if (!db->isActive(iech)) continue;
    if (rother[iech] < 0) continue;
    VectorInt vech = {iech};
    c00            = model->evalCovMat(db, db, -1, -1, vech, vech).getValues();
    s              = model->evalCovMat(db, db, -1, -1, rutil, vech).getValues();

    tutil.prodVecMatInPlace(s, aux3);
    double estim   = VH::innerProduct(aux2, aux3);
    data_est[iech] = estim + model->getMean(0);

    if (flag_abs)
    {
      double true_value = db->getZVariable(iech, 0);
      if (FFFF(true_value))
        data_est[iech] = TEST;
      else
        data_est[iech] = ABS(data_est[iech] - true_value);
    }

    invsig.prodVecMatInPlace(aux3, aux4);
    double variance = VH::innerProduct(aux3, aux4);
    data_var[iech]  = c00[0] - variance;
  }
  return 0;
}

/****************************************************************************/
/*!
 **  Evaluate the improvement in adding a new pivot on the global score
 **
 ** \return  Error retun code
 **
 ** \param[in]  db         Db structure
 ** \param[in]  model      Model structure
 ** \param[in]  ranks1     Ranks of exact pivots
 ** \param[in]  rother     Ranks of the idle samples
 **
 ** \param[out] crit       Array of criterion
 **
 *****************************************************************************/
Id st_crit_global(Db* db,
                  Model* model,
                  VectorInt& ranks1,
                  VectorInt& rother,
                  double* crit)
{
  Id ecr;
  double estim, sigma, value;

  /* Initializations */

  Id nsize1 = static_cast<Id>(ranks1.size());
  Id ndat   = db->getNSample(true);
  Id nutil  = ndat - nsize1;
  if (nsize1 <= 0) return 1;

  /* Core allocation */

  VectorDouble datm(ndat);
  VectorDouble olderr(nutil);
  VectorDouble olddiv(nutil);
  MatrixDense temp(nsize1, nutil);
  VectorDouble aux1(nutil);
  MatrixSymmetric invc;
  VectorDouble c00;
  VectorDouble cs;
  VectorDouble cs1;
  VectorDouble temp_loc;

  /* Establish the Kriging matrix on the pivot samples */

  invc = model->evalCovMatSym(db, ranks1, -1);
  invc.invert();

  /* Set the data vector (corrected by the mean */

  VectorDouble data = db->getColumnByLocator(ELoc::Z);
  for (Id i = 0; i < nsize1; i++)
    datm[i] = data[ranks1[i]] - model->getMean(0);

  /* Loop on the non-pivots */

  for (Id iech = ecr = 0; iech < ndat; iech++)
  {
    if (!db->isActive(iech)) continue;
    if (rother[iech] < 0) continue;

    VectorInt vech = {iech};

    c00 = model->evalCovMat(db, db, -1, -1, vech, vech).getValues();
    cs  = model->evalCovMat(db, db, -1, -1, ranks1, vech).getValues();
    invc.prodMatVecInPlace(cs, temp_loc);
    temp.setColumn(ecr, temp_loc);

    estim       = VH::innerProduct(datm, temp_loc);
    olderr[ecr] = estim + model->getMean(0) - db->getZVariable(iech, 0);

    sigma       = VH::innerProduct(cs, temp_loc);
    olddiv[ecr] = olderr[ecr] / (c00[0] - sigma);
    ecr++;
  }

  /* Loop on the candidates */

  for (Id iech = ecr = 0; iech < ndat; iech++)
  {
    crit[iech] = TEST;
    if (!db->isActive(iech)) continue;
    if (rother[iech] < 0) continue;

    VectorInt vech = {iech};

    cs  = model->evalCovMat(db, db, -1, -1, vech, ranks1).getValues();
    cs1 = model->evalCovMat(db, db, -1, -1, vech, rother).getValues();

    temp.prodVecMatInPlace(cs, aux1);
    VH::linearCombinationInPlace(1., cs1, -1, aux1, cs1);
    VH::linearCombinationInPlace(1., olderr, -olddiv[ecr], cs1, cs1);

    value = 0.;
    for (Id i = 0; i < nutil; i++)
      value += cs1[i] * cs1[i];
    crit[iech] = value / nutil;

    ecr++;
  }
  return 0;
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
 ** \param[in]  ranks1     Ranks of exact pivots
 ** \param[in]  method2    Criterion for choosing ACP pivots
 **                        1 : Local evaluation
 **                        2 : Global evaluation
 ** \param[in]  nsize2_max Maximum number of ACP pivots
 ** \param[in]  ranks2     Ranks of ACP pivots
 ** \param[in]  verbose    1 for a verbose output
 **
 *****************************************************************************/
Id sampling_f(Db* db,
              Model* model,
              double beta,
              Id method1,
              Id nsize1_max,
              VectorInt& ranks1,
              Id method2,
              Id nsize2_max,
              VectorInt& ranks2,
              Id verbose)
{
  Id best_rank;
  double best_ecart;

  /* Initializations */

  Id nsize1 = static_cast<Id>(ranks1.size());
  Id nsize2 = static_cast<Id>(ranks2.size());
  Id nech   = db->getNSample();

  /* Preliminary checks */

  if (method2 != 1)
  {
    messerr("The Global Evaluation method for choosing ACP pivots");
    messerr("has not been programmed yet");
    return 1;
  }
  if (nsize1_max > 0 && nsize1 == 0)
  {
    messerr("The sampling requires a first sample to be defined 'ranks1'");
    return 1;
  }

  /* Core allocation */

  VectorDouble data_est(nech);
  VectorDouble data_var(nech);
  VectorInt rother = st_ranks_other(nech, ranks1, ranks2);

  /* Sample the exact pivots */

  while (nsize1 < nsize1_max)
  {
    if (method1 == 1)
    {
      if (st_krige_data(db, model, beta, ranks1, ranks2, rother,
                        1, data_est.data(), data_var.data())) return 1;
      best_rank  = VectorHelper::whereMaximum(VectorHelper::initVDouble(data_est.data(), nech));
      best_ecart = data_est[best_rank];
    }
    else
    {
      if (st_crit_global(db, model, ranks1, rother, data_est.data())) return 1;
      best_rank  = VectorHelper::whereMinimum(VectorHelper::initVDouble(data_est.data(), nech));
      best_ecart = data_est[best_rank];
    }
    if (verbose)
      message("Exact Pivots (%3d/%3d): Rank = %3d - value = %lf\n", nsize1 + 1,
              nsize1_max, best_rank + 1, best_ecart);
    ranks1[nsize1]    = best_rank;
    rother[best_rank] = -1;
    nsize1++;
  }

  /* Sample the ACP pivots */

  while (nsize2 < nsize2_max)
  {
    if (st_krige_data(db, model, beta, ranks1, ranks2, rother,
                      1, data_est.data(), data_var.data())) return 1;
    best_rank  = VectorHelper::whereMaximum(VectorHelper::initVDouble(data_est.data(), nech));
    best_ecart = data_est[best_rank];
    if (verbose)
      message("ACP   Pivots (%3d/%3d): Rank = %3d - value = %lf\n", nsize2 + 1,
              nsize2_max, best_rank + 1, best_ecart);
    ranks2[nsize2]    = best_rank;
    rother[best_rank] = -1;
    nsize2++;
  }

  /* Calculation of statistics on reproduction errors */

  if (verbose)
  {
    if (st_krige_data(db, model, beta, ranks1, ranks2, rother,
                      1, data_est.data(), data_var.data())) return 1;
    StatResults stats = ut_statistics(nech, data_est.data());
    mestitle(1, "Statistics on estimation errors");
    message("Count   = %d \n", stats.nvalid);
    message("Minimum = %lf\n", stats.mini);
    message("Mean    = %lf\n", stats.mean);
    message("St. Dev.= %lf\n", stats.stdv);
    message("Maximum = %lf\n", stats.maxi);
  }

  return 0;
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
 ** \param[in]  ranks1     Ranks of exact pivots
 ** \param[in]  ranks2     Ranks of ACP pivots
 ** \param[in]  flag_std   Option for storing the standard deviation
 ** \param[in]  verbose    Verbose flag
 **
 *****************************************************************************/
Id krigsampling_f(Db* dbin,
                  Db* dbout,
                  Model* model,
                  double beta,
                  VectorInt& ranks1,
                  VectorInt& ranks2,
                  bool flag_std,
                  Id verbose)
{
  Id ntot, nutil, i;
  double estim;
  VectorInt rutil;
  MatrixDense tutil;
  MatrixSquare invsig;

  /* Preliminary checks */

  Id nsize1    = static_cast<Id>(ranks1.size());
  Id nsize2    = static_cast<Id>(ranks2.size());
  double sigma = 0.;
  st_global_init(dbin, dbout);
  FLAG_EST = true;
  FLAG_STD = flag_std;
  if (st_check_environment(1, 1, model)) return 1;
  Id nvar = model->getNVar();
  Id nech = dbin->getNSample();

  /* Preliminary checks */

  if (nvar != 1)
  {
    messerr("This method is only programmed for monovariate case");
    return 1;
  }
  if (nsize1 + nsize2 <= 0)
  {
    messerr("You must specify some pivots in 'ranks1' or 'ranks2'");
    return 1;
  }

  /* Add the attributes for storing the results */

  if (FLAG_EST)
  {
    IPTR_EST = dbout->addColumnsByConstant(nvar, 0.);
    if (IPTR_EST < 0) return 1;
  }
  if (FLAG_STD)
  {
    IPTR_STD = dbout->addColumnsByConstant(nvar, 0.);
    if (IPTR_STD < 0) return 1;
  }

  /* Core allocation */

  VectorInt rother = st_ranks_other(nech, ranks1, ranks2);

  /* Perform local sampling */

  if (st_sampling_krige_data(dbin, model, beta, ranks1, ranks2,
                             rother, &ntot, &nutil, rutil, tutil, invsig))
    return 1;

  /* Optional printout */

  if (verbose)
  {
    message("Printout of intermediate arrays\n");
    print_imatrix("Pivot ranks", 0, 1, 1, ntot, NULL, rutil.data());
    print_matrix("Inv-Sigma", 0, 1, ntot, ntot, NULL, invsig.getValues().data());
    print_matrix("U", 0, 1, ntot, nutil, NULL, tutil.getValues().data());
  }

  /* Second core allocation */

  VectorDouble datm(nutil, 0);
  VectorDouble aux1(ntot, 0);
  VectorDouble aux2(ntot, 0);
  VectorDouble aux3(ntot, 0);
  VectorDouble aux4;
  VectorDouble s;
  VectorDouble c00;
  if (FLAG_STD)
    aux4.resize(ntot, 0);

  /* Get the vector of active data and subtract the mean */

  VectorDouble data = dbin->getColumnByLocator(ELoc::Z);
  for (i = 0; i < nutil; i++)
    datm[i] = data[rutil[i]] - model->getMean(0);
  tutil.prodVecMatInPlace(datm, aux1);
  invsig.prodVecMatInPlace(aux1, aux2);

  /* Loop on the target samples */

  for (IECH_OUT = 0; IECH_OUT < DBOUT->getNSample(); IECH_OUT++)
  {
    mes_process("Kriging sample", DBOUT->getNSample(), IECH_OUT);
    OptDbg::setCurrentIndex(IECH_OUT + 1);
    if (!dbout->isActive(IECH_OUT)) continue;
    if (OptDbg::query(EDbg::KRIGING) || OptDbg::query(EDbg::NBGH) || OptDbg::query(EDbg::RESULTS))
    {
      mestitle(1, "Target location");
      db_sample_print(dbout, IECH_OUT, 1, 0, 0, 0);
    }

    VectorInt vech = {IECH_OUT};
    s              = model->evalCovMat(dbin, dbout, -1, -1, rutil, vech).getValues();
    if (FLAG_STD)
      c00 = model->evalCovMat(dbout, dbout, -1, -1, vech, vech).getValues();

    tutil.prodVecMatInPlace(s, aux3);
    estim = VH::innerProduct(aux2, aux3) + model->getMean(0);
    DBOUT->setArray(IECH_OUT, IPTR_EST, estim);

    if (FLAG_STD)
    {
      invsig.prodVecMatInPlace(aux3, aux4);
      sigma = VH::innerProduct(aux3, aux4);
      sigma = c00[0] - sigma;
      sigma = (sigma > 0) ? sqrt(sigma) : 0.;
      DBOUT->setArray(IECH_OUT, IPTR_STD, sigma);
    }

    /* Optional printout */

    if (OptDbg::query(EDbg::RESULTS))
    {
      tab_printg(" - Estimate  = ", estim);
      message("\n");
      if (FLAG_STD)
      {
        tab_printg(" - Std. Dev. = ", sigma);
        message("\n");
      }
    }
  }
  return 0;
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
static void st_declustering_stats(Id mode, Id method, Db* db, Id iptr)
{
  double mean, var, zval, coeff, mini, maxi, sumwgt;

  mean = var = sumwgt = 0.;
  mini                = MAXIMUM_BIG;
  maxi                = MINIMUM_BIG;
  for (Id iech = 0; iech < db->getNSample(); iech++)
  {
    if (!db->isActive(iech)) continue;
    zval = db->getZVariable(iech, 0);
    if (FFFF(zval)) continue;
    coeff = (mode == 0) ? 1. : db->getArray(iech, iptr);
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
 **
 *****************************************************************************/
static void st_declustering_truncate_and_rescale(Db* db, Id iptr)
{
  double total, coeff;

  /* Truncate the negative weights */

  total = 0;
  for (Id iech = 0; iech < db->getNSample(); iech++)
  {
    if (!db->isActive(iech)) continue;
    if (FFFF(db->getZVariable(iech, 0))) continue;
    coeff = db->getArray(iech, iptr);
    if (coeff < 0)
      db->setArray(iech, iptr, 0.);
    else
      total += coeff;
  }

  /* Rescale */

  for (Id iech = 0; iech < db->getNSample(); iech++)
  {
    if (!db->isActive(iech)) continue;
    if (FFFF(db->getZVariable(iech, 0))) continue;
    db->updArray(iech, iptr, EOperator::DIVIDE, total);
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
static Id st_declustering_1(Db* db, Id iptr, const VectorDouble& radius)
{
  Id ndim = db->getNDim();
  VectorDouble vect(ndim);

  if (radius.empty())
  {
    messerr("This method requires the definition of 'radius'");
    return 1;
  }

  /* Loop on the target sample */

  for (Id iech = 0; iech < db->getNSample(); iech++)
  {
    if (!db->isActive(iech)) continue;
    if (FFFF(db->getZVariable(iech, 0))) continue;

    /* Loop on the second sample */

    for (Id jech = 0; jech < db->getNSample(); jech++)
    {
      if (!db->isActive(jech)) continue;
      double value = db->getZVariable(iech, 0);
      if (FFFF(value)) continue;
      (void)distance_intra(db, iech, jech, vect.data());

      /* Normalize the distance */

      double dist = 0.;
      for (Id idim = 0; idim < db->getNDim(); idim++)
      {
        vect[idim] /= radius[idim];
        dist += vect[idim] * vect[idim];
      }
      if (dist > 1) continue;
      db->updArray(iech, iptr, EOperator::ADD, 1);
    }
  }

  /* Normalization step */

  double total = 0.;
  for (Id iech = 0; iech < db->getNSample(); iech++)
  {
    if (!db->isActive(iech)) continue;
    if (FFFF(db->getZVariable(iech, 0))) continue;
    total += 1. / db->getArray(iech, iptr);
  }
  for (Id iech = 0; iech < db->getNSample(); iech++)
  {
    if (!db->isActive(iech)) continue;
    if (FFFF(db->getZVariable(iech, 0))) continue;
    db->setArray(iech, iptr, 1. / db->getArray(iech, iptr) / total);
  }
  return 0;
}

/****************************************************************************/
/*!
 **  Perform the Declustering task as weight of the Mean Kriging (Unique Neigh)
 **
 ** \return  Error return code
 **
 ** \param[in]  db         input Db structure
 ** \param[in]  model      Model structure
 ** \param[in]  neigh      ANeigh structure (should be Unique)
 ** \param[in]  iptr       Rank of the declustering weight
 **
 *****************************************************************************/
static Id st_declustering_2(Db* db,
                            Model* model,
                            ANeigh* neigh,
                            Id iptr)
{
  KrigOpt krigopt;
  krigopt.setOptionCalcul(EKrigOpt::DRIFT);
  KrigingSystem ksys(db, db, model, neigh, krigopt);
  if (ksys.setKrigOptDataWeights(iptr, true)) return 1;
  if (!ksys.isReady()) return 1;
  if (ksys.estimate(0)) return 1;
  ksys.conclusion();

  /* Truncate the negative weights */

  st_declustering_truncate_and_rescale(db, iptr);

  return 0;
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
 ** \param[in]  model      Model structure
 ** \param[in]  neigh      ANeigh structure
 ** \param[in]  ndiscs     Array of discretization counts
 ** \param[in]  iptr       Rank of the declustering weight
 **
 *****************************************************************************/
static Id st_declustering_3(Db* db,
                            Db* dbgrid,
                            Model* model,
                            ANeigh* neigh,
                            const VectorInt& ndiscs,
                            Id iptr)
{
  // Preliminary checks

  if (neigh == nullptr)
  {
    messerr("This function requires a Neighborhood");
    return 1;
  }
  if (neigh->getType() == ENeigh::IMAGE)
  {
    messerr("This tool cannot function with an IMAGE neighborhood");
    return 1;
  }
  if (ndiscs.empty())
  {
    messerr("The Cell discretization must be provided");
    return 1;
  }

  /* Setting options */

  KrigOpt krigopt;
  krigopt.setOptionCalcul(EKrigOpt::BLOCK, ndiscs);
  KrigingSystem ksys(db, dbgrid, model, neigh, krigopt);
  if (ksys.setKrigOptDataWeights(iptr, false)) return 1;
  if (!ksys.isReady()) return 1;

  /* Loop on the targets to be processed */

  for (Id iech_out = 0; iech_out < dbgrid->getNSample(); iech_out++)
  {
    mes_process("Kriging sample", dbgrid->getNSample(), iech_out);
    if (ksys.estimate(iech_out)) return 1;
  }

  ksys.conclusion();

  /* Truncate the negative weights */

  st_declustering_truncate_and_rescale(db, iptr);

  return 0;
}

/****************************************************************************/
/*!
 **  Perform the Declustering task
 **
 ** \return  Error return code
 **
 ** \param[in]  dbin       input Db structure
 ** \param[in]  model      Model structure
 ** \param[in]  method     Method for declustering
 ** \param[in]  neigh      ANeigh structure
 ** \param[in]  dbgrid     Grid auxiliary Db structure
 ** \param[in]  radius     Array of neighborhood radius
 ** \param[in]  ndiscs     Array of discretization
 ** \param[in]  flag_sel   1 to mask off samples with zero weight
 ** \param[in]  verbose    Verbose option
 **
 *****************************************************************************/
Id declustering(Db* dbin,
                Model* model,
                Id method,
                ANeigh* neigh,
                DbGrid* dbgrid,
                const VectorDouble& radius,
                const VectorInt& ndiscs,
                Id flag_sel,
                bool verbose)
{
  if (!dbin->isNVarComparedTo(1)) return 1;

  /* Add the kriging weight as a new variable */

  Id iptr = dbin->addColumnsByConstant(1, 0.);
  if (iptr < 0) return 1;

  /* Produce statistics on the target variable before declustering */

  if (verbose) st_declustering_stats(0, method, dbin, iptr);

  /* Dispatch */

  switch (method)
  {
    case 1: /* Weight proportional to number of samples */
    {
      if (st_declustering_1(dbin, iptr, radius)) return 1;
      break;
    }

    case 2: /* Weight of the Mean */
    {
      if (model == nullptr)
      {
        messerr("A Model is needed for this declustering method");
        return 1;
      }
      if (st_declustering_2(dbin, model, neigh, iptr)) return 1;
      break;
    }

    case 3: /* Average weight of the Block Kriging */
    {
      if (model == nullptr)
      {
        messerr("A Model is needed for this declustering method");
        return 1;
      }
      if (st_declustering_3(dbin, dbgrid, model, neigh, ndiscs, iptr))
        return 1;
      break;
    }

    default:
      messerr("Not yet implemented");
      return 1;
  }

  /* Store the selection (optional) */

  if (flag_sel)
  {
    Id iptr_sel = dbin->addColumnsByConstant(1, 0.);
    if (iptr_sel < 0) return 1;
    for (Id iech = 0; iech < dbin->getNSample(); iech++)
    {
      dbin->setArray(iech, iptr_sel, 0.);
      if (!dbin->isActive(iech)) continue;
      double indic = (dbin->getArray(iech, iptr) > 0.);
      dbin->setArray(iech, iptr_sel, indic);
    }
  }

  /* Produce statistics on the target variable after declustering */

  if (verbose) st_declustering_stats(1, method, dbin, iptr);

  return 0;
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
static VectorDouble st_calcul_covmat(const char* title,
                                     Db* db1,
                                     Id test_def1,
                                     Db* db2,
                                     Id test_def2,
                                     Model* model)
{
  Id n1, n2, i1, i2;
  VectorDouble covgen;

  /* Initializations */

  n1 = (test_def1) ? db1->getNSampleActiveAndDefined(0) : db1->getNSample(true);
  n2 = (test_def2) ? db2->getNSampleActiveAndDefined(0) : db2->getNSample(true);

  /* Core allocation */

  covgen.resize(n1 * n2);

  for (Id ii1 = i1 = 0; ii1 < db1->getNSample(); ii1++)
  {
    if (test_def1)
    {
      if (!db1->isActiveAndDefined(ii1, 0)) continue;
    }
    else
    {
      if (!db1->isActive(ii1)) continue;
    }

    for (Id ii2 = i2 = 0; ii2 < db2->getNSample(); ii2++)
    {
      if (test_def2)
      {
        if (!db2->isActiveAndDefined(ii2, 0)) continue;
      }
      else
      {
        if (!db2->isActive(ii2)) continue;
      }

      for (Id idim = 0; idim < db1->getNDim(); idim++)
        d1_global[idim] = db1->getDistance1D(ii1, ii2, idim);

      COVGEN(i1, i2) = model->evaluateOneGeneric(nullptr, d1_global);
      i2++;
    }
    i1++;
  }

  /* Optional printout */

  if (INH_FLAG_VERBOSE)
    print_matrix(title, INH_FLAG_LIMIT, 1, n2, n1, NULL, covgen.data());

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
static VectorDouble st_calcul_drfmat(const char* title,
                                     Db* db1,
                                     Id test_def1,
                                     Model* model)
{
  Id i1, n1, nbfl;
  VectorDouble drftab;

  /* Initializations */

  n1   = (test_def1) ? db1->getNSampleActiveAndDefined(0) : db1->getNSample(true);
  nbfl = model->getNDrift();

  /* Core allocation */

  drftab.resize(n1 * nbfl);

  /* Loop on the samples */

  i1 = 0;
  for (Id ii1 = 0; ii1 < db1->getNSample(); ii1++)
  {
    if (test_def1)
    {
      if (!db1->isActiveAndDefined(ii1, 0)) continue;
    }
    else
    {
      if (!db1->isActive(ii1)) continue;
    }

    VectorDouble drfloc = model->evalDriftBySample(db1, ii1, ECalcMember::LHS);
    (void)memcpy(&drftab[i1 * nbfl], drfloc.data(), nbfl * sizeof(double));
    i1++;
  }

  /* Optional printout */

  if (INH_FLAG_VERBOSE)
    print_matrix(title, INH_FLAG_LIMIT, 1, nbfl, n1, NULL, drftab.data());

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
static VectorDouble st_calcul_distmat(const char* title,
                                      Db* db1,
                                      Id test_def1,
                                      Db* db2,
                                      Id test_def2,
                                      double power)
{
  Id n1, ns, i1, is, ndim;
  VectorDouble distgen;
  double dist;

  /* Initializations */

  n1   = (test_def1) ? db1->getNSampleActiveAndDefined(0) : db1->getNSample(true);
  ns   = (test_def2) ? db2->getNSampleActiveAndDefined(0) : db2->getNSample(true);
  ndim = db1->getNDim();

  /* Core allocation */

  distgen.resize(n1 * ns);

  for (Id ii1 = i1 = 0; ii1 < db1->getNSample(); ii1++)
  {
    if (test_def1)
    {
      if (!db1->isActiveAndDefined(ii1, 0)) continue;
    }
    else
    {
      if (!db1->isActive(ii1)) continue;
    }

    for (Id iis = is = 0; iis < db2->getNSample(); iis++)
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
      for (Id idim = 0; idim < ndim; idim++)
      {
        d1_global[idim] = db1->getDistance1D(ii1, iis, idim);
        dist += d1_global[idim] * d1_global[idim];
      }

      DISTGEN(i1, is) = 1. / pow(dist, power / 2.);
      is++;
    }
    i1++;
  }

  /* Optional printout */

  if (INH_FLAG_VERBOSE)
    print_matrix(title, INH_FLAG_LIMIT, 1, ns, n1, NULL, distgen.data());

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
static VectorDouble st_calcul_product(const char* title,
                                      Id n1,
                                      Id ns,
                                      const VectorDouble& covss,
                                      const VectorDouble& distgen)
{
  VectorDouble prodgen;

  prodgen.resize(n1 * ns);

  for (Id i1 = 0; i1 < n1; i1++)
    for (Id is = 0; is < ns; is++)
    {
      PRODGEN(i1, is) = 0.;
      for (Id js = 0; js < ns; js++)
        PRODGEN(i1, is) += COVSS(is, js) * DISTGEN(i1, js);
    }

  /* Optional printout */

  if (INH_FLAG_VERBOSE)
    print_matrix(title, INH_FLAG_LIMIT, 1, ns, n1, NULL, prodgen.data());

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
static VectorDouble st_inhomogeneous_covpp(Db* dbdat,
                                           Db* dbsrc,
                                           Model* model_dat,
                                           const VectorDouble& distps,
                                           const VectorDouble& prodps)
{
  Id np, ns;

  /* Initializations */

  np = dbdat->getNSampleActiveAndDefined(0);
  ns = dbsrc->getNSample(true);

  /* Covariance matrix between Mesures */

  auto covpp = st_calcul_covmat("Covariance P-P", dbdat, 1, dbdat, 1, model_dat);

  /* Calculate the LHS matrix */

  for (Id ip = 0; ip < np; ip++)
    for (Id jp = ip; jp < np; jp++)
      for (Id is = 0; is < ns; is++)
      {
        COVPP(ip, jp) += DISTPS(ip, is) * PRODPS(jp, is);
        if (jp > ip) COVPP(jp, ip) = COVPP(ip, jp);
      }

  /* Set the error return code */

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
static VectorDouble st_inhomogeneous_covgp(Db* dbdat,
                                           Db* dbsrc,
                                           Db* dbout,
                                           Id flag_source,
                                           Model* model_dat,
                                           const double* distps,
                                           const double* prodps,
                                           const double* prodgs)
{
  Id np, ns, ng;

  /* Initializations */

  np = dbdat->getNSampleActiveAndDefined(0);
  ns = dbsrc->getNSample(true);
  ng = dbout->getNSample(true);

  /* Covariance matrix between Mesures and Target */

  auto covgp = st_calcul_covmat("Covariance G-P", dbout, 0, dbdat, 1, model_dat);

  /* Add the contribution of the source */

  if (!flag_source)
  {
    for (Id ig = 0; ig < ng; ig++)
      for (Id ip = 0; ip < np; ip++)
        for (Id is = 0; is < ns; is++)
          COVGP(ig, ip) += DISTPS(ip, is) * PRODGS(ig, is);
  }
  else
  {
    for (Id ig = 0; ig < ng; ig++)
      for (Id ip = 0; ip < np; ip++)
        COVGP(ig, ip) = PRODPS(ip, ig);
  }

  /* Set the error return code */

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
static VectorDouble st_inhomogeneous_covgg(Db* dbsrc,
                                           Db* dbout,
                                           Id flag_source,
                                           Model* model_dat,
                                           const double* distgs,
                                           const VectorDouble& prodgs)
{
  Id ns = dbsrc->getNSample(true);
  Id ng = dbout->getNSample(true);
  VectorDouble covgg(ng, 0);

  /* Calculate the variance term (for a zero-distance) */

  double c00 = model_dat->evaluateOneGeneric(nullptr, VectorDouble());

  /* Calculate the variance vector */

  if (!flag_source)
  {
    for (Id ig = 0; ig < ng; ig++)
    {
      covgg[ig] = c00;
      for (Id is = 0; is < ns; is++)
        covgg[ig] += DISTGS(ig, is) * PRODGS(ig, is);
    }
  }
  else
  {
    for (Id ig = 0; ig < ng; ig++)
      covgg[ig] = c00;
  }
  return covgg;
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
 ** \param[out] ymat        Array: t(F) %*% C^-1
 ** \param[out] zmat        Array: (t(F) %*% C^-1 %*% F)^-1
 **
 ** \remarks The returned arrays 'yloc' and 'zloc' must be freed by the
 ** \remarks calling function
 **
 *****************************************************************************/
static Id st_drift_prepar(Id np,
                          Id nbfl,
                          const double* covpp,
                          const double* drftab,
                          VectorDouble& ymat,
                          VectorDouble& zmat)
{
  double value;
  Id error, ecr;

  /* Initialization */

  error = 1;

  /* First returned array */

  ymat.resize(nbfl * np);

  ecr = 0;
  for (Id il = 0; il < nbfl; il++)
    for (Id ip = 0; ip < np; ip++)
    {
      value = 0.;
      for (Id jp = 0; jp < np; jp++)
        value += COVPP(ip, jp) * DRFTAB(jp, il);
      ymat[ecr++] = value;
    }

  /* Second retrned array */

  zmat.resize(nbfl * nbfl);

  ecr = 0;
  for (Id il = 0; il < nbfl; il++)
    for (Id jl = 0; jl < nbfl; jl++)
    {
      value = 0.;
      for (Id ip = 0; ip < np; ip++)
        value += YMAT(ip, il) * DRFTAB(ip, jl);
      zmat[ecr++] = value;
    }

  /* Invert 'zmat' */

  if (matrix_invert(zmat.data(), nbfl, -1)) goto label_end;

  /* Set the error return code */

  error = 0;

label_end:
  if (error)
  {
    ymat.clear();
    zmat.clear();
  }
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
static void st_drift_update(Id np,
                            Id nbfl,
                            const double* covgp,
                            const double* driftg,
                            const double* ymat,
                            double* zmat,
                            double* maux,
                            double* lambda,
                            double* mu)
{
  double value;

  /* Calculate the Lagrange vector */

  for (Id il = 0; il < nbfl; il++)
  {
    value = 0.;
    for (Id ip = 0; ip < np; ip++)
      value = YMAT(ip, il) * covgp[ip] - driftg[il];
    maux[il] = value;
  }
  matrix_product_safe(nbfl, nbfl, 1, zmat, maux, mu);

  /* Update the vector of kriging weights */

  for (Id ip = 0; ip < np; ip++)
    for (Id il = 0; il < nbfl; il++)
      lambda[ip] -= YMAT(ip, il) * mu[il];
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
Id inhomogeneous_kriging(Db* dbdat,
                         Db* dbsrc,
                         Db* dbout,
                         double power,
                         Id flag_source,
                         Model* model_dat,
                         Model* model_src)
{
  Id error, np, ip, ns, ng, nvar, neq, nred, nfeq, nbfl;
  VectorDouble covss, distps, distgs, covpp, covgp, prodps, prodgs;
  VectorDouble driftp, ymat, zmat;
  double* rhs;
  double estim, stdev, auxval;
  VectorInt nbgh_ranks;
  VectorDouble driftg;
  VectorDouble covgg;
  VectorDouble lambda;
  VectorDouble data;
  VectorDouble mu;
  VectorDouble maux;

  /* Preliminary checks */

  error = nvar        = 1;
  NeighUnique* neighU = NeighUnique::create(false);
  st_global_init(dbdat, dbout);
  FLAG_EST = true;
  FLAG_STD = true;
  if (st_check_environment(1, 1, model_dat)) goto label_end;

  /* Preliminary checks */

  if (model_dat->getNVar() != nvar)
  {
    messerr("The Model for the Data must be Monovariate");
    goto label_end;
  }
  if (model_src->getNVar() != nvar)
  {
    messerr("The Model for the Sources must be Monovariate");
    goto label_end;
  }

  /* Add the attributes for storing the results */

  if (FLAG_EST)
  {
    IPTR_EST = dbout->addColumnsByConstant(nvar, 0.);
    if (IPTR_EST < 0) goto label_end;
  }
  if (FLAG_STD)
  {
    IPTR_STD = dbout->addColumnsByConstant(nvar, 0.);
    if (IPTR_STD < 0) goto label_end;
  }
  nred = neq = np = dbdat->getNSampleActiveAndDefined(0);
  nfeq            = 0;
  ns              = dbsrc->getNSample(true);
  ng              = dbout->getNSample(true);
  nbfl            = model_dat->getNDrift();

  /* Core allocation */

  lambda.resize(np);
  data.resize(np);

  /* Pre-calculations */

  if (st_model_manage(1, model_dat)) goto label_end;
  if (st_krige_manage(1, nvar, model_dat, neighU)) goto label_end;
  if (krige_koption_manage(1, 1, EKrigOpt::POINT, 1, VectorInt())) goto label_end;

  /* Constitute the Data vector */

  for (Id iip = ip = 0; iip < dbdat->getNSample(); iip++)
  {
    if (!dbdat->isActiveAndDefined(iip, 0)) continue;
    data[ip] = dbdat->getZVariable(iip, 0);
    ip++;
  }

  /* Establish the covariance matrix between Sources */

  covss = st_calcul_covmat("Covarance S_S", dbsrc, 0, dbsrc, 0, model_src);

  /* Establish the distance matrix between Data and Sources */

  distps = st_calcul_distmat("Distance P-S", dbdat, 1, dbsrc, 0, power);

  /* Establish the distance matrix between Target and Sources */

  if (!flag_source)
  {
    distgs = st_calcul_distmat("Distance G-S", dbout, 0, dbsrc, 0, power);
  }

  /* Establish the Data-Source Product matrix */

  prodps = st_calcul_product("Convolve P-S", np, ns, covss, distps);

  /* Establish the complete kriging matrix */

  covpp = st_inhomogeneous_covpp(dbdat, dbsrc, model_dat, distps, prodps);
  if (OptDbg::query(EDbg::KRIGING) || OptDbg::isReferenceDefined())
    krige_lhs_print(np, neq, nred, NULL, covpp.data());

  /* Invert the Kriging Matrix */

  if (matrix_invert(covpp.data(), np, -1)) goto label_end;

  /* Establish the drift at Data */

  if (nbfl > 0)
  {
    mu.resize(nbfl);
    maux.resize(nbfl);

    driftp = st_calcul_drfmat("Drift P", dbdat, 1, model_dat);

    /* Prepare auxiliary arrays */

    if (st_drift_prepar(np, nbfl, covpp.data(), driftp.data(), ymat, zmat)) goto label_end;
  }

  /* Establish the Target-Source Product matrix */

  if (!flag_source)
  {
    prodgs = st_calcul_product("Convolve G-S", ng, ns, covss, distgs);
  }

  /* Establish the COVGP */

  covgp = st_inhomogeneous_covgp(dbdat, dbsrc, dbout, flag_source, model_dat,
                                 distps.data(), prodps.data(), prodgs.data());

  /* Establish the drift at Target */

  if (nbfl > 0) driftg.resize(nbfl);

  /* Establish the variance at targets */

  covgg = st_inhomogeneous_covgg(dbsrc, dbout, flag_source, model_dat, distgs.data(),
                                 prodgs);

  /* Loop on the targets to be processed */

  for (IECH_OUT = 0; IECH_OUT < DBOUT->getNSample(); IECH_OUT++)
  {
    mes_process("Kriging sample", DBOUT->getNSample(), IECH_OUT);
    OptDbg::setCurrentIndex(IECH_OUT + 1);
    if (!dbout->isActive(IECH_OUT)) continue;
    if (OptDbg::query(EDbg::KRIGING) || OptDbg::query(EDbg::NBGH) || OptDbg::query(EDbg::RESULTS))
    {
      mestitle(1, "Target location");
      db_sample_print(dbout, IECH_OUT, 1, 0, 0, 0);
    }

    // Neighborhood search

    neighU->select(IECH_OUT, nbgh_ranks);
    rhs = &COVGP(IECH_OUT, 0);

    /* Optional printout of the R.H.S */

    if (OptDbg::force()) krige_rhs_print(nvar, np, neq, nred, NULL, rhs);

    /* Fill the drift at Target point (optional) */

    model_dat->evalDriftBySampleInPlace(dbout, IECH_OUT, ECalcMember::LHS, driftg);

    /* Calculate the Kriging weights */

    matrix_product_safe(np, np, 1, covpp.data(), rhs, lambda.data());
    if (OptDbg::force())
      st_krige_wgt_print(0, nvar, nvar, nfeq, nbgh_ranks, nred, -1, NULL,
                         lambda.data());

    /* Update vector of weights in presence of drift */

    if (nbfl > 0)
    {

      /* Evaluate the drift at Target */

      model_dat->evalDriftBySampleInPlace(dbout, IECH_OUT, ECalcMember::LHS, driftg);

      /* Update the kriging weights */

      st_drift_update(np, nbfl, rhs, driftg.data(), ymat.data(), zmat.data(),
                      maux.data(), lambda.data(), mu.data());
    }

    /* Perform the estimation */

    matrix_product_safe(1, np, 1, data.data(), lambda.data(), &estim);
    matrix_product_safe(1, np, 1, rhs, lambda.data(), &stdev);

    /* Update the variance in presence of drift */

    if (nbfl > 0)
    {
      matrix_product_safe(1, nbfl, 1, mu.data(), maux.data(), &auxval);
      stdev += auxval;
    }

    /* Update the variance calculation */

    VAR0(0, 0) = covgg[IECH_OUT];
    stdev      = covgg[IECH_OUT] - stdev;
    stdev      = (stdev > 0) ? sqrt(stdev) : 0.;

    /* Store the result */

    dbout->setArray(IECH_OUT, IPTR_EST, estim);
    dbout->setArray(IECH_OUT, IPTR_STD, stdev);

    /* Optional printout */

    if (OptDbg::query(EDbg::KRIGING) || OptDbg::force())
      st_result_kriging_print(0, nvar, 0);
  }

  /* Set the error return flag */

  error = 0;

label_end:
  OptDbg::setCurrentIndex(0);
  (void)st_model_manage(-1, model_dat);
  (void)st_krige_manage(-1, 1, model_dat, neighU);
  (void)krige_koption_manage(-1, 1, EKrigOpt::POINT, 1, VectorInt());
  delete neighU;
  return (error);
}

} // namespace gstlrn
