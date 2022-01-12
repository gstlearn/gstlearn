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
#include "geoslib_f_private.h"
#include "Variogram/Vario.hpp"
#include "Basic/Utilities.hpp"
#include "Basic/Law.hpp"
#include "Basic/MathFunc.hpp"
#include "Stats/Classical.hpp"
#include "Basic/AException.hpp"
#include "LithoRule/Rule.hpp"
#include "LithoRule/RuleShift.hpp"
#include "LithoRule/RuleProp.hpp"
#include "LithoRule/PropDef.hpp"
#include "Db/Db.hpp"
#include "Model/Model.hpp"

#include <math.h>
#include <string.h>

/*! \cond */
typedef struct
{
  int opt_correl;
  int npar;
  int flag_rho;
  double rho;
  double params[4];
  double modif[16];
} Local_CorPgs;

typedef struct
{
  int flag_trace;
  int idir;
  int ipas;
  int nrow;
  int ncol;
  VectorDouble trace;
} Local_TracePgs;

typedef struct
{
  Db *db;
  mutable const Rule *rule;
  PropDef *propdef;
  int flag_stat;
  int flag_facies;
  ECalcVario calcul_type;
  int igrfcur;
  int idircur;
  int ipascur;
  int ngrf;
  int npair;
  int nfacies;
  int ifirst;
  int ilast;
  VectorDouble d0;
  VectorDouble d1;
  VectorDouble memint;
  VectorDouble stat_proba;
  VectorDouble stat_thresh;
  Local_CorPgs corpgs;
  Local_TracePgs tracepgs;
  Model *model;
  Vario *vario;
  Vario *varioind;
  Vario_Order *vorder;
} Local_Pgs;

#define VARS(ivar,jvar) (vario->vars[(ivar) * vario->getNVar() + (jvar)])
#define COVS(ivar,jvar) (covs[(ivar) * nfacies + (jvar)])
#define MEMINT(ipair)   (local_pgs->memint[ipair])
#define STAT_PROBA(i,j) (M_R(local_pgs->stat_proba,local_pgs->nfacies,i,j))
#define STAT_THRESH(ifac,igrf,rank) (local_pgs->stat_thresh[2*(nfacies * (igrf) + (ifac))+(rank)])
#define LAG_USED(idir,ipas) (vario->getSwByIndex(idir,vario->getLagNumber(idir) + ipas + 1) > 0 && \
                         vario->getUtilizeByIndex(idir,vario->getLagNumber(idir) + ipas + 1))
#define TABOUT(i,j)      tabout[(j)*neq+(i)]
#define EIGVEC(i,j)      eigvec[(i)*neq+(j)]
#define RULES(ir,i)     (rules[(ir)  * NRULE  + (i)])
#define RULES1(ir,i)    (rules1[(ir) * NRULE  + (i)])
#define RULES2(ir,i)    (rules2[(ir) * NRULE  + (i)])
#define DIVS(is,i)      (divs[(is)   * ncur   + (i)])
#define FIPOSAD(ir,i)   ((ir) * NCOLOR + (i))

#define QUANT_DIR 10000
#define F(i,j) (st_index(i,j))

static double EPS = 1.e-05;
static double GS_TOLSTOP = 5.e-02;
static double GS_TOLSTOP_RHO = 1.e-01;
static CTables *CTABLES = NULL;
static bool TEST_DISCRET = false;
static int NCOLOR = 0;
static int NGRF = 0;
static int NRULE = 0;
static int BASE = 0;

// Needed declarations due to intricated recursions
static Relem* st_relem_free(Relem *relem);
static void st_relem_explore(Relem *relem, int verbose);
/*! \endcond */

/****************************************************************************
 **
 ** PURPOSE: Define if the calculations must be performed using
 ** PURPOSE: the Discretized version or not
 **
 ** IN_ARGS:  flag_discret  : Flag for the discrete option
 **
 *****************************************************************************/
void set_test_discrete(bool flag_discret)
{
  TEST_DISCRET = flag_discret;
}

/****************************************************************************
 **
 ** FUNCTION: st_relem_alloc
 **
 ** PURPOSE:  Allocate a new Relem structure
 **
 ** IN_ARGS:  old_split  : Pointer to the calling Split
 **
 *****************************************************************************/
static Relem* st_relem_alloc(Split *old_split)

{
  Relem *relem;

  /* Initializations of the New Relem structure */

  relem = new Relem;
  relem->nsplit = 0;
  relem->nrule = 0;
  relem->nbyrule = 0;
  relem->facies = VectorInt();
  relem->Rrules = nullptr;
  relem->Rfipos = nullptr;
  relem->old_split = old_split;
  return (relem);
}

/****************************************************************************
 **
 ** FUNCTION: st_split_alloc
 **
 ** PURPOSE:  Allocate a new Split structure
 **
 ** IN_ARGS:  old_relem : Pointer to the calling Relem
 **
 *****************************************************************************/
static Split* st_split_alloc(Relem *old_relem)

{
  Split *split;

  /* Initializations of the New Split structure */

  split = new Split;
  split->oper = 0;
  split->nrule = 0;
  split->nbyrule = 0;
  split->old_relem = old_relem;
  split->Srules = nullptr;
  split->Sfipos = nullptr;
  split->relems.resize(2);
  for (int i = 0; i < 2; i++)
    split->relems[i] = nullptr;
  return (split);
}

/****************************************************************************
 **
 ** FUNCTION: st_define_fipos
 **
 ** PURPOSE:  Define the position as a result of the operator and the position
 **
 ** IN_ARGS:  oper  : Rank of the operator (starting from 1)
 ** IN_ARGS:  side  : Side of the operand (1 for left and 0 for right)
 **
 *****************************************************************************/
static int st_define_fipos(int oper, int side)
{
  int reponse;

  if (IFFFF(side))
    reponse = 1;
  else
    reponse = 2 * (oper - 1) + (1 - side) + 1;

  return (reponse);
}

/****************************************************************************
 **
 ** FUNCTION: st_relem_define
 **
 ** PURPOSE:  Define the list of facies in the current Relem structure
 ** 
 *****************************************************************************/
static void st_relem_define(Relem *relem,
                            int nfacies,
                            const VectorInt &facies,
                            int side,
                            int *poss)
{
  int ecr, number;

  if (relem == (Relem*) NULL) return;

  if (poss == nullptr)
    number = nfacies;
  else
  {
    number = 0;
    for (int i = 0; i < nfacies; i++)
      if (poss[i] == side) number++;
  }

  relem->facies.resize(number, 0);
  relem->Rfipos = (int*) mem_alloc(sizeof(int) * NCOLOR, 1);
  for (int i = 0; i < NCOLOR; i++)
    relem->Rfipos[i] = 0;

  ecr = 0;
  for (int i = 0; i < nfacies; i++)
  {
    if (poss == nullptr || poss[i] == side) relem->facies[ecr++] = facies[i];
  }

  if (number == 1) relem->Rfipos[relem->facies[0] - 1] = 1;
}

/****************************************************************************
 **
 ** FUNCTION: st_rule_print
 **
 *****************************************************************************/
static void st_rule_print(int rank,
                          int nbyrule,
                          int *rules,
                          int *fipos,
                          bool flag_rank,
                          int flag_similar,
                          int flag_igrf,
                          double score)
{
  int value, iscore, loc0, loc1;

  // Print the Rank (optional)

  if (flag_rank) message("%4d:", rank + 1);

  // Print the Rule

  for (int ic = 0; ic < nbyrule; ic++)
  {
    value = RULES(rank, ic);
    if (value == 1001)
      message("  S");
    else if (value == 1002)
      message("  T");
    else
      message(" %2d", value);
  }

  // Print the score (if available)

  if (!FFFF(score))
  {
    iscore = (int) score;
    message(" -> %d", iscore);
  }

  // Print the Facies rank

  message(" (");
  for (int ic = 0; ic < NCOLOR; ic++)
    message(" %3d", fipos[FIPOSAD(rank, ic)]);
  message(" )");

  // Print the similar score

  if (flag_similar >= 0)
  {
    message(" [Same as: %3d (", flag_similar + 1);

    loc0 = flag_igrf;
    for (int igrf = 0; igrf < NGRF; igrf++)
    {
      loc1 = loc0 / 2;
      if (loc0 - 2 * loc1 > 0) message(" by Symmetry around G%d", igrf + 1);
      loc0 = loc1;
    }
    message(" )]");
  }

  // Print the end-of-line

  message("\n");
}

/****************************************************************************
 **
 ** FUNCTION: st_rules_print
 **
 *****************************************************************************/
static void st_rules_print(const char *title,
                           int nrule,
                           int nbyrule,
                           int *rules,
                           int *fipos)
{
  if (nrule <= 0) return;
  message("%s (Nrule=%d, Nbyrule=%d):\n", title, nrule, nbyrule);
  for (int ir = 0; ir < nrule; ir++)
    st_rule_print(ir, nbyrule, rules, fipos, false, -1, -1, TEST);
}

/****************************************************************************
 **
 ** FUNCTION: st_relem_subdivide
 **
 ** PURPOSE:  Define and store all the possibilities of splitting
 ** PURPOSE:  the list of facies of the current Relem
 **
 ** IN_ARGS: relem0  Pointer on the current Relem structure
 ** IN_ARGS: half    1 to take the symmetry into account
 ** IN_ARGS: noper   Number of underlying GRF
 **
 *****************************************************************************/
static void st_relem_subdivide(Relem *relem0, int half, int noper)
{
  Split *split;
  int *divs, ndiv, number, ncur, previous_oper, verbose;

  verbose = 0;
  ncur = static_cast<int>(relem0->facies.size());
  if (ncur <= 1) return;

  previous_oper = 1;
  if (relem0->old_split != NULL) previous_oper = relem0->old_split->oper;

  divs = ut_split_into_two(ncur, half, verbose, &ndiv);
  number = ndiv * noper;
  divs = (int*) mem_free((char* ) divs);

  relem0->splits.resize(number);

  number = 0;
  for (int oper = 1; oper <= noper; oper++)
  {
    half = (oper == previous_oper);
    divs = ut_split_into_two(ncur, half, verbose, &ndiv);
    for (int is = 0; is < ndiv; is++, number++)
    {
      relem0->splits[number] = split = st_split_alloc(relem0);
      split->oper = oper;
      for (int i = 0; i < 2; i++)
      {
        split->relems[i] = st_relem_alloc(split);
        st_relem_define(split->relems[i], ncur, relem0->facies, 1 - i,
                        &DIVS(is, 0));
        st_relem_subdivide(split->relems[i], 0, NGRF);
      }
    }
    divs = (int*) mem_free((char* ) divs);
  }

  relem0->splits.resize(number);
  relem0->nsplit = number;
}

/****************************************************************************
 **
 ** FUNCTION: st_split_free
 **
 ** PURPOSE:  Deallocate the Split structure
 **
 ** IN_ARGS:  split : Address of the calling Split
 **
 *****************************************************************************/
static Split* st_split_free(Split *split)
{
  if (split == (Split*) NULL) return (split);

  /* Free the descending substructures (relem) */

  for (int i = 0; i < 2; i++)
    split->relems[i] = st_relem_free(split->relems[i]);

  /* Free the local arrays */

  split->Srules = (int*) mem_free((char* ) split->Srules);
  split->Sfipos = (int*) mem_free((char* ) split->Sfipos);

  /* Free the Split structure itself */

  delete split;
  split = nullptr;

  return (split);
}

/****************************************************************************
 **
 ** FUNCTION: st_relem_free
 **
 ** PURPOSE:  Deallocate the Relem structure
 **
 ** IN_ARGS:  relem : Address of the newly deallocated Relem
 **
 *****************************************************************************/
static Relem* st_relem_free(Relem *relem)
{
  if (relem == NULL) return (relem);

  /* Free the descending substructures (Split) */

  for (int is = 0; is < relem->nsplit; is++)
    relem->splits[is] = st_split_free(relem->splits[is]);

  /* Free the local arrays */

  relem->Rrules = (int*) mem_free((char* ) relem->Rrules);
  relem->Rfipos = (int*) mem_free((char* ) relem->Rfipos);

  /* Free the Relem structure itself */

  delete relem;
  relem = nullptr;

  return (relem);
}

/****************************************************************************/
/*!
 **  Establish the array of variances
 **
 ** \param[in,out]  vario        Vario structure for the GRFs to be filled
 ** \param[in]      rule         Lithotype Rule definition
 ** \param[in]      ngrf         Number of GRFs
 **
 *****************************************************************************/
static void st_variogram_define_vars(Vario *vario, const Rule *rule, int ngrf)
{
  int igrf, jgrf;

  for (igrf = 0; igrf < ngrf; igrf++)
    for (jgrf = 0; jgrf < ngrf; jgrf++)
    {
      if (igrf == jgrf)
        vario->setVar(igrf, jgrf, 1.);
      else
        vario->setVar(igrf, jgrf, rule->getRho());
    }
}

/****************************************************************************/
/*!
 **  Define the bounds for Gaussian integrals and store them in relevant variables
 **
 *****************************************************************************/
static void st_set_bounds(Db *db,
                          int flag_one,
                          int ngrf,
                          int nfacies,
                          int ifac,
                          int iech,
                          double t1min,
                          double t1max,
                          double t2min,
                          double t2max)
{
  int jfac;
  if (!TEST_DISCRET)
  {
    jfac = (flag_one) ? 0 :
                        ifac;
    db->setBounds(iech, jfac, t1min, t1max);
    if (ngrf > 1)
    {
      jfac = (flag_one) ? 1 :
                          nfacies + ifac;
      db->setBounds(iech, jfac, t2min, t2max);
    }
  }
  else
  {
    jfac = (flag_one) ? 0 :
                        ifac;
    db->setIntervals(iech, jfac,
                     (double) ct_tableone_getrank_from_proba(CTABLES, t1min),
                     (double) ct_tableone_getrank_from_proba(CTABLES, t1max));
    if (ngrf > 1)
    {
      jfac = (flag_one) ? 1 :
                          nfacies + ifac;
      db->setIntervals(iech, jfac,
                       (double) ct_tableone_getrank_from_proba(CTABLES, t2min),
                       (double) ct_tableone_getrank_from_proba(CTABLES, t2max));
    }
  }
}

/****************************************************************************/
/*!
 **  Modify rho where it is needed
 **
 ** \param[in]  rho         Rho value
 **
 ** \param[out] local_pgs   Local_Pgs structure
 **
 *****************************************************************************/
static void st_set_rho(double rho, Local_Pgs *local_pgs)
{
  double rho2;
  int iech, ifac, ngrf;
  Db *db = local_pgs->db;
  PropDef *propdef = local_pgs->propdef;
  const Rule *rule = local_pgs->rule;
  int flag_stat = local_pgs->flag_stat;
  double t1min, t1max, t2min, t2max;

  local_pgs->corpgs.rho = rho;
  local_pgs->rule->setRho(rho);
  rho2 = rho * rho;
  st_variogram_define_vars(local_pgs->vario, local_pgs->rule, local_pgs->ngrf);

  /* Define the Thresholds */

  if (flag_stat)
  {
    rule->setProportions(propdef->proploc);
  }
  else
  {
    ngrf = local_pgs->ngrf;
    for (iech = 0; iech < db->getSampleNumber(); iech++)
    {
      if (!db->isActive(iech)) continue;
      ifac = (int) db->getVariable(iech, 0);
      if (rule_thresh_define(propdef, db, rule, ifac, iech, 0, 0, 0, &t1min,
                             &t1max, &t2min, &t2max)) return;
      st_set_bounds(db, 1, ngrf, local_pgs->nfacies, ifac, iech, t1min, t1max,
                    t2min, t2max);
    }
  }

  /* Update the modif matrix if necessary */

  if (local_pgs->corpgs.opt_correl == 2)
  {
    M_R(local_pgs->corpgs.modif,4,0,1) = rho;
    M_R(local_pgs->corpgs.modif,4,0,2) = rho;
    M_R(local_pgs->corpgs.modif,4,0,3) = rho2;
    M_R(local_pgs->corpgs.modif,4,1,3) = 1 - rho2;
  }
}

/****************************************************************************/
/*!
 **  Calculate the probability for two independent GRFs
 **
 ** \return  Evaluation value
 **
 ** \param[in] correl    Correlation between covariance
 ** \param[in] low       Array of lower thresholds or lower indices
 **                      (Dimension: 2)
 ** \param[in] up        Array of upper thresholds or upper indices
 **                      (Dimension: 2)
 ** \param[in] iconf     Rank for the discrete calculation
 **
 *****************************************************************************/
static double st_get_proba_ind(double correl,
                               double *low,
                               double *up,
                               int iconf)
{
  int ier, infin[2];
  double p1min, p1max, p2min, p2max, proba, err;

  double releps = 0.;
  double abseps = EPS;
  int maxpts = 8000;

  proba = TEST;

  if (!TEST_DISCRET)
  {
    if (correl == 0.)
    {
      p1min = law_cdf_gaussian(low[0]);
      p1max = law_cdf_gaussian(up[0]);
      p2min = law_cdf_gaussian(low[1]);
      p2max = law_cdf_gaussian(up[1]);
      proba = (p1max - p1min) * (p2max - p2min);
    }
    else
    {
      infin[0] = mvndst_infin(low[0], up[0]);
      infin[1] = mvndst_infin(low[1], up[1]);
      mvndst(2, low, up, infin, &correl, maxpts, abseps, releps, &err, &proba,
             &ier);
    }
  }
  else
  {
    proba = ct_tableone_calculate_by_rank(CTABLES, iconf, low, up);
  }
  return proba;
}

/****************************************************************************/
/*!
 **  Calculate the thresholds in the stationary case
 **
 ** \returns Error return code
 **
 ** \param[in]  local_pgs     Local_Pgs structure
 **
 *****************************************************************************/
static int st_calculate_thresh_stat(Local_Pgs *local_pgs)

{
  double t1min, t1max, t2min, t2max, cround;

  int nfacies = local_pgs->nfacies;
  int ngrf = local_pgs->ngrf;
  int iconf0 = 0;

  for (int ifac = 0; ifac < nfacies; ifac++)
  {
    if (rule_thresh_define(local_pgs->propdef, local_pgs->db, local_pgs->rule,
                           ifac + 1, 0, 0, 0, 0, &t1min, &t1max, &t2min,
                           &t2max)) return (1);
    if (!TEST_DISCRET)
    {
      STAT_THRESH(ifac,0,0) = t1min;
      STAT_THRESH(ifac,0,1) = t1max;
      STAT_THRESH(ifac,1,0) = t2min;
      STAT_THRESH(ifac,1,1) = t2max;
    }
    else
    {
      STAT_THRESH(ifac,0,0) = (double) ct_tableone_getrank_from_proba(CTABLES,
                                                                      t1min);
      STAT_THRESH(ifac,0,1) = (double) ct_tableone_getrank_from_proba(CTABLES,
                                                                      t1max);
      if (ngrf > 1)
      {
        STAT_THRESH(ifac,1,0) = (double) ct_tableone_getrank_from_proba(CTABLES,
                                                                        t2min);
        STAT_THRESH(ifac,1,1) = (double) ct_tableone_getrank_from_proba(CTABLES,
                                                                        t2max);
      }
      else
      {
        STAT_THRESH(ifac,1,0) = (double) ct_tableone_getrank_from_proba(CTABLES,
                                                                        -10.);
        STAT_THRESH(ifac,1,1) = (double) ct_tableone_getrank_from_proba(CTABLES,
                                                                        +10.);
      }
    }
  }

  // Verification

  double total = 0.;
  double correl = local_pgs->corpgs.rho;
  double low[2], up[2];
  if (TEST_DISCRET) iconf0 = ct_tableone_covrank(CTABLES, correl, &cround);

  for (int ifac = 0; ifac < nfacies; ifac++)
  {
    low[0] = STAT_THRESH(ifac, 0, 0);
    up[0] = STAT_THRESH(ifac, 0, 1);
    low[1] = STAT_THRESH(ifac, 1, 0);
    up[1] = STAT_THRESH(ifac, 1, 1);
    double proba = st_get_proba_ind(correl, low, up, iconf0);
    total += proba;
  }
  if (ABS(total - 1.) > EPSILON3)
    messerr(
        "In st_calculate_thresh_stat, the sum of Probabilities (%lf) is not close to 1.",
        total);
  return (0);
}

/****************************************************************************/
/*!
 **  Manage local variables for variopgs calculation (Non-stationary case)
 **
 ** \return  Error return code
 **
 ** \param[in]  mode          Type of usage
 **                           1 for allocation
 **                           0 for valuation (rule dependent)
 **                          -1 for deallocation
 ** \param[in]  ngrf          Number of grfs
 ** \param[in]  nfacies       Number of facies
 ** \param[in]  flag_one      1 for considering only the Facies at data point
 **                           0 for considering all facies
 ** \param[in]  flag_prop     1 for allocating variable for proportions
 ** \param[in]  db            Db structure
 ** \param[in]  propdef       PropDef structure
 ** \param[in]  rule          Lithotype Rule definition
 **
 *****************************************************************************/
static int st_vario_pgs_variable(int mode,
                                 int ngrf,
                                 int nfacies,
                                 int flag_one,
                                 int flag_prop,
                                 Db *db,
                                 PropDef *propdef,
                                 const Rule *rule)
{
  int number, ifac, jfac, nloop, iptr;
  double t1min, t1max, t2min, t2max;
  static bool is_prop_defined;

  // Dispatch

  number = (flag_one) ? ngrf :
                        ngrf * nfacies;
  if (db == nullptr) return 0;
  switch (mode)
  {
    case 1:

      /* The proportions (if not already present in correct number) */

      is_prop_defined = false;
      if (flag_prop && db->getProportionNumber() != nfacies)
      {
        iptr = db->addFieldsByConstant(nfacies, 0., String(), ELoc::P);
        if (iptr < 0) return (1);
        is_prop_defined = true;
      }

      /* The bounds */

      if (!TEST_DISCRET)
      {
        iptr = db->addFieldsByConstant(number, 0., "Lower", ELoc::L);
        if (iptr < 0) return (1);

        iptr = db->addFieldsByConstant(number, 0., "Upper", ELoc::U);
        if (iptr < 0) return (1);
      }
      else
      {
        iptr = db->addFieldsByConstant(number, 0., "Lower Rank", ELoc::RKLOW);
        if (iptr < 0) return (1);

        iptr = db->addFieldsByConstant(number, 0., "Upper Rank", ELoc::RKUP);
        if (iptr < 0) return (1);
      }
      break;

    case 0:

      /* Evaluate the bounds */
      /* Use dummy rho value in order to avoid discarding pairs in geometry */

      nloop = (flag_one) ? 1 :
                           nfacies;
      for (int iech = 0; iech < db->getSampleNumber(); iech++)
      {
        if (!db->isActive(iech)) continue;

        for (int i = 0; i < nloop; i++)
        {
          ifac = (flag_one) ? (int) db->getVariable(iech, 0) :
                              i;
          jfac = (flag_one) ? ifac :
                              ifac + 1;
          if (rule_thresh_define(propdef, db, rule, jfac, iech, 0, 0, 0, &t1min,
                                 &t1max, &t2min, &t2max)) return (1);

          /* Define the proportions */

          if (flag_prop) db->setProportion(iech, ifac, propdef->propmem[ifac]);

          /* Define the bounds */

          st_set_bounds(db, flag_one, ngrf, nfacies, ifac, iech, t1min, t1max,
                        t2min, t2max);
        }
      }
      break;

    case -1:

      /* Deallocation */

      if (flag_prop && is_prop_defined)
      {
        db->deleteFieldByLocator(ELoc::P);
      }
      if (!TEST_DISCRET)
      {
        db->deleteFieldByLocator(ELoc::L);
        db->deleteFieldByLocator(ELoc::U);
      }
      else
      {
        db->deleteFieldByLocator(ELoc::RKLOW);
        db->deleteFieldByLocator(ELoc::RKUP);
      }
      break;
  }
  return (0);
}

/****************************************************************************
 **
 ** FUNCTION: st_rule_encode
 **
 ** PURPOSE: Encode the set of numerical rule items into a rule sentence
 ** PURPOSE: that can be used to create a Rule structure
 **
 *****************************************************************************/
static Rule* st_rule_encode(int *string)
{
  VectorInt n_type = VectorInt(NRULE);
  VectorInt n_facs = VectorInt(NRULE);

  for (int i = 0; i < NRULE; i++)
  {
    if (string[i] <= 1000)
    {
      n_type[i] = 0;
      n_facs[i] = string[i];
    }
    else if (string[i] == 1001)
    {
      n_type[i] = 1;
      n_facs[i] = 0;
    }
    else if (string[i] == 1002)
    {
      n_type[i] = 2;
      n_facs[i] = 0;
    }
    else
      messageAbort("Unexpected rule number");
  }
  return Rule::createFromNumericalCoding(n_type, n_facs);
}

/****************************************************************************/
/*!
 **  Extract the information of the Trace
 **
 ** \return The calculated score
 **
 ** \param[in]   local_pgs   Local_Pgs structure
 **
 *****************************************************************************/
static double st_extract_trace(Local_Pgs *local_pgs)

{
  Local_TracePgs *tracepgs = &local_pgs->tracepgs;
  int nrow = tracepgs->nrow;
  int ncol = tracepgs->ncol;
  if (nrow <= 0 || ncol <= 0) return TEST;

  /* Evaluate the sum of the score */

  double totsum = 0.;
  for (int irow = 0; irow < nrow; irow++)
  {
    totsum += tracepgs->trace[ncol * irow + 2];
    if (ncol >= 5) totsum += tracepgs->trace[ncol * irow + 4];
  }
  set_keypair("vario.pgs_score", 1, 1, 1, &totsum);

  if (tracepgs->flag_trace)
    set_keypair("vario.pgs_stack", 1, nrow, ncol, tracepgs->trace.data());

  return (totsum);
}

/****************************************************************************/
/*!
 **  Patch the central value (dist=0) of the covariances
 **
 ** \param[in]  local_pgs   Local_Pgs structure
 ** \param[in]  vario       Vario structure for the GRFs to be filled
 ** \param[in]  idir        Rank of the direction
 ** \param[in]  rho         Correlation coefficient
 **
 *****************************************************************************/
static void st_variogram_patch_C00(Local_Pgs *local_pgs,
                                   Vario *vario,
                                   int idir,
                                   double rho)
{
  Db *db = local_pgs->db;
  int nech = (db == nullptr) ? 0 :
                               db->getActiveSampleNumber();
  vario->patchCenter(idir, nech, rho);
}

/****************************************************************************/
/*!
 **  Add a new row to the trace
 **
 ** \param[in]  local_pgs  Local_Pgs structure
 **
 *****************************************************************************/
static void trace_add_row(Local_Pgs *local_pgs)
{
  Local_TracePgs *tracepgs;
  int nrow, ncol, iad;

  /* Initializations */

  tracepgs = &local_pgs->tracepgs;
  if (!tracepgs->flag_trace) return;
  ncol = tracepgs->ncol;
  nrow = tracepgs->nrow;
  iad = ncol * nrow;

  nrow++;
  tracepgs->trace.resize(nrow * ncol);
  for (int icol = 0; icol < ncol; icol++)
    tracepgs->trace[iad + icol] = TEST;
  tracepgs->nrow = nrow;
}

/****************************************************************************/
/*!
 **  Local searching function
 **
 ** \return  Evaluation value
 **
 ** \param[in]  correl     Correlation parameter
 ** \param[in]  user_data  User Data
 **
 *****************************************************************************/
static double st_func_search_stat(double correl, void *user_data)
{
  double low[2], up[2], cround;
  Local_Pgs *local_pgs;

  /* Initializations */

  local_pgs = (Local_Pgs*) user_data;
  int iconf0 = 0;

  int nfacies = local_pgs->nfacies;
  int ipas = local_pgs->ipascur;
  int idir = local_pgs->idircur;
  int igrf = local_pgs->igrfcur;
  Vario *vario = local_pgs->varioind;

  if (TEST_DISCRET) iconf0 = ct_tableone_covrank(CTABLES, correl, &cround);

  double sum = 0.;
  for (int ifac1 = 0; ifac1 < nfacies; ifac1++)
    for (int ifac2 = 0; ifac2 < nfacies; ifac2++)
    {
      low[0] = STAT_THRESH(ifac1, igrf, 0);
      up[0] = STAT_THRESH(ifac1, igrf, 1);
      low[1] = STAT_THRESH(ifac2, igrf, 0);
      up[1] = STAT_THRESH(ifac2, igrf, 1);
      double proba = st_get_proba_ind(correl, low, up, iconf0);

      double logp = (proba <= 0.) ? -1.e30 :
                                    log(proba);
      int iad = vario->getDirAddress(idir, ifac1, ifac2, ipas, false, 1);
      double sw = vario->getSwByIndex(idir, iad);
      double gg = vario->getGgByIndex(idir, iad);
      iad = vario->getDirAddress(idir, ifac1, ifac2, ipas, false, -1);
      gg += vario->getGgByIndex(idir, iad);
      sum -= logp * gg * sw / 2.;
    }

  return (0.5 * sum);
}

/****************************************************************************/
/*!
 **  Local searching function
 **
 ** \return  Evaluation value
 **
 ** \param[in]  correl     Correlation parameter
 ** \param[in]  user_data  User Data
 **
 *****************************************************************************/
static double st_func_search_nostat(double correl, void *user_data)
{
  double low[2], up[2], proba, sum, dist, w1, w2, logp, cround;
  int ipair, i, i1, i2, ifac1, ifac2, iconf0;
  Local_Pgs *local_pgs;

  /* Initializations */

  local_pgs = (Local_Pgs*) user_data;
  iconf0 = 0;

  if (TEST_DISCRET) iconf0 = ct_tableone_covrank(CTABLES, correl, &cround);

  /* Reset the pre-calculation array (only if flag_stat) */
  if (local_pgs->flag_stat)
    for (i = 0; i < (int) local_pgs->stat_proba.size(); i++)
      local_pgs->stat_proba[i] = TEST;

  sum = 0.;
  for (ipair = local_pgs->ifirst; ipair < local_pgs->ilast; ipair++)
  {
    vario_order_get_indices(local_pgs->vorder, ipair, &i1, &i2, &dist);
    w1 = local_pgs->db->getWeight(i1);
    w2 = local_pgs->db->getWeight(i2);
    ifac1 = ifac2 = -1;
    proba = TEST;
    if (local_pgs->flag_stat)
    {

      /* In the stationary case, search in the lookup table first */

      ifac1 = (int) local_pgs->db->getVariable(i1, 0) - 1;
      if (ifac1 < 0 || ifac1 >= local_pgs->nfacies) continue;
      ifac2 = (int) local_pgs->db->getVariable(i2, 0) - 1;
      if (ifac2 < 0 || ifac2 >= local_pgs->nfacies) continue;
      proba = STAT_PROBA(ifac1, ifac2);
    }

    if (FFFF(proba))
    {
      if (!TEST_DISCRET)
      {
        low[0] = local_pgs->db->getLowerBound(i1, local_pgs->igrfcur);
        up[0] = local_pgs->db->getUpperBound(i1, local_pgs->igrfcur);
        low[1] = local_pgs->db->getLowerBound(i2, local_pgs->igrfcur);
        up[1] = local_pgs->db->getUpperBound(i2, local_pgs->igrfcur);
      }
      else
      {
        low[0] = local_pgs->db->getLowerInterval(i1, local_pgs->igrfcur);
        up[0] = local_pgs->db->getUpperInterval(i1, local_pgs->igrfcur);
        low[1] = local_pgs->db->getLowerInterval(i2, local_pgs->igrfcur);
        up[1] = local_pgs->db->getUpperInterval(i2, local_pgs->igrfcur);
      }
      proba = st_get_proba_ind(correl, low, up, iconf0);
      if (local_pgs->flag_stat)
      STAT_PROBA(ifac1,ifac2) = STAT_PROBA(ifac2,ifac1) = proba;
    }
    logp = (proba <= 0.) ? -1.e30 :
                           log(proba);
    sum -= w1 * w2 * logp;
  }
  return (0.5 * sum);
}

/****************************************************************************/
/*!
 **  Update the Trace array
 **
 ** \param[in]  local_pgs  Local_Pgs structure
 ** \param[in]  value0     First value in a Trace row
 ** \param[in]  value1     Second value in a Trace row
 ** \param[in]  origin     Origin for values in record (after 2 heading values)
 ** \param[in]  number     Number of values assigned
 ** \param[in]  values     Array of values assigned
 **
 *****************************************************************************/
static void trace_define(Local_Pgs *local_pgs,
                         double value0,
                         double value1,
                         int origin,
                         int number,
                         double *values)
{
  Local_TracePgs *tracepgs;
  int i, ncol, nrow, iad;

  /* Initializations */

  tracepgs = &local_pgs->tracepgs;
  if (!tracepgs->flag_trace) return;
  nrow = tracepgs->nrow;
  ncol = tracepgs->ncol;
  iad = ncol * (nrow - 1);
  if (2 + origin + number > ncol)
    messageAbort("Error in Trace dimension (ncol=%d origin=%d number=%d)", ncol,
                 origin, number);

  /* Store the information */

  tracepgs->trace[iad] = value0;
  tracepgs->trace[iad + 1] = value1;

  for (i = 0; i < number; i++)
    tracepgs->trace[iad + 2 + origin + i] = values[i];
}

/****************************************************************************/
/*!
 **  Performing the variogram calculations (stationary case)
 **
 ** \return  Error return code
 **
 ** \param[in]  vario         Vario structure for the GRFs to be filled
 ** \param[in]  local_pgs     Local_Pgs structure
 ** \param[in]  ngrf          Number of GRFs
 **
 *****************************************************************************/
static int st_varcalc_from_vario_stat(Vario *vario,
                                      Local_Pgs *local_pgs,
                                      int ngrf)
{
  int iad;
  double result, testval, varloc, niter;

  /* Initializations */

  st_set_rho(0., local_pgs);

  /* Loop on the directions */

  for (int idir = 0; idir < vario->getDirectionNumber(); idir++)
  {
    local_pgs->idircur = idir;

    /* Set the value of C(0) */

    st_variogram_patch_C00(local_pgs, vario, idir, 0.);

    /* Loop on the lags */

    for (int ipas = 0; ipas < vario->getLagNumber(idir); ipas++)
    {
      mes_process("Inverting Variogram Lag", vario->getLagNumber(idir), ipas);
      local_pgs->ipascur = ipas;
      trace_add_row(local_pgs);

      /* Loop on the GRFs */

      for (int igrf = 0; igrf < ngrf; igrf++)
      {
        local_pgs->igrfcur = igrf;
        result = golden_search(st_func_search_stat, (void*) local_pgs,
                               GS_TOLSTOP, -1., 1., &testval, &niter);
        trace_define(local_pgs, idir + 1, ipas + 1, 2 * igrf, 1, &testval);
        trace_define(local_pgs, idir + 1, ipas + 1, 2 * igrf + 1, 1, &niter);

        for (int jgrf = 0; jgrf <= igrf; jgrf++)
        {
          varloc = (igrf == jgrf) ? result :
                                    0.;
          iad = vario->getDirAddress(idir, igrf, jgrf, ipas, false, 1);
          vario->setGgByIndex(idir, iad, varloc);
          iad = vario->getDirAddress(idir, igrf, jgrf, ipas, false, -1);
          vario->setGgByIndex(idir, iad, varloc);

          if (debug_query("converge"))
            message("Lag:%d - Grf:%d - Variogram(%d) = %lf\n", ipas, igrf, iad,
                    varloc);
        }
      }
    }
  }
  return (0);
}

/****************************************************************************/
/*!
 **  Define the trace
 **
 ** \param[in]  flag_rho    1 if rho has to be calculated, 0 otherwise
 ** \param[in]  flag_correl 1 for the correlated case; 0 otherwise
 ** \param[in,out]  local_pgs   Local_Pgs structure
 **
 *****************************************************************************/
static void st_define_trace(int flag_rho, int flag_correl, Local_Pgs *local_pgs)
{
  Local_TracePgs *tracepgs;

  tracepgs = &local_pgs->tracepgs;
  tracepgs->flag_trace = !flag_rho;
  if (!tracepgs->flag_trace) return;

  if (!flag_correl)
    tracepgs->ncol = 2 + 2 * local_pgs->ngrf;
  else
    tracepgs->ncol = 2 + 2 + local_pgs->corpgs.npar;

  tracepgs->nrow = 0;
  tracepgs->trace.clear();
}

/****************************************************************************/
/*!
 **  Reset the Local_TracePgs structure
 **
 ** \param[in,out]  local_pgs  Local_TracePgs structure
 **
 *****************************************************************************/
static void st_retrace_define(Local_Pgs *local_pgs)

{
  Local_TracePgs *tracepgs;

  tracepgs = &local_pgs->tracepgs;
  if (!tracepgs->flag_trace) return;

  // Initialize the new trace (which clears any previously defined trace)

  st_define_trace(0, 0, local_pgs);
}

/****************************************************************************/
/*!
 **  Evaluate the variogram of one underlying GRF
 **
 ** \param[in]  local_pgs  Local_Pgs structure
 ** \param[in]   idir      Rank of the direction
 **
 *****************************************************************************/
static void st_varcalc_uncorrelated_grf(Local_Pgs *local_pgs, int idir)
{
  int ipas, iad, igrf, jgrf, ngrf;
  double result, testval, niter, varloc;
  Vario *vario;

  vario = local_pgs->vario;
  ngrf = local_pgs->ngrf;

  /* Loop on the lags */

  for (ipas = 0; ipas < vario->getLagNumber(idir); ipas++)
  {
    mes_process("Inverting Variogram Lag", vario->getLagNumber(idir), ipas);
    local_pgs->ipascur = ipas;
    trace_add_row(local_pgs);
    if (!LAG_USED(idir, ipas)) continue;
    vario_order_get_bounds(local_pgs->vorder, idir, ipas, &local_pgs->ifirst,
                           &local_pgs->ilast);
    if (local_pgs->ifirst >= local_pgs->ilast) continue;

    for (igrf = 0; igrf < ngrf; igrf++)
    {
      local_pgs->igrfcur = igrf;
      result = golden_search(st_func_search_nostat, (void*) local_pgs,
                             GS_TOLSTOP, -1., 1., &testval, &niter);
      trace_define(local_pgs, idir + 1, ipas + 1, 2 * igrf, 1, &testval);
      trace_define(local_pgs, idir + 1, ipas + 1, 2 * igrf + 1, 1, &niter);
      for (jgrf = 0; jgrf <= igrf; jgrf++)
      {
        varloc = (igrf == jgrf) ? result :
                                  0.;
        iad = vario->getDirAddress(idir, igrf, jgrf, ipas, false, 1);
        vario->setGgByIndex(idir, iad, varloc);
        iad = vario->getDirAddress(idir, igrf, jgrf, ipas, false, -1);
        vario->setGgByIndex(idir, iad, varloc);

        if (debug_query("converge"))
          message("Lag:%d - Grf:%d - Variogram(%d) = %lf\n", ipas, igrf, iad,
                  varloc);
      }
    }
  }

  return;
}

/****************************************************************************
 **
 ** FUNCTION: st_rule_calcul
 **
 ** \param[in] localpgs    Local_Pgs structure
 ** \param[in] string      Rule string
 **
 *****************************************************************************/
static double st_rule_calcul(Local_Pgs *local_pgs, int *string)
{
  double score;

  /* Preliminary assignments */

  score = 0.;
  local_pgs->rule = st_rule_encode(string);
  local_pgs->ngrf = local_pgs->rule->getGRFNumber();
  local_pgs->vario->setNVar(local_pgs->ngrf);
  // TODO : ngrf now is 1 (but vario had nvar = 2)
  // The following instruction provoques a message from Vario.cpp :
  //     Invalid dimension for 'means' (2)
  //     It should match the number of variables in 'Db' (1)
  //local_pgs->vario->internalVariableResize();
  //local_pgs->vario->internalDirectionResize();
  st_retrace_define(local_pgs);

  if (local_pgs->flag_stat)
  {
    st_set_rho(0., local_pgs);
    (void) st_calculate_thresh_stat(local_pgs);
    st_varcalc_from_vario_stat(local_pgs->vario, local_pgs, local_pgs->ngrf);
  }
  else
  {
    (void) st_vario_pgs_variable(0, local_pgs->ngrf, local_pgs->nfacies, 1, 0,
                                 local_pgs->db, local_pgs->propdef,
                                 local_pgs->rule);
    st_set_rho(0., local_pgs);
    for (int idir = 0; idir < local_pgs->vario->getDirectionNumber(); idir++)
    {
      local_pgs->idircur = idir;
      st_variogram_patch_C00(local_pgs, local_pgs->vario, idir,
                             local_pgs->rule->getRho());
      st_varcalc_uncorrelated_grf(local_pgs, idir);
    }
  }

  /* Deallocation of the Rule */

  local_pgs->rule = rule_free(local_pgs->rule);

  score = st_extract_trace(local_pgs);
  return (score);
}

/****************************************************************************
 **
 ** FUNCTION: st_permut
 **
 *****************************************************************************/
static int st_permut(int value, int igrf)
{
  if (igrf == 0)
  {
    if (value == 1) return (2);
    if (value == 2) return (1);
  }
  else if (igrf == 1)
  {
    if (value == 3) return (4);
    if (value == 4) return (3);
  }
  else if (igrf == 2)
  {
    if (value == 5) return (6);
    if (value == 6) return (5);
  }
  else
  {
    messageAbort("Function st_permut has been programmed up to 3 GRFs");
  }
  return (value);
}

/****************************************************************************
 **
 ** FUNCTION: st_fipos_encode
 **
 ** \return The encoded value
 **
 ** \param[in] fgrf      Array of codes (skipped if < 0)
 **
 *****************************************************************************/
static int st_fipos_encode(VectorInt &fgrf)
{
  int nmax, found, fipos;

  nmax = 1 + NGRF;

  found = fipos = 0;
  for (int i = nmax - 1; i >= 0; i--)
  {
    if (fgrf[i] > 0) found++;
    if (found == 1) fipos = 1;
    fipos = fipos * BASE + fgrf[i];
  }
  return (fipos);
}

/****************************************************************************
 **
 ** FUNCTION: st_fipos_decode
 **
 ** \param[in] fipos    Value to be decoded
 **
 ** \param[out] fgrf    Array of codes
 **
 *****************************************************************************/
static void st_fipos_decode(int fipos, VectorInt &fgrf)
{
  int nmax, div;

  nmax = 1 + NGRF;
  for (int i = 0; i < nmax; i++)
    fgrf[i] = 0;
  for (int i = 0; i < nmax; i++)
  {
    fipos = fipos - 1;
    div = fipos / BASE;
    if (div > 0) fgrf[i] = fipos - BASE * div + 1;
    fipos = div;
  }
}

/****************************************************************************
 **
 ** FUNCTION: st_update_orientation
 **
 *****************************************************************************/
static int st_update_orientation(int fac0, int igrf_cas, VectorInt &fgrf)
{
  int fac, facp, nmax, loc0, loc1;

  /* Preliminary check */

  nmax = 1 + NGRF;
  fac = fac0;
  if (fac0 < 0) return (fac0);

  /* Decomposition */

  st_fipos_decode(fac, fgrf);

  /* Loop on the GRF permutations */

  loc0 = igrf_cas;
  for (int igrf = 0; igrf < NGRF; igrf++)
  {
    loc1 = loc0 / 2;
    if (loc0 - 2 * loc1 > 0)
    {

      /* Update the orientation of 'grf' */

      for (int i = 0; i < nmax; i++)
        fgrf[i] = st_permut(fgrf[i], igrf);
    }
    loc0 = loc1;
  }

  /* Recomposition */

  facp = st_fipos_encode(fgrf);

  return (facp);
}

/****************************************************************************
 **
 ** FUNCTION: st_same_score
 **
 ** REMARKS In this function, we change the orientation of one or several GRFs
 **
 *****************************************************************************/
static int st_same_score(Relem *relem,
                         int ir0,
                         int igrf_cas,
                         VectorInt &fgrf,
                         VectorInt &fcmp)
{
  int *fipos, flag_same;

  fipos = relem->Rfipos;
  if (ir0 <= 0) return (-1);

  // Modify the orientation of 'grf' for the current 'fipos'

  for (int ic = 0; ic < NCOLOR; ic++)
    fcmp[ic] = st_update_orientation(fipos[FIPOSAD(ir0, ic)], igrf_cas, fgrf);

  // Look if the same 'fipos' has already been calculated

  for (int ir = 0; ir < ir0; ir++)
  {
    flag_same = 1;
    for (int ic = 0; ic < NCOLOR && flag_same; ic++)
    {
      if (fipos[FIPOSAD(ir, ic)] != fcmp[ic]) flag_same = 0;
    }
    if (flag_same) return (ir);
  }

  return (-1);
}

/****************************************************************************
 **
 ** FUNCTION: st_relem_evaluate
 **
 *****************************************************************************/
static VectorDouble st_relem_evaluate(Relem *relem,
                                      int verbose,
                                      VectorInt &fgrf,
                                      VectorInt &fcmp,
                                      Local_Pgs *local_pgs,
                                      int *nscore,
                                      int *r_opt)
{
  int *rules, *fipos;
  int nrule, indice, nmax, flag_check, igrf_cas, number, igrf_opt;
  double score_ref;
  VectorDouble scores;

  /* Initializations */

  flag_check = (int) get_keypone("Multi_Score_Check", 0.);
  nmax = (int) pow(2., (double) NGRF);
  nrule = relem->nrule;
  rules = relem->Rrules;
  fipos = relem->Rfipos;
  *nscore = nrule;

  /* Core allocation */

  scores.resize(nrule);

  *r_opt = number = 0;
  for (int ir = 0; ir < nrule; ir++)
  {

    // Check if the same flag has already been found

    indice = igrf_opt = -1;
    for (igrf_cas = 1; igrf_cas < nmax && igrf_opt < 0; igrf_cas++)
    {
      indice = st_same_score(relem, ir, igrf_cas, fgrf, fcmp);
      if (indice >= 0) igrf_opt = igrf_cas;
    }

    // Set the score 

    if (indice >= 0)
      scores[ir] = scores[indice];
    else
    {
      number++;
      scores[ir] = st_rule_calcul(local_pgs, &RULES(ir, 0));
      propdef_reset(local_pgs->propdef);
    }

    // When Multi_Score_Check, calculate the score even if already defined
    // and compare both results

    if (flag_check && indice >= 0)
    {
      score_ref = st_rule_calcul(local_pgs, &RULES(ir, 0));
      if (ABS(scores[ir] - score_ref) > 1.e-10 * score_ref)
      {
        messerr("Warning: Difference between score stored and re-evaluated:");
        messerr("- as already stored = %lf", scores[ir]);
        messerr("- as re-evaluated   = %lf", score_ref);
      }
    }

    // Optional printout

    if (verbose)
    {
      if (indice < 0)
        st_rule_print(ir, NRULE, rules, fipos, true, indice, igrf_opt,
                      scores[ir]);
    }

    // Ranking the Minimum score

    if (scores[ir] < scores[*r_opt]) *r_opt = ir;
  }

  // Store the different rules as well as the scores in keypair mechanism

  set_keypair("rule_auto_scores", 1, 1, nrule, scores.data());
  set_keypair_int("rule_auto_allrules", 1, nrule, NRULE, rules);
  set_keypair_int("rule_auto_best_rule", 1, 1, NRULE, &RULES(*r_opt, 0));

  return (scores);
}

/****************************************************************************
 **
 ** FUNCTION: st_rule_glue
 **
 *****************************************************************************/
static void st_rule_glue(Relem *relem,
                         int nrule1,
                         int nbyrule1,
                         int *rules1,
                         int *fipos1)
{
  int *rules, *fipos, nrule, ir, nnew;

  if (relem == (Relem*) NULL) return;
  if (nrule1 <= 0) return;

  nrule = ir = relem->nrule;
  nnew = nrule + nrule1;

  relem->Rrules = rules = (int*) mem_realloc((char* ) relem->Rrules,
                                             sizeof(int) * NRULE * nnew, 1);
  relem->Rfipos = fipos = (int*) mem_realloc((char* ) relem->Rfipos,
                                             sizeof(int) * NCOLOR * nnew, 1);

  for (int i1 = 0; i1 < nrule1; i1++, ir++)
  {
    for (int ic = 0; ic < nbyrule1; ic++)
      RULES(ir,ic) = RULES1(i1, ic);
    for (int ic = 0; ic < NCOLOR; ic++)
      fipos[FIPOSAD(ir, ic)] = fipos1[FIPOSAD(i1, ic)];
  }

  relem->nrule = nnew;
  relem->nbyrule = nbyrule1;
}

/****************************************************************************
 **
 ** FUNCTION: st_rule_product
 **
 *****************************************************************************/
static void st_rule_product(Split *split,
                            int nprod,
                            int nrule1,
                            int nbyrule1,
                            int *rules1,
                            int *fipos1,
                            int nrule2,
                            int nbyrule2,
                            int *rules2,
                            int *fipos2)
{
  int *rules, *fipos, ir, ic, oper;
  int flag_debug = 0;

  split->Srules = rules = (int*) mem_alloc(sizeof(int) * NRULE * nprod, 1);
  for (int i = 0; i < NRULE * nprod; i++)
    rules[i] = 0;
  split->Sfipos = fipos = (int*) mem_alloc(sizeof(int) * NCOLOR * nprod, 1);
  for (int i = 0; i < NCOLOR * nprod; i++)
    fipos[i] = 0;

  ir = 0;
  for (int i1 = 0; i1 < nrule1; i1++)
    for (int i2 = 0; i2 < nrule2; i2++, ir++)
    {
      ic = 0;
      oper = split->oper;
      RULES(ir,ic++) = 1000 + oper;

      if (flag_debug)
      {
        message("Rule Product (with operator %d)\n", oper);
        st_rule_print(i1, nbyrule1, rules1, fipos1, false, -1, -1, TEST);
        st_rule_print(i2, nbyrule2, rules2, fipos2, false, -1, -1, TEST);
      }

      for (int i = 0; i < nbyrule1; i++)
        RULES(ir,ic++) = RULES1(i1, i);
      for (int i = 0; i < nbyrule2; i++)
        RULES(ir,ic++) = RULES2(i2, i);
      for (int i = 0; i < NCOLOR; i++)
      {
        if (fipos1[FIPOSAD(i1, i)] > 0)
          fipos[FIPOSAD(ir, i)] = fipos1[FIPOSAD(i1, i)] * BASE
              + st_define_fipos(oper, 1);
        if (fipos2[FIPOSAD(i2, i)] > 0)
          fipos[FIPOSAD(ir, i)] = fipos2[FIPOSAD(i2, i)] * BASE
              + st_define_fipos(oper, 0);
      }

      if (flag_debug)
      {
        message("Product result=");
        st_rule_print(ir, nbyrule1 + nbyrule2 + 1, rules, fipos, false, -1, -1,
                      TEST);
      }
    }
  split->nrule = nprod;
  split->nbyrule = nbyrule1 + nbyrule2 + 1;
}

/****************************************************************************
 **
 ** FUNCTION: st_split_collapse
 **
 *****************************************************************************/
static void st_split_collapse(Split *split, int verbose)

{
  Relem *relem;
  int num[2], nby[2], *ptr[2], nprod;

  if (split == (Split*) NULL) return;

  // Explore the two Relems

  for (int i = 0; i < 2; i++)
    st_relem_explore(split->relems[i], verbose);

  // Prepare collapsing

  if (split->nrule <= 0)
  {
    for (int i = 0; i < 2; i++)
    {
      relem = split->relems[i];
      if (relem->facies.size() <= 1)
      {
        num[i] = 1;
        nby[i] = 1;
        ptr[i] = &relem->facies[0];
      }
      else
      {
        num[i] = relem->nrule;
        nby[i] = relem->nbyrule;
        ptr[i] = relem->Rrules;
      }
    }

    // Merge the rules of the two Relem (by product)

    nprod = num[0] * num[1];
    split->nbyrule = nby[0] + nby[1] + 1;
    if (nprod > 0)
    {
      st_rule_product(split, nprod, num[0], nby[0], ptr[0],
                      split->relems[0]->Rfipos, num[1], nby[1], ptr[1],
                      split->relems[1]->Rfipos);
      if (verbose)
        st_rules_print("Split", split->nrule, split->nbyrule, split->Srules,
                       split->Sfipos);
    }
  }
}

/****************************************************************************
 **
 ** FUNCTION: st_relem_explore
 **
 *****************************************************************************/
static void st_relem_explore(Relem *relem, int verbose)
{
  Split *split;

  if (relem == (Relem*) NULL) return;

  for (int is = 0; is < relem->nsplit; is++)
  {
    split = relem->splits[is];
    st_split_collapse(split, verbose);
    st_rule_glue(relem, split->nrule, split->nbyrule, split->Srules,
                 split->Sfipos);
    if (verbose)
      st_rules_print("Relem", relem->nrule, relem->nbyrule, relem->Rrules,
                     relem->Rfipos);
  }
}

/****************************************************************************/
/*!
 **  Manage the Variogram Order structure
 **
 ** \return Pointer to the Vario_Order structure
 **
 ** \param[in]  mode        Usage:
 ** \li                     1 : to initialize
 ** \li                     0 : to clean the geometry
 ** \li                    -1 : to delete
 ** \param[in]  flag_dist   1 if distances are stored; 0 otherwise
 ** \param[in]  size_aux    Size (in bytes) of the auxiliary array
 ** \param[in]  vorder      Vario_Order structure
 **
 *****************************************************************************/
Vario_Order* vario_order_manage(int mode,
                                                int flag_dist,
                                                int size_aux,
                                                Vario_Order *vorder)
{
  Vario_Order *vorder_loc;

  /* Dispatch */

  vorder_loc = (Vario_Order*) NULL;
  switch (mode)
  {
    case 1:
      vorder_loc = (Vario_Order*) mem_alloc(sizeof(Vario_Order), 0);
      if (vorder_loc == (Vario_Order*) NULL) return (vorder_loc);
      vorder_loc->npair = 0;
      vorder_loc->nalloc = 0;
      vorder_loc->size_aux = size_aux;
      vorder_loc->flag_dist = flag_dist;
      vorder_loc->tab_iech = nullptr;
      vorder_loc->tab_jech = nullptr;
      vorder_loc->tab_ipas = nullptr;
      vorder_loc->tab_sort = nullptr;
      vorder_loc->tab_aux_iech = nullptr;
      vorder_loc->tab_aux_jech = nullptr;
      vorder_loc->tab_dist = nullptr;
      break;

    case 0:
      vorder_loc = vorder;
      if (vorder == (Vario_Order*) NULL) return (vorder_loc);
      vorder_loc->tab_iech = (int*) mem_free((char* ) vorder_loc->tab_iech);
      vorder_loc->tab_jech = (int*) mem_free((char* ) vorder_loc->tab_jech);
      vorder_loc->tab_ipas = (int*) mem_free((char* ) vorder_loc->tab_ipas);
      vorder_loc->tab_sort = (int*) mem_free((char* ) vorder_loc->tab_sort);
      vorder_loc->tab_aux_iech = (char*) mem_free(
          (char* ) vorder_loc->tab_aux_iech);
      vorder_loc->tab_aux_jech = (char*) mem_free(
          (char* ) vorder_loc->tab_aux_jech);
      vorder_loc->tab_dist = (double*) mem_free((char* ) vorder_loc->tab_dist);
      break;

    case -1:
      vorder_loc = vorder;
      if (vorder == (Vario_Order*) NULL) return (vorder_loc);
      if (vorder_loc != nullptr)
      {
        vorder_loc->tab_iech = (int*) mem_free((char* ) vorder_loc->tab_iech);
        vorder_loc->tab_jech = (int*) mem_free((char* ) vorder_loc->tab_jech);
        vorder_loc->tab_ipas = (int*) mem_free((char* ) vorder_loc->tab_ipas);
        vorder_loc->tab_sort = (int*) mem_free((char* ) vorder_loc->tab_sort);
        vorder_loc->tab_aux_iech = (char*) mem_free(
            (char* ) vorder_loc->tab_aux_iech);
        vorder_loc->tab_aux_jech = (char*) mem_free(
            (char* ) vorder_loc->tab_aux_jech);
        vorder_loc->tab_dist = (double*) mem_free(
            (char* ) vorder_loc->tab_dist);
        vorder_loc = (Vario_Order*) mem_free((char* ) vorder_loc);
      }
      break;
  }
  return (vorder_loc);
}

/****************************************************************************/
/*!
 **  Add a record to the Variogram Order structure
 **
 ** \return Error return code
 **
 ** \param[in]  vorder      Vario_Order structure
 ** \param[in]  iech        Rank of the first sample
 ** \param[in]  jech        Rank of the second sample
 ** \param[in]  aux_iech    Auxiliary array for sample 'iech' (or NULL)
 ** \param[in]  aux_jech    Auxiliary array for sample 'jech' (or NULL)
 ** \param[in]  ipas        Rank of the lag
 ** \param[in]  idir        Rank of the direction (or 0)
 ** \param[in]  dist        Calculated distance (only stored if flag_dist == 1)
 **
 *****************************************************************************/
int vario_order_add(Vario_Order *vorder,
                                    int iech,
                                    int jech,
                                    void *aux_iech,
                                    void *aux_jech,
                                    int ipas,
                                    int idir,
                                    double dist)
{
  int iad;
  static int VARIO_ORDER_QUANT = 1000;

  if (vorder == (Vario_Order*) NULL) return (0);

  /* Resize the array */

  if (vorder->npair >= vorder->nalloc)
  {
    vorder->nalloc += VARIO_ORDER_QUANT;
    vorder->tab_iech = (int*) mem_realloc((char* ) vorder->tab_iech,
                                          vorder->nalloc * sizeof(int), 0);
    if (vorder->tab_iech == nullptr) return (1);
    vorder->tab_jech = (int*) mem_realloc((char* ) vorder->tab_jech,
                                          vorder->nalloc * sizeof(int), 0);
    if (vorder->tab_jech == nullptr) return (1);
    vorder->tab_ipas = (int*) mem_realloc((char* ) vorder->tab_ipas,
                                          vorder->nalloc * sizeof(int), 0);
    if (vorder->tab_ipas == nullptr) return (1);
    vorder->tab_sort = (int*) mem_realloc((char* ) vorder->tab_sort,
                                          vorder->nalloc * sizeof(int), 0);
    if (vorder->tab_sort == nullptr) return (1);
    if (vorder->size_aux > 0)
    {
      vorder->tab_aux_iech = (char*) mem_realloc(
          (char* ) vorder->tab_aux_iech, vorder->nalloc * vorder->size_aux, 0);
      if (vorder->tab_aux_iech == nullptr) return (1);
      vorder->tab_aux_jech = (char*) mem_realloc(
          (char* ) vorder->tab_aux_jech, vorder->nalloc * vorder->size_aux, 0);
      if (vorder->tab_aux_jech == nullptr) return (1);
    }
    if (vorder->flag_dist)
    {
      vorder->tab_dist = (double*) mem_realloc((char* ) vorder->tab_dist,
                                               vorder->nalloc * sizeof(double),
                                               0);
      if (vorder->tab_dist == nullptr) return (1);
    }
  }

  /* Add the new information */

  vorder->tab_iech[vorder->npair] = (dist > 0) ? iech :
                                                 jech;
  vorder->tab_jech[vorder->npair] = (dist > 0) ? jech :
                                                 iech;
  vorder->tab_ipas[vorder->npair] = ipas + idir * QUANT_DIR;
  if (vorder->flag_dist) vorder->tab_dist[vorder->npair] = dist;
  if (vorder->size_aux > 0)
  {
    iad = vorder->npair * vorder->size_aux;
    if (aux_iech != NULL)
      (void) memcpy(&vorder->tab_aux_iech[iad], aux_iech, vorder->size_aux);
    if (aux_jech != NULL)
      (void) memcpy(&vorder->tab_aux_jech[iad], aux_jech, vorder->size_aux);
  }
  vorder->npair++;
  return (0);
}

/****************************************************************************/
/*!
 **  Print the Vario_Order structure
 **
 ** \param[in]  vorder      Vario_Order structure
 ** \param[in]  idir_target Rank of the target direction (starting from 0) or -1
 ** \param[in]  ipas_target Rank of the target lag (starting from 0) or -1
 ** \param[in]  verbose     1 for a complete printout
 **
 *****************************************************************************/
void vario_order_print(Vario_Order *vorder,
                                       int idir_target,
                                       int ipas_target,
                                       int verbose)
{
  int i, j, ipas, idir, flag_first;

  if (vorder == (Vario_Order*) NULL) return;

  mestitle(0, "Variogram Order structure");
  message("Allocated size    = %d\n", vorder->nalloc);
  message("Number of pairs   = %d\n", vorder->npair);
  if (!verbose) return;
  flag_first = 1;

  for (i = 0; i < vorder->npair; i++)
  {
    j = (vorder->tab_sort == nullptr) ? i :
                                        vorder->tab_sort[i];
    ipas = vorder->tab_ipas[j];
    idir = ipas / QUANT_DIR;
    ipas = ipas - QUANT_DIR * idir;
    if (idir_target >= 0 && idir != idir_target) continue;
    if (ipas_target >= 0 && ipas != ipas_target) continue;

    if (flag_first)
    {
      if (!vorder->flag_dist)
        message("Rank - Dir - Lag - I - J\n");
      else
        message("Rank - Dir - Lag - I - J - Dist\n");
      flag_first = 0;
    }

    message("%5d", i + 1);
    message(" %5d", idir + 1);
    message(" %5d", ipas + 1);
    message(" %5d", vorder->tab_iech[j] + 1);
    message(" %5d", vorder->tab_jech[j] + 1);
    if (vorder->flag_dist) message(" %lf", vorder->tab_dist[j]);
    message("\n");
  }
}

/****************************************************************************/
/*!
 **  Resize the array and sort it
 **
 ** \return Pointer to the Vario_Order structure
 **
 ** \param[in]  vorder      Vario_Order structure
 ** \param[in]  npair       Final number of pairs
 **
 *****************************************************************************/
Vario_Order* vario_order_final(Vario_Order *vorder, int *npair)
{
  int i, error;

  *npair = 0;
  if (vorder == (Vario_Order*) NULL) return (vorder);

  error = 0;
  if (vorder->npair > 0)
  {
    vorder->tab_iech = (int*) mem_realloc((char* ) vorder->tab_iech,
                                          vorder->npair * sizeof(int), 0);
    if (vorder->tab_iech == nullptr) error = 1;
    vorder->tab_jech = (int*) mem_realloc((char* ) vorder->tab_jech,
                                          vorder->npair * sizeof(int), 0);
    if (vorder->tab_jech == nullptr) error = 1;
    vorder->tab_ipas = (int*) mem_realloc((char* ) vorder->tab_ipas,
                                          vorder->npair * sizeof(int), 0);
    if (vorder->tab_ipas == nullptr) error = 1;
    vorder->tab_sort = (int*) mem_realloc((char* ) vorder->tab_sort,
                                          vorder->npair * sizeof(int), 0);
    if (vorder->tab_sort == nullptr) error = 1;
    if (vorder->flag_dist)
    {
      vorder->tab_dist = (double*) mem_realloc((char* ) vorder->tab_dist,
                                               vorder->npair * sizeof(double),
                                               0);
      if (vorder->tab_dist == nullptr) error = 1;
    }
    if (vorder->size_aux > 0)
    {
      vorder->tab_aux_iech = (char*) mem_realloc(
          (char* ) vorder->tab_aux_iech, vorder->npair * vorder->size_aux, 0);
      if (vorder->tab_aux_iech == nullptr) error = 1;
      vorder->tab_aux_jech = (char*) mem_realloc(
          (char* ) vorder->tab_aux_jech, vorder->npair * vorder->size_aux, 0);
      if (vorder->tab_aux_iech == nullptr) error = 1;
    }
  }
  vorder->nalloc = vorder->npair;

  if (error)
  {
    vorder = vario_order_manage(-1, vorder->flag_dist, vorder->size_aux,
                                vorder);
    *npair = 0;
  }
  else if (vorder->npair > 0)
  {
    for (i = 0; i < vorder->npair; i++)
      vorder->tab_sort[i] = i;
    ut_sort_int(1, vorder->npair, vorder->tab_sort, vorder->tab_ipas);
    *npair = vorder->npair;
  }
  return (vorder);
}

/****************************************************************************/
/*!
 **  Returns the two samples for a given (ordered) pair
 **
 ** \param[in]  vorder      Vario_Order structure
 ** \param[in]  ipair       Rank of the sorted pair
 **
 ** \param[out] iech        Rank of the first sample
 ** \param[out] jech        Rank of the second sample
 ** \param[out] dist        Calculated distance or TEST (if flag_dist == 0)
 **
 *****************************************************************************/
void vario_order_get_indices(Vario_Order *vorder,
                                             int ipair,
                                             int *iech,
                                             int *jech,
                                             double *dist)
{
  int jpair;

  if (vorder->tab_sort == nullptr) messageAbort("vario_order_get_indices");
  jpair = vorder->tab_sort[ipair];
  *iech = vorder->tab_iech[jpair];
  *jech = vorder->tab_jech[jpair];
  *dist = (vorder->flag_dist) ? vorder->tab_dist[jpair] :
                                TEST;
}

/****************************************************************************/
/*!
 **  Returns the two auxiliary arrays for a given (ordered) pair
 **
 ** \param[in]  vorder      Vario_Order structure
 ** \param[in]  ipair       Rank of the sorted pair
 **
 ** \param[out] aux_iech    Array to auxiliary information for sample 'iech'
 ** \param[out] aux_jech    Array to auxiliary information for sample 'jech'
 **
 *****************************************************************************/
void vario_order_get_auxiliary(Vario_Order *vorder,
                                               int ipair,
                                               char *aux_iech,
                                               char *aux_jech)
{
  int jpair, iad;

  if (vorder->tab_sort == nullptr) messageAbort("vario_order_get_auxiliary");
  jpair = vorder->tab_sort[ipair];
  iad = vorder->size_aux * jpair;
  (void) memcpy(aux_iech, &vorder->tab_aux_iech[iad], vorder->size_aux);
  (void) memcpy(aux_jech, &vorder->tab_aux_jech[iad], vorder->size_aux);
}

/****************************************************************************/
/*!
 **  Returns the first and last indices matching a target lag
 **
 ** \param[in]  vorder      Vario_Order structure
 ** \param[in]  idir        Rank of the target direction
 ** \param[in]  ipas        Rank of the target lag
 **
 ** \param[out] ifirst      Rank of the first sample of the lag (included)
 ** \param[out] ilast       Rank of the last sample of the lag (excluded)
 **
 *****************************************************************************/
void vario_order_get_bounds(Vario_Order *vorder,
                                            int idir,
                                            int ipas,
                                            int *ifirst,
                                            int *ilast)
{
  int ipair, jpair, ival;

  ival = ipas + idir * QUANT_DIR;
  if (vorder->npair > 0 && vorder->tab_sort == nullptr)
    messageAbort("vario_order_get_bounds");
  *ifirst = vorder->npair;
  *ilast = -1;
  for (ipair = 0; ipair < vorder->npair; ipair++)
  {
    jpair = vorder->tab_sort[ipair];
    if (vorder->tab_ipas[jpair] == ival)
    {
      if (ipair < *ifirst) *ifirst = ipair;
    }
    else
    {
      if (ipair > *ifirst)
      {
        *ilast = ipair;
        return;
      }
    }
  }

  /* Particular case of the last lag */

  if (*ifirst < vorder->npair)
  {
    *ilast = vorder->npair;
    return;
  }
  return;
}

/****************************************************************************/
/*!
 **  Calculate the generalized inverse of a square symmetric matrix
 **
 ** \return  Error returned code
 **
 ** \param[in]  a         Matrix to be inverted
 ** \param[in]  neq       Number of equations
 **
 ** \param[out] tabout    Inverted matrix
 **
 *****************************************************************************/
static int invgen(double *a, int neq, double *tabout)
{
  double *eigvec, *eigval, value;
  int i, j, k, error;

  /* Initializations */

  error = 1;
  eigvec = eigval = nullptr;

  /* Core allocation */

  eigval = (double*) mem_alloc(sizeof(double) * neq, 0);
  if (eigval == nullptr) goto label_end;
  eigvec = (double*) mem_alloc(sizeof(double) * neq * neq, 0);
  if (eigvec == nullptr) goto label_end;

  /* Calculate the eigen vectors */

  if (matrix_eigen(a, neq, eigval, eigvec)) goto label_end;

  /* Calculate the generalized inverse */

  for (i = 0; i < neq; i++)
    for (j = 0; j < neq; j++)
    {
      value = 0.;
      for (k = 0; k < neq; k++)
      {
        if (ABS(eigval[k]) > 1e-10) value += EIGVEC(k,i)* EIGVEC(k,j) / eigval[k];
      }
      TABOUT(i,j)= value;
    }

    /* Set the error returned code */

  error = 0;

  label_end: eigval = (double*) mem_free((char* ) eigval);
  eigvec = (double*) mem_free((char* ) eigvec);
  return (error);
}

/****************************************************************************/
/*!
 ** Calculate the indexes of each parameter
 **
 ** \param[in]  i Index
 ** \param[in]  j Index
 **
 *****************************************************************************/

static int st_index(int i, int j)
{
  int value;

  value = 0;
  switch (i)
  {
    case 0:
      value = (j == 0) ? 1 :
                         0;
      break;

    case 1:
      value = (j == 0) ? 3 :
                         0;
      break;

    case 2:
      value = (j == 0) ? 2 :
                         1;
      break;

    case 3:
      value = (j == 0) ? 3 :
                         2;
      break;

  }
  return (value);
}

/****************************************************************************/
/*!
 **  Establish the total vector C1(h), C12(h), C21(h), C2(h)
 ** \param[in]  corpgs      Local_CorPgs structure
 ** \param[in]  params_in   Parameters (Dimension corpgs.npar)
 **
 ** \param[out] params      Parameters (Dimension = 4)
 **
 *****************************************************************************/
static void st_compute_params(Local_CorPgs *corpgs,
                              double *params_in,
                              double *params)
{
  double rho, rho2;
  rho = corpgs->rho;

  switch (corpgs->opt_correl)
  {
    case 0:
      params[0] = params_in[0];
      params[1] = params_in[1];
      params[2] = params_in[2];
      params[3] = params_in[3];
      break;

    case 1: /* Symmetrical case */
      params[0] = params_in[0];
      params[1] = params[2] = params_in[1];
      params[3] = params_in[2];
      break;

    case 2: /* Residual case */
      rho2 = rho * rho;
      params[0] = params_in[0];
      params[1] = params[2] = rho * params_in[0];
      params[3] = rho2 * params_in[0] + (1. - rho2) * params_in[1];
      break;
  }
}

/****************************************************************************/
/*!
 **  Establish the correlation for C1(h), C12(h), C21(h), C2(h)
 **
 ** \param[in]  corpgs      Local_CorPgs structure
 ** \param[in]  params_in   Array of parameters
 **
 ** \param[out] correl      Correlation matrix (Dimension = 4*4)
 **
 *****************************************************************************/
static void st_build_correl(Local_CorPgs *corpgs,
                            double *params_in,
                            double *correl)
{

  int i;
  double rho;
  double params[4];

  for (i = 0; i < 4; i++)
    params[i] = 0.;
  st_compute_params(corpgs, params_in, params);
  rho = corpgs->rho;

  for (i = 0; i < 4; i++)
    M_R(correl,4,i,i) = 1.;

  M_R(correl,4,2,0) = rho;
  M_R(correl,4,3,1) = rho;

  M_R(correl,4,1,0) = params[0];
  M_R(correl,4,2,1) = params[2];
  M_R(correl,4,3,0) = params[1];
  M_R(correl,4,3,2) = params[3];

  matrix_fill_symmetry(4, correl);
}

/****************************************************************************/
/*!
 **  Update the following matrices according to constraints on model
 **
 ** \param[in]  corpgs      Local_CorPgs structure
 ** \param[in,out]  Grad        Vector of gradients (Dimension = npar)
 ** \param[in,out]  Hess        Matrix of Hessian (Dimension = npar * npar)
 ** \param[in,out]  JJ          Matrix of t(JJ) * JJ (Dimension = npar * npar)
 **
 *****************************************************************************/
static void st_update_constraints(Local_CorPgs *corpgs,
                                  double *Grad,
                                  double *Hess,
                                  double *JJ)
{
  double v[4], m[16], *modif;
  int i, j, k, l, npar;

  /* Initializations */

  modif = corpgs->modif;
  npar = corpgs->npar;

  /* Update the Grad */

  matrix_combine(4, 1., Grad, 0., NULL, v);
  matrix_combine(4, 0., NULL, 0., NULL, Grad);
  for (i = 0; i < npar; i++)
    for (j = 0; j < 4; j++)
      Grad[i] += v[j] * M_R(modif, 4, i, j);

  /* Update the Hessian */

  matrix_combine(16, 1., Hess, 0., NULL, m);
  matrix_combine(16, 0., NULL, 0., NULL, Hess);
  for (i = 0; i < npar; i++)
    for (j = 0; j < npar; j++)
      for (k = 0; k < 4; k++)
        for (l = 0; l < 4; l++)
          M_R(Hess,npar,i,j) +=
          M_R(modif,4,i,k) * M_R(m, 4, k, l) * M_R(modif, 4, j, l);

  /* Update JJ */

  if (JJ != NULL)
  {
    matrix_combine(16, 1., JJ, 0., NULL, m);
    matrix_combine(16, 0., NULL, 0., NULL, JJ);
    for (i = 0; i < npar; i++)
      for (j = 0; j < npar; j++)
        for (k = 0; k < 4; k++)
          for (l = 0; l < 4; l++)
            M_R(JJ,npar,i,j) +=
            M_R(modif,4,i,k) * M_R(m, 4, l, k) * M_R(modif, 4, j, l);
  }
}

/****************************************************************************/
/*!
 **  Compute the derivatives (first and second) of the smallest eigenvalue
 **
 ** \param[in] corpgs  Local_CorPgs structure
 ** \param[in] eigval  Current eigen value
 **
 ** \param[out] ev     Output array
 ** \param[out] d1     First order derivative
 ** \param[out] d2     Second order derivative
 **
 *****************************************************************************/
static void st_deriv_eigen(Local_CorPgs *corpgs,
                           double eigval,
                           double *ev,
                           double *d1,
                           double *d2)
{
  double temp[16], invGn[16];
  int i, j;

  matrix_combine(16, 0., NULL, 0., NULL, d2);
  matrix_combine(16, 0., NULL, 0., NULL, temp);
  st_build_correl(corpgs, corpgs->params, temp);
  matrix_combine(16, -1., temp, 0., NULL, temp);

  for (i = 0; i < 4; i++)
    M_R(temp,4,i,i) += eigval;

  invgen(temp, 4, invGn);

  for (i = 0; i < 4; i++)
    d1[i] = 2 * ev[12 + F(i, 0)] * ev[12 + F(i, 1)];

  for (i = 0; i < 4; i++)
    for (j = 0; j < i; j++)
    {
      M_R(d2,4,i,j) = ev[12 + F(i, 0)] * ev[12 + F(j, 0)]
                      * M_R(invGn, 4, F(i,1), F(j,1));
      M_R(d2,4,i,j) += ev[12 + F(i, 1)] * ev[12 + F(j, 0)]
                       * M_R(invGn, 4, F(i,0), F(j,1));
      M_R(d2,4,i,j) += ev[12 + F(i, 0)] * ev[12 + F(j, 1)]
                       * M_R(invGn, 4, F(i,1), F(j,0));
      M_R(d2,4,i,j) += ev[12 + F(i, 1)] * ev[12 + F(j, 1)]
                       * M_R(invGn, 4, F(i,0), F(j,0));
      M_R(d2,4,i,j) *= 2;
    }

  matrix_fill_symmetry(4, d2);
  st_update_constraints(corpgs, d1, d2, NULL);

}

/****************************************************************************/
/*!
 **  Expand the vector of parameters into C1, C12, C21 and C2
 **  according to the constraints
 **
 ** \return  Returned parameter
 **
 ** \param[in]  local_pgs  Local_Pgs structure
 ** \param[in]  igrf       Rank of the first variable
 ** \param[in]  jgrf       Rank of the second variable
 ** \param[in]  idir       positive (1) or negative (-1) distance
 **
 *****************************************************************************/
static double st_param_expand(Local_Pgs *local_pgs,
                              int igrf,
                              int jgrf,
                              int idir)
{
  double rho, rho2;
  Local_CorPgs *corpgs;

  corpgs = &local_pgs->corpgs;
  rho = corpgs->rho;

  switch (corpgs->opt_correl)
  {
    case 0:
      if (igrf == 0 && jgrf == 0)
        return (corpgs->params[0]);
      else if (igrf == 1 && jgrf == 1)
        return (corpgs->params[3]);
      else
      {
        if (idir > 0)
          return (corpgs->params[1]);
        else
          return (corpgs->params[2]);
      }
      break;

    case 1:
      if (igrf == 0 && jgrf == 0)
        return (corpgs->params[0]);
      else if (igrf == 1 && jgrf == 1)
        return (corpgs->params[2]);
      else
        return (corpgs->params[1]);
      break;

    case 2:
      rho2 = rho * rho;
      if (igrf == 0 && jgrf == 0)
        return (corpgs->params[0]);
      else if (igrf == 1 && jgrf == 1)
        return (corpgs->params[0] * rho2 + corpgs->params[1] * (1. - rho2));
      else
        return (corpgs->params[0] * rho);
      break;
  }
  return (0.);
}

/****************************************************************************/
/*!
 **  Compute the modif matrix
 **
 ** \param[out]  corpgs   Local_CorPgs structure
 **
 *****************************************************************************/
static void st_set_modif(Local_CorPgs *corpgs)

{
  double *modif = corpgs->modif;
  double rho, rho2;
  int i;

  rho = corpgs->rho;
  matrix_combine(16, 0., NULL, 0., NULL, modif);

  switch (corpgs->opt_correl)
  {
    case 0: /* Full parameters */
      corpgs->npar = 4;
      for (i = 0; i < 4; i++)
        M_R(modif,4,i,i) = 1.;
      break;

    case 1: /* Symmetrical case */
      corpgs->npar = 3;
      M_R(modif,4,0,0) = 1;
      M_R(modif,4,1,1) = 1;
      M_R(modif,4,1,2) = 1;
      M_R(modif,4,2,3) = 1;
      break;

    case 2: /* Residual case */
      rho2 = rho * rho;
      corpgs->npar = 2;
      M_R(modif,4,0,0) = 1;
      M_R(modif,4,0,1) = rho;
      M_R(modif,4,0,2) = rho;
      M_R(modif,4,0,3) = rho2;
      M_R(modif,4,1,3) = 1 - rho2;
      break;
  }
}

/****************************************************************************/
/*!
 **  Define the correlation option
 **
 ** \param[in]  option      Correlation option
 ** \li                     0 : full parameters (C1, C12, C21, C2)
 ** \li                     1 : symmetrical case (C12 = C21)
 ** \li                     2 : residual case
 ** \param[in]  flag_rho    1 if rho has to be calculated, 0 otherwise
 ** \param[in]  rho         Correlation between GRFs
 ** \param[in,out]  local_pgs   Local_Pgs structure
 **
 *****************************************************************************/
static void st_define_corpgs(int option,
                             int flag_rho,
                             double rho,
                             Local_Pgs *local_pgs)
{
  Local_CorPgs *corpgs;
  int i;

  /* Initializations */

  corpgs = &local_pgs->corpgs;
  corpgs->rho = rho;
  corpgs->opt_correl = option;
  st_set_modif(corpgs);

  for (i = 0; i < 4; i++)
    corpgs->params[i] = 0.;

  corpgs->flag_rho = flag_rho;
}

/****************************************************************************/
/*!
 **  Count the number of pairs with the target facies
 **
 ** \param[in]  local_pgs    Local_Pgs structure
 ** \param[in]  ifac1        First target facies (starting from 1)
 ** \param[in]  ifac2        Second target facies (starting from 1)
 **
 *****************************************************************************/
static int st_get_count(Local_Pgs *local_pgs, int ifac1, int ifac2)
{
  int ipair, i1, i2, number;
  double dist, w1, w2;

  /* Initializations */

  number = 0;

  for (ipair = local_pgs->ifirst; ipair < local_pgs->ilast; ipair++)
  {
    vario_order_get_indices(local_pgs->vorder, ipair, &i1, &i2, &dist);
    if (ifac1 != local_pgs->db->getVariable(i1, 0)) continue;
    if (ifac2 != local_pgs->db->getVariable(i2, 0)) continue;
    w1 = local_pgs->db->getWeight(i1);
    w2 = local_pgs->db->getWeight(i2);
    number += (int) (w1 * w2);
  }
  return (number);
}

/****************************************************************************/
/*!
 ** PRUPOSE: Internal function
 **
 *****************************************************************************/
static double st_rkl(int maxpts,
                     double x,
                     double y,
                     double *lower,
                     double *upper,
                     double *corr1,
                     double *covar,
                     double *temp)
{
  double cste[2], vect[2], mean[2], v1, v2, S, error;
  int inform;
  static double abseps = 1.e-12;
  static double releps = 0.;

  cste[0] = 0.;
  cste[1] = 0.;
  vect[0] = x;
  vect[1] = y;
  matrix_product(2, 2, 1, temp, vect, mean);
  v1 = law_df_bigaussian(vect, cste, corr1);
  mvndst2n(lower, upper, mean, covar, maxpts, abseps, releps, &error, &v2,
           &inform);
  S = v1 * v2;
  return (S);
}

/****************************************************************************/
/*!
 ** \return  Calculate
 ** \return    d/(d xk) * d/(d xl)
 ** \return    int_lower1^upper1
 ** \return    int_lower2^upper2
 ** \return    int_lower3^upper3
 ** \return    int_lower4^upper4
 ** \return    gaussian density(x1,x2,x3,x4,rho) dx1 dx2 dx3 dx4
 **
 ** \param[in]  maxpts       Maximum number of evaluations in mvndst
 ** \param[in]  index1       First derivative index
 ** \param[in]  index2       Second derivative index
 ** \param[in]  lower        Array of lower bounds
 ** \param[in]  upper        Array of upper bounds
 ** \param[in]  correl       Correlation matrix (Dimension = 4*4)
 **
 *****************************************************************************/
static double st_ikl(int maxpts,
                     int index1,
                     int index2,
                     double *lower,
                     double *upper,
                     double *correl)
{
  double low[2], upp[2], value, x, y, S;
  double corr1[4], corr2[4], corrc[4], covar[4], inv_corr1[4], temp[4];
  int i, j, k, index[2];

  // Initializations 
  index[0] = index1;
  index[1] = index2;

  // Build submatrices
  matrix_manage(4, 1, -2, 0, index, NULL, lower, low);
  matrix_manage(4, 1, -2, 0, index, NULL, upper, upp);
  matrix_manage(4, 4, 2, 2, index, index, correl, corr1);
  matrix_manage(4, 4, -2, 2, index, index, correl, corrc);
  matrix_manage(4, 4, -2, -2, index, index, correl, corr2);
  if (matrix_invert_copy(corr1, 2, inv_corr1)) messageAbort("st_ikl #1");
  matrix_product(2, 2, 2, corrc, inv_corr1, temp);

  // Derive covar
  for (i = 0; i < 2; i++)
    for (j = 0; j < 2; j++)
    {
      value = 0;
      for (k = 0; k < 2; k++)
        value += M_R(temp,2,k,i) * M_R(corrc, 2, k, j);
      M_R(covar,2,i,j) = M_R(corr2,2,i,j) - value;
    }

  S = 0.;
  x = upper[index1];
  if (IS_GAUSS_DEF(x))
  {
    y = upper[index2];
    if (IS_GAUSS_DEF(y))
      S += st_rkl(maxpts, x, y, low, upp, corr1, covar, temp);
    y = lower[index2];
    if (IS_GAUSS_DEF(y))
      S -= st_rkl(maxpts, x, y, low, upp, corr1, covar, temp);
  }
  x = lower[index1];
  if (IS_GAUSS_DEF(x))
  {
    y = lower[index2];
    if (IS_GAUSS_DEF(y))
      S += st_rkl(maxpts, x, y, low, upp, corr1, covar, temp);
    y = upper[index2];
    if (IS_GAUSS_DEF(y))
      S -= st_rkl(maxpts, x, y, low, upp, corr1, covar, temp);
  }
  return (S / 2.);
}

/****************************************************************************/
/*!
 ** \return  Calculate the other derivatives
 **
 *****************************************************************************/
static double st_nkl(double *u,
                     double lower,
                     double upper,
                     double *invvari,
                     int index2,
                     double meanj,
                     double varj,
                     double stdj)
{
  double dflow, dfupp, cdflow, cdfupp, invval, invpart[3], total, S;

  dfupp = law_dnorm(upper, meanj, stdj);
  dflow = law_dnorm(lower, meanj, stdj);
  cdfupp = law_cdf_gaussian((upper - meanj) / stdj);
  cdflow = law_cdf_gaussian((lower - meanj) / stdj);
  invval = invvari[index2];
  matrix_manage(4, 1, -1, 0, &index2, NULL, invvari, invpart);
  matrix_product(1, 3, 1, invpart, u, &total);

  S = (dfupp - dflow) * varj * invval - (cdfupp - cdflow)
      * (invval * meanj + total);
  return (S);
}

/****************************************************************************/
/*!
 ** \return  Calculate second-order derivatives
 **
 ** \param[in]  maxpts       Maximum number of evaluations in mvndst
 ** \param[in]  index1       First derivative index
 ** \param[in]  index2       Second derivative index
 ** \param[in]  lower        Array of lower bounds (Dimension = 4)
 ** \param[in]  upper        Array of upper bounds (Dimension = 4)
 ** \param[in]  correl       Correlation matrix (Dimension = 4*4)
 **
 *****************************************************************************/
static double st_d2_dkldkl(int maxpts,
                           int index1,
                           int index2,
                           double *lower,
                           double *upper,
                           double *correl)
{
  double corri[16], v1, v2, S;
  double deltaparam = 1.e-6;
  int i;

  for (i = 0; i < 16; i++)
    corri[i] = correl[i];
  M_R(corri,4,index1,index2) += deltaparam;
  M_R(corri,4,index2,index1) += deltaparam;
  v1 = st_ikl(maxpts, index1, index2, lower, upper, corri);

  for (i = 0; i < 16; i++)
    corri[i] = correl[i];
  M_R(corri,4,index1,index2) -= deltaparam;
  M_R(corri,4,index2,index1) -= deltaparam;
  v2 = st_ikl(maxpts, index1, index2, lower, upper, corri);

  S = (v1 - v2) / (2. * deltaparam);
  return (S / 2.);
}

/****************************************************************************/
/*!
 **  Calculate the gradients for a pair of facies and return it
 **
 ** \return  Calculate d²S/dC12dC21 and d²S/dC1dC2
 **
 ** \param[in]  lower        Array of lower bounds (Dimension = 4)
 ** \param[in]  upper        Array of upper bounds (Dimension = 4)
 ** \param[in]  correl       Correlation matrix (Dimension = 4*4)
 **
 *****************************************************************************/
static double st_d2_dkldij(double *lower, double *upper, double *correl)
{
  int i, i1, i2, i3, i4, grid[4], flag_out;
  double u[4], S;

  S = 0.;
  for (i4 = 0; i4 < 2; i4++)
    for (i3 = 0; i3 < 2; i3++)
      for (i2 = 0; i2 < 2; i2++)
        for (i1 = 0; i1 < 2; i1++)
        {
          grid[0] = i1;
          grid[1] = i2;
          grid[2] = i3;
          grid[3] = i4;

          flag_out = 0;
          for (i = 0; i < 4 && flag_out == 0; i++)
          {
            u[i] = (grid[i]) ? upper[i] :
                               lower[i];
            flag_out = !IS_GAUSS_DEF(u[i]);
          }
          if (!flag_out)
            S += pow(-1., i1 + i2 + i3 + i4) * law_df_quadgaussian(u, correl);
        }
  return (S / 2.);
}

/****************************************************************************/
/*!
 ** \return  Calculate the first-order derivatives
 **
 ** \param[in]  index1      First derivative index
 ** \param[in]  index2      Second derivative index
 ** \param[in]  lower       Array of lower bounds (Dimension = 4)
 ** \param[in]  upper       Array of upper bounds (Dimension = 4)
 ** \param[in]  correl      Correlation matrix (Dimension = 4*4)
 **
 *****************************************************************************/
static double st_d2_dkldkj(int index1,
                           int index2,
                           double *lower,
                           double *upper,
                           double *correl)
{
  double varcori[9], invvarcor[16], invvarcori[4], corr1[9], crosscor[3];
  double invcorr1[9], temp[3], lowi[3], uppi[3], u[3];
  double lowj, uppj, corr2, covar, sdcovar, S, mu, random;
  int i, i1, i2, i3, grid[3], flag_out;

  matrix_manage(4, 4, -1, -1, &index2, &index2, correl, varcori);
  if (matrix_invert_copy(correl, 4, invvarcor)) messageAbort("st_d2_dkldkj #1");
  matrix_manage(4, 4, 1, 0, &index1, NULL, invvarcor, invvarcori);
  matrix_manage(4, 4, -1, -1, &index2, &index2, correl, corr1);
  matrix_manage(4, 4, 1, -1, &index2, &index2, correl, crosscor);
  matrix_manage(4, 4, 1, 1, &index2, &index2, correl, &corr2);
  if (matrix_invert_copy(corr1, 3, invcorr1)) messageAbort("st_d2_dkldkj #2");
  matrix_product(1, 3, 3, crosscor, invcorr1, temp);
  matrix_product(1, 3, 1, crosscor, temp, &covar);
  covar = corr2 - covar;
  sdcovar = sqrt(covar);
  matrix_manage(4, 1, -1, 0, &index2, NULL, lower, lowi);
  matrix_manage(4, 1, -1, 0, &index2, NULL, upper, uppi);
  matrix_manage(4, 1, 1, 0, &index2, NULL, lower, &lowj);
  matrix_manage(4, 1, 1, 0, &index2, NULL, upper, &uppj);

  S = 0.;
  for (i3 = 0; i3 < 2; i3++)
    for (i2 = 0; i2 < 2; i2++)
      for (i1 = 0; i1 < 2; i1++)
      {
        grid[0] = i1;
        grid[1] = i2;
        grid[2] = i3;

        for (i = flag_out = 0; i < 3 && flag_out == 0; i++)
        {
          u[i] = (grid[i]) ? uppi[i] :
                             lowi[i];
          flag_out = !IS_GAUSS_DEF(u[i]);
        }
        if (flag_out) continue;

        matrix_product(1, 3, 1, temp, u, &mu);
        random = law_df_multigaussian(3, u, varcori);

        S += pow(-1., 3 - i1 + i2 + i3) * random
             * st_nkl(u, lowj, uppj, invvarcori, index2, mu, covar, sdcovar);
      }
  return (S / 2);
}

/****************************************************************************/
/*!
 **  Global calculation in the stationary case
 **
 ** \return  The global score
 **
 ** \param[in]  local_pgs    Local_Pgs structure
 ** \param[in]  flag_deriv   1 if the derivatives must be calculated
 ** \param[in]  flag_reset   1 to update the probability calculations
 ** \param[in]  correl       Correlation matrix updated
 **
 ** \param[out] Grad         Vector of cumulated gradients (Dimension= 4)
 ** \param[out] Hess         Matrix of cumulated Hessian (Dimension= 4*4)
 ** \param[out] JJ           Matrix of cumulated JJ (Dimension= 4*4)
 **
 *****************************************************************************/
static double st_calcul_stat(Local_Pgs *local_pgs,
                             int flag_deriv,
                             int flag_reset,
                             double *correl,
                             double *Grad,
                             double *Hess,
                             double *JJ)
{
  double grad[4], lower[4], upper[4], hess[16], gradgrad[16];
  double s, S, rj2, erval, ggval;
  int ifac1, ifac2, nfifj, i, j, inform;
  static double abseps = 1.e-6;
  static double releps = 0.;
  static int maxpts4 = 10000;
  static int maxpts2 = 4000;

  S = 0.;
  matrix_combine(4, 0., NULL, 0., NULL, grad);
  matrix_combine(16, 0., NULL, 0., NULL, hess);

  for (ifac1 = 0; ifac1 < local_pgs->nfacies; ifac1++)
  {
    for (ifac2 = 0; ifac2 < local_pgs->nfacies; ifac2++)
    {
      nfifj = st_get_count(local_pgs, ifac1 + 1, ifac2 + 1);
      if (nfifj <= 0) continue;

      /* Get the bounds */

      VectorDouble bounds;
      bounds = local_pgs->rule->getThresh(ifac1 + 1);
      lower[0] = bounds[0];
      upper[0] = bounds[1];
      lower[2] = bounds[2];
      upper[2] = bounds[3];
      bounds = local_pgs->rule->getThresh(ifac2 + 1);
      lower[1] = bounds[0];
      upper[1] = bounds[1];
      lower[2] = bounds[2];
      lower[3] = bounds[3];
      if (flag_reset)
      {
        mvndst4(lower, upper, correl, maxpts4, abseps, releps, &erval, &s,
                &inform);
        STAT_PROBA(ifac1,ifac2) = s;
      }
      else
        s = STAT_PROBA(ifac1, ifac2);

      rj2 = -2. * log(s);
      S += (double) nfifj * rj2;

      /* Calculate the derivative */

      if (!flag_deriv) continue;
      grad[0] = -st_ikl(maxpts2, 0, 1, lower, upper, correl) / s;
      grad[1] = -st_ikl(maxpts2, 0, 3, lower, upper, correl) / s;
      grad[2] = -st_ikl(maxpts2, 1, 2, lower, upper, correl) / s;
      grad[3] = -st_ikl(maxpts2, 2, 3, lower, upper, correl) / s;
      matrix_product(4, 1, 4, grad, grad, gradgrad);

      M_R(hess,4,3,0) = M_R(hess,4,2,1) = -st_d2_dkldij(lower, upper, correl);
      M_R(hess,4,1,0) = st_d2_dkldkj(0, 2, lower, upper, correl);
      M_R(hess,4,2,0) = st_d2_dkldkj(1, 3, lower, upper, correl);
      M_R(hess,4,3,1) = st_d2_dkldkj(3, 1, lower, upper, correl);
      M_R(hess,4,3,2) = st_d2_dkldkj(2, 0, lower, upper, correl);
      M_R(hess,4,0,0) = st_d2_dkldkl(maxpts4, 0, 1, lower, upper, correl);
      M_R(hess,4,1,1) = st_d2_dkldkl(maxpts4, 0, 3, lower, upper, correl);
      M_R(hess,4,2,2) = st_d2_dkldkl(maxpts4, 1, 2, lower, upper, correl);
      M_R(hess,4,3,3) = st_d2_dkldkl(maxpts4, 2, 3, lower, upper, correl);
      matrix_fill_symmetry(4, hess);

      for (i = 0; i < 4; i++)
      {
        Grad[i] += nfifj * grad[i];
        for (j = 0; j < 4; j++)
        {
          ggval = M_R(gradgrad, 4, i, j);
          M_R(Hess,4,i,j) -= nfifj * (M_R(hess,4,i,j) / s - ggval);
          M_R(JJ,4,i,j) += nfifj * ggval / rj2;
        }
      }
    }
  }
  return (S / 2.);
}

/****************************************************************************/
/*!
 **  Global calculation in the non-stationary case
 **
 ** \return  The global score
 **
 ** \param[in]  local_pgs    Local_Pgs structure
 ** \param[in]  flag_deriv   1 if the derivatives must be calculated
 ** \param[in]  flag_reset   1 to update the probability calculations
 ** \param[in]  correl       Correlation matrix updated
 **
 ** \param[out] Grad         Vector of cumulated gradients (Dimension= 4)
 ** \param[out] Hess         Matrix of cumulated Hessian (Dimension= 4*4)
 ** \param[out] JJ           Matrix of cumulated JJ (Dimension= 4*4)
 **
 *****************************************************************************/
static double st_calcul_nostat(Local_Pgs *local_pgs,
                               int flag_deriv,
                               int flag_reset,
                               double *correl,
                               double *Grad,
                               double *Hess,
                               double *JJ)
{
  double grad[4], lower[4], upper[4], hess[16], gradgrad[16];
  double s, S, rj2, erval, ggval, dist, w1, w2;
  int i1, i2, ifac1, ifac2, i, j, inform, ipair;
  static double abseps = 1.e-6;
  static double releps = 0.;
  static int maxpts4 = 10000;
  static int maxpts2 = 4000;

  S = 0.;
  matrix_combine(4, 0., NULL, 0., NULL, grad);
  matrix_combine(16, 0., NULL, 0., NULL, hess);

  for (ipair = local_pgs->ifirst; ipair < local_pgs->ilast; ipair++)
  {
    vario_order_get_indices(local_pgs->vorder, ipair, &i1, &i2, &dist);
    ifac1 = (int) local_pgs->db->getVariable(i1, 0);
    ifac2 = (int) local_pgs->db->getVariable(i2, 0);
    w1 = local_pgs->db->getWeight(i1);
    w2 = local_pgs->db->getWeight(i2);
    /* Get the bounds */

    (void) rule_thresh_define(local_pgs->propdef, local_pgs->db,
                              local_pgs->rule, ifac1, i1, 0, 0, 1, &lower[0],
                              &upper[0], &lower[2], &upper[2]);
    (void) rule_thresh_define(local_pgs->propdef, local_pgs->db,
                              local_pgs->rule, ifac2, i2, 0, 0, 1, &lower[1],
                              &upper[1], &lower[3], &upper[3]);

    if (flag_reset)
    {
      mvndst4(lower, upper, correl, maxpts4, abseps, releps, &erval, &s,
              &inform);
      MEMINT(ipair) = s;
    }
    else
      s = MEMINT(ipair);
    rj2 = -2. * log(s);
    S += w1 * w2 * rj2;

    /* Calculate the derivative */

    if (!flag_deriv) continue;
    grad[0] = -st_ikl(maxpts2, 0, 1, lower, upper, correl) / s;
    grad[1] = -st_ikl(maxpts2, 0, 3, lower, upper, correl) / s;
    grad[2] = -st_ikl(maxpts2, 1, 2, lower, upper, correl) / s;
    grad[3] = -st_ikl(maxpts2, 2, 3, lower, upper, correl) / s;
    matrix_product(4, 1, 4, grad, grad, gradgrad);

    M_R(hess,4,3,0) = M_R(hess,4,2,1) = -st_d2_dkldij(lower, upper, correl);
    M_R(hess,4,1,0) = st_d2_dkldkj(0, 2, lower, upper, correl);
    M_R(hess,4,2,0) = st_d2_dkldkj(1, 3, lower, upper, correl);
    M_R(hess,4,3,1) = st_d2_dkldkj(3, 1, lower, upper, correl);
    M_R(hess,4,3,2) = st_d2_dkldkj(2, 0, lower, upper, correl);
    M_R(hess,4,0,0) = st_d2_dkldkl(maxpts4, 0, 1, lower, upper, correl);
    M_R(hess,4,1,1) = st_d2_dkldkl(maxpts4, 0, 3, lower, upper, correl);
    M_R(hess,4,2,2) = st_d2_dkldkl(maxpts4, 1, 2, lower, upper, correl);
    M_R(hess,4,3,3) = st_d2_dkldkl(maxpts4, 2, 3, lower, upper, correl);
    matrix_fill_symmetry(4, hess);

    for (i = 0; i < 4; i++)
    {
      Grad[i] += w1 * w2 * grad[i];
      for (j = 0; j < 4; j++)
      {
        ggval = M_R(gradgrad, 4, i, j);
        M_R(Hess,4,i,j) -= w1 * w2 * (M_R(hess,4,i,j) / s - ggval);
        M_R(JJ,4,i,j) += w1 * w2 * ggval / rj2;
      }
    }
  }
  return (S / 2.);
}

/****************************************************************************/
/*!
 **  Global calculation
 **
 ** \return  The global score
 **
 ** \param[in]  local_pgs    Local_Pgs structure
 ** \param[in]  flag_deriv   1 if the derivatives must be calculated
 ** \param[in]  flag_reset   1 to update the probability calculations
 ** \param[in]  params       Array of parameters
 **
 ** \param[out] Grad         Vector of cumulated gradients (Dimension= 4)
 ** \param[out] Hess         Matrix of cumulated Hessian (Dimension= 4*4)
 ** \param[out] JJ           Matrix of cumulated JJ (Dimension= 4*4)
 **
 *****************************************************************************/
static double st_calcul(Local_Pgs *local_pgs,
                        int flag_deriv,
                        int flag_reset,
                        double *params,
                        double *Grad,
                        double *Hess,
                        double *JJ)
{
  double S, correl[16];

  S = 0.;

  if (flag_deriv)
  {
    matrix_combine(4, 0., NULL, 0., NULL, Grad);
    matrix_combine(16, 0., NULL, 0., NULL, Hess);
    matrix_combine(16, 0., NULL, 0., NULL, JJ);
  }

  st_build_correl(&local_pgs->corpgs, params, correl);

  if (local_pgs->flag_stat)
    S = st_calcul_stat(local_pgs, flag_deriv, flag_reset, correl, Grad, Hess,
                       JJ);
  else
    S = st_calcul_nostat(local_pgs, flag_deriv, flag_reset, correl, Grad, Hess,
                         JJ);

  /* Modify the results due to the constraints on parameters */

  if (flag_deriv) st_update_constraints(&local_pgs->corpgs, Grad, Hess, JJ);

  return (S / 2.);
}

/****************************************************************************/
/*!
 **  Initialize the parameters
 **
 ** \param[out]  corpgs  Local_CorPgs Structure
 **
 *****************************************************************************/
static void st_initialize_params(Local_CorPgs *corpgs)
{
  double r = corpgs->rho;

  switch (corpgs->opt_correl)
  {
    case 0:
      corpgs->params[0] = corpgs->params[3] = fabs(r);
      corpgs->params[1] = corpgs->params[2] = fabs(r) * r;
      break;

    case 1:
      corpgs->params[0] = corpgs->params[2] = fabs(r);
      corpgs->params[1] = fabs(r) * r;
      break;

    case 2:
      corpgs->params[0] = corpgs->params[1] = fabs(r);
      break;
  }

}

/****************************************************************************/
/*!
 **  Optimize the lag
 **
 ** \param[in]  local_pgs    Local_Pgs structure
 ** \param[in]  tolsort      Tolerance value
 ** \param[in]  new_val      Flag indicating if parameters must be initialized
 **
 *****************************************************************************/
static double st_optim_onelag_pgs(Local_Pgs *local_pgs,
                                  double tolsort,
                                  int new_val)
{
  Local_CorPgs *corpgs;
  double Hess[16], JJ[16], eigval[4], eigvec[16], Gn[16], invGn[16], correl[16],
      d2[16];
  double hsd[4], hgn[4], step[4], a[4], gr[4], param_temp[4], hgna[4], Grad[4],
      d1[4];
  double normgrad, normgrad2, alpha, delta2, c, a2, beta, hgna2, niter;
  double Sr, Snew, delta, mdiminution, mdiminution_pred, stepgr, rval;
  int flag_sortie, flag_moved, npar, npar2, i;

  static double maxiter = 100;
  static double delta0 = 1;

  double Spen = 0.;
  double Srpen = 0.;
  double penalize = 1000;
  int barrier = 0;

  /* Initializations */

  corpgs = &local_pgs->corpgs;
  npar = corpgs->npar;
  npar2 = npar * npar;
  delta = delta0;
  mdiminution = 0.;

  if (new_val) st_initialize_params(corpgs);

  matrix_combine(npar, 1., corpgs->params, 0., NULL, param_temp);
  matrix_combine(npar, 0., NULL, 0., NULL, Grad);
  matrix_combine(npar2, 0., NULL, 0., NULL, Hess);
  matrix_combine(npar2, 0., NULL, 0., NULL, JJ);

  /* Calculate the score and the derivatives */

  Sr = st_calcul(local_pgs, 1, 1, corpgs->params, Grad, Hess, JJ);

  niter = 0.;
  flag_sortie = 0;
  flag_moved = 1;

  while (!flag_sortie)
  {
    if (barrier)
    {
      st_build_correl(corpgs, param_temp, correl);
      matrix_eigen(correl, 4, eigval, eigvec);
      st_deriv_eigen(corpgs, eigval[3], eigvec, d1, d2);
      Srpen = Sr - penalize * log(eigval[3]);
      matrix_combine(npar, 1, Grad, -penalize / eigval[3], d1, Grad);
      matrix_combine(npar, npar, Hess, penalize / (eigval[3] * eigval[3]), d2,
                     Hess);
      matrix_combine(npar, npar, JJ, penalize / (eigval[3] * eigval[3]), d2,
                     JJ);
      penalize *= 0.5;
    }
    niter++;
    delta2 = delta * delta;
    if (flag_moved)
    {
      matrix_combine(npar, 1., Grad, 0., NULL, gr);
      matrix_combine(npar2, 1., Hess, 0., NULL, Gn);
      if (!is_matrix_definite_positive(npar, Gn, eigval, eigvec, 0))
        matrix_combine(npar2, 1., JJ, 0., NULL, Gn);
      matrix_combine(npar, -1., gr, 0., NULL, hsd);
      if (invgen(Gn, npar, invGn)) messageAbort("st_optim_lag");
      matrix_product(npar, npar, 1, invGn, hsd, hgn);
    }

    /* Determine the lag (hgn, alpha*hsd) or a convex combinaison of both */

    if (matrix_norm(hgn, npar) <= delta2)
    {
      matrix_combine(npar, 1., hgn, 0., NULL, step);
    }
    else
    {
      normgrad2 = matrix_norm(gr, npar);
      alpha = normgrad2 / matrix_normA(gr, Gn, npar, npar);
      normgrad = sqrt(normgrad2);
      if (normgrad > (delta / alpha))
      {
        matrix_combine(npar, delta / normgrad, hsd, 0., NULL, step);
      }
      else
      {
        matrix_combine(npar, alpha, hsd, 0., NULL, a);
        matrix_combine(npar, 1., hgn, -1., a, hgna);
        matrix_product(1, npar, 1, a, hgn, &c);
        a2 = matrix_norm(a, npar);
        hgna2 = matrix_norm(hgna, npar);
        if (c <= 0.)
          beta = (-c + sqrt(c * c + hgna2 * (delta2 - a2))) / hgna2;
        else
          beta = (delta2 - a2) / (c + sqrt(c * c + hgna2 * (delta2 - a2)));
        matrix_combine(npar, beta, hgn, (1. - beta), a, step);
      }
    }

    matrix_combine(npar, 1., step, 1., corpgs->params, param_temp);
    st_build_correl(corpgs, param_temp, correl);
    while (!is_matrix_definite_positive(4, correl, eigval, eigvec, 0))
    {
      matrix_combine(npar, 0.9, step, 0., NULL, step);
      matrix_combine(npar, 1., step, 1., corpgs->params, param_temp);
      st_build_correl(corpgs, param_temp, correl);
    }

    Snew = st_calcul(local_pgs, 0, 1, param_temp, Grad, Hess, JJ);

    if (barrier) Spen = Snew - penalize * log(eigval[3]);
    if (!FFFF(Snew))
    {

      mdiminution = Snew - Sr;
      if (barrier) mdiminution = Spen - Srpen;
      matrix_product(1, npar, 1, step, gr, &stepgr);
      mdiminution_pred = stepgr + 0.5 * matrix_normA(step, Gn, npar, npar);
      rval = mdiminution / mdiminution_pred;
      flag_moved = (mdiminution < 0);
    }
    else
    {
      flag_moved = 0;
      rval = 0.;
    }

    if (flag_moved)
    {
      Sr = Snew;
      Srpen = Spen;
      matrix_combine(npar, 1, param_temp, 0, NULL, corpgs->params);
      Snew = st_calcul(local_pgs, 1, 0, corpgs->params, Grad, Hess, JJ);
      if (barrier)
      {
        st_deriv_eigen(corpgs, eigval[3], eigvec, d1, d2);
        matrix_combine(npar, 1, Grad, penalize / eigval[3], d1, Grad);
        matrix_combine(npar, npar, Hess, -penalize / (eigval[3] * eigval[3]),
                       d2, Hess);
        matrix_combine(npar, npar, JJ, -penalize / (eigval[3] * eigval[3]), d2,
                       JJ);
        penalize /= 2.;
      }
      if (rval > 0.75) delta = MAX(delta, 3. * sqrt(matrix_norm(step, npar)));
    }
    if (rval < 0.25) delta /= 2.;

    flag_sortie = (matrix_norminf(npar, step) < tolsort || niter == maxiter
                   || matrix_norminf(npar, Grad) < 0.05
                   || (fabs(mdiminution) < tolsort && flag_moved));
  }

  /* Returning arguments */

  if (debug_query("converge"))
  {
    message("Lag %d - S = %lf Parameters =", local_pgs->ipascur, Sr);
    for (i = 0; i < corpgs->npar; i++)
      message(" %lf", corpgs->params[i]);
    message("\n");
  }

  /* Store the trace */

  trace_add_row(local_pgs);
  trace_define(local_pgs, niter, Snew, 0, npar, Grad);

  return (Sr);
}

/****************************************************************************/
/*!
 **  Discard a data if:
 **  - its facies is unknown
 **  - its thresholds are so close that it leads to a zero probability
 **
 ** \return  Error return code
 **
 ** \param[in]  local_pgs   Local_Pgs structure
 ** \param[in]  iech        Rank of the sample
 **
 *****************************************************************************/
static int st_discard_point(Local_Pgs *local_pgs, int iech)
{
  int ifac;
  double low, up;

  /* This function is bypassed in the stationary case */

  if (local_pgs->flag_stat) return (0);

  /* The following checks must not be performed if not on facies */

  if (!local_pgs->flag_facies) return (0);

  /* Check on the facies */

  ifac = (int) local_pgs->db->getVariable(iech, 0);
  if (ifac < 1 || ifac > local_pgs->nfacies) return (1);

  /* Check on the thresholds */

  if (!TEST_DISCRET)
  {
    if (local_pgs->db->getIntervalNumber() <= 0) return (0);
    low = local_pgs->db->getLowerBound(iech, local_pgs->igrfcur);
    up = local_pgs->db->getUpperBound(iech, local_pgs->igrfcur);
  }
  else
  {
    if (get_LOCATOR_NITEM(local_pgs->db, ELoc::RKLOW) <= 0 && get_LOCATOR_NITEM(
        local_pgs->db, ELoc::RKUP)
                                                              <= 0) return (0);
    low = local_pgs->db->getLowerInterval(iech, local_pgs->igrfcur);
    up = local_pgs->db->getUpperInterval(iech, local_pgs->igrfcur);
  }
  if (up <= low) return (1);

  return (0);
}

/****************************************************************************/
/*!
 **  Compress the Geometry of all pairs
 **
 ** \return  Error return code
 **
 ** \param[in]  local_pgs   Local_Pgs structure
 **
 *****************************************************************************/
static int st_variogram_geometry_pgs_final(Local_Pgs *local_pgs)
{
  int npair;

  local_pgs->vorder = vario_order_final(local_pgs->vorder, &npair);
  if (local_pgs->vorder == (Vario_Order*) NULL) return (1);
  if (npair > 0 && !local_pgs->flag_stat)
  {
    local_pgs->memint.resize(npair);
  }

  return (0);
}

/****************************************************************************/
/*!
 **  Correct the experimental variogram for GRFs
 **
 ** \param[in]  local_pgs   Local_Pgs structure
 ** \param[in]  vario       Vario structure
 ** \param[in]  idir        Rank of the direction
 **
 *****************************************************************************/
static void st_variogram_geometry_pgs_correct(Local_Pgs *local_pgs,
                                              Vario *vario,
                                              int idir)
{
  int ipas, igrf, jgrf, iad, ngrf;

  ngrf = local_pgs->ngrf;

  for (ipas = 0; ipas < vario->getLagNumber(idir); ipas++)
    for (igrf = 0; igrf < ngrf; igrf++)
      for (jgrf = 0; jgrf <= igrf; jgrf++)
      {
        iad = vario->getDirAddress(idir, igrf, jgrf, ipas, false, 1);
        vario->setGgByIndex(idir, iad, st_param_expand(local_pgs, igrf, jgrf, 1));
        if (vario->getSwByIndex(idir, iad) > 0.)
          vario->setHhByIndex(idir, iad,
                       vario->getHhByIndex(idir, iad) / vario->getSwByIndex(idir, iad));
        iad = vario->getDirAddress(idir, igrf, jgrf, ipas, false, -1);
        vario->setGgByIndex(idir, iad, st_param_expand(local_pgs, igrf, jgrf, -1));
        if (vario->getSwByIndex(idir, iad) > 0.)
          vario->setHhByIndex(idir, iad,
                       vario->getHhByIndex(idir, iad) / vario->getSwByIndex(idir, iad));
      }
}

/****************************************************************************/
/*!
 **  Determine the Geometry of all pairs
 **
 ** \return  Error return code
 **
 ** \param[in]  local_pgs   Local_Pgs structure
 ** \param[in]  vario       Vario structure
 ** \param[in]  idir        Rank of the direction
 **
 *****************************************************************************/
static int st_variogram_geometry_pgs_calcul(Local_Pgs *local_pgs,
                                            Vario *vario,
                                            int idir)
{
  int iech, jech, iiech, jjech, nech, ipas, iad, ivar, jvar, nvar, error;
  Db *db;
  double psmin, ps, dist, maxdist;
  VectorInt rindex;

  /* Retrieve information from Local_pgs structure */

  error = 1;
  db = local_pgs->db;
  nech = db->getSampleNumber();
  nvar = vario->getVariableNumber();
  maxdist = vario->getMaximumDistance(idir);
  const DirParam &dirparam = vario->getDirParam(idir);

  /* Initializations */

  ps = 0.;
  psmin = _variogram_convert_angular_tolerance(dirparam.getTolAngle());

  /* Sort the data */

  rindex = db->getSortArray();

  /* Loop on the first point */

  for (iiech = 0; iiech < nech - 1; iiech++)
  {
    iech = rindex[iiech];
    if (!db->isActive(iech)) continue;
    if (FFFF(db->getWeight(iech))) continue;
    if (st_discard_point(local_pgs, iech)) continue;
    mes_process("Calculating Variogram Geometry", nech, iech);

    for (jjech = iiech + 1; jjech < nech; jjech++)
    {
      jech = rindex[jjech];
      if (variogram_maximum_dist1D_reached(db, iech, jech, maxdist)) break;
      if (!db->isActive(jech)) continue;
      if (FFFF(db->getWeight(jech))) continue;
      if (st_discard_point(local_pgs, jech)) continue;

      /* Check if the pair must be kept (Code criterion) */

      if (code_comparable(db, db, iech, jech, dirparam.getOptionCode(),
                          (int) dirparam.getTolCode())) continue;

      /* Check if the pair must be kept */

      dist = distance_intra(db, iech, jech, NULL);
      if (variogram_reject_pair(db, iech, jech, dist, psmin,
                                dirparam.getBench(), dirparam.getCylRad(),
                                dirparam.getCodir(), &ps)) continue;

      /* Get the rank of the lag */

      ipas = variogram_get_lag(vario, idir, ps, psmin, &dist);
      if (IFFFF(ipas)) continue;

      /* Add the sample (only positive lags are of interest) */

      if (ipas < 0) ipas = -ipas;
      if (vario_order_add(local_pgs->vorder, iech, jech, NULL, NULL, ipas, idir,
                          dist)) goto label_end;
      dist = ABS(dist);

      /* Update the distance and weight for all GRFs */

      for (ivar = 0; ivar < nvar; ivar++)
        for (jvar = 0; jvar <= ivar; jvar++)
        {
          if (vario->getFlagAsym())
          {
            iad = vario->getDirAddress(idir, ivar, jvar, ipas, false, 1);
            vario->setGgByIndex(idir, iad, 0.);
            vario->setHhByIndex(idir, iad, vario->getHhByIndex(idir, iad) - dist);
            vario->setSwByIndex(idir, iad, vario->getSwByIndex(idir, iad) + 1);
            iad = vario->getDirAddress(idir, ivar, jvar, ipas, false, -1);
            vario->setGgByIndex(idir, iad, 0.);
            vario->setHhByIndex(idir, iad, vario->getHhByIndex(idir, iad) + dist);
            vario->setSwByIndex(idir, iad, vario->getSwByIndex(idir, iad) + 1);
          }
          else
          {
            iad = vario->getDirAddress(idir, ivar, jvar, ipas, false, 0);
            vario->setGgByIndex(idir, iad, 0.);
            vario->setHhByIndex(idir, iad, vario->getHhByIndex(idir, iad) + dist);
            vario->setSwByIndex(idir, iad, vario->getSwByIndex(idir, iad) + 1);
          }
        }
    }
  }

  /* Set the error return code */

  error = 0;

  label_end: return (error);
}

/****************************************************************************/
/*!
 **  Set the model-type (opt_correl)
 **
 ** \param[in]  opt     The model-type to set
 **
 ** \param[out] corpgs  Local_CorPgs structure
 **
 ****************************************************************************/
static void st_set_opt_correl(int opt, Local_CorPgs *corpgs)

{
  double params[4];
  int i = 0;

  for (i = 0; i < 4; i++)
    params[i] = 0;

  st_compute_params(corpgs, corpgs->params, params);

  switch (opt)
  {
    case 0:
      corpgs->params[0] = params[0];
      corpgs->params[1] = params[1];
      corpgs->params[2] = params[2];
      corpgs->params[3] = params[3];
      break;

    case 1:
      corpgs->params[0] = params[0];
      corpgs->params[1] = corpgs->params[2] = (params[1] + params[2]) / 2.;
      corpgs->params[2] = params[3];
      break;

    case 2:
      corpgs->params[0] = params[0];
      corpgs->params[1] = params[3];
      break;

  }

  corpgs->opt_correl = opt;
  st_set_modif(corpgs);

}

/****************************************************************************/
/*!
 **  Evaluate the variogram of the underlying GRFs (assuming the two GRFs
 **  of the PGS model are correlated)
 **
 ** \param[in]  local_pgs  Local_Pgs structure
 ** \param[in]  idir      Rank of the direction
 **
 *****************************************************************************/
static double st_varcalc_correlated_grf(Local_Pgs *local_pgs, int idir)
{
  double value;
  int ipas, iad, igrf, jgrf, opt_temp;
  Vario *vario;

  opt_temp = local_pgs->corpgs.opt_correl;
  value = 0.;
  vario = local_pgs->vario;

  for (ipas = 0; ipas < vario->getLagNumber(idir); ipas++)
  {
    mes_process("Inverting Variogram Lag", vario->getLagNumber(idir), ipas);
    local_pgs->ipascur = ipas;
    trace_add_row(local_pgs);
    if (!LAG_USED(idir, ipas)) continue;
    vario_order_get_bounds(local_pgs->vorder, idir, ipas, &local_pgs->ifirst,
                           &local_pgs->ilast);
    if (local_pgs->ifirst >= local_pgs->ilast) continue;

    if (opt_temp != 2) st_set_opt_correl(2, &local_pgs->corpgs);

    st_optim_onelag_pgs(local_pgs, 1e-3, 1);
    st_set_opt_correl(opt_temp, &local_pgs->corpgs);
    value += (vario->getUtilizeByIndex(idir, vario->getLagNumber(idir) + ipas)
        * st_optim_onelag_pgs(local_pgs, 1e-3, 0));

    for (igrf = 0; igrf < local_pgs->ngrf; igrf++)
      for (jgrf = 0; jgrf <= igrf; jgrf++)
      {
        iad = vario->getDirAddress(idir, igrf, jgrf, ipas, false, 1);
        vario->setGgByIndex(idir, iad, st_param_expand(local_pgs, igrf, jgrf, 1));
        iad = vario->getDirAddress(idir, igrf, jgrf, ipas, false, -1);
        vario->setGgByIndex(idir, iad, st_param_expand(local_pgs, igrf, jgrf, -1));
      }
  }
  return (value);
}

/****************************************************************************/
/*!
 **  Manage the Local_CorPgs structure
 **
 ** \param[in,out]  local_corpgs  Local_CorPgs structure
 **
 *****************************************************************************/
static void st_manage_corpgs(Local_CorPgs *local_corpgs)

{
  int i;

  local_corpgs->opt_correl = 0;
  local_corpgs->npar = 0;
  local_corpgs->flag_rho = 0;
  local_corpgs->rho = 0.;

  for (i = 0; i < 4; i++)
    local_corpgs->params[i] = 0.;
  for (i = 0; i < 16; i++)
    local_corpgs->modif[i] = 0.;
}

/****************************************************************************/
/*!
 **  Manage the Local_TracePgs structure
 **
 ** \param[in,out]  local_tracepgs  Local_TracePgs structure
 **
 *****************************************************************************/
static void st_manage_trace(Local_TracePgs *local_tracepgs)
{
  local_tracepgs->flag_trace = 0;
  local_tracepgs->idir = 0;
  local_tracepgs->ipas = 0;
  local_tracepgs->nrow = 0;
  local_tracepgs->ncol = 0;
  local_tracepgs->trace = VectorDouble();
}

/****************************************************************************/
/*!
 **  Manage the Local_Pgs structure
 **
 ** \param[in]  mode         0 initialization; 1 allocation; -1 deallocation
 ** \param[in,out] local_pgs Local_Pgs structure
 ** \param[in]  db           Db structure
 ** \param[in]  rule         Lithotype Rule definition
 ** \param[in]  vario        Vario structure
 ** \param[in]  varioind     Indicator Vario structure
 ** \param[in]  model        Model structure
 ** \param[in]  propdef      PropDef structure
 ** \param[in]  flag_stat    1 for stationary; 0 otherwise
 ** \param[in]  flag_facies  1 when processed on facies; 0 otherwise
 ** \param[in]  flag_dist    1 if distances are stored; 0 otherwise
 ** \param[in]  ngrf         Number of GRFs
 ** \param[in]  nfacies      Number of facies
 ** \param[in]  calcul_type  Type of the calculation (covariance, variogram, ...)
 **
 *****************************************************************************/
static void st_manage_pgs(int mode,
                          Local_Pgs *local_pgs,
                          Db *db = nullptr,
                          const Rule *rule = nullptr,
                          Vario *vario = nullptr,
                          Vario *varioind = nullptr,
                          Model *model = nullptr,
                          PropDef *propdef = nullptr,
                          int flag_stat = 0,
                          int flag_facies = 0,
                          int flag_dist = 0,
                          int ngrf = 0,
                          int nfacies = 0,
                          const ECalcVario &calcul_type = ECalcVario::UNDEFINED)
{
  /* Dispatch */

  switch (mode)
  {
    case 0:
      local_pgs->db = nullptr;
      local_pgs->rule = nullptr;
      local_pgs->propdef = nullptr;
      local_pgs->flag_stat = 0;
      local_pgs->flag_facies = 0;
      local_pgs->calcul_type = ECalcVario::UNDEFINED;
      local_pgs->igrfcur = 0;
      local_pgs->idircur = 0;
      local_pgs->ipascur = 0;
      local_pgs->ngrf = 0;
      local_pgs->npair = 0;
      local_pgs->nfacies = 0;
      local_pgs->ifirst = 0;
      local_pgs->ilast = 0;
      local_pgs->d0 = VectorDouble();
      local_pgs->d1 = VectorDouble();
      local_pgs->memint = VectorDouble();
      local_pgs->stat_proba = VectorDouble();
      local_pgs->stat_thresh = VectorDouble();
      local_pgs->model = nullptr;
      local_pgs->vario = nullptr;
      local_pgs->varioind = nullptr;
      local_pgs->vorder = (Vario_Order*) NULL;
      break;

    case 1:
      local_pgs->db = db;
      local_pgs->rule = rule;
      local_pgs->propdef = propdef;
      local_pgs->flag_stat = flag_stat;
      local_pgs->flag_facies = flag_facies;
      local_pgs->calcul_type = calcul_type;
      local_pgs->igrfcur = 0;
      local_pgs->ipascur = 0;
      local_pgs->ngrf = ngrf;
      local_pgs->npair = 0;
      local_pgs->nfacies = nfacies;
      local_pgs->vario = vario;
      local_pgs->varioind = varioind;
      local_pgs->model = model;
      if (model != nullptr)
      {
        int ndim = model->getDimensionNumber();
        local_pgs->d0.resize(ndim);
        local_pgs->d1.resize(ndim);
      }
      local_pgs->vorder = vario_order_manage(1, flag_dist, 0, NULL);
      if (flag_stat)
      {
        local_pgs->stat_proba.resize(nfacies * nfacies, 0.);
        local_pgs->stat_thresh.resize(nfacies * 2 * 2, 0.); // Do not use ngrf, use 2 instead
      }
      st_manage_corpgs(&local_pgs->corpgs);
      st_manage_trace(&local_pgs->tracepgs);
      break;

    case -1:
      local_pgs->vorder = vario_order_manage(-1, 0, 0, local_pgs->vorder);
      break;
  }

  return;
}

/****************************************************************************/
/*!
 **  Performing the variogram calculations
 **
 ** \return  Error return code
 **
 ** \param[in]  vario         Vario structure for the GRFs to be filled
 ** \param[in]  rule          Lithotype Rule definition
 ** \param[in]  local_pgs     Local_Pgs structure
 ** \param[in]  ngrf          Number of GRFs
 ** \param[in]  opt_correl    0 full model; 1 symetrical; 2 residuals
 ** \param[out] flag_geometry 1 if Geometry must be established per direction
 **                           0 if Geometry is already calculated before
 **                             calling this function
 **
 *****************************************************************************/
static int st_variopgs_calcul_norho(Vario *vario,
                                    const Rule *rule,
                                    Local_Pgs *local_pgs,
                                    int ngrf,
                                    int opt_correl,
                                    int flag_geometry)
{

  int idir;

  st_set_rho(rule->getRho(), local_pgs);

  /* Loop on the directions */

  for (idir = 0; idir < vario->getDirectionNumber(); idir++)
  {
    local_pgs->idircur = idir;

    /* Establish the geometry */

    if (flag_geometry)
    {
      if (st_variogram_geometry_pgs_calcul(local_pgs, vario, idir)) return (1);
      st_variogram_geometry_pgs_correct(local_pgs, vario, idir);
      if (st_variogram_geometry_pgs_final(local_pgs)) return (1);
    }

    /* Set the value of C(0) */

    st_variogram_patch_C00(local_pgs, vario, idir, rule->getRho());

    if (ngrf > 1 && (opt_correl != 2 || rule->getRho() != 0))
      st_varcalc_correlated_grf(local_pgs, idir);
    else
      st_varcalc_uncorrelated_grf(local_pgs, idir);

    /* Clear the geometry */

    if (flag_geometry)
      local_pgs->vorder = vario_order_manage(0, 0, 0, local_pgs->vorder);
  }
  return (0);
}

/****************************************************************************/
/*!
 **  Make some lags inactive
 **
 ** \param[out]  vario    Vario structure
 **
 *****************************************************************************/
static void st_make_some_lags_inactive(Vario *vario)
{
  for (int idir = 0; idir < vario->getDirectionNumber(); idir++)
    for (int ipas = 0; ipas < vario->getLagNumber(idir); ipas++)
      vario->setUtilizeByIndex(idir, vario->getLagNumber(idir) + ipas, 1.);
}

/****************************************************************************/
/*!
 **  Make all lags active
 **
 ** \param[out]  vario    Vario structure
 **
 *****************************************************************************/
static void st_make_all_lags_active(Vario *vario)
{
  for (int idir = 0; idir < vario->getDirectionNumber(); idir++)
    for (int ipas = 0; ipas < vario->getLagNumber(idir); ipas++)
      vario->setUtilizeByIndex(idir, vario->getLagNumber(idir) + ipas, 1.);
}

/****************************************************************************/
/*!
 **  Local searching function for rho
 **
 ** \return  Evaluation value
 **
 ** \param[in]  rho        rho parameter
 ** \param[in]  user_data  User Data
 **
 *****************************************************************************/
static double st_rho_search(double rho, void *user_data)
{

  double sum;
  int idir, ndir;
  Local_Pgs *local_pgs;

  /* Initializations */

  sum = 0;
  local_pgs = (Local_Pgs*) user_data;
  ndir = local_pgs->vario->getDirectionNumber();
  st_set_rho(rho, local_pgs);

  /* Evaluation of the global cost-function */

  for (idir = 0; idir < ndir; idir++)
    sum += st_varcalc_correlated_grf(local_pgs, idir);

  if (debug_query("converge"))
    message("Value of the evaluating function = %lf - rho value %lf\n", sum,
            rho);

  return (sum);
}

/****************************************************************************/
/*!
 **  Performing the variogram calculations (in the case of flag.rho)
 **
 ** \return  Error return code
 **
 ** \param[in]  vario         Vario structure for the GRFs to be filled
 ** \param[in]  rule          Lithotype Rule definition
 ** \param[in]  local_pgs     Local_Pgs structure
 ** \param[in]  ngrf          Number of GRFs
 ** \param[in]  opt_correl    0 full model; 1 symetrical; 2 residuals
 **
 *****************************************************************************/
static int st_variopgs_calcul_rho(Vario *vario,
                                  const Rule *rule,
                                  Local_Pgs *local_pgs,
                                  int ngrf,
                                  int opt_correl)
{
  int idir;
  double testval, niter;

  /* Calculate the geometry */

  for (idir = 0; idir < vario->getDirectionNumber(); idir++)
  {
    if (st_variogram_geometry_pgs_calcul(local_pgs, vario, idir)) return (1);
    st_variogram_geometry_pgs_correct(local_pgs, vario, idir);
  }
  if (st_variogram_geometry_pgs_final(local_pgs)) return (1);

  st_make_some_lags_inactive(vario);
  st_set_rho(
      golden_search(st_rho_search, (void*) local_pgs, GS_TOLSTOP_RHO, -1., 1.,
                    &testval, &niter),
      local_pgs);
  st_make_all_lags_active(vario);

  /* Perform the calculations with fixed rho */

  if (st_variopgs_calcul_norho(vario, rule, local_pgs, ngrf, opt_correl, 0))
    return (1);

  /* Clean the geometry */

  local_pgs->vorder = vario_order_manage(0, 0, 0, local_pgs->vorder);

  return (0);
}

/****************************************************************************/
/*!
 **  Check if the Discrete Calculations make sens or not
 **
 ** \return  Error code
 **
 ** \param[in]  mode        Lithotype mode (ENUM_RULES)
 ** \param[in]  flag_rho    1 if rho has to be calculated, 0 otherwise
 **
 *****************************************************************************/
static int st_check_test_discret(const ERule &mode, int flag_rho)
{
  // Avoiding Discretizing calculation is always valid
  if (!TEST_DISCRET) return 0;

  // Case where the Discrete calculation option has been switch ON
  if (mode == ERule::STD)
  {
    // Only authorized when GRFs are not correlated
    if (flag_rho)
    {
      messerr("Calculations may not be perfored using Discretized Version");
      messerr("when underlying GRFs are correlated");
      return 1;
    }
  }
  else
  {
    messerr("Calculations may not be performed using Discretized version");
    messerr("when the Rule is not Standard (ERule::STD)");
    return 1;
  }
  return 0;
}

/****************************************************************************/
/*!
 **  Check that the arguments are correct
 **
 ** \return  Error return code
 **
 ** \param[in]  flag_db       1 if the input Db must be provided
 **                          -1 if it is optional
 **                           0 if the Db is not tested
 ** \param[in]  flag_rule     1 if the Rule must be defined
 ** \param[in]  flag_varioind 1 if the Indicator Variogram must be defined
 ** \param[in]  db            Db structure
 ** \param[in]  dbprop        Db Grid used for proportions (non-stationary)
 ** \param[in]  vario         Vario structure for the GRFs to be filled
 ** \param[in]  varioind      Vario structure for Indicator
 ** \param[in]  rule          Lithotype Rule definition
 **
 *****************************************************************************/
static int st_vario_pgs_check(int flag_db,
                              int flag_rule,
                              int flag_varioind,
                              Db *db,
                              const Db *dbprop,
                              const Vario *vario,
                              Vario *varioind,
                              const Rule *rule)
{

  /* Experimental variogram (compulsory) */

  if (vario == nullptr)
  {
    messerr("You must define the Input Variogram for the GRFs");
    return (1);
  }
  if (vario->getCalculType() != ECalcVario::COVARIANCE && vario->getCalculType()
      != ECalcVario::COVARIANCE_NC
      && vario->getCalculType() != ECalcVario::VARIOGRAM)
  {
    messerr("Only the Variogram is calculated here");
    return (1);
  }

  /* Input Db file (optional) */

  if (flag_db != 0)
  {
    if (flag_db > 0 && db == nullptr)
    {
      messerr("You must define the Input Db");
      return (1);
    }
    if (db != nullptr)
    {
      if (!db->isVariableNumberComparedTo(1)) return 1;
      if (db->getNDim() != vario->getDimensionNumber())
      {
        messerr("Space Dimension inconsistency between Input Db and Vario");
        return (1);
      }
    }
  }

  /* Rule (optional) */

  if (flag_rule)
  {
    if (rule == nullptr)
    {
      messerr("You must define the Rule");
      return (1);
    }
    if (rule->getModeRule() != ERule::STD)
    {
      messerr("This function is only programmed for standard rule");
      return (1);
    }
  }

  /* Optional Proportion File */

  if (dbprop != nullptr && dbprop->getNDim() != vario->getDimensionNumber())
  {
    messerr("Space Dimension inconsistency between Dbprop and Vario");
    return (1);
  }

  /* Indicator variogram (optional) */

  if (flag_varioind)
  {
    if (varioind == nullptr)
    {
      messerr("You must define the Indicator Variogram (stationary case)");
      return (1);
    }
  }
  return (0);
}

/****************************************************************************/
/*!
 **  Calculate the gaussian variograms
 **
 ** \return  Error return code
 **
 ** \param[in]  db           Db structure
 ** \param[in]  dbprop       Db Grid used for proportions (non-stationary)
 ** \param[in]  vario        Vario structure for the GRFs to be filled
 ** \param[in]  rule         Lithotype Rule definition
 ** \param[in]  propcst      Array of proportions for the facies
 ** \param[in]  flag_rho     1 if the correlation coefficient must be regressed
 ** \param[in]  opt_correl   0 full model; 1 symmetrical; 2 residuals
 **
 *****************************************************************************/
static int st_variogram_pgs_nostat(Db *db,
                                   const Db *dbprop,
                                   Vario *vario,
                                   const Rule *rule,
                                   const VectorDouble &propcst,
                                   int flag_rho,
                                   int opt_correl)
{
  Local_Pgs local_pgs;
  int flag_correl, flag_stat;
  int error, nfacies, ngrf;
  PropDef *propdef;

  /* Initializations */

  error = 1;
  ngrf = 0;
  flag_stat = nfacies = 0;
  propdef = nullptr;
  st_manage_pgs(0, &local_pgs);

  /* Preliminary checks */

  if (st_vario_pgs_check(1, 1, 0, db, dbprop, vario, NULL, rule))
    goto label_end;

  /*******************/
  /* Core allocation */
  /*******************/

  ngrf = rule->getGRFNumber();
  nfacies = rule->getFaciesNumber();
  propdef = proportion_manage(1, 1, flag_stat, ngrf, 0, nfacies, 0, db, dbprop,
                              propcst, propdef);
  if (propdef == nullptr) goto label_end;
  flag_correl = ngrf > 1 && (opt_correl != 2 || rule->getRho() != 0);
  if (rule->particularities(db, dbprop, NULL, 1, flag_stat)) goto label_end;
  proportion_rule_process(propdef, EProcessOper::COPY);

  /**************************/
  /* Allocate the variables */
  /**************************/

  if (st_vario_pgs_variable(1, ngrf, nfacies, 1, 0, db, propdef, rule))
    goto label_end;

  /****************************/
  /* Perform the calculations */
  /****************************/

  /* Initialize the Local_Pgs structure */

  st_manage_pgs(1, &local_pgs, db, rule, vario, nullptr, nullptr, propdef,
                flag_stat, 1, 0, ngrf, nfacies, vario->getCalculType());
  st_define_corpgs(opt_correl, flag_rho, rule->getRho(), &local_pgs);
  st_define_trace(flag_rho, flag_correl, &local_pgs);

  /* Infer the variogram of PGS */

  if (st_vario_pgs_variable(0, ngrf, nfacies, 1, 0, db, propdef, rule))
    goto label_end;
  if (!flag_rho)
  {
    st_set_rho(rule->getRho(), &local_pgs);
    if (st_variopgs_calcul_norho(vario, rule, &local_pgs, ngrf, opt_correl, 1))
      goto label_end;
  }
  else
  {
    if (st_variopgs_calcul_rho(vario, rule, &local_pgs, ngrf, opt_correl))
      goto label_end;
  }

  /* Set the error return flag */

  error = 0;

  label_end: (void) st_extract_trace(&local_pgs);
  st_manage_pgs(-1, &local_pgs, db, rule, vario, nullptr, nullptr, propdef,
                flag_stat, 1, 0, ngrf, nfacies, vario->getCalculType());
  (void) st_vario_pgs_variable(-1, ngrf, nfacies, 1, 0, db, propdef, rule);
  propdef = proportion_manage(-1, 1, flag_stat, ngrf, 0, nfacies, 0, db, dbprop,
                              propcst, propdef);
  return (error);
}

/****************************************************************************/
/*!
 **  Calculate the covariance matrix
 **
 ** \param[in]  local_pgs  Local_Pgs structure
 **
 ** \param[out] flag_ind   1 if the two GRF are independent
 ** \param[out] iconf      Array of ranks of the discretized covariance
 ** \param[out] cov        Matrix of covariance
 **
 *****************************************************************************/
static void st_calcul_covmatrix(Local_Pgs *local_pgs,
                                int *flag_ind,
                                int *iconf,
                                double *cov)
{
  double cov0[4], covh[4], cround;
  CovCalcMode mode;

  const Rule *rule = local_pgs->rule;
  int nvar = local_pgs->model->getVariableNumber();
  int ngrf = local_pgs->rule->getGRFNumber();

  /* Calculate the covariance for the zero distance */
  for (int i = 0; i < local_pgs->model->getDimensionNumber(); i++)
    local_pgs->d0[i] = 0.;
  model_calcul_cov(NULL,local_pgs->model, mode, 1, 1., local_pgs->d0, cov0);

  /* Calculate the covariance for the given shift */
  model_calcul_cov(NULL,local_pgs->model, mode, 1, 1., local_pgs->d1, covh);

  if (rule->getModeRule() == ERule::STD)
  {
    cov[0] = covh[0]; /* C11(h)  */
    if (ngrf > 1)
    {
      cov[1] = cov0[1]; /* C21(0)  */
      cov[2] = covh[2]; /* C21(-h) */
      cov[3] = covh[1]; /* C21(h)  */
      cov[4] = cov0[2]; /* C21(0)  */
      cov[5] = covh[3]; /* C22(h)  */
    }
  }
  else if (rule->getModeRule() == ERule::SHIFT)
  {
    RuleShift *ruleshift = (RuleShift*) rule;
    cov[0] = covh[0]; /* C11(h)  */
    cov[5] = (nvar == 1) ? covh[0] :
                           covh[3]; /* C22(h)  */

    for (int i = 0; i < local_pgs->model->getDimensionNumber(); i++)
      local_pgs->d0[i] = ruleshift->getShift(i);

    model_calcul_cov(NULL,local_pgs->model, mode, 1, 1., local_pgs->d0, covh);
    cov[1] = (nvar == 1) ? covh[0] :
                           covh[1]; /* C21(s)  */
    cov[4] = (nvar == 1) ? covh[0] :
                           covh[1]; /* C21(s)  */

    for (int i = 0; i < local_pgs->model->getDimensionNumber(); i++)
      local_pgs->d0[i] = local_pgs->d1[i] - ruleshift->getShift(i);
    model_calcul_cov(NULL,local_pgs->model, mode, 1, 1., local_pgs->d0, covh);
    cov[2] = (nvar == 1) ? covh[0] :
                           covh[1]; /* C21(h-s) */

    for (int i = 0; i < local_pgs->model->getDimensionNumber(); i++)
      local_pgs->d0[i] = local_pgs->d1[i] + ruleshift->getShift(i);
    model_calcul_cov(NULL,local_pgs->model, mode, 1, 1., local_pgs->d0, covh);
    cov[3] = (nvar == 1) ? covh[0] :
                           covh[1]; /* C21(h+s)  */
  }
  else
    messageAbort("This rule is not expected in st_calcul_covmatrix");

  /* Check if the two GRFs are spatially independent */

  (*flag_ind) = 1;
  if (ngrf > 1)
  {
    for (int i = 1; i <= 4; i++)
      if (ABS(cov[i]) > 1.e-8) (*flag_ind) = 0;
  }

  /* In TEST_DISCRET case, identify the ranks of the discretized covariance */

  if (TEST_DISCRET)
  {
    iconf[0] = ct_tableone_covrank(CTABLES, cov[0], &cround);
    if (local_pgs->ngrf > 1)
      iconf[1] = ct_tableone_covrank(CTABLES, cov[5], &cround);
  }
  return;
}

/****************************************************************************/
/*!
 **  Calculate the probability
 **
 ** \return  Evaluation value
 **
 ** \param[in] local_pgs Local_Pgs structure
 ** \param[in] flag_ind  1 if the GRFs are independent
 ** \param[in] low       Array of lower thresholds (Dimension: 4)
 ** \param[in] up        Array of upper thresholds (Dimension: 4)
 ** \param[in] iconf     Array of ranks of the discretized covariance
 ** \param[in] cov       Covariance matrix
 **
 *****************************************************************************/
static double st_get_proba(Local_Pgs *local_pgs,
                           int flag_ind,
                           double *low,
                           double *up,
                           int *iconf,
                           double *cov)
{
  double proba = TEST;
  int ngrf = local_pgs->ngrf;

  if (ngrf == 1)
  {
    proba = st_get_proba_ind(cov[0], low, up, iconf[0]); // TODO: to be checked
  }
  else
  {
    if (flag_ind)
    {
      proba = st_get_proba_ind(cov[0], &low[0], &up[0], iconf[0])
          * st_get_proba_ind(cov[5], &low[2], &up[2], iconf[1]);
    }
    else
    {
      // Case when the two GRFs are correlated (presence of RHO)

      if (!TEST_DISCRET)
      {
        int ier;
        double err;
        double releps = 0.;
        double abseps = EPS;
        int maxpts = 8000;

        int infin[4];
        infin[0] = mvndst_infin(low[0], up[0]);
        infin[1] = mvndst_infin(low[1], up[1]);
        infin[2] = mvndst_infin(low[2], up[2]);
        infin[3] = mvndst_infin(low[3], up[3]);
        mvndst(4, low, up, infin, cov, maxpts, abseps, releps, &err, &proba,
               &ier);
      }
      else
      {
        my_throw(
            "Discrete calculation for correlated GRFs is not performed yet");
      }
    }
  }

  return (proba);
}

/****************************************************************************/
/*!
 **  Define the thresholds of the GRF(s)
 **
 ** \param[in]  local_pgs     Local_Pgs structure
 ** \param[in]  iech1         Rank of the first sample
 ** \param[in]  iech2         Rank of the second sample
 ** \param[in]  ifac1         Rank of the first facies
 ** \param[in]  ifac2         Rank of the second facies
 **
 ** \param[out] low           Array of lower thresholds (Dimension: 4)
 ** \param[out] up            Array of upper thresholds (Dimension: 4)
 ** \param[out] ploc          Array of proportions (Dimension: 2)
 **
 ** REMARKS: Warning: in the case of TEST_DISCRET, the returned arguments
 ** REMARKS: 'low' and 'up' return the ranks in the discretized covariance
 **
 *****************************************************************************/
static void st_define_bounds(Local_Pgs *local_pgs,
                             int iech1,
                             int iech2,
                             int ifac1,
                             int ifac2,
                             double *low,
                             double *up,
                             double *ploc)
{
  int nfacies;

  nfacies = local_pgs->nfacies;
  if (local_pgs->flag_stat)
  {
    ploc[0] = local_pgs->propdef->propfix[ifac1];
    ploc[1] = local_pgs->propdef->propfix[ifac2];
    low[0] = STAT_THRESH(ifac1, 0, 0);
    up[0] = STAT_THRESH(ifac1, 0, 1);
    low[1] = STAT_THRESH(ifac2, 0, 0);
    up[1] = STAT_THRESH(ifac2, 0, 1);
    if (local_pgs->ngrf > 1)
    {
      low[2] = STAT_THRESH(ifac1, 1, 0);
      up[2] = STAT_THRESH(ifac1, 1, 1);
      low[3] = STAT_THRESH(ifac2, 1, 0);
      up[3] = STAT_THRESH(ifac2, 1, 1);
    }
  }
  else
  {
    ploc[0] = local_pgs->db->getProportion(iech1, ifac1);
    ploc[1] = local_pgs->db->getProportion(iech2, ifac2);

    if (!TEST_DISCRET)
    {
      low[0] = local_pgs->db->getLowerBound(iech1, ifac1);
      up[0] = local_pgs->db->getUpperBound(iech1, ifac1);
      low[1] = local_pgs->db->getLowerBound(iech2, ifac2);
      up[1] = local_pgs->db->getUpperBound(iech2, ifac2);
      if (local_pgs->ngrf > 1)
      {
        low[2] = local_pgs->db->getLowerBound(iech1, nfacies + ifac1);
        up[2] = local_pgs->db->getUpperBound(iech1, nfacies + ifac1);
        low[3] = local_pgs->db->getLowerBound(iech2, nfacies + ifac2);
        up[3] = local_pgs->db->getUpperBound(iech2, nfacies + ifac2);
      }
    }
    else
    {
      low[0] = local_pgs->db->getLowerInterval(iech1, ifac1);
      up[0] = local_pgs->db->getUpperInterval(iech1, ifac1);
      low[1] = local_pgs->db->getLowerInterval(iech2, ifac2);
      up[1] = local_pgs->db->getUpperInterval(iech2, ifac2);
      if (local_pgs->ngrf > 1)
      {
        low[2] = local_pgs->db->getLowerInterval(iech1, nfacies + ifac1);
        up[2] = local_pgs->db->getUpperInterval(iech1, nfacies + ifac1);
        low[3] = local_pgs->db->getLowerInterval(iech2, nfacies + ifac2);
        up[3] = local_pgs->db->getUpperInterval(iech2, nfacies + ifac2);
      }
    }
  }
  return;
}

/****************************************************************************/
/*!
 **  Calculate the variogram or covariance contribution
 **
 ** \return  Evaluation value
 **
 ** \param[in] local_pgs Local_Pgs structure
 ** \param[in] flag_ind  1 if the GRFs are independent
 ** \param[in] iech1     Rank of the first sample
 ** \param[in] iech2     Rank of the second sample
 ** \param[in] ifac1     Rank of the first facies
 ** \param[in] ifac2     Rank of the second facies
 ** \param[in] iconf     Array of ranks of the discretized covariance
 ** \param[in] cov       Covariance matrix
 **
 *****************************************************************************/
static double st_get_value(Local_Pgs *local_pgs,
                           int flag_ind,
                           int iech1,
                           int iech2,
                           int ifac1,
                           int ifac2,
                           int *iconf,
                           double *cov)
{
  double value, g1, g2, ploc[2], low[4], up[4];

  if (local_pgs->calcul_type == ECalcVario::VARIOGRAM)
  {
    if (ifac1 == ifac2)
    {
      st_define_bounds(local_pgs, iech1, iech2, ifac1, ifac2, low, up, ploc);
      g1 = st_get_proba(local_pgs, flag_ind, low, up, iconf, cov);
      value = (ploc[0] + ploc[1]) * 0.5 - g1;
    }
    else
    {
      st_define_bounds(local_pgs, iech1, iech2, ifac1, ifac2, low, up, ploc);
      g1 = st_get_proba(local_pgs, flag_ind, low, up, iconf, cov);
      st_define_bounds(local_pgs, iech2, iech1, ifac1, ifac2, low, up, ploc);
      g2 = st_get_proba(local_pgs, flag_ind, low, up, iconf, cov);
      value = -0.5 * (g1 + g2);
    }
  }
  else
  {
    st_define_bounds(local_pgs, iech1, iech2, ifac1, ifac2, low, up, ploc);
    value = st_get_proba(local_pgs, flag_ind, low, up, iconf, cov);
  }
  return (value);
}

/****************************************************************************/
/*!
 **  Performing the variogram calculations in the non-stationary case
 **
 ** \return  Error return code
 **
 ** \param[in]  local_pgs     Local_Pgs structure
 **
 *****************************************************************************/
static int st_vario_indic_model_nostat(Local_Pgs *local_pgs)

{
  double dist, cov[6];
  int ipas, ifac, jfac, nfacies, ipair, iech, jech, i, idir, flag_ind, iconf[2];
  Vario *vario;

  /* Initializations */

  nfacies = local_pgs->nfacies;
  vario = local_pgs->vario;

  /* Loop on the directions */

  for (idir = 0; idir < vario->getDirectionNumber(); idir++)
  {
    /* Establish the geometry */

    if (st_variogram_geometry_pgs_calcul(local_pgs, vario, idir)) return (1);
    if (st_variogram_geometry_pgs_final(local_pgs)) return (1);

    /* Loop on the lags */

    for (ipas = 0; ipas < vario->getLagNumber(idir); ipas++)
    {
      vario_order_get_bounds(local_pgs->vorder, idir, ipas, &local_pgs->ifirst,
                             &local_pgs->ilast);
      if (local_pgs->ifirst >= local_pgs->ilast) continue;

      /* Loop on the pairs of the lag */

      for (ipair = local_pgs->ifirst; ipair < local_pgs->ilast; ipair++)
      {
        vario_order_get_indices(local_pgs->vorder, ipair, &iech, &jech, &dist);

        /* Calculate the distance vector */

        dist = distance_intra(local_pgs->db, iech, jech, local_pgs->d1.data());
        st_calcul_covmatrix(local_pgs, &flag_ind, iconf, cov);

        /* Loops on the facies */

        for (ifac = 0; ifac < nfacies; ifac++)
          for (jfac = 0; jfac <= ifac; jfac++)
          {
            if (local_pgs->vario->getFlagAsym())
            {
              i = vario->getDirAddress(idir, ifac, jfac, ipas, false, 1);
              vario->setGgByIndex(
                  idir,
                  i,
                  vario->getGgByIndex(idir, i) + st_get_value(local_pgs, flag_ind,
                                                       iech, jech, ifac, jfac,
                                                       iconf, cov));
              i = vario->getDirAddress(idir, ifac, jfac, ipas, false, -1);
              vario->setGgByIndex(
                  idir,
                  i,
                  vario->getGgByIndex(idir, i) + st_get_value(local_pgs, flag_ind,
                                                       jech, iech, ifac, jfac,
                                                       iconf, cov));
            }
            else
            {
              i = vario->getDirAddress(idir, ifac, jfac, ipas, false, 0);
              vario->setGgByIndex(
                  idir,
                  i,
                  vario->getGgByIndex(idir, i) + st_get_value(local_pgs, flag_ind,
                                                       iech, jech, ifac, jfac,
                                                       iconf, cov));
            }
          }
      }
    }

    /* Clear the geometry */

    local_pgs->vorder = vario_order_manage(0, 0, 0, local_pgs->vorder);

    /* Scale the variogram */

    variogram_scale(vario, idir);
  }
  return (0);
}

/****************************************************************************/
/*!
 **  Duplicate the information from the input variogram to the output variogram
 **  This operation considers (Sw, Hh, Gg) and replicates this information for
 **  all directions.
 **  Per direction, the information of the first simple input variogram is
 **  replicated to all simple and cross-variograms output variograms
 **
 ** \return  Error return code
 **
 ** \param[in]  vario1       Input Variogram
 ** \param[out] vario2       Output Variogram
 ** \param[in]  flagSw       True if Sw must be replicated
 ** \param[in]  flagHh       True if Hh must be replicated
 ** \param[in]  flagGg       True if Gg must be replicated
 **
 *****************************************************************************/
static int st_copy_swhh(const Vario *vario1,
                        Vario *vario2,
                        bool flagSw,
                        bool flagHh,
                        bool flagGg)
{
  if (vario1->getDirectionNumber() != vario2->getDirectionNumber())
  {
    messerr("Both variograms should share the same number of Directions");
    return 1;
  }
  int nvar = vario2->getVariableNumber();

  for (int idir = 0; idir < vario2->getDirectionNumber(); idir++)
  {
    for (int ipas = 0; ipas < vario1->getLagTotalNumber(idir); ipas++)
    {
      int iad1 = vario1->getDirAddress(idir, 0, 0, ipas, true, 0);
      for (int ivar = 0; ivar < nvar; ivar++)
        for (int jvar = 0; jvar < nvar; jvar++)
        {
          int iad2 = vario2->getDirAddress(idir, ivar, jvar, ipas, true, 0);
          if (flagSw) vario2->setSwByIndex(idir, iad2, vario1->getSwByIndex(idir, iad1));
          if (flagHh) vario2->setHhByIndex(idir, iad2, vario1->getHhByIndex(idir, iad1));
          if (flagGg) vario2->setGgByIndex(idir, iad2, vario1->getGgByIndex(idir, iad1));
        }
    }
  }
  return 0;
}

/****************************************************************************/
/*!
 **  Performing the variogram calculations in the stationary case
 **
 ** \return  Error return code
 **
 ** \param[in]  local_pgs     Local_Pgs structure
 **
 *****************************************************************************/
static int st_vario_indic_model_stat(Local_Pgs *local_pgs)

{
  double cov[6];
  int ii, flag_ind, iconf[2];

  /* Initializations */

  int nfacies = local_pgs->nfacies;
  Vario *vario = local_pgs->vario;
  for (int i = 0; i < 6; i++)
    cov[i] = 0.;

  // Duplicate Number and Distance for all lags (from the first simple variogram)

  if (st_copy_swhh(local_pgs->varioind, local_pgs->vario, true, true, false))
    return 1;

  /* Loop on the directions */

  for (int idir = 0; idir < vario->getDirectionNumber(); idir++)
  {

    /* Loop on the lags */

    for (int ipas = 0; ipas < vario->getLagNumber(idir); ipas++)
    {

      /* Calculate the distance vector */

      for (int i = 0; i < vario->getDimensionNumber(); i++)
      {
        int jpas = vario->getDirAddress(idir, 0, 0, ipas, false, 1);
        local_pgs->d1[i] = vario->getHhByIndex(idir, jpas) * vario->getCodir(idir, i);
      }
      st_calcul_covmatrix(local_pgs, &flag_ind, iconf, cov);

      /* Loops on the facies */

      for (int ifac = 0; ifac < nfacies; ifac++)
        for (int jfac = 0; jfac <= ifac; jfac++)
        {
          if (local_pgs->vario->getFlagAsym())
          {
            ii = vario->getDirAddress(idir, ifac, jfac, ipas, false, 1);
            vario->setGgByIndex(
                idir,
                ii,
                st_get_value(local_pgs, flag_ind, 0, 0, ifac, jfac, iconf,
                             cov));
            ii = vario->getDirAddress(idir, ifac, jfac, ipas, false, -1);
            vario->setGgByIndex(
                idir,
                ii,
                st_get_value(local_pgs, flag_ind, 0, 0, jfac, ifac, iconf,
                             cov));
          }
          else
          {
            ii = vario->getDirAddress(idir, ifac, jfac, ipas, false, 0);
            vario->setGgByIndex(
                idir,
                ii,
                st_get_value(local_pgs, flag_ind, 0, 0, ifac, jfac, iconf,
                             cov));
          }
        }
    }
  }
  return (0);
}

/****************************************************************************/
/*!
 **  Establish the theoretical variance of the simple and
 **  cross-variograms of the indicators in the stationary case
 **
 ** \param[in]  local_pgs  Local_Pgs structure
 **
 *****************************************************************************/
static void st_update_variance_stat(Local_Pgs *local_pgs)

{
  int ivar, jvar, idir, iad, nfacies;
  double pivar, pjvar;
  Vario *vario;

  /* Initializations */

  vario = local_pgs->vario;
  nfacies = local_pgs->nfacies;

  /* Evaluate the theoretical variances */

  for (ivar = 0; ivar < nfacies; ivar++)
    for (jvar = 0; jvar < nfacies; jvar++)
    {
      pivar = local_pgs->propdef->propfix[ivar];
      pjvar = local_pgs->propdef->propfix[jvar];
      if (ivar == jvar)
        vario->setVar(ivar, jvar, pivar * (1. - pivar));
      else
        vario->setVar(ivar, jvar, -pivar * pjvar);
      if (!vario->getFlagAsym()) continue;

      for (idir = 0; idir < vario->getDirectionNumber(); idir++)
      {
        iad = vario->getDirAddress(idir, ivar, jvar, 0, false, 0);
        vario->setSwByIndex(idir, iad, 1);
        vario->setHhByIndex(idir, iad, 0);

        switch (local_pgs->calcul_type.toEnum())
        {
          case ECalcVario::E_VARIOGRAM:
            break;

          case ECalcVario::E_COVARIANCE:
            vario->setGgByIndex(idir, iad, vario->getVar(ivar, jvar));
            break;

          case ECalcVario::E_COVARIANCE_NC:
            vario->setGgByIndex(idir, iad, (ivar == jvar) ? pivar :
                                                     0.);
            break;
          default:
            break;
        }
      }
    }
  return;
}

/****************************************************************************/
/*!
 **  Establish the theoretical variance of the simple and
 **  cross-variograms of the indicators in the non-stationary case
 **
 ** \return  Error return code
 **
 ** \param[in]  local_pgs  Local_Pgs structure
 **
 *****************************************************************************/
static int st_update_variance_nostat(Local_Pgs *local_pgs)

{
  double *mean, *covs;

  /* Initializations */

  int error = 1;
  int number = 0;
  Db *dbin = local_pgs->db;
  Vario *vario = local_pgs->vario;
  int nfacies = local_pgs->nfacies;
  mean = covs = nullptr;

  /* Core allocation */

  mean = (double*) mem_alloc(nfacies * sizeof(double), 0);
  if (mean == nullptr) goto label_end;
  for (int i = 0; i < nfacies; i++)
    mean[i] = 0.;
  covs = (double*) mem_alloc(nfacies * nfacies * sizeof(double), 0);
  if (covs == nullptr) goto label_end;
  for (int i = 0; i < nfacies * nfacies; i++)
    covs[i] = 0.;

  /* Loop on the samples */

  for (int iech = 0; iech < dbin->getSampleNumber(); iech++)
  {
    if (!dbin->isActive(iech)) continue;

    /* Loop on the variables */

    for (int ivar = 0; ivar < nfacies; ivar++)
    {
      double p1 = dbin->getProportion(iech, ivar);
      mean[ivar] += p1;

      for (int jvar = 0; jvar < nfacies; jvar++)
      {
        double p2 = dbin->getProportion(iech, jvar);
        COVS(ivar,jvar) += p1 * p2;
      }
    }
    number++;
  }

  /* Normalization */

  for (int ivar = 0; ivar < nfacies; ivar++)
    mean[ivar] /= (double) number;
  for (int ivar = 0; ivar < nfacies; ivar++)
    for (int jvar = 0; jvar < nfacies; jvar++)
      COVS(ivar,jvar) = COVS(ivar,jvar) / (double) number
          - mean[ivar] * mean[jvar];

  /* Store the results */

  for (int ivar = 0; ivar < nfacies; ivar++)
    for (int jvar = 0; jvar < nfacies; jvar++)
    {
      vario->setVar(ivar, jvar, COVS(ivar, jvar));

      if (!vario->getFlagAsym()) continue;

      // Set the C00 term
      for (int idir = 0; idir < vario->getDirectionNumber(); idir++)
      {
        int iad = vario->getDirAddress(idir, ivar, jvar, 0, false, 0);
        vario->setSwByIndex(idir, iad, dbin->getSampleNumber());
        vario->setHhByIndex(idir, iad, 0);

        switch (local_pgs->calcul_type.toEnum())
        {
          case ECalcVario::E_VARIOGRAM:
            break;

          case ECalcVario::E_COVARIANCE:
            vario->setGgByIndex(idir, iad, vario->getVar(ivar, jvar));
            break;

          case ECalcVario::E_COVARIANCE_NC:
            vario->setGgByIndex(idir, iad, (ivar == jvar) ? mean[ivar] :
                                                     0.);
            break;
          default:
            break;
        }
      }
    }

  /* Set the error return code */

  error = 0;

  label_end: mean = (double*) mem_free((char* ) mean);
  covs = (double*) mem_free((char* ) covs);
  return (error);
}

/****************************************************************************/
/*!
 **  Evaluate the experimental variogram of indicators in PluriGaussian case
 **
 ** \return  Error return code
 **
 ** \param[in]  db         Db descriptor
 ** \param[in]  varioparam VarioParam structure
 ** \param[in]  ruleprop   RuleProp structure
 ** \param[in]  model1     First Model structure
 ** \param[in]  model2     Second Model structure (optional)
 **
 *****************************************************************************/
Vario* model_pgs(Db *db,
                                 const VarioParam *varioparam,
                                 const RuleProp *ruleprop,
                                 const Model *model1,
                                 const Model *model2)
{
  Vario *vario = nullptr;
  Vario *varioind = nullptr;

  if (varioparam == nullptr)
  {
    messerr("The VarioParam must be provided");
    return nullptr;
  }
  if (ruleprop == nullptr)
  {
    messerr("RuleProp must be defined");
    return nullptr;
  }
  int flag_stat = ruleprop->isFlagStat();
  const Rule *rule = ruleprop->getRule();
  const VectorDouble &propcst = ruleprop->getPropCst();
  const Db *dbprop = ruleprop->getDbprop();

  Local_Pgs local_pgs;
  PropDef *propdef;
  Model *new_model;

  /*******************/
  /* Initializations */
  /*******************/

  int error = 1;
  int nfacies = rule->getFaciesNumber();
  int ngrf = rule->getGRFNumber();
  if (rule->getModeRule() == ERule::SHIFT) ngrf++;

  new_model = nullptr;
  propdef = nullptr;
  st_manage_pgs(0, &local_pgs);
  if (st_check_test_discret(rule->getModeRule(), 0)) goto label_end;

  /* Merge the models */

  new_model = model_rule_combine(model1, model2, rule);
  if (new_model == nullptr)
  {
    messerr("The Model(s) must be defined");
    return nullptr;
  }
  if (new_model->getVariableNumber() != ngrf)
  {
    messerr("The number of GRF is not equal to the number of variables");
    messerr("defined in the combined Model");
    return nullptr;
  }

  if (!flag_stat)
  {
    if (db == nullptr)
    {
      messerr("You must define the Input Db");
      return nullptr;
    }
    if (db->getNDim() != varioparam->getDimensionNumber())
    {
      messerr("Inconsistent parameters:");
      messerr("Input DB : NDIM=%d", db->getNDim());
      messerr("Variogram: NDIM=%d", varioparam->getDimensionNumber());
      return nullptr;
    }
    if (dbprop != nullptr && dbprop->getNDim()
        != varioparam->getDimensionNumber())
    {
      messerr("Space Dimension inconsistency between Dbprop and Vario");
      return nullptr;
    }
    if (new_model->getDimensionNumber() != db->getNDim())
    {
      messerr("The Space Dimension of the Db structure (%d)", db->getNDim());
      messerr("Does not correspond to the Space Dimension of the model (%d)",
              new_model->getDimensionNumber());
      return nullptr;
    }
  }

  // Initiate the output class

  vario = new Vario(varioparam, db);
  vario->setCalculName("vg");
  vario->setNVar(nfacies);
  vario->internalVariableResize();
  vario->internalDirectionResize();

  // Calculate the variogram of Indicators
  if (flag_stat)
  {
    varioind = new Vario(varioparam, db);
    if (varioind->computeIndic("vg")) return nullptr;
  }

  /* Preliminary checks */

  if (st_vario_pgs_check(-1, 1, flag_stat, db, dbprop, vario, varioind, rule))
    goto label_end;

  /*******************/
  /* Core allocation */
  /*******************/

  propdef = proportion_manage(1, 1, flag_stat, ngrf, 0, nfacies, 0, db, dbprop,
                              propcst, propdef);
  if (propdef == nullptr) goto label_end;

  if (rule->particularities(db, dbprop, new_model, 0, flag_stat))
    goto label_end;

  proportion_rule_process(propdef, EProcessOper::COPY);

  /* Pre-calculation of integrals: Define the structure */

  if (TEST_DISCRET)
    CTABLES = ct_tables_manage(1, 0, 1, 200, 100, -1., 1., NULL);

  st_manage_pgs(1, &local_pgs, db, rule, vario, varioind, new_model, propdef,
                flag_stat, 0, 1, ngrf, nfacies, vario->getCalculType());

  /* Calculate the variogram and the variance matrix */

  if (flag_stat)
  {
    if (st_calculate_thresh_stat(&local_pgs)) goto label_end;
    if (st_vario_indic_model_stat(&local_pgs)) goto label_end;
    st_update_variance_stat(&local_pgs);
  }
  else
  {
    if (st_vario_pgs_variable(1, ngrf, nfacies, 0, 1, db, propdef, rule))
      goto label_end;
    if (st_vario_pgs_variable(0, ngrf, nfacies, 0, 1, db, propdef, rule))
      goto label_end;
    if (st_vario_indic_model_nostat(&local_pgs)) goto label_end;
    if (st_update_variance_nostat(&local_pgs)) goto label_end;
  }

  /* Set the error return code */

  error = 0;

  label_end: if (TEST_DISCRET)
    CTABLES = ct_tables_manage(-1, 0, 1, 200, 100, -1., 1., CTABLES);
  st_manage_pgs(-1, &local_pgs, db, rule, vario, varioind, new_model, propdef,
                flag_stat, 0, 1, ngrf, nfacies, vario->getCalculType());
  new_model = model_free(new_model);
  (void) st_vario_pgs_variable(-1, ngrf, nfacies, 0, 1, db, propdef, rule);
  propdef = proportion_manage(-1, 1, flag_stat, ngrf, 0, nfacies, 0, db, dbprop,
                              propcst, propdef);
  if (error)
  {
    delete vario;
    vario = nullptr;
  }
  return vario;
}

/****************************************************************************/
/*!
 **  Calculate the gaussian variograms in the stationary case
 **
 ** \return  Error return code
 **
 ** \param[in]  db           Db structure (for Space dimension)
 ** \param[in]  vario        Vario structure for the GRFs to be filled
 ** \param[in]  varioind     Indicator Vario structure
 ** \param[in]  rule         Lithotype Rule definition
 ** \param[in]  propcst      Array of proportions for the facies
 **
 *****************************************************************************/
static int st_variogram_pgs_stat(Db *db,
                                 Vario *vario,
                                 Vario *varioind,
                                 const Rule *rule,
                                 const VectorDouble &propcst)
{
  Local_Pgs local_pgs;
  int node_tot, nmax_tot, ny1, ny2, error, nfacies, ngrf, flag_stat;
  double prop_tot;
  PropDef *propdef;

  /* Initializations */

  error = 1;
  ngrf = nfacies = 0;
  flag_stat = 1;
  propdef = nullptr;
  st_manage_pgs(0, &local_pgs);

  /* Preliminary checks */

  if (st_vario_pgs_check(0, 1, 1, db, NULL, vario, varioind, rule))
    goto label_end;
  ngrf = rule->getGRFNumber();
  nfacies = rule->getFaciesNumber();

  /*******************/
  /* Core allocation */
  /*******************/

  rule->statistics(0, &node_tot, &nfacies, &nmax_tot, &ny1, &ny2, &prop_tot);
  propdef = proportion_manage(1, 1, flag_stat, ngrf, 0, nfacies, 0,
  NULL,
                              NULL, propcst, propdef);
  if (propdef == nullptr) goto label_end;
  if (rule->particularities(NULL, NULL, NULL, 1, flag_stat)) goto label_end;
  proportion_rule_process(propdef, EProcessOper::COPY);

  /****************************/
  /* Perform the calculations */
  /****************************/

  /* Initialize the Local_Pgs structure */

  st_manage_pgs(1, &local_pgs, db, rule, vario, varioind, nullptr, propdef,
                flag_stat, 1, 0, ngrf, nfacies, vario->getCalculType());
  st_define_corpgs(0, 0, rule->getRho(), &local_pgs);
  st_define_trace(0, 0, &local_pgs);
  st_set_rho(0., &local_pgs);
  if (st_calculate_thresh_stat(&local_pgs)) goto label_end;

  // Copy the count and distance information from indicator variogram into calculation one

  if (st_copy_swhh(local_pgs.varioind, local_pgs.vario, true, true, false))
    goto label_end;

  /* Infer the variogram of PGS */

  st_varcalc_from_vario_stat(vario, &local_pgs, ngrf);

  /* Set the error return flag */

  error = 0;

  label_end: (void) st_extract_trace(&local_pgs);
  st_manage_pgs(-1, &local_pgs, db, rule, vario, varioind, NULL, propdef,
                flag_stat, 1, 0, ngrf, nfacies, vario->getCalculType());
  propdef = proportion_manage(-1, 1, 1, ngrf, 0, nfacies, 0,
  NULL,
                              NULL, propcst, propdef);
  return (error);
}

/****************************************************************************/
/*!
 **  Calculate the gaussian variograms
 **
 ** \return  Error return code
 **
 ** \param[in]  db           Db structure
 ** \param[in]  varioparam   VarioParam structure for the GRFs
 ** \param[in]  ruleprop     RuleProp structure
 ** \param[in]  flag_rho     1 if the correlation coefficient must be regressed
 ** \param[in]  opt_correl   0 full model; 1 symmetrical; 2 residuals
 **
 ** \remarks This is simply a routine dispatching between the stationary function
 ** \remarks and the non-stationary one
 **
 *****************************************************************************/
Vario* variogram_pgs(Db *db,
                                     const VarioParam *varioparam,
                                     const RuleProp *ruleprop,
                                     int flag_rho,
                                     int opt_correl)
{
  Vario *vario = nullptr;
  Vario *varioind = nullptr;

  if (db == NULL)
  {
    messerr("The Db must be provided");
    return nullptr;
  }
  if (varioparam == nullptr)
  {
    messerr("The VarioParam must be provided");
    return nullptr;
  }
  if (ruleprop == nullptr)
  {
    messerr("RuleProp must be defined");
    return nullptr;
  }
  int flag_stat = ruleprop->isFlagStat();
  const Rule *rule = ruleprop->getRule();
  const VectorDouble &propcst = ruleprop->getPropCst();
  const Db *dbprop = ruleprop->getDbprop();

  int error;

  // Preliminary checks

  if (st_check_test_discret(rule->getModeRule(), flag_rho)) return nullptr;
  if (varioparam->getDirectionNumber() <= 0)
  {
    messerr("The variogram must contain at least one calculation Direction");
    return nullptr;
  }
  if (db->getVariableNumber() != 1)
  {
    messerr("The number of variables (%d) must be equal to 1",
            db->getVariableNumber());
    return nullptr;
  }
  int nclass = rule->getFaciesNumber();
  if (nclass <= 0)
  {
    messerr("No Facies class have been found");
    return nullptr;
  }

  // In Stationary case, create the variogram of indicators to speed up calculations

  VectorDouble props;
  if (flag_stat)
  {
    if (!propcst.empty())
    {
      if ((int) propcst.size() != nclass)
      {
        messerr(
            "Number of proportions in 'propcst' (%d) should match Number of Facies in 'rule' (%d)",
            (int) propcst.size(), rule->getFaciesNumber());
        return nullptr;
      }
      props = propcst;
    }
    else
    {
      // Calculate the number of Facies in 'Db'
      props = dbStatisticsFacies(db);
      if ((int) props.size() != nclass)
      {
        messerr(
            "Number of Facies in 'db' (%d) should match Number of facies in 'rule' (%d)",
            (int) props.size(), rule->getFaciesNumber());
        return nullptr;
      }
    }

    // Calculate the variogram of Indicators
    varioind = new Vario(varioparam, db);
    if (varioind->computeIndic("covnc")) return nullptr;
  }

  /* Pre-calculation of integrals: Define the structure */

  if (TEST_DISCRET)
    CTABLES = ct_tables_manage(1, 0, 1, 200, 100, -1., 1., NULL);

  // Initiate the output class

  vario = new Vario(varioparam, db);
  vario->setCalculName("covnc");
  vario->setNVar(rule->getGRFNumber());
  vario->internalVariableResize();
  vario->internalDirectionResize();

  /* Perform the calculations */

  if (flag_stat)
    error = st_variogram_pgs_stat(db, vario, varioind, rule, props);
  else
    error = st_variogram_pgs_nostat(db, dbprop, vario, rule, props, flag_rho,
                                    opt_correl);

  /* Final operations */

  if (TEST_DISCRET)
    CTABLES = ct_tables_manage(-1, 0, 1, 200, 100, -1., 1., CTABLES);
  if (varioind != nullptr) delete varioind;
  if (error) delete vario;
  return vario;
}

/****************************************************************************/
/*!
 **  Find the optimal Truncation Scheme from Variopgs score
 **
 ** \return  The newly created Rule structure
 **
 ** \param[in]  db           Db structure
 ** \param[in]  varioparam   VarioParam structure for the GRFs
 ** \param[in]  ruleprop     RuleProp structure
 ** \param[in]  ngrfmax      Maximum number of underlying GRFs (1 or 2)
 ** \param[in]  verbose      Verbose flag
 **
 *****************************************************************************/
Rule* _rule_auto(Db *db,
                 const VarioParam *varioparam,
                 const RuleProp *ruleprop,
                 int ngrfmax,
                 int verbose)
{
  if (ruleprop == nullptr)
  {
    messerr("RuleProp must be defined");
    return nullptr;
  }
  int flag_stat = ruleprop->isFlagStat();
  const VectorDouble &propcst = ruleprop->getPropCst();
  const Db *dbprop = ruleprop->getDbprop();

  int nscore, r_opt;
  int *rules, flag_rho, flag_correl, opt_correl;
  Vario *vario = nullptr;
  Vario *varioind = nullptr;
  Local_Pgs local_pgs;
  VectorInt facies;
  VectorInt fcmp;
  VectorInt fgrf;
  VectorDouble scores;

  /* Initializations */

  int error = 1;
  Rule *rule = nullptr;
  Relem *Pile_Relem = (Relem*) NULL;
  PropDef *propdef = nullptr;
  String calcName("covnc");

  NCOLOR = db->getFaciesNumber();
  NGRF = ngrfmax;
  NRULE = 2 * NCOLOR - 1;
  BASE = 2 * NGRF;
  flag_rho = 0;
  flag_correl = 0;
  opt_correl = 0;

  /* Core allocation */

  facies.resize(NCOLOR);
  for (int i = 0; i < NCOLOR; i++)
    facies[i] = i + 1;

  /* Preliminary tasks (as in variogram.pgs) */

  if (flag_stat)
  {
    // Calculate the variogram of Indicators
    varioind = new Vario(varioparam, db);
    if (varioind->computeIndic(calcName)) goto label_end;
  }

  if (st_check_test_discret(ERule::STD, 0)) goto label_end;
  st_manage_pgs(0, &local_pgs);

  vario = new Vario(varioparam, db);
  vario->setCalculName(calcName);
  vario->setNVar(NGRF);
  vario->internalVariableResize();
  vario->internalDirectionResize();

  if (st_vario_pgs_check(0, 0, flag_stat, db, NULL, vario, varioind, NULL))
    goto label_end;

  propdef = proportion_manage(1, 1, flag_stat, NGRF, 0, NCOLOR, 0, db, dbprop,
                              propcst, propdef);
  if (propdef == nullptr) goto label_end;
  proportion_rule_process(propdef, EProcessOper::COPY);

  /* Pre-calculation of integrals: Define the structure */

  if (TEST_DISCRET)
    CTABLES = ct_tables_manage(1, 0, 1, 200, 100, -1., 1., NULL);

  /* Allocation */

  st_manage_pgs(1, &local_pgs, db, nullptr, vario, varioind, nullptr, propdef,
                flag_stat, 1, 0, NGRF, NCOLOR, vario->getCalculType());

  if (flag_stat)
  {
    st_define_corpgs(0, 0, 0, &local_pgs);
    st_define_trace(0, 0, &local_pgs);
  }
  else
  {
    st_define_corpgs(opt_correl, flag_rho, 0., &local_pgs);
    st_define_trace(flag_rho, flag_correl, &local_pgs);

    // Prepare the geometry 

    for (int idir = 0; idir < vario->getDirectionNumber(); idir++)
    {
      local_pgs.idircur = idir;
      if (st_variogram_geometry_pgs_calcul(&local_pgs, vario, idir))
        goto label_end;
      st_variogram_geometry_pgs_correct(&local_pgs, vario, idir);
    }
    if (st_variogram_geometry_pgs_final(&local_pgs)) goto label_end;

    // The thresholds are added lately in order to allow calculation of 
    // geometry (without checking the threshold interval (not defined yet)
    if (st_vario_pgs_variable(1, NGRF, NCOLOR, 1, 0, db, propdef, NULL))
      goto label_end;
  }

  /* Elaborate the whole tree of possible Lithotype Rules */

  if (verbose)
    mestitle(1, "Construction of the Tree of candidate Lithotype Rules:");
  Pile_Relem = st_relem_alloc(NULL);
  st_relem_define(Pile_Relem, NCOLOR, facies, ITEST, NULL);
  st_relem_subdivide(Pile_Relem, 1, 1);
  st_relem_explore(Pile_Relem, verbose && debug_query("converge"));

  // Evaluate all possibilities

  fcmp.resize(NCOLOR);
  fgrf.resize(1 + NGRF);
  if (verbose)
  {
    mestitle(1, "List of Rules and corresponding scores:");
    if (flag_stat)
      message("Stationary case");
    else
      message("Non-stationary case");
    if (TEST_DISCRET)
      message(" (Discrete Integration)\n");
    else
      message("\n");
  }
  scores = st_relem_evaluate(Pile_Relem, verbose, fgrf, fcmp, &local_pgs,
                             &nscore, &r_opt);

  /* Get the resulting optimal Rule */

  if (verbose)
  {
    mestitle(1, "Optimal Lithotype Rule:");
    st_rule_print(r_opt, NRULE, Pile_Relem->Rrules, Pile_Relem->Rfipos, false,
                  -1, -1, TEST);
  }
  rules = Pile_Relem->Rrules;
  rule = st_rule_encode(&RULES(r_opt, 0));

  /* Clean the geometry (non-stationary case) */

  if (!flag_stat)
    local_pgs.vorder = vario_order_manage(0, 0, 0, local_pgs.vorder);

  /* Set the error return code */

  error = 0;

  label_end: Pile_Relem = st_relem_free(Pile_Relem);
  if (TEST_DISCRET)
    CTABLES = ct_tables_manage(-1, 0, 1, 200, 100, -1., 1., CTABLES);
  st_manage_pgs(-1, &local_pgs, db, nullptr, vario, varioind, nullptr, propdef,
                flag_stat, 1, 0, NGRF, NCOLOR, vario->getCalculType());
  (void) st_vario_pgs_variable(-1, NGRF, NCOLOR, 1, 0, db, propdef, NULL);

  propdef = proportion_manage(-1, 1, flag_stat, NGRF, 0, NCOLOR, 0, db, dbprop,
                              propcst, propdef);
  if (varioind != nullptr) delete varioind;
  if (vario != nullptr) delete vario;
  if (error) rule = rule_free(rule);
  return (rule);
}
