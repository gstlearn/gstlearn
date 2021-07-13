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
#include "Basic/Utilities.hpp"
#include "Stats/Classical.hpp"

/*! \cond */
typedef struct {
  int    opt_correl;
  int    npar;
  int    flag_rho;
  double rho;
  double params[4];
  double modif[16];
} Local_CorPgs;

typedef struct {
  int     flag_trace;
  int     idir;
  int     ipas;
  int     nrow;
  int     ncol;
  VectorDouble trace;
} Local_TracePgs;

typedef struct {
  Db     *db;
  Rule   *rule;
  Props  *propdef;
  int     flag_stat;
  int     flag_facies;
  int     covtype;
  int     igrfcur;
  int     idircur;
  int     ipascur;
  int     ngrf;
  int     npair;
  int     nfacies;
  int     nfac2;
  int     ifirst;
  int     ilast;
  VectorDouble   d0;
  VectorDouble   d1;
  VectorDouble   memint;
  VectorDouble   stat_proba;
  VectorDouble   stat_thresh;
  Local_CorPgs   corpgs;
  Local_TracePgs tracepgs;
  Model         *model;
  Vario         *vario;
  Vario         *varioind;
  Vario_Order   *vorder; 
} Local_Pgs;

#define VARS(ivar,jvar) (vario->vars[(ivar) * vario->getNVar() + (jvar)])
#define COVS(ivar,jvar) (covs[(ivar) * nfacies + (jvar)])
#define MEMINT(ipair)   (local_pgs->memint[ipair])
#define STAT_PROBA(i,j) (M_R(local_pgs->stat_proba,local_pgs->nfacies,i,j))
#define STAT_THRESH(ifac,igrf,rank) (local_pgs->stat_thresh[2*(nfacies * (igrf) + (ifac))+(rank)])
#define LAG_USED(ipas)  (dir.getSw(dir.getNPas() + ipas + 1) > 0 && \
                         dir.getUtilize(dir.getNPas() + ipas + 1))
#define TABOUT(i,j)      tabout[(j)*neq+(i)]
#define EIGVEC(i,j)      eigvec[(i)*neq+(j)]
#define RULES(ir,i)     (rules[(ir)  * NRULE  + (i)])
#define RULES1(ir,i)    (rules1[(ir) * NRULE  + (i)])
#define RULES2(ir,i)    (rules2[(ir) * NRULE  + (i)])
#define DIVS(is,i)      (divs[(is)   * ncur   + (i)])
#define FIPOS(ir,i)     (fipos[(ir)  * NCOLOR + (i)])
#define FIPOS1(ir,i)    (fipos1[(ir) * NCOLOR + (i)])
#define FIPOS2(ir,i)    (fipos2[(ir) * NCOLOR + (i)])

#define QUANT_DIR 10000
#define F(i,j) (st_index(i,j))

static double EPS = 1.e-05;
static double GS_TOLSTOP     = 5.e-02;
static double GS_TOLSTOP_RHO = 1.e-01;
static CTables *CTABLES      = NULL;
static int TEST_DISCRET      = 0;
static int TEST_TIME         = 0;
static int NCOLOR            = 0;
static int NGRF              = 0;
static int NRULE             = 0;
static int BASE              = 0;
static double clock_start    = 0;

// Needed declarations due to intricated recursions
static Relem *st_relem_free(Relem *relem);
static void   st_relem_explore(Relem *relem,int verbose);
/*! \endcond */

/****************************************************************************
 **
 ** FUNCTION: st_relem_alloc
 **
 ** PURPOSE:  Allocate a new Relem structure
 **
 ** IN_ARGS:  old_split  : Pointer to the calling Split
 **
 *****************************************************************************/
static Relem *st_relem_alloc(Split *old_split)

{
  Relem *relem;

  /* Initializations of the New Relem structure */

  relem = (Relem *) mem_alloc(sizeof(Relem),1);
  relem->nfacies   = 0;
  relem->nsplit    = 0;
  relem->nrule     = 0;
  relem->nbyrule   = 0;
  relem->facies    = (int *) NULL;
  relem->rules     = (int *) NULL;
  relem->fipos     = (int *) NULL;
  relem->old_split = old_split;
  return(relem);
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
static Split *st_split_alloc(Relem *old_relem)

{
  Split *split;

  /* Initializations of the New Split structure */

  split = (Split *) mem_alloc(sizeof(Split),1);
  split->oper      = 0;
  split->nrule     = 0;
  split->nbyrule   = 0;
  split->old_relem = old_relem;
  split->rules     = (int *) NULL;
  split->fipos     = (int *) NULL;
  for (int i=0; i<2; i++) split->relems[i] = (Relem *) NULL;
  return(split);
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
static int st_define_fipos(int oper,
                           int side)
{
  int reponse;

  if (IFFFF(side))
    reponse = 1;
  else
    reponse = 2*(oper-1) + (1-side) + 1;

  return(reponse);
}

/****************************************************************************
 **
 ** FUNCTION: st_relem_define
 **
 ** PURPOSE:  Define the list of facies in the current Relem structure
 ** 
 *****************************************************************************/
static void st_relem_define(Relem *relem,
                            int    nfacies,
                            int   *facies,
                            int    side,
                            int   *poss)
{
  int ecr,number;

  if (relem == (Relem *) NULL) return;

  if (poss == (int *) NULL) 
    number = nfacies;
  else
  {
    number = 0;
    for (int i=0; i<nfacies; i++)
      if (poss[i] == side) number++;
  }

  relem->nfacies = number;
  relem->facies  = (int *) mem_alloc(sizeof(int) * number,1);
  relem->fipos   = (int *) mem_alloc(sizeof(int) * NCOLOR,1);
  for (int i=0; i<NCOLOR; i++) relem->fipos[i] = 0;

  ecr = 0;
  for (int i=0; i<nfacies; i++) 
  {
    if (poss == (int *) NULL || poss[i] == side)
      relem->facies[ecr++] = facies[i];
  }

  if (number == 1)
    relem->fipos[relem->facies[0]-1] = 1;
}

/****************************************************************************
 **
 ** FUNCTION: st_rule_print
 **
 *****************************************************************************/
static void st_rule_print(int    rank,
                          int    nbyrule,
                          int   *rules,
                          int   *fipos,
                          int    flag_rank,
                          int    flag_similar,
                          int    flag_igrf,
                          double score)
{
  int value,iscore,loc0,loc1;

  // Print the Rank (optional)

  if (flag_rank)
    message("%4d:",rank+1);

  // Print the Rule

  for (int ic=0; ic<nbyrule; ic++)
  {
    value = RULES(rank,ic);
    if (value == 1001)
      message("  S");
    else if (value == 1002)
      message("  T");
    else
      message(" %2d",value);
  }

  // Print the score (if available)

  if (! FFFF(score)) 
  {
    iscore = (int) score;
    message(" -> %d",iscore);
  }

  // Print the Facies rank

  message(" (");
  for (int ic=0; ic<NCOLOR; ic++) 
    message(" %3d",FIPOS(rank,ic));
  message(" )");

  // Print the similar score

  if (flag_similar >= 0) 
  {
    message(" [Same as: %3d (",flag_similar+1);
    
    loc0 = flag_igrf;
    for (int igrf=0; igrf<NGRF; igrf++)
    {
      loc1 = loc0 / 2;
      if (loc0 - 2 * loc1 > 0) message(" G%d",igrf+1);
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
                           int  nrule,
                           int  nbyrule,
                           int *rules,
                           int *fipos)
{
  if (nrule <= 0) return;
  message("%s (Nrule=%d, Nbyrule=%d):\n",title,nrule,nbyrule);
  for (int ir=0; ir<nrule; ir++)
    st_rule_print(ir,nbyrule,rules,fipos,0,-1,-1,TEST);
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
static void st_relem_subdivide(Relem *relem0,
                               int    half,
                               int    noper)
{
  Split *split;
  int   *divs,ndiv,number,ncur,previous_oper,verbose;

  verbose = 0;
  ncur    = relem0->nfacies;
  if (ncur <= 1) return;

  previous_oper = 1;
  if (relem0->old_split != NULL)
    previous_oper = relem0->old_split->oper;

  divs   = ut_split_into_two(ncur,half,verbose,&ndiv);
  number = ndiv * noper;
  divs = (int *) mem_free((char *) divs);
  
  relem0->splits = (Local_Split **) mem_alloc(sizeof(Local_Split *) * number,1);
  
  number = 0;
  for (int oper=1; oper<=noper; oper++)
  {
    half = (oper == previous_oper);
    divs = ut_split_into_two(ncur,half,verbose,&ndiv);
    for (int is=0; is<ndiv; is++,number++)
    {
      relem0->splits[number] = split = st_split_alloc(relem0);
      split->oper = oper;
      for (int i=0; i<2; i++)
      {
        split->relems[i] = st_relem_alloc(split);
        st_relem_define(split->relems[i],ncur,relem0->facies,
                        1-i,&DIVS(is,0));
        st_relem_subdivide(split->relems[i],0,NGRF);
      }
    }
    divs = (int *) mem_free((char *) divs);
  }

  relem0->splits = (Local_Split **) mem_alloc(sizeof(Local_Split *) * number,1);
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
static Split *st_split_free(Split *split)
{
  if (split == (Split *) NULL) return(split);

  /* Free the descending substructures (relem) */

  for (int i=0; i<2; i++)
    split->relems[i] = st_relem_free(split->relems[i]);

  /* Free the local arrays */

  split->rules = (int *) mem_free((char *) split->rules);
  split->fipos = (int *) mem_free((char *) split->fipos);

  /* Free the Split structure itself */

  split = (Split *) mem_free((char *) split);

  return(split);
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
static Relem *st_relem_free(Relem *relem)
{
  if (relem == NULL) return(relem);

  /* Free the descending substructures (Split) */

  for (int is=0; is<relem->nsplit; is++)
    relem->splits[is] = st_split_free(relem->splits[is]);

  /* Free the local arrays */

  relem->facies = (int *) mem_free((char *) relem->facies);
  relem->rules  = (int *) mem_free((char *) relem->rules);
  relem->fipos  = (int *) mem_free((char *) relem->fipos);

  /* Free the Relem structure itself */

  relem = (Relem *) mem_free((char *) relem);

  return(relem);
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
static void st_variogram_define_vars(Vario *vario,
                                     Rule  *rule,
                                     int    ngrf)
{
  int igrf,jgrf;
  
  for (igrf=0; igrf<ngrf; igrf++)
    for (jgrf=0; jgrf<ngrf; jgrf++)
    {
      if (igrf == jgrf)
        vario->setVars(igrf,jgrf,1.);
      else
        vario->setVars(igrf,jgrf,rule->getRho());
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
static void st_set_rho(double rho,
                       Local_Pgs *local_pgs)
{
  double rho2;
  int iech,ifac,ngrf;
  Db *db = local_pgs->db;
  Props *propdef = local_pgs->propdef;
  Rule * rule    = local_pgs->rule;
  int flag_stat  = local_pgs->flag_stat;
  double t1min, t1max, t2min, t2max;

  local_pgs->corpgs.rho = rho;
  local_pgs->rule->setRho(rho);
  rho2 = rho * rho;
  st_variogram_define_vars(local_pgs->vario,local_pgs->rule,local_pgs->ngrf);

  /* Define the Thresholds */
  
  if (flag_stat)
  {
    rule->setProportions(propdef->proploc);
  }
  else
  {
    ngrf = local_pgs->ngrf;
    for (iech=0; iech<get_NECH(db); iech++)
    {
      if (! db->isActive(iech)) continue;
      ifac = (int) db->getVariable(iech,0);
      if (rule_thresh_define(propdef,db,rule,ifac,iech,0,0,0,
                             &t1min,&t1max,&t2min,&t2max)) return;
      db->setLowerBound(iech,0,t1min);
      db->setUpperBound(iech,0,t1max);
      if (ngrf > 1)
      {
        db->setLowerBound(iech,1,t2min);
        db->setUpperBound(iech,1,t2max);
      }
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
**  Calculate the thresholds in the stationary case
**
** \returns Error return code
**
** \param[in]  local_pgs     Local_Pgs structure
**
*****************************************************************************/
static int st_calculate_thresh_stat(Local_Pgs *local_pgs)

{
  int nfacies,ifac,ngrf;
  double t1min,t1max,t2min,t2max;

  nfacies = local_pgs->nfacies;
  ngrf    = local_pgs->ngrf;

  for (ifac=0; ifac<nfacies; ifac++)
  {
    if (rule_thresh_define(local_pgs->propdef,local_pgs->db,local_pgs->rule,
                           ifac+1,0,0,0,0,&t1min,&t1max,&t2min,&t2max)) return(1);
    if (! TEST_DISCRET)
    {
      STAT_THRESH(ifac,0,0) = t1min;
      STAT_THRESH(ifac,0,1) = t1max;
      if (ngrf > 1) 
      {
        STAT_THRESH(ifac,1,0) = t2min;
        STAT_THRESH(ifac,1,1) = t2max;
      }
    }
    else
    {
      STAT_THRESH(ifac,0,0) = (double) 
        ct_tableone_getrank_from_proba(CTABLES,t1min);
      STAT_THRESH(ifac,0,1) = (double) 
        ct_tableone_getrank_from_proba(CTABLES,t1max);
      if (ngrf > 1)
      {
        STAT_THRESH(ifac,1,0) = (double)
          ct_tableone_getrank_from_proba(CTABLES,t2min);
        STAT_THRESH(ifac,1,1) = (double)
          ct_tableone_getrank_from_proba(CTABLES,t2max);
      }
    }
  }
  return(0);
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
** \param[in]  propdef       Props structure
** \param[in]  rule          Lithotype Rule definition
**
** \param[out] iptr_p        Pointer to the proportions
** \param[out] iptr_l        Pointer to the lower bound(s)
** \param[out] iptr_u        Pointer to the upper bound(s)
** \param[out] iptr_rl       Pointer to the discretized lower bound indice(s)
** \param[out] iptr_ru       Pointer to the discretized upper bound indice(s)
**
*****************************************************************************/
static int st_vario_pgs_variable(int    mode,
                                 int    ngrf, 
                                 int    nfacies,
                                 int    flag_one,
                                 int    flag_prop,
                                 Db    *db,
                                 Props *propdef,
                                 Rule  *rule,
                                 int   *iptr_p,
                                 int   *iptr_l,
                                 int   *iptr_u,
                                 int   *iptr_rl,
                                 int   *iptr_ru)
{
  int    number,ifac,jfac,nloop;
  double t1min,t1max,t2min,t2max;

  // Dispatch

  number = (flag_one) ? ngrf : ngrf * nfacies;
  switch (mode)
  {
    case 1:
      
      /* Allocation */
      
      (*iptr_p) = (*iptr_l) = (*iptr_u) = (*iptr_rl) = (*iptr_ru) = -1;
      
      /* The proportions */
      
      if (flag_prop && db->getProportionNumber() == 0)
      {
        (*iptr_p) = db->addFields(nfacies,0.);
        if ((*iptr_p) < 0) return(1);
        db->setLocatorsByAttribute(nfacies,*iptr_p,LOC_P);
      }
      
      /* The bounds */
      
      if (! TEST_DISCRET)
      {
        (*iptr_l) = db->addFields(number,0.);
        if ((*iptr_l) < 0) return(1);
        db->setLocatorsByAttribute(number,*iptr_l,LOC_L);
        
        (*iptr_u) = db->addFields(number,0.);
        if ((*iptr_u) < 0) return(1);
        db->setLocatorsByAttribute(number,*iptr_u,LOC_U);
      }
      else
      {
        (*iptr_rl) = db->addFields(number,0.);
        if ((*iptr_rl) < 0) return(1);
        db->setLocatorsByAttribute(number,*iptr_rl,LOC_RKLOW);
        
        (*iptr_ru) = db->addFields(number,0.);
        if ((*iptr_ru) < 0) return(1);
        db->setLocatorsByAttribute(number,*iptr_ru,LOC_RKUP);
      }
      break;

    case 0:

      /* Evaluate the bounds */
      /* Use dummy rho value in order to avoid discarding pairs in geometry */
    
      nloop = (flag_one) ? 1 : nfacies;
      for (int iech=0; iech<get_NECH(db); iech++)
      {
        if (! db->isActive(iech)) continue;
        
        for (int i=0; i<nloop; i++)
        {
          ifac = (flag_one) ? (int) db->getVariable(iech,0) : i;
          
          jfac = (flag_one) ? ifac : ifac+1;
          if (rule_thresh_define(propdef,db,rule,jfac,iech,0,0,0,
                                 &t1min,&t1max,&t2min,&t2max)) return(1);
          
          /* Define the proportions */
          
          if ((*iptr_p) >= 0) db->setProportion(iech,ifac,propdef->propmem[ifac]);
          
          /* Define the bounds */
          
          if (! TEST_DISCRET)
          {
            jfac = (flag_one) ? 0 : ifac;
            db->setLowerBound(iech,jfac,t1min);
            db->setUpperBound(iech,jfac,t1max);
            if (ngrf > 1)
            {
              jfac = (flag_one) ? 1 : nfacies + ifac;
              db->setLowerBound(iech,jfac,t2min);
              db->setUpperBound(iech,jfac,t2max);
            }
          }
          else
          {
            jfac = (flag_one) ? 0 : ifac;
            db->setLowerInterval(iech,jfac,
                      (double) ct_tableone_getrank_from_proba(CTABLES,t1min));
            db->setUpperInterval(iech,jfac,
                      (double) ct_tableone_getrank_from_proba(CTABLES,t1max));
            if (ngrf > 1)
            {
              jfac = (flag_one) ? 1 : nfacies + ifac;
              db->setLowerInterval(iech,jfac,
                        (double) ct_tableone_getrank_from_proba(CTABLES,t2min));
              db->setUpperInterval(iech,jfac,
                        (double) ct_tableone_getrank_from_proba(CTABLES,t2max));
            }
          }
        }
      }
      break;

    case -1:
      
      /* Deallocation */
      
      if ((*iptr_p) >= 0) 
        (void) db_attribute_del_mult(db,*iptr_p,nfacies);
      if ((*iptr_l) >= 0) 
        (void) db_attribute_del_mult(db,*iptr_l,number);
      if ((*iptr_u) >= 0) 
        (void) db_attribute_del_mult(db,*iptr_u,number);
      if ((*iptr_rl) >= 0) 
        (void) db_attribute_del_mult(db,*iptr_rl,number);
      if ((*iptr_ru) >= 0) 
        (void) db_attribute_del_mult(db,*iptr_ru,number);
      break;
  }
  return(0);
}

/****************************************************************************
 **
 ** FUNCTION: st_rule_encode
 **
 ** PURPOSE: Encode the set of numerical rule items into a rule sentence
 ** PURPOSE: that can be used to create a Rule structure
 **
 *****************************************************************************/
static Rule *st_rule_encode(int *string)
{
  Rule *rule;

  rule   = (Rule *) NULL;
  VectorInt n_type = VectorInt(NRULE);
  VectorInt n_facs = VectorInt(NRULE);

  for (int i=0; i<NRULE; i++)
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

  rule = new Rule(n_type,n_facs);
  return(rule);
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
  Local_TracePgs *tracepgs;
  double total[2],totsum;
  int    nrow,ncol;

  tracepgs = &local_pgs->tracepgs;
  nrow = tracepgs->nrow;
  ncol = tracepgs->ncol;
  if (nrow <= 0 || ncol <= 0) return TEST;

  /* Evaluate the sum of the score */

  total[0] = total[1] = 0.;
  for (int irow=0; irow<nrow; irow++)
  {
    total[0] += tracepgs->trace[ncol * irow + 2];
    if (ncol >= 5) total[1] += tracepgs->trace[ncol * irow + 4];
  }
  totsum = total[0] + total[1];
  set_keypair("vario.pgs_score",1,1,1,&totsum);

  if (tracepgs->flag_trace)
    set_keypair("vario.pgs_stack",1,nrow,ncol,tracepgs->trace.data());
  
  return(totsum);
}

/****************************************************************************/
/*!
**  Patch the central value (dist=0) of the covariances
**
** \param[in]  local_pgs   Local_Pgs structure
** \param[in]  vario       Vario structure for the GRFs to be filled
** \param[in]  idir        Rank of the direction
** \param[in]  ngrf        Number of GRFs
** \param[in]  rho         Correlation coefficient
**
*****************************************************************************/
static void st_variogram_patch_C00(Local_Pgs *local_pgs,
                                   Vario *vario,
                                   int    idir,
                                   int    ngrf,
                                   double rho)
{
  Db* db  = local_pgs->db;
  int nech = 0.;
  if (db != nullptr) nech = db->getActiveSampleNumber();

  const Dir& dir = vario->getDirs(idir);
  for (int igrf=0; igrf<ngrf; igrf++)
    for (int jgrf=0; jgrf<=igrf; jgrf++)
    {      
      // Get the central address
      int iad = dir.getAddress(igrf,jgrf,0,false,0);
      vario->setSw(idir,iad,nech);
      vario->setHh(idir,iad,0.);
      if (igrf == jgrf)
        vario->setGg(idir,iad,1.);
      else
        vario->setGg(idir,iad,rho);
    }
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
  int nrow,ncol,iad;

  /* Initializations */

  tracepgs = &local_pgs->tracepgs;
  if (! tracepgs->flag_trace) return;
  ncol = tracepgs->ncol;
  nrow = tracepgs->nrow;
  iad  = ncol * nrow;

  nrow++;
  tracepgs->trace.resize(nrow * ncol);
  for (int icol=0; icol<ncol; icol++) tracepgs->trace[iad+icol] = -1.11;
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
static double st_func_search_stat(double  correl,
                                  void *user_data)
{
  double  low[2],up[2],err,cround;
  int     infin[2],ier;
  Local_Pgs *local_pgs;

  /* Initializations */

  local_pgs         = (Local_Pgs *) user_data;
  int iconf0        = 0;
  double releps     = 0.;
  double abseps     = EPS;
  double maxpts     = 8000;

  int nfacies    = local_pgs->nfacies;
  int ipas       = local_pgs->ipascur;
  int idir       = local_pgs->idircur;
  int igrf       = local_pgs->igrfcur;
  Dir dir        = local_pgs->varioind->getDirs(idir);

  if (TEST_DISCRET)
    iconf0 = ct_tableone_covrank(CTABLES,correl,&cround);
  
  double sum = 0.;
  for (int ifac1=0; ifac1<nfacies; ifac1++)
    for (int ifac2=0; ifac2<nfacies; ifac2++)
    {
      low[0]  = STAT_THRESH(ifac1,igrf,0);
      up[0]   = STAT_THRESH(ifac1,igrf,1);
      low[1]  = STAT_THRESH(ifac2,igrf,0);
      up[1]   = STAT_THRESH(ifac2,igrf,1);
      double proba;
      if (! TEST_DISCRET)
      {
        if (correl == 0.)
        {
          double p1min = law_cdf_gaussian(low[0]);
          double p1max = law_cdf_gaussian( up[0]);
          double p2min = law_cdf_gaussian(low[1]);
          double p2max = law_cdf_gaussian( up[1]);
          proba = (p1max - p1min) * (p2max - p2min);
        }
        else
        {
          infin[0] = mvndst_infin(low[0],up[0]);
          infin[1] = mvndst_infin(low[1],up[1]);
          mvndst(2,low,up,infin,&correl,maxpts,abseps,releps,&err,&proba,&ier);
        }
      }
      else
      {
        proba = ct_tableone_calculate_by_rank(CTABLES,iconf0,low,up);
      }

      double logp = (proba <= 0.) ? -1.e30 : log(proba);
      int iad     = dir.getAddress(ifac1,ifac2,ipas,false,1);
      double sw   = dir.getSw(iad);
      double gg   = dir.getGg(iad);
      iad         = dir.getAddress(ifac1,ifac2,ipas,false,-1);
      gg         += dir.getGg(iad);
      sum -= logp * gg * sw / 2.;
    }
  return(0.5 * sum);
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
static double st_func_search_nostat(double  correl,
                                    void *user_data)
{
  double  low[2],up[2],rlow[2],rup[2],abseps,releps,err,proba,sum;
  double  p1min,p1max,p2min,p2max,dist,w1,w2,logp,cround;
  int     ipair,i,i1,i2,infin[2],maxpts,ier,ifac1,ifac2,iconf0;
  Local_Pgs *local_pgs;

  /* Initializations */

  local_pgs  = (Local_Pgs *) user_data;
  iconf0     = 0;
  releps     = 0.;
  abseps     = EPS;
  maxpts     = 8000;

  if (TEST_DISCRET)
    iconf0 = ct_tableone_covrank(CTABLES,correl,&cround); 

  /* Reset the pre-calculation array (only if flag_stat) */
  if (local_pgs->flag_stat)
    for (i=0; i<local_pgs->nfac2; i++) local_pgs->stat_proba[i] = TEST;

  sum = 0.;
  for (ipair=local_pgs->ifirst; ipair<local_pgs->ilast; ipair++)
  {
    vario_order_get_indices(local_pgs->vorder,ipair,&i1,&i2,&dist);
    w1 = local_pgs->db->getWeight(i1);
    w2 = local_pgs->db->getWeight(i2);
    ifac1 = ifac2 = -1;
    proba = TEST;
    if (local_pgs->flag_stat)
    {

      /* In the stationary case, search in the lookup table first */
 
      ifac1 = (int) local_pgs->db->getVariable(i1,0) - 1;
      if (ifac1 < 0 || ifac1 >= local_pgs->nfacies) continue;
      ifac2 = (int) local_pgs->db->getVariable(i2,0) - 1;
      if (ifac2 < 0 || ifac2 >= local_pgs->nfacies) continue;
      proba = STAT_PROBA(ifac1,ifac2);
    }

    if (FFFF(proba))
    {
      if (! TEST_DISCRET)
      {
        low[0] = local_pgs->db->getLowerBound(i1,local_pgs->igrfcur);
        up[0]  = local_pgs->db->getUpperBound(i1,local_pgs->igrfcur);
        low[1] = local_pgs->db->getLowerBound(i2,local_pgs->igrfcur);
        up[1]  = local_pgs->db->getUpperBound(i2,local_pgs->igrfcur);
        
        if (correl == 0.)
        {
          p1min = law_cdf_gaussian(low[0]);
          p1max = law_cdf_gaussian( up[0]);
          p2min = law_cdf_gaussian(low[1]);
          p2max = law_cdf_gaussian( up[1]);
          proba = (p1max - p1min) * (p2max - p2min);
        }
        else
        {
          infin[0] = mvndst_infin(low[0],up[0]);
          infin[1] = mvndst_infin(low[1],up[1]);
          mvndst(2,low,up,infin,&correl,maxpts,abseps,releps,&err,&proba,&ier);
        }
      }
      else
      {
        rlow[0] = local_pgs->db->getLowerInterval(i1,local_pgs->igrfcur);
        rup[0]  = local_pgs->db->getUpperInterval(i1,local_pgs->igrfcur);
        rlow[1] = local_pgs->db->getLowerInterval(i2,local_pgs->igrfcur);
        rup[1]  = local_pgs->db->getUpperInterval(i2,local_pgs->igrfcur);
        proba = ct_tableone_calculate_by_rank(CTABLES,iconf0,rlow,rup);
      }
      if (local_pgs->flag_stat) 
        STAT_PROBA(ifac1,ifac2) = STAT_PROBA(ifac2,ifac1) = proba;
    }
    logp = (proba <= 0.) ? -1.e30 : log(proba);
    sum -= w1 * w2 * logp;
  }
  return(0.5 * sum);
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
                         double     value0,
                         double     value1,
                         int        origin,
                         int        number,
                         double    *values)
{
  Local_TracePgs *tracepgs;
  int i,ncol,nrow,iad;

  /* Initializations */

  tracepgs = &local_pgs->tracepgs;
  if (! tracepgs->flag_trace) return;
  nrow  = tracepgs->nrow;
  ncol  = tracepgs->ncol;
  iad   = ncol * (nrow-1);
  if (2+origin+number > ncol) 
    messageAbort("Error in Trace dimension (ncol=%d origin=%d number=%d)",
              ncol,origin,number);

  /* Store the information */

  tracepgs->trace[iad]   = value0;
  tracepgs->trace[iad+1] = value1;

  for (i=0; i<number; i++)
    tracepgs->trace[iad+2+origin+i] = values[i];
}

/****************************************************************************/
/*!
**  Performing the variogram calculations (stationary case)
**
** \return  Error return code
**
** \param[in]  vario         Vario structure for the GRFs to be filled
** \param[in]  rule          Lithotype Rule definition
** \param[in]  propdef       Props structure
** \param[in]  local_pgs     Local_Pgs structure
** \param[in]  ngrf          Number of GRFs
**
*****************************************************************************/
static int st_varcalc_from_vario_stat(Vario *vario,
                                      Rule  *rule,
                                      Props *propdef,
                                      Local_Pgs *local_pgs,
                                      int ngrf)
{
  int  iad;
  double result,testval,varloc,niter;
  
  /* Initializations */

  st_set_rho(0.,local_pgs);
    
  /* Loop on the directions */

  for (int idir=0; idir<vario->getDirectionNumber(); idir++)
  {
    const Dir& dir = vario->getDirs(idir);
    local_pgs->idircur = idir;

    /* Set the value of C(0) */

    st_variogram_patch_C00(local_pgs,vario,idir,ngrf,0.);
    
    /* Loop on the lags */

    for (int ipas=0; ipas<dir.getNPas(); ipas++)
    {
      mes_process("Inverting Variogram Lag",dir.getNPas(),ipas);
      local_pgs->ipascur = ipas;
      trace_add_row(local_pgs);

      /* Loop on the GRFs */

      for (int igrf=0; igrf<ngrf; igrf++)
      { 
        local_pgs->igrfcur = igrf; 
        result = golden_search(st_func_search_stat,(void *) local_pgs,
                               GS_TOLSTOP,-1.,1.,&testval,&niter);
        trace_define(local_pgs,idir+1,ipas+1,2*igrf  ,1,&testval);
        trace_define(local_pgs,idir+1,ipas+1,2*igrf+1,1,&niter);

        for (int jgrf=0; jgrf<=igrf; jgrf++)
        {
          varloc = (igrf == jgrf) ? result : 0.;
          iad = dir.getAddress(igrf,jgrf,ipas,false,1);
          vario->setGg(idir,iad,varloc);
          iad = dir.getAddress(igrf,jgrf,ipas,false,-1);
          vario->setGg(idir,iad,varloc);
          
          if (debug_query("converge"))
            message("Lag:%d - Grf:%d - Variogram(%d) = %lf\n",
                    ipas,igrf,iad,varloc);
        }
      }
    }
  }
  return(0);
}
  
/****************************************************************************/
/*!
**  Define the trace
**
** \param[in]  flag_rho    1 if rho has to be calculated, 0 otherwise
** \param[in]  flag_correl 1 for the correlated case; 0 otherwise
** \param[in]  local_pgs   Local_Pgs structure
**
** \param[out]  local_pgs   Local_Pgs structure
**
*****************************************************************************/
static void st_define_trace(int        flag_rho,
                            int        flag_correl,
                            Local_Pgs *local_pgs)
{
  Local_TracePgs *tracepgs;

  tracepgs = &local_pgs->tracepgs;
  tracepgs->flag_trace = ! flag_rho;
  if (! tracepgs->flag_trace) return;

  if (! flag_correl)
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
** \param[in]  local_pgs  Local_TracePgs structure
**
** \param[out]  local_pgs  Local_Pgs structure
**
*****************************************************************************/
static void st_retrace_define(Local_Pgs *local_pgs)

{
  Local_TracePgs *tracepgs;

  tracepgs = &local_pgs->tracepgs;
  if (! tracepgs->flag_trace) return;

  // Clean the previous trace 

  tracepgs->trace.clear();

  // Initialize the new trace

  st_define_trace(0,0,local_pgs);
}

/****************************************************************************/
/*!
**  Evaluate the variogram of one underlying GRF 
**
** \param[in]  local_pgs  Local_Pgs structure
** \param[in]   idir      Rank of the direction
**
*****************************************************************************/
static void st_varcalc_uncorrelated_grf(Local_Pgs *local_pgs,
                                        int        idir)
{
  int    ipas,iad,igrf,jgrf,ngrf;
  double result,testval,niter,varloc;
  Vario * vario;

  vario = local_pgs->vario;
  const Dir& dir = vario->getDirs(idir);
  ngrf  = local_pgs->ngrf;

  /* Loop on the lags */

  for (ipas=0; ipas<dir.getNPas(); ipas++)
  {
    mes_process("Inverting Variogram Lag",dir.getNPas(),ipas);
    local_pgs->ipascur = ipas;
    trace_add_row(local_pgs);
    if (! LAG_USED(ipas)) continue;
    vario_order_get_bounds(local_pgs->vorder,idir,ipas,
                           &local_pgs->ifirst,&local_pgs->ilast);
    if (local_pgs->ifirst >= local_pgs->ilast) continue;

    for (igrf=0; igrf<ngrf; igrf++)
    { 
      local_pgs->igrfcur = igrf; 
      result = golden_search(st_func_search_nostat,(void *) local_pgs,
                             GS_TOLSTOP,-1.,1.,&testval,&niter);
      trace_define(local_pgs,idir+1,ipas+1,2*igrf,  1,&testval);
      trace_define(local_pgs,idir+1,ipas+1,2*igrf+1,1,&niter);
      for (jgrf=0; jgrf<=igrf; jgrf++)
      {
        varloc = (igrf == jgrf) ? result : 0.;
        iad = dir.getAddress(igrf,jgrf,ipas,false,1);
        vario->setGg(idir,iad,varloc);
        iad = dir.getAddress(igrf,jgrf,ipas,false,-1);
        vario->setGg(idir,iad,varloc);
        
        if (debug_query("converge"))
          message("Lag:%d - Grf:%d - Variogram(%d) = %lf\n",
                  ipas,igrf,iad,varloc);
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
static double st_rule_calcul(Local_Pgs *local_pgs,
                             int       *string)
{
  double score;
  int iptr_p,iptr_l,iptr_u,iptr_rl,iptr_ru;

  /* Preliminary assignments */

  score = 0.;
  local_pgs->rule = st_rule_encode(string);
  local_pgs->ngrf = local_pgs->rule->getGRFNumber();
  st_retrace_define(local_pgs);

  if (local_pgs->flag_stat)
  {
    st_set_rho(0.,local_pgs);
    (void) st_calculate_thresh_stat(local_pgs);
    st_varcalc_from_vario_stat(local_pgs->vario,local_pgs->rule,
                               local_pgs->propdef,local_pgs,
                               local_pgs->ngrf);
  }
  else
  {
    (void) st_vario_pgs_variable(0,local_pgs->ngrf,local_pgs->nfacies,
                                 1,0,local_pgs->db,local_pgs->propdef,
                                 local_pgs->rule,
                                 &iptr_p,&iptr_l,&iptr_u,&iptr_rl,&iptr_ru);
    st_set_rho(0.,local_pgs);
    for (int idir=0; idir<local_pgs->vario->getDirectionNumber(); idir++)
    {
      local_pgs->idircur = idir;
      st_variogram_patch_C00(local_pgs,local_pgs->vario,idir,
                             local_pgs->ngrf,local_pgs->rule->getRho());
      st_varcalc_uncorrelated_grf(local_pgs,idir);
    }
  }

  /* Deallocation of the Rule */

  local_pgs->rule = rule_free(local_pgs->rule);

  score = st_extract_trace(local_pgs);
  return(score);
}

/****************************************************************************
 **
 ** FUNCTION: st_permut
 **
 *****************************************************************************/
static int st_permut(int value,
                     int igrf)
{
  if (igrf == 0)
  {
    if (value == 1) return(2);
    if (value == 2) return(1);
  }
  else if (igrf == 1)
  {
    if (value == 3) return(4);
    if (value == 4) return(3);
  }
  else if (igrf == 2)
  {
    if (value == 5) return(6);
    if (value == 6) return(5);
  }
  else
  {
    messageAbort("Function st_permut has been programmed up to 3 GRFs");
  }
  return(value);
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
static int st_fipos_encode(int *fgrf)
{
  int nmax,found,fipos;

  nmax = 1 + NGRF;

  found = fipos = 0;
  for (int i=nmax-1; i>=0; i--)
  {
    if (fgrf[i] > 0) found++;
    if (found == 1) fipos = 1;
    fipos = fipos * BASE + fgrf[i];
  }
  return(fipos);
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
static void st_fipos_decode(int  fipos,
                            int *fgrf)
{
  int nmax,div;

  nmax = 1 + NGRF;
  for (int i=0; i<nmax; i++) fgrf[i] = 0;
  for (int i=0; i<nmax; i++)
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
static int st_update_orientation(int  fac0,
                                 int  igrf_cas,
                                 int *fgrf)
{
  int fac,facp,nmax,loc0,loc1;

  /* Preliminary check */

  nmax = 1 + NGRF;
  fac  = fac0;
  if (fac0 < 0) return(fac0);

  /* Decomposition */

  st_fipos_decode(fac,fgrf);

  /* Loop on the GRF permutations */

  loc0 = igrf_cas;
  for (int igrf=0; igrf<NGRF; igrf++)
  {
    loc1 = loc0 / 2;
    if (loc0 - 2 * loc1 > 0)
    {
      
      /* Update the orientation of 'grf' */
      
      for (int i=0; i<nmax; i++)
        fgrf[i] = st_permut(fgrf[i],igrf);
    }
    loc0 = loc1;
  }
  
  /* Recomposition */
  
  facp = st_fipos_encode(fgrf);
  
  return(facp);
}

/****************************************************************************
 **
 ** FUNCTION: st_same_score
 **
 ** REMARKS In this function, we change the orientation of one or several GRFs
 **
 *****************************************************************************/
static int st_same_score(Relem *relem,
                         int    ir0,
                         int    igrf_cas,
                         int   *fgrf,
                         int   *fcmp)
{
  int *fipos,flag_same;

  fipos = relem->fipos;
  if (ir0 <= 0) return(-1);

  // Modify the orientation of 'grf' for the current 'fipos'

  for (int ic=0; ic<NCOLOR; ic++)
    fcmp[ic] = st_update_orientation(FIPOS(ir0,ic),igrf_cas,fgrf);

  // Look if the same 'fipos' has already been calculated

  for (int ir=0; ir<ir0; ir++)
  {
    flag_same = 1;
    for (int ic=0; ic<NCOLOR && flag_same; ic++)
    {
      if (FIPOS(ir,ic) != fcmp[ic]) flag_same = 0;
    }
    if (flag_same) return(ir);
  }

  return(-1);
}

/****************************************************************************
 **
 ** FUNCTION: st_relem_evaluate
 **
 *****************************************************************************/
static double *st_relem_evaluate(Relem *relem,
                                 int    verbose,
                                 int   *fgrf,
                                 int   *fcmp,
                                 Local_Pgs *local_pgs,
                                 int   *nscore,
                                 int   *r_opt)
{
  int    *rules,*fipos;
  int     nrule,indice,nmax,flag_check,flag_skip,igrf_cas,number,igrf_opt;
  double *scores,score_ref;

  /* Initializations */

  flag_check = (int) get_keypone("Multi_Score_Check",0.);
  flag_skip  = (int) get_keypone("Multi_Score_Skip_Print",0.);
  nmax       = (int) pow(2.,(double) NGRF);
  nrule      = relem->nrule;
  rules      = relem->rules;
  fipos      = relem->fipos;
  *nscore    = nrule;
  
  /* Core allocation */

  scores = (double *) mem_alloc(sizeof(double) * nrule,1);
  
  *r_opt = number = 0;
  for (int ir=0; ir<nrule; ir++)
  {

    // Check if the same flag has already been found

    indice = igrf_opt = -1;
    for (igrf_cas=1; igrf_cas<nmax && igrf_opt < 0; igrf_cas++)
    {
      indice = st_same_score(relem,ir,igrf_cas,fgrf,fcmp);
      if (indice >= 0) igrf_opt = igrf_cas;
    }

    // Set the score 

    if (indice >= 0)
      scores[ir] = scores[indice];
    else
    {
      number++;
      scores[ir] = st_rule_calcul(local_pgs,&RULES(ir,0));
    }

    // When Multi_Score_Check, calculate the score even if already defined
    // and compare both results

    if (flag_check && indice >= 0)
    {
      score_ref = st_rule_calcul(local_pgs,&RULES(ir,0));
      if (ABS(scores[ir] - score_ref) > 1.e-10 * score_ref)
      {
        messerr("Warning: Difference between score stored and re-evaluated:");
        messerr("- as already stored = %lf",scores[ir]);
        messerr("- as re-evaluated   = %lf",score_ref);
      }
    }

    // Optional printout

    if (verbose) 
    {
      if (! flag_skip || indice < 0)
        st_rule_print(ir,NRULE,rules,fipos,1,indice,igrf_opt,scores[ir]);
    }

    // Ranking the scores

    if (scores[ir] < scores[*r_opt]) *r_opt = ir;
  }

  // Store the different rules as well as the scores in keypair mechanism

  set_keypair("rule_auto_scores",1,1,nrule,scores);
  set_keypair_int("rule_auto_allrules",1,nrule,NRULE,rules);
  set_keypair_int("rule_auto_best_rule",1,1,NRULE,&RULES(*r_opt,0));

  return(scores);
}

/****************************************************************************
 **
 ** FUNCTION: st_rule_glue
 **
 *****************************************************************************/
static void st_rule_glue(Relem *relem,
                         int    nrule1,
                         int    nbyrule1,
                         int   *rules1,
                         int   *fipos1)
{
  int *rules,*fipos,nrule,ir,nnew;

  if (relem == (Relem *) NULL) return;
  if (nrule1 <= 0) return;

  nrule = ir = relem->nrule;
  nnew  = nrule + nrule1;

  relem->rules = rules = (int *)
    mem_realloc((char *) relem->rules,sizeof(int) * NRULE  * nnew,1);
  relem->fipos = fipos = (int *)
    mem_realloc((char *) relem->fipos,sizeof(int) * NCOLOR * nnew,1);

  for (int i1=0; i1<nrule1; i1++,ir++)
  {
    for (int ic=0; ic<nbyrule1; ic++) RULES(ir,ic) = RULES1(i1,ic);
    for (int ic=0; ic<NCOLOR;   ic++) FIPOS(ir,ic) = FIPOS1(i1,ic);
  }

  relem->nrule   = nnew;
  relem->nbyrule = nbyrule1;
}

/****************************************************************************
 **
 ** FUNCTION: st_rule_product
 **
 *****************************************************************************/
static void st_rule_product(Split *split,
                            int    nprod,
                            int    nrule1,
                            int    nbyrule1,
                            int   *rules1,
                            int   *fipos1,
                            int    nrule2,
                            int    nbyrule2,
                            int   *rules2,
                            int   *fipos2)
{
  int *rules,*fipos,ir,ic,oper;
  int flag_debug = 0;

  split->rules = rules = (int *) mem_alloc(sizeof(int) * NRULE  * nprod,1);
  for (int i=0; i<NRULE  * nprod; i++) rules[i] = 0;
  split->fipos = fipos = (int *) mem_alloc(sizeof(int) * NCOLOR * nprod,1);
  for (int i=0; i<NCOLOR * nprod; i++) fipos[i] = 0;

  ir = 0;
  for (int i1=0; i1<nrule1; i1++)
    for (int i2=0; i2<nrule2; i2++,ir++)
    {
      ic   = 0;
      oper = split->oper;
      RULES(ir,ic++) = 1000 + oper;

      if (flag_debug)
      {
        message("Rule Product (with operator %d)\n",oper);
        st_rule_print(i1,nbyrule1,rules1,fipos1,0,-1,-1,TEST);
        st_rule_print(i2,nbyrule2,rules2,fipos2,0,-1,-1,TEST);
      }

      for (int i=0; i<nbyrule1; i++) RULES(ir,ic++) = RULES1(i1,i);
      for (int i=0; i<nbyrule2; i++) RULES(ir,ic++) = RULES2(i2,i);
      for (int i=0; i<NCOLOR;   i++)
      {
        if (FIPOS1(i1,i) > 0)
          FIPOS(ir,i) = FIPOS1(i1,i) * BASE + st_define_fipos(oper,1);
        if (FIPOS2(i2,i) > 0)
          FIPOS(ir,i) = FIPOS2(i2,i) * BASE + st_define_fipos(oper,0);
      }

      if (flag_debug)
      {
        message("Product result=");
        st_rule_print(ir,nbyrule1+nbyrule2+1,rules,fipos,0,-1,-1,TEST);
      }
    }
  split->nrule   = nprod;
  split->nbyrule = nbyrule1 + nbyrule2 + 1;
}

/****************************************************************************
 **
 ** FUNCTION: st_split_collapse
 **
 *****************************************************************************/
static void st_split_collapse(Split *split,
                              int    verbose)

{
  Relem *relem;
  int num[2],nby[2],*ptr[2],nprod;

  if (split == (Split *) NULL) return;

  // Explore the two Relems

  for (int i=0; i<2; i++)
    st_relem_explore(split->relems[i],verbose);

  // Prepare collapsing

  if (split->nrule <= 0)
  {
    for (int i=0; i<2; i++)
    {
      relem = split->relems[i];
      if (relem->nfacies <= 1)
      {
        num[i]   = 1;
        nby[i]   = 1;
        ptr[i]   = &relem->facies[0];
      }
      else
      {
        num[i]   = relem->nrule;
        nby[i]   = relem->nbyrule;
        ptr[i]   = relem->rules;
      }
    }
    
    // Merge the rules of the two Relem (by product)
    
    nprod = num[0] * num[1];
    split->nbyrule = nby[0] + nby[1] + 1;
    if (nprod > 0)
    {
      st_rule_product(split,nprod,
                      num[0],nby[0],ptr[0],split->relems[0]->fipos,
                      num[1],nby[1],ptr[1],split->relems[1]->fipos);
      if (verbose)
        st_rules_print("Split",split->nrule,split->nbyrule,
                       split->rules,split->fipos);
    }
  }
}

/****************************************************************************
 **
 ** FUNCTION: st_relem_explore
 **
 *****************************************************************************/
static void st_relem_explore(Relem *relem,
                             int verbose)
{
  Split *split;

  if (relem == (Relem *) NULL) return;

  for (int is=0; is<relem->nsplit; is++)
  {
    split = relem->splits[is];
    st_split_collapse(split,verbose);
    st_rule_glue(relem,split->nrule,split->nbyrule,
                 split->rules,split->fipos);
    if (verbose)
      st_rules_print("Relem",relem->nrule,relem->nbyrule,
                     relem->rules,relem->fipos);
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
GEOSLIB_API Vario_Order *vario_order_manage(int mode,
                                            int flag_dist,
                                            int size_aux,
                                            Vario_Order *vorder)
{
  Vario_Order *vorder_loc;

  /* Dispatch */

  vorder_loc = (Vario_Order *) NULL;
  switch (mode)
  {
    case 1:
      vorder_loc = (Vario_Order *) mem_alloc(sizeof(Vario_Order),0);
      if (vorder_loc == (Vario_Order *) NULL) return(vorder_loc);
      vorder_loc->npair        = 0;
      vorder_loc->nalloc       = 0;
      vorder_loc->size_aux     = size_aux;
      vorder_loc->flag_dist    = flag_dist;
      vorder_loc->tab_iech     = (int    *) NULL;
      vorder_loc->tab_jech     = (int    *) NULL;
      vorder_loc->tab_ipas     = (int    *) NULL;
      vorder_loc->tab_sort     = (int    *) NULL;
      vorder_loc->tab_aux_iech = (char   *) NULL;
      vorder_loc->tab_aux_jech = (char   *) NULL;
      vorder_loc->tab_dist     = (double *) NULL;
      break;
    
    case 0:
      vorder_loc = vorder;
      if (vorder == (Vario_Order *) NULL) return(vorder_loc);
      vorder_loc->tab_iech     = (int    *) 
        mem_free((char *) vorder_loc->tab_iech);
      vorder_loc->tab_jech     = (int    *) 
        mem_free((char *) vorder_loc->tab_jech);
      vorder_loc->tab_ipas     = (int    *) 
        mem_free((char *) vorder_loc->tab_ipas);
      vorder_loc->tab_sort     = (int    *) 
        mem_free((char *) vorder_loc->tab_sort);
      vorder_loc->tab_aux_iech = (char   *) 
        mem_free((char *) vorder_loc->tab_aux_iech);
      vorder_loc->tab_aux_jech = (char   *) 
        mem_free((char *) vorder_loc->tab_aux_jech);
      vorder_loc->tab_dist     = (double *) 
        mem_free((char *) vorder_loc->tab_dist);
      break;

    case -1:
      vorder_loc = vorder;
      if (vorder == (Vario_Order *) NULL) return(vorder_loc);
      vorder_loc->tab_iech     = (int    *) 
        mem_free((char *) vorder_loc->tab_iech);
      vorder_loc->tab_jech     = (int    *) 
        mem_free((char *) vorder_loc->tab_jech);
      vorder_loc->tab_ipas     = (int    *) 
        mem_free((char *) vorder_loc->tab_ipas);
      vorder_loc->tab_sort     = (int    *) 
        mem_free((char *) vorder_loc->tab_sort);
      vorder_loc->tab_aux_iech = (char   *) 
        mem_free((char *) vorder_loc->tab_aux_iech);
      vorder_loc->tab_aux_jech = (char   *) 
        mem_free((char *) vorder_loc->tab_aux_jech);
      vorder_loc->tab_dist     = (double *) 
        mem_free((char *) vorder_loc->tab_dist);
      vorder_loc = (Vario_Order *) mem_free((char *) vorder_loc);
      break;
  }
  return(vorder_loc);
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
GEOSLIB_API int vario_order_add(Vario_Order *vorder,
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

  if (vorder == (Vario_Order *) NULL) return(0);

  /* Resize the array */

  if (vorder->npair >= vorder->nalloc)
  {
    vorder->nalloc += VARIO_ORDER_QUANT;
    vorder->tab_iech = (int *)
      mem_realloc((char *) vorder->tab_iech, 
                  vorder->nalloc * sizeof(int),0);
    if (vorder->tab_iech == (int *) NULL) return(1);
    vorder->tab_jech = (int *)
      mem_realloc((char *) vorder->tab_jech, 
                  vorder->nalloc * sizeof(int),0);
    if (vorder->tab_jech == (int *) NULL) return(1);
    vorder->tab_ipas = (int *)
      mem_realloc((char *) vorder->tab_ipas, 
                  vorder->nalloc * sizeof(int),0);
    if (vorder->tab_ipas == (int *) NULL) return(1);
    vorder->tab_sort = (int *)
      mem_realloc((char *) vorder->tab_sort, 
                  vorder->nalloc * sizeof(int),0);
    if (vorder->tab_sort == (int *) NULL) return(1);
    if (vorder->size_aux > 0)
    {
      vorder->tab_aux_iech = (char *)
        mem_realloc((char *) vorder->tab_aux_iech, 
                    vorder->nalloc * vorder->size_aux,0);
      if (vorder->tab_aux_iech == (char *) NULL) return(1);
      vorder->tab_aux_jech = (char *)
        mem_realloc((char *) vorder->tab_aux_jech, 
                    vorder->nalloc * vorder->size_aux,0);
      if (vorder->tab_aux_jech == (char *) NULL) return(1);
    }
    if (vorder->flag_dist)
    {
      vorder->tab_dist = (double *)
        mem_realloc((char *) vorder->tab_dist, 
                    vorder->nalloc * sizeof(double),0);
      if (vorder->tab_dist == (double *) NULL) return(1);
    }
  }

  /* Add the new information */

  vorder->tab_iech[vorder->npair] = (dist > 0) ? iech : jech;
  vorder->tab_jech[vorder->npair] = (dist > 0) ? jech : iech;
  vorder->tab_ipas[vorder->npair] = ipas + idir * QUANT_DIR;
  if (vorder->flag_dist) vorder->tab_dist[vorder->npair] = dist;
  if (vorder->size_aux > 0)
  {
    iad = vorder->npair * vorder->size_aux;
    if (aux_iech != NULL)
      (void) memcpy(&vorder->tab_aux_iech[iad],aux_iech,vorder->size_aux);
    if (aux_jech != NULL)
      (void) memcpy(&vorder->tab_aux_jech[iad],aux_jech,vorder->size_aux);
  }
  vorder->npair++;
  return(0);
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
GEOSLIB_API void vario_order_print(Vario_Order *vorder,
                                   int idir_target,
                                   int ipas_target,
                                   int verbose)
{
  int i,j,ipas,idir,flag_first;

  if (vorder == (Vario_Order *) NULL) return;

  mestitle(0,"Variogram Order structure");
  message("Allocated size    = %d\n",vorder->nalloc);
  message("Number of pairs   = %d\n",vorder->npair);
  if (! verbose) return;
  flag_first = 1;

  for (i=0; i<vorder->npair; i++)
  {    
    j = (vorder->tab_sort == (int *) NULL) ? i : vorder->tab_sort[i];
    ipas = vorder->tab_ipas[j];
    idir = ipas / QUANT_DIR;
    ipas = ipas - QUANT_DIR * idir;
    if (idir_target >= 0 && idir != idir_target) continue;
    if (ipas_target >= 0 && ipas != ipas_target) continue;
    
    if (flag_first)
    {
      if (! vorder->flag_dist) 
        message("Rank - Dir - Lag - I - J\n");
      else
        message("Rank - Dir - Lag - I - J - Dist\n");
      flag_first = 0;
    }

    message("%5d",i+1);
    message(" %5d",idir+1);
    message(" %5d",ipas+1);
    message(" %5d",vorder->tab_iech[j]+1);
    message(" %5d",vorder->tab_jech[j]+1);
    if (vorder->flag_dist) message(" %lf",vorder->tab_dist[j]);
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
GEOSLIB_API Vario_Order *vario_order_final(Vario_Order *vorder,
                                           int *npair)
{
  int i,error;

  *npair = 0;
  if (vorder == (Vario_Order *) NULL) return(vorder);

  error = 0;
  if (vorder->npair > 0)
  {
    vorder->tab_iech = (int *) 
      mem_realloc((char *) vorder->tab_iech,vorder->npair * sizeof(int),0);
    if (vorder->tab_iech == (int *) NULL) error = 1;
    vorder->tab_jech = (int *) 
      mem_realloc((char *) vorder->tab_jech,vorder->npair * sizeof(int),0);
    if (vorder->tab_jech == (int *) NULL) error = 1;
    vorder->tab_ipas = (int *) 
      mem_realloc((char *) vorder->tab_ipas,vorder->npair * sizeof(int),0);
    if (vorder->tab_ipas == (int *) NULL) error = 1;
    vorder->tab_sort = (int *) 
      mem_realloc((char *) vorder->tab_sort,vorder->npair * sizeof(int),0);
    if (vorder->tab_sort == (int *) NULL) error = 1;
    if (vorder->flag_dist)
    {
      vorder->tab_dist = (double *) 
        mem_realloc((char *) vorder->tab_dist,vorder->npair *sizeof(double),0);
      if (vorder->tab_dist == (double *) NULL) error = 1;
    }
    if (vorder->size_aux > 0)
    {
      vorder->tab_aux_iech = (char *) 
        mem_realloc((char *) vorder->tab_aux_iech,
                    vorder->npair * vorder->size_aux,0);
      if (vorder->tab_aux_iech == (char *) NULL) error = 1;
      vorder->tab_aux_jech = (char *) 
        mem_realloc((char *) vorder->tab_aux_jech,
                    vorder->npair * vorder->size_aux,0);
      if (vorder->tab_aux_iech == (char *) NULL) error = 1;
    }
  }
  vorder->nalloc = vorder->npair;

  if (error) 
  {
    vorder = vario_order_manage(-1,vorder->flag_dist,vorder->size_aux,vorder);
    *npair = 0;
  }
  else if (vorder->npair > 0)
  {
    for (i=0; i<vorder->npair; i++) vorder->tab_sort[i] = i;
    ut_sort_int(1,vorder->npair,vorder->tab_sort,vorder->tab_ipas);
    *npair = vorder->npair;
  }
  return(vorder);
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
GEOSLIB_API void vario_order_get_indices(Vario_Order *vorder,
                                         int     ipair,
                                         int    *iech,
                                         int    *jech,
                                         double *dist)
{
  int jpair;
  
  if (vorder->tab_sort == (int *) NULL) messageAbort("vario_order_get_indices");
  jpair = vorder->tab_sort[ipair];
  *iech = vorder->tab_iech[jpair];
  *jech = vorder->tab_jech[jpair];
  *dist = (vorder->flag_dist) ? vorder->tab_dist[jpair] : TEST;
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
GEOSLIB_API void vario_order_get_auxiliary(Vario_Order *vorder,
                                           int     ipair,
                                           char   *aux_iech,
                                           char   *aux_jech)
{
  int jpair,iad;
  
  if (vorder->tab_sort == (int *) NULL) messageAbort("vario_order_get_auxiliary");
  jpair = vorder->tab_sort[ipair];
  iad   = vorder->size_aux * jpair;
  (void) memcpy(aux_iech,&vorder->tab_aux_iech[iad],vorder->size_aux);
  (void) memcpy(aux_jech,&vorder->tab_aux_jech[iad],vorder->size_aux);
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
GEOSLIB_API void vario_order_get_bounds(Vario_Order *vorder,
                                        int  idir,
                                        int  ipas,
                                        int *ifirst,
                                        int *ilast)
{
  int ipair,jpair,ival;

  ival = ipas + idir * QUANT_DIR;
  if (vorder->npair > 0 && 
      vorder->tab_sort == (int *) NULL) messageAbort("vario_order_get_bounds");
  *ifirst = vorder->npair;
  *ilast  = -1;
  for (ipair=0; ipair<vorder->npair; ipair++)
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
static int invgen(double *a,
                  int     neq,
                  double *tabout)
{
  double *eigvec,*eigval,value;
  int     i,j,k,error;
  
  /* Initializations */
  
  error  = 1;
  eigvec = eigval = (double *) NULL;
  
  /* Core allocation */
  
  eigval = (double *) mem_alloc(sizeof(double) * neq,0);
  if (eigval == (double *) NULL) goto label_end;
  eigvec = (double *) mem_alloc(sizeof(double) * neq * neq,0);
  if (eigvec == (double *) NULL) goto label_end;
  
  /* Calculate the eigen vectors */
  
  if (matrix_eigen(a,neq,eigval,eigvec)) goto label_end;

  /* Calculate the generalized inverse */

  for (i=0; i<neq; i++)
    for (j=0; j<neq; j++)
    {
      value = 0.;
      for (k=0; k<neq; k++)
      {
        if (ABS(eigval[k]) > 1e-10)
          value += EIGVEC(k,i) * EIGVEC(k,j) / eigval[k];
      }
      TABOUT(i,j) = value;
    }

  /* Set the error returned code */
  
  error = 0;
  
label_end:
  eigval = (double *) mem_free((char *) eigval);
  eigvec = (double *) mem_free((char *) eigvec);
  return(error);
}


/****************************************************************************/
/*!
** Calculate the indexes of each parameter 
**
** \param[in]  i Index
** \param[in]  j Index
** 
*****************************************************************************/



static int st_index(int i,
                    int j)
{
  int value;  

  value = 0;
  switch(i)
  {
    case 0:
      value = (j==0) ? 1 : 0;
      break;
      
    case 1:
      value = (j==0) ? 3 : 0;
      break;
      
    case 2:
      value = (j==0) ? 2 : 1;
      break;
      
    case 3:
      value = (j==0) ? 3 : 2;
      break;
      
  }
  return(value);
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
static void st_compute_params(Local_CorPgs * corpgs,
                              double *params_in,
                              double *params)
{
  double rho,rho2;
  rho = corpgs->rho;
  
  switch (corpgs->opt_correl)
  {
    case 0:
      params[0] = params_in[0];
      params[1] = params_in[1];
      params[2] = params_in[2]; 
      params[3] = params_in[3]; 
      break;
      
    case 1:			/* Symmetrical case */
      params[0] = params_in[0];
      params[1] = params[2] = params_in[1];
      params[3] = params_in[2]; 
      break;
      
    case 2:			/* Residual case */
      rho2 = rho * rho;
      params[0] = params_in[0];
      params[1] = params[2] = rho * params_in[0];
      params[3] = rho2 * params_in[0] +  (1. - rho2) * params_in[1];
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

  for (i=0; i<4; i++) params[i] = 0.;
  st_compute_params(corpgs,params_in,params);
  rho = corpgs->rho;

  for (i=0; i<4; i++) 
    M_R(correl,4,i,i) = 1.;
  
  M_R(correl,4,2,0) = rho;
  M_R(correl,4,3,1) = rho;
  
  M_R(correl,4,1,0) = params[0];
  M_R(correl,4,2,1) = params[2];
  M_R(correl,4,3,0) = params[1];
  M_R(correl,4,3,2) = params[3];
  
  matrix_fill_symmetry(4,correl);
}

/****************************************************************************/
/*!
**  Update the following matrices according to constraints on model
**
** \param[in]  corpgs      Local_CorPgs structure
** \param[in]  Grad        Vector of gradients (Dimension = 4)
** \param[in]  Hess        Matrix of Hessian (Dimension = 4*4)
** \param[in]  JJ          Matrix of t(JJ) * JJ (Dimension = 4*4)
**
** \param[out]  Grad        Vector of gradients (Dimension = npar)
** \param[out]  Hess        Matrix of Hessian (Dimension = npar * npar)
** \param[out]  JJ          Matrix of t(JJ) * JJ (Dimension = npar * npar)
**
*****************************************************************************/
static void st_update_constraints(Local_CorPgs *corpgs,
                                  double *Grad,
                                  double *Hess,
                                  double *JJ)
{
  double v[4],m[16],*modif;
  int i,j,k,l,npar;
 
  /* Initializations */

  modif = corpgs->modif;
  npar  = corpgs->npar;

  /* Update the Grad */
  
  matrix_combine(4,1.,Grad,0.,NULL,v);
  matrix_combine(4,0.,NULL,0.,NULL,Grad);
  for (i=0; i<npar; i++)
    for (j=0; j<4; j++)
      Grad[i] += v[j] * M_R(modif,4,i,j);

  /* Update the Hessian */

  matrix_combine(16,1.,Hess,0.,NULL,m);
  matrix_combine(16,0.,NULL,0.,NULL,Hess);
  for (i=0; i<npar; i++)
    for (j=0; j<npar; j++)
      for (k=0; k<4; k++)
        for (l=0; l<4; l++)
          M_R(Hess,npar,i,j) += 
            M_R(modif,4,i,k) * M_R(m,4,k,l) * M_R(modif,4,j,l);

  /* Update JJ */
  
  if(JJ!=NULL)
  {
    matrix_combine(16,1.,JJ,0.,NULL,m);
    matrix_combine(16,0.,NULL,0.,NULL,JJ);
    for (i=0; i<npar; i++)
      for (j=0; j<npar; j++)
        for (k=0; k<4; k++)
          for (l=0; l<4; l++)
            M_R(JJ,npar,i,j) += 
              M_R(modif,4,i,k) * M_R(m,4,l,k) * M_R(modif,4,j,l);
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
                           double * ev,
                           double * d1, 
                           double * d2)
{
  double temp[16], invGn[16];
  int i,j;
  
  matrix_combine(16,0.,NULL,0.,NULL,d2);
  matrix_combine(16,0.,NULL,0.,NULL,temp);
  st_build_correl(corpgs,corpgs->params,temp);
  matrix_combine(16,-1.,temp,0.,NULL,temp);
  
  for(i=0;i<4;i++)
    M_R(temp,4,i,i) += eigval;
  
 
  invgen(temp,4,invGn);
  
  for(i=0;i<4;i++)
    d1[i] = 2 * ev[12+F(i,0)]*ev[12+F(i,1)];

  for(i=0;i<4;i++)
    for(j=0;j<i;j++)
    {M_R(d2,4,i,j)  = ev[12+F(i,0)]*ev[12+F(j,0)] * M_R(invGn,4,F(i,1),F(j,1));
      M_R(d2,4,i,j) += ev[12+F(i,1)]*ev[12+F(j,0)] * M_R(invGn,4,F(i,0),F(j,1));
      M_R(d2,4,i,j) += ev[12+F(i,0)]*ev[12+F(j,1)] * M_R(invGn,4,F(i,1),F(j,0));
      M_R(d2,4,i,j) += ev[12+F(i,1)]*ev[12+F(j,1)] * M_R(invGn,4,F(i,0),F(j,0));
      M_R(d2,4,i,j) *= 2;
    }
  
  matrix_fill_symmetry(4,d2);
  st_update_constraints(corpgs,d1,d2,NULL);

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
  double rho,rho2;
  Local_CorPgs *corpgs;
  
  corpgs = &local_pgs->corpgs;
  rho    = corpgs->rho;

  switch (corpgs->opt_correl)
  {
    case 0:
      if (igrf == 0 && jgrf == 0)
        return(corpgs->params[0]);
      else if (igrf == 1 && jgrf == 1) 
        return(corpgs->params[3]);
      else
      {
        if (idir > 0) 
          return(corpgs->params[1]);
        else
          return(corpgs->params[2]);
      }
      break;
    
    case 1:
      if (igrf == 0 && jgrf == 0)
        return(corpgs->params[0]);
      else if (igrf == 1 && jgrf == 1) 
        return(corpgs->params[2]);
      else
        return(corpgs->params[1]);
      break;
      
    case 2:
      rho2 = rho * rho;
      if (igrf == 0 && jgrf == 0)
        return(corpgs->params[0]);
      else if (igrf == 1 && jgrf == 1) 
        return(corpgs->params[0] * rho2 + corpgs->params[1] * (1. - rho2));
      else
        return(corpgs->params[0] * rho);
      break;
  }
  return(0.);
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
  matrix_combine(16,0.,NULL,0.,NULL,modif);
  
  switch (corpgs->opt_correl)
  {
    case 0:			/* Full parameters */
      corpgs->npar = 4;
      for (i=0; i<4; i++) M_R(modif,4,i,i) = 1.;
      break;
    
    case 1:			/* Symmetrical case */
      corpgs->npar = 3;
      M_R(modif,4,0,0) = 1;
      M_R(modif,4,1,1) = 1;
      M_R(modif,4,1,2) = 1;
      M_R(modif,4,2,3) = 1;
      break;
    
    case 2:			/* Residual case */
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
** \param[in]  local_pgs   Local_Pgs structure
**
** \param[out]  local_pgs   Local_Pgs structure
**
*****************************************************************************/
static void st_define_corpgs(int option, 
                             int flag_rho,
                             double rho,
                             Local_Pgs *local_pgs)
{
  Local_CorPgs *corpgs;
  int    i;

  /* Initializations */

  corpgs = &local_pgs->corpgs;
  corpgs->rho         = rho;
  corpgs->opt_correl  = option;
  st_set_modif(corpgs);

  for (i=0; i<4; i++) corpgs->params[i] = 0.;

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
static int st_get_count(Local_Pgs *local_pgs,
                        int     ifac1,
                        int     ifac2)
{
  int ipair,i1,i2,number;
  double dist,w1,w2;

  /* Initializations */

  number = 0;

  for (ipair=local_pgs->ifirst; ipair<local_pgs->ilast; ipair++)
  {
    vario_order_get_indices(local_pgs->vorder,ipair,&i1,&i2,&dist);
    if (ifac1 != local_pgs->db->getVariable(i1,0)) continue;
    if (ifac2 != local_pgs->db->getVariable(i2,0)) continue;
    w1 = local_pgs->db->getWeight(i1);
    w2 = local_pgs->db->getWeight(i2);
    number += (int) (w1 * w2);
  }
  return(number);
}

/****************************************************************************/
/*!
** PRUPOSE: Internal function
**
*****************************************************************************/
static double st_rkl(int    maxpts,
                     double x,
                     double y,
                     double *lower,
                     double *upper,
                     double *corr1,
                     double *covar,
                     double *temp)
{
  double cste[2],vect[2],mean[2],v1,v2,S,error;
  int    inform;
  static double abseps = 1.e-12;
  static double releps = 0.;

  cste[0] = 0.;
  cste[1] = 0.;
  vect[0] = x;
  vect[1] = y;
  matrix_product(2,2,1,temp,vect,mean);
  v1 = law_df_bigaussian(vect,cste,corr1);
  mvndst2n(lower,upper,mean,covar,maxpts,abseps,releps,&error,&v2,&inform);
  S  = v1 * v2;
  return(S);
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
static double st_ikl(int     maxpts,
                     int     index1,
                     int     index2,
                     double *lower,
                     double *upper,
                     double *correl)
{
  double low[2],upp[2],value,x,y,S;
  double corr1[4],corr2[4],corrc[4],covar[4],inv_corr1[4],temp[4];
  int    i,j,k,index[2];

  // Initializations 
  index[0] = index1;
  index[1] = index2;

  // Build submatrices
  matrix_manage(4,1,-2, 0,index,NULL ,lower ,low);
  matrix_manage(4,1,-2, 0,index,NULL ,upper ,upp);
  matrix_manage(4,4, 2, 2,index,index,correl,corr1);
  matrix_manage(4,4,-2, 2,index,index,correl,corrc);
  matrix_manage(4,4,-2,-2,index,index,correl,corr2);
  if (matrix_invert_copy(corr1,2,inv_corr1)) messageAbort("st_ikl #1");
  matrix_product(2,2,2,corrc,inv_corr1,temp);

  // Derive covar
  for (i=0; i<2; i++)
    for (j=0; j<2; j++)
    {
      value  = 0;
      for (k=0; k<2; k++) value += M_R(temp,2,k,i) * M_R(corrc,2,k,j);
      M_R(covar,2,i,j) = M_R(corr2,2,i,j) - value;
    }

  S = 0.;
  x = upper[index1];
  if (IS_GAUSS_DEF(x))
  {
    y = upper[index2];
    if (IS_GAUSS_DEF(y)) S += st_rkl(maxpts,x,y,low,upp,corr1,covar,temp);
    y = lower[index2];
    if (IS_GAUSS_DEF(y)) S -= st_rkl(maxpts,x,y,low,upp,corr1,covar,temp);
  }
  x = lower[index1];
  if (IS_GAUSS_DEF(x))
  {
    y = lower[index2];
    if (IS_GAUSS_DEF(y)) S += st_rkl(maxpts,x,y,low,upp,corr1,covar,temp);
    y = upper[index2];
    if (IS_GAUSS_DEF(y)) S -= st_rkl(maxpts,x,y,low,upp,corr1,covar,temp);
  }
  return(S/2.);
}

/****************************************************************************/
/*!
** \return  Calculate the other derivatives
**
*****************************************************************************/
static double st_nkl(double *u,
                     double  lower,
                     double  upper,
                     double *invvari,
                     int     index2,
                     double  meanj,
                     double  varj,
                     double  stdj)
{
  double dflow,dfupp,cdflow,cdfupp,invval,invpart[3],total,S;

  dfupp  = law_dnorm(upper,meanj,stdj);
  dflow  = law_dnorm(lower,meanj,stdj);
  cdfupp = law_cdf_gaussian((upper - meanj) / stdj);
  cdflow = law_cdf_gaussian((lower - meanj) / stdj);
  invval = invvari[index2];
  matrix_manage(4,1,-1,0,&index2,NULL,invvari,invpart);
  matrix_product(1,3,1,invpart,u,&total);

  S = (dfupp - dflow) * varj * invval -
    (cdfupp - cdflow) * (invval * meanj + total);
  return(S);
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
static double st_d2_dkldkl(int     maxpts,
                           int     index1,
                           int     index2,
                           double *lower,
                           double *upper,
                           double *correl)
{
  double corri[16],v1,v2,S;
  double deltaparam = 1.e-6;
  int i;

  for (i=0; i<16; i++) corri[i] = correl[i];
  M_R(corri,4,index1,index2) += deltaparam;
  M_R(corri,4,index2,index1) += deltaparam;
  v1 = st_ikl(maxpts,index1,index2,lower,upper,corri);

  for (i=0; i<16; i++) corri[i] = correl[i];
  M_R(corri,4,index1,index2) -= deltaparam;
  M_R(corri,4,index2,index1) -= deltaparam;
  v2 = st_ikl(maxpts,index1,index2,lower,upper,corri);

  S = (v1 - v2) / (2. * deltaparam);
  return(S/2.);
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
static double st_d2_dkldij(double *lower,
                           double *upper,
                           double *correl)
{
  int    i,i1,i2,i3,i4,grid[4],flag_out;
  double u[4],S;

  S = 0.;
  for (i4=0; i4<2; i4++)
    for (i3=0; i3<2; i3++)
      for (i2=0; i2<2; i2++)
        for (i1=0; i1<2; i1++)
        {
          grid[0] = i1;
          grid[1] = i2;
          grid[2] = i3;
          grid[3] = i4;

          flag_out = 0;
          for (i=0; i<4 && flag_out==0; i++)
          {
            u[i] = (grid[i]) ? upper[i] : lower[i];
            flag_out = ! IS_GAUSS_DEF(u[i]);
          }
          if (! flag_out) 
            S += pow(-1.,i1+i2+i3+i4) * law_df_quadgaussian(u,correl);
        }
  return(S/2.);
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
static double st_d2_dkldkj(int     index1,
                           int     index2,
                           double *lower,
                           double *upper,
                           double *correl)
{
  double varcori[9],invvarcor[16],invvarcori[4],corr1[9],crosscor[3];
  double invcorr1[9],temp[3],lowi[3],uppi[3],u[3];
  double lowj,uppj,corr2,covar,sdcovar,S,mu,random;
  int i,i1,i2,i3,grid[3],flag_out;

  matrix_manage(4,4,-1,-1,&index2,&index2,correl,varcori);
  if (matrix_invert_copy(correl,4,invvarcor)) messageAbort("st_d2_dkldkj #1");
  matrix_manage(4,4, 1, 0,&index1,NULL,invvarcor,invvarcori);
  matrix_manage(4,4,-1,-1,&index2,&index2,correl,corr1);
  matrix_manage(4,4, 1,-1,&index2,&index2,correl,crosscor);
  matrix_manage(4,4, 1, 1,&index2,&index2,correl,&corr2);
  if (matrix_invert_copy(corr1,3,invcorr1)) messageAbort("st_d2_dkldkj #2");
  matrix_product(1,3,3,crosscor,invcorr1,temp);
  matrix_product(1,3,1,crosscor,temp,&covar);
  covar   = corr2 - covar ;
  sdcovar = sqrt(covar);
  matrix_manage(4,1,-1, 0,&index2,NULL,lower,lowi);
  matrix_manage(4,1,-1, 0,&index2,NULL,upper,uppi);
  matrix_manage(4,1, 1, 0,&index2,NULL,lower,&lowj);
  matrix_manage(4,1, 1, 0,&index2,NULL,upper,&uppj);
  
  S = 0.;
  for (i3=0; i3<2; i3++)
    for (i2=0; i2<2; i2++)
      for (i1=0; i1<2; i1++)
      {
        grid[0] = i1;
        grid[1] = i2;
        grid[2] = i3;
        
        for (i=flag_out=0; i<3 && flag_out==0; i++)
        {
          u[i] = (grid[i]) ? uppi[i] : lowi[i];
          flag_out = ! IS_GAUSS_DEF(u[i]);
        }
        if (flag_out) continue;
        
        matrix_product(1,3,1,temp,u,&mu);
        random = law_df_multigaussian(3,u,varcori);
        
        S += pow(-1.,3-i1+i2+i3) * random * 
          st_nkl(u,lowj,uppj,invvarcori,index2,mu,covar,sdcovar);
      }
  return(S/2);
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
** \param[in]  params       Array of parameters
** \param[in]  correl       Correlation matrix updated
**
** \param[out] Grad         Vector of cumulated gradients (Dimension= 4)
** \param[out] Hess         Matrix of cumulated Hessian (Dimension= 4*4)
** \param[out] JJ           Matrix of cumulated JJ (Dimension= 4*4)
**
*****************************************************************************/
static double st_calcul_stat(Local_Pgs *local_pgs,
                             int        flag_deriv,
                             int        flag_reset,
                             double    *params,
                             double    *correl,
                             double    *Grad,
                             double    *Hess,
                             double    *JJ)
{
  double grad[4],lower[4],upper[4],hess[16],gradgrad[16];
  double s,S,rj2,erval,ggval;
  int    ifac1,ifac2,nfifj,i,j,inform;
  static double abseps = 1.e-6;
  static double releps = 0.;
  static int maxpts4   = 10000;
  static int maxpts2   = 4000;

  S = 0.;
  matrix_combine(4 ,0.,NULL,0.,NULL,grad);
  matrix_combine(16,0.,NULL,0.,NULL,hess);

  for (ifac1=0; ifac1<local_pgs->nfacies; ifac1++)
  {
    for (ifac2=0; ifac2<local_pgs->nfacies; ifac2++)
    {
      nfifj = st_get_count(local_pgs,ifac1+1,ifac2+1);
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
        mvndst4(lower,upper,correl,maxpts4,abseps,releps,&erval,&s,&inform);
        STAT_PROBA(ifac1,ifac2) = s;
      }
      else
        s = STAT_PROBA(ifac1,ifac2);
      
      rj2 = -2. * log(s);
      S  += (double) nfifj * rj2;
      
      /* Calculate the derivative */
      
      if (! flag_deriv) continue;
      grad[0] = -st_ikl(maxpts2,0,1,lower,upper,correl) / s;
      grad[1] = -st_ikl(maxpts2,0,3,lower,upper,correl) / s;
      grad[2] = -st_ikl(maxpts2,1,2,lower,upper,correl) / s;
      grad[3] = -st_ikl(maxpts2,2,3,lower,upper,correl) / s;
      matrix_product(4,1,4,grad,grad,gradgrad);
      
      M_R(hess,4,3,0) = M_R(hess,4,2,1) = -st_d2_dkldij(lower,upper,correl);
      M_R(hess,4,1,0) = st_d2_dkldkj(0,2,lower,upper,correl);
      M_R(hess,4,2,0) = st_d2_dkldkj(1,3,lower,upper,correl);
      M_R(hess,4,3,1) = st_d2_dkldkj(3,1,lower,upper,correl);
      M_R(hess,4,3,2) = st_d2_dkldkj(2,0,lower,upper,correl);
      M_R(hess,4,0,0) = st_d2_dkldkl(maxpts4,0,1,lower,upper,correl);
      M_R(hess,4,1,1) = st_d2_dkldkl(maxpts4,0,3,lower,upper,correl);
      M_R(hess,4,2,2) = st_d2_dkldkl(maxpts4,1,2,lower,upper,correl);
      M_R(hess,4,3,3) = st_d2_dkldkl(maxpts4,2,3,lower,upper,correl);
      matrix_fill_symmetry(4,hess);

      for (i=0; i<4; i++)
      {
        Grad[i] += nfifj * grad[i];
        for (j=0; j<4; j++)
        {
          ggval = M_R(gradgrad,4,i,j);
          M_R(Hess,4,i,j) -= nfifj * (M_R(hess,4,i,j) / s - ggval);
          M_R(JJ,4,i,j)   += nfifj * ggval / rj2;
        }
      }
    }
  }
  return(S / 2.);
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
** \param[in]  params       Array of parameters
** \param[in]  correl       Correlation matrix updated
**
** \param[out] Grad         Vector of cumulated gradients (Dimension= 4)
** \param[out] Hess         Matrix of cumulated Hessian (Dimension= 4*4)
** \param[out] JJ           Matrix of cumulated JJ (Dimension= 4*4)
**
*****************************************************************************/
static double st_calcul_nostat(Local_Pgs *local_pgs,
                               int        flag_deriv,
                               int        flag_reset,
                               double    *params,
                               double    *correl,
                               double    *Grad,
                               double    *Hess,
                               double    *JJ)
{
  double  grad[4],lower[4],upper[4],hess[16],gradgrad[16];
  double  s,S,rj2,erval,ggval,dist,w1,w2;
  int     i1,i2,ifac1,ifac2,i,j,inform,ipair;
  static  double abseps = 1.e-6;
  static  double releps = 0.;
  static  int maxpts4   = 10000;
  static  int maxpts2   = 4000;

  S = 0.;
  matrix_combine(4 ,0.,NULL,0.,NULL,grad);
  matrix_combine(16,0.,NULL,0.,NULL,hess);

  for (ipair=local_pgs->ifirst; ipair<local_pgs->ilast; ipair++)
  {
    vario_order_get_indices(local_pgs->vorder,ipair,&i1,&i2,&dist);
    ifac1 = (int) local_pgs->db->getVariable(i1,0);
    ifac2 = (int) local_pgs->db->getVariable(i2,0);
    w1 = local_pgs->db->getWeight(i1);
    w2 = local_pgs->db->getWeight(i2);
    /* Get the bounds */

    (void) rule_thresh_define(local_pgs->propdef,local_pgs->db,local_pgs->rule,
                              ifac1,i1,0,0,1,
                              &lower[0],&upper[0],&lower[2],&upper[2]);
    (void) rule_thresh_define(local_pgs->propdef,local_pgs->db,local_pgs->rule,
                              ifac2,i2,0,0,1,
                              &lower[1],&upper[1],&lower[3],&upper[3]);
    
    if (flag_reset)
    {
      mvndst4(lower,upper,correl,maxpts4,abseps,releps,&erval,&s,&inform);
      MEMINT(ipair) = s;
    }
    else
      s = MEMINT(ipair);
    rj2 = -2. * log(s);
    S  += w1 * w2 * rj2;

    /* Calculate the derivative */

    if (! flag_deriv) continue;
    grad[0] = -st_ikl(maxpts2,0,1,lower,upper,correl) / s;
    grad[1] = -st_ikl(maxpts2,0,3,lower,upper,correl) / s;
    grad[2] = -st_ikl(maxpts2,1,2,lower,upper,correl) / s;
    grad[3] = -st_ikl(maxpts2,2,3,lower,upper,correl) / s;
    matrix_product(4,1,4,grad,grad,gradgrad);

    M_R(hess,4,3,0) = M_R(hess,4,2,1) = -st_d2_dkldij(lower,upper,correl);
    M_R(hess,4,1,0) = st_d2_dkldkj(0,2,lower,upper,correl);
    M_R(hess,4,2,0) = st_d2_dkldkj(1,3,lower,upper,correl);
    M_R(hess,4,3,1) = st_d2_dkldkj(3,1,lower,upper,correl);
    M_R(hess,4,3,2) = st_d2_dkldkj(2,0,lower,upper,correl);
    M_R(hess,4,0,0) = st_d2_dkldkl(maxpts4,0,1,lower,upper,correl);
    M_R(hess,4,1,1) = st_d2_dkldkl(maxpts4,0,3,lower,upper,correl);
    M_R(hess,4,2,2) = st_d2_dkldkl(maxpts4,1,2,lower,upper,correl);
    M_R(hess,4,3,3) = st_d2_dkldkl(maxpts4,2,3,lower,upper,correl);
    matrix_fill_symmetry(4,hess);
  
    for (i=0; i<4; i++)
    {
      Grad[i] += w1 * w2 * grad[i];
      for (j=0; j<4; j++)
      {
        ggval = M_R(gradgrad,4,i,j);
        M_R(Hess,4,i,j) -= w1 * w2 * (M_R(hess,4,i,j) / s - ggval);
        M_R(JJ,4,i,j)   +=  w1 * w2 * ggval / rj2;
      }
    }
  }
  return(S / 2.);
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
                        int        flag_deriv,
                        int        flag_reset,
                        double    *params,
                        double    *Grad,
                        double    *Hess,
                        double    *JJ)
{
  double S,correl[16];

  S = 0.;
 
  if (flag_deriv)
  {
    matrix_combine(4 ,0.,NULL,0.,NULL,Grad);
    matrix_combine(16,0.,NULL,0.,NULL,Hess);
    matrix_combine(16,0.,NULL,0.,NULL,JJ);
  }

  st_build_correl(&local_pgs->corpgs,params,correl);
    
  if (local_pgs->flag_stat)
    S = st_calcul_stat(local_pgs,flag_deriv,flag_reset,
                       params,correl,Grad,Hess,JJ);
  else
    S = st_calcul_nostat(local_pgs,flag_deriv,flag_reset,
                         params,correl,Grad,Hess,JJ);
  
  /* Modify the results due to the constraints on parameters */
  
  if (flag_deriv)
    st_update_constraints(&local_pgs->corpgs,Grad,Hess,JJ);
  
  return(S / 2.);
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
  double Hess[16],JJ[16],eigval[4],eigvec[16],Gn[16],invGn[16],correl[16],d2[16];
  double hsd[4],hgn[4],step[4],a[4],gr[4],param_temp[4],hgna[4],Grad[4],d1[4];
  double normgrad,normgrad2,alpha,delta2,c,a2,beta,hgna2,niter;
  double Sr,Snew,delta,mdiminution,mdiminution_pred,stepgr,rval;
  int    flag_sortie,flag_moved,npar,npar2,i;

  static double maxiter = 100;
  static double delta0  = 1;
  
  double Spen,Srpen;
  double penalize = 1000;
  int barrier = 0;

  /* Initializations */

  corpgs      = &local_pgs->corpgs;
  npar        = corpgs->npar;
  npar2       = npar * npar;
  delta       = delta0;
  mdiminution = 0.;
 
  if(new_val)
    st_initialize_params(corpgs);

  matrix_combine(npar ,1.,corpgs->params,0.,NULL,param_temp);
  matrix_combine(npar ,0.,NULL,0.,NULL,Grad);
  matrix_combine(npar2,0.,NULL,0.,NULL,Hess);
  matrix_combine(npar2,0.,NULL,0.,NULL,JJ);

  /* Calculate the score and the derivatives */

  Sr = st_calcul(local_pgs,1,1,corpgs->params,Grad,Hess,JJ);

  niter       = 0.;
  flag_sortie = 0;
  flag_moved  = 1;

  while (! flag_sortie)
  {
    if (barrier)
    {
      st_build_correl(corpgs,param_temp,correl);
      matrix_eigen(correl,4,eigval,eigvec);
      st_deriv_eigen(corpgs,eigval[3],eigvec,d1,d2);
      Srpen = Sr - penalize * log(eigval[3]);
      matrix_combine(npar,1,Grad,-penalize/eigval[3],d1,Grad);
      matrix_combine(npar,npar,Hess,penalize/(eigval[3]*eigval[3]),d2,Hess);
      matrix_combine(npar,npar,JJ,penalize/(eigval[3]*eigval[3]),d2,JJ);
      penalize *= 0.5;
    }
    niter  ++;
    delta2 = delta * delta;
    if (flag_moved)
    {
      matrix_combine(npar ,1.,Grad,0.,NULL,gr);
      matrix_combine(npar2,1.,Hess,0.,NULL,Gn);
      if (! is_matrix_definite_positive(npar,Gn,eigval,eigvec,0))
        matrix_combine(npar2,1.,JJ,0.,NULL,Gn);
      matrix_combine(npar,-1.,gr,0.,NULL,hsd);
      if (invgen(Gn,npar,invGn)) messageAbort("st_optim_lag");
      matrix_product(npar,npar,1,invGn,hsd,hgn);
    }

    /* Determine the lag (hgn, alpha*hsd) or a convex combinaison of both */
    
    if (matrix_norm(hgn,npar) <= delta2)
    {
      matrix_combine(npar,1.,hgn,0.,NULL,step);
    }
    else
    {
      normgrad2 = matrix_norm(gr,npar);
      alpha = normgrad2 / matrix_normA(gr,Gn,npar,npar);
      normgrad = sqrt(normgrad2);
      if (normgrad > (delta / alpha))
      {
        matrix_combine(npar,delta / normgrad,hsd,0.,NULL,step);
      }
      else
      {
        matrix_combine(npar,alpha,hsd,0.,NULL,a);
        matrix_combine(npar,1.,hgn,-1.,a,hgna);
        matrix_product(1,npar,1,a,hgn,&c);
        a2    = matrix_norm(a,npar);
        hgna2 = matrix_norm(hgna,npar);
        if (c <= 0.)
          beta = (-c + sqrt(c*c + hgna2 * (delta2 - a2))) / hgna2;
        else
          beta = (delta2 - a2) / (c + sqrt(c*c + hgna2 * (delta2 - a2)));
        matrix_combine(npar,beta,hgn,(1.-beta),a,step);
      }
    }
    
    matrix_combine(npar,1.,step,1.,corpgs->params,param_temp);
    st_build_correl(corpgs,param_temp,correl);
    while (! is_matrix_definite_positive(4,correl,eigval,eigvec,0))
    {
      matrix_combine(npar,0.9,step,0.,NULL,step);
      matrix_combine(npar,1.,step,1.,corpgs->params,param_temp);
      st_build_correl(corpgs,param_temp,correl);
    }

    Snew = st_calcul(local_pgs,0,1,param_temp,Grad,Hess,JJ);
    
    if(barrier)
      Spen = Snew - penalize * log(eigval[3]);
    if (! FFFF(Snew))
    {
      
      mdiminution = Snew - Sr;
      if(barrier)
        mdiminution = Spen -Srpen;
      matrix_product(1,npar,1,step,gr,&stepgr);
      mdiminution_pred = stepgr + 0.5 * matrix_normA(step,Gn,npar,npar);
      rval = mdiminution / mdiminution_pred;
      flag_moved = (mdiminution < 0);
    }
    else
    {
      flag_moved= 0;
      rval = 0.;
    }
    
    if (flag_moved)
    {
      Sr = Snew;
      Srpen = Spen;
      matrix_combine(npar,1,param_temp,0,NULL,corpgs->params);
      Snew = st_calcul(local_pgs,1,0,corpgs->params,Grad,Hess,JJ);
      if (barrier)
      { 
        st_deriv_eigen(corpgs,eigval[3],eigvec,d1,d2);
        matrix_combine(npar,1,Grad,penalize/eigval[3],d1,Grad);
        matrix_combine(npar,npar,Hess,-penalize/(eigval[3]*eigval[3]),d2,Hess);
        matrix_combine(npar,npar,JJ,-penalize/(eigval[3]*eigval[3]),d2,JJ);
        penalize/= 2.;
      }
      if (rval > 0.75) delta = MAX(delta, 3. * sqrt(matrix_norm(step,npar)));
    }
    if (rval < 0.25) delta /= 2.;
    
    flag_sortie = 
      (matrix_norminf(npar,step) < tolsort ||
       niter == maxiter                    ||
       matrix_norminf(npar,Grad) < 0.05    ||
       (fabs(mdiminution) < tolsort && flag_moved));
  }
  
  /* Returning arguments */
  
  if (debug_query("converge"))
  {
    message("Lag %d - S = %lf Parameters =",local_pgs->ipascur,Sr);
    for (i=0; i<corpgs->npar; i++) message(" %lf",corpgs->params[i]);
    message("\n");
  }

  /* Store the trace */

  trace_add_row(local_pgs);
  trace_define(local_pgs,niter,Snew,0,npar,Grad);

  return(Sr);
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
static int st_discard_point(Local_Pgs *local_pgs,
                            int iech)
{
  int    ifac;
  double low,up;

  /* This function is bypassed in the stationary case */

  if (local_pgs->flag_stat) return(0);

  /* The following checks must not be performed if not on facies */

  if (! local_pgs->flag_facies) return(0);

  /* Check on the facies */

  ifac = (int) local_pgs->db->getVariable(iech,0);
  if (ifac < 1 || ifac > local_pgs->nfacies) return(1);

  /* Check on the thresholds */

  if (! TEST_DISCRET)
  {
    if (local_pgs->db->getIntervalNumber() <= 0) return(0);
    low = local_pgs->db->getLowerBound(iech,local_pgs->igrfcur);
    up  = local_pgs->db->getUpperBound(iech,local_pgs->igrfcur);
  }
  else
  {
    if (get_LOCATOR_NITEM(local_pgs->db,LOC_RKLOW) <= 0 &&
        get_LOCATOR_NITEM(local_pgs->db,LOC_RKUP)  <= 0) return(0);
    low = local_pgs->db->getLowerInterval(iech,local_pgs->igrfcur);
    up  = local_pgs->db->getUpperInterval(iech,local_pgs->igrfcur);
  }
  if (up <= low) return(1);
  
  return(0);
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

  local_pgs->vorder = vario_order_final(local_pgs->vorder,&npair);
  if (local_pgs->vorder == (Vario_Order *) NULL) return(1);
  if (npair > 0 && ! local_pgs->flag_stat)
  {
    local_pgs->memint.resize(npair);
  }

  return(0);
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
                                              int    idir)
{
  int ipas,igrf,jgrf,iad,ngrf;

  const Dir& dir  = vario->getDirs(idir);
  ngrf = local_pgs->ngrf;

  for (ipas=0; ipas<dir.getNPas(); ipas++)
    for (igrf=0; igrf<ngrf; igrf++)
      for (jgrf=0; jgrf<=igrf; jgrf++)
      {
        iad = dir.getAddress(igrf,jgrf,ipas,false,1);
        vario->setGg(idir,iad,st_param_expand(local_pgs,igrf,jgrf,1));
        if (dir.getSw(iad) > 0.)
          vario->setHh(idir,iad, dir.getHh(iad) / dir.getSw(iad));
        iad = dir.getAddress(igrf,jgrf,ipas,false,-1);
        vario->setGg(idir,iad,st_param_expand(local_pgs,igrf,jgrf,-1));
        if (dir.getSw(iad) > 0.)
          vario->setHh(idir,iad, dir.getHh(iad) / dir.getSw(iad));
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
                                            int    idir)
{
  int *rindex,iech,jech,iiech,jjech,nech,ipas,iad,ivar,jvar,nvar,error;
  Db  *db;
  double psmin,ps,dist,maxdist;

  /* Retrieve information from Local_pgs structure */

  error   = 1;
  rindex  = (int *) NULL;
  db      = local_pgs->db;
  nech    = get_NECH(db);
  nvar    = vario->getVariableNumber();
  const Dir& dir = vario->getDirs(idir);
  maxdist = variogram_maximum_distance(dir);

  /* Initializations */

  ps    = 0.;
  psmin = _variogram_convert_angular_tolerance(dir.getTolAngle());

  /* Sort the data */

  rindex = variogram_sort(db);

  /* Loop on the first point */

  for (iiech=0; iiech<nech-1; iiech++)
  {
    iech = rindex[iiech];
    if (! db->isActive(iech)) continue;
    if (FFFF(db->getWeight(iech))) continue;
    if (st_discard_point(local_pgs,iech)) continue;
    mes_process("Calculating Variogram Geometry",nech,iech);

    for (jjech=iiech+1; jjech<nech; jjech++)
    {
      jech = rindex[jjech];
      if (variogram_maximum_dist1D_reached(db,iech,jech,maxdist)) break;
      if (! db->isActive(jech)) continue;
      if (FFFF(db->getWeight(jech))) continue;
      if (st_discard_point(local_pgs,jech)) continue;
    
      /* Check if the pair must be kept (Code criterion) */

      if (code_comparable(db,db,iech,jech,dir.getOptionCode(),
                          (int) dir.getTolCode())) continue;

      /* Check if the pair must be kept */

      dist = distance_intra(db,iech,jech,NULL);
      if (variogram_reject_pair(db,iech,jech,dist,psmin,
                                dir.getBench(),dir.getCylRad(),
                                dir.getCodir(),&ps)) continue;

      /* Get the rank of the lag */

      ipas = variogram_get_lag(vario,dir,ps,psmin,&dist);
      if (IFFFF(ipas)) continue;

      /* Add the sample (only positive lags are of interest) */

      if (ipas < 0) ipas = -ipas;
      if (vario_order_add(local_pgs->vorder,
                          iech,jech,NULL,NULL,ipas,idir,dist)) goto label_end;
      dist = ABS(dist);

      /* Update the distance and weight for all GRFs */
      
      for (ivar=0; ivar<nvar; ivar++)
        for (jvar=0; jvar<=ivar; jvar++)
        {
          if (vario->getFlagAsym())
          {
            iad = dir.getAddress(ivar,jvar,ipas,false,1);
            vario->setGg(idir, iad, TEST);
            vario->setHh(idir, iad, dir.getHh(iad) - dist);
            vario->setSw(idir, iad, dir.getSw(iad) + 1);
            iad = dir.getAddress(ivar,jvar,ipas,false,-1);
            vario->setGg(idir, iad, TEST);
            vario->setHh(idir, iad, dir.getHh(iad) + dist);
            vario->setSw(idir, iad, dir.getSw(iad) + 1);
          }      
          else
          {
            iad = dir.getAddress(ivar,jvar,ipas,false,0);
            vario->setGg(idir, iad, TEST);
            vario->setHh(idir, iad, dir.getHh(iad) + dist);
            vario->setSw(idir, iad, dir.getSw(iad) + 1);
          }
        }
    }
  }
  
  /* Set the error return code */

  error = 0;

label_end:
  rindex = (int *) mem_free((char *) rindex);
  return(error);
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
static void st_set_opt_correl(int opt,Local_CorPgs *corpgs)

{
  double params[4];
  int i = 0;
  
  for(i=0;i<4;i++)
    params[i] = 0;

  st_compute_params(corpgs,corpgs->params,params);

  switch(opt)
  {
    case 0:
      corpgs->params[0] = params[0];
      corpgs->params[1] = params[1];
      corpgs->params[2] = params[2];
      corpgs->params[3] = params[3];
      break;

    case 1:
      corpgs->params[0] = params[0];
      corpgs->params[1] = corpgs->params[2] = (params[1] + params[2])/2.;
      corpgs->params[2] = params[3];
      break;
  
    case 2:
      corpgs->params[0] = params[0];
      corpgs->params[1] = params[3];
      break;

  }
  
  corpgs->opt_correl=opt;
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
static double st_varcalc_correlated_grf(Local_Pgs *local_pgs,
                                        int        idir)
{
  double value;
  int    ipas,iad,igrf,jgrf,opt_temp;
  Vario  *vario;

  opt_temp = local_pgs->corpgs.opt_correl;
  value = 0.; 
  vario = local_pgs->vario;
  const Dir& dir = vario->getDirs(idir);

  for (ipas=0; ipas<dir.getNPas(); ipas++)
  {
    mes_process("Inverting Variogram Lag",dir.getNPas(),ipas);
    local_pgs->ipascur = ipas;
    trace_add_row(local_pgs);
    if (! LAG_USED(ipas)) continue;
    vario_order_get_bounds(local_pgs->vorder,idir,ipas,
                           &local_pgs->ifirst,&local_pgs->ilast);
    if (local_pgs->ifirst >= local_pgs->ilast) continue;
    
    if (opt_temp!=2)
      st_set_opt_correl(2,&local_pgs->corpgs);
    
    st_optim_onelag_pgs(local_pgs,1e-3,1);
    st_set_opt_correl(opt_temp,&local_pgs->corpgs);
    value += (dir.getUtilize(dir.getNPas() + ipas) *
              st_optim_onelag_pgs(local_pgs,1e-3,0));
  
    for (igrf=0; igrf<local_pgs->ngrf; igrf++)
      for (jgrf=0; jgrf<=igrf; jgrf++)
      {
        iad = dir.getAddress(igrf,jgrf,ipas,false,1);
        vario->setGg(idir, iad,st_param_expand(local_pgs,igrf,jgrf,1));
        iad = dir.getAddress(igrf,jgrf,ipas,false,-1);
        vario->setGg(idir,iad,st_param_expand(local_pgs,igrf,jgrf,-1));
      }
  }
  return(value);
}

/****************************************************************************/
/*!
**  Manage the Local_CorPgs structure
**
** \param[in]  local_corpgs  Local_CorPgs structure
**
** \param[out]  local_corpgs  Local_CorPgs structure
**
*****************************************************************************/
static void st_manage_corpgs(Local_CorPgs *local_corpgs)

{
  int i;

  local_corpgs->opt_correl = 0;
  local_corpgs->npar       = 0;
  local_corpgs->flag_rho   = 0;
  local_corpgs->rho        = 0.;
  
  for (i=0; i<4; i++)
    local_corpgs->params[i] = 0.;
  for (i=0; i<16; i++)
    local_corpgs->modif[i] = 0.;
}

/****************************************************************************/
/*!
**  Manage the Local_TracePgs structure
**
** \param[in]  local_tracepgs  Local_TracePgs structure
**
** \param[out]  local_tracepgs  Local_TracePgs structure
**
*****************************************************************************/
static void st_manage_trace(Local_TracePgs *local_tracepgs)
{
  local_tracepgs->flag_trace = 0;
  local_tracepgs->idir       = 0;
  local_tracepgs->ipas       = 0;
  local_tracepgs->nrow       = 0;
  local_tracepgs->ncol       = 0;
  local_tracepgs->trace      = VectorDouble();
}

/****************************************************************************/
/*!
**  Manage the Local_Pgs structure
**
** \param[in]  mode         0 initialization; 1 allocation; -1 deallocation
** \param[in]  local_pgs    Local_Pgs structure
** \param[in]  db           Db structure
** \param[in]  rule         Lithotype Rule definition
** \param[in]  vario        Vario structure
** \param[in]  varioind     Indicator Vario structure
** \param[in]  model        Model structure
** \param[in]  propdef      Props structure
** \param[in]  flag_stat    1 for stationary; 0 otherwise
** \param[in]  flag_facies  1 when processed on facies; 0 otherwise
** \param[in]  flag_dist    1 if distances are stored; 0 otherwise
** \param[in]  ngrf         Number of GRFs
** \param[in]  nfacies      Number of facies
** \param[in]  covtype      Type of the calculation (covariance, variogram, ...)
**
** \param[out] local_pgs    Local_Pgs structure
**
*****************************************************************************/
static void st_manage_pgs(int        mode,
                          Local_Pgs *local_pgs,
                          Db        *db,
                          Rule      *rule,
                          Vario     *vario,
                          Vario     *varioind,
                          Model     *model,
                          Props     *propdef,
                          int        flag_stat,
                          int        flag_facies,
                          int        flag_dist,
                          int        ngrf,
                          int        nfacies,
                          int        covtype)
{
  int ndim,ncova;

  /* Initializations */
  
  ndim = ncova = -1;

  /* Dispatch */

  switch (mode)
  {
    case 0:
      local_pgs->db          = (Db    *) NULL;
      local_pgs->rule        = (Rule  *) NULL;
      local_pgs->propdef     = (Props *) NULL;
      local_pgs->flag_stat   = 0;
      local_pgs->flag_facies = 0;
      local_pgs->covtype     = 0;
      local_pgs->igrfcur     = 0;
      local_pgs->idircur     = 0;
      local_pgs->ipascur     = 0;
      local_pgs->ngrf        = 0;
      local_pgs->npair       = 0;
      local_pgs->nfacies     = 0;
      local_pgs->nfac2       = 0;
      local_pgs->ifirst      = 0;
      local_pgs->ilast       = 0;
      local_pgs->d0          = VectorDouble();
      local_pgs->d1          = VectorDouble();
      local_pgs->memint      = VectorDouble();
      local_pgs->stat_proba  = VectorDouble();
      local_pgs->stat_thresh = VectorDouble();
      local_pgs->model       = (Model          *) NULL;
      local_pgs->vario       = (Vario          *) NULL;
      local_pgs->varioind    = (Vario          *) NULL;
      local_pgs->vorder      = (Vario_Order    *) NULL;
      break;

    case 1:
      local_pgs->db          = db;
      local_pgs->rule        = rule;
      local_pgs->propdef     = propdef;
      local_pgs->flag_stat   = flag_stat;
      local_pgs->flag_facies = flag_facies;
      local_pgs->covtype     = covtype;
      local_pgs->igrfcur     = 0;
      local_pgs->ipascur     = 0;
      local_pgs->ngrf        = ngrf;
      local_pgs->npair       = 0;
      local_pgs->nfacies     = nfacies;
      local_pgs->nfac2       = nfacies * nfacies;
      local_pgs->vario       = vario;
      local_pgs->varioind    = varioind;
      local_pgs->model       = model;
      if (model != (Model *) NULL)
      {
        ndim  = model->getDimensionNumber();
        ncova = model->getCovaNumber();
        local_pgs->d0.resize(ndim);
        local_pgs->d1.resize(ndim);
      }
      local_pgs->vorder      = vario_order_manage(1,flag_dist,0,NULL);
      if (flag_stat)
      {
        local_pgs->stat_proba.resize(local_pgs->nfac2,0.);
        local_pgs->stat_thresh.resize(nfacies * ngrf * 2,0.);
      }
      st_manage_corpgs(&local_pgs->corpgs);
      st_manage_trace (&local_pgs->tracepgs);
      break;
      
    case -1:
      local_pgs->vorder  = vario_order_manage(-1,0,0,local_pgs->vorder);
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
** \param[in]  db            Db structure
** \param[in]  vario         Vario structure for the GRFs to be filled
** \param[in]  rule          Lithotype Rule definition
** \param[in]  propdef       Props structure
** \param[in]  local_pgs     Local_Pgs structure
** \param[in]  ngrf          Number of GRFs
** \param[in]  opt_correl    0 full model; 1 symetrical; 2 residuals
** \param[out] flag_geometry 1 if Geometry must be established per direction
**                           0 if Geometry is already calculated before
**                             calling this function
**
*****************************************************************************/
static int st_variopgs_calcul_norho(Db    *db,
                                    Vario *vario,
                                    Rule  *rule,
                                    Props *propdef,
                                    Local_Pgs *local_pgs,
                                    int ngrf,
                                    int opt_correl,
                                    int flag_geometry)
{
 
  int idir;
  
  st_set_rho(rule->getRho(),local_pgs);
    
  /* Loop on the directions */

  for (idir=0; idir<vario->getDirectionNumber(); idir++)
  {
    local_pgs->idircur = idir;

    /* Establish the geometry */

    if (flag_geometry)
    {
      if (st_variogram_geometry_pgs_calcul(local_pgs,vario,idir)) return(1);
      st_variogram_geometry_pgs_correct(local_pgs,vario,idir);
      if (st_variogram_geometry_pgs_final(local_pgs)) return(1);
    }
    
    /* Set the value of C(0) */

    st_variogram_patch_C00(local_pgs,vario,idir,ngrf,rule->getRho());
    
    if (ngrf > 1 && (opt_correl != 2 || rule->getRho() != 0))
      st_varcalc_correlated_grf(local_pgs,idir);
    else
      st_varcalc_uncorrelated_grf(local_pgs,idir);

    /* Clear the geometry */

    if (flag_geometry)
      local_pgs->vorder = vario_order_manage(0,0,0,local_pgs->vorder);
  }
  return(0);
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
  int idir, ipas;

  for (idir=0; idir<vario->getDirectionNumber(); idir++)
  {
    const Dir& dir = vario->getDirs(idir);
    for (ipas=0; ipas<dir.getNPas(); ipas++)
      vario->setUtilize(idir, dir.getNPas() + ipas,1.);
  }  
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
  int idir, ipas;

  for (idir=0; idir<vario->getDirectionNumber(); idir++)
  {
    const Dir& dir = vario->getDirs(idir);
    for (ipas=0; ipas<dir.getNPas(); ipas++)
      vario->setUtilize(idir, dir.getNPas() + ipas,1.);
  }  
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
static double st_rho_search(double  rho,
                            void *user_data)
{

  double sum; 
  int  idir,ndir;   
  Local_Pgs *local_pgs;
  
  /* Initializations */

  sum = 0;
  local_pgs  = (Local_Pgs *) user_data;
  ndir = local_pgs->vario->getDirectionNumber();
  st_set_rho(rho,local_pgs);
 
  /* Evaluation of the global cost-function */

  for (idir=0; idir<ndir; idir++)
    sum += st_varcalc_correlated_grf(local_pgs,idir);

  if (debug_query("converge"))
    message("Value of the evaluating function = %lf - rho value %lf\n",
            sum,rho);

  return(sum);
}

/****************************************************************************/
/*!
**  Performing the variogram calculations (in the case of flag.rho)
**
** \return  Error return code
**
** \param[in]  db            Db structure
** \param[in]  vario         Vario structure for the GRFs to be filled
** \param[in]  rule          Lithotype Rule definition
** \param[in]  propdef       Props structure
** \param[in]  local_pgs     Local_Pgs structure
** \param[in]  ngrf          Number of GRFs
** \param[in]  opt_correl    0 full model; 1 symetrical; 2 residuals
**
*****************************************************************************/
static int st_variopgs_calcul_rho(Db    *db,
                                  Vario *vario,
                                  Rule  *rule,
                                  Props *propdef,
                                  Local_Pgs *local_pgs,
                                  int ngrf,
                                  int opt_correl)
{
  int idir;
  double testval,niter;

  /* Calculate the geometry */

  for (idir=0; idir<vario->getDirectionNumber(); idir++)
  {
    if (st_variogram_geometry_pgs_calcul(local_pgs,vario,idir)) return(1);
    st_variogram_geometry_pgs_correct(local_pgs,vario,idir);
  }
  if (st_variogram_geometry_pgs_final(local_pgs)) return(1);
    
  st_make_some_lags_inactive(vario);
  st_set_rho(golden_search(st_rho_search,(void *) local_pgs,
                           GS_TOLSTOP_RHO,-1.,1.,&testval,&niter),
             local_pgs);
  st_make_all_lags_active(vario);
  
  /* Perform the calculations with fixed rho */

  if (st_variopgs_calcul_norho(db,vario,rule,propdef,local_pgs,
                               ngrf,opt_correl,0)) return(1);
  
  /* Clean the geometry */

  local_pgs->vorder = vario_order_manage(0,0,0,local_pgs->vorder);

  return(0);
}
  
/****************************************************************************/
/*!
**  Initialize the computing time
**
** \param[in] where   Name of the calling function
**
*****************************************************************************/
static void st_timer_start(const char *where)
{
  if (! TEST_TIME) return;
  clock_start = clock();
  message("CPU time started in %s\n",where);
}

/****************************************************************************/
/*!
**  Print the computing time
**
*****************************************************************************/
static void st_timer_print(void)
{
  double cpu_process,new_clock;

  if (! TEST_TIME) return;
  new_clock = clock();
  cpu_process = 1000. * (new_clock - clock_start);
  message("CPU elapsed (Discret=%d): %lf (ms)\n",
          TEST_DISCRET,cpu_process  / CLOCKS_PER_SEC);
  clock_start = new_clock;
}

/****************************************************************************/
/*!
**  Returns the defaulted value for TEST_DISCRET parameters
**
** \return  Default value for TEST_DISCRET parameter
**
** \param[in]  rule        Lithotype Rule definition
** \param[in]  flag_rho    1 if rho has to be calculated, 0 otherwise
**
*****************************************************************************/
static int st_def_test_discret(Rule *rule,
                               int   flag_rho)
{
  if (rule->getModeRule() == RULE_SHIFT)
    return(0);
  else if (rule->getModeRule() == RULE_STD)
  {
    if (rule->getRho() <= 0. && ! flag_rho)
      return(1);
    else
      return(0);
  }
  else
    messageAbort("This rule is not expected in st_def_test_discret");
  return(0);
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
static int st_vario_pgs_check(int    flag_db,
                              int    flag_rule,
                              int    flag_varioind,
                              Db    *db,
                              Db    *dbprop,
                              Vario *vario,
                              Vario *varioind,
                              Rule  *rule)
{

  /* Experimental variogram (compulsory) */

  if (vario == (Vario *) NULL)
  {
    messerr("You must define the Input Variogram for the GRFs");
    return(1);
  }
  if (vario->getCalculType() != CALCUL_COVARIANCE    &&
      vario->getCalculType() != CALCUL_COVARIANCE_NC &&
      vario->getCalculType() != CALCUL_VARIOGRAM)
  {
    messerr("Only the Variogram is calculated here");
    return(1);
  }

  // Resize to the number of Underlying GRF
  vario->internalResize(db->getNDim(), rule->getGRFNumber(), "cov");

  /* Input Db file (optional) */

  if (flag_db != 0)
  {
    if (flag_db > 0 && db == (Db *) NULL) 
    {
      messerr("You must define the Input Db");
      return(1);
    }
    if (db != (Db *) NULL)
    {
      if (! db->isVariableNumberComparedTo(1)) return 1;
      if (db->getNDim() != vario->getDimensionNumber())
      {
        messerr("Space Dimension inconsistency between Input Db and Vario");
        return(1);
      }
    }
  }

  /* Rule (optional) */

  if (flag_rule)
  {
    if (rule == (Rule  *) NULL)
    {
      messerr("You must define the Rule");
      return(1);
    }
    if (rule->getModeRule() != RULE_STD)
    {
      messerr("This function is only programmed for standard rule");
      return(1);
    }
  }

  /* Optional Proportion File */

  if (dbprop != (Db *) NULL && dbprop->getNDim() != vario->getDimensionNumber())
  {
    messerr("Space Dimension inconsistency between Dbprop and Vario");
    return(1);
  }

  /* Indicator variogram (optional) */

  if (flag_varioind)
  {
    if (varioind == (Vario *) NULL)
    {
      messerr("You must define the Indicator Variogram (stationary case)");
      return(1);
    }
  }
  return(0);
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
static int st_variogram_pgs_nostat(Db      *db,
                                   Db      *dbprop,
                                   Vario   *vario,
                                   Rule    *rule,
                                   const    VectorDouble& propcst,
                                   int      flag_rho,
                                   int      opt_correl)
{
  Local_Pgs local_pgs;
  int    flag_correl,flag_stat,iptr_p,iptr_l,iptr_u,iptr_rl,iptr_ru;
  int    node_tot,nmax_tot,ny1,ny2,error,nfacies,ngrf;
  double prop_tot;
  Props *propdef;
  
  /* Initializations */

  error  = 1;
  ngrf   = 0;
  flag_stat = nfacies = 0;
  iptr_p = iptr_l = iptr_u = iptr_rl = iptr_ru = -1;
  propdef = (Props *) NULL;
  st_manage_pgs(0,&local_pgs,NULL,NULL,NULL,NULL,NULL,NULL,0,0,0,0,0,0);

  /* Preliminary checks */

  if (st_vario_pgs_check(1,1,0,db,dbprop,vario,NULL,rule)) goto label_end;

  /*******************/
  /* Core allocation */
  /*******************/

  ngrf = rule->getGRFNumber();
  rule->statistics(0,&node_tot,&nfacies,&nmax_tot,&ny1,&ny2,&prop_tot);

  propdef = proportion_manage(1,1,flag_stat,ngrf,0,nfacies,0,
                              db,dbprop,propcst,propdef);
  if (propdef == (Props *) NULL) goto label_end;
  flag_correl = ngrf > 1 && (opt_correl != 2 || rule->getRho() != 0);
  if (rule->particularities(db,dbprop,NULL,1,flag_stat)) goto label_end;
  proportion_rule_process(propdef,0);

  /**************************/
  /* Allocate the variables */
  /**************************/

  if (st_vario_pgs_variable(1,ngrf,nfacies,1,0,db,propdef,rule,
                            &iptr_p,&iptr_l,&iptr_u,&iptr_rl,&iptr_ru))
    goto label_end;

  /****************************/
  /* Perform the calculations */
  /****************************/
  
  /* Initialize the Local_Pgs structure */

  st_manage_pgs(1,&local_pgs,db,rule,vario,NULL,NULL,propdef,
                flag_stat,1,0,ngrf,nfacies,vario->getCalculType());
  st_define_corpgs(opt_correl,flag_rho,rule->getRho(),&local_pgs);
  st_define_trace(flag_rho,flag_correl,&local_pgs);

  /* Infer the variogram of PGS */

  if (st_vario_pgs_variable(0,ngrf,nfacies,1,0,db,propdef,rule,
                            &iptr_p,&iptr_l,&iptr_u,&iptr_rl,&iptr_ru))
    goto label_end;
  if (! flag_rho)
  {
    st_set_rho(rule->getRho(),&local_pgs);
    if (st_variopgs_calcul_norho(db,vario,rule,propdef,&local_pgs,ngrf,
                                 opt_correl,1)) goto label_end;
  }
  else
  {
    if (st_variopgs_calcul_rho(db,vario,rule,propdef,&local_pgs,ngrf,
                               opt_correl)) goto label_end;
  }

  /* Set the error return flag */

  error = 0;

label_end:
  (void) st_extract_trace(&local_pgs);
  st_manage_pgs(-1,&local_pgs,db,rule,vario,NULL,NULL,propdef,
                flag_stat,1,0,ngrf,nfacies,vario->getCalculType());
  (void) st_vario_pgs_variable(-1,ngrf,nfacies,1,0,db,propdef,rule,
                               &iptr_p,&iptr_l,&iptr_u,&iptr_rl,&iptr_ru);
  propdef = proportion_manage(-1,1,flag_stat,ngrf,0,nfacies,0,
                              db,dbprop,propcst,propdef);
  return(error);
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
                                int       *flag_ind,
                                int       *iconf,
                                double    *cov)
{
  Rule *rule;
  double covtab[4],cij00,cround;
  int    i;
  int nvar;
  CovCalcMode mode;
  
  nvar = local_pgs->model->getVariableNumber();
  rule = local_pgs->rule;

  /* Calculate the covariance for the zero distance */
  for (i=0; i<local_pgs->model->getDimensionNumber(); i++) local_pgs->d0[i] = 0.;
  model_calcul_cov(local_pgs->model,mode,1,1.,local_pgs->d0,covtab);
  cij00 = covtab[1];

  /* Calculate the covariance for the given shift */
  model_calcul_cov(local_pgs->model,mode,1,1.,local_pgs->d1,covtab);

  if (rule->getModeRule() == RULE_STD)
  {
    cov[0] = covtab[0];         /* C11(h)  */
    cov[1] = cij00;             /* C21(0)  */
    cov[2] = covtab[2];         /* C21(-h) */
    cov[3] = covtab[1];         /* C21(h)  */
    cov[4] = cij00;             /* C21(0)  */
    cov[5] = covtab[3];         /* C22(h)  */
  }
  else if (rule->getModeRule() == RULE_SHIFT)
  {
    cov[0] = covtab[0];                          /* C11(h)  */
    cov[5] = (nvar==1) ? covtab[0] : covtab[3];  /* C22(h)  */
    
    for (i=0; i<local_pgs->model->getDimensionNumber(); i++)
      local_pgs->d0[i] = rule->getShift(i);
    
    model_calcul_cov(local_pgs->model,mode,1,1.,local_pgs->d0,covtab);
    cov[1] = (nvar==1) ? covtab[0] : covtab[1];		/* C21(s)  */
    cov[4] = (nvar==1) ? covtab[0] : covtab[1];		/* C21(s)  */
    
    for (i=0; i<local_pgs->model->getDimensionNumber(); i++)
      local_pgs->d0[i] = local_pgs->d1[i] - rule->getShift(i);
    model_calcul_cov(local_pgs->model,mode,1,1.,local_pgs->d0,covtab);
    cov[2] = (nvar==1) ? covtab[0] : covtab[1];		/* C21(h-s) */
    
    for (i=0; i<local_pgs->model->getDimensionNumber(); i++)
      local_pgs->d0[i] = local_pgs->d1[i] + rule->getShift(i);
    model_calcul_cov(local_pgs->model,mode,1,1.,local_pgs->d0,covtab);
    cov[3] = (nvar==1) ? covtab[0] : covtab[1];		/* C21(h+s)  */
  }
  else
    messageAbort("This rule is not expected in st_calcul_covmatrix");

  /* Check if the two GRFs are independent */

  (*flag_ind) = 1;
  for (i=1; i<=4; i++)
    if (ABS(cov[i]) > 1.e-8) (*flag_ind) = 0;
  
  /* In TEST_DISCRET case, identify the ranks of the discretized covariance */

  if (TEST_DISCRET)
  {
    iconf[0] = ct_tableone_covrank(CTABLES,cov[0],&cround); 
    if (local_pgs->ngrf > 1)
      iconf[1] = ct_tableone_covrank(CTABLES,cov[5],&cround); 
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
                           int        flag_ind,
                           double    *low,
                           double    *up,
                           int       *iconf,
                           double    *cov)
{
  double  abseps,releps,err,proba1,proba2,proba;
  double  p1min,p1max,p2min,p2max;
  int     infin[4],maxpts,ngrf,ier;

  /* Initializations */

  releps  = 0.;
  abseps  = EPS;
  maxpts  = 8000;
  ngrf    = local_pgs->ngrf;
  proba   = TEST;
  
  /* First GRF */
  
  if (ngrf == 1 || flag_ind)
  {
    if (! TEST_DISCRET)
    {
      if (ABS(cov[0]) <= 0.)
      {
        
        /* Case where the covariance of first GRF between samples is zero */
        
        p1min  = law_cdf_gaussian(low[0]);
        p1max  = law_cdf_gaussian( up[0]);
        p2min  = law_cdf_gaussian(low[1]);
        p2max  = law_cdf_gaussian( up[1]);
        proba1 = (p1max - p1min) * (p2max - p2min);
      }
      else
      {
        infin[0] = mvndst_infin(low[0],up[0]);
        infin[1] = mvndst_infin(low[1],up[1]);
        mvndst(2,&low[0],&up[0],&infin[0],&cov[0],
               maxpts,abseps,releps,&err,&proba1,&ier);
      }
    }
    else
    {
      proba1 = ct_tableone_calculate_by_rank(CTABLES,iconf[0],low,up);
    }
    proba = proba1;
  }
  
  /* Second GRF */
  
  if (ngrf == 2)
  {
    if (flag_ind)
    {
      if (! TEST_DISCRET)
      {
        if (ABS(cov[5]) <= 0.)
        {
          
          /* Case where the covariance of second GRF between samples is zero */
          
          p1min  = law_cdf_gaussian(low[2]);
          p1max  = law_cdf_gaussian( up[2]);
          p2min  = law_cdf_gaussian(low[3]);
          p2max  = law_cdf_gaussian( up[3]);
          proba2 = (p1max - p1min) * (p2max - p2min);
        }
        else
        {
          infin[2] = mvndst_infin(low[2],up[2]);
          infin[3] = mvndst_infin(low[3],up[3]);
          mvndst(2,&low[2],&up[2],&infin[2],&cov[5],
                 maxpts,abseps,releps,&err,&proba2,&ier);
        }
      }
      else
      {
        proba2 = ct_tableone_calculate_by_rank(CTABLES,
                                               iconf[1],&low[2],&up[2]);
      }
      proba *= proba2;
    }
    else
    {
      infin[0] = mvndst_infin(low[0],up[0]);
      infin[1] = mvndst_infin(low[1],up[1]);
      infin[2] = mvndst_infin(low[2],up[2]);
      infin[3] = mvndst_infin(low[3],up[3]);
      mvndst(4,low,up,infin,cov,maxpts,abseps,releps,&err,&proba,&ier);
    }
  }
 
  return(proba);
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
    low[0]  = STAT_THRESH(ifac1,0,0);
    up[0]   = STAT_THRESH(ifac1,0,1);
    low[1]  = STAT_THRESH(ifac2,0,0);
    up[1]   = STAT_THRESH(ifac2,0,1);
    if (local_pgs->ngrf > 1)
    {
      low[2]  = STAT_THRESH(ifac1,1,0);
      up[2]   = STAT_THRESH(ifac1,1,1);
      low[3]  = STAT_THRESH(ifac2,1,0);
      up[3]   = STAT_THRESH(ifac2,1,1);
    }
  }
  else
  {
    ploc[0] = local_pgs->db->getProportion(iech1,ifac1);
    ploc[1] = local_pgs->db->getProportion(iech2,ifac2);

    if (! TEST_DISCRET)
    {
      low[0] = local_pgs->db->getLowerBound(iech1,ifac1);
      up[0]  = local_pgs->db->getUpperBound(iech1,ifac1);
      low[1] = local_pgs->db->getLowerBound(iech2,ifac2);
      up[1]  = local_pgs->db->getUpperBound(iech2,ifac2);
      if (local_pgs->ngrf > 1) 
      {
        low[2] = local_pgs->db->getLowerBound(iech1,nfacies+ifac1);
        up[2]  = local_pgs->db->getUpperBound(iech1,nfacies+ifac1);
        low[3] = local_pgs->db->getLowerBound(iech2,nfacies+ifac2);
        up[3]  = local_pgs->db->getUpperBound(iech2,nfacies+ifac2);
      }
    }
    else
    {
      low[0] = local_pgs->db->getLowerInterval(iech1,ifac1);
      up[0]  = local_pgs->db->getUpperInterval(iech1,ifac1);
      low[1] = local_pgs->db->getLowerInterval(iech2,ifac2);
      up[1]  = local_pgs->db->getUpperInterval(iech2,ifac2);
      if (local_pgs->ngrf > 1) 
      {
        low[2] = local_pgs->db->getLowerInterval(iech1,nfacies+ifac1);
        up[2]  = local_pgs->db->getUpperInterval(iech1,nfacies+ifac1);
        low[3] = local_pgs->db->getLowerInterval(iech2,nfacies+ifac2);
        up[3]  = local_pgs->db->getUpperInterval(iech2,nfacies+ifac2);
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
                           int        flag_ind,
                           int        iech1,
                           int        iech2,
                           int        ifac1,
                           int        ifac2,
                           int       *iconf,
                           double    *cov)
{
  double value,g1,g2,ploc[2],low[4],up[4];

  if (local_pgs->covtype == CALCUL_VARIOGRAM)
  {
    if (ifac1 == ifac2)
    {
      st_define_bounds(local_pgs,iech1,iech2,ifac1,ifac2,low,up,ploc);
      g1 = st_get_proba(local_pgs,flag_ind,low,up,iconf,cov);
      value = (ploc[0] + ploc[1]) * 0.5 - g1;
    }
    else
    {
      st_define_bounds(local_pgs,iech1,iech2,ifac1,ifac2,low,up,ploc);
      g1 = st_get_proba(local_pgs,flag_ind,low,up,iconf,cov);
      st_define_bounds(local_pgs,iech2,iech1,ifac1,ifac2,low,up,ploc);
      g2 = st_get_proba(local_pgs,flag_ind,low,up,iconf,cov);
      value = -0.5 * (g1 + g2);
    }
  }
  else
  {
    st_define_bounds(local_pgs,iech1,iech2,ifac1,ifac2,low,up,ploc);
    value = st_get_proba(local_pgs,flag_ind,low,up,iconf,cov);
  }
  return(value);
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
  double  dist,cov[6];
  int     ipas,ifac,jfac,nfacies,ipair,iech,jech,i,idir,flag_ind,iconf[2];
  Vario  *vario;
 
  /* Initializations */

  nfacies = local_pgs->nfacies;
  vario   = local_pgs->vario;

  /* Loop on the directions */

  for (idir=0; idir<vario->getDirectionNumber(); idir++)
  {
    const Dir& dir = vario->getDirs(idir);

    /* Establish the geometry */

    if (st_variogram_geometry_pgs_calcul(local_pgs,vario,idir)) return(1);
    if (st_variogram_geometry_pgs_final(local_pgs)) return(1);

    /* Clean the distance and variogram */

    for (i=0; i<dir.getSize(); i++)
    {
      vario->setSw(idir, i,0.);
      vario->setHh(idir, i,0.);
      vario->setGg(idir, i,0.);
    }

    /* Loop on the lags */

    for (ipas=0; ipas<dir.getNPas(); ipas++)
    {
      vario_order_get_bounds(local_pgs->vorder,idir,ipas,
                             &local_pgs->ifirst,&local_pgs->ilast);
      if (local_pgs->ifirst >= local_pgs->ilast) continue;
      
      /* Loop on the pairs of the lag */

      for (ipair=local_pgs->ifirst; ipair<local_pgs->ilast; ipair++)
      {
        vario_order_get_indices(local_pgs->vorder,ipair,&iech,&jech,&dist);
        
        /* Calculate the distance vector */
        
        dist = distance_intra(local_pgs->db,iech,jech,local_pgs->d1.data());
        st_calcul_covmatrix(local_pgs,&flag_ind,iconf,cov);
        
        /* Loops on the facies */
        
        for (ifac=0; ifac<nfacies; ifac++)
          for (jfac=0; jfac<=ifac; jfac++)
          {
            if (local_pgs->vario->getFlagAsym())
            {
              i = dir.getAddress(ifac,jfac,ipas,false,1);
              vario->setSw(idir, i, dir.getSw(i) + 1.);
              vario->setHh(idir, i, dir.getHh(i) + dist);
              vario->setGg(idir, i, dir.getGg(i) +
                st_get_value(local_pgs,flag_ind,iech,jech,ifac,jfac,iconf,cov));
              i = dir.getAddress(ifac,jfac,ipas,false,-1);
              vario->setSw(idir, i, dir.getSw(i) + 1.);
              vario->setHh(idir, i, dir.getHh(i) - dist);
              vario->setGg(idir, i, dir.getGg(i) +
                st_get_value(local_pgs,flag_ind,jech,iech,ifac,jfac,iconf,cov));
            }        
            else
            {
              i = dir.getAddress(ifac,jfac,ipas,false,0);
              vario->setSw(idir, i, dir.getSw(i) + 1.);
              vario->setHh(idir, i, dir.getHh(i) + dist);
              vario->setGg(idir, i, dir.getGg(i) +
                st_get_value(local_pgs,flag_ind,iech,jech,ifac,jfac,iconf,cov));
            }
          }
      }
    }
    
    /* Clear the geometry */
    
    local_pgs->vorder = vario_order_manage(0,0,0,local_pgs->vorder);

    /* Scale the variogram */

    variogram_scale(vario,idir);
  }
  return(0);
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
  int    ipas,jpas,ifac,jfac,nfacies,i,idir,flag_ind,iconf[2];
  Vario *vario;
 
  /* Initializations */

  nfacies = local_pgs->nfacies;
  vario   = local_pgs->vario;
  
  /* Loop on the directions */
  
  for (idir=0; idir<vario->getDirectionNumber(); idir++)
  {
    const Dir& dir = vario->getDirs(idir);
    
    /* Loop on the lags */
    
    for (ipas=0; ipas<dir.getNPas(); ipas++)
    {
      
      /* Calculate the distance vector */
      
      for (i=0; i<vario->getDimensionNumber(); i++)
      {
        jpas = (local_pgs->vario->getFlagAsym()) ? dir.getNPas() + ipas : ipas;
        local_pgs->d1[i] = dir.getHh(jpas) * dir.getCodir(i);
      }
      st_calcul_covmatrix(local_pgs,&flag_ind,iconf,cov);
      
      /* Loops on the facies */
      
      for (ifac=0; ifac<nfacies; ifac++)
        for (jfac=0; jfac<=ifac; jfac++)
        {
          if (local_pgs->vario->getFlagAsym())
          {
            i = dir.getAddress(ifac,jfac,ipas,false,1);
            vario->setGg(idir, i,
              st_get_value(local_pgs,flag_ind,0,0,ifac,jfac,iconf,cov));
            i = dir.getAddress(ifac,jfac,ipas,false,-1);
            vario->setGg(idir, i,
              st_get_value(local_pgs,flag_ind,0,0,jfac,ifac,iconf,cov));
          }        
          else
          {
            i = dir.getAddress(ifac,jfac,ipas,false,0);
            vario->setGg(idir, i,
              st_get_value(local_pgs,flag_ind,0,0,ifac,jfac,iconf,cov));
          }        
        }
    }
  }
  return(0);
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
  int ivar,jvar,idir,iad,nfacies;
  double pivar,pjvar;
  Vario *vario;
  
  /* Initializations */

  vario   = local_pgs->vario;
  nfacies = local_pgs->nfacies;

  /* Evaluate the theoretical variances */

  for (ivar=0; ivar<nfacies; ivar++)
    for (jvar=0; jvar<nfacies; jvar++)
    {
      pivar = local_pgs->propdef->propfix[ivar];
      pjvar = local_pgs->propdef->propfix[jvar];
      if (ivar == jvar)
        vario->setVars(ivar,jvar,pivar * (1. - pivar));
      else
        vario->setVars(ivar,jvar,-pivar * pjvar);
      if(!vario->getFlagAsym()) continue;
      
      for (idir=0; idir<vario->getDirectionNumber(); idir++)
      {
        const Dir& dir = vario->getDirs(idir);
        iad = dir.getAddress(ivar,jvar,0,false,0);
        vario->setSw(idir,iad,1);
        vario->setHh(idir,iad,0);
        
        switch (local_pgs->covtype)
        {   
          case CALCUL_VARIOGRAM:
            break;
            
          case CALCUL_COVARIANCE:
            vario->setGg(idir,iad,vario->getVars(ivar,jvar));
            break;
            
          case CALCUL_COVARIANCE_NC:
            vario->setGg(idir,iad,(ivar == jvar) ? pivar : 0.);
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
  int     i,nfacies,ivar,jvar,number,error,iech,idir,iad;
  double *mean,*covs,p1,p2;
  Vario *vario;
  Db  *dbin;
  
  /* Initializations */

  error   = 1;
  dbin    = local_pgs->db;
  vario   = local_pgs->vario;
  nfacies = local_pgs->nfacies;
  mean    = covs = (double *) NULL;

  /* Core allocation */

  mean = (double *) mem_alloc(nfacies * sizeof(double), 0);
  if (mean == (double *) NULL) goto label_end;
  for (i=0; i<nfacies; i++) mean[i] = 0.;
  covs = (double *) mem_alloc(nfacies * nfacies * sizeof(double), 0);
  if (covs == (double *) NULL) goto label_end;
  for (i=0; i<nfacies * nfacies; i++) covs[i] = 0.;

  /* Loop on the samples */
  
  for (iech=number=0; iech<get_NECH(dbin); iech++)
  {
    if (! dbin->isActive(iech)) continue;
    
    /* Loop on the variables */

    for (ivar=0; ivar<nfacies; ivar++)
    {
      p1 = dbin->getProportion(iech,ivar);
      mean[ivar] += p1;

      for (jvar=0; jvar<nfacies; jvar++)
      {        
        p2 = dbin->getProportion(iech,jvar);
        COVS(ivar,jvar) += p1 * p2;
      }        
    }
    number++;
  }

  /* Normalization */

  for (ivar=0; ivar<nfacies; ivar++)
  {
    mean[ivar] /= number;
    for (jvar=0; jvar<nfacies; jvar++)
      COVS(ivar,jvar) /= number;
  }

  /* Evaluate the variances */

  for (ivar=0; ivar<nfacies; ivar++)
    for (jvar=0; jvar<nfacies; jvar++)
    {
      if (ivar == jvar)
        vario->setVars(ivar,jvar,mean[ivar] - COVS(ivar,ivar));
      else
        vario->setVars(ivar,jvar,-COVS(ivar,jvar));

      if(! vario->getFlagAsym()) continue;
      for (idir=0; idir<vario->getDirectionNumber(); idir++)
      {
        const Dir& dir = vario->getDirs(idir);
        iad = dir.getAddress(ivar,jvar,0,false,0);
        vario->setSw(idir,iad,get_NECH(dbin));
        vario->setHh(idir,iad,0);
        
        switch (local_pgs->covtype)
        {
          case CALCUL_VARIOGRAM:
            break;
            
          case CALCUL_COVARIANCE:
            vario->setGg(idir,iad,vario->getVars(ivar,jvar));
            break;
            
          case CALCUL_COVARIANCE_NC:
            vario->setGg(idir,iad,(ivar == jvar) ? mean[ivar] : 0.);
              break;
              
        } 
      }
    }
  
  /* Set the error return code */
  
  error = 0;

label_end:
  mean = (double *) mem_free((char *) mean);
  covs = (double *) mem_free((char *) covs);
  return(error);
}

/****************************************************************************/
/*!
**  Evaluate the experimental variogram of indicators in Plurigaussian case
**
** \return  Error return code
**
** \param[in]  db         Db descriptor
** \param[in]  dbprop     Db descriptor for the grid of proportions
** \param[in]  vario      Vario structure
** \param[in]  rule       Lithotype Rule definition for first point
** \param[in]  propcst    Array giving the constant proportions
** \param[in]  flag_stat  1 for stationary; 0 otherwise
** \param[in]  model1     First Model structure
** \param[in]  model2     Second Model structure (optional)
**
** \remark  At this stage, the number of variables is equal to the number
** \remark  of indicators.
** \remark  However, the models(s) are defined for a single variable
**
*****************************************************************************/
GEOSLIB_API int model_pgs(Db     *db,
                          Db     *dbprop,
                          Vario  *vario,
                          Rule   *rule,
                          const   VectorDouble& propcst,
                          int     flag_stat,
                          Model  *model1,
                          Model  *model2)
{
  Local_Pgs local_pgs;
  int     error,nfacies,ngrf,node_tot,nmax_tot;
  int     ny1,ny2,iptr_l,iptr_u,iptr_p,iptr_rl,iptr_ru;
  double  prop_tot;
  Props  *propdef;
  Model  *new_model;
  
  /*******************/
  /* Initializations */
  /*******************/

  st_timer_start("model_pgs");
  TEST_DISCRET = (int) get_keypone("TEST_DISCRET",
                                   st_def_test_discret(rule,0));
  TEST_TIME    = (int) get_keypone("TEST_TIME",0);
  error  = 1;
  ngrf   = nfacies = 0;
  iptr_l = iptr_u = iptr_p = iptr_rl = iptr_ru = -1;
  new_model = (Model *) NULL;
  propdef   = (Props *) NULL;
  st_manage_pgs(0,&local_pgs,NULL,NULL,NULL,NULL,NULL,NULL,0,0,0,0,0,0);

  /* Preliminary checks */

  if (st_vario_pgs_check(-1,1,0,db,dbprop,vario,NULL,rule)) goto label_end;
  
  /* Merge the models */

  new_model = model_rule_combine(model1,model2,rule);
  if (new_model == (Model *) NULL) goto label_end;

  ngrf = rule->getGRFNumber();
  rule->statistics(0,&node_tot,&nfacies,&nmax_tot,&ny1,&ny2,&prop_tot);
  if(rule->getModeRule() == RULE_SHIFT) ngrf++;

  if (nfacies != vario->getVariableNumber())
  {
    messerr("Inconsistency between the Variogram and the Rule");
    messerr("- Number of variables in the Variogram = %d",vario->getVariableNumber());
    messerr("- Number of facies in the Rule = %d",nfacies);
    return(1);
  }
  if (new_model == (Model *) NULL)
  {
    messerr("The Model(s) must be defined");
    return(1);
  }
  if (new_model->getVariableNumber() != ngrf)
  {
    messerr("The number of GRF is not equal to the number of variables");
    messerr("defined in the combined Model");
    return(1);
  }
  if (! flag_stat)
  {
    if (db == (Db *) NULL)
    {
      messerr("You must define the Input Db");
      return(1);
    }
    if (db->getNDim() != vario->getDimensionNumber())
    {
      messerr("Inconsistent parameters:");
      messerr("Input DB : NDIM=%d",db->getNDim());
      messerr("Variogram: NDIM=%d",vario->getDimensionNumber());
      return(1);
    }
    if (dbprop != (Db *) NULL && dbprop->getNDim() != vario->getDimensionNumber())
    {
      messerr("Space Dimension inconsistency between Dbprop and Vario");
      return(1);
    }
    if (new_model->getDimensionNumber() != db->getNDim())
    {
      messerr("The Space Dimension of the Db structure (%d)",db->getNDim());
      messerr("Does not correspond to the Space Dimension of the model (%d)",
              new_model->getDimensionNumber());
      goto label_end;
    }
  }

  /*******************/
  /* Core allocation */
  /*******************/
  
  propdef = proportion_manage(1,1,flag_stat,ngrf,0,nfacies,0,
                              db,dbprop,propcst,propdef);
  if (propdef == (Props *) NULL) goto label_end;

  if (rule->particularities(db,dbprop,new_model,0,flag_stat)) goto label_end;

  proportion_rule_process(propdef,0);

  /* Pre-calculation of integrals: Define the structure */

  if (TEST_DISCRET)
    CTABLES = ct_tables_manage(1,0,1,2,200,100,-1.,1.,NULL);

  st_manage_pgs(1,&local_pgs,db,rule,vario,NULL,new_model,propdef,
                flag_stat,0,1,ngrf,nfacies,vario->getCalculType());

  if (! flag_stat)
  {
    if (st_vario_pgs_variable(1,ngrf,nfacies,0,1,db,propdef,rule,
                              &iptr_p,&iptr_l,&iptr_u,&iptr_rl,&iptr_ru))
      goto label_end;
  }
  else
  {
    if (st_calculate_thresh_stat(&local_pgs)) goto label_end;
  }

  /* Calculate the variance matrix and the variogram */

  if (flag_stat)
  {
    if (st_vario_indic_model_stat(&local_pgs)) goto label_end;
    st_update_variance_stat(&local_pgs);
  }
  else
  {
    if (st_vario_pgs_variable(0,ngrf,nfacies,0,1,db,propdef,rule,
                              &iptr_p,&iptr_l,&iptr_u,&iptr_rl,&iptr_ru))
      goto label_end;
    if (st_vario_indic_model_nostat(&local_pgs)) goto label_end;
    if (st_update_variance_nostat(&local_pgs)) goto label_end;
  }

  /* Set the error return code */

  error = 0;

label_end:
  if (TEST_DISCRET)
    CTABLES = ct_tables_manage(-1,0,1,2,200,100,-1.,1.,CTABLES);
  st_manage_pgs(-1,&local_pgs,db,rule,vario,NULL,new_model,propdef,
                flag_stat,0,1,ngrf,nfacies,vario->getCalculType());
  new_model = model_free(new_model);
  (void) st_vario_pgs_variable(-1,ngrf,nfacies,0,1,db,propdef,rule,
                               &iptr_p,&iptr_l,&iptr_u,&iptr_rl,&iptr_ru);
  propdef = proportion_manage(-1,1,flag_stat,ngrf,0,nfacies,0,
                              db,dbprop,propcst,propdef);
  st_timer_print();
  return(error);
}

static void st_copy_swhh(Local_Pgs* local_pgs)
{
  Vario* vario    = local_pgs->vario;
  Vario* varioind = local_pgs->varioind;
  int ngrf = local_pgs->ngrf;

  for (int idir = 0; idir < (int) vario->getDirectionNumber(); idir++)
  {
    const Dir& dir1 = vario->getDirs(idir);
    const Dir& dir2 = varioind->getDirs(idir);
    for (int ipas = 0; ipas < (int) dir1.getLagTotalNumber(); ipas++)
    {
      int iad2 = dir2.getAddress(0,0,ipas,true,0);
      for (int igrf = 0; igrf < ngrf; igrf++)
        for (int jgrf = 0; jgrf < ngrf; jgrf++)
        {
          int iad1 = dir1.getAddress(igrf,jgrf,ipas,true,0);
          vario->setSw(idir,iad1, dir2.getSw(iad2));
          vario->setHh(idir,iad1, dir2.getHh(iad2));
        }
    }
  }
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
static int st_variogram_pgs_stat(Db     *db,
                                 Vario  *vario,
                                 Vario  *varioind,
                                 Rule   *rule,
                                 const   VectorDouble& propcst)
{
  Local_Pgs local_pgs;
  int    node_tot,nmax_tot,ny1,ny2,error,nfacies,ngrf,flag_stat;
  double prop_tot;
  Props *propdef;
  
  /* Initializations */

  error  = 1;
  ngrf   = nfacies = 0;
  flag_stat = 1;
  propdef = (Props *) NULL;
  st_manage_pgs(0,&local_pgs,NULL,NULL,NULL,NULL,NULL,NULL,0,0,0,0,0,0);

  /* Preliminary checks */

  if (st_vario_pgs_check(0,1,1,db,NULL,vario,varioind,rule)) goto label_end;

  /*******************/
  /* Core allocation */
  /*******************/

  ngrf = rule->getGRFNumber();
  rule->statistics(0,&node_tot,&nfacies,&nmax_tot,&ny1,&ny2,&prop_tot);
  propdef = proportion_manage(1,1,flag_stat,ngrf,0,nfacies,0,
                              NULL,NULL,propcst,propdef);
  if (propdef == (Props *) NULL) goto label_end;
  if (rule->particularities(NULL,NULL,NULL,1,flag_stat)) goto label_end;
  proportion_rule_process(propdef,0);

  /****************************/
  /* Perform the calculations */
  /****************************/
  
  /* Initialize the Local_Pgs structure */

  st_manage_pgs(1,&local_pgs,NULL,rule,vario,varioind,NULL,propdef,
                flag_stat,1,0,ngrf,nfacies,vario->getCalculType());
  st_define_corpgs(0,0,rule->getRho(),&local_pgs);
  st_define_trace(0,0,&local_pgs);
  st_set_rho(0.,&local_pgs);
  if (st_calculate_thresh_stat(&local_pgs)) goto label_end;

  // Copy the count and distance information from indicator variogram into calculation one

  st_copy_swhh(&local_pgs);

  /* Infer the variogram of PGS */

  st_varcalc_from_vario_stat(vario,rule,propdef,&local_pgs,ngrf);

  /* Bring the calculations back from a covariance to a variogram */


  /* Set the error return flag */

  error = 0;

label_end:
  (void) st_extract_trace(&local_pgs);
  st_manage_pgs(-1,&local_pgs,NULL,rule,vario,varioind,NULL,propdef,
                flag_stat,1,0,ngrf,nfacies,vario->getCalculType());
  propdef = proportion_manage(-1,1,1,ngrf,0,nfacies,0,
                              NULL,NULL,propcst,propdef);
  return(error);
}

/****************************************************************************/
/*!
**  Calculate the gaussian variograms 
**
** \return  Error return code
**
** \param[in]  db           Db structure
** \param[in]  vario        Vario structure for the GRFs to be filled
** \param[in]  rule         Lithotype Rule definition
** \param[in]  propcst      Array of proportions for the facies
** \param[in]  dbprop       Db Grid used for proportions (non-stationary)
** \param[in]  flag_stat    1 for stationary and 0 otherwise
** \param[in]  flag_rho     1 if the correlation coefficient must be regressed
** \param[in]  opt_correl   0 full model; 1 symmetrical; 2 residuals
**
** \remarks This is simply a routine dispatching between the stationary function
** \remarks and the non-stationary one
**
*****************************************************************************/
GEOSLIB_API int variogram_pgs(Db     *db,
                              Vario*  vario,
                              Rule*   rule,
                              const   VectorDouble& propcst,
                              Db     *dbprop,
                              int     flag_stat,
                              int     flag_rho,
                              int     opt_correl)
{
  Vario* varioind = nullptr;
  VectorDouble props;
  int error;

  /* Initializations */

  st_timer_start("variogram_pgs");
  TEST_DISCRET = (int) get_keypone("TEST_DISCRET",
                                   st_def_test_discret(rule,0));
  TEST_TIME    = (int) get_keypone("TEST_TIME",0);

  // Preliminary checks

  if (vario == NULL)
  {
    messerr("The output Variogram must be provided (empty)");
    return 1;
  }
  if (vario->getDirectionNumber() <= 0)
  {
    messerr("The variogram must contain at least one calculation Direction");
    return 1;
  }
  if (db->getVariableNumber() != 1)
  {
    messerr("The number of variables (%d) must be equal to 1",db->getVariableNumber());
    return 1;
  }
  vario->setCalculName("covnc");
  int iatt = db->getAttribute(LOC_Z,0);
  int nclass = rule->getFaciesNumber();
  if (nclass <= 0)
  {
    messerr("No Facies class have been found");
    return 1;
  }

  // In Stationary case, create the variogram of indicators to speed up calculations

  if (flag_stat)
  {
    if (! propcst.empty())
    {
      if ((int) propcst.size() != nclass)
      {
        messerr("Number of proportions in 'propcst' (%d) should match Number of Facies in 'rule' (%d)",
                (int) propcst.size(),rule->getFaciesNumber());
        return 1;
      }
      props.resize(nclass);
      props = propcst;
    }
    else
    {
      // Calculate the number of Facies in 'Db'
      props = dbStatisticsFacies(db);
      if ((int) props.size() != nclass)
      {
        messerr("Number of Facies in 'db' (%d) should match Number of facies in 'rule' (%d)",
                (int) props.size(), rule->getFaciesNumber());
        return 1;
      }
    }

    // Translate the 'Facies' into 'categories'
    Limits limits = Limits(nclass);
    if (limits.toIndicator(db,iatt))
    {
      messerr("Problem when translating Facies into Categories");
      return 1;
    }

    // Calculate the variogram of Indicators
    varioind = (Vario*) vario->clone();
    VectorDouble vars(nclass * nclass,0.);
    int ecr = 0;
    for (int iclass = 0; iclass < nclass; iclass++)
      for (int jclass = 0 ; jclass < nclass; jclass++)
        vars[ecr++] = (iclass != jclass) ? 0. : props[iclass] * (1. - props[iclass]);
    if (varioind->compute(db,"covnc",props,vars))
    {
      messerr("Error when calculating the Variogram of Indicators");
      return 1;
    }
  }

  /* Pre-calculation of integrals: Define the structure */

  if (TEST_DISCRET)
    CTABLES = ct_tables_manage(1,0,1,2,200,100,-1.,1.,NULL);

  /* Perform the calculations */

  if (flag_stat)
    error = st_variogram_pgs_stat(db,vario,varioind,rule,props);
  else
    error = st_variogram_pgs_nostat(db,dbprop,vario,rule,props,flag_rho,opt_correl);

  /* Final operations */

  if (TEST_DISCRET)
    CTABLES = ct_tables_manage(-1,0,1,2,200,100,-1.,1.,CTABLES);
  st_timer_print();
  return(error);
}
  
/****************************************************************************/
/*!
**  Find the optimal Truncation Scheme from Variopgs score
**
** \return  The newly created Rule structure
**
** \param[in]  db           Db structure
** \param[in]  dbprop       Db Grid used for proportions (non-stationary)
** \param[in]  vario        Vario structure for the GRFs to be filled
** \param[in]  varioind     Indicator Vario structure 
** \param[in]  propcst      Array of proportions for the facies
** \param[in]  ncolor       Number of different facies
** \param[in]  ngrf         Number of underlying GRFs (1 or 2)
** \param[in]  flag_stat    1 for stationary and 0 otherwise
** \param[in]  verbose      Verbosity flag
**
*****************************************************************************/
GEOSLIB_API Rule *rule_auto(Db     *db,
                            Db     *dbprop,
                            Vario  *vario,
                            Vario  *varioind,
                            const   VectorDouble& propcst,
                            int     ncolor,
                            int     ngrf,
                            int     flag_stat,
                            int     verbose)
{
  int    *facies,*fcmp,*fgrf,*string,error,nscore,r_opt,nfacies;
  int     iptr_p,iptr_l,iptr_u,iptr_rl,iptr_ru;
  int    *rules,flag_rho,flag_correl,opt_correl;
  Rule   *rule;
  Relem  *Pile_Relem;
  double *scores;
  Props  *propdef;
  Local_Pgs local_pgs;

  /* Initializations */

  error        = 1;
  rule         = (Rule *) NULL;
  scores       = (double *) NULL;
  facies       = fcmp = fgrf = string = (int *) NULL;
  Pile_Relem   = (Relem *) NULL;
  propdef      = (Props *) NULL;
  st_timer_start("rule_auto");
  TEST_DISCRET = 1;
  TEST_TIME    = (int) get_keypone("TEST_TIME",0);
  NCOLOR       = nfacies = ncolor;
  NGRF         = ngrf;
  NRULE        = 2 * NCOLOR - 1;
  BASE         = 2 * NGRF;
  flag_rho     = 0;
  flag_correl  = 0;
  opt_correl   = 0;
  iptr_p       = iptr_l = iptr_u = iptr_rl = iptr_ru = -1;

  /* Core allocation */

  facies = (int *) mem_alloc(sizeof(int) * NCOLOR,1);
  for (int i=0; i<NCOLOR; i++) facies[i] = i+1;

  /* Preliminary tasks (as in variogram.pgs) */

  st_manage_pgs(0,&local_pgs,NULL,NULL,NULL,NULL,NULL,NULL,0,0,0,0,0,0);
  if (st_vario_pgs_check(0,0,flag_stat,db,NULL,vario,varioind,NULL))
    goto label_end;
  propdef = proportion_manage(1,1,flag_stat,ngrf,0,nfacies,0,
                              db,dbprop,propcst,propdef);
  if (propdef == (Props *) NULL) goto label_end;
  proportion_rule_process(propdef,0);
  
  /* Pre-calculation of integrals: Define the structure */

  CTABLES = ct_tables_manage(1,verbose,1,2,200,100,-1.,1.,NULL);

  /* Allocation */

  st_manage_pgs(1,&local_pgs,db,NULL,vario,varioind,NULL,propdef,
                flag_stat,1,0,ngrf,nfacies,vario->getCalculType());

  if (flag_stat)
  {
    st_define_corpgs(0,0,0,&local_pgs);
    st_define_trace(0,0,&local_pgs);
  }
  else
  {
    st_define_corpgs(opt_correl,flag_rho,0.,&local_pgs);
    st_define_trace(flag_rho,flag_correl,&local_pgs);

    // Prepare the geometry 
    
    for (int idir=0; idir<vario->getDirectionNumber(); idir++)
    {
      local_pgs.idircur = idir;
      if (st_variogram_geometry_pgs_calcul(&local_pgs,vario,idir)) 
        goto label_end;
      st_variogram_geometry_pgs_correct(&local_pgs,vario,idir);
    }
    if (st_variogram_geometry_pgs_final(&local_pgs)) goto label_end;

    // The thresholds are added lately in order to allow calculation of 
    // geometry (without checking the threshold interval (not defined yet)
    if (st_vario_pgs_variable(1,ngrf,nfacies,1,0,db,propdef,NULL,
                              &iptr_p,&iptr_l,&iptr_u,&iptr_rl,&iptr_ru))
      goto label_end;
  }

  /* Elaborate the whole tree */

  Pile_Relem = st_relem_alloc(NULL);
  st_relem_define(Pile_Relem,NCOLOR,facies,ITEST,NULL);
  st_relem_subdivide(Pile_Relem,1,1);
  st_relem_explore(Pile_Relem,0);

  // Evaluate all possibilities

  fcmp = (int *) mem_alloc(sizeof(int) * NCOLOR,1);
  fgrf = (int *) mem_alloc(sizeof(int) * (1+NGRF),1);
  scores = st_relem_evaluate(Pile_Relem,verbose,fgrf,fcmp,
                             &local_pgs,&nscore,&r_opt);
  fcmp = (int *) mem_free((char *) fcmp);
  fgrf = (int *) mem_free((char *) fgrf);

  /* Get the resulting optimal Rule */

  st_rule_print(r_opt,NRULE,Pile_Relem->rules,Pile_Relem->fipos,0,-1,-1,TEST);
  rules = Pile_Relem->rules;
  string = &RULES(r_opt,0);
  rule = st_rule_encode(string);

  /* Clean the geometry (non-stationary case) */

  if (! flag_stat)
    local_pgs.vorder = vario_order_manage(0,0,0,local_pgs.vorder);

  /* Set the error return code */

  error = 0;

label_end:
  Pile_Relem = st_relem_free(Pile_Relem);
  facies = (int    *) mem_free((char *) facies);
  scores = (double *) mem_free((char *) scores);
  CTABLES = ct_tables_manage(-1,verbose,1,2,200,100,-1.,1.,CTABLES);
  st_manage_pgs(-1,&local_pgs,db,NULL,vario,varioind,NULL,propdef,
                flag_stat,1,0,ngrf,nfacies,vario->getCalculType());
  propdef = proportion_manage(-1,1,flag_stat,ngrf,0,nfacies,0,
                              db,dbprop,propcst,propdef);
  if (error) rule = rule_free(rule);
  st_timer_print();
  return(rule);
}
