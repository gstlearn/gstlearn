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
#include "LithoRule/Rule.hpp"
#include "Basic/Utilities.hpp"

/*! \cond */
#define INDLOC(ifac1,ifac2)  ((ifac2)*propdef->nfac[0]+(ifac1))
#define PROPFIX(ifac1,ifac2) (propdef->propfix[INDLOC(ifac1,ifac2)])
#define PROPWRK(ifac1,ifac2) (propdef->propwrk[INDLOC(ifac1,ifac2)])
/*! \endcond */

/****************************************************************************/
/*!
**  Free a Rule structure
**
** \return  Pointer to the newly freed Rule structure
**
** \param[in]  rule Rule structure to be freed
**
*****************************************************************************/
GEOSLIB_API Rule *rule_free(Rule *rule)

{
  if (rule == (Rule *) NULL) return(rule);
  delete rule;
  rule = nullptr;
  return(rule);
}

/****************************************************************************/
/*!
**  Print the Gaussian to Facies translation
**
** \param[in]  string     Type of set on which translation is carried on
** \param[in]  isimu      Simulation rank
** \param[in]  iech       Sample rank
** \param[in]  y          Gaussian values
** \param[in]  facies     Translated facies
**
*****************************************************************************/
static void st_print_gaus2fac(const char *string,
                              int         isimu,
                              int         iech,
                              double     *y,
                              double      facies)
{
  if (iech == 0)
    mestitle(1,"%s: Gaussian -> Facies (Simulation=%d)",string,isimu+1);
  
  message("Sample (%6d) - ",iech+1);
  if (! FFFF(y[0])) message(" Y1 = %8.5lf",y[0]);
  if (! FFFF(y[1])) message(" Y2 = %8.5lf",y[1]);
  message(" -> Facies = %d\n",(int) facies);
}

/****************************************************************************/
/*!
**  Combine the underlying GRF into a facies value at data points (Shadow)
**
** \return  Error return code
**
** \param[in]  propdef    Props structure
** \param[in]  dbin       Db input structure
** \param[in]  dbout      Db output structure
** \param[in]  rule       Rule structure
** \param[in]  flag_used  1 if the gaussian is used; 0 otherwise
** \param[in]  ipgs       Rank of the PGS
** \param[in]  isimu      Rank of the simulation
** \param[in]  nbsimu     Number of simulations
**
** \remark Attributes LOC_GAUSFAC are mandatory
** \remark Attributes LOC_FACIES are mandatory
**
*****************************************************************************/
GEOSLIB_API int rule_gaus2fac_data_shadow(Props  *propdef,
                                          Db     *dbin,
                                          Db     *dbout,
                                          Rule   *rule,
                                          int    *flag_used,
                                          int     ipgs,
                                          int     isimu,
                                          int     nbsimu)
{
  int    iech,error,igrf,icase;
  double y[2],facies,t1min,t1max,t2min,t2max,sh_dsup,sh_down;

  /* Initializations */

  error = 1;
  check_mandatory_attribute("rule_gaus2fac_data_shadow",dbin,LOC_GAUSFAC);

  /* Processing the translation */

  for (iech=0; iech<get_NECH(dbin); iech++)
  {
    if (! dbin->isActive(iech)) continue;
    
    /* Initializations */
    
    facies = TEST;
    for (igrf=0; igrf<2; igrf++) y[igrf] = TEST;
    
    if (rule_thresh_define_shadow(propdef,dbin,rule,ITEST,
                                  iech,isimu,nbsimu,1,
                                  &t1min,&t1max,&t2min,&t2max,
                                  &sh_dsup,&sh_down)) goto label_end;
    
    for (igrf=0; igrf<2; igrf++)
    {
      icase = get_rank_from_propdef(propdef,ipgs,igrf);
      y[igrf] = (flag_used[igrf]) ?
        dbin->getSimvar(LOC_GAUSFAC,iech,isimu,0,icase,nbsimu,1) : 0.;
    }
    facies = rule->getFaciesFromGaussian(y[0],y[1]);
    
    /* Combine the underlying GRFs to derive Facies*/
    
    dbin->setSimvar(LOC_FACIES,iech,isimu,0,ipgs,nbsimu,1,facies);
    
    /* Conditional printout of the facies at data */
    
    if (debug_query("condexp"))
      st_print_gaus2fac("Data",isimu,iech,y,facies);
  }
  
  /* Set the error return code */
  
  error = 0;
  
label_end:
  return(error);
}

/****************************************************************************/
/*!
**  Combine the underlying GRF into a facies value at data points
**
** \return  Error return code
**
** \param[in]  propdef    Props structure
** \param[in]  dbin       Db input structure
** \param[in]  dbout      Db output structure
** \param[in]  rule       Rule structure
** \param[in]  flag_used  1 if the gaussian is used; 0 otherwise
** \param[in]  ipgs       Rank of the PGS
** \param[in]  isimu      Rank of the simulation
** \param[in]  nbsimu     Number of simulations
**
** \remark Attributes LOC_GAUSFAC are mandatory
** \remark Attributes LOC_FACIES are mandatory
**
*****************************************************************************/
GEOSLIB_API int rule_gaus2fac_data(Props  *propdef,
                                   Db     *dbin,
                                   Db     *dbout,
                                   Rule   *rule,
                                   int    *flag_used,
                                   int     ipgs,
                                   int     isimu,
                                   int     nbsimu)
{
  int    iech,error,igrf,icase;
  double y[2],facies,t1min,t1max,t2min,t2max;

  /* Initializations */

  error = 1;
  check_mandatory_attribute("rule_gaus2fac_data",dbin,LOC_GAUSFAC);

  /* Processing the translation */

  for (iech=0; iech<get_NECH(dbin); iech++)
  {
    if (! dbin->isActive(iech)) continue;
    
    /* Initializations */
    
    facies = TEST;
    for (igrf=0; igrf<2; igrf++) y[igrf] = TEST;
    if (rule_thresh_define(propdef,dbin,rule,ITEST,
                           iech,isimu,nbsimu,1,
                           &t1min,&t1max,&t2min,&t2max)) goto label_end;
    
    for (igrf=0; igrf<2; igrf++)
    {
      icase = get_rank_from_propdef(propdef,ipgs,igrf);
      y[igrf] = (flag_used[igrf]) ?
        dbin->getSimvar(LOC_GAUSFAC,iech,isimu,0,icase,nbsimu,1) : 0.;
    }
    facies = rule->getFaciesFromGaussian(y[0],y[1]);
    
    /* Combine the underlying GRFs to derive Facies */
    
    dbin->setSimvar(LOC_FACIES,iech,isimu,0,ipgs,nbsimu,1,facies);
    
    /* Conditional printout of the facies at data */
    
    if (debug_query("condexp"))
      st_print_gaus2fac("Data",isimu,iech,y,facies);
  }
  
  /* Set the error return code */
  
  error = 0;
  
label_end:
  return(error);
}

/****************************************************************************/
/*!
**  Combine the underlying GRF into a facies value (Shadow case)
**
** \return  Error return code
**
** \param[in]  propdef    Props structure
** \param[in]  dbout      Db output structure
** \param[in]  rule       Rule structure
** \param[in]  flag_used  1 if the gaussian is used; 0 otherwise
** \param[in]  ipgs       Rank of the PGS
** \param[in]  isimu      Rank of the simulation
** \param[in]  nbsimu     Number of simulations
**
** \remark Attributes LOC_FACIES and LOC_SIMU are mandatory
**
*****************************************************************************/
GEOSLIB_API int rule_gaus2fac_result_shadow(Props  *propdef,
                                            Db     *dbout,
                                            Rule   *rule,
                                            int    *flag_used,
                                            int     ipgs,
                                            int     isimu,
                                            int     nbsimu)
{
  int ndim,iech,jech,error,idim,nstep,istep,flag,flag_shadow,igrf,icase;
  double *del,y[2],facies,dinc,dy,ys,yc_dsup,yc_down;
  double  t1min,t1max,t2min,t2max,s1min,s1max,s2min,s2max,sh_dsup,sh_down,seuil;

  /* Initializations */

  check_mandatory_attribute("rule_gaus2fac_result_shadow",dbout,LOC_FACIES);
  check_mandatory_attribute("rule_gaus2fac_result_shadow",dbout,LOC_SIMU);
  error = 1;
  del   = (double *) NULL;
  dy    = 0.;
  nstep = 0;
  ndim  = dbout->getNDim();
  icase = get_rank_from_propdef(propdef,ipgs,0);
  VectorDouble xyz(ndim);
  VectorInt ind1(ndim);
  VectorInt ind2(ndim);

  /* Initializations */

  del = db_vector_alloc(dbout);
  if (del == (double *) NULL) goto label_end;
  dinc  = rule->getIncr();
  nstep = (int) floor(rule->getDMax() / dinc);
  dy    = dinc * rule->getTgte();
  for (idim=0; idim<ndim; idim++) del[idim] = dinc * rule->getShift(idim);

  /* Processing the translation */

  for (iech=0; iech<get_NECH(dbout); iech++)
  {
    if (! dbout->isActive(iech)) continue;
    
    /* Initializations */
    
    facies = TEST;
    for (igrf=0; igrf<2; igrf++) y[igrf] = TEST;
    
    y[0] = dbout->getSimvar(LOC_SIMU,iech,isimu,0,icase,nbsimu,1);
    if (FFFF(y[0])) break;
    if (rule_thresh_define_shadow(propdef,dbout,rule,SHADOW_WATER,
                                  iech,isimu,nbsimu,1,
                                  &t1min,&t1max,&t2min,&t2max,
                                  &yc_dsup,&yc_down)) goto label_end;
    db_index_sample_to_grid(dbout,iech,ind2.data());
    grid_to_point(dbout,ind2.data(),NULL,xyz.data());
    
    if (y[0] >= t1max)
      facies = SHADOW_ISLAND;
    else
    {
      flag_shadow = 0;
      db_index_sample_to_grid(dbout,iech,ind2.data());
      grid_to_point(dbout,ind2.data(),NULL,xyz.data());
      for (istep=1; istep<=nstep && flag_shadow==0; istep++)
      {
        for (idim=0; idim<ndim; idim++) xyz[idim] -= del[idim];
        flag = point_to_grid(dbout,xyz.data(),0,ind2.data());
        if (flag > 0) break;
        if (flag < 0) continue;
        ys = rule->st_grid_eval(dbout,isimu,icase,nbsimu,xyz);
        if (FFFF(ys)) continue;
        jech  = db_index_grid_to_sample(dbout,ind2.data());
        if (rule_thresh_define_shadow(propdef,dbout,rule,SHADOW_WATER,
                                      jech,isimu,nbsimu,1,
                                      &s1min,&s1max,&s2min,&s2max,
                                      &sh_dsup,&sh_down)) return(1);
        if (ys < s1max) continue;  /* Upstream point not in island */
        seuil = t1max - yc_down + dy * istep;
        flag_shadow = (MIN(ys,s1max + sh_dsup) > seuil);
      }
      facies = (flag_shadow) ? SHADOW_SHADOW : SHADOW_WATER;
    }
    
    /* Combine the underlying GRFs to derive Facies */
    
    dbout->setSimvar(LOC_FACIES,iech,isimu,0,ipgs,nbsimu,1,facies);
    
    /* Optional printout */
    
    if (debug_query("condexp")) 
      st_print_gaus2fac("Results",isimu,iech,y,facies);
  }
  
  /* Set the error return code */
  
  error = 0;

label_end:
  del  = db_vector_free(del);
  return(error);
}

/****************************************************************************/
/*!
**  Combine the underlying GRF into a facies value
**
** \return  Error return code
**
** \param[in]  propdef    Props structure
** \param[in]  dbout      Db output structure
** \param[in]  rule       Rule structure
** \param[in]  flag_used  1 if the gaussian is used; 0 otherwise
** \param[in]  ipgs       Rank of the PGS
** \param[in]  isimu      Rank of the simulation
** \param[in]  nbsimu     Number of simulations
**
** \remark Attributes LOC_FACIES and LOC_SIMU are mandatory
**
*****************************************************************************/
GEOSLIB_API int rule_gaus2fac_result(Props  *propdef,
                                     Db     *dbout,
                                     Rule   *rule,
                                     int    *flag_used,
                                     int     ipgs,
                                     int     isimu,
                                     int     nbsimu)
{
  int    ndim,iech,jech,error,idim,igrf,icase;
  double t1min,t1max,t2min,t2max,facies,y[2];

  /* Initializations */

  check_mandatory_attribute("rule_gaus2fac_result",dbout,LOC_FACIES);
  check_mandatory_attribute("rule_gaus2fac_result",dbout,LOC_SIMU);
  error  = 1;
  ndim   = dbout->getNDim();
  VectorDouble xyz(ndim);
  VectorInt ind1(ndim);
  VectorInt ind2(ndim);

  /* Processing the translation */

  for (iech=0; iech<get_NECH(dbout); iech++)
  {
    if (! dbout->isActive(iech)) continue;
    
    /* Initializations */
    
    facies = TEST;
    for (igrf=0; igrf<2; igrf++) y[igrf] = TEST;
    
    switch (rule->getModeRule())
    {
      case RULE_STD:
        if (rule_thresh_define(propdef,dbout,rule,ITEST,
                               iech,isimu,nbsimu,1,
                               &t1min,&t1max,&t2min,&t2max)) goto label_end;
        for (igrf=0; igrf<2; igrf++)
        {
          icase = get_rank_from_propdef(propdef,ipgs,igrf);
          y[igrf] = (flag_used[igrf]) ?
            dbout->getSimvar(LOC_SIMU,iech,isimu,0,icase,nbsimu,1) : 0.;
        }
        facies = rule->getFaciesFromGaussian(y[0],y[1]);
        break;
      
      case RULE_SHIFT:
        icase = get_rank_from_propdef(propdef,ipgs,0);
        y[0] = dbout->getSimvar(LOC_SIMU,iech,isimu,0,icase,nbsimu,1);
        if (FFFF(y[0])) break;
      
        if (rule_thresh_define(propdef,dbout,rule,ITEST,
                               iech,isimu,nbsimu,1,
                               &t1min,&t1max,&t2min,&t2max)) goto label_end;
        db_index_sample_to_grid(dbout,iech,ind2.data());
        for (idim=0; idim<ndim; idim++) ind2[idim] -= ind1[idim];
        jech = db_index_grid_to_sample(dbout,ind2.data());
        if (jech >= 0)
          y[1] = dbout->getSimvar(LOC_SIMU,jech,isimu,0,icase,nbsimu,1);
        else
          y[1] = TEST;
        facies = rule->getFaciesFromGaussian(y[0],y[1]);
        break;
      
      default:
        messageAbort("Other type of rule mode is forbidden");
        break;
    }
    
    /* Combine the underlying GRFs to derive Facies */
    
    dbout->setSimvar(LOC_FACIES,iech,isimu,0,ipgs,nbsimu,1,facies);
    
    /* Optional printout */
    
    if (debug_query("condexp")) 
      st_print_gaus2fac("Results",isimu,iech,y,facies);
  }
  
  /* Set the error return code */
  
  error = 0;

label_end:
  return(error);
}

/****************************************************************************/
/*!
**  Locate the current proportions
**
** \param[in]  propdef    Props structure
** \param[in]  ifac_ref   Conditional (first variable) facies
**                        (Only used for PROCESS_CONDITIONAL)
**
****************************************************************************/
static int st_proportion_locate(Props  *propdef,
                                int     ifac_ref)
{
  int ifac;
  
  switch (propdef->mode)
  {
    case PROCESS_COPY:
    case PROCESS_MARGINAL:
      for (ifac=0; ifac<propdef->nfaccur; ifac++)
        propdef->proploc[ifac] = propdef->propwrk[ifac];
      break;
      
    case PROCESS_CONDITIONAL:
      for (ifac=0; ifac<propdef->nfaccur; ifac++)
        propdef->proploc[ifac] = PROPWRK(ifac_ref-1,ifac);
      break;
  }
  return(0);
}

/****************************************************************************/
/*!
**  Transform the proportions (from CST to WRK)
**
** \return  -1 if the proportion is not defined; 0 otherwise
**
** \param[in]  propdef   Props structure
**
****************************************************************************/
static int st_proportion_transform(Props  *propdef)

{
  double total,pp;
  int    ifac1,ifac2;

  /* Dispatch */
  
  switch (propdef->mode)
  {
    case PROCESS_COPY:
      for (ifac1=0; ifac1<propdef->nfac[0]; ifac1++)
      {
        pp = PROPFIX(ifac1,0);
        if (FFFF(pp)) return(-1);
        PROPWRK(ifac1,0) = pp;
      }
      break;
      
    case PROCESS_MARGINAL:
      for (ifac1=0; ifac1<propdef->nfac[0]; ifac1++)
      {
        PROPWRK(ifac1,0) = 0.;
        for (ifac2=0; ifac2<propdef->nfac[1]; ifac2++)
        {
          pp = PROPFIX(ifac1,ifac2);
          if (FFFF(pp)) return(-1);
          PROPWRK(ifac1,0) += pp;
        }
      }
      break;
      
    case PROCESS_CONDITIONAL:
      for (ifac1=0; ifac1<propdef->nfac[0]; ifac1++)
        for (ifac2=0; ifac2<propdef->nfac[1]; ifac2++)
        {
          pp = PROPFIX(ifac1,ifac2);
          if (FFFF(pp)) return(-1);
          PROPWRK(ifac1,ifac2) = pp;
        }
      
      for (ifac1=0; ifac1<propdef->nfac[0]; ifac1++)
      {
        total = 0.;
        for (ifac2=0; ifac2<propdef->nfac[1]; ifac2++)
          total += PROPWRK(ifac1,ifac2);
        
        for (ifac2=0; ifac2<propdef->nfac[1]; ifac2++)
          PROPWRK(ifac1,ifac2) = (total <= 0.) ? 
            1. / propdef->nfac[1] : PROPWRK(ifac1,ifac2) / total;
      }
      break;
      
    default:
      messageAbort("This should never happen in st_proportion_transform");
      break;
  }
  return(0);
}

/****************************************************************************/
/*!
**  Set the method to compute Proportions
**
** \param[in]  propdef  Props structure
** \param[in]  mode     Type of operation 
** \li                  PROCESS_COPY
** \li                  PROCESS_MARGINAL
** \li                  PROCESS_CONDITIONAL
**
****************************************************************************/
GEOSLIB_API void proportion_rule_process(Props *propdef, int mode)
{
  /* Assignments */

  propdef->mode = mode;

  /* Assign the current value for the number of facies */

  if (mode == PROCESS_COPY ||
      mode == PROCESS_MARGINAL)    propdef->nfaccur = propdef->nfac[0];
  if (mode == PROCESS_CONDITIONAL) propdef->nfaccur = propdef->nfac[1];

  /* In the stationary case, transform the proportions (from CST to WRK) */

  if (propdef->case_stat) st_proportion_transform(propdef);
  
  return;
}

/****************************************************************************/
/*!
**  Print the (non-stationary) proportions
**
** \param[in]  propdef   Props structure
**
*****************************************************************************/
GEOSLIB_API void proportion_print(Props  *propdef)

{
  if (propdef == (Props *) NULL) return;
  mestitle(0,"Proportions");

  print_matrix("Initial :",0,1,propdef->nfac[1],propdef->nfac[0],NULL,
               propdef->propfix.data());
  
  print_matrix("Working :",0,1,propdef->nfac[1],propdef->nfac[0],NULL,
               propdef->propwrk.data());

  print_matrix("Current :",0,1,propdef->nfaccur,1,NULL,
               propdef->proploc.data());
}

/****************************************************************************/
/*!
**  Check if the proportion has changed since the previous usage
**  and store the current proportions for future comparison
**
** \return  1 if the proportions are unchanged; 0 otherwise
**
** \param[in]  propdef Props structure
**
****************************************************************************/
static int st_proportion_changed(Props *propdef)

{
  int ifac,modify;
  
  /* Compare with the memory proportion array */

  modify = 0;
  for (ifac=0; ifac<propdef->nfaccur; ifac++)
    if (propdef->proploc[ifac] != propdef->propmem[ifac]) modify = 1;
  if (! modify) return(1);

  /* Print the proportions (optional) */

  if (debug_query("props")) proportion_print(propdef);

  for (ifac=0; ifac<propdef->nfaccur; ifac++)
    propdef->propmem[ifac] = propdef->proploc[ifac];

  return(0);
}

/****************************************************************************/
/*!
**  Set the (non-stationary) proportions
**
** \return  Error return code
** \return  - the target point does not lie within the proportion grid
** \return  - in conditional processing, the reference facies does not exist
**
** \param[in]  propdef    Props structure
** \param[in]  db         Db input structure
** \param[in]  iech       Rank of the data in the input Db
** \param[in]  isimu      Rank of the simulation (PROCESS_CONDITIONAL)
** \param[in]  nbsimu     Number of simulations
**
** \param[out] jech       Rank of the auxiliary data in the input Db
**
** \remark  At the end of this function, the local proportions are stored
** \remark  in the array proploc of the structure Props
** \remark  The argument 'isimu' is only used for
** \remark            propdef->mode == PROCESS_CONDITIONAL (simbipgs)
**
*****************************************************************************/
static int st_proportion_define(Props  *propdef,
                                Db     *db,
                                int     iech,
                                int     isimu,
                                int     nbsimu,
                                int    *jech)
{
  int ifac,ifac_ref;

  /* Non-stationary case : Load the proportions in propcst */

  (*jech) = 0;
  if (! propdef->case_stat)
  {
    if (propdef->case_prop_interp)
    {
    
      /* Case where the proportions must be interpolated */

      (*jech) = index_point_to_grid(db,iech,1,propdef->dbprop,
                                    propdef->coor.data());
      if ((*jech) < 0)
      {
        messerr("At the data #%d, the proportion matrix is undefined",iech+1);
        return(1);
      }
      
      /* Load the proportions (into CST) */
      
      for (ifac=0; ifac<propdef->nfacprod; ifac++)
        propdef->propfix[ifac] = propdef->dbprop->getProportion(*jech,ifac);
    }
    else
    {

      /* The proportions are already available from the dbin */

      for (ifac=0; ifac<propdef->nfacprod; ifac++)
        propdef->propfix[ifac] = db->getProportion(iech,ifac);
    }
      
    /* Transform proportions (from CST to WRK) */
    
    st_proportion_transform(propdef);
  }

  /* Locate the current proportions (from WRK to LOC) */

  ifac_ref = -1;
  if (propdef->mode == PROCESS_CONDITIONAL)
  {
    ifac_ref = (int) db->getSimvar(LOC_FACIES,iech,isimu,0,0,nbsimu,1);
    if (ifac_ref < 1 || ifac_ref > propdef->nfac[0]) return(1);
  }
  st_proportion_locate(propdef,ifac_ref);

  return(0);
}

/****************************************************************************/
/*!
**  Set the (non-stationary) proportions and define thresholds (for shadow only)
**
** \return  Error return code
**
** \param[in]  propdef    Props structure
** \param[in]  db         Db input structure
** \param[in]  rule       Rule structure
** \param[in]  facies     Facies of interest (or GV_ITEST)
** \param[in]  iech       Rank of the data in the input Db
** \param[in]  isimu      Rank of the simulation (PROCESS_CONDITIONAL)
** \param[in]  nbsimu     Number of simulations (PROCESS_CONDITIONAL)
** \param[in]  flag_check 1 if the consistency check with the actual
**                        proportion of the current facies must be done
**
** \param[out] t1min      Minimum threshold for Y1
** \param[out] t1max      Maximum threshold for Y1
** \param[out] t2min      Minimum threshold for Y2
** \param[out] t2max      Maximum threshold for Y2
** \param[out] sh_dsup    Local or global upwards shift (shadow)
** \param[out] sh_down    Local or global downwards shift (shadow)
**
*****************************************************************************/
GEOSLIB_API int rule_thresh_define_shadow(Props  *propdef,
                                          Db     *db,
                                          Rule   *rule,
                                          int     facies,
                                          int     iech,
                                          int     isimu,
                                          int     nbsimu,
                                          int     flag_check,
                                          double *t1min,
                                          double *t1max,
                                          double *t2min,
                                          double *t2max,
                                          double *sh_dsup,
                                          double *sh_down)
{
  int    unmodify,facloc,jech;

  /* Set the debugging information */

  debug_index(iech+1);

  /* Processing an "unknown" facies */

  if (! IFFFF(facies) && (facies < 1 || facies > propdef->nfaccur))
  {
    *t1min = *t2min = get_rule_extreme(-1);
    *t1max = *t2max = get_rule_extreme(+1);
    return(0);
  }

  /* Define the proportions */

  if (st_proportion_define(propdef,db,iech,isimu,nbsimu,&jech)) 
  {
    *t1min = *t2min = get_rule_extreme(-1);
    *t1max = *t2max = get_rule_extreme(+1);
    return(0);
  }

  /* Check if the proportions have been changed */

  unmodify = st_proportion_changed(propdef);

  /* In case of Shadow, return the upwards and downwards values */

  *sh_dsup = (propdef->case_stat) ? rule->getShDsup() : propdef->proploc[1];
  *sh_down = (propdef->case_stat) ? rule->getShDown() : propdef->proploc[2];

  /* In the special cases, only the first proportion is significant */

  propdef->proploc[1] = (1 - propdef->proploc[0]) / 2;
  propdef->proploc[2] = (1 - propdef->proploc[0]) / 2;

  /* Set the proportions and translate proportions into thresholds */

  if (! unmodify)
  {
    if (rule->setProportions(propdef->proploc)) return(1);
  }

  /* Convert the proportions into thresholds */

  facloc = (IFFFF(facies)) ? 1 : facies;
  VectorDouble bounds = rule->getThresh(facloc);
  *t1min = bounds[0];
  *t1max = bounds[1];
  *t2min = bounds[2];
  *t2max = bounds[3];

  return(0);
}

/****************************************************************************/
/*!
**  Set the (non-stationary) proportions and define thresholds
**
** \return  Error return code
**
** \param[in]  propdef    Props structure
** \param[in]  db         Db input structure
** \param[in]  rule       Rule structure
** \param[in]  facies     Facies of interest (or ITEST)
** \param[in]  iech       Rank of the data in the input Db
** \param[in]  isimu      Rank of the simulation (PROCESS_CONDITIONAL)
** \param[in]  nbsimu     Number of simulations
** \param[in]  flag_check 1 if the consistency check with the actual
**                        proportion of the current facies must be done
**
** \param[out] t1min      Minimum threshold for Y1
** \param[out] t1max      Maximum threshold for Y1
** \param[out] t2min      Minimum threshold for Y2
** \param[out] t2max      Maximum threshold for Y2
**
*****************************************************************************/
GEOSLIB_API int rule_thresh_define(Props  *propdef,
                                   Db     *db,
                                   Rule   *rule,
                                   int     facies,
                                   int     iech,
                                   int     isimu,
                                   int     nbsimu,
                                   int     flag_check,
                                   double *t1min,
                                   double *t1max,
                                   double *t2min,
                                   double *t2max)
{
  int    unmodify,facloc,jech;

  /* Set the debugging information */

  debug_index(iech+1);

  /* Processing an "unknown" facies */

  if (! IFFFF(facies) && (facies < 1 || facies > propdef->nfaccur))
  {
    *t1min = *t2min = get_rule_extreme(-1);
    *t1max = *t2max = get_rule_extreme(+1);
    return(0);
  }

  /* Define the proportions */

  if (st_proportion_define(propdef,db,iech,isimu,nbsimu,&jech)) 
  {
    *t1min = *t2min = get_rule_extreme(-1);
    *t1max = *t2max = get_rule_extreme(+1);
    return(0);
  }

  /* Check if the proportions have been changed */

  unmodify = st_proportion_changed(propdef);

  /* Check that the facies is compatible with the proportions */

  if (flag_check && ! IFFFF(facies) && rule->getModeRule() == RULE_STD)
  {
    if (propdef->proploc[facies-1] <= 0.) 
    {
      messerr("The presence of facies (%d) at sample (%d) is not consistent with the zero proportion",
              facies,iech+1);
      if (! propdef->case_stat)
        messerr("Check the proportions in the cell (%d) of the Proportion Db",
                jech+1);
      return(1);
    }
  }

  /* Set the proportions and translate proportions into thresholds */

  if (! unmodify)
  {
    if (rule->setProportions(propdef->proploc)) return(1);

    /* In the case of SHIFT, update the thresholds */

    if (rule->getModeRule() == RULE_SHIFT && 0)
      rule->updateShift();
  }

  /* Convert the proportions into thresholds */

  facloc = (IFFFF(facies)) ? 1 : facies;
  VectorDouble bounds = rule->getThresh(facloc);
  *t1min = bounds[0];
  *t1max = bounds[1];
  *t2min = bounds[2];
  *t2max = bounds[3];

  return(0);
}

/****************************************************************************/
/*!
**  Apply the Rule transformation to the GRFs of a Db
**  (Shadow case)
**
** \return  Error return code
**
** \param[in]  db        Output Db structure
** \param[in]  dbprop    Db structure used for proportions (non-stationary case)
** \param[in]  rule      Lithotype Rule definition
** \param[in]  model     First Model structure (only for SHIFT)
** \param[in]  props     Array of proportions for the facies
** \param[in]  flag_stat 1 for stationary; 0 otherwise
** \param[in]  nfacies   Number of facies
**
** \remark The input variable must be locatorized as Z or LOC_SIMU
** \remark It will be changed in this function to locator LOC_SIMU
**
*****************************************************************************/
GEOSLIB_API int db_rule_shadow(Db     *db,
                               Db     *dbprop,
                               Rule   *rule,
                               Model  *model,
                               const VectorDouble& props,
                               int     flag_stat,
                               int     nfacies)
{
  int    iptr,error,flag_used[2],nbsimu,igrf,ngrf;
  Props *propdef;

  /* Initializations */

  error   = 1;
  nbsimu  = 1;
  iptr    = -1;
  propdef = (Props *) NULL;

  /* Preliminary checks */

  ngrf = rule->getGRFNumber();
  for (igrf=0; igrf<2; igrf++) 
    flag_used[igrf] = rule->isYUsed(igrf);

  propdef = proportion_manage(1,1,flag_stat,ngrf,0,nfacies,0,
                              db,dbprop,props,propdef);
  if (propdef == (Props *) NULL) goto label_end;

  /* General setting for lithotype */

  rule->particularities_shadow(db, dbprop, model,1,flag_stat);
  proportion_rule_process(propdef,PROCESS_COPY);

  /**********************/
  /* Add the attributes */
  /**********************/

  /* Storage of the simulations in the output file */
  iptr = db->addFields(nbsimu,0.);
  if (iptr < 0) goto label_end;
  db->setLocatorsByAttribute(nbsimu,iptr,LOC_FACIES);

  /* Identify the Non conditional simulations at target points */
  for (igrf=0; igrf<2; igrf++)
  {
    if (! flag_used[igrf]) continue;
    iptr = db_attribute_identify(db,LOC_SIMU,igrf);
    if (iptr < 0)
    {
      iptr = db_attribute_identify(db,LOC_Z,igrf);
      if (iptr < 0)
      {
        messerr("The variable containing the simulation of the GRF %d is missing in the Db",igrf+1);
        goto label_end;
      }
      db->setLocatorByAttribute(iptr,LOC_SIMU,igrf+1);
    }
  }

  /* Combine the conditional simulation for each GRF */

  for (int isimu=0; isimu<nbsimu; isimu++)
    if (rule_gaus2fac_result_shadow(propdef,db,rule,flag_used,
                                    0,isimu,nbsimu)) goto label_end;

  /* Set the error return flag */

  error = 0;

label_end:
  propdef = proportion_manage(-1,1,flag_stat,ngrf,0,nfacies,0,
                              db,dbprop,props,propdef);
  return(error);
}

/****************************************************************************/
/*!
**  Apply the Rule transformation to convert a set of Gaussian vectors
**  into the corresponding Facies in a Db
**
** \return  Error return code
**
** \param[in]  db        Output Db structure
** \param[in]  dbprop    Db structure used for proportions (non-stationary case)
** \param[in]  rule      Lithotype Rule definition
** \param[in]  props     Array of proportions for the facies
** \param[in]  flag_stat 1 for stationary; 0 otherwise
** \param[in]  model     First Model structure (only for SHIFT)
** \param[in]  namconv   Naming convention
**
** \remark The input variable must be locatorized as Z or LOC_SIMU
**
*****************************************************************************/
GEOSLIB_API int db_rule(Db     *db,
                        Rule   *rule,
                        const   VectorDouble& props,
                        Db     *dbprop,
                        int     flag_stat,
                        Model  *model,
                        NamingConvention namconv)
{
  int error = 1;
  int iptr    = -1;
  Props* propdef = (Props *) NULL;
  int ngrf    = rule->getGRFNumber();
  VectorInt flagUsed = rule->whichGRFUsed();
  int nfacies = rule->getFaciesNumber();
  bool flagReturn = false;

  /* Preliminary checks */

  if (db->getLocatorNumber(LOC_SIMU) != ngrf &&
      db->getLocatorNumber(LOC_Z) != ngrf)
  {
    messerr("The Rule specifies the use of %d underlying GRF(s)",ngrf);
    messerr("The input 'db' should have one variable per GRF with locator 'SIMU' or 'Z'");
    goto label_end;
  }

  propdef = proportion_manage(1,1,flag_stat,ngrf,0,nfacies,0,
                              db,dbprop,props,propdef);
  if (propdef == (Props *) NULL) goto label_end;
  if (rule->particularities(db,dbprop,model,1,flag_stat)) goto label_end;
  proportion_rule_process(propdef,PROCESS_COPY);

  /**********************/
  /* Add the attributes */
  /**********************/

  /* Storage of the simulations in the output file */
  iptr = db->addFields(1,0.,"Facies",LOC_FACIES);
  if (iptr < 0) goto label_end;

  /* Identify the Non conditional simulations at target points */

  for (int igrf=0; igrf<2; igrf++)
  {
    if (! flagUsed[igrf]) continue;
    if (db->getLocatorNumber(LOC_SIMU) == ngrf)
    {
      iptr = db_attribute_identify(db,LOC_SIMU,igrf);
      flagReturn = false;
    }
    else
    {
      iptr = db_attribute_identify(db,LOC_Z,igrf);
      flagReturn = true;
    }
    db->setLocatorByAttribute(iptr,LOC_SIMU,igrf+1);
  }

  /* Translate Gaussian into Facies */

  if (rule_gaus2fac_result(propdef,db,rule,flagUsed.data(),0,0,1)) goto label_end;

  // Returning to the initial locators (if the initial variable
  // had a LOC_Z locator which has been temporarily modified into LOC_SIMU)

  if (flagReturn)
  {
    for (int igrf=0; igrf<2; igrf++)
    {
      if (! flagUsed[igrf]) continue;
      iptr = db_attribute_identify(db,LOC_SIMU,igrf);
      db->setLocatorByAttribute(iptr,LOC_SIMU,igrf+1);
    }
  }

  // Naming convention

  namconv.setNamesAndLocators(nullptr, VectorInt(), db, iptr);
  error = 0;

label_end:
  propdef = proportion_manage(-1,1,flag_stat,ngrf,0,nfacies,0,
                              db,dbprop,props,propdef);
  return(error);
}

/****************************************************************************/
/*!
**  Apply the Rule transformation to derive the bounds variables for a Db
**  (Shadow case)
**
** \return  Error return code
**
** \param[in]  db        Db structure
** \param[in]  dbprop    Db structure used for proportions (non-stationary case)
** \param[in]  rule      Lithotype Rule definition
** \param[in]  model     First Model structure (only for SHIFT)
** \param[in]  props     Array of proportions for the facies
** \param[in]  flag_stat 1 for stationary; 0 otherwise
** \param[in]  nfacies   Number of facies
**
*****************************************************************************/
GEOSLIB_API int db_bounds_shadow(Db     *db,
                                 Db     *dbprop,
                                 Rule   *rule,
                                 Model  *model,
                                 const   VectorDouble& props,
                                 int     flag_stat,
                                 int     nfacies)
{
  int     flag_used[2],ngrf,error,iptr,igrf;
  double *coor;
  Props  *propdef;

  /* Initializations */

  error   = 1;
  ngrf    = 0;
  coor    = (double *) NULL;
  propdef = (Props  *)NULL;

  /**********************/
  /* Preliminary checks */
  /**********************/

  /* Input Db */

  if (db == (Db *) NULL)
  {
    messerr("The Db is not defined");
    goto label_end;
  }
  if (! db->isVariableNumberComparedTo(1)) goto label_end;

  /* Rule */

  if (rule == (Rule *) NULL)
  {
    messerr("The Rule is not defined");
    goto label_end;
  }
  ngrf = rule->getGRFNumber();
  for (igrf=0; igrf<2; igrf++) 
    flag_used[igrf] = rule->isYUsed(igrf);

  /*******************/
  /* Core allocation */
  /*******************/

  coor = db_sample_alloc(db,LOC_X);
  if (coor == (double *) NULL) goto label_end;

  propdef = proportion_manage(1,1,flag_stat,ngrf,0,nfacies,0,
                              db,dbprop,props,propdef);
  if (propdef == (Props *) NULL) goto label_end;

  /* General setting for lithotype */

  rule->particularities_shadow(db, dbprop, model,1,flag_stat);
  proportion_rule_process(propdef,PROCESS_COPY);

  /**********************/
  /* Add the attributes */
  /**********************/

  /* Lower bound at input data points */
  if (db_locator_attribute_add(db,LOC_L,ngrf,0,0.,&iptr))
    goto label_end;
  
  /* Upper bound at input data points */
  if (db_locator_attribute_add(db,LOC_U,ngrf,0,0.,&iptr))
    goto label_end;

  /* Calculate the thresholds and store them in the Db file */

  for (igrf=0; igrf<ngrf; igrf++)
  {
    if (! flag_used[igrf]) continue;
    if (rule_evaluate_bounds_shadow(propdef,db,db,rule,0,igrf,0,0,0.)) 
      goto label_end;
  }

  /* Set the error return flag */

  error = 0;

label_end:
  propdef = proportion_manage(-1,1,flag_stat,ngrf,0,nfacies,0,
                              db,dbprop,props,propdef);
  coor = db_sample_free(coor);
  return(error);
}

/****************************************************************************/
/*!
**  Apply the Rule transformation to derive the bounds variables for a Db
**
** \return  Error return code
**
** \param[in]  db        Db structure
** \param[in]  rule      Lithotype Rule definition
** \param[in]  props     Array of proportions for the facies
** \param[in]  dbprop    Db structure used for proportions (non-stationary case)
** \param[in]  flag_stat 1 for stationary; 0 otherwise
** \param[in]  model     First Model structure (only for SHIFT)
** \param[in]  namconv   Naming convention
**
*****************************************************************************/
GEOSLIB_API int db_bounds(Db     *db,
                          Rule   *rule,
                          const   VectorDouble& props,
                          Db     *dbprop,
                          int     flag_stat,
                          Model  *model,
                          NamingConvention namconv)
{
  int iptrl, iptru, nvar;
  VectorInt flagUsed;
  int nfacies = 0;
  int ngrf = 0;
  Props* propdef = nullptr;
  int error = 1;

  /* Input Db */

  if (db == (Db *) NULL)
  {
    messerr("The Db is not defined");
    goto label_end;
  }
  nvar = db->getVariableNumber();
  if (! db->isVariableNumberComparedTo(1)) goto label_end;

  /* Rule */

  if (rule == (Rule *) NULL)
  {
    messerr("The Rule is not defined");
    goto label_end;
  }
  ngrf = rule->getGRFNumber();
  flagUsed = rule->whichGRFUsed();
  nfacies = rule->getFaciesNumber();

  /* Model */

  if (rule->getModeRule() == RULE_SHIFT)
  {
    if (model == (Model *) NULL)
    {
      messerr("No Model is provided");
      goto label_end;
    }
    if (model->getVariableNumber() != nvar)
    {
      messerr("The number of variables in the Model (%d) does not match",
              model->getVariableNumber());
      messerr(" the number of variables in the Db (%d)",nvar);
      goto label_end;
    }
  }

  /*******************/
  /* Core allocation */
  /*******************/

  propdef = proportion_manage(1, 1, flag_stat, ngrf, 0, nfacies, 0, db, dbprop,
                              props, propdef);
  if (propdef == (Props *) NULL) goto label_end;

  /* General setting for lithotype */

  if (rule->particularities(db, dbprop, model, 1, flag_stat)) goto label_end;
  proportion_rule_process(propdef, PROCESS_COPY);

  /**********************/
  /* Add the attributes */
  /**********************/

  /* Lower bound at input data points */
  if (db_locator_attribute_add(db,LOC_L,ngrf,0,0.,&iptrl))
    goto label_end;
  
  /* Upper bound at input data points */
  if (db_locator_attribute_add(db,LOC_U,ngrf,0,0.,&iptru))
    goto label_end;

  /* Calculate the thresholds and store them in the Db file */

  for (int igrf=0; igrf<ngrf; igrf++)
  {
    if (! flagUsed[igrf]) continue;
    if (rule_evaluate_bounds(propdef,db,db,rule,0,igrf,0,0)) goto label_end;
  }

  // Naming convention

  namconv.setLocatorOutType(LOC_L);
  namconv.setNamesAndLocators(nullptr, VectorInt(), db, iptrl, "Lower", ngrf);
  namconv.setLocatorOutType(LOC_U);
  namconv.setNamesAndLocators(nullptr, VectorInt(), db, iptru, "Upper", ngrf);
  error = 0;

 label_end:
  propdef = proportion_manage(-1,1,flag_stat,ngrf,0,nfacies,0,
                              db,dbprop,props,propdef);
  return(error);
}

/****************************************************************************/
/*!
**  Allocate or deallocate a proportion array
**
** \return  Pointer on the returned Props structure
**
** \param[in]  mode        1 for allocation; -1 for deallocation
** \param[in]  flag_facies 1 if Gibbs is used for facies
** \param[in]  flag_stat   1 if the proportions are stationary
** \param[in]  ngrf1       Number of GRFs for the first PGS
** \param[in]  ngrf2       Number of GRFs for the second PGS
** \param[in]  nfac1       Number of facies for the first PGS
** \param[in]  nfac2       Number of facies for the second PGS
** \param[in]  db          Db structure containing the data
** \param[in]  dbprop      Db structure containing the proportions
**                         (only used in the non-stationary case)
** \param[in]  propcst     Constant set of proportions (used if flag_stat)
** \param[in]  proploc     Props structure (used for mode<0)
**
****************************************************************************/
GEOSLIB_API Props *proportion_manage(int     mode,
                                     int     flag_facies,
                                     int     flag_stat,
                                     int     ngrf1,
                                     int     ngrf2,
                                     int     nfac1,
                                     int     nfac2,
                                     Db     *db,
                                     Db     *dbprop,
                                     const   VectorDouble& propcst,
                                     Props  *proploc)
{
  int ifac,error,nfacmax,nfacprod;
  Db *db_loc;
  Props *propdef;

  /* Initializations */

  error    = 1;
  nfacmax  = MAX(nfac1,nfac2);
  nfacprod = nfac1;
  if (nfac2 > 0) nfacprod *= nfac2;

  /* Dispatch */

  if (mode > 0)
  {
    propdef = new Props;
    propdef->case_facies = flag_facies;
    propdef->case_stat   = flag_stat;
    propdef->case_prop_interp = (dbprop != (Db *) NULL && is_grid(dbprop));
    propdef->ngrf[0]     = ngrf1;
    propdef->ngrf[1]     = ngrf2;
    propdef->nfac[0]     = nfac1;
    propdef->nfac[1]     = nfac2;
    propdef->nfaccur     = nfac1;
    propdef->nfacprod    = nfacprod;
    propdef->mode        = PROCESS_UNDEFINED;
    if (propdef->nfaccur <= 0)
    {
      messerr(" The number of facies may not be zero");
      goto label_end;
    }
    propdef->propfix.resize(nfacprod,0.);
    propdef->propwrk.resize(nfacprod,0.);
    propdef->proploc.resize(nfacprod,0.);
    propdef->propmem.resize(nfacprod,0.);

    if (flag_facies)
    {
      
      // Case of facies: Use of the proportions
      
      if (!flag_stat)
      {
        // Non-stationary case

        db_loc = (propdef->case_prop_interp) ? dbprop : db;
        if (db == nullptr)
        {
          messerr("You have requested Non-stationary proportions");
          messerr("No file is provided containing Proportion variables");
          messerr("Please provide variables with 'proportion' locators");
          messerr("either in the input 'Db' or in 'dbprop'");
          goto label_end;
        }
        if (db_loc->getProportionNumber() != nfacprod)
        {
          messerr(
              "In the non-stationary case, the number of proportion variables (%d)",
              db_loc->getProportionNumber());
          messerr(
              "must be equal to the number of facies (%d) in the Lithotype Rule",
              nfacprod);
          goto label_end;
        }
        propdef->dbprop = db_loc;
        propdef->coor.resize(db_loc->getNDim());
      }
      else
      {
        
        // Stationary case
        
        double pref = 1. / (double) nfacprod;
        for (ifac=0; ifac<nfacprod; ifac++) 
        {
          propdef->propfix[ifac] = (propcst.empty()) ? pref : propcst[ifac];
          propdef->propwrk[ifac] = (propcst.empty()) ? pref : propcst[ifac];
          propdef->proploc[ifac] = (propcst.empty()) ? pref : propcst[ifac];
          propdef->propmem[ifac] = (propcst.empty()) ? pref : propcst[ifac];
        }
      }    
      
      /* Set memory proportion so as to provoke the update at first usage */
      
      for (ifac=0; ifac<nfacmax; ifac++) propdef->propmem[ifac] = -1;
    }
  }
  else
  {
    propdef = proploc;
    if (propdef == (Props *) NULL) return(propdef);

    /* Deallocation */

    propdef->nfaccur  = 0;
    propdef->nfacprod = 0;
    propdef->dbprop   = (Db *) NULL; 

    delete propdef;
    propdef = nullptr;
  }

  /* Set the error return code */

  error = 0;

label_end:
  if (error)
  {
    if (propdef != nullptr) delete propdef;
    propdef = nullptr;
  }
  return(propdef);
}

/****************************************************************************/
/*!
**  Calculate all the thresholds at each sample of a Db
**
** \return  Error return code
**
** \param[in]  db        Db structure
** \param[in]  rule      Lithotype Rule definition
** \param[in]  propcst   Array of proportions for the facies
** \param[in]  dbprop    Db structure used for proportions (non-stationary case)
** \param[in]  flag_stat 1 for stationary; 0 otherwise
** \param[in]  model     First Model structure (only for SHIFT)
** \param[in]  namconv   Naming Convention
**
*****************************************************************************/
GEOSLIB_API int db_threshold(Db *db,
                             Rule *rule,
                             const VectorDouble& propcst,
                             Db *dbprop,
                             int flag_stat,
                             Model *model,
                             NamingConvention namconv)
{
  int    rank, iptr;
  double t1min,t1max,t2min,t2max;

  /* Initializations */

  int error   = 1;
  int ngrf = 0;
  int nfacies = 0;
  Props* propdef = (Props *) NULL;

  /**********************/
  /* Preliminary checks */
  /**********************/

  /* Input Db */

  if (db == (Db *) NULL)
  {
    messerr("The Db is not defined");
    goto label_end;
  }

  /* Rule */

  if (rule == (Rule *) NULL)
  {
    messerr("The Rule is not defined");
    goto label_end;
  }
  if (rule->getModeRule() != RULE_STD)
  {
    messerr("This function is only programmed for standard rule");
    goto label_end;
  }
  ngrf = rule->getGRFNumber();

  /* Model */

  if (rule->getModeRule() == RULE_SHIFT)
  {
    if (model == (Model *) NULL)
    {
      messerr("No Model is provided");
      goto label_end;
    }
  }

  /*******************/
  /* Core allocation */
  /*******************/

  nfacies = rule->getFaciesNumber();
  propdef = proportion_manage(1,1,flag_stat,ngrf,0,nfacies,0,
                              db,dbprop,propcst,propdef);
  if (propdef == (Props *) NULL) goto label_end;

  if (rule->particularities(db,dbprop,model,1,flag_stat)) goto label_end;
  proportion_rule_process(propdef,PROCESS_COPY);

  /**********************/
  /* Add the attributes */
  /**********************/

  iptr = db->addFields(2 * ngrf * nfacies,0.);
  if (iptr < 0) goto label_end;

  /* Calculate the thresholds and store them in the Db file */

  for (int iech=0; iech<get_NECH(db); iech++)
  {
    if (! db->isActive(iech)) continue;
    rank = 0;
    for (int ifac=0; ifac<nfacies; ifac++)
    {
      if (rule_thresh_define(propdef,db,rule,ifac+1,iech,0,0,0,
                             &t1min,&t1max,&t2min,&t2max)) goto label_end;
      db->setArray(iech,iptr+rank,t1min);
      rank++;
      db->setArray(iech,iptr+rank,t1max);
      rank++;
      if (ngrf == 1) continue;
      db->setArray(iech,iptr+rank,t2min);
      rank++;
      db->setArray(iech,iptr+rank,t2max);
      rank++;
    }
  }

  // Naming convention

  rank = 0;
  for (int ifac = 0; ifac < nfacies; ifac++)
  {
    namconv.setNamesAndLocators(db, iptr + rank,
        concatenateStrings("Thresh-F", intToString(ifac + 1), "-Y1-Low"));
    rank++;
    namconv.setNamesAndLocators(db, iptr + rank,
        concatenateStrings("Thresh-F", intToString(ifac + 1), "-Y1-Up"));
    rank++;
    if (ngrf == 1) continue;
    namconv.setNamesAndLocators(db, iptr + rank,
        concatenateStrings("Thresh-F", intToString(ifac + 1), "-Y2-Low"));
    rank++;
    namconv.setNamesAndLocators(db, iptr + rank,
        concatenateStrings("Thresh-F", intToString(ifac + 1), "-Y2-Up"));
    rank++;
  }
  error = 0;

label_end:
  propdef = proportion_manage(-1,1,flag_stat,ngrf,0,nfacies,0,
                              db,dbprop,propcst,propdef);
  return(error);
}

/****************************************************************************/
/*!
**  Combine two basic models into a bivariate model (residuals model)
**
** \return  The newly Model structure
**
** \param[in]  model1      First input Model
** \param[in]  model2      Second input Model
** \param[in]  rule        Rule
**
** \remarks: The drift is not copied into the new model
**
*****************************************************************************/
GEOSLIB_API Model *model_rule_combine(Model *model1,
                                      Model *model2,
                                      Rule  *rule)
{
  Model *new_model;
  int ngrf;
  double rho;

  /* Initializations */

  ngrf = 0;
  new_model = (Model *) NULL;

  /* Preliminary checks */

  if (rule == (Rule *) NULL)
  {
    messerr("This function requires a valid rule.");
    return(new_model);
  }
  if (model1 == (Model *) NULL)
  {
    messerr("This function requires the first model to be defined");
    return(new_model);
  }
  ngrf = rule->getGRFNumber();

  /* Case of a bivariate input model or monogaussian: simply duplicate */

  if (model1->getVariableNumber() == 2 || ngrf == 1)
  {
    new_model = model_duplicate(model1,0.,0);
    return(new_model);
  }

  /* If model2 is not defined, consider model1 */

  if (model2 == (Model *) NULL) 
  {
    if (rule->getModeRule() == RULE_SHIFT)
    { 
      new_model = model_duplicate(model1,0.,0);
      return(new_model);
    } 

    model2 = model1;
  }

  /* Subsequent checks */

  if (model1->getVariableNumber() != 1 || model2->getVariableNumber() != 1)
  {
    messerr("This function can only combine monovariate models");
    return(new_model);
  }
  if (model1->getDimensionNumber() != model2->getDimensionNumber())
  {
    messerr("The two models to be combined must share the space dimension");
    return(new_model);
  }
  if (model1->isFlagLinked() || model2->isFlagLinked())
  {
    messerr("This function cannot combine models with linked drifts");
    return(new_model);
  }

  /* Calculate the correlation coefficient */

  rho = 0.;
  if (rule->getModeRule() == RULE_STD) rho = rule->getRho();

  new_model = model_combine(model1,model2,rho);
  return(new_model);
}
