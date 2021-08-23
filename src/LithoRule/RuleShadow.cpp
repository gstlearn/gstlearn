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
#include "Basic/Utilities.hpp"
#include "Basic/String.hpp"
#include "Basic/AException.hpp"
#include "LithoRule/RuleShadow.hpp"
#include "LithoRule/Rule.hpp"
#include "LithoRule/Node.hpp"
#include "LithoRule/PropDef.hpp"
#include "Model/Model.hpp"
#include "Basic/Law.hpp"
#include "geoslib_f.h"
#include "geoslib_enum.h"

/**
 * Defining the Lithotype Rule in the case of Shadow
 * @param slope   Slope definition
 * @param sh_dsup Maximum threshold
 * @param sh_down Minimum threshold
 * @param shift   Shift orientation
 */
RuleShadow::RuleShadow(double slope,
           double sh_dsup,
           double sh_down,
           const VectorDouble& shift)
    : Rule(),
      _shDsup(sh_dsup),
      _shDown(sh_down),
      _slope(slope),
      _dMax(TEST),
      _tgte(TEST),
      _incr(TEST),
      _shift(shift),
      _xyz(),
      _ind1(),
      _ind2()
{
  setModeRule(RULE_SHADOW);
  VectorString nodnames = {"S", "T", "F", "F", "F"};
  setMainNodeFromNodNames(nodnames);

//  rule->main         = st_node_alloc("S1",THRESH_Y1,0);
//  rule->main->r1     = st_node_alloc("T1",THRESH_Y2,0);
//  rule->main->r1->r1 = st_node_alloc("F3",THRESH_IDLE,SHADOW_SHADOW);
//  rule->main->r1->r2 = st_node_alloc("F2",THRESH_IDLE,SHADOW_WATER);
//  rule->main->r2     = st_node_alloc("F1",THRESH_IDLE,SHADOW_ISLAND);
}

RuleShadow::RuleShadow(const RuleShadow& m)
    : _shDsup(m._shDsup),
      _shDown(m._shDown),
      _slope(m._slope),
      _dMax(m._dMax),
      _tgte(m._tgte),
      _incr(m._incr),
      _shift(m._shift)
{
}

RuleShadow& RuleShadow::operator=(const RuleShadow& m)
{
  if (this != &m)
  {
    _shDsup = m._shDsup;
    _shDown = m._shDown;
    _slope = m._slope;
    _dMax = m._dMax;
    _tgte = m._tgte;
    _incr = m._incr;
    _shift = m._shift;
  }
  return *this;
}

RuleShadow::~RuleShadow()
{
}

/****************************************************************************/
/*!
**  Define the particularities of the PGS model (for Shadow)
**
** \return  Error return code
**
** \param[in]  db              Db structure
** \param[in]  dbprop          Db structure used for proportions
** \param[in]  model           Model structure (only used for shift option)
** \param[in]  flag_grid_check 1 if grid is compulsory; 0 otherwise
**                             (only for SHIFT)
** \param[in]  flag_stat       1 for stationary; 0 otherwise
**
*****************************************************************************/
int RuleShadow::particularities(Db *db,
                                const Db *dbprop,
                                Model *model,
                                int flag_grid_check,
                                int flag_stat)
{
  double sh_dsup_max,sh_down_max;
  int ndim = (model != (Model *) NULL) ? model->getDimensionNumber() : 0;

  double norme = 0.;
  for (int idim=0; idim<ndim; idim++)
    norme += _shift[idim] * _shift[idim];
  norme = sqrt(norme);
  if (norme <= 0) return(1);
  for (int idim=0; idim<ndim; idim++) _shift[idim] /= norme;

  _incr = 1.e30;
  for (int idim=0; idim<ndim; idim++)
    if (_shift[idim] != 0) _incr = MIN(_incr, db->getDX(idim));

  /* Calculate the maximum distance */

  _tgte = tan(ut_deg2rad(_slope));
  _st_shadow_max(dbprop,flag_stat,&sh_dsup_max,&sh_down_max);
  _dMax = (_tgte > 0) ? (sh_dsup_max + sh_down_max) / _tgte : 0.;

  return(0);
}

void RuleShadow::_st_shadow_max(const Db *dbprop,
                                int flag_stat,
                                double *sh_dsup_max,
                                double *sh_down_max)
{
  int iech;
  double val2,val3;

  if (flag_stat || dbprop == (Db *) NULL)
  {
    /* Stationary case */

    *sh_dsup_max = getShDsup();
    *sh_down_max = getShDown();
  }
  else
  {
    *sh_dsup_max = *sh_down_max = 0.;
    for (iech=0; iech<dbprop->getSampleNumber(); iech++)
    {
      val2 = dbprop->getProportion(iech,1);
      if (val2 > (*sh_dsup_max)) (*sh_dsup_max) = val2;
      val3 = dbprop->getProportion(iech,2);
      if (val3 > (*sh_down_max)) (*sh_down_max) = val3;
    }
  }
  return;
}

/****************************************************************************/
/*!
**  Evaluation of the value from a grid by inverse square distance
**  interpolation from the 4 surrounding grid nodes
**
** \return  Interpolated value or FFFF if out of grid
**
** \param[in]  dbgrid Db structure
** \param[in]  isimu  Rank of the simulation
** \param[in]  icase  Rank of the Simulation storage
** \param[in]  nbsimu Number of simulations
**
** \param[out]  xyz0  Working array
**
*****************************************************************************/
double RuleShadow::st_grid_eval(Db *dbgrid,
                                int isimu,
                                int icase,
                                int nbsimu,
                                VectorDouble& xyz0)
{
  double top = 0.;
  double bot = 0.;
  int ndim = dbgrid->getNDim();

  /* First point */
  int iech = db_index_grid_to_sample(dbgrid,_ind2.data());
  double z = dbgrid->getSimvar(LOC_SIMU,iech,isimu,0,icase,nbsimu,1);
  if (! FFFF(z))
  {
    double d2 = 0.;
    grid_to_point(dbgrid,_ind2.data(),NULL,xyz0.data());
    for (int idim=0; idim<ndim; idim++)
    {
      double delta = _xyz[idim] - xyz0[idim];
      d2 += delta * delta;
    }
    if (d2 <= 0.) return(z);
    top += z / d2;
    bot += 1./ d2;
  }

  /* Second point */
  _ind2[0] += 1;
  iech = db_index_grid_to_sample(dbgrid,_ind2.data());
  z = dbgrid->getSimvar(LOC_SIMU,iech,isimu,0,icase,nbsimu,1);
  if (! FFFF(z))
  {
    double d2 = 0.;
    grid_to_point(dbgrid,_ind2.data(),NULL,xyz0.data());
    for (int idim=0; idim<ndim; idim++)
    {
      double delta = _xyz[idim] - xyz0[idim];
      d2 += delta * delta;
    }
    if (d2 <= 0.) return(z);
    top += z / d2;
    bot += 1./ d2;
  }

  /* Third point */
  _ind2[1] += 1;
  iech = db_index_grid_to_sample(dbgrid,_ind2.data());
  z = dbgrid->getSimvar(LOC_SIMU,iech,isimu,0,icase,nbsimu,1);
  if (! FFFF(z))
  {
    double d2 = 0.;
    grid_to_point(dbgrid,_ind2.data(),NULL,xyz0.data());
    for (int idim=0; idim<ndim; idim++)
    {
      double delta = _xyz[idim] - xyz0[idim];
      d2 += delta * delta;
    }
    if (d2 <= 0.) return(z);
    top += z / d2;
    bot += 1./ d2;
  }

  /* Fourth point */
  _ind2[0] -= 1;
  iech = db_index_grid_to_sample(dbgrid,_ind2.data());
  z = dbgrid->getSimvar(LOC_SIMU,iech,isimu,0,icase,nbsimu,1);
  if (! FFFF(z))
  {
    double d2 = 0.;
    grid_to_point(dbgrid,_ind2.data(),NULL,xyz0.data());
    for (int idim=0; idim<ndim; idim++)
    {
      double delta = _xyz[idim] - xyz0[idim];
      d2 += delta * delta;
    }
    if (d2 <= 0.) return(z);
    top += z / d2;
    bot += 1./ d2;
  }

  /* Final interpolation */
  _ind2[1] -= 1;
  z = (bot != 0) ? top / bot : TEST;
  return(z);
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
int RuleShadow::gaus2facData(PropDef *propdef,
                             Db *dbin,
                             Db *dbout,
                             int *flag_used,
                             int ipgs,
                             int isimu,
                             int nbsimu)
{
  double y[2],facies,t1min,t1max,t2min,t2max,sh_dsup,sh_down;

  /* Initializations */

  check_mandatory_attribute("rule_gaus2fac_data_shadow",dbin,LOC_GAUSFAC);

  /* Processing the translation */

  for (int iech=0; iech<dbin->getSampleNumber(); iech++)
  {
    if (! dbin->isActive(iech)) continue;

    /* Initializations */

    facies = TEST;
    for (int igrf=0; igrf<2; igrf++) y[igrf] = TEST;

    if (rule_thresh_define_shadow(propdef,dbin,this,ITEST,
                                  iech,isimu,nbsimu,1,
                                  &t1min,&t1max,&t2min,&t2max,
                                  &sh_dsup,&sh_down)) return 1;

    for (int igrf=0; igrf<2; igrf++)
    {
      int icase = get_rank_from_propdef(propdef,ipgs,igrf);
      y[igrf] = (flag_used[igrf]) ?
        dbin->getSimvar(LOC_GAUSFAC,iech,isimu,0,icase,nbsimu,1) : 0.;
    }
    facies = getFaciesFromGaussian(y[0],y[1]);

    /* Combine the underlying GRFs to derive Facies*/

    dbin->setSimvar(LOC_FACIES,iech,isimu,0,ipgs,nbsimu,1,facies);
  }
  return 0;
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
int RuleShadow::gaus2facResult(PropDef *propdef,
                               Db *dbout,
                               int *flag_used,
                               int ipgs,
                               int isimu,
                               int nbsimu)
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
  dinc  = getIncr();
  nstep = (int) floor(getDMax() / dinc);
  dy    = dinc * getTgte();
  for (idim=0; idim<ndim; idim++) del[idim] = dinc * getShift(idim);

  /* Processing the translation */

  for (iech=0; iech<dbout->getSampleNumber(); iech++)
  {
    if (! dbout->isActive(iech)) continue;

    /* Initializations */

    facies = TEST;
    for (igrf=0; igrf<2; igrf++) y[igrf] = TEST;

    y[0] = dbout->getSimvar(LOC_SIMU,iech,isimu,0,icase,nbsimu,1);
    if (FFFF(y[0])) break;
    if (rule_thresh_define_shadow(propdef,dbout,this,SHADOW_WATER,
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
        ys = st_grid_eval(dbout,isimu,icase,nbsimu,xyz);
        if (FFFF(ys)) continue;
        jech  = db_index_grid_to_sample(dbout,ind2.data());
        if (rule_thresh_define_shadow(propdef,dbout,this,SHADOW_WATER,
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
  }

  /* Set the error return code */

  error = 0;

label_end:
  del  = db_vector_free(del);
  return(error);
}

/****************************************************************************/
/*!
**  Set the bounds and possibly add replicates (Shadow)
**
** \return  Error return code
**
** \param[in]  propdef    PropDef structure
** \param[in]  dbin       Db structure
** \param[in]  dbout      Db grid structure
** \param[in]  isimu      Rank of the simulation (PROCESS_CONDITIONAL)
** \param[in]  igrf       Rank of the GRF
** \param[in]  ipgs       Rank of the GS
** \param[in]  nbsimu     Number of simulations (PROCESS_CONDITIONAL)
**
*****************************************************************************/
int RuleShadow::evaluateBounds(PropDef *propdef,
                               Db *dbin,
                               Db *dbout,
                               int isimu,
                               int igrf,
                               int ipgs,
                               int nbsimu)
{
  int    iech,jech,nadd,nech,idim,facies,nstep,istep,valid;
  double dist,t1min,t1max,t2min,t2max,s1min,s1max,s2min,s2max;
  double dinc,seuil,alea,sh_dsup,sh_down,yc_down,dval;

  /* Initializations */

  if (dbin == (Db *) NULL) return(0);
  nadd = 0;
  nech = dbin->getSampleNumber();
  dist = 0.;
  dinc  = getIncr();
  nstep = (int)floor(getDMax() / dinc);

  /* Case of the shadow */

  if (igrf == 1) return(0);

  /* Loop on the data */
  for (iech=0; iech<nech; iech++)
  {
    /* Convert the proportions into thresholds for data point */
    if (! dbin->isActive(iech)) continue;
    if (! point_inside_grid(dbin,iech,dbout)) continue;
    facies = (int) dbin->getVariable(iech,0);
    if (rule_thresh_define_shadow(propdef,dbin,this,facies,iech,isimu,nbsimu,1,
                                  &t1min,&t1max,&t2min,&t2max,
                                  &sh_dsup,&sh_down)) return(1);
    yc_down = sh_down;
    dbin->setLowerBound(iech,get_rank_from_propdef(propdef,ipgs,igrf),t1min);
    dbin->setUpperBound(iech,get_rank_from_propdef(propdef,ipgs,igrf),t1max);

    /* The data belongs to the island, no replicate */

    if (facies == SHADOW_ISLAND) continue;

    /* In the case of data belonging to the SHADOW */
    /* Generate one replicate in ISLAND at uniform distance upstream */

    if (facies == SHADOW_SHADOW)
    {
      /* Add one replicate */
      jech = dbin->addSamples(1,0.);
      if (jech < 0) return(1);

      /* Set the coordinates of the replicate */
      /* - at a point where proportions are known */
      /* - after truncation, the point can create shadow at target */

      alea  = 1;
      valid = 0;
      while (! valid)
      {
        dist = 0.;
        alea = law_uniform(0.,1.);
        for (idim=0; idim<dbin->getNDim(); idim++)
        {
          dval = alea * getDMax() * getShift(idim);
          dbin->setCoordinate(jech,idim,dbin->getCoordinate(iech,idim) - dval);
          dist += dval * dval;
        }
        dist = sqrt(dist);

        /* Can the replicate be added */

        if (replicateInvalid(dbin,dbout,jech))
        {
          dbin->deleteSample(jech);
          continue;
        }

        /* Get proportion at the tentative replicate */
        if (rule_thresh_define_shadow(propdef,dbin,this,facies,
                                      jech,isimu,nbsimu,1,
                                      &s1min,&s1max,&s2min,&s2max,
                                      &sh_dsup,&sh_down))
        {
          dbin->deleteSample(jech);
          return(1);
        }
        seuil = t1max - yc_down + dist * getTgte();
        if (seuil > s1max + sh_dsup) continue;
        valid = 1;
      }

      /* Set the attributes of the replicate */
      dbin->setVariable(jech,0,SHADOW_ISLAND);
      dbin->setLowerBound(jech,get_rank_from_propdef(propdef,ipgs,igrf),
               MAX(seuil,s1max));
      dbin->setUpperBound(jech,get_rank_from_propdef(propdef,ipgs,igrf),
               THRESH_SUP);
      nadd++;
    }

    /* In the case of data belonging to the WATER            */
    /* Generate series of replicates (maximum number: nstep) */
    /* whose "elevation" which will not create any shadow    */

    if (facies == SHADOW_WATER)
    {
      /* Loop on the replicates */
      for (istep=1; istep<=nstep; istep++)
      {
        jech = dbin->addSamples(1,0.);
        if (jech < 0) return(1);

        /* Set the coordinates of the replicate */
        dist = 0.;
        for (idim=0; idim<dbin->getNDim(); idim++)
        {
          dval = dinc * getShift(idim) * istep;
          dbin->setCoordinate(jech,idim,dbin->getCoordinate(iech,idim) - dval);
          dist += dval * dval;
        }
        dist = sqrt(dist);

        /* Can the replicate be added */
        if (replicateInvalid(dbin,dbout,jech))
        {
          dbin->deleteSample(jech);
          continue;
        }

        /* Get proportion at the tentative replicate */
        if (rule_thresh_define_shadow(propdef,dbin,this,facies,
                                      jech,isimu,nbsimu,1,
                                      &s1min,&s1max,&s2min,&s2max,
                                      &sh_dsup,&sh_down))
        {
          dbin->deleteSample(jech);
          return(1);
        }

        /* Set the attributes of the replicate */
        seuil = t1max - yc_down + dist * getTgte();
        if (seuil > s1max + sh_dsup)
        {
          /* The replicate is not necessary */
          dbin->deleteSample(jech);
          continue;
        }

        dbin->setVariable(jech,0,SHADOW_IDLE);
        dbin->setLowerBound(jech,get_rank_from_propdef(propdef,ipgs,igrf),
                 THRESH_INF);
        dbin->setUpperBound(jech,get_rank_from_propdef(propdef,ipgs,igrf),
                 MAX(seuil,s1max));
        nadd++;
      }
    }
    jech++;
  }

  if (igrf == 0 && nadd > 0)
  {
    message("Initial count of data = %d\n",nech);
    message("Number of replicates  = %d\n",nadd);
  }
  return(0);
}
