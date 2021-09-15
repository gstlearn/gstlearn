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
#include "Gibbs/AGibbs.hpp"
#include "Basic/Utilities.hpp"
#include "Basic/AStringable.hpp"
#include "geoslib_old_f.h"
#include "geoslib_define.h"
#include "geoslib_enum.h"

AGibbs::AGibbs()
    : _npgs(1),
      _ngrf(1),
      _nbsimu(1),
      _nburn(1),
      _niter(1),
      _flagOrder(true),
      _flagCategory(true),
      _rho(1.),
      _sqr(0.),
      _eps(EPSILON3)
{
}

AGibbs::AGibbs(int npgs, int ngrf, int nbsimu, int nburn, int niter,
               int flag_order, bool flag_category, double rho, double eps)
    : _npgs(1),
      _ngrf(1),
      _nbsimu(1),
      _nburn(1),
      _niter(1),
      _flagOrder(false),
      _flagCategory(false),
      _rho(1.),
      _sqr(0.),
      _eps(eps)
{
  init(npgs, ngrf, nbsimu, nburn, niter, flag_order, flag_category, rho, eps);
  _sqr = sqrt(1. - _rho * rho);
}

AGibbs::AGibbs(const AGibbs &r)
    : _npgs(r._npgs),
      _ngrf(r._ngrf),
      _nbsimu(r._nbsimu),
      _nburn(r._nburn),
      _niter(r._niter),
      _flagOrder(r._flagOrder),
      _flagCategory(r._flagCategory),
      _rho(r._rho),
      _sqr(r._sqr),
      _eps(r._eps)
{
}

AGibbs& AGibbs::operator=(const AGibbs &r)
{
  if (this != &r)
  {
    _npgs = r._npgs;
    _ngrf = r._ngrf;
    _nbsimu = r._nbsimu;
    _nburn = r._nburn;
    _niter = r._niter;
    _flagOrder = r._flagOrder;
    _flagCategory = r._flagCategory;
    _rho = r._rho;
    _sqr = r._sqr;
    _eps = r._eps;
  }
  return *this;
}

AGibbs::~AGibbs()
{
}

void AGibbs::init(int npgs,
                  int ngrf,
                  int nbsimu,
                  int nburn,
                  int niter,
                  int flag_order,
                  bool flag_category,
                  double rho,
                  double eps)
{
  _npgs = npgs;
  _ngrf = ngrf;
  _nbsimu = nbsimu;
  _nburn = nburn;
  _niter = niter;
  _flagOrder = flag_order;
  _flagCategory = flag_category;
  _rho = rho;
  _eps = eps;

  _sqr = sqrt(1. - _rho * rho);
}

int AGibbs::getRank(int ipgs, int igrf) const
{
  if (ipgs <= 0)
    return(igrf);
  else
    return(getNgrf() + igrf);
}

/****************************************************************************/
/*!
**  Check for the presence of mandatory attributes
**
** \param[in]  method  Name of the method
** \param[in]  db      Db structure
** \param[in]  locatorType  Mandatory attribute type
**
*****************************************************************************/
int AGibbs::_checkMandatoryAttribute(const String& method,
                                       Db *db,
                                       ENUM_LOCS locatorType)
{
  if (get_LOCATOR_NITEM(db,locatorType) <= 0)
  {
    messerr("%s : Attributes %d are mandatory",method.c_str(),locatorType);
    return 1;
  }
  return 0;
}

/****************************************************************************/
/*!
**  Check/Show the facies against gaussian at wells
**
** \return Error return code
**
** \param[in]  db         Db structure
** \param[in]  model      Model structure
** \param[in]  isimu      Rank of the simulation
** \param[in]  ipgs       Rank of the GS
** \param[in]  igrf       Rank of the bounds (starting from 0)
**
** \remark Attributes LOC_GAUSFAC are mandatory
**
*****************************************************************************/
int AGibbs::checkGibbs(Db *db, Model *model, int isimu, int ipgs, int igrf)
{
  int nech,iech,number,icase,icase0;
  double gaus,vmin,vmax;

  /* Initializations */

  _checkMandatoryAttribute("st_check_gibbs",db,LOC_GAUSFAC);
  number = 0;
  nech   = db->getSampleNumber();
  icase  = getRank(ipgs,igrf);
  icase0 = getRank(ipgs,0);
  mestitle(1,"Checking gaussian values from Gibbs vs. bounds (PGS=%d GRF=%d Simu=%d)",
           ipgs+1,igrf+1,isimu+1);

  /* Loop on the data */

  for (iech=0; iech<nech; iech++)
  {
    if (! db->isActive(iech)) continue;

    /* Read the bounds */

    vmin = db->getLowerBound(iech,icase);
    vmax = db->getUpperBound(iech,icase);
    if (FFFF(vmin)) vmin = -1.e30;
    if (FFFF(vmax)) vmax =  1.e30;

    /* Read the gaussian value */

    gaus = db->getSimvar(LOC_GAUSFAC,iech,isimu,0,icase,_nbsimu,1);
    if (igrf == 1)
      gaus = _sqr * gaus + _rho *
        db->getSimvar(LOC_GAUSFAC,iech,isimu,0,icase0,_nbsimu,1);

    /* Check inconsistency */

    if ((! FFFF(vmin) && gaus < vmin) ||
        (! FFFF(vmax) && gaus > vmax))
    {
      message("- Sample (#%d):",iech+1);
      message(" Simu#%d of Y%d=%lf",isimu+1,igrf+1,gaus);
      message(" does not lie within [");
      if (FFFF(vmin))
        message("NA,");
      else
        message("%lf",vmin);
      message(";");
      if (FFFF(vmax))
        message("NA");
      else
        message("%lf",vmax);
      message("]\n");

      number++;
    }
  }

  if (number <= 0) message("No problem found\n");

  return number;
}

/****************************************************************************/
/*!
**  Check if bounds are corrects
**
** \return  Error code: 1 if bounds are incorrect
**
** \param[in]     db          Db structure
** \param[in]     iech0       Rank of the sample
** \param[in]     data        Data value (of TEST)
** \param[in]     vmin        Minimum bound (or TEST)
** \param[in]     vmax        Maximum bound (or TEST)
** \param[in]     iech        Rank of the previous sample
** \param[in]     value       Previous value (or TEST)
** \param[in]     vemin       Minimum bound for previous sample (or TEST)
** \param[in]     vemax       Maximum bound for previous sample (or TEST)
**
** \remarks If an error occurs, the error message is printed.
** \remarks If the previous sample ('iech') is undefined, no comparison
** \remarks is printed.
** \remarks Bounds 'vmin' and 'vmax' are modified in presence of hard data
**
*****************************************************************************/
int AGibbs::_boundsCheck(Db *db,
                         int iech0,
                         double data,
                         double *vmin,
                         double *vmax,
                         int iech,
                         double value,
                         double vemin,
                         double vemax)
{
  int flag_err_bnd,flag_err_min,flag_err_max;
  double vlmin,vlmax;

  /* Initializations */

  flag_err_min = 0;
  flag_err_max = 0;
  flag_err_bnd = 0;

  // Check the bound validity

  flag_err_bnd = (! FFFF(*vmin) && ! FFFF(*vmax) && (*vmin) > (*vmax));

  // Check against the data

  vlmin = *vmin;
  vlmax = *vmax;
  if (! FFFF(data))
  {

    /* Case where the data value is defined */

    flag_err_min = (! FFFF(*vmin) && data < (*vmin));
    flag_err_max = (! FFFF(*vmax) && data > (*vmax));
    db->setLowerBound(iech0,0,data);
    db->setUpperBound(iech0,0,data);
    *vmin = *vmax = data;
  }

  if (! (flag_err_min || flag_err_max || flag_err_bnd)) return(0);

  // Print the error message

  messerr("Sample %d",iech0+1);
  if (flag_err_bnd)
    messerr("Bounds are wrongly ordered: Vmin(%lf) > Vmax(%lf)",
            vlmin,vlmax);
  if (flag_err_min || flag_err_max)
    messerr("Data (%lf) does not lie in [%lf ; %lf]",data,vlmin,vlmax);

  if (! IFFFF(iech))
  {
    messerr("Compared to sample %d",iech+1);
    if (!FFFF(value))
      messerr("- Value = %lf",value);
    if (!FFFF(vemin))
      messerr("- Lower Bound = %lf",vemin);
    if (!FFFF(vemax))
      messerr("- Upper Bound = %lf",vemax);
  }
  if (! FFFF(data))
  {
    messerr("... Data information superseeds the bounds");
    return(0);
  }

  return(1);
}

/****************************************************************************/
/*!
**  Correct the bounds according to the order relationship
**
** \return  Error code: 1 if there is no solution; 0 otherwise
**
** \param[in]  flag_category 1 for categorical; 0 for continuous
** \param[in]  flag_order    Order relationship
** \li                        1 if the ascending order must be honored
** \li                       -1 if the descending order must be honored
** \li                        0 if no order relationship must be honored
** \param[in]  db            Db structure
** \param[in]  iech0         Rank of the sample
** \param[in]  ivar          Rank of the variable
** \param[in]  icase         Rank of the GRF / PGS
** \param[in]  nvar          Number of variables
** \param[in]  vlmin_arg     Input minimum bound
** \param[in]  vlmax_arg     Input maximum bound
**
** \param[out]  vlmin_arg   Output minimum bound
** \param[out]  vlmax_arg   Output maximum bound
**
** \remark Attributes LOC_GAUSFAC are mandatory
**
*****************************************************************************/
int AGibbs::_correctBoundsOrder(int flag_category,
                                int flag_order,
                                Db *db,
                                int iech0,
                                int ivar,
                                int icase,
                                int nvar,
                                double *vlmin_arg,
                                double *vlmax_arg)
{
  int    iech;
  double vlmin,vlmax,vemin,vemax,vimin,vimax,value,data;

  /* Preliminary check */

  check_mandatory_attribute("st_correct_bounds_order",db,LOC_GAUSFAC);
  iech  = -1;
  value = 0.;
  data  = TEST;
  if (! flag_category) data = db->getVariable(iech0,0);
  vimin = db->getLowerBound(iech0,icase);
  vimax = db->getUpperBound(iech0,icase);
  if (_boundsCheck(db, iech0, data, &vimin, &vimax,
                   ITEST, TEST, TEST, TEST)) return (1);
  vlmin = vimin;
  vlmax = vimax;

  /* Dispatch */

  switch(flag_order)
  {
    case -1:      /* Descending order */
      for (iech=0; iech<iech0; iech++)
      {
        value = db->getSimvar(LOC_GAUSFAC,iech,0,ivar,icase,1,nvar);
        if (!FFFF(value) && (FFFF(vlmax) || value < vlmax)) vlmax = value;

        vemin = db->getLowerBound(iech,icase);
        vemax = db->getUpperBound(iech,icase);
        if (!FFFF(vemax) && (FFFF(vlmax) || vemax < vlmax)) vlmax = vemax;

        if (_boundsCheck(db,iech0,data,&vlmin,&vlmax,
                         iech,value,vemin,vemax)) return(1);
      }
      for (iech=iech0+1; iech<db->getSampleNumber(); iech++)
      {
        value = db->getSimvar(LOC_GAUSFAC,iech,0,ivar,icase,1,nvar);
        if (!FFFF(value) && (FFFF(vlmin) || value > vlmin)) vlmin = value;

        vemin = db->getLowerBound(iech,icase);
        vemax = db->getUpperBound(iech,icase);
        if (!FFFF(vemin) && (FFFF(vlmin) || vemin > vlmin)) vlmin = vemin;

        if (_boundsCheck(db,iech0,data,&vlmin,&vlmax,
                         iech,value,vemin,vemax)) return(1);
      }
      break;

    case 1:     /* Ascending order */
      for (iech=0; iech<iech0; iech++)
      {
        value = db->getSimvar(LOC_GAUSFAC,iech,0,ivar,icase,1,nvar);
        if (!FFFF(value) && (FFFF(vlmin) || value > vlmin)) vlmin = value;

        vemin = db->getLowerBound(iech,icase);
        vemax = db->getUpperBound(iech,icase);
        if (!FFFF(vemin) && (FFFF(vlmin) || vemin > vlmin)) vlmin = vemin;

        if (_boundsCheck(db,iech0,data,&vlmin,&vlmax,
                         iech,value,vemin,vemax)) return(1);
      }
      for (iech=iech0+1; iech<db->getSampleNumber(); iech++)
      {
        value = db->getSimvar(LOC_GAUSFAC,iech,0,ivar,icase,1,nvar);
        if (!FFFF(value) && (FFFF(vlmax) || value < vlmax)) vlmax = value;

        vemin = db->getLowerBound(iech,icase);
        vemax = db->getUpperBound(iech,icase);
        if (!FFFF(vemax) && (FFFF(vlmax) || vemax < vlmax)) vlmax = vemax;

        if (_boundsCheck(db,iech0,data,&vlmin,&vlmax,
                         iech,value,vemin,vemax)) return(1);
      }
      break;

    default:      /* No order relationship */
      break;
  }

  *vlmin_arg = vlmin;
  *vlmax_arg = vlmax;
  return(0);
}

/****************************************************************************/
/*!
**  Print the inequality
**
** \param[in]  db        Db structure
** \param[in]  ifirst    If 0, print the header
** \param[in]  iech      Rank of the sample
** \param[in]  ivar      Rank of the variable
** \param[in]  nfois     Rank of the iteration (<=0 for bootstrap)
** \param[in]  flag_cv   1 to print the convergence criterion
** \param[in]  simval    Simulated value
** \param[in]  vmin      Lower threshold
** \param[in]  vmax      Upper threshold
** \param[in]  mean      Current mean value
** \param[in]  delta     Convergence increment
**
*****************************************************************************/
void AGibbs::_printInequalities(Db *db,
                                int ifirst,
                                int iech,
                                int ivar,
                                int nfois,
                                int flag_cv,
                                double simval,
                                double vmin,
                                double vmax,
                                double mean,
                                double delta)
{
  int flag_min,flag_max,idim;

  /* Initializations */

  flag_min = flag_max = 1;
  if (FFFF(vmin)) flag_min = 0;
  if (FFFF(vmax)) flag_max = 0;

  /* Print the title (first sample) */

  if (ifirst == 0 && nfois >= 0)
    message("Iteration = %d\n",nfois);

  /* Print the simulated value */

  message("Sample (%3d) - Variable (%3d) = %8.4lf in ",iech+1,ivar+1,simval);

  /* Print the bounds */

  if (! flag_min)
    message("[      NA,");
  else
    message("[%8.4lf,",vmin);
  if (! flag_max)
    message("      NA]");
  else
    message("%8.4lf]",vmax);

  /* Print the coordinates */

  message(" at point (");
  for (idim=0; idim<db->getNDim(); idim++)
  {
    if (idim != 0) message(",");
    message("%8.4lf",db->getCoordinate(iech,idim));
  }
  message(")");

  /* Print the mean and the evolution increment */

  if (flag_cv)
    message(" - Mean = %8.4lf (Delta = %5.2lf (percent))",mean,delta);

  message("\n");
}

/****************************************************************************/
/*!
 **  Print the initial status for Gibbs iterations
 **
 ** \param[in]  title       Title of the printout
 ** \param[in]  dbin        Db structure
 ** \param[in]  nvar        Number of variables
 ** \param[in]  nbsimu      Number of simulations
 ** \param[in]  isimu       Rank of the simulation
 ** \param[in]  icase       Case for PGS or GRF
 **
 *****************************************************************************/
void AGibbs::_gibbsInitPrint(const char *title,
                             Db *dbin,
                             int nvar,
                             int nbsimu,
                             int isimu,
                             int icase)
{
  int nech, iech, ivar;
  double simval, vmin, vmax;

  nech = dbin->getSampleNumber();
  mestitle(1, "%s (Simu:%d)", title, isimu + 1);
  for (ivar = 0; ivar < nvar; ivar++)
    for (iech = 0; iech < nech; iech++)
    {
      if (!dbin->isActive(iech)) continue;
      vmin = dbin->getLowerBound(iech, icase);
      vmax = dbin->getUpperBound(iech, icase);
      simval = dbin->getSimvar(LOC_GAUSFAC, iech, isimu, ivar, icase, nbsimu, nvar);
      _printInequalities(dbin, 0, iech, ivar, -1, 0, simval, vmin, vmax, TEST, TEST);
    }
}

/****************************************************************************/
/*!
**  Print the final scores for Gibbs iterations
**
** \param[in]  title       Title of the printout
** \param[in]  dbin        Db structure
** \param[in]  nvar        Number of variables
** \param[in]  isimu       Rank of the simulation
** \param[in]  niter       Total number of iterations
** \param[in]  icase       Case for PGS or GRF
**
*****************************************************************************/
void AGibbs::_gibbsIterPrint(const char *title,
                             Db *dbin,
                             int nvar,
                             int isimu,
                             int niter,
                             int icase)
{
  int    nech,iech,ivar,iecr;
  double simval,vmin,vmax;

  nech = dbin->getSampleNumber();
  mestitle(1,"%s (Simu:%d)",title,isimu+1);
  message("Number of samples              = %d\n",nech);
  message("Number of bootstrap iterations = %d\n",_nburn);
  message("Maximum number of iterations   = %d\n",_niter);
  message("Total number of iterations     = %d\n",niter);
  message("Relative stability criterion   = %5.2lf per cent\n",_eps);

  mestitle(1,"Assignment of Values at data points ");
  iecr = 0;
  for (ivar=0; ivar<nvar; ivar++)
    for (iech=0; iech<nech; iech++)
    {
      if (! dbin->isActive(iech)) continue;
      vmin = dbin->getLowerBound(iech,icase);
      vmax = dbin->getUpperBound(iech,icase);
      simval = dbin->getSimvar(LOC_GAUSFAC,iech,isimu,ivar,icase,
                               getNbsimu(),nvar);
      _printInequalities(dbin,iecr,iech,ivar,-1,0,simval,vmin,vmax,TEST,TEST);
      iecr++;
    }
}
