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
#include "Covariances/CovAniso.hpp"
#include "Basic/Law.hpp"
#include "geoslib_e.h"
#include "geoslib_old_f.h"

/*! \cond */
#define QUANTUM 1000
#define DISCRET(idisc)     (GV_PI * (0.5 + (idisc)) / ((double) ndisc))
#define IPTR(ix,iy)        ((iy) * nx + (ix))
/*! \endcond */

/*****************************************************************************/
/*!
**  Simulates the discrete distribution on 0,1,...,n
**
** \param[in]  nfreq  Number of frequencies
** \param[in]  freqs  Array of frequencies (which add up to 1)
**
*****************************************************************************/
static int st_gdiscrete(int     nfreq,
                        double *freqs)
{
  double partvec,u;

  u = law_uniform(0.,1.);

  partvec = 0.;
  for (int ifreq=0; ifreq<nfreq; ifreq++)
  {
    partvec += freqs[ifreq];
    if (u < partvec) return(ifreq);
  }
  return(nfreq-1);
}

/*****************************************************************************/
/*!
**  Allocate or Reallocate memory for the frequencies of the spectrum
**
** \returns The allocated array
**
** \param[in]  freqs   Array provided in input
** \param[in]  nfreq   Current number of frequencies
**
** \param[out] nquloc  Number of quantums
**
*****************************************************************************/
static double *st_freq_array_realloc(double *freqs,
                                     int     nfreq,
                                     int    *nquloc)
{
  int nquant,newsize;

  nquant = *nquloc;
  if (nfreq >= nquant * QUANTUM)
  {
    nquant++;
    newsize = nquant * QUANTUM;
    if (freqs == nullptr)
      freqs = (double *) mem_alloc(sizeof(double) * newsize,0);
    else {
      freqs = (double *) mem_realloc((char *) freqs,sizeof(double) * newsize,0);
    }
    if (freqs == nullptr) return(freqs);
    *nquloc = nquant;
  }
  return(freqs);
}

/*****************************************************************************/
/*!
**  Generate the spectrum for Chentsov
**
** \returns The array 'freqs' or NULL
**
** \param[in]  tol     Tolerance on the sum of probabilities
** \param[in]  nfmax   Maximum number of frequencies (or < 0)
**
** \param[out] nfreq   Number of frequencies
**
*****************************************************************************/
static double *st_spectrum_chentsov(double  tol,
                                    int     nfmax,
                                    int    *nfreq)
{
  double *freqs,ratio,total;
  int     nquant,ifreq;

  /* Initializations */
  
  nquant = ifreq = 0;
  total  = 0.;
  freqs  = nullptr;
  
  /* Core allocation */
  
  freqs = st_freq_array_realloc(freqs,ifreq,&nquant);
  if (freqs == nullptr) return(freqs);
  
  /* Loop on the spectrum items */
  
  freqs[ifreq] = 0.;            // ifreq = 0
  total += freqs[ifreq];
  ifreq++;

  freqs[ifreq] = 0.75;          // ifreq = 1
  total += freqs[ifreq];
  ifreq++;
  
  while (1)
  {
    freqs = st_freq_array_realloc(freqs,ifreq,&nquant);
    if (freqs == nullptr) return(freqs);

    freqs[ifreq] = 0.;
    ifreq++;

    freqs = st_freq_array_realloc(freqs,ifreq,&nquant);
    if (freqs == nullptr) return(freqs);

    ratio = ((double) (ifreq - 2.)) / ((double) (ifreq + 1.));
    freqs[ifreq]  = freqs[ifreq-2] * (2.*ifreq+1.) / (2.*ifreq-3.);
    freqs[ifreq] *= ratio * ratio;
    total += freqs[ifreq];
    ifreq++;

    if (ABS(1. - total) < tol) break;
    if (nfmax > 0 && ifreq >= nfmax) break;
  }    

  *nfreq = ifreq;
  freqs = (double *) mem_realloc((char *) freqs,sizeof(double) * ifreq,0);
  return(freqs);
}

/*****************************************************************************/
/*!
**  Generate the spectrum for Exponential
**
** \returns The array 'freqs' or NULL
**
** \param[in]  model   Model (used for its range)
** \param[in]  tol     Tolerance on the sum of probabilities
** \param[in]  nfmax   Maximum number of frequencies (or < 0)
**
** \param[out] nfreq   Number of frequencies
**
*****************************************************************************/
static double *st_spectrum_exponential(Model  *model,
                                       double  tol,
                                       int     nfmax,
                                       int    *nfreq)
{
  double *freqs,fcs,fcs2,r1,r2,expfc,total;
  int     nquant,ifreq;
  
  /* Initializations */

  nquant = ifreq = 0;
  total  = 0.;
  freqs  = nullptr;
  fcs    = 1. / model->getCova(0)->getTheoretical();
  fcs2   = fcs * fcs;
  expfc  = exp(-fcs * GV_PI);
  
  /* Core allocation */
  
  freqs = st_freq_array_realloc(freqs,ifreq,&nquant);
  if (freqs == nullptr) return(freqs);
  
  freqs[ifreq] = (1. + expfc) / (fcs2 + 1.) / 2.;
  if (freqs[ifreq] < 0.) freqs[ifreq] = 0.;
  total += freqs[ifreq];
  ifreq++;
  
  freqs[ifreq] = 3. * (1. - expfc) / (fcs2 + 4.) / 2.;
  if (freqs[ifreq] < 0.) freqs[ifreq] = 0.;
  total += freqs[ifreq];
  ifreq++;

  /* Loop on the spectrum items */
  
  while (1)
  {
    freqs = st_freq_array_realloc(freqs,ifreq,&nquant);
    if (freqs == nullptr) return(freqs);

    r1 = ifreq + 1.;
    r2 = ifreq - 2.;
    freqs[ifreq] = freqs[ifreq-2] * (2.*ifreq+1.) / (2.*ifreq-3.);
    freqs[ifreq] *=(fcs2 + r2*r2) / (fcs2 + r1*r1);
    total += freqs[ifreq];
    ifreq++;

    if (ABS(1. - total) < tol) break;
    if (nfmax > 0 && ifreq >= nfmax) break;
  }
  
  *nfreq = ifreq;
  freqs = (double *) mem_realloc((char *) freqs,sizeof(double) * ifreq,0);
  return(freqs);
}

/*****************************************************************************/
/*!
**  Check the degrees and orders
**
** \return Error code
**
** \param[in]  nfreq   Number of frequencies
** \param[in]  nbf     Number of basic functions
** \param[in]  degree  Array of degrees
** \param[in]  order   Array of orders
** \param[in]  verbose Verbose flag
**
*****************************************************************************/
static int st_check_degree_order(int  nfreq,
                                 int  nbf,
                                 int *degree,
                                 int *order,
                                 int  verbose)
{
  int degmax,ordmin,ordmax;

  degmax =  0;
  ordmin =  nfreq;
  ordmax = -nfreq;

  for (int ibf=0; ibf<nbf; ibf++)
  {
    if (degree[ibf] > degmax) degmax = degree[ibf];
    if (order[ibf]  < ordmin) ordmin = order[ibf];
    if (order[ibf]  > ordmax) ordmax = order[ibf];
    if (order[ibf]  < -degree[ibf] ||
        order[ibf]  > +degree[ibf])
    {
      messerr("Order(%d) must lie in [-degree;+degree] where degree=%d",
              order[ibf],degree[ibf]);
      return(1);
    }
  }

  if (verbose)
  {
    message("Maximum degree            = %d\n",degmax);
    message("Minimum order             = %d\n",ordmin);
    message("Maximum order             = %d\n",ordmax);
  }
  return(0);
}

/*****************************************************************************/
/*!
**  Normalize the spectrum
**
** \param[in]  verbose Verbose flag
** \param[in]  nfreq   Number of frequencies
**
** \param[out] freqs   Array of frequencies
**
*****************************************************************************/
static void st_spectrum_normalize(int     verbose,
                                  int     nfreq,
                                  double *freqs)
{
  double totpos,totneg;
  
  totpos = totneg = 0.;
  for (int ifreq=0; ifreq<nfreq; ifreq++)
  {
    if (freqs[ifreq] < 0)
    {
      totneg -= freqs[ifreq];
      freqs[ifreq] = 0.;
    }
    else
    {
      totpos += freqs[ifreq];
    }
  }

  for (int ifreq=0; ifreq<nfreq; ifreq++) freqs[ifreq] /= totpos;

  /* Printout (optional) */

  if (verbose)
  {
    message("Cumulated Spectrum        = %lf\n",totpos);
    message("Sum of negative weights   = %lf\n",totneg);
  }
}

/*****************************************************************************/
/*!
**  Generate the spectrum for any covariance
**
** \returns The array 'freqs' or NULL
**
** \param[in]  model   Model (defined in Euclidean space) to be used
** \param[in]  ndisc   Discretization of the intergral
** \param[in]  tol     Tolerance on the sum of probabilities
** \param[in]  nfmax   Maximum number of frequencies (or < 0)
**
** \param[out] nfreq   Number of frequencies
**
*****************************************************************************/
static double *st_spectrum_any(Model  *model,
                               int     ndisc,
                               double  tol,
                               int     nfmax,
                               int    *nfreq)
{
  double *freqs,*covs,an,alpha,ca,sina,cosa,total,dincr;
  int     nquant,ifreq;
  VectorDouble dd;

  /* Initializations */

  dd.resize(2);
  dd[0]  = dd[1] = 0.;
  nquant = ifreq = 0;
  freqs  = covs = nullptr;
  dincr  = GV_PI / ((double) ndisc);

  /* Core allocation */
  
  covs  = (double *) mem_alloc(sizeof(double) * ndisc,0);
  if (covs  == nullptr) return(freqs);
  
  /* Calculate the discretized covariance values */

  for (int idisc=0; idisc<ndisc; idisc++)
  {
    alpha = DISCRET(idisc);
    dd[0] = 2. * sin(alpha / 2.);
    ca    = 0.;
    for (int icova=0; icova<model->getCovaNumber(); icova++)
      ca += model_calcul_basic(model,icova,ECalcMember::LHS,dd);
    covs[idisc] = ca;
  }
  
  /* Loop on the frequencies */

  total = 0.;
  while(1)
  {
    freqs = st_freq_array_realloc(freqs,ifreq,&nquant);
    if (freqs == nullptr) return(freqs);

    /* Discretization of the frequency item */
    
    an = 0.;
    for (int idisc=0; idisc<ndisc; idisc++)
    {
      alpha = DISCRET(idisc);
      cosa  = cos(alpha);
      sina  = sin(alpha);
      an   += covs[idisc] * sina * ut_legendre(1,ifreq,cosa);
    }
    freqs[ifreq] = an * dincr * sqrt((2. * ifreq + 1) / 2.);
    total += freqs[ifreq];
    ifreq++;
    
    if (total > 1.)
    {
      ifreq--;
      break;
    }
    if (ABS(1. - total) < tol) break;
    if (nfmax > 0 && ifreq >= nfmax) break;
  }

  *nfreq = ifreq;
  freqs = (double *) mem_realloc((char *) freqs,sizeof(double) * ifreq,0);

  // Core deallocation
  
  covs = (double *) mem_free((char *) covs);

  return(freqs);
}

/*****************************************************************************/
/*!
**  Simulates the random function on the sphere
**
** \param[in]  db          Data base containing the coordinates of target points
**                         These coordinates must be expressed in long/lat
** \param[in]  model       Model (defined in Euclidean space) to be used
** \param[in]  seed        Seed for random number generation
** \param[in]  special     Special type of Model
**                         0 : Standard
**                         1 : Chentsov
**                         2 : Exponential
** \param[in]  nbf         Number of basic functions
** \param[in]  nfmax       Maximum number of frequencies (or <0)
** \param[in]  verbose     Verbose flag
** \param[in]  flag_test   1 to perform the test of basic functions
** \param[in]  test_degree Testing Degree
** \param[in]  test_order  Testing Order
** \param[in]  test_phase  Testing Phase
**
** \remarks: Some variables can be defined using the keypair mechanism
** \remarks: Simsph_Shunt is set:
** \remarks:              0 for complete performance
** \remarks:              1 for shortcut after normalized Spectrum calculation
** \remarks:              2 for shortcut after drawing degree and order arrays
** \remarks: Simsph_Maximum_Degree: to define the Maximum Degree 
** \remarks: Simsph_Fixed_Degree: to define the Maximum Degree 
**
*****************************************************************************/
GEOSLIB_API int simsph_f(Db    *db,
                       Model *model,
                       int    seed,
                       int    special,
                       int    nbf,
                       int    nfmax,
                       int    verbose,
                       int    flag_test,
                       int    test_degree,
                       int    test_order,
                       double test_phase)
{
  int    *degree,*order,flag_sphere;
  int     iptr,error,ndisc,nfreq,degmax,shunt,nech,nx,ny,ecr,ntot;
  double *freqs,*phase;
  double  theta,phi,t1,t2,val,tol,degree_loc,order_loc,phase_loc;

  /* Initializations */

  error  = 1;
  iptr   = nfreq = 0;
  degree = order = nullptr;
  freqs  = phase = nullptr;
  ndisc  = (int) get_keypone("Simsph_Ndisc",360);
  shunt  = (int) get_keypone("Simsph_Shunt",0);
  tol    = get_keypone("Simsph_Spectrum_Tolerance",1.e-5);
  nx     = db->getNX(0);
  ny     = db->getNX(1);
  nech   = db->getSampleNumber();
  if (flag_test) nbf = 1;

  /* Preliminary checks */

  law_set_random_seed(seed);
  variety_query(&flag_sphere);
  if (! flag_sphere)
  {
    messerr("The Spherical Simulation is restricted to Spherical coordinates");
    goto label_end;
  }
  if (! ut_is_legendre_defined()) goto label_end;
  if (db->getNDim() != 2)
  {
    messerr("The Simulation on Sphere is restricted to 2-D case");
    goto label_end;
  }
  for (int icova=0; icova<model->getCovaNumber(); icova++)
  {
    if (model->getCova(icova)->getFlagAniso())
    {
      messerr("Only Isotropic Models may be used for Spherical Simulations");
      goto label_end;
    }
  }

  /* Create the new variable in the Data base */

  iptr = db->addFields(1,0.,String(),ELoc::SIMU);

  /* Core allocation */

  phase  = (double *) mem_alloc(sizeof(double) * nbf,0);
  if (phase  == nullptr) goto label_end;
  degree = (int    *) mem_alloc(sizeof(int)    * nbf,0);
  if (degree == nullptr) goto label_end;
  order  = (int    *) mem_alloc(sizeof(int)    * nbf,0);
  if (order  == nullptr) goto label_end;

  /* Define the spectrum */

  if (special == 1)
    freqs = st_spectrum_chentsov(tol,nfmax,&nfreq);
  else if (special == 2)
    freqs = st_spectrum_exponential(model,tol,nfmax,&nfreq);
  else
    freqs = st_spectrum_any(model,ndisc,tol,nfmax,&nfreq);
  if (freqs == nullptr) goto label_end;
  set_keypair("Simsph_Spectrum_Frequencies",1,nfreq,1,freqs);

  /* Optional printout */

  if (verbose)
  {
    mestitle(1,"Generation of Covariance spectrum in Spherical Coordinates");
    if (special == 1)
      message("Using Chentsov construction\n");
    else if (special == 2)
      message("Using Exponential construction\n");
    else
      message("Number of discretization  = %d\n",ndisc);
    message("Spectrum Tolerance        = %lg\n",tol);
    message("Random generation seed    = %d\n",law_get_random_seed());
    message("Number of frequencies     = %d\n",nfreq);
    message("Number of basic functions = %d\n",nbf);
  }
  st_spectrum_normalize(verbose,nfreq,freqs);
  if (shunt == 1) goto label_short;

  /* Get the model ingredients (generated by separate flows in order to  */
  /* avoid intermeshing dependencies) */

  if (flag_test)
  {
    degree[0] = test_degree;
    order[0]  = test_order;
    phase[0]  = test_phase;
  }
  else
  {
    degmax = (int) get_keypone("Simsph_Maximum_Degree",-1.);
    for (int ibf=0; ibf<nbf; ibf++)
    {
      degree[ibf] = st_gdiscrete(nfreq,freqs);
      if (degmax > 0) degree[ibf] = MIN(degmax,degree[ibf]);
    }
    for (int ibf=0; ibf<nbf; ibf++)
      order[ibf]  = law_int_uniform(-degree[ibf],degree[ibf]);
    for (int ibf=0; ibf<nbf; ibf++)
      phase[ibf]  = law_uniform(0.,2. * GV_PI);
  }
  if (st_check_degree_order(nfreq,nbf,degree,order,verbose)) goto label_end;

  /* Saving options */
  
  set_keypair_int("Simsph_Array_Degree",1,nbf,1,degree);
  set_keypair_int("Simsph_Array_Order" ,1,nbf,1,order);
  set_keypair    ("Simsph_Array_Phase" ,1,nbf,1,phase);
  if (shunt == 2) goto label_short;
    
  /* Loop on the samples of the Data Base */
  /* We benefit in writing this as a double loop on coordinates */
  /* as some calculations can be factorized as they only concern latitude */

  ecr  = 0;
  ntot = nbf * nech;
  for (int iy=0; iy<ny; iy++)
  {
    theta  = ut_deg2rad(db->getCoordinate(IPTR(0,iy),1) + 90.); // Latitude[-90,90]
    for (int ibf=0; ibf<nbf; ibf++)
    {
      degree_loc = degree[ibf];
      order_loc  = order[ibf];
      phase_loc  = phase[ibf];
      t1 = ut_flegendre(1,(int) degree_loc,(int) order_loc,theta);
      for (int ix=0; ix<nx; ix++,ecr++)
      {
        int jech = IPTR(ix,iy);
        mes_process("Simulation on Sphere",ntot,ecr);
        if (! db->isActive(jech)) continue;
        phi = ut_deg2rad(db->getCoordinate(jech,0));       // Longitude [0,360]
        t2  = cos(phi * order_loc + phase_loc);
        db->updArray(jech,iptr,0,t1*t2);
      }
    }
  }
  
  /* Final normalization */
    
  val = 2. / sqrt((double) nbf);
  for (int iech=0; iech<nech; iech++)
  {
    if (db->isActive(iech)) db->updArray(iech,iptr,3,val);
  }

  /* Set the error return code */

label_short:
  error = 0;

label_end:
  if (error) db->deleteFieldByLocator(ELoc::SIMU);
  degree = (int    *) mem_free((char *) degree);
  order  = (int    *) mem_free((char *) order);
  freqs  = (double *) mem_free((char *) freqs);
  phase  = (double *) mem_free((char *) phase);
  return(error);
}
