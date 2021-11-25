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
#include "Basic/EJustify.hpp"
#include "geoslib_e.h"
#include "geoslib_old_f.h"

static int     FLAG_SIMU = 0;
static int     FLAG_NO_VAR_CHECK = 0;
static int     NBGH_INITIALIZED = 0;
static int    *NBGH_ind   = nullptr;
static int    *NBGH_isect = nullptr;
static int    *NBGH_nsect = nullptr;
static double *NBGH_x1    = nullptr;
static double *NBGH_x2    = nullptr;
static double *NBGH_dst   = nullptr;

/****************************************************************************/
/*!
**  Print the information selected in the neighborhood
**
** \param[in]  dbin      input Db structure
** \param[in]  neigh     Neigh structure
** \param[in]  rank      Array of the data ranks
** \li                   -1 if not selected
** \li                   >=0 gives the angular sector in ENeigh::MOVING
** 
*****************************************************************************/
static void st_neigh_print(Db* dbin,
                           Neigh* neigh,
                           int* rank)
{
  int iech,idim,ndim,sel,flag_ext;
  String string;

  /* Initializations */

  ndim = neigh->getNDim();
  flag_ext = is_flag_data_disc_defined();

  /* Neighborhood data */

  mestitle(1,"Data selected in neighborhood");
  tab_prints(NULL,1,EJustify::RIGHT,"Rank");
  tab_prints(NULL,1,EJustify::RIGHT,"Sample");
  if (dbin->hasCode()) tab_prints(NULL,1,EJustify::RIGHT,"Code");
  for (idim=0; idim<ndim; idim++)
  {
    string = getLocatorName(ELoc::X,idim);
    tab_prints(NULL,1,EJustify::RIGHT,string.c_str());
  }
  if (flag_ext)
    for (idim=0; idim<ndim; idim++)
    {
      string = getLocatorName(ELoc::BLEX,idim);
      tab_prints(NULL,1,EJustify::RIGHT,string.c_str());
    }
  if (neigh->getType() == ENeigh::MOVING)
    tab_prints(NULL,1,EJustify::RIGHT,"Sector");
  message("\n");

  /* Loop on the sample points */

  for (iech=sel=0; iech<dbin->getSampleNumber(); iech++)
  {
    if (rank[iech] < 0) continue;
    
    tab_printi(NULL,1,EJustify::RIGHT,sel+1);
    tab_printi(NULL,1,EJustify::RIGHT,iech+1);
    if (dbin->hasCode())
      tab_printi(NULL,1,EJustify::RIGHT,static_cast<int> (dbin->getCode(iech)));
    for (idim=0; idim<ndim; idim++)
      tab_printg(NULL,1,EJustify::RIGHT,dbin->getCoordinate(iech,idim));
    if (flag_ext)
    {
      for (idim=0; idim<ndim; idim++)
        tab_printg(NULL,1,EJustify::RIGHT,dbin->getBlockExtension(iech,idim));
    }
    if (neigh->getType() == ENeigh::MOVING)
      tab_printi(NULL,1,EJustify::RIGHT,rank[iech]+1);
    message("\n"); 
    sel++;
  }

  return;
}

/****************************************************************************/
/*!
**  Discard a sample for which all variables are undefined
**
**  Returns 1 if all variables are undefined; 0 otherwise
**
** \param[in]  dbin      Input Db structure
** \param[in]  iech      Rank of the sample
**
** \remarks When the Neighborhood is performed in the case of Simulations
** \remarks checking for all variables being undefined is performed
** \remarks on ELoc::SIMU rather than on ELoc::Z
**
*****************************************************************************/
static int st_discard_undefined(Db  *dbin,
                                int  iech)
{
  int ivar;

  /* Dispatch */

  if (FLAG_NO_VAR_CHECK)
  {
    if (dbin->getVariableNumber() <= 0)
      return 0;
    else
    {
      for (ivar=0; ivar<dbin->getVariableNumber(); ivar++)
        if (! FFFF(dbin->getVariable(iech,ivar))) return 0;
      return 1;
    }
  }

  if (! FLAG_SIMU)
  {
    for (ivar=0; ivar<dbin->getVariableNumber(); ivar++)
      if (! FFFF(dbin->getVariable(iech,ivar))) return 0;
  }
  else
  {
    // In the case of simulations, the test is performed on the 
    // simulation error for the first variable and first simulation
    if (! FFFF(dbin->getSimvar(ELoc::SIMU,iech,0,0,0,1,0))) return 0;
  }

  return 1;
}

/****************************************************************************/
/*!
**  Mask the data sample in the case of cross-validation
**
** \return  1 if the sample is masked; 0 otherwise
**
** \param[in]  dbin     intput Db structure
** \param[in]  dbout    output Db structure
** \param[in]  iech_in  Rank in the input Db structure
** \param[in]  iech_out Rank in the output Db structure
** \param[in]  neigh    Neigh structure
**
*****************************************************************************/
static int st_xvalid(Db    *dbin,
                     Db    *dbout,
                     int    iech_in,
                     int    iech_out,
                     Neigh *neigh)
{
  static double eps = 1.e-8;

  if (neigh->getFlagXvalid() == 0)
    return(0);
  else if (neigh->getFlagXvalid() > 0)
  {
    if (distance_inter(dbin,dbout,iech_in,iech_out,NULL) < eps) return(1);
  }
  else
  {
    if (! dbin->hasCode()) return(0);
    if (dbin->getCode(iech_in) == dbout->getCode(iech_out)) return(1);
  }
  return(0);
}
                            
/****************************************************************************/
/*!
**  Select the unique neighborhood (or Image Neighborhood)
**  Pay attention to the cross-validation flag
**
** \param[in]  dbin      input Db structure
** \param[in]  dbout     output Db structure
** \param[in]  iech_out  rank of the output sample
** \param[in]  neigh     Neigh structure
**
** \param[out]  rank     Neighborhood selection array
**
*****************************************************************************/
static void st_unique(Db     *dbin,
                      Db     *dbout,
                      int     iech_out,
                      Neigh  *neigh,
                      int    *rank)
{
  int iech;

  /* Loop on samples */

  for (iech=0; iech<dbin->getSampleNumber(); iech++)
  {
    /* Discard the masked input sample */

    if (! dbin->isActive(iech)) continue;

    /* Discard samples where all variables are undefined */

    if (st_discard_undefined(dbin,iech)) continue;

    /* Discard the target sample for the cross-validation option */

    if (neigh->getFlagXvalid() != 0)
    {
      if (st_xvalid(dbin,dbout,iech,iech_out,neigh)) continue;
    }
    rank[iech] = 0;
  }

  return;
}

/****************************************************************************/
/*!
**  Search for the bench neighborhood, according to the last
**  coordinate
**
** \param[in]  dbin      input Db structure
** \param[in]  dbout     output Db structure
** \param[in]  iech_out  rank of the output sample
** \param[in]  neigh     Neigh structure
**
** \param[out]  rank     Neighborhood selection array
**
*****************************************************************************/
static void st_bench(Db     *dbin,
                     Db     *dbout,
                     int     iech_out,
                     Neigh  *neigh,
                     int    *rank)
{
  double z0;
  int iech,idim_bench;

  /* Initializations */

  idim_bench = dbin->getNDim() - 1;
  z0 = dbout->getCoordinate(iech_out,idim_bench);

  /* Loop on samples */

  for (iech=0; iech<dbin->getSampleNumber(); iech++)
  {
    /* Discard the masked input sample */

    if (! dbin->isActive(iech)) continue;

    /* Discard samples where all variables are undefined */

    if (st_discard_undefined(dbin,iech)) continue;

    /* Discard the target sample for the cross-validation option */

    if (neigh->getFlagXvalid() != 0)
    {
      if (st_xvalid(dbin,dbout,iech,iech_out,neigh)) continue;
    }

    /* Discard sample located outside the bench */

    if (ABS(dbin->getCoordinate(iech,idim_bench) - z0) <= neigh->getWidth())
      rank[iech] = 0;
  }

  return;
}

/****************************************************************************/
/*!
**  Returns the additional variance for continuous moving neighborhood
**
** \return  Additional variance or 0
**
** \param[in]  neigh    Neigh structure
** \param[in]  db1      First Db
** \param[in]  rank1    Rank of the sample in the first Db
** \param[in]  db2      Second Db
** \param[in]  rank2    Rank of the sample in the second Db
**
** \remarks In the case of a neighborhood which is not MOVING (or undefined),
** \remarks this function systematically returns 0.
**
*****************************************************************************/
GSTLEARN_EXPORT double neigh_continuous_variance(Neigh *neigh,
                                             Db    *db1,
                                             int    rank1,
                                             Db    *db2,
                                             int    rank2)
{
  double *dd,dist,var;
  int     ndim;
  static double eps_dist = 1.e-4;
    
  /* Initializations */

  dd = nullptr;
  if (neigh == nullptr) return(0.);
  if (neigh->getType() != ENeigh::MOVING) return(0.);
  ndim = neigh->getNDim();

  /* Core allocation */
  
  dd = (double *) mem_alloc(sizeof(double) * ndim,1);

  /* Calculate the distance increment */

  for (int idim=0; idim<ndim; idim++)
    dd[idim] = db1->getCoordinate(rank1,idim) - db2->getCoordinate(rank2,idim);

  /* Anisotropic neighborhood */

  if (neigh->getFlagAniso())
  {
      
    /* Rotated anisotropy ellipsoid */

    if (neigh->getFlagRotation())
      matrix_product_safe(1,ndim,ndim,dd,neigh->getAnisoRotMat().data(),dd);
    for (int idim=0; idim<ndim; idim++) dd[idim] /= neigh->getAnisoCoeff(idim);
  }

  /* Calculate the distance */

  matrix_product(1,ndim,1,dd,dd,&dist);
  dist = sqrt(dist) / neigh->getRadius();
  var  = 0.;
  if (dist > neigh->getDistCont())
  {
    if (ABS(1. - dist) < eps_dist) dist = 1. - eps_dist;
    var  = (dist - neigh->getDistCont()) / (1. - dist);
    var  = var * var;
  }

  /* Core deallocation */
  
  dd = (double *) mem_free((char *) dd);

  return(var);
}

/****************************************************************************/
/*!
**  Calculates the distance between an input sample and the target
**
** \return  Distance
**
** \param[in]  dbin     input Db structure
** \param[in]  dbout    output Db structure
** \param[in]  iech_in  Rank of the sample in the input Db structure
** \param[in]  iech_out Rank of the sample in the output Db structure
** \param[in]  neigh    Neigh structure
**
*****************************************************************************/
static double st_moving_dist(Db     *dbin,
                             Db     *dbout,
                             int     iech_in,
                             int     iech_out,
                             Neigh  *neigh)
{
  double dist;
  int    idim,ndim;
    
  /* Initializations */

  ndim = dbin->getNDim();

  /* Calculate the distance to the target */

  for (idim=0; idim<ndim; idim++)
    NBGH_x1[idim] = 
      dbout->getCoordinate(iech_out,idim) - dbin->getCoordinate(iech_in,idim);

  /* Anisotropic neighborhood */

  if (neigh->getFlagAniso())
  {
      
    /* Rotated anisotropy ellipsoid */

    if (neigh->getFlagRotation())
    {
      matrix_product(1,ndim,ndim,NBGH_x1,neigh->getAnisoRotMat().data(),NBGH_x2);
      (void) memcpy((void *) NBGH_x1,(void *) NBGH_x2,
                    sizeof(double) * ndim);
    }
    for (idim=0; idim<ndim; idim++) NBGH_x1[idim] /= neigh->getAnisoCoeff(idim);
  }

  /* Calculate the distance */

  matrix_product(1,ndim,1,NBGH_x1,NBGH_x1,&dist);
  dist = sqrt(dist);
  return(dist);
}

/****************************************************************************/
/*!
**  Returns the rank of the sector (using the first two coordinates)
**
** \return  Sector index
**
** \param[in]  neigh Neigh structure
** \param[in]  dx    increment along X
** \param[in]  dy    increment along Y
**
*****************************************************************************/
static int st_moving_sector_define(Neigh  *neigh,
                                   double  dx,
                                   double  dy)
{
  double angle;
  int    isect;

  isect = 0;
  if (neigh->getNSect() > 1)
  {
    if (dx == 0.)
    {
      if (dy >= 0.)
        angle = GV_PI/2.;
      else
        angle = 1.5*GV_PI;
    }
    else if (dx > 0.)
    {
      if (dy >= 0.)
        angle = atan(dy/dx);
      else
        angle = 2.*GV_PI - atan(-dy/dx);
    }
    else
    {
      if (dy > 0.)
        angle = GV_PI/2. + atan(-dx/dy);
      else
        angle = GV_PI    + atan( dy/dx);
    }
    isect = (int)(neigh->getNSect() * angle / (2.*GV_PI));
  }
  return(isect);
}

/****************************************************************************/
/*!
**  For each angular sector, select the first samples until
**  the maximum number of samples is reached
**
** \param[in]  neigh  Neigh structure
** \param[in]  nsel   Number of selected samples
**
** \param[out]  rank  Array of active data point ranks
**
*****************************************************************************/
static void st_moving_sector_nsmax(Neigh *neigh,
                                   int    nsel,
                                   int   *rank)
{
  int i,j,isect,n_end,n_ang;

  for (isect=n_end=0; isect<neigh->getNSect(); isect++)
  {
    for (i=n_ang=0; i<nsel; i++)
    {
      j = NBGH_ind[i];
      if (rank[j] != isect) continue;
      if (n_ang < neigh->getNSMax())
        n_ang++;
      else
        rank[j] = -1;
    }
    n_end += n_ang;
  }
  return;
}

/****************************************************************************/
/*!
**  Select the closest data per sector, until the maximum number
**  neighbors is reached
**
** \param[in]  neigh Neigh structure
** \param[in]  nsel  Number of ellegible data points
** \param[in]  rank  Rank of the ellegible samples
**
** \remark  The samples beyond the maximum number of neighbors have their
** \remark  rank turned into -1
**
*****************************************************************************/
static void st_moving_select(Neigh *neigh,
                             int    nsel,
                             int   *rank)
{
  int i,j,isect,jsect,number;

  /* Initializations */

  if (neigh->getNMaxi() <= 0) return;
  for (isect=0; isect<neigh->getNSect(); isect++)
    NBGH_nsect[isect] = NBGH_isect[isect] = 0;

  /* Count the number of samples per sector */

  for (i=number=0; i<nsel; i++)
  {
    j = NBGH_ind[i];
    isect = rank[j];
    if (isect < 0) continue;
    NBGH_nsect[isect]++;
    number++;
  }
  if (number < neigh->getNMaxi()) return;

  /* Find the rank of the admissible data per sector */

  number = 0;
  while (number < neigh->getNMaxi())
  {
    for (isect=0; isect<neigh->getNSect(); isect++)
    {
      if (NBGH_isect[isect] >= NBGH_nsect[isect]) continue;
      NBGH_isect[isect]++;
      number++;
      if (number >= neigh->getNMaxi()) break;
    }
  }

  /* Discard the data beyond the admissible rank per sector */

  for (isect=0; isect<neigh->getNSect(); isect++)
  {
    if (NBGH_isect[isect] >= NBGH_nsect[isect]) continue;
    for (i=number=0; i<nsel; i++)
    {
      j = NBGH_ind[i];
      jsect = rank[j];
      if (jsect < 0) continue;
      if (isect != jsect) continue;
      number++;
      if (number > NBGH_isect[isect]) rank[j] = -1;
    }
  }

  return;
}

/****************************************************************************/
/*!
**  Moving neighborhood search
**
** \return  Error return code
** \return  0 : No error
** \return  1 : The number of data is smaller than the minimum number of
** \return      data in Moving Neighborhood
**
** \param[in]  dbin      Input Db structure
** \param[in]  dbout     Output Db structure
** \param[in]  iech_out  Rank of the target in the output Db structure
** \param[in]  neigh     Neigh structure
**
** \param[out]  rank     Array of active data point ranks
**
*****************************************************************************/
static int st_moving(Db     *dbin,
                     Db     *dbout,
                     int     iech_out,
                     Neigh  *neigh,
                     int    *rank)
{
  int    iech,isect,nsel,isel;
  double dist,distmax;
  static double eps = 1.e-9;

  /* Initializations */

  isect = 0;
  if (dbin->getSampleNumber() < neigh->getNMini()) return(1);

  /* Loop on the data points */

  distmax = 0.;
  for (iech=nsel=0; iech<dbin->getSampleNumber(); iech++)
  {
    
    /* Discard the masked input sample */

    if (! dbin->isActive(iech)) continue;

    /* Discard samples where all variables are undefined */

    if (st_discard_undefined(dbin,iech)) continue;

    /* Discard the target sample for the cross-validation option */

    if (neigh->getFlagXvalid() != 0)
    {
      if (st_xvalid(dbin,dbout,iech,iech_out,neigh)) continue;
    }

    /* Calculate the distance between data and target */

    dist = st_moving_dist(dbin,dbout,iech,iech_out,neigh);
    if (! FFFF(neigh->getRadius()) && dist > neigh->getRadius()) continue;
    if (dist > distmax) distmax = dist;

    /* Calculate the angular sector to which the sample belongs */

    if (neigh->getFlagSector())
      isect = st_moving_sector_define(neigh,NBGH_x1[0],NBGH_x1[1]);

    /* The sample may be selected */

    NBGH_ind[nsel]  = iech;
    NBGH_dst[nsel]  = dist;
    rank[iech]      = isect;
    nsel++;
  }
  if (nsel < neigh->getNMini()) return(1);

  /* Slightly modify the distances in order to ensure the sorting results */
  /* In the case of equal distances                                       */

  for (isel=0; isel<nsel; isel++)
    NBGH_dst[isel] += distmax * isel * eps;

  /* Sort the selected samples according to the distance */

  ut_sort_double(0,nsel,NBGH_ind,NBGH_dst);
  
  /* For each angular sector, select the first sample up to the maximum */

  if (neigh->getFlagSector() && neigh->getNSMax() > 0)
  {
    st_moving_sector_nsmax(neigh,nsel,rank);
    if (nsel < neigh->getNMini()) return(1);
  }

  /* Select the first data samples */

  st_moving_select(neigh,nsel,rank);

  return(0);
}

/****************************************************************************/
/*!
**  Frees the Neigh structure
**
** \return  Pointer on the freed Neigh structure
**
** \param[in]  neigh Neigh structure to be freed
**
*****************************************************************************/
GSTLEARN_EXPORT Neigh *neigh_free(Neigh *neigh)

{
  delete neigh;
  neigh = NULL;
  return(neigh);
}

/****************************************************************************/
/*!
**  Initialize a Bench Neighborhood
**
** \return  Pointer to the newly created Neigh structure
**
** \param[in]  ndim         Space dimension
** \param[in]  flag_xvalid  !=0 if the cross-validation option is required
** \param[in]  width        width of the bench 
**                          (centered on the target vertically)
**
*****************************************************************************/
GSTLEARN_EXPORT Neigh *neigh_init_bench(int    ndim,
                                    int    flag_xvalid,
                                    double width)
{
  Neigh *neigh;

  neigh = neigh_init(ndim,ENeigh::BENCH,flag_xvalid,0,0,0,0,0,0,0,0,0,
                     width,0.,0.,VectorDouble(),VectorDouble(),VectorInt());

  return(neigh);
}

/****************************************************************************/
/*!
**  Initialize an Image Neighborhood
**
** \return  Pointer to the newly created Neigh structure
**
** \param[in]  ndim         Space dimension
** \param[in]  flag_xvalid  !=0 if the cross-validation option is required
** \param[in]  skip         skipping factor
** \param[in]  nbgh_image   Vector of image neighborhood radius
**
*****************************************************************************/
GSTLEARN_EXPORT Neigh *neigh_init_image(int     ndim,
                                    int     flag_xvalid,
                                    int     skip,
                                    const VectorInt& nbgh_image)
{
  Neigh *neigh;

  neigh = neigh_init(ndim,ENeigh::IMAGE,flag_xvalid,0,0,0,0,0,0,0,0,skip,
                     0.,0.,0.,VectorDouble(),VectorDouble(),nbgh_image);

  return(neigh);
}

/****************************************************************************/
/*!
**  Initialize a Unique Neighborhood
**
** \return  Pointer to the newly created Neigh structure
**
** \param[in]  ndim        Space dimension
**
** \remark  The Unique Neighborhood is incompatible with the
** \remark  Cross-Validation option
**
*****************************************************************************/
GSTLEARN_EXPORT Neigh *neigh_init_unique(int ndim)

{
  Neigh *neigh;

  neigh = neigh_init(ndim,ENeigh::UNIQUE,0,0,0,0,0,0,0,0,0,0,
                     0.,0.,0.,VectorDouble(),VectorDouble(),VectorInt());

  return(neigh);
}

/****************************************************************************/
/*!
**  Initialize a Neighborhood
**
** \return  Pointer to the newly created Neigh structure
**
*****************************************************************************/
static Neigh *st_neigh_alloc(void)

{
  Neigh *neigh;

  /* Core allocation */

  neigh = new Neigh();
/*
  neigh->setNDim(0);
  neigh->setType(0);
  neigh->setNMini(0); // 1 in default constructor
  neigh->setNMaxi(0);
  neigh->setNSect(1);
  neigh->setNSMax(0);
  neigh->setWidth(0);
  neigh->setRadius(0);
  neigh->setDistCont(0);
  neigh->setSkip(0);
  neigh->setFlagXvalid(0);
  neigh->setFlagSector(0);
  neigh->setFlagAniso(0);
  neigh->setFlagRotation(0);
  neigh->setFlagContinuous(0);
  neigh->setAnisoCoeff(VectorDouble());
  neigh->setAnisoRotMat(VectorDouble());
  neigh->setImageRadius(VectorInt());
*/
  return(neigh);
}

/****************************************************************************/
/*!
**  Load anisotropy argument from the Neigh structure
**
** \param[in]  neigh       Neigh structure
**
** \param[out] nbgh_radius Array of anisotropic search radii
**
*****************************************************************************/
static void st_get_neigh_anisotropy(Neigh *neigh,
                                    VectorDouble& nbgh_radius)
{
  int ndim = neigh->getNDim();
  nbgh_radius.resize(ndim);

  if (! neigh->getAnisoCoeff().empty())
    for (int i=0; i<ndim; i++)
      nbgh_radius[i] = neigh->getRadius() * neigh->getAnisoCoeff(i);
  else
    for (int i=0; i<ndim; i++)
      nbgh_radius[i] = neigh->getRadius();
}

/****************************************************************************/
/*!
**  Allocate the Neigh structure
**
** \return  Pointer on the Neigh structure allocated
**
** \param[in]  ndim        Space dimension
** \param[in]  type        Neighborhood type (\sa ENeigh)
** \param[in]  flag_xvalid Option
** \li                     0 no option
** \li                     >0 discard the target from the neighborhood
** \li                     <0 discard all data with same code as target
** \param[in]  flag_sector    1 if the MOVING neighborhood uses sectors
** \param[in]  flag_aniso     1 if the MOVING neighborhood is anisotropic
** \param[in]  flag_rotation  1 if the anisotropic MOVING neighborhood
**                            is rotated
** \param[in]  flag_continuous 1 for continuous mving neighborhood
** \param[in]  nmini       Minimum number of points in the neighborhood (or 0)
**                         (only used for ENeigh::MOVING)
** \param[in]  nmaxi       Maximum number of points in the neighborhood (or 0)
**                         (only used for ENeigh::MOVING)
** \param[in]  nsect       Number of angular sectors (minimum=1)
**                         (only used for ENeigh::MOVING with flag_sector)
** \param[in]  nsmax       Maximum number of samples per angular sector (or 0)
**                         (only used for ENeigh::MOVING with flag_sector)
** \param[in]  skip        Skipping factor 
**                         (only used for ENeigh::IMAGE)
** \param[in]  width       Width of the bench (centered on target)
**                         (used for ENeigh::BENCH)
** \param[in]  radius      Search radius or TEST
**                         (only used for ENeigh::MOVING)
** \param[in]  dist_cont   Normalized Distance where to apply continuous
**                         moving neighborhood
** \param[in]  nbgh_radius Array of anisotropic search radii
**                         (only used for ENeigh::MOVING)
** \param[in]  nbgh_rotmat Rotation matrix for the anisotropy
**                         (only used for ENeigh::MOVING)
** \param[in]  nbgh_image  Vector of image neighborhood radius
**                         (only used for ENeigh::IMAGE)
**
*****************************************************************************/
GSTLEARN_EXPORT Neigh *neigh_init(int ndim,
                              ENeigh type,
                              int flag_xvalid,
                              int flag_sector,
                              int flag_aniso,
                              int flag_rotation,
                              int flag_continuous,
                              int nmini,
                              int nmaxi,
                              int nsect,
                              int nsmax,
                              int skip,
                              double  width,
                              double  radius,
                              double  dist_cont,
                              const VectorDouble& nbgh_radius,
                              const VectorDouble& nbgh_rotmat,
                              const VectorInt& nbgh_image)
{
  Neigh *neigh;

  /* Initializations */

  neigh = nullptr;
  
  /* Allocation */

  neigh = st_neigh_alloc();
  if (neigh == nullptr) return(neigh);

  /* Load the parameters */

  neigh->setNDim(ndim);
  neigh->setType(type);
  neigh->setNMini(nmini);
  neigh->setNMaxi(nmaxi);
  neigh->setNSect((flag_sector) ? MAX(nsect,1) : 1);
  neigh->setNSMax(nsmax);
  neigh->setWidth(width);
  neigh->setRadius(radius);
  neigh->setDistCont(dist_cont);
  neigh->setSkip(skip);
  neigh->setFlagXvalid(flag_xvalid);
  neigh->setFlagSector(flag_sector && ndim >= 2);
  neigh->setFlagAniso(flag_aniso && ! nbgh_radius.empty());
  neigh->setFlagRotation(flag_rotation && flag_aniso && ! nbgh_rotmat.empty());
  neigh->setFlagContinuous((! IFFFF(flag_continuous)) ? flag_continuous : 0);

  /* Core allocation */

  if (neigh->getFlagAniso() && ! nbgh_radius.empty())
  {
    neigh->setRadius(0.);
    for (int i=0; i<neigh->getNDim(); i++)
      neigh->setRadius(MAX(neigh->getRadius(),nbgh_radius[i]));
    neigh->setAnisoCoeff(nbgh_radius);
    neigh->anisoRescale();
  }
  if (neigh->getFlagRotation() && ! nbgh_rotmat.empty())
    neigh->setAnisoRotMat(nbgh_rotmat);
  if (type == ENeigh::IMAGE && ! nbgh_image.empty())
    neigh->setImageRadius(nbgh_image);

  return(neigh);
}

/****************************************************************************/
/*!
**  Print the characteristics of the Neigh stucture
**
** \param[in]  neigh Neigh structure
**
*****************************************************************************/
GSTLEARN_EXPORT void neigh_print(const Neigh *neigh)

{
  int ndim,idim;
  double *ranges;

  /* Initializations */

  if (neigh == nullptr) return;
  ndim   = neigh->getNDim();
  ranges = nullptr;

  /* Neighborhood options */

  mestitle(0,"Neighborhood characteristics");

  switch (neigh->getType().toEnum())
  {
    case ENeigh::E_UNIQUE:
      message("Unique neighborhood option\n");
      message("Space dimension = %d\n",neigh->getNDim());
      break;

    case ENeigh::E_BENCH:
      message("Bench neighborhood option\n");
      message("Space dimension = %d\n",neigh->getNDim());
      message("Bench width     = %lf\n",neigh->getWidth());
      break;

    case ENeigh::E_MOVING:
      message("Moving neighborhood option\n");
      message("Space dimension                     = %d\n",neigh->getNDim());
      if (neigh->getNMini() > 0)
        message("Minimum number of samples           = %d\n",neigh->getNMini());
      if (neigh->getNMaxi() > 0)
        message("Maximum number of samples           = %d\n",neigh->getNMaxi());
      if (neigh->getFlagSector())
      {
        message("Number of angular sectors           = %d\n",neigh->getNSect());
        if (neigh->getNSMax() > 0)
          message("Maximum number of points per sector = %d\n",neigh->getNSMax());
      }
      if (neigh->getFlagContinuous())
      {
        message("Norm. dist. for continuous Neigh.   = %lf\n",neigh->getDistCont());
      }

      if (! neigh->getFlagAniso())
      {
        if (! FFFF(neigh->getRadius()))
          message("Maximum horizontal distance         = %lf\n",neigh->getRadius());
      }
      else
      {
        ranges = (double *) mem_alloc(sizeof(double) * ndim,1);
        for (idim=0; idim<ndim; idim++)
          ranges[idim] = neigh->getRadius() * neigh->getAnisoCoeff(idim);
        print_matrix("Anisotropic Ranges :",0,1,ndim,1,nullptr,ranges);
        ranges = (double *) mem_free((char *) ranges);

        if (neigh->getFlagRotation())
        {
          print_matrix("Anisotropy Rotation :",0,1,ndim,ndim,
                       nullptr,neigh->getAnisoRotMat().data());
        }
      }
      break;

    case ENeigh::E_IMAGE:
      message("Image neighborhood option\n");
      message("Skipping factor = %d\n",neigh->getSkip());
      print_imatrix("Image radius :",0,1,ndim,1,nullptr,
                   neigh->getImageRadius().data());
      break;
  }
  
  /* Cross-validation option */

  if (neigh->getFlagXvalid() != 0)
    message("The Cross-Validation Option is switched ON\n");
  
  return;
}

/****************************************************************************/
/*!
**  Calculate the Neighborhood Parameters
**
** \param[in]  dbin      input Db structure
** \param[in]  neigh     Neigh structure
** \param[in]  rank      Neighborhood selection array
** \param[in]  nsel      Number of selected samples
**
** \param[out] tab       Resulting array (Dimension = 5)
** \li                    0 : Number of active samples
** \li                    1 : Minimum distance
** \li                    2 : Maximum distance
** \li                    3 : Number of non-empty sectors
** \li                    4 : Number of consecutive empty sectors
**
*****************************************************************************/
GSTLEARN_EXPORT void neigh_echo(Db     *dbin,
                            Neigh  *neigh,
                            int    *rank,
                            int     nsel,
                            double *tab)
{
  int iech,jech,isect,n_empty,number;
  double dist,dmin,dmax;

  /* Number of selected samples */

  number = 0;
  for (iech=0; iech<dbin->getSampleNumber(); iech++)
  {
    if (rank[iech] < 0) continue;
    number++;
  }
  tab[0] = (double) number;

  /* Maximum distance */

  dmax = TEST;
  for (jech=0; jech<nsel; jech++)
  {
    dist = NBGH_dst[jech];
    if (FFFF(dmax) || dist > dmax) dmax = dist;
  }
  tab[1] = dmax;
  
  /* Minimum distance */

  dmin = TEST;
  for (jech=0; jech<nsel; jech++)
  {
    dist = NBGH_dst[jech];
    if (FFFF(dmin) || dist < dmin) dmin = dist;
  }
  tab[2] = dmin;

  /* Number of sectors containing neighborhood information */

  number = 0;
  for (isect=0; isect<neigh->getNSect(); isect++)
  {
    if (NBGH_nsect[isect] > 0) number++;
  }
  tab[3] = (double) number;
  
  /* Number of consecutive empty sectors */

  number = n_empty = 0;
  for (isect=0; isect<2 * neigh->getNSect(); isect++)
  {
    if (NBGH_nsect[isect] > 0) 
      n_empty = 0;
    else
    {
      n_empty++;
      if (n_empty > number) number = n_empty;
    }
  }
  tab[4] = (double) number;

  return;
}

/****************************************************************************/
/*!
**  Select the neighborhood
**
** \return  Error return code
**
** \param[in]  dbin          input Db structure
** \param[in]  dbout         output Db structure
** \param[in]  iech_out      rank of the output sample
** \param[in]  neigh         Neigh structure
** \param[in]  flag_simu     1 if used for Simulation
** \param[in]  flag_no_var_check 1 if no check on variable definition
**
** \param[out]  nsel  Number of active points in the neighborhood
** \param[out]  rank  Array of active data point ranks
**
** \remarks When the Neighborhood is performed in the case of Simulations
** \remarks checking for all variables being undefined is performed
** \remarks on ELoc::SIMU rather than on ELoc::Z
**
*****************************************************************************/
GSTLEARN_EXPORT int neigh_select(Db     *dbin,
                             Db     *dbout,
                             int     iech_out,
                             Neigh  *neigh,
                             int     flag_simu,
                             int     flag_no_var_check,
                             int    *nsel,
                             int    *rank)
{
  int iech;

  /* Initializations */

  if (! NBGH_INITIALIZED) return(1);
  for (iech=0; iech<dbin->getSampleNumber(); iech++) rank[iech] = -1;
  FLAG_SIMU = flag_simu;
  FLAG_NO_VAR_CHECK = flag_no_var_check;

  /* Select the active data points */

  switch (neigh->getType().toEnum())
  {
    case ENeigh::E_IMAGE:
    case ENeigh::E_UNIQUE:
      st_unique(dbin,dbout,iech_out,neigh,rank);
      break;

    case ENeigh::E_BENCH:
      st_bench(dbin,dbout,iech_out,neigh,rank);
      break;

    case ENeigh::E_MOVING:
      if (st_moving(dbin,dbout,iech_out,neigh,rank)) return(1);
      break;
  }

  /* Print the Neighborhood search result */

  if (debug_query("nbgh")) st_neigh_print(dbin,neigh,rank);

  /* Count and rank the active points */

  for (iech=(*nsel)=0; iech<dbin->getSampleNumber(); iech++)
    if (rank[iech] >= 0) rank[(*nsel)++] = iech;

  /* Blank out the rest of the array */

  for (iech=(*nsel); iech<dbin->getSampleNumber(); iech++) rank[iech] = -1;

  return(0);
}

/****************************************************************************/
/*!
**  Perform the core allocations before using the neighborhood
**
** \return  Error return code
**
** \param[in]  dbin   input Db structure
** \param[in]  neigh  Neigh structure
**
*****************************************************************************/
GSTLEARN_EXPORT int neigh_start(Db    *dbin,
                            Neigh *neigh)

{
  /* Initializations */

  if (dbin == nullptr) return(1);
  if (NBGH_INITIALIZED) return(1);

  /* Core allocation */

  NBGH_ind   = (int    *) mem_alloc(sizeof(int)    * dbin->getSampleNumber(),0);
  if (NBGH_ind   == nullptr) return(1);
  NBGH_dst   = (double *) mem_alloc(sizeof(double) * dbin->getSampleNumber(),0);
  if (NBGH_dst   == nullptr) return(1);
  NBGH_isect = (int    *) mem_alloc(sizeof(int)    * neigh->getNSect(),0);
  if (NBGH_isect == nullptr) return(1);
  NBGH_nsect = (int    *) mem_alloc(sizeof(int)    * neigh->getNSect(),0);
  if (NBGH_nsect == nullptr) return(1);
  NBGH_x1    = db_sample_alloc(dbin,ELoc::X);
  if (NBGH_x1    == nullptr) return(1);
  NBGH_x2    = db_sample_alloc(dbin,ELoc::X);
  if (NBGH_x2    == nullptr) return(1);

  NBGH_INITIALIZED = 1;
  return(0);
}

/****************************************************************************/
/*!
**  Performs the core deallocation after using neighborhood
**
*****************************************************************************/
GSTLEARN_EXPORT void neigh_stop(void)

{
  /* Initialization */

  if (! NBGH_INITIALIZED) return;

  /* Core deallocation */

  NBGH_ind   = (int    *) mem_free((char *) NBGH_ind);
  NBGH_dst   = (double *) mem_free((char *) NBGH_dst);
  NBGH_isect = (int    *) mem_free((char *) NBGH_isect);
  NBGH_nsect = (int    *) mem_free((char *) NBGH_nsect);
  NBGH_x1    = (double *) mem_free((char *) NBGH_x1);
  NBGH_x2    = (double *) mem_free((char *) NBGH_x2);

  NBGH_INITIALIZED = 0;
  return;
}

/****************************************************************************/
/*!
**  Ask the characteristics of the Neigh structure
**
** \return  Error returned code
**
** \param[in]  neigh  Neigh structure
**
** \param[out]  type        Neighborhood type
** \param[out]  nmini       Minimum number of points in the neighborhood (or 0)
**                          (only used for ENeigh::MOVING)
** \param[out]  nmaxi       Maximum number of points in the neighborhood (or 0)
**                          (only used for ENeigh::MOVING)
** \param[out]  nsect       Number of angular sectors (minimum=1)
**                          (only used for ENeigh::MOVING with flag_sector)
** \param[out]  nsmax       Maximum number of samples per angular sector (or 0)
**                          (only used for ENeigh::MOVING with flag_sector)
** \param[out]  skip        Skipping factor 
**                          (only used for ENeigh::IMAGE)
** \param[out]  flag_sector    1 if the MOVING neighborhood uses sectors
** \param[out]  flag_aniso     1 if the MOVING neighborhood is anisotropic
** \param[out]  flag_rotation  1 if the anisotropic MOVING neighborhood
**                             is rotated
** \param[out]  flag_continuous 1 for continuous mving neighborhood
** \param[out]  width       Width of the bench (centered on target)
**                          (used for ENeigh::BENCH)
** \param[out]  radius      Maximum search radius or TEST
**                          (only used for ENeigh::MOVING)
** \param[out]  dist_cont   Normalized Distance where to apply continuous
**                          moving neighborhood
** \param[out]  nbgh_rotmat Rotation matrix for the anisotropy
**                          (only used for ENeigh::MOVING)
** \param[out]  nbgh_radius Anisotropy coefficients
**                          (only used for ENeigh::MOVING)
** \param[out]  nbgh_image  Vector of image neighborhood radius
**                          (only used for ENeigh::IMAGE)
**
*****************************************************************************/
GSTLEARN_EXPORT int neigh_extract(Neigh  *neigh,
                              ENeigh *type,
                              int    *nmini,
                              int    *nmaxi,
                              int    *nsect,
                              int    *nsmax,
                              int    *skip,
                              int    *flag_sector,
                              int    *flag_aniso,
                              int    *flag_rotation,
                              int    *flag_continuous,
                              double *width,
                              double *radius,
                              double *dist_cont,
                              VectorDouble& nbgh_rotmat,
                              VectorDouble& nbgh_radius,
                              VectorInt& nbgh_image)
{
  int i,j,ecr,ndim;

  /* Initializations */

  ndim = neigh->getNDim();

  /* Returning arguments */

  *type  = neigh->getType();
  *nmini = neigh->getNMini();
  *nmaxi = neigh->getNMaxi();
  *nsect = neigh->getNSect();
  *nsmax = neigh->getNSMax();
  *skip  = neigh->getSkip();
  *flag_sector     = neigh->getFlagSector();
  *flag_aniso      = neigh->getFlagAniso();
  *flag_rotation   = neigh->getFlagRotation();
  *flag_continuous = neigh->getFlagContinuous();
  *width     = neigh->getWidth();
  *radius    = neigh->getRadius();
  *dist_cont = neigh->getDistCont();
  
  if (! neigh->getAnisoRotMat().empty())
    nbgh_rotmat = neigh->getAnisoRotMat();
  else
    for (i=ecr=0; i<ndim; i++)
      for (j=0; j<ndim; j++,ecr++)
        nbgh_rotmat[ecr] = (i == j);

  st_get_neigh_anisotropy(neigh,nbgh_radius);

  if (! neigh->getImageRadius().empty())
    nbgh_image = neigh->getImageRadius();
  else
    nbgh_image.resize(ndim);

  return(0);
}
