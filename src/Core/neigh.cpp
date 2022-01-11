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
#include "Basic/Utilities.hpp"
#include "Basic/EJustify.hpp"
#include "Db/Db.hpp"
#include "Neigh/Neigh.hpp"
#include "Neigh/NeighWork.hpp"

#include <math.h>
#include <string.h>

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
double neigh_continuous_variance(Neigh *neigh,
                                 Db *db1,
                                 int rank1,
                                 Db *db2,
                                 int rank2)
{
  double *dd, dist, var;
  int ndim;
  static double eps_dist = 1.e-4;

  /* Initializations */

  dd = nullptr;
  if (neigh == nullptr) return (0.);
  if (neigh->getType() != ENeigh::MOVING) return (0.);
  ndim = neigh->getNDim();

  /* Core allocation */

  dd = (double*) mem_alloc(sizeof(double) * ndim, 1);

  /* Calculate the distance increment */

  for (int idim = 0; idim < ndim; idim++)
    dd[idim] = db1->getCoordinate(rank1, idim)
        - db2->getCoordinate(rank2, idim);

  /* Anisotropic neighborhood */

  if (neigh->getFlagAniso())
  {

    /* Rotated anisotropy ellipsoid */

    if (neigh->getFlagRotation())
      matrix_product_safe(1, ndim, ndim, dd, neigh->getAnisoRotMats().data(),
                          dd);
    for (int idim = 0; idim < ndim; idim++)
      dd[idim] /= neigh->getAnisoCoeff(idim);
  }

  /* Calculate the distance */

  matrix_product(1, ndim, 1, dd, dd, &dist);
  dist = sqrt(dist) / neigh->getRadius();
  var = 0.;
  if (dist > neigh->getDistCont())
  {
    if (ABS(1. - dist) < eps_dist) dist = 1. - eps_dist;
    var = (dist - neigh->getDistCont()) / (1. - dist);
    var = var * var;
  }

  /* Core deallocation */

  dd = (double*) mem_free((char* ) dd);

  return (var);
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
Neigh* neigh_free(Neigh *neigh)

{
  delete neigh;
  neigh = NULL;
  return (neigh);
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
Neigh* neigh_init_bench(int ndim, int flag_xvalid, double width)
{
  Neigh *neigh;

  neigh = neigh_init(ndim, ENeigh::BENCH, flag_xvalid, 0, 0, 0, 0, 0, 0, 0, 0,
                     0, width, 0., 0., VectorDouble(), VectorDouble(),
                     VectorInt());

  return (neigh);
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
Neigh* neigh_init_image(int ndim,
                        int flag_xvalid,
                        int skip,
                        const VectorInt &nbgh_image)
{
  Neigh *neigh;

  neigh = neigh_init(ndim, ENeigh::IMAGE, flag_xvalid, 0, 0, 0, 0, 0, 0, 0, 0,
                     skip, 0., 0., 0., VectorDouble(), VectorDouble(),
                     nbgh_image);

  return (neigh);
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
Neigh* neigh_init_unique(int ndim)

{
  Neigh *neigh;

  neigh = neigh_init(ndim, ENeigh::UNIQUE, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0., 0.,
                     0., VectorDouble(), VectorDouble(), VectorInt());

  return (neigh);
}

/****************************************************************************/
/*!
 **  Initialize a Neighborhood
 **
 ** \return  Pointer to the newly created Neigh structure
 **
 *****************************************************************************/
static Neigh* st_neigh_alloc(void)
{
  Neigh *neigh;

  /* Core allocation */

  neigh = new Neigh();
  return (neigh);
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
static void st_get_neigh_anisotropy(Neigh *neigh, VectorDouble &nbgh_radius)
{
  int ndim = neigh->getNDim();
  nbgh_radius.resize(ndim);

  if (!neigh->getAnisoCoeffs().empty())
    for (int i = 0; i < ndim; i++)
      nbgh_radius[i] = neigh->getRadius() * neigh->getAnisoCoeff(i);
  else
    for (int i = 0; i < ndim; i++)
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
Neigh* neigh_init(int ndim,
                  const ENeigh &type,
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
                  double width,
                  double radius,
                  double dist_cont,
                  const VectorDouble &nbgh_radius,
                  const VectorDouble &nbgh_rotmat,
                  const VectorInt &nbgh_image)
{
  Neigh *neigh;

  /* Initializations */

  neigh = nullptr;

  /* Allocation */

  neigh = st_neigh_alloc();
  if (neigh == nullptr) return (neigh);

  /* Load the parameters */

  neigh->setNDim(ndim);
  neigh->setType(type);
  neigh->setNMini(nmini);
  neigh->setNMaxi(nmaxi);
  neigh->setNSect((flag_sector) ? MAX(nsect, 1) : 1);
  neigh->setNSMax(nsmax);
  neigh->setWidth(width);
  neigh->setRadius(radius);
  neigh->setDistCont(dist_cont);
  neigh->setSkip(skip);
  neigh->setFlagXvalid(flag_xvalid);
  neigh->setFlagSector(flag_sector && ndim >= 2);
  neigh->setFlagAniso(flag_aniso && !nbgh_radius.empty());
  neigh->setFlagRotation(flag_rotation && flag_aniso && !nbgh_rotmat.empty());
  neigh->setFlagContinuous((!IFFFF(flag_continuous)) ? flag_continuous :
                                                       0);

  /* Core allocation */

  if (neigh->getFlagAniso() && !nbgh_radius.empty())
  {
    neigh->setRadius(0.);
    for (int i = 0; i < neigh->getNDim(); i++)
      neigh->setRadius(MAX(neigh->getRadius(), nbgh_radius[i]));
    neigh->setAnisoCoeffs(nbgh_radius);
    neigh->anisoRescale();
  }
  if (neigh->getFlagRotation() && !nbgh_rotmat.empty())
    neigh->setAnisoRotMat(nbgh_rotmat);
  if (type == ENeigh::IMAGE && !nbgh_image.empty())
    neigh->setImageRadius(nbgh_image);

  return (neigh);
}

/****************************************************************************/
/*!
 **  Print the characteristics of the Neigh stucture
 **
 ** \param[in]  neigh Neigh structure
 **
 *****************************************************************************/
void neigh_print(const Neigh *neigh)

{
  int ndim, idim;
  double *ranges;

  /* Initializations */

  if (neigh == nullptr) return;
  ndim = neigh->getNDim();
  ranges = nullptr;

  /* Neighborhood options */

  mestitle(0, "Neighborhood characteristics");

  switch (neigh->getType().toEnum())
  {
    case ENeigh::E_UNIQUE:
      message("Unique neighborhood option\n");
      message("Space dimension = %d\n", neigh->getNDim());
      break;

    case ENeigh::E_BENCH:
      message("Bench neighborhood option\n");
      message("Space dimension = %d\n", neigh->getNDim());
      message("Bench width     = %lf\n", neigh->getWidth());
      break;

    case ENeigh::E_MOVING:
      message("Moving neighborhood option\n");
      message("Space dimension                     = %d\n", neigh->getNDim());
      if (neigh->getNMini() > 0)
        message("Minimum number of samples           = %d\n",
                neigh->getNMini());
      if (neigh->getNMaxi() > 0)
        message("Maximum number of samples           = %d\n",
                neigh->getNMaxi());
      if (neigh->getFlagSector())
      {
        message("Number of angular sectors           = %d\n",
                neigh->getNSect());
        if (neigh->getNSMax() > 0)
          message("Maximum number of points per sector = %d\n",
                  neigh->getNSMax());
      }
      if (neigh->getFlagContinuous())
      {
        message("Norm. dist. for continuous Neigh.   = %lf\n",
                neigh->getDistCont());
      }

      if (!neigh->getFlagAniso())
      {
        if (!FFFF(neigh->getRadius()))
          message("Maximum horizontal distance         = %lf\n",
                  neigh->getRadius());
      }
      else
      {
        ranges = (double*) mem_alloc(sizeof(double) * ndim, 1);
        for (idim = 0; idim < ndim; idim++)
          ranges[idim] = neigh->getRadius() * neigh->getAnisoCoeff(idim);
        print_matrix("Anisotropic Ranges :", 0, 1, ndim, 1, nullptr, ranges);
        ranges = (double*) mem_free((char* ) ranges);

        if (neigh->getFlagRotation())
        {
          print_matrix("Anisotropy Rotation :", 0, 1, ndim, ndim, nullptr,
                       neigh->getAnisoRotMats().data());
        }
      }
      break;

    case ENeigh::E_IMAGE:
      message("Image neighborhood option\n");
      message("Skipping factor = %d\n", neigh->getSkip());
      print_imatrix("Image radius :", 0, 1, ndim, 1, nullptr,
                    neigh->getAllImageRadius().data());
      break;
  }

  /* Cross-validation option */

  if (neigh->getFlagXvalid() != 0)
    message("The Cross-Validation Option is switched ON\n");

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
int neigh_extract(Neigh *neigh,
                  ENeigh *type,
                  int *nmini,
                  int *nmaxi,
                  int *nsect,
                  int *nsmax,
                  int *skip,
                  int *flag_sector,
                  int *flag_aniso,
                  int *flag_rotation,
                  int *flag_continuous,
                  double *width,
                  double *radius,
                  double *dist_cont,
                  VectorDouble &nbgh_rotmat,
                  VectorDouble &nbgh_radius,
                  VectorInt &nbgh_image)
{
  int i, j, ecr, ndim;

  /* Initializations */

  ndim = neigh->getNDim();

  /* Returning arguments */

  *type = neigh->getType();
  *nmini = neigh->getNMini();
  *nmaxi = neigh->getNMaxi();
  *nsect = neigh->getNSect();
  *nsmax = neigh->getNSMax();
  *skip = neigh->getSkip();
  *flag_sector = neigh->getFlagSector();
  *flag_aniso = neigh->getFlagAniso();
  *flag_rotation = neigh->getFlagRotation();
  *flag_continuous = neigh->getFlagContinuous();
  *width = neigh->getWidth();
  *radius = neigh->getRadius();
  *dist_cont = neigh->getDistCont();

  if (!neigh->getAnisoRotMats().empty())
    nbgh_rotmat = neigh->getAnisoRotMats();
  else
    for (i = ecr = 0; i < ndim; i++)
      for (j = 0; j < ndim; j++, ecr++)
        nbgh_rotmat[ecr] = (i == j);

  st_get_neigh_anisotropy(neigh, nbgh_radius);

  if (!neigh->getAllImageRadius().empty())
    nbgh_image = neigh->getAllImageRadius();
  else
    nbgh_image.resize(ndim);

  return (0);
}
