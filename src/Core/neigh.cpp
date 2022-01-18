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
  Neigh* neigh = new Neigh();
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
  neigh->setFlagContinuous((!IFFFF(flag_continuous)) ? flag_continuous : 0);

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
