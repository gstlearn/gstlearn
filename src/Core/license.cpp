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
#include "geoslib_old_f.h"
#include "License/LicenseKey.hpp"

/****************************************************************************/
/*!
 *  Register the license from the License File
 *
 * \return Error return code
 *
 * \param[in]  file_name     File Name
 * \param[in]  target_name   Target Name
 *
 ****************************************************************************/
int register_license_file(const char *file_name,
                                          const char *target_name)
{
  if (! LicenseKey::registerLicenseFromFile(target_name, file_name)) return(1);
  return(0);
}

