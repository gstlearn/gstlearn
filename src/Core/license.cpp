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

/****************************************************************************/
/*!
 *  Check if the Library Geoslib is authorized and setup the License Id.
 *
 * \return Error return code (error = 1)
 *
 * \param[in]  target_name   Target Name
 *
 * \remarks The License File must have been defined beforehand. 
 * \remarks This may be done using register_license_file() where
 * \remarks the License File is checked
 * \remarks The name of the License File can also be passed using the
 * \remarks Environment Variable GEOSLIB_LICENSE (on LINUX system only)
 *
 ****************************************************************************/
GEOSLIB_API int setup_license(const char *target_name)
{
  if (! LicenseKey::registerLicense(target_name)) return(1);
  return(0);
}

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
GEOSLIB_API int register_license_file(const char *file_name,
                                      const char *target_name)
{
  if (! LicenseKey::registerLicenseFromFile(target_name, file_name)) return(1);
  return(0);
}

