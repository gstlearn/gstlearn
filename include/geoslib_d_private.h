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
#ifndef GEOSLIB_D_PRIVATE_H
#define GEOSLIB_D_PRIVATE_H

/* External function definition */

#if defined(_WIN32) || defined(_WIN64)
#define GEOSLIB_API __declspec(dllexport)
#if !defined(strcasecmp)
#define strcasecmp _stricmp
#endif
#if !defined(strncasecmp)
#define strncasecmp _strnicmp
#endif
#else
#define GEOSLIB_API extern
#include <dirent.h>
#endif

#endif
