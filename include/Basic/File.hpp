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
#pragma once

#include "gstlearn_export.hpp"
#include "geoslib_define.h"

#include <fstream>

/**
 * This file contains all the OLD-STYLE declarations causing warnings on Windows
 * They should gradually be replaced by modern statements
 * However, they are kept there to keep track on these statements
 */


// TODO : File manipulation class

// Skips the Byte Order Mark (BOM) that defines UTF-8 in some text files.
//https://stackoverflow.com/a/17219495
GSTLEARN_EXPORT void skipBOM(std::ifstream &in);
GSTLEARN_EXPORT FILE* gslFopen(const char *path, const char* mode);
GSTLEARN_EXPORT FILE* gslFopen(const String& path, const String& mode);
GSTLEARN_EXPORT bool gslFileExist(const char *path, const char* mode);
GSTLEARN_EXPORT bool gslFileExist(const String& path, const String& mode);

GSTLEARN_EXPORT String gslGetEnv(const String& name);

