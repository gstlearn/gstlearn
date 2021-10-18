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

#include "Basic/String.hpp"


/**
 * This file contains all the OLD-STYLE declarations causing warnings on Windows
 * They should gradually be replaced by modern statements
 * However, they are kept there to keep track on these statements
 */
template<typename T, size_t N> size_t gslArraySize(T (&array)[N]) { return N; }

FILE* gslFopen(const char *path, const char* mode);
FILE* gslFopen(const String& path, const String& mode);
bool gslFileExist(const char *path, const char* mode);
bool gslFileExist(const String& path, const String& mode);

char* gslStrcpy(char* dst, int dst_size, const char* src);
void gslStrcpy(String& dst, const String& src);

char* gslStrcat(char* dst, int dst_size, const char* src);
void gslStrcat(String& dst, const String& src);

int gslSPrintf(char* dst, int dst_size, const char* format, ...);
int gslSPrintf(String& dst, String format, ...);
