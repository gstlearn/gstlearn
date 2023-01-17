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

#include "Enum/AEnum.hpp"

/**
 * Define a set of constant values.
 * Their initial value is provided in the file OptCst.cpp
 * Note: 'LOCMOD' is still in use in RGeostats (not to be deleted)
 */
#define ENUM_CST ECst, NTCAR, \
                     NTCAR,   1, "Number of characters in printout", \
                     NTDEC,   2, "Number of decimal digits in printout", \
                     NTROW,   3, "Maximum number of rows in table printout", \
                     NTCOL,   4, "Maximum number of columns in table printout", \
                     NTBATCH, 5, "Number of elements per line for display", \
                     NTNAME,  6, "Maximum number of characters for Names", \
                     NTRANK,  7, "Maximum Number of characters for Ranks", \
                     NPROC,   8, "Percentage for Display the Progress Bar", \
                     LOCMOD,  9, "Default for setting Locator"

ENUM_DECLARE(ENUM_CST)
