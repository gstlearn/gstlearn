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

#define ENUM_CST ECst, NTCAR, \
                     NTCAR,   1, "Number of characters in printout", \
                     NTDEC,   2, "Number of decimal digits in printout", \
                     NTROW,   3, "Maximum number of rows in table printout", \
                     NTCOL,   4, "Maximum number of columns in table printout", \
                     NTBATCH, 5, "Number of elements per line for display", \
                     NTNAME,  6, "Maximum number of characters for Names", \
                     NTRANK,  7, "Maximum Number of characters for Ranks", \
                     NPROC,   8, "Percentage for Display the Progress Bar", \
                     LOCMOD,  9, "Option for updating locator of new variable", \
                     LOCNEW, 10, "When defining new locator, option for old ones", \
                     TOLINV, 11, "Tolerance for matrix inversion", \
                     TOLGEN, 12, "Tolerance for matrix generalized inversion", \
                     EPSMAT, 13, "Tolerance value for Matrix calculations", \
                     EPSSVD, 14, "Tolerance value for SVD Matrix calculations", \
                     ASP,    15, "Graphic Aspect Ratio"

ENUM_DECLARE(ENUM_CST)
