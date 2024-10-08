/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES Paris / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://gstlearn.org                                              */
/* License: BSD 3-clause                                                      */
/*                                                                            */
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
                     NPROC,   8, "Number of Quantiles for Display the Progress Bar", \
                     LOCMOD,  9, "Default for setting Locator"

ENUM_DECLARE(ENUM_CST)
