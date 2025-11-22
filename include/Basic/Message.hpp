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

#include "Basic/String.hpp"
#include "Basic/VectorNumT.hpp"
#include "geoslib_define.h"
#include "gstlearn_export.hpp"

namespace gstlrn
{
// Set of functions regarding the printout
GSTLEARN_EXPORT void messageFlush(const String& string);
GSTLEARN_EXPORT void messerrFlush(const String& string);
GSTLEARN_EXPORT void messerr(const char* format, ...);
GSTLEARN_EXPORT void message(const char* format, ...);
GSTLEARN_EXPORT void messageNoDiff(const char* format, ...);
GSTLEARN_EXPORT void mesArg(const char* title, Id current, Id nmax);
GSTLEARN_EXPORT bool checkArg(const char* title, Id current, Id nmax);
GSTLEARN_EXPORT void messageAbort(const char* format, ...);
GSTLEARN_EXPORT void mestitle(Id level, const char* format, ...);
GSTLEARN_EXPORT void mes_process(const char* string, Id ntot, Id iech);

// Old-fashion printing formats
GSTLEARN_EXPORT void printElement(const String& title,
                                  const String& string,
                                  Id ncol          = 1,
                                  Id justification = 1);
GSTLEARN_EXPORT void printElement(const String& title,
                                  double value,
                                  Id ncol             = 1,
                                  Id justification    = 1,
                                  bool roundZero      = true,
                                  bool flagScientific = false);
GSTLEARN_EXPORT void printElement(const String& title,
                                  Id value,
                                  Id ncol          = 1,
                                  Id justification = 1);

GSTLEARN_EXPORT void printVector(const String& title,
                                 const VectorDouble& tab,
                                 bool flagIgnoreMaxNCols = false,
                                 bool newLineAfterTitle  = false);
GSTLEARN_EXPORT void printVector(const String& title,
                                 const VectorInt& tab,
                                 bool flagIgnoreMaxNCols = false,
                                 bool newLineAfterTitle  = false);

GSTLEARN_EXPORT void printMatrix(const String& title,
                                 Id flag_limit,
                                 Id bycol,
                                 Id nx,
                                 Id ny,
                                 const VectorDouble& tab);
GSTLEARN_EXPORT void printMatrix(const String& title,
                                 Id flag_limit,
                                 Id bycol,
                                 Id nx,
                                 Id ny,
                                 const VectorInt& tab);
GSTLEARN_EXPORT void printMatrix(const String& title,
                                 Id flag_limit,
                                 const AMatrix& mat);
GSTLEARN_EXPORT void printTriMat(const String& title,
                                 Id mode,
                                 Id neq,
                                 const VectorDouble& tl);

} // namespace gstlrn