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
  GSTLEARN_EXPORT bool checkArg(const char* title, Id current, Id nmax);
  GSTLEARN_EXPORT void messageAbort(const char* format, ...);
  GSTLEARN_EXPORT void mestitle(Id level, const char* format, ...);
  GSTLEARN_EXPORT void mes_process(const char* string, Id ntot, Id iech);
  GSTLEARN_EXPORT void
    messerrNotImplemented(const String& className, const String& methodName);

  // Old-fashion printing formats
  GSTLEARN_EXPORT void printElement(
    const String& string,
    const String& title = String(),
    Id ncol = 1,
    Id justification = 1);
  GSTLEARN_EXPORT void printElement(
    double value,
    const String& title = String(),
    Id ncol = 1,
    Id justification = 1,
    bool roundZero = true,
    bool flagScientific = false);
  GSTLEARN_EXPORT void printElement(
    Id value,
    const String& title = String(),
    Id ncol = 1,
    Id justification = 1);

  GSTLEARN_EXPORT void printVector(
    const VectorDouble& tab,
    const String& title = String(),
    bool flagIgnoreMaxNCols = false,
    bool newLineAfterTitle = false);
  GSTLEARN_EXPORT void printVector(
    const VectorInt& tab,
    const String& title = String(),
    bool flagIgnoreMaxNCols = false,
    bool newLineAfterTitle = false);

  GSTLEARN_EXPORT void printMatrix(
    const VectorDouble& tab,
    Id nx,
    Id ny,
    const String& title = String(),
    Id flag_limit = 0,
    Id bycol = 0);
  GSTLEARN_EXPORT void printMatrix(
    const VectorInt& tab,
    Id nx,
    Id ny,
    const String& title = String(),
    Id flag_limit = 0,
    Id bycol = 0);
  GSTLEARN_EXPORT void printMatrix(
    const AMatrix& mat,
    const String& title = String(),
    Id flag_limit = 0);

  GSTLEARN_EXPORT void printTriMat(
    const VectorDouble& tl,
    Id neq,
    Id mode = 1,
    const String& title = String());

} // namespace gstlrn
