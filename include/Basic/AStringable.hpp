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

#include "Basic/AStringFormat.hpp"
#include "Basic/String.hpp"
#include "Basic/VectorNumT.hpp"
#include "geoslib_define.h"
#include "gstlearn_export.hpp"

namespace gstlrn
{
class AMatrix;

class GSTLEARN_EXPORT AStringable
{
public:
  AStringable();
  AStringable(const AStringable& r);
  AStringable& operator=(const AStringable& r);
  virtual ~AStringable();

  virtual String toString(const AStringFormat* strfmt = nullptr) const;

  virtual void display(const AStringFormat* strfmt = nullptr) const final;
#ifndef SWIG // TODO : overload not available in customized SWIG 4.2.3 and more
  virtual void display(Id level) const final;
#endif
};

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
GSTLEARN_EXPORT void tab_prints(const String& title,
                                const String& string,
                                Id ncol          = 1,
                                Id justification = 1);
GSTLEARN_EXPORT void tab_printg(const String& title,
                                double value,
                                Id ncol          = 1,
                                Id justification = 1);
GSTLEARN_EXPORT void tab_printd(const String& title,
                                double value,
                                Id ncol          = 1,
                                Id justification = 1);
GSTLEARN_EXPORT void tab_printi(const String& title,
                                Id value,
                                Id ncol          = 1,
                                Id justification = 1);

GSTLEARN_EXPORT void tab_print_rc(const String& title,
                                  Id mode,
                                  Id value,
                                  Id ncol          = 1,
                                  Id justification = 1);
GSTLEARN_EXPORT void tab_print_rowname(const char* string, Id taille);
GSTLEARN_EXPORT void print_matrix(const String& title,
                                  Id flag_limit,
                                  Id bycol,
                                  Id nx,
                                  Id ny,
                                  const double* sel,
                                  const double* tab);
GSTLEARN_EXPORT void print_matrix(const String& title,
                                  Id flag_limit,
                                  const AMatrix& mat);
GSTLEARN_EXPORT void print_trimat(const String& title,
                                  Id mode,
                                  Id neq,
                                  const double* tl);
GSTLEARN_EXPORT void print_imatrix(const String& title,
                                   Id flag_limit,
                                   Id bycol,
                                   Id nx,
                                   Id ny,
                                   const double* sel,
                                   const Id* tab);
GSTLEARN_EXPORT void print_vector(const String& title,
                                  Id flag_limit,
                                  Id ntab,
                                  const double* tab);
GSTLEARN_EXPORT void print_vector(const String& title,
                                  Id flag_limit,
                                  Id ntab,
                                  const VectorDouble& tab);
GSTLEARN_EXPORT void print_ivector(const String& title,
                                   Id flag_limit,
                                   Id ntab,
                                   const Id* itab);
GSTLEARN_EXPORT void print_ivector(const String& title,
                                   Id flag_limit,
                                   Id ntab,
                                   const VectorInt& itab);
} // namespace gstlrn