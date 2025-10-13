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

#include "geoslib_define.h"
#include "gstlearn_export.hpp"

#include "Enum/EJustify.hpp"
#include "Basic/AStringFormat.hpp"

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
GSTLEARN_EXPORT void   messageFlush(const String& string);
GSTLEARN_EXPORT void   messerrFlush(const String& string);
GSTLEARN_EXPORT void   messerr(const char *format,...);
GSTLEARN_EXPORT void   message(const char *format,...);
GSTLEARN_EXPORT void   messageNoDiff(const char *format,...);
GSTLEARN_EXPORT void   mesArg(const char* title, Id current, Id nmax);
GSTLEARN_EXPORT bool   checkArg(const char* title, Id current, Id nmax);
GSTLEARN_EXPORT void   messageAbort(const char* format, ...);
GSTLEARN_EXPORT void   mestitle(Id level,const char *format,...);
GSTLEARN_EXPORT void   mes_process(const char *string, Id ntot, Id iech);
GSTLEARN_EXPORT String toTitle(Id level, const char* format, ...);
GSTLEARN_EXPORT String toMatrix(const String &title,
                                const AMatrix &mat,
                                bool flagOverride = false,
                                bool flagSkipZero = false);
GSTLEARN_EXPORT String toMatrix(const String& title,
                                const VectorString& colnames,
                                const VectorString& rownames,
                                bool bycol,
                                Id nrows,
                                Id ncols,
                                const VectorDouble &tab,
                                bool flagOverride = false,
                                bool flagSkipZero = false);
GSTLEARN_EXPORT String toMatrix(const String& title,
                                const VectorString& colnames,
                                const VectorString& rownames,
                                bool bycol,
                                Id nrows,
                                Id ncols,
                                const double* tab,
                                bool flagOverride = false,
                                bool flagSkipZero = false);
GSTLEARN_EXPORT String toMatrix(const String& title,
                                const VectorString& colnames,
                                const VectorString& rownames,
                                bool bycol,
                                Id nrows,
                                Id ncols,
                                const VectorInt &tab,
                                bool flagOverride = false,
                                bool flagSkipZero = false);
GSTLEARN_EXPORT String toVector(const String& title,
                                const VectorDouble& tab,
                                bool flagOverride = true);
GSTLEARN_EXPORT String toVector(const String& title,
                                const VectorVectorDouble& tab,
                                bool flagOverride = true);
GSTLEARN_EXPORT String toVector(const String& title,
                                const VectorVectorInt& tab,
                                bool flagOverride = true);
GSTLEARN_EXPORT String toVector(const String& title, const VectorInt& tab, bool flagOverride = true);
GSTLEARN_EXPORT String toVector(const String& title,
                                const VectorString& tab,
                                bool flagOverride = true);
GSTLEARN_EXPORT String toVector(const String& title, 
                                constvect tab, 
                                bool flagOverride = true);

GSTLEARN_EXPORT String toStr(const String& string,
                             const EJustify& justify = EJustify::fromKey("RIGHT"),
                             Id localSize = 0);
GSTLEARN_EXPORT String toDouble(double value,
                                const EJustify& justify = EJustify::fromKey("RIGHT"));
GSTLEARN_EXPORT String toInt(Id value,
                             const EJustify& justify = EJustify::fromKey("RIGHT"));
GSTLEARN_EXPORT String toInterval(double zmin, double zmax);
GSTLEARN_EXPORT VectorString toVectorDouble(const VectorDouble& values,
                                            const EJustify& justify = EJustify::fromKey("RIGHT"));

// Old-fashion printing formats
GSTLEARN_EXPORT void tab_prints(const char* title,
                                const char* string,
                                Id ncol = 1,
                                const EJustify &justify = EJustify::fromKey("RIGHT"));
GSTLEARN_EXPORT void tab_printg(const char *title,
                                double value,
                                Id ncol = 1,
                                const EJustify &justify = EJustify::fromKey("RIGHT"));
GSTLEARN_EXPORT void tab_printd(const char *title,
                                double value,
                                Id ncol = 1,
                                const EJustify &justify = EJustify::fromKey("RIGHT"));
GSTLEARN_EXPORT void tab_printi(const char *title,
                                Id value,
                                Id ncol = 1,
                                const EJustify &justify = EJustify::fromKey("RIGHT"));
GSTLEARN_EXPORT void tab_print_rc(const char *title,
                                  Id mode,
                                  Id value,
                                  Id ncol = 1,
                                  const EJustify &justify = EJustify::fromKey("RIGHT"));
GSTLEARN_EXPORT void tab_print_rowname(const char *string, Id taille);
GSTLEARN_EXPORT void print_matrix(const char *title,
                                  Id flag_limit,
                                  Id bycol,
                                  Id nx,
                                  Id ny,
                                  const double *sel,
                                  const double *tab);
GSTLEARN_EXPORT void print_matrix(const char *title,
                                  Id flag_limit,
                                  const AMatrix& mat);
GSTLEARN_EXPORT void print_trimat(const char *title,
                                  Id mode,
                                  Id neq,
                                  const double *tl);
GSTLEARN_EXPORT void print_imatrix(const char *title,
                                   Id flag_limit,
                                   Id bycol,
                                   Id nx,
                                   Id ny,
                                   const double *sel,
                                   const Id *tab);
GSTLEARN_EXPORT void print_vector(const char *title,
                                  Id flag_limit,
                                  Id ntab,
                                  const double *tab);
GSTLEARN_EXPORT void print_vector(const char *title,
                                  Id flag_limit,
                                  Id ntab,
                                  const VectorDouble &tab);
GSTLEARN_EXPORT void print_ivector(const char *title,
                                   Id flag_limit,
                                   Id ntab,
                                   const Id *itab);
GSTLEARN_EXPORT void print_ivector(const char *title,
                                   Id flag_limit,
                                   Id ntab,
                                   const VectorInt &itab);
}