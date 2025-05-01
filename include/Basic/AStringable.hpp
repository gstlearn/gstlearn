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
#ifndef SWIG // TODO : overload not available in SWIG 4.2.0b
  virtual void display(int level) const final;
#endif
};

// Set of functions regarding the printout
GSTLEARN_EXPORT void   messageFlush(const String& string);
GSTLEARN_EXPORT void   messerrFlush(const String& string);
GSTLEARN_EXPORT void   messerr(const char *format,...);
GSTLEARN_EXPORT void   message(const char *format,...);
GSTLEARN_EXPORT void   messageNoDiff(const char *format,...);
GSTLEARN_EXPORT void   mesArg(const char* title, int current, int nmax);
GSTLEARN_EXPORT bool   checkArg(const char* title, int current, int nmax);
GSTLEARN_EXPORT void   messageAbort(const char* format, ...);
GSTLEARN_EXPORT void   mestitle(int level,const char *format,...);
GSTLEARN_EXPORT void   mes_process(const char *string, int ntot, int iech);
GSTLEARN_EXPORT String toTitle(int level, const char* format, ...);
GSTLEARN_EXPORT String toMatrix(const String &title,
                                const AMatrix &mat,
                                bool flagOverride = false,
                                bool flagSkipZero = false);
GSTLEARN_EXPORT String toMatrix(const String& title,
                                const VectorString& colnames,
                                const VectorString& rownames,
                                bool bycol,
                                int nrows,
                                int ncols,
                                const VectorDouble &tab,
                                bool flagOverride = false,
                                bool flagSkipZero = false);
GSTLEARN_EXPORT String toMatrix(const String& title,
                                const VectorString& colnames,
                                const VectorString& rownames,
                                bool bycol,
                                int nrows,
                                int ncols,
                                const double* tab,
                                bool flagOverride = false,
                                bool flagSkipZero = false);
GSTLEARN_EXPORT String toMatrix(const String& title,
                                const VectorString& colnames,
                                const VectorString& rownames,
                                bool bycol,
                                int nrows,
                                int ncols,
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
                             int localSize = 0);
GSTLEARN_EXPORT String toDouble(double value,
                                const EJustify& justify = EJustify::fromKey("RIGHT"));
GSTLEARN_EXPORT String toInt(int value,
                             const EJustify& justify = EJustify::fromKey("RIGHT"));
GSTLEARN_EXPORT String toInterval(double zmin, double zmax);
GSTLEARN_EXPORT VectorString toVectorDouble(const VectorDouble& values,
                                            const EJustify& justify = EJustify::fromKey("RIGHT"));

// Old-fashion printing formats
GSTLEARN_EXPORT void tab_prints(const char* title,
                                const char* string,
                                int ncol = 1,
                                const EJustify &justify = EJustify::fromKey("RIGHT"));
GSTLEARN_EXPORT void tab_printg(const char *title,
                                double value,
                                int ncol = 1,
                                const EJustify &justify = EJustify::fromKey("RIGHT"));
GSTLEARN_EXPORT void tab_printd(const char *title,
                                double value,
                                int ncol = 1,
                                const EJustify &justify = EJustify::fromKey("RIGHT"));
GSTLEARN_EXPORT void tab_printi(const char *title,
                                int value,
                                int ncol = 1,
                                const EJustify &justify = EJustify::fromKey("RIGHT"));
GSTLEARN_EXPORT void tab_print_rc(const char *title,
                                  int mode,
                                  int value,
                                  int ncol = 1,
                                  const EJustify &justify = EJustify::fromKey("RIGHT"));
GSTLEARN_EXPORT void tab_print_rowname(const char *string, int taille);
GSTLEARN_EXPORT void print_matrix(const char *title,
                                  int flag_limit,
                                  int bycol,
                                  int nx,
                                  int ny,
                                  const double *sel,
                                  const double *tab);
GSTLEARN_EXPORT void print_matrix(const char *title,
                                  int flag_limit,
                                  const AMatrix& mat);
GSTLEARN_EXPORT void print_trimat(const char *title,
                                  int mode,
                                  int neq,
                                  const double *tl);
GSTLEARN_EXPORT void print_imatrix(const char *title,
                                   int flag_limit,
                                   int bycol,
                                   int nx,
                                   int ny,
                                   const double *sel,
                                   const int *tab);
GSTLEARN_EXPORT void print_vector(const char *title,
                                  int flag_limit,
                                  int ntab,
                                  const double *tab);
GSTLEARN_EXPORT void print_vector(const char *title,
                                  int flag_limit,
                                  int ntab,
                                  const VectorDouble &tab);
GSTLEARN_EXPORT void print_ivector(const char *title,
                                   int flag_limit,
                                   int ntab,
                                   const int *itab);
GSTLEARN_EXPORT void print_ivector(const char *title,
                                   int flag_limit,
                                   int ntab,
                                   const VectorInt &itab);
