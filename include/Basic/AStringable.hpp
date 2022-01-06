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
#include "Basic/AStringFormat.hpp"

#include "geoslib_define.h"
// Put it in the header because inherited objects will need it
#include <sstream>

class cs;

class GSTLEARN_EXPORT AStringable
{
public:
  AStringable();
  AStringable(const AStringable& r);
  AStringable& operator=(const AStringable& r);
  virtual ~AStringable();

  virtual String toString(const AStringFormat* strfmt = nullptr) const;

  virtual void display(const AStringFormat* strfmt = nullptr) const final;
  virtual void display(int level) const final;
};

// Set of functions regarding the printout
GSTLEARN_EXPORT void   messageFlush(const String& string);
GSTLEARN_EXPORT void   messerrFlush(const String& string);
GSTLEARN_EXPORT void   messerr(const char *format,...);
GSTLEARN_EXPORT void   message(const char *format,...);
GSTLEARN_EXPORT void   mesArg(const char *title, int current, int nmax, bool flagStartOne = false);
GSTLEARN_EXPORT void   messageAbort(const char *format,...);
GSTLEARN_EXPORT void   mestitle(int level,const char *format,...);
GSTLEARN_EXPORT String toTitle(int level, const char* format, ...);
GSTLEARN_EXPORT String toMatrix(const String& title,
                                const VectorString& colnames,
                                const VectorString& rownames,
                                bool bycol,
                                int ncols,
                                int nrows,
                                const VectorDouble &tab,
                                bool flagOverride = false);
GSTLEARN_EXPORT String toMatrix(const String& title,
                                const VectorString& colnames,
                                const VectorString& rownames,
                                bool bycol,
                                int ncols,
                                int nrows,
                                const VectorInt &tab,
                                bool flagOverride = false);
GSTLEARN_EXPORT String toMatrixSymmetric(const String& title,
                                         const VectorString& colnames,
                                         const VectorString& rownames,
                                         bool bycol,
                                         int ncols,
                                         const VectorDouble &tab,
                                         bool flagOverride = false);
GSTLEARN_EXPORT String toMatrixDiagonal(const String& title,
                                        const VectorString& colnames,
                                        const VectorString& rownames,
                                        int ncols,
                                        const VectorDouble &tab,
                                        bool flagOverride = false);
GSTLEARN_EXPORT String toMatrixDiagCst(const String& title,
                                       const VectorString& colnames,
                                       const VectorString& rownames,
                                       int ncols,
                                       const VectorDouble &tab,
                                       bool flagOverride = false);
GSTLEARN_EXPORT String toMatrix(const String& title,
                                const cs* A,
                                bool  flagOverride = false);
GSTLEARN_EXPORT String toVector(const String& title,
                                const VectorDouble& tab);
GSTLEARN_EXPORT String toVector(const String& title,
                                const VectorVectorDouble& tab);
GSTLEARN_EXPORT String toVector(const String& title,
                                const VectorInt& tab);
GSTLEARN_EXPORT String toStr(const String& string, int justify = 1);
GSTLEARN_EXPORT String toDouble(double value, int justify = 1);
GSTLEARN_EXPORT String toInt(int value, int justify = 1);
GSTLEARN_EXPORT String toInterval(double zmin, double zmax);
GSTLEARN_EXPORT void   setFormatDecimalNumber(int number = 3);
GSTLEARN_EXPORT void   setFormatColumnSize(int column = 10);
GSTLEARN_EXPORT void   setFormatColnameSize(int column = 10);
GSTLEARN_EXPORT void   setFormatMaxNRows(int maxnrows = 7);
GSTLEARN_EXPORT void   setFormatMaxNCols(int maxncols = 7);
GSTLEARN_EXPORT void   setFormatRCSize(int column = 3);
GSTLEARN_EXPORT void   setFormatBatchNumber(int nbatch = 7);
