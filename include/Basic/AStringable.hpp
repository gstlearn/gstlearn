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

#include "geoslib_define.h"
#include "csparse_f.h" // Cannot use forward declaration for cs and don't know why!
/// TODO : include geoslib_f.h here makes the compilation crash !!!

class AStringable
{
public:
  AStringable() {};
  virtual ~AStringable() {};

  virtual String toString(int level = ITEST) const;

  virtual void display(int level = ITEST) const final;
};

// Set of functions regarding the printout
void   messageFlush(const String& string);
void   messerrFlush(const String& string);
void   messerr(const char *format,...);
void   message(const char *format,...);
String stringCompose(const char *format,...);
void   mesArg(const char *title, int current, int nmax, bool flagStartOne = false);
void   messageAbort(const char *format,...);
void   mestitle(int level,const char *format,...);
String toTitle(int level, const char* format, ...);
String toMatrix(const String& title,
                const VectorString& colnames,
                const VectorString& rownames,
                bool bycol,
                int ncols,
                int nrows,
                const VectorDouble &tab,
                bool flagOverride = false);
String toMatrix(const String& title,
                const VectorString& colnames,
                const VectorString& rownames,
                bool bycol,
                int ncols,
                int nrows,
                const VectorInt &tab,
                bool flagOverride = false);
String toMatrixSymmetric(const String& title,
                         const VectorString& colnames,
                         const VectorString& rownames,
                         bool bycol,
                         int ncols,
                         const VectorDouble &tab,
                         bool flagOverride = false);
String toMatrixDiagonal(const String& title,
                        const VectorString& colnames,
                        const VectorString& rownames,
                        int ncols,
                        const VectorDouble &tab,
                        bool flagOverride = false);
String toMatrixDiagCst(const String& title,
                       const VectorString& colnames,
                       const VectorString& rownames,
                       int ncols,
                       const VectorDouble &tab,
                       bool flagOverride = false);
String toMatrix(const String& title,
                const cs* A,
                bool  flagOverride = false);
String toVector(const String& title,
                const VectorDouble& tab);
String toVector(const String& title,
                const VectorInt& tab);
String toStr(const String& string, int justify = 1);
String toDouble(double value, int justify = 1);
String toInt(int value, int justify = 1);
String toInterval(double zmin, double zmax);
void   setFormatDecimalNumber(int number = 3);
void   setFormatColumnSize(int column = 10);
void   setFormatColnameSize(int column = 10);
void   setFormatMaxNRows(int maxnrows = 7);
void   setFormatMaxNCols(int maxncols = 7);
void   setFormatRCSize(int column = 3);
void   setFormatBatchNumber(int nbatch = 7);
