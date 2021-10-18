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
#include "Basic/AStringable.hpp"
#include "Basic/Vector.hpp"
#include "Basic/String.hpp"
#include "Basic/Utilities.hpp"
#include "geoslib_f.h"
#include "geoslib_old_f.h"

#include <string>
#include <iostream>
#include <sstream>
#include <typeinfo>
#include <iomanip>
#include <stdio.h>
#include <stdarg.h>

#define JUSTIFY_LEFT  0
#define JUSTIFY_RIGHT 1

#define CASE_COL 0
#define CASE_ROW 1

// TODO : move this as AStringable static members
static int _columnSize = 10;
static int _colnameSize = 12;
static int _nDec = 3;
static int _nRC = 3;
static int _maxNCols = 7;
static int _maxNRows = 7;
static int _nBatch = 7;

String AStringable::toString(int level) const
{
  std::stringstream sstr;
  sstr << "toString Not yet implemented for " << typeid(*this).name() << std::endl;
  return sstr.str();
}

std::stringstream _formatColumn(int justify, int localSize = 0)
{
  std::stringstream sstr;
  int size = (localSize > 0) ? localSize : _columnSize;
  sstr << std::fixed << std::setw(size) << std::setprecision(_nDec);
  if (justify == JUSTIFY_LEFT)
    sstr << std::left;
  else
    sstr << std::right;
  return sstr;
}

String _tabPrintString(const String& string, int justify, int localSize = 0)
{
  std::stringstream sstr = _formatColumn(justify, localSize);
  int size = static_cast<int> (string.size());
  int truncSize = (localSize > 0) ? localSize : _columnSize;
  if (size > truncSize)
  {
    // String must be truncated

    String strloc = string;
    strloc.erase(0, size - truncSize);
    strloc.replace(0, 2, " *");
    sstr << strloc;
  }
  else
  {
    sstr << string;
  }
  return sstr.str();
}

String _tabPrintDouble(double value, int justify, int localSize = 0)
{
  std::stringstream sstr = _formatColumn(justify, localSize);
  if (FFFF(value))
    sstr << "N/A";
  else
    sstr << value;

  return sstr.str();
}

String _tabPrintInt(int value, int justify, int localSize = 0)
{
  std::stringstream sstr = _formatColumn(justify, localSize);
  if (IFFFF(value))
    sstr << "N/A";
  else
    sstr << value;

  return sstr.str();
}

String _tabPrintRowColumn(int icase, int value, int flagAdd)
{
  std::stringstream sstr;
  sstr << std::setw(_columnSize - _nRC - 1) << std::right;
  if (icase == CASE_ROW)
  {
    if (!flagAdd)
      sstr << "[" << std::setw(_nRC) << value << ",]";
    else
      sstr << "[" << std::setw(_nRC) << value << "+]";
  }
  else
  {
    if (!flagAdd)
      sstr << "[," << std::setw(_nRC) << value << "]";
    else
      sstr << "[ " << std::setw(_nRC) << value << "]";
  }
  return sstr.str();
}

String _printColumnHeader(const VectorString& colnames,
                          int colfrom,
                          int colto,
                          int rowSize = _columnSize,
                          int colSize = _columnSize)
{
  std::stringstream sstr;
  if (!colnames.empty())
  {
    // By Names
    sstr << _tabPrintString(" ", JUSTIFY_RIGHT) << " ";
    for (int ix = colfrom; ix < colto; ix++)
      sstr << _tabPrintString(colnames[ix], JUSTIFY_RIGHT, colSize);
    sstr << std::endl;
  }
  else
  {
    // By Numbers
    sstr << _tabPrintString(" ", JUSTIFY_RIGHT) << " ";
    for (int ix = colfrom; ix < colto; ix++)
      sstr << _tabPrintRowColumn(CASE_COL, ix + 1, false);
    sstr << std::endl;
  }
  return sstr.str();
}

String _printRowHeader(const VectorString& rownames, int iy, int rowSize = _columnSize)
{
  std::stringstream sstr;
  if (!rownames.empty())
    sstr << _tabPrintString(rownames[iy], JUSTIFY_LEFT, rowSize);
  else
    sstr << _tabPrintRowColumn(CASE_ROW, iy + 1, false);
  return sstr.str();
}

String _printTrailer(int ncols, int nrows, int ncols_util, int nrows_util)
{
  std::stringstream sstr;
  if (ncols != ncols_util || nrows != nrows_util)
  {
    if (ncols == ncols_util)
      sstr << "(Ncols=" << ncols;
    else
      sstr << "(Ncols=" << ncols_util << "[from " << ncols << "]";

    if (nrows == nrows_util)
      sstr << ",Nrows=" << nrows << ")";
    else
      sstr << ",Nrows=" << nrows_util << "[from " << nrows << "])";
    sstr << std::endl;
  }
  return sstr.str();
}

/**
 * Setting the number of decimals for printing real value
 * @param number Number of decimals
 */
void setFormatDecimalNumber(int number)
{
  _nDec = number;
}

/**
 * Setting the number of characters per Column for printing an element
 * @param column Number of characters
 */
void setFormatColumnSize(int column)
{
  _columnSize = column;
}

/**
 * Setting the maximum number of characters for printing column name
 * @param column Number of characters
 */
void setFormatColnameSize(int column)
{
  _colnameSize = column;
}

/**
 * Define the number of columns for encoding the Row/Column rank
 * @param column Number of columns
 */
void setFormatRCSize(int column)
{
  _nRC = column;
}

/**
 * Set the maximum number of Rows in matrix printout
 * @param maxnrows Maximum number of -1
 */
void setFormatMaxNRows(int maxnrows)
{
  _maxNRows = maxnrows;
}
/**
 * Set the maximum number of Columns in matrix printout
 * @param maxncols Maximum number or -1
 */
void setFormatMaxNCols(int maxncols)
{
  _maxNCols = maxncols;
}
/**
 * Set the number of columns per Batch
 * @param nbatch Number of columns per batch
 */
void setFormatBatchNumber(int nbatch)
{
  _nBatch = nbatch;
}

/**
 * Print a formatted message
 * @param format Output format
 * @param ...    Additional arguments
 */
void message(const char *format, ...)
{
  char STRING[1000];
  va_list ap;

  va_start(ap, format);
  (void) vsprintf(STRING, format, ap);
  va_end(ap);
  message_extern(STRING);

  return;
}

String stringCompose(const char *format,...)
{
  std::stringstream sstr;
  char STRING[1000];
  va_list ap;

  va_start(ap, format);
  (void) vsprintf(STRING, format, ap);
  sstr << STRING;
  return sstr.str();
}

/**
 * When message has been collected as a String, this function produces it out
 * without passing through useless internal buffering
 * @param string String to be printed out
 */
void messageFlush(const String& string)
{
  message_extern(string.c_str());
}

/**
 * When the error message has been collected as a String, this function produces it out
 * without passing through useless internal buffering
 * @param string String to be produced
 * @remark This function is similar to messageFlush but dedicated to Errors
 */
void messerrFlush(const String& string)
{
  message_extern(string.c_str());
}

/**
 * Print Error message
 * @param format Output format
 * @param ...    Additional arguments
 */
void messerr(const char *format, ...)
{
  char STRING[1000];
  va_list ap;

  va_start(ap, format);
  (void) vsprintf(STRING, format, ap);
  va_end(ap);
  message_extern(STRING);
  message_extern("\n");

  return;
}

/**
 * Print a standard Error Message if an argument does not lie in Interval
 * @param title   Title to be printed
 * @param current Current value of the argument
 * @param nmax    Maximum (inclusive) possible value
 * @param flagStartOne True if the count starts at 1 (instead of 0)
 */
void mesArg(const char *title, int current, int nmax, bool flagStartOne)
{
  if (nmax <= 0)
    messerr("Error in %s (%d). No element of this type is recorded yet",title,current);
  else
  {
    if (flagStartOne)
      messerr("Error in %s (%d). Argument should lie within [1,%d]", title, current,
              nmax);
    else
      messerr("Error in %s (%d). Argument should lie within [0,%d[", title, current,
              nmax);
  }
}

/**
 * Print a message and underlines it with various formats
 * @param level  Level of the title
 * @param format Output format
 * @param ...    Additional arguments
 */
void mestitle(int level, const char *format, ...)
{
  char STRING[1000];
  va_list ap;

  message_extern("\n");
  va_start(ap, format);
  (void) vsprintf(STRING, format, ap);
  va_end(ap);
  int size = strlen(STRING);

  (void) strcat(STRING, "\n");
  message_extern(STRING);

  /* Underline the string */

  (void) strcpy(STRING, "");
  for (int i = 0; i < size; i++)
  {
    switch (level)
    {
      case 0:
        (void) strcat(STRING, "=");
        break;

      case 1:
        (void) strcat(STRING, "-");
        break;

      case 2:
        (void) strcat(STRING, ".");
        break;
    }
  }
  (void) strcat(STRING, "\n");
  message_extern(STRING);

  return;
}

/**
 * Print a message and underlines it with various formats
 * @param level  Level of the title
 * @param format Output format
 * @param ...    Additional arguments
 */
String toTitle(int level, const char* format, ...)
{
  std::stringstream sstr;
  char STRING[1000];
  va_list ap;

  sstr << std::endl;
  va_start(ap, format);
  (void) vsprintf(STRING, format, ap);
  va_end(ap);
  sstr << STRING << std::endl;

  /* Underline the string */

  int size = strlen(STRING);
  (void) strcpy(STRING, "");
  for (int i = 0; i < size; i++)
  {
    switch (level)
    {
      case 0:
        (void) strcat(STRING, "=");
        break;

      case 1:
        (void) strcat(STRING, "-");
        break;

      case 2:
        (void) strcat(STRING, ".");
        break;
    }
  }
  sstr << STRING << std::endl;

  return sstr.str();
}

/**
 * Function for aborting the API
 * @param format Fatal error format
 * @param ...    Additional arguments
 */
void messageAbort(const char *format, ...)
{
  char STRING[1000];
  va_list ap;

  va_start(ap, format);
  (void) vsprintf(STRING, format, ap);
  va_end(ap);
  message_extern("Abort : ");
  message_extern(STRING);
  message_extern("\n");
  exit_extern();
}

/**
 * Send the String to the display function
 */
void AStringable::display(int level) const
{
  message_extern(toString(level).c_str());
}

/**
 * Print the contents of a VectorDouble in a Matrix Form
 * @param title        Title of the printout
 * @param colnames     Names of the columns (optional)
 * @param rownames     Names of the rows (optional)
 * @param bycol        true if values as sorted by column; false otherwise
 * @param ncols        Number of columns
 * @param nrows        Number of rows
 * @param tab          VectorDouble containing the values
 * @param flagOverride true to override printout limitations
 */
String toMatrix(const String& title,
                const VectorString& colnames,
                const VectorString& rownames,
                bool bycol,
                int ncols,
                int nrows,
                const VectorDouble &tab,
                bool flagOverride)
{
  std::stringstream sstr;
  if (tab.empty() || ncols <= 0 || nrows <= 0) return sstr.str();

  /* Initializations */

  int ncutil = ncols;
  int nrutil = nrows;
  if (_maxNCols > 0 && ncutil > _maxNCols && !flagOverride) ncutil = _maxNCols;
  if (_maxNRows > 0 && nrutil > _maxNRows && !flagOverride) nrutil = _maxNRows;
  int npass = (int) ceil((double) ncutil / (double) _nBatch);
  bool multi_row = nrutil > 1 || npass > 1;

  int colSize = 0;
  if (colnames.empty())
    colSize = _columnSize;
  else
  {
    colSize = MIN(_colnameSize, getMaxStringSize(colnames) + 1);
    colSize = MAX(colSize, _columnSize);
  }
  int rowSize = 0;
  if (rownames.empty())
    rowSize = _columnSize;
  else
    rowSize = MAX(getMaxStringSize(rownames) + 1, _columnSize);

  /* Print the title (optional) */

  if (! title.empty())
  {
    sstr << title;
    if (multi_row) sstr << std::endl;
  }

  // Loop on the batches

  for (int ipass = 0; ipass < npass; ipass++)
  {
    int jdeb = ipass * _nBatch;
    int jfin = MIN(jdeb + _nBatch, ncutil);

    /* Print the names of the columns and the column numbers */

    if (multi_row)
      sstr << _printColumnHeader(colnames, jdeb, jfin, rowSize, colSize);

    /* Loop on the rows */

    for (int iy = 0; iy < nrutil; iy++)
    {
      if (multi_row) sstr << _printRowHeader(rownames, iy, rowSize);

      /* Loop on the columns */
      for (int ix = jdeb; ix < jfin; ix++)
      {
        int iad = (bycol) ? iy + nrows * ix : ix + ncols * iy;
        sstr << _tabPrintDouble(tab[iad], JUSTIFY_RIGHT, colSize);
      }
      sstr << std::endl;
    }
  }

  /* Print the trailer */

  sstr << _printTrailer(ncols, nrows, ncutil, nrutil);
  return sstr.str();
}

/**
 * Specific printout dedicated to square symmetrical matrices
 * Only the Upper left side is printed
 * @param title        Title of the printout
 * @param colnames     Names of the columns (optional)
 * @param rownames     Names of the rows (optional)
 * @param bycol        true if values as sorted by column; false otherwise
 * @param ncols        Number of columns = Number of rows
 * @param tab          VectorDouble containing the values (Dimension: n * (n+1) /2)
 * @param flagOverride Override the printout limitations
 * @return
 */
String toMatrixSymmetric(const String& title,
                const VectorString& colnames,
                const VectorString& rownames,
                bool  bycol,
                int   ncols,
                const VectorDouble &tab,
                bool  flagOverride)
{
  std::stringstream sstr;
  int nrows = ncols;
  if (tab.empty() || ncols <= 0 || nrows <= 0) return sstr.str();

  /* Initializations */

  int ncutil = ncols;
  int nrutil = nrows;
  if (_maxNCols > 0 && ncutil > _maxNCols && !flagOverride) ncutil = _maxNCols;
  if (_maxNRows > 0 && nrutil > _maxNRows && !flagOverride) nrutil = _maxNRows;
  int npass = (int) ceil((double) ncutil / (double) _nBatch);
  bool multi_row = nrutil > 1 || npass > 1;

  int colSize = 0;
  if (colnames.empty())
    colSize = _columnSize;
  else
  {
    colSize = MIN(_colnameSize, getMaxStringSize(colnames) + 1);
    colSize = MAX(colSize, _columnSize);
  }
  int rowSize = 0;
  if (rownames.empty())
    rowSize = _columnSize;
  else
    rowSize = MAX(getMaxStringSize(rownames) + 1, _columnSize);

  /* Print the title (optional) */

  if (! title.empty())
  {
    sstr << title;
    if (multi_row) sstr << std::endl;
  }

  // Loop on the batches

  for (int ipass = 0; ipass < npass; ipass++)
  {
    int jdeb = ipass * _nBatch;
    int jfin = MIN(jdeb + _nBatch, ncutil);

    /* Print the names of the columns and the column numbers */

    if (multi_row)
      sstr << _printColumnHeader(colnames, jdeb, jfin, rowSize, colSize);

    /* Loop on the rows */

    for (int iy = 0; iy < nrutil; iy++)
    {
      if (multi_row) sstr << _printRowHeader(rownames, iy, rowSize);

      /* Loop on the columns */
      for (int ix = jdeb; ix < jfin; ix++)
      {
        if (ix <= iy)
         {
           int iad = (bycol) ? iy + nrows * ix : ix + ncols * iy;
           sstr << _tabPrintDouble(tab[iad], JUSTIFY_RIGHT, colSize);
         }
         else
         {
           sstr << _tabPrintString(" ", JUSTIFY_RIGHT);
         }
      }
      sstr << std::endl;
    }
  }

  /* Print the trailer */

  sstr << _printTrailer(ncols, nrows, ncutil, nrutil);
  return sstr.str();
}

/**
 * Specific printout dedicated to square diagonal matrices
 * @param title        Title of the printout
 * @param colnames     Names of the columns (optional)
 * @param rownames     Names of the rows (optional)
 * @param ncols        Number of columns = Number of rows
 * @param tab          VectorDouble containing the values (Dimension: ncols)
 * @param flagOverride Override the printout limitations
 * @return
 */
String toMatrixDiagonal(const String& title,
                        const VectorString& colnames,
                        const VectorString& rownames,
                        int ncols,
                        const VectorDouble &tab,
                        bool flagOverride)
{
  std::stringstream sstr;
  int nrows = ncols;
  if (tab.empty() || ncols <= 0 || nrows <= 0) return sstr.str();

  /* Initializations */

  int ncutil = ncols;
  int nrutil = nrows;
  if (_maxNCols > 0 && ncutil > _maxNCols && !flagOverride) ncutil = _maxNCols;
  if (_maxNRows > 0 && nrutil > _maxNRows && !flagOverride) nrutil = _maxNRows;
  int npass = (int) ceil((double) ncutil / (double) _nBatch);
  bool multi_row = nrutil > 1 || npass > 1;

  int colSize = 0;
  if (colnames.empty())
    colSize = _columnSize;
  else
  {
    colSize = MIN(_colnameSize, getMaxStringSize(colnames) + 1);
    colSize = MAX(colSize, _columnSize);
  }
  int rowSize = 0;
  if (rownames.empty())
    rowSize = _columnSize;
  else
    rowSize = MAX(getMaxStringSize(rownames) + 1, _columnSize);

  /* Print the title (optional) */

  if (! title.empty())
  {
    sstr << title;
    if (multi_row) sstr << std::endl;
  }

  // Loop on the batches

  for (int ipass = 0; ipass < npass; ipass++)
  {
    int jdeb = ipass * _nBatch;
    int jfin = MIN(jdeb + _nBatch, ncutil);

    /* Print the names of the columns and the column numbers */

    if (multi_row)
      sstr << _printColumnHeader(colnames, jdeb, jfin, rowSize, colSize);

    /* Loop on the rows */

    for (int iy = 0; iy < nrutil; iy++)
    {
      if (multi_row) sstr << _printRowHeader(rownames, iy, rowSize);

      /* Loop on the columns */
      for (int ix = jdeb; ix < jfin; ix++)
      {
        if (ix == iy)
         {
           sstr << _tabPrintDouble(tab[ix], JUSTIFY_RIGHT, colSize);
         }
         else
         {
           sstr << _tabPrintString(" ", JUSTIFY_RIGHT);
         }
      }
      sstr << std::endl;
    }
  }

  /* Print the trailer */

  sstr << _printTrailer(ncols, nrows, ncutil, nrutil);
  return sstr.str();
}

/**
 * Specific printout dedicated to square diagonal matrices
 * @param title        Title of the printout
 * @param colnames     Names of the columns (optional)
 * @param rownames     Names of the rows (optional)
 * @param ncols        Number of columns = Number of rows
 * @param tab          VectorDouble containing the values (Dimension: 1)
 * @param flagOverride Override the printout limitations
 * @return
 */
String toMatrixDiagCst(const String& title,
                       const VectorString& colnames,
                       const VectorString& rownames,
                       int ncols,
                       const VectorDouble &tab,
                       bool flagOverride)
{
  std::stringstream sstr;
  int nrows = ncols;
  if (tab.empty() || ncols <= 0 || nrows <= 0) return sstr.str();

  /* Initializations */

  int ncutil = ncols;
  int nrutil = nrows;
  if (_maxNCols > 0 && ncutil > _maxNCols && !flagOverride) ncutil = _maxNCols;
  if (_maxNRows > 0 && nrutil > _maxNRows && !flagOverride) nrutil = _maxNRows;
  int npass = (int) ceil((double) ncutil / (double) _nBatch);
  bool multi_row = nrutil > 1 || npass > 1;

  int colSize = 0;
  if (colnames.empty())
    colSize = _columnSize;
  else
  {
    colSize = MIN(_colnameSize, getMaxStringSize(colnames) + 1);
    colSize = MAX(colSize, _columnSize);
  }
  int rowSize = 0;
  if (rownames.empty())
    rowSize = _columnSize;
  else
    rowSize = MAX(getMaxStringSize(rownames) + 1, _columnSize);

  /* Print the title (optional) */

  if (! title.empty())
  {
    sstr << title;
    if (multi_row) sstr << std::endl;
  }

  // Loop on the batches

  for (int ipass = 0; ipass < npass; ipass++)
  {
    int jdeb = ipass * _nBatch;
    int jfin = MIN(jdeb + _nBatch, ncutil);

    /* Print the names of the columns and the column numbers */

    if (multi_row)
      sstr << _printColumnHeader(colnames, jdeb, jfin, rowSize, colSize);

    /* Loop on the rows */

    for (int iy = 0; iy < nrutil; iy++)
    {
      if (multi_row) sstr << _printRowHeader(rownames, iy, rowSize);

      /* Loop on the columns */
      for (int ix = jdeb; ix < jfin; ix++)
      {
        if (ix == iy)
         {
           sstr << _tabPrintDouble(tab[0], JUSTIFY_RIGHT, colSize);
         }
         else
         {
           sstr << _tabPrintString(" ", JUSTIFY_RIGHT);
         }
      }
      sstr << std::endl;
    }
  }

  /* Print the trailer */

  sstr << _printTrailer(ncols, nrows, ncutil, nrutil);
  return sstr.str();
}

/**
 * Print the contents of a VectorDouble in a Matrix Form
 * @param title        Title of the printout
 * @param colnames     Names of the columns (optional)
 * @param rownames     Names of the rows (optional)
 * @param bycol        true if values as sorted by column; false otherwise
 * @param ncols        Number of columns
 * @param nrows        Number of rows
 * @param tab          VectorInt containing the values
 * @param flagOverride true to override printout limitations
 */
String toMatrix(const String& title,
                const VectorString& colnames,
                const VectorString& rownames,
                bool bycol,
                int ncols,
                int nrows,
                const VectorInt &tab,
                bool flagOverride)
{
  std::stringstream sstr;
  if (tab.empty() || ncols <= 0 || nrows <= 0) return sstr.str();

  /* Initializations */

  int ncutil = ncols;
  int nrutil = nrows;
  if (_maxNCols > 0 && ncutil > _maxNCols && !flagOverride) ncutil = _maxNCols;
  if (_maxNRows > 0 && nrutil > _maxNRows && !flagOverride) nrutil = _maxNRows;
  int npass = (int) ceil((double) ncutil / (double) _nBatch);
  bool multi_row = nrutil > 1 || npass > 1;

  int colSize = 0;
  if (colnames.empty())
    colSize = _columnSize;
  else
  {
    colSize = MIN(_colnameSize, getMaxStringSize(colnames) + 1);
    colSize = MAX(colSize, _columnSize);
  }
  int rowSize = 0;
  if (rownames.empty())
    rowSize = _columnSize;
  else
    rowSize = MAX(getMaxStringSize(rownames) + 1, _columnSize);

  /* Print the title (optional) */

  if (! title.empty())
  {
    sstr << title;
    if (multi_row) sstr << std::endl;
  }

  // Loop on the batches

  for (int ipass = 0; ipass < npass; ipass++)
  {
    int jdeb = ipass * _nBatch;
    int jfin = MIN(jdeb + _nBatch, ncutil);

    /* Print the names of the columns and the column numbers */

    if (multi_row)
      sstr << _printColumnHeader(colnames, jdeb, jfin, rowSize, colSize);

    /* Loop on the rows */

    for (int iy = 0; iy < nrutil; iy++)
    {
      if (multi_row) sstr << _printRowHeader(rownames, iy, rowSize);

      /* Loop on the columns */

      for (int ix = jdeb; ix < jfin; ix++)
      {
        int iad = (bycol) ? iy + nrows * ix : ix + ncols * iy;
        sstr << _tabPrintInt(tab[iad], JUSTIFY_RIGHT, colSize);
      }
      sstr << std::endl;
    }
  }

  /* Print the trailer */

  sstr << _printTrailer(ncols, nrows, ncutil, nrutil);
  return sstr.str();
}

String toMatrix(const String& title, const cs* A, bool flagOverride)
{
  std::stringstream sstr;
  int nrows, ncols;

  if (!A) return sstr.str();
  if (A->nz >= 0) return sstr.str();

  int nrutil = nrows = A->m;
  int ncutil = ncols = A->n;
  if (_maxNCols > 0 && ncutil > _maxNCols && !flagOverride) ncutil = _maxNCols;
  if (_maxNRows > 0 && nrutil > _maxNRows && !flagOverride) nrutil = _maxNRows;
  int npass = (int) ceil((double) ncutil / (double) _nBatch);
  bool multi_row = nrutil > 1 || npass > 1;

  int* Ap = A->p;
  int* Ai = A->i;
  double* Ax = A->x;

  /* Print the title (optional) */

  if (! title.empty())
  {
    sstr << title;
    if (multi_row) sstr << std::endl;
  }

  /* Loop on the passes */

  for (int ipass = 0; ipass < npass; ipass++)
  {
    int jdeb = ipass * _nBatch;
    int jfin = MIN(jdeb + _nBatch, ncutil);

    /* Title of the columns */

    if (multi_row) sstr << _printColumnHeader(VectorString(), jdeb, jfin);

    /* Loop on the lines */

    for (int iy = 0; iy < nrutil; iy++)
    {
      if (multi_row) sstr << _printRowHeader(VectorString(), iy);

      /* Loop on the columns */

      for (int ix = jdeb; ix < jfin; ix++)
      {

        /* Search for the correct line number */

        int found = -1;
        for (int p = Ap[ix]; p < Ap[ix + 1] && found < 0; p++)
        {
          if (Ai[p] == iy) found = p;
        }

        if (found < 0)
          sstr << _tabPrintString(".", JUSTIFY_RIGHT, _columnSize);
        else
          sstr << _tabPrintDouble(Ax[found], JUSTIFY_RIGHT, _columnSize);
      }
      sstr << std::endl;
    }
    sstr << std::endl;
  }

  // Print the trailer
  sstr << _printTrailer(ncols, nrows, ncutil, nrutil);
  return sstr.str();
}

/**
 * Printout a vector in a formatted manner
 * @param title Title of the printout (or empty string)
 * @param tab   Vector (real values) to be printed
 * @return The string (terminated with a newline)
 */
String toVector(const String& title, const VectorDouble& tab)
{
  std::stringstream sstr;
  if (tab.empty()) return sstr.str();

  int ncutil = static_cast<int> (tab.size());
  bool multi_row = ncutil > _nBatch;

  /* Print the title (optional) */

  if (! title.empty())
  {
    sstr << title;
    if (multi_row) sstr << std::endl;
  }

  int lec = 0;
  if (multi_row) sstr << _printColumnHeader(VectorString(), 0, _nBatch);

  for (int i = 0; i < ncutil; i += _nBatch)
  {
    if (multi_row) sstr << _printRowHeader(VectorString(), i);

    for (int j = 0; j < _nBatch; j++)
    {
      if (lec >= ncutil) continue;
      sstr << toDouble(tab[lec]);
      lec++;
    }
    sstr << std::endl;
  }
  return sstr.str();
}

/**
 * Printout a vector in a formatted manner
 * @param title Title of the printout (or empty string)
 * @param tab   Vector (integer values) to be printed
 * @return The string (terminated with a newline)
 */
String toVector(const String& title, const VectorInt& tab)
{
  std::stringstream sstr;
  if (tab.empty()) return sstr.str();

  int ncutil = static_cast<int> (tab.size());
  bool multi_row = ncutil > _nBatch;

  /* Print the title (optional) */

  if (! title.empty())
  {
    sstr << title;
    if (multi_row) sstr << std::endl;
  }

  int lec = 0;
  if (multi_row) sstr << _printColumnHeader(VectorString(), 0, _nBatch);

  for (int i = 0; i < ncutil; i += _nBatch)
  {
    if (multi_row) sstr << _printRowHeader(VectorString(), i);

    for (int j = 0; j < _nBatch; j++)
    {
      if (lec >= ncutil) continue;
      sstr << toInt(tab[lec]);
      lec++;
    }
    sstr << std::endl;
  }
  return sstr.str();
}

String toStr(const String& string, int justify)
{
  std::stringstream sstr;
  sstr << _tabPrintString(string, justify);
  return sstr.str();
}

String toDouble(double value, int justify)
{
  std::stringstream sstr;
  sstr << _tabPrintDouble(value, justify);
  return sstr.str();
}

String toInt(int value, int justify)
{
  std::stringstream sstr;
  sstr << _tabPrintInt(value, justify);
  return sstr.str();
}

String toInterval(double zmin, double zmax)
{
  std::stringstream sstr;

  sstr << "Bounded in [";
  if (FFFF(zmin))
    sstr << "N/A";
  else
    sstr << zmin;
  message(" ; ");
  if (FFFF(zmax))
    sstr << "N/A";
  else
    sstr << zmax;
  sstr << "]" << std::endl;

  return sstr.str();
}
