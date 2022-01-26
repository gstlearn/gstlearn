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
#include "geoslib_f.h"
#include "geoslib_old_f.h"
#include "Basic/AStringable.hpp"
#include "Basic/Vector.hpp"
#include "Basic/String.hpp"
#include "Basic/Utilities.hpp"
#include "Basic/File.hpp"
#include "Basic/OptCst.hpp"
#include "csparse_d.h"

#include <string>
#include <iostream>
#include <sstream>
#include <typeinfo>
#include <iomanip>
#include <stdio.h>
#include <stdarg.h>
#include <math.h>
#include <string.h>

#define JUSTIFY_LEFT  0
#define JUSTIFY_RIGHT 1

#define CASE_COL 0
#define CASE_ROW 1

static int _getColumnRank()
{
  return (int) OptCst::query(ECst::NTRANK);
}
static int _getColumnName()
{
  return (int) OptCst::query(ECst::NTNAME);
}
static int _getColumnSize()
{
  return (int) OptCst::query(ECst::NTCAR);
}
static int _getDecimalNumber()
{
  return (int) OptCst::query(ECst::NTDEC);
}
static double _getThresh()
{
  int ndec = (int) OptCst::query(ECst::NTDEC);
  // Recalculate threshold under which any small value must be displayed has 0.0
  double thresh = (0.5 * pow(10, - ndec));
  return thresh;
}
static int _getMaxNCols()
{
  return (int) OptCst::query(ECst::NTCOL);
}
static int _getMaxNRows()
{
  return (int) OptCst::query(ECst::NTROW);
}
static int _getNBatch()
{
  return (int) OptCst::query(ECst::NTBATCH);
}

AStringable::AStringable()
{
}

/**
 * Copy constructor: don't copy temporary file info
 */
AStringable::AStringable(const AStringable& /*r*/)
{
}
/**
 * Assignment operator: don't copy temporary file info
 */
AStringable& AStringable::operator=(const AStringable& /*r*/)
{
  return *this;
}


AStringable::~AStringable()
{
}

String AStringable::toString(const AStringFormat* strfmt) const
{
  std::stringstream sstr;
  sstr << "toString Not yet implemented for " << typeid(*this).name() << std::endl;
  return sstr.str();
}

std::stringstream _formatColumn(int justify, int localSize = 0)
{
  std::stringstream sstr;
  int size = (localSize > 0) ? localSize : _getColumnSize();
  sstr << std::fixed << std::setw(size) << std::setprecision(_getDecimalNumber());
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
  int truncSize = (localSize > 0) ? localSize : _getColumnSize();
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
  {
    // Prevent -0.00 : https://stackoverflow.com/a/12536500/3952924
    value = (ABS(value) < _getThresh()) ? 0. : value;
    sstr << value;
  }

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
  sstr << std::setw(_getColumnSize() - _getColumnRank() - 1) << std::right;
  if (icase == CASE_ROW)
  {
    if (!flagAdd)
      sstr << "[" << std::setw(_getColumnRank()) << value << ",]";
    else
      sstr << "[" << std::setw(_getColumnRank()) << value << "+]";
  }
  else
  {
    if (!flagAdd)
      sstr << "[," << std::setw(_getColumnRank()) << value << "]";
    else
      sstr << "[ " << std::setw(_getColumnRank()) << value << "]";
  }
  return sstr.str();
}

String _printColumnHeader(const VectorString& colnames,
                          int colfrom,
                          int colto,
                          int colSize = _getColumnSize())
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

String _printRowHeader(const VectorString& rownames, int iy, int rowSize = _getColumnSize())
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
 * Print a formatted message
 * @param format Output format
 * @param ...    Additional arguments
 */
void message(const char *format, ...)
{
  char str[LONG_SIZE];
  va_list ap;

  va_start(ap, format);
  // TODO : use non old_style functions
  (void) vsprintf(str, format, ap);
  va_end(ap);
  message_extern(str);

  return;
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
  char str[1000];
  va_list ap;

  va_start(ap, format);
  // TODO : use non old_style functions
  (void) vsprintf(str, format, ap);
  va_end(ap);

  message_extern(str);
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

  (void) gslStrcat(STRING, "\n");
  message_extern(STRING);

  /* Underline the string */

  (void) gslStrcpy(STRING, "");
  for (int i = 0; i < size; i++)
  {
    switch (level)
    {
      case 0:
        (void) gslStrcat(STRING, "=");
        break;

      case 1:
        (void) gslStrcat(STRING, "-");
        break;

      case 2:
        (void) gslStrcat(STRING, ".");
        break;
    }
  }
  (void) gslStrcat(STRING, "\n");
  message_extern(STRING);

  return;
}

/**
 * Conditionally print the progress of a procedure
 * @param string String to be printed
 * @param ntot   Total number of samples
 * @param iech   Rank of the current sample
 */
void mes_process(const char *string, int ntot, int iech)
{
  static int memo = 0;
  double ratio;
  int nproc, jech, percent;

  nproc = (int) OptCst::query(ECst::NPROC);
  if (nproc <= 0) return;
  jech = iech + 1;

  /* Calculate the current percentage */

  ratio = 100. * (double) jech / (double) ntot;
  percent = (int) (ratio / (double) nproc) * nproc;

  /* Conditional printout */

  if (percent != memo) message("%s : %d (percent)\n", string, percent);
  memo = percent;

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
  (void) gslStrcpy(STRING, "");
  for (int i = 0; i < size; i++)
  {
    switch (level)
    {
      case 0:
        (void) gslStrcat(STRING, "=");
        break;

      case 1:
        (void) gslStrcat(STRING, "-");
        break;

      case 2:
        (void) gslStrcat(STRING, ".");
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
void AStringable::display(const AStringFormat* strfmt) const
{
  message_extern(toString(strfmt).c_str());
}

void AStringable::display(int level) const
{
  AStringFormat sf(level);
  display(&sf);
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
  if (_getMaxNCols() > 0 && ncutil > _getMaxNCols() && !flagOverride) ncutil = _getMaxNCols();
  if (_getMaxNRows() > 0 && nrutil > _getMaxNRows() && !flagOverride) nrutil = _getMaxNRows();
  int npass = (int) ceil((double) ncutil / (double) _getNBatch());
  bool multi_row = nrutil > 1 || npass > 1;

  int colSize = 0;
  if (colnames.empty())
    colSize = _getColumnSize();
  else
  {
    colSize = MIN(_getColumnName(), getMaxStringSize(colnames) + 1);
    colSize = MAX(colSize, _getColumnSize());
  }
  int rowSize = 0;
  if (rownames.empty())
    rowSize = _getColumnSize();
  else
    rowSize = MAX(getMaxStringSize(rownames) + 1, _getColumnSize());

  /* Print the title (optional) */

  if (! title.empty())
  {
    sstr << title;
    if (multi_row) sstr << std::endl;
  }

  // Loop on the batches

  for (int ipass = 0; ipass < npass; ipass++)
  {
    int jdeb = ipass * _getNBatch();
    int jfin = MIN(jdeb + _getNBatch(), ncutil);

    /* Print the names of the columns and the column numbers */

    if (multi_row)
      sstr << _printColumnHeader(colnames, jdeb, jfin, colSize);

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
  if (_getMaxNCols() > 0 && ncutil > _getMaxNCols() && !flagOverride) ncutil = _getMaxNCols();
  if (_getMaxNRows() > 0 && nrutil > _getMaxNRows() && !flagOverride) nrutil = _getMaxNRows();
  int npass = (int) ceil((double) ncutil / (double) _getNBatch());
  bool multi_row = nrutil > 1 || npass > 1;

  int colSize = 0;
  if (colnames.empty())
    colSize = _getColumnSize();
  else
  {
    colSize = MIN(_getColumnName(), getMaxStringSize(colnames) + 1);
    colSize = MAX(colSize, _getColumnSize());
  }
  int rowSize = 0;
  if (rownames.empty())
    rowSize = _getColumnSize();
  else
    rowSize = MAX(getMaxStringSize(rownames) + 1, _getColumnSize());

  /* Print the title (optional) */

  if (! title.empty())
  {
    sstr << title;
    if (multi_row) sstr << std::endl;
  }

  // Loop on the batches

  for (int ipass = 0; ipass < npass; ipass++)
  {
    int jdeb = ipass * _getNBatch();
    int jfin = MIN(jdeb + _getNBatch(), ncutil);

    /* Print the names of the columns and the column numbers */

    if (multi_row)
      sstr << _printColumnHeader(colnames, jdeb, jfin, colSize);

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
  if (_getMaxNCols() > 0 && ncutil > _getMaxNCols() && !flagOverride) ncutil = _getMaxNCols();
  if (_getMaxNRows() > 0 && nrutil > _getMaxNRows() && !flagOverride) nrutil = _getMaxNRows();
  int npass = (int) ceil((double) ncutil / (double) _getNBatch());
  bool multi_row = nrutil > 1 || npass > 1;

  int colSize = 0;
  if (colnames.empty())
    colSize = _getColumnSize();
  else
  {
    colSize = MIN(_getColumnName(), getMaxStringSize(colnames) + 1);
    colSize = MAX(colSize, _getColumnSize());
  }
  int rowSize = 0;
  if (rownames.empty())
    rowSize = _getColumnSize();
  else
    rowSize = MAX(getMaxStringSize(rownames) + 1, _getColumnSize());

  /* Print the title (optional) */

  if (! title.empty())
  {
    sstr << title;
    if (multi_row) sstr << std::endl;
  }

  // Loop on the batches

  for (int ipass = 0; ipass < npass; ipass++)
  {
    int jdeb = ipass * _getNBatch();
    int jfin = MIN(jdeb + _getNBatch(), ncutil);

    /* Print the names of the columns and the column numbers */

    if (multi_row)
      sstr << _printColumnHeader(colnames, jdeb, jfin, colSize);

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
  if (_getMaxNCols() > 0 && ncutil > _getMaxNCols() && !flagOverride) ncutil = _getMaxNCols();
  if (_getMaxNRows() > 0 && nrutil > _getMaxNRows() && !flagOverride) nrutil = _getMaxNRows();
  int npass = (int) ceil((double) ncutil / (double) _getNBatch());
  bool multi_row = nrutil > 1 || npass > 1;

  int colSize = 0;
  if (colnames.empty())
    colSize = _getColumnSize();
  else
  {
    colSize = MIN(_getColumnName(), getMaxStringSize(colnames) + 1);
    colSize = MAX(colSize, _getColumnSize());
  }
  int rowSize = 0;
  if (rownames.empty())
    rowSize = _getColumnSize();
  else
    rowSize = MAX(getMaxStringSize(rownames) + 1, _getColumnSize());

  /* Print the title (optional) */

  if (! title.empty())
  {
    sstr << title;
    if (multi_row) sstr << std::endl;
  }

  // Loop on the batches

  for (int ipass = 0; ipass < npass; ipass++)
  {
    int jdeb = ipass * _getNBatch();
    int jfin = MIN(jdeb + _getNBatch(), ncutil);

    /* Print the names of the columns and the column numbers */

    if (multi_row)
      sstr << _printColumnHeader(colnames, jdeb, jfin, colSize);

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
  if (_getMaxNCols() > 0 && ncutil > _getMaxNCols() && !flagOverride) ncutil = _getMaxNCols();
  if (_getMaxNRows() > 0 && nrutil > _getMaxNRows() && !flagOverride) nrutil = _getMaxNRows();
  int npass = (int) ceil((double) ncutil / (double) _getNBatch());
  bool multi_row = nrutil > 1 || npass > 1;

  int colSize = 0;
  if (colnames.empty())
    colSize = _getColumnSize();
  else
  {
    colSize = MIN(_getColumnName(), getMaxStringSize(colnames) + 1);
    colSize = MAX(colSize, _getColumnSize());
  }
  int rowSize = 0;
  if (rownames.empty())
    rowSize = _getColumnSize();
  else
    rowSize = MAX(getMaxStringSize(rownames) + 1, _getColumnSize());

  /* Print the title (optional) */

  if (! title.empty())
  {
    sstr << title;
    if (multi_row) sstr << std::endl;
  }

  // Loop on the batches

  for (int ipass = 0; ipass < npass; ipass++)
  {
    int jdeb = ipass * _getNBatch();
    int jfin = MIN(jdeb + _getNBatch(), ncutil);

    /* Print the names of the columns and the column numbers */

    if (multi_row)
      sstr << _printColumnHeader(colnames, jdeb, jfin, colSize);

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
  if (_getMaxNCols() > 0 && ncutil > _getMaxNCols() && !flagOverride) ncutil = _getMaxNCols();
  if (_getMaxNRows() > 0 && nrutil > _getMaxNRows() && !flagOverride) nrutil = _getMaxNRows();
  int npass = (int) ceil((double) ncutil / (double) _getNBatch());
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
    int jdeb = ipass * _getNBatch();
    int jfin = MIN(jdeb + _getNBatch(), ncutil);

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
          sstr << _tabPrintString(".", JUSTIFY_RIGHT, _getColumnSize());
        else
          sstr << _tabPrintDouble(Ax[found], JUSTIFY_RIGHT, _getColumnSize());
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
  bool multi_row = ncutil > _getNBatch();

  /* Print the title (optional) */

  if (! title.empty())
  {
    sstr << title;
    if (multi_row) sstr << std::endl;
  }

  int lec = 0;
  if (multi_row) sstr << _printColumnHeader(VectorString(), 0, _getNBatch());

  for (int i = 0; i < ncutil; i += _getNBatch())
  {
    if (multi_row) sstr << _printRowHeader(VectorString(), i);

    for (int j = 0; j < _getNBatch(); j++)
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
 * Printout a list of vectors in a formatted manner
 * @param title Title of the printout (or empty string)
 * @param tab   Vector of vectors (real values) to be printed
 * @return The string (terminated with a newline)
 */
String toVector(const String& title, const VectorVectorDouble& tab)
{
  std::stringstream sstr;
  if (tab.empty()) return sstr.str();

  if (! title.empty())
    sstr << title << std::endl;

  for (int i = 0; i < (int) tab.size(); i++)
    sstr << toVector(String(), tab[i]);

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
  bool multi_row = ncutil > _getNBatch();

  /* Print the title (optional) */

  if (! title.empty())
  {
    sstr << title;
    if (multi_row) sstr << std::endl;
  }

  int lec = 0;
  if (multi_row) sstr << _printColumnHeader(VectorString(), 0, _getNBatch());

  for (int i = 0; i < ncutil; i += _getNBatch())
  {
    if (multi_row) sstr << _printRowHeader(VectorString(), i);

    for (int j = 0; j < _getNBatch(); j++)
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
