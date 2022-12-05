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
#include "geoslib_old_f.h"
#include "csparse_d.h"

#include "Enum/EJustify.hpp"

#include "Basic/AStringable.hpp"
#include "Basic/VectorNumT.hpp"
#include "Basic/String.hpp"
#include "Basic/Utilities.hpp"
#include "Basic/File.hpp"
#include "Basic/OptCst.hpp"

#include <string>
#include <iostream>
#include <sstream>
#include <typeinfo>
#include <iomanip>
#include <stdio.h>
#include <stdarg.h>
#include <math.h>
#include <string.h>

#define CASE_DOUBLE 0
#define CASE_REAL   1
#define CASE_INT    2
#define CASE_COL    3
#define CASE_ROW    4

static char FORMAT[100];
static char DECODE[100];
static char TABSTR[100];

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
static void _buildFormat(int mode)
{
  switch (mode)
  {
    case CASE_INT:
      (void) gslSPrintf(FORMAT, "%%%dd", (int) OptCst::query(ECst::NTCAR));
      break;

    case CASE_REAL:
      (void) gslSPrintf(FORMAT, "%%%d.%dlf", (int) OptCst::query(ECst::NTCAR),
                        (int) OptCst::query(ECst::NTDEC));
      break;

    case CASE_DOUBLE:
      (void) gslSPrintf(FORMAT, "%%%d.%dlg", (int) OptCst::query(ECst::NTCAR),
                        (int) OptCst::query(ECst::NTDEC));
      break;

    case CASE_COL:
      (void) gslSPrintf(FORMAT, "[,%%%dd]", (int) OptCst::query(ECst::NTCAR) - 3);
      break;

    case CASE_ROW:
      (void) gslSPrintf(FORMAT, "[%%%dd,]", (int) OptCst::query(ECst::NTCAR) - 3);
      break;
  }
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

String AStringable::toString(const AStringFormat* /*strfmt*/) const
{
  std::stringstream sstr;
  sstr << "toString Not yet implemented for " << typeid(*this).name() << std::endl;
  return sstr.str();
}

std::stringstream _formatColumn(const EJustify& justify, int localSize = 0)
{
  std::stringstream sstr;
  int size = (localSize > 0) ? localSize : _getColumnSize();
  sstr << std::fixed << std::setw(size) << std::setprecision(_getDecimalNumber());
  if (justify == EJustify::LEFT)
    sstr << std::left;
  else
    sstr << std::right;
  return sstr;
}

String _tabPrintString(const String& string,
                       const EJustify& justify,
                       int localSize = 0)
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

String _tabPrintDouble(double value, const EJustify& justify, int localSize = 0)
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

String _tabPrintInt(int value, const EJustify& justify, int localSize = 0)
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
    sstr << _tabPrintString(" ", EJustify::RIGHT) << " ";
    for (int ix = colfrom; ix < colto; ix++)
      sstr << _tabPrintString(colnames[ix], EJustify::RIGHT, colSize);
    sstr << std::endl;
  }
  else
  {
    // By Numbers
    sstr << _tabPrintString(" ", EJustify::RIGHT) << " ";
    for (int ix = colfrom; ix < colto; ix++)
      sstr << _tabPrintRowColumn(CASE_COL, ix, false);
    sstr << std::endl;
  }
  return sstr.str();
}

String _printRowHeader(const VectorString& rownames, int iy, int rowSize = _getColumnSize())
{
  std::stringstream sstr;
  if (!rownames.empty())
    sstr << _tabPrintString(rownames[iy], EJustify::LEFT, rowSize);
  else
    sstr << _tabPrintRowColumn(CASE_ROW, iy, false);
  return sstr.str();
}

String _printTrailer(int ncols, int nrows, int ncols_util, int nrows_util)
{
  std::stringstream sstr;

  bool used = false;
  if (ncols != ncols_util)
  {
    used = true;
    if (ncols == ncols_util)
      sstr << "(Ncols=" << ncols;
    else
      sstr << "(Ncols=" << ncols_util << "[from " << ncols << "]";
  }

  if (nrows != nrows_util)
  {
    used = true;
    if (nrows == nrows_util)
      sstr << ",Nrows=" << nrows << ")";
    else
      sstr << ",Nrows=" << nrows_util << "[from " << nrows << "])";
  }

  if (used) sstr << std::endl;
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
  int size = (int) strlen(STRING);

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

  int size = (int) strlen(STRING);
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
        sstr << _tabPrintDouble(tab[iad], EJustify::RIGHT, colSize);
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
String toMatrixSymmetric(const String &title,
                         const VectorString &colnames,
                         const VectorString &rownames,
                         bool bycol,
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
        if (ix <= iy)
         {
           int iad = (bycol) ? iy + nrows * ix : ix + ncols * iy;
           sstr << _tabPrintDouble(tab[iad], EJustify::RIGHT, colSize);
         }
         else
         {
           sstr << _tabPrintString(" ", EJustify::RIGHT);
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
           sstr << _tabPrintDouble(tab[ix], EJustify::RIGHT, colSize);
         }
         else
         {
           sstr << _tabPrintString(" ", EJustify::RIGHT);
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
           sstr << _tabPrintDouble(tab[0], EJustify::RIGHT, colSize);
         }
         else
         {
           sstr << _tabPrintString(" ", EJustify::RIGHT);
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
        sstr << _tabPrintInt(tab[iad], EJustify::RIGHT, colSize);
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
          sstr << _tabPrintString(".", EJustify::RIGHT, _getColumnSize());
        else
          sstr << _tabPrintDouble(Ax[found], EJustify::RIGHT, _getColumnSize());
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
 * @param flagOverride true to override printout limitations
 * @return The string (terminated with a newline)
 */
String toVector(const String& title, const VectorDouble& tab, bool flagOverride)
{
  std::stringstream sstr;
  if (tab.empty()) return sstr.str();

  int ncols = static_cast<int> (tab.size());
  int ncutil = ncols;
  if (_getMaxNCols() > 0 && ncutil > _getMaxNCols() && !flagOverride) ncutil = _getMaxNCols();
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

  // Print the trailer
  sstr << _printTrailer(ncols, 0, ncutil, 0);

  return sstr.str();
}

/**
 * Printout a list of vectors in a formatted manner
 * @param title Title of the printout (or empty string)
 * @param tab   Vector of vectors (real values) to be printed
 * @param flagOverride true to override printout limitations
 * @return The string (terminated with a newline)
 */
String toVector(const String& title, const VectorVectorDouble& tab, bool flagOverride)
{
  std::stringstream sstr;
  if (tab.empty()) return sstr.str();

  if (! title.empty())
    sstr << title << std::endl;

  int nrows = (int) tab.size();
  int nrutil = nrows;
  if (_getMaxNRows() > 0 && nrutil > _getMaxNRows() && !flagOverride) nrutil = _getMaxNRows();

  for (int i = 0; i < nrutil; i++)
    sstr << toVector(String(), tab[i], flagOverride);

  // Print the trailer
  sstr << _printTrailer(0, nrows, 0, nrutil);

  return sstr.str();
}

String toVector(const String& title, const VectorString& tab, bool flagOverride)
{
  std::stringstream sstr;
  if (tab.empty()) return sstr.str();

  int ncols = static_cast<int> (tab.size());
  int ncutil = ncols;
  if (_getMaxNCols() > 0 && ncutil > _getMaxNCols() && !flagOverride) ncutil = _getMaxNCols();
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
      sstr << tab[lec];
      lec++;
    }
    sstr << std::endl;
  }

  // Print the trailer
  sstr << _printTrailer(ncols, 0, ncutil, 0);

  return sstr.str();
}

/**
 * Printout a vector in a formatted manner
 * @param title Title of the printout (or empty string)
 * @param tab   Vector (integer values) to be printed
 * @param flagOverride true to override printout limitations
 * @return The string (terminated with a newline)
 */
String toVector(const String& title, const VectorInt& tab, bool flagOverride)
{
  std::stringstream sstr;
  if (tab.empty()) return sstr.str();

  int ncols = static_cast<int> (tab.size());
  int ncutil = ncols;
  if (_getMaxNCols() > 0 && ncutil > _getMaxNCols() && !flagOverride) ncutil = _getMaxNCols();
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

  // Print the trailer
  sstr << _printTrailer(ncols, 0, ncutil, 0);

  return sstr.str();
}

String toStr(const String& string, const EJustify& justify)
{
  std::stringstream sstr;
  sstr << _tabPrintString(string, justify);
  return sstr.str();
}

String toDouble(double value, const EJustify& justify)
{
  std::stringstream sstr;
  sstr << _tabPrintDouble(value, justify);
  return sstr.str();
}

VectorString toVectorDouble(const VectorDouble& values, const EJustify& justify)
{
  VectorString strings;
  for (int i = 0; i < (int) values.size(); i++)
    strings.push_back(toDouble(values[i], justify));
  return strings;
}

String toInt(int value, const EJustify& justify)
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

/****************************************************************************/
/*!
 **  Tabulated printout of a string
 **
 ** \param[in]  title    optional title (NULL if not defined)
 ** \param[in]  string   String to be written
 ** \param[in]  ncol     number of columns for the printout
 ** \param[in]  justify  justification flag
 **                      (EJustify::LEFT, EJustify::CENTER or EJustify::RIGHT)
 **
 *****************************************************************************/
void tab_prints(const char *title,
                const char *string,
                int ncol,
                const EJustify &justify)
{
  int taille = (1 + (int) OptCst::query(ECst::NTCAR)) * ncol;
  int size = static_cast<int>(strlen(string));
  int neff = MIN(taille, size);
  int nrst = taille - neff;
  int n1 = nrst / 2;
  int n2 = taille - size - n1;

  /* Encode the title (if defined) */

  if (title != NULL) message("%s", title);

  /* Blank the string out */

  (void) gslStrcpy(TABSTR, "");

  /* Switch according to the justification */

  switch (justify.toEnum())
  {
    case EJustify::E_LEFT:
      (void) gslStrncpy(TABSTR, string, neff);
      TABSTR[neff] = '\0';
      for (int i = 0; i < nrst; i++)
        (void) gslStrcat(TABSTR, " ");
      break;

    case EJustify::E_CENTER:
      for (int i = 0; i < n1; i++)
        (void) gslStrcat(TABSTR, " ");
      (void) gslStrncpy(&TABSTR[n1], string, neff);
      TABSTR[n1 + neff] = '\0';
      for (int i = 0; i < n2; i++)
        (void) gslStrcat(TABSTR, " ");
      break;

    case EJustify::E_RIGHT:
      for (int i = 0; i < nrst; i++)
        (void) gslStrcat(TABSTR, " ");
      (void) gslStrncpy(&TABSTR[nrst], string, neff);
      TABSTR[nrst + neff] = '\0';
      break;
  }
  message(TABSTR);
}

/****************************************************************************/
/*!
 **  Tabulated printout of a real value
 **
 ** \param[in]  title    optional title (NULL if not defined)
 ** \param[in]  value    Value to be written
 ** \param[in]  ncol     number of columns for the printout
 ** \param[in]  justify  justification flag
 **                      (EJustify::LEFT, EJustify::CENTER or EJustify::RIGHT)
 **
 *****************************************************************************/
void tab_printg(const char *title,
                double value,
                int ncol,
                const EJustify &justify)
{
  _buildFormat(CASE_REAL);

  if (FFFF(value))
    (void) gslStrcpy(DECODE, "N/A");
  else
  {
    // Prevent -0.00 : https://stackoverflow.com/a/12536500/3952924
    value = (ABS(value) < _getThresh()) ? 0. : value;
    (void) gslSPrintf(DECODE, FORMAT, value);
  }
  tab_prints(title, DECODE, ncol, justify);
}

/****************************************************************************/
/*!
 **  Tabulated printout of a double value
 **
 ** \param[in]  title    optional title (NULL if not defined)
 ** \param[in]  value    Value to be written
 ** \param[in]  ncol     number of columns for the printout
 ** \param[in]  justify  justification flag
 **                      (EJustify::LEFT, EJustify::CENTER or EJustify::RIGHT)
 **
 *****************************************************************************/
void tab_printd(const char *title,
                double value,
                int ncol,
                const EJustify &justify)
{
  _buildFormat(CASE_DOUBLE);

  if (FFFF(value))
    (void) gslStrcpy(DECODE, "N/A");
  else
    (void) gslSPrintf(DECODE, FORMAT, value);

  tab_prints(title, DECODE, ncol, justify);
}

/****************************************************************************/
/*!
 **  Tabulated printout of an integer value
 **
 ** \param[in]  title    optional title (NULL if not defined)
 ** \param[in]  value    Value to be written
 ** \param[in]  ncol     number of columns for the printout
 ** \param[in]  justify  justification flag
 **                      (EJustify::LEFT, EJustify::CENTER or EJustify::RIGHT)
 **
 *****************************************************************************/
void tab_printi(const char *title, int value, int ncol, const EJustify &justify)
{
  _buildFormat(CASE_INT);

  if (IFFFF(value))
    (void) gslStrcpy(DECODE, "N/A");
  else
    (void) gslSPrintf(DECODE, FORMAT, value);

  tab_prints(title, DECODE, ncol, justify);
}

/****************************************************************************/
/*!
 **  Tabulated printout of a row or column value
 **
 ** \param[in]  title    optional title (NULL if not defined)
 ** \param[in]  mode     CASE_ROW or CASE_COL
 ** \param[in]  value    Value to be written
 ** \param[in]  ncol     number of columns for the printout
 ** \param[in]  justify  justification flag
 **                      (EJustify::LEFT, EJustify::CENTER or EJustify::RIGHT)
 **
 *****************************************************************************/
void tab_print_rc(const char *title,
                  int mode,
                  int value,
                  int ncol,
                  const EJustify &justify)
{
  _buildFormat(mode);

  (void) gslSPrintf(DECODE, FORMAT, value);
  string_strip_blanks(DECODE, 0);
  tab_prints(title, DECODE, ncol, justify);
}

/****************************************************************************/
/*!
 **  Tabulated printout of a string (character size provided)
 **
 ** \param[in]  string   String to be written
 ** \param[in]  taille   Number of characters
 **
 ** \remarks The string is printed (left-adjusted) on 'taille' characters
 **
 *****************************************************************************/
void tab_print_rowname(const char *string, int taille)
{
  int size = static_cast<int>(strlen(string));
  int neff = MIN(taille, size);
  int nrst = taille - neff;

  /* Blank the string out */

  (void) gslStrcpy(TABSTR, "");
  (void) gslStrncpy(TABSTR, string, neff);
  TABSTR[neff] = '\0';
  for (int i = 0; i < nrst; i++)
    (void) gslStrcat(TABSTR, " ");
  message(TABSTR);
}

/****************************************************************************/
/*!
 **  Tabulated printout of a matrix
 **
 ** \param[in]  title  Title (Optional)
 ** \param[in]  flag_limit  option for the limits
 ** \li                      1 if limits must be applied
 ** \li                      0 if the whole matrix is printed
 ** \param[in]  bycol  1 if values in 'tab' are sorted by column, 0 otherwise
 ** \param[in]  nx     number of columns in the matrix
 ** \param[in]  ny     number of rows in the matrix
 ** \param[in]  sel    array of selection or NULL
 ** \param[in]  tab    array containing the matrix
 **
 ** \remarks The order of the dimension (nx,ny) is opposite
 ** \remarks of the one used in R-packages where dim[1]=nrow and dim[2]=ncol
 **
 *****************************************************************************/
void print_matrix(const char *title,
                  int flag_limit,
                  int bycol,
                  int nx,
                  int ny,
                  const double *sel,
                  const double *tab)
{
  if (tab == nullptr || nx <= 0 || ny <= 0) return;
  int nx_util = (flag_limit && (int) OptCst::query(ECst::NTCOL) > 0) ?
      MIN((int) OptCst::query(ECst::NTCOL), nx) : nx;
  int ny_util = (flag_limit && (int) OptCst::query(ECst::NTROW) > 0) ?
      MIN((int) OptCst::query(ECst::NTROW), ny) : ny;
  int multi_row = (ny > 1 || title == NULL);

  /* Print the title (optional) */

  if (title != NULL)
  {
    if (multi_row)
      message("%s\n", title);
    else
      message("%s ", title);
  }

  /* Print the header */

  if (multi_row)
  {
    tab_prints(NULL, " ");
    for (int ix = 0; ix < nx_util; ix++)
      tab_print_rc(NULL, CASE_COL, ix + 1);
    message("\n");
  }

  /* Print the contents of the array */

  int ny_done = 0;
  for (int iy = 0; iy < ny; iy++)
  {
    if (sel != nullptr && !sel[iy]) continue;
    ny_done++;
    if (ny_done > ny_util) break;
    if (multi_row) tab_print_rc(NULL, CASE_ROW, iy + 1);
    for (int ix = 0; ix < nx_util; ix++)
    {
      int iad = (bycol) ? iy + ny * ix : ix + nx * iy;
      tab_printg(NULL, tab[iad]);
    }
    message("\n");
  }

  /* Print the trailor */

  if (nx != nx_util || ny != ny_util)
  {
    if (nx == nx_util)
      message("(Ncol=%d", nx);
    else
      message("(Ncol=%d[from %d]", nx_util, nx);

    if (ny == ny_util)
      message(",Nrow=%d)", ny);
    else
      message(",Nrow=%d[from %d])", ny_util, ny);
    message("\n");
  }
}

/****************************************************************************/
/*!
 **  Tabulated printout of a upper triangular matrix
 **
 ** \param[in]  title  Title (Optional)
 ** \param[in]  mode   1 if the matrix is stored linewise
 **                    2 if the matrix is stored columnwise
 ** \param[in]  neq    size of the matrix
 ** \param[in]  tl     array containing the upper triangular matrix
 **
 ** \remarks The ordering (compatible with matrix_solve is mode==2)
 **
 *****************************************************************************/
void print_trimat(const char *title, int mode, int neq, const double *tl)
{
#define TRI(i)        (((i) * ((i) + 1)) / 2)
#define TL1(i,j)      (tl[(j)*neq+(i)-TRI(j)])  /* only for i >= j */
#define TL2(i,j)      (tl[TRI(i)+(j)])          /* only for i >= j */

  /* Initializations */

  if (tl == nullptr || neq <= 0) return;

  /* Print the title (optional) */

  if (title != NULL) message("%s\n", title);

  /* Print the header */

  tab_prints(NULL, " ");
  for (int ix = 0; ix < neq; ix++)
    tab_print_rc(NULL, CASE_COL, ix + 1);
  message("\n");

  /* Print the contents of the array */

  for (int iy = 0; iy < neq; iy++)
  {
    tab_print_rc(NULL, CASE_ROW, iy + 1);
    for (int ix = 0; ix < neq; ix++)
    {
      if (ix >= iy)
      {
        if (mode == 1)
          tab_printg(NULL, TL1(ix, iy));
        else
          tab_printg(NULL, TL2(ix, iy));
      }
      else
        tab_prints(NULL, " ");
    }
    message("\n");
  }
#undef TRI
#undef TL1
#undef TL2
}

/****************************************************************************/
/*!
 **  Tabulated printout of a matrix (integer version)
 **
 ** \param[in]  title  Title (Optional)
 ** \param[in]  flag_limit  option for the limits
 ** \li                      1 if limits must be applied
 ** \li                      0 if the whole matrix is printed
 ** \param[in]  bycol  1 if values in 'tab' are sorted by column, 0 otherwise
 ** \param[in]  nx     number of columns in the matrix
 ** \param[in]  ny     number of rows in the matrix
 ** \param[in]  sel    array of selection or NULL
 ** \param[in]  tab    array containing the matrix
 **
 *****************************************************************************/
void print_imatrix(const char *title,
                   int flag_limit,
                   int bycol,
                   int nx,
                   int ny,
                   const double *sel,
                   const int *tab)
{
  if (tab == nullptr || nx <= 0 || ny <= 0) return;
  int nx_util = (flag_limit && (int) OptCst::query(ECst::NTCOL) > 0) ?
      MIN((int) OptCst::query(ECst::NTCOL), nx) : nx;
  int ny_util = (flag_limit && (int) OptCst::query(ECst::NTROW) > 0) ?
      MIN((int) OptCst::query(ECst::NTROW), ny) : ny;
  int multi_row = (ny > 1 || title == NULL);

  /* Print the title (optional) */

  if (title != NULL)
  {
    if (multi_row)
      message("%s\n", title);
    else
      message("%s ", title);
  }

  /* Print the header */

  if (multi_row)
  {
    tab_prints(NULL, " ");
    for (int ix = 0; ix < nx_util; ix++)
      tab_print_rc(NULL, CASE_COL, ix + 1);
    message("\n");
  }

  /* Print the contents of the array */

  int ny_done = 0;
  for (int iy = 0; iy < ny; iy++)
  {
    if (sel != nullptr && !sel[iy]) continue;
    ny_done++;
    if (ny_done > ny_util) break;
    if (multi_row) tab_print_rc(NULL, CASE_ROW, iy + 1);
    for (int ix = 0; ix < nx_util; ix++)
    {
      int iad = (bycol) ? iy + ny * ix : ix + nx * iy;
      tab_printi(NULL, tab[iad]);
    }
    message("\n");
  }

  /* Print the trailing part */

  if (nx != nx_util || ny != ny_util)
  {
    if (nx == nx_util)
      message("(Ncol=%d", nx);
    else
      message("(Ncol=%d[from %d]", nx_util, nx);

    if (ny == ny_util)
      message(",Nrow=%d)", ny);
    else
      message(",Nrow=%d[from %d])", ny_util, ny);
    message("\n");
  }
}

/****************************************************************************/
/*!
 **  Print a vector of real values in a matrix form
 **
 ** \param[in]  title      Title (Optional)
 ** \param[in]  flag_limit 1 if NTCOL is used; 0 otherwise
 ** \param[in]  ntab       Number of elements in the array
 ** \param[in]  tab        Array to be printed
 **
 *****************************************************************************/
void print_vector(const char *title,
                  int flag_limit,
                  int ntab,
                  const double *tab)
{
  static int nby_def = 5;

  /* Initializations */

  if (ntab <= 0) return;
  int nby = (flag_limit && (int) OptCst::query(ECst::NTCOL) >= 0) ?
      (int) OptCst::query(ECst::NTCOL) : nby_def;
  bool flag_many = (ntab > nby);

  if (title != NULL)
  {
    message("%s", title);
    if (flag_many) message("\n");
  }
  int lec = 0;
  for (int i = 0; i < ntab; i += nby)
  {
    if (flag_many) message(" %2d+  ", i);
    for (int j = 0; j < nby; j++)
    {
      if (lec >= ntab) continue;
      message(" %10f", tab[lec]);
      lec++;
    }
    message("\n");
  }
}

void print_vector(const char *title,
                  int flag_limit,
                  int ntab,
                  const VectorDouble &tab)
{
  print_vector(title, flag_limit, ntab, tab.data());
}

/****************************************************************************/
/*!
 **  Print a vector of integer values in a matrix form
 **
 ** \param[in]  title      Title (Optional)
 ** \param[in]  flag_limit 1 if NTCOL is used; 0 otherwise
 ** \param[in]  ntab       Number of elements in the array
 ** \param[in]  itab       Array to be printed
 **
 *****************************************************************************/
void print_ivector(const char *title, int flag_limit, int ntab, const int *itab)
{
  static int nby_def = 5;

  /* Initializations */

  if (ntab <= 0) return;
  int nby = (flag_limit && (int) OptCst::query(ECst::NTCOL) >= 0) ?
      (int) OptCst::query(ECst::NTCOL) : nby_def;
  bool flag_many = (ntab > nby);

  if (title != NULL)
  {
    message("%s", title);
    if (flag_many) message("\n");
  }
  int lec = 0;
  for (int i = 0; i < ntab; i += nby)
  {
    if (flag_many) message(" %2d+  ", i);
    for (int j = 0; j < nby; j++)
    {
      if (lec >= ntab) continue;
      message(" %10d", itab[lec]);
      lec++;
    }
    message("\n");
  }
}

void print_ivector(const char *title,
                   int flag_limit,
                   int ntab,
                   const VectorInt &itab)
{
  print_ivector(title, flag_limit, ntab, itab.data());
}
