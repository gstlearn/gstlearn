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
#include "Basic/AStringable.hpp"
#include "Basic/OptCst.hpp"
#include "Basic/String.hpp"
#include "Basic/Utilities.hpp"
#include "Basic/VectorNumT.hpp"

#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <typeinfo>

#define CASE_DOUBLE 0
#define CASE_REAL   1
#define CASE_INT    2
#define CASE_COL    3
#define CASE_ROW    4

namespace gstlrn
{
class EJustify;

static String FORMAT;
static String DECODE;
static String TABSTR;

static Id _getColumnRank()
{
  return static_cast<Id>(OptCst::query(ECst::NTRANK));
}
static Id _getColumnName()
{
  return static_cast<Id>(OptCst::query(ECst::NTNAME));
}
static Id _getColumnSize()
{
  return static_cast<Id>(OptCst::query(ECst::NTCAR));
}
static Id _getDecimalNumber()
{
  return static_cast<Id>(OptCst::query(ECst::NTDEC));
}
static double _getThresh()
{
  Id ndec       = static_cast<Id>(OptCst::query(ECst::NTDEC));
  double thresh = (0.5 * pow(10, -ndec));
  return thresh;
}
static Id _getMaxNCols()
{
  return static_cast<Id>(OptCst::query(ECst::NTCOL));
}
static Id _getMaxNRows()
{
  return static_cast<Id>(OptCst::query(ECst::NTROW));
}
static Id _getNBatch()
{
  return static_cast<Id>(OptCst::query(ECst::NTBATCH));
}
static void _buildFormat(Id mode)
{
  switch (mode)
  {
    case CASE_INT:
      (void)gslSPrintf(FORMAT, "%%%dd", static_cast<Id>(OptCst::query(ECst::NTCAR)));
      break;

    case CASE_REAL:
      (void)gslSPrintf(FORMAT, "%%%d.%dlf", static_cast<Id>(OptCst::query(ECst::NTCAR)),
                       static_cast<Id>(OptCst::query(ECst::NTDEC)));
      break;

    case CASE_DOUBLE:
      (void)gslSPrintf(FORMAT, "%%%d.%dlg", static_cast<Id>(OptCst::query(ECst::NTCAR)),
                       static_cast<Id>(OptCst::query(ECst::NTDEC)));
      break;

    case CASE_COL:
      (void)gslSPrintf(FORMAT, "[,%%%dd]", static_cast<Id>(OptCst::query(ECst::NTCAR)) - 3);
      break;

    case CASE_ROW:
      (void)gslSPrintf(FORMAT, "[%%%dd,]", static_cast<Id>(OptCst::query(ECst::NTCAR)) - 3);
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
  sstr << "toString is not yet implemented for " << typeid(*this).name() << std::endl;
  return sstr.str();
}

std::stringstream _formatColumn(const EJustify& justify, Id localSize = 0)
{
  std::stringstream sstr;
  auto size = static_cast<I32>((localSize > 0) ? localSize : _getColumnSize());
  auto prec = static_cast<I32>(_getDecimalNumber());
  sstr << std::fixed << std::setw(size) << std::setprecision(prec);
  if (justify == EJustify::LEFT)
    sstr << std::left;
  else
    sstr << std::right;
  return sstr;
}

String _tabPrintString(const String& string,
                       const EJustify& justify,
                       Id localSize = 0)
{
  std::stringstream sstr = _formatColumn(justify, localSize);
  Id size                = static_cast<Id>(string.size());
  Id truncSize           = (localSize > 0) ? localSize : _getColumnSize();
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

String _tabPrintDouble(double value, const EJustify& justify, Id localSize = 0)
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

String _tabPrintInt(Id value, const EJustify& justify, Id localSize = 0)
{
  std::stringstream sstr = _formatColumn(justify, localSize);
  if (IFFFF(value))
    sstr << "N/A";
  else
    sstr << value;

  return sstr.str();
}

String _tabPrintRowColumn(Id icase, Id value, Id flagAdd)
{
  std::stringstream sstr;
  I32 rank  = static_cast<I32>(_getColumnRank());
  I32 width = static_cast<I32>(_getColumnSize() - _getColumnRank() - 1);
  sstr << std::setw(width) << std::right;
  if (icase == CASE_ROW)
  {
    if (!flagAdd)
      sstr << "[" << std::setw(rank) << value << ",]";
    else
      sstr << "[" << std::setw(rank) << value << "+]";
  }
  else
  {
    if (!flagAdd)
      sstr << "[," << std::setw(rank) << value << "]";
    else
      sstr << "[ " << std::setw(rank) << value << "]";
  }
  return sstr.str();
}

String _printColumnHeader(const VectorString& colnames,
                          Id colfrom,
                          Id colto,
                          Id colSize = _getColumnSize())
{
  std::stringstream sstr;
  if (!colnames.empty())
  {
    // By Names
    sstr << _tabPrintString(" ", EJustify::RIGHT) << " ";
    for (Id ix = colfrom; ix < colto; ix++)
      sstr << _tabPrintString(colnames[ix], EJustify::RIGHT, colSize);
    sstr << std::endl;
  }
  else
  {
    // By Numbers
    sstr << _tabPrintString(" ", EJustify::RIGHT) << " ";
    for (Id ix = colfrom; ix < colto; ix++)
      sstr << _tabPrintRowColumn(CASE_COL, ix, false);
    sstr << std::endl;
  }
  return sstr.str();
}

String _printRowHeader(const VectorString& rownames, Id iy, Id rowSize = _getColumnSize())
{
  std::stringstream sstr;
  if (!rownames.empty())
    sstr << _tabPrintString(rownames[iy], EJustify::LEFT, rowSize);
  else
    sstr << _tabPrintRowColumn(CASE_ROW, iy, false);
  return sstr.str();
}

String _printTrailer(Id ncols, Id nrows, Id ncols_util, Id nrows_util)
{
  std::stringstream sstr;

  bool one_used = (ncols != ncols_util || nrows != nrows_util);
  bool all_used = (ncols != ncols_util && nrows != nrows_util);

  if (one_used) sstr << "(";

  if (ncols != ncols_util)
  {
    if (ncols == ncols_util)
      sstr << "Ncols=" << ncols;
    else
      sstr << "Ncols=" << ncols_util << "[from " << ncols << "]";
  }

  if (all_used) sstr << ",";

  if (nrows != nrows_util)
  {
    if (nrows == nrows_util)
      sstr << "Nrows=" << nrows;
    else
      sstr << "Nrows=" << nrows_util << "[from " << nrows << "]";
  }

  if (one_used) sstr << ")" << std::endl;
  return sstr.str();
}

/**
 * Print a formatted message
 * @param format Output format
 * @param ...    Additional arguments
 */
void message(const char* format, ...)
{
  char str[LONG_SIZE];

  va_list ap;
  va_start(ap, format);
  (void)vsnprintf(str, sizeof(str), format, ap);
  va_end(ap);
  message_extern(str);
}

/**
 * Print a formatted message (with "#NO_DIFF#" prefix)
 * @param format Output format
 * @param ...    Additional arguments
 */
void messageNoDiff(const char* format, ...)
{
  char str[LONG_SIZE];
  va_list ap;

  va_start(ap, format);
  (void)vsnprintf(str, sizeof(str), format, ap);
  va_end(ap);
  std::stringstream sstr;
  sstr << "#NO_DIFF# " << str;
  message_extern(sstr.str().c_str());
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
 * Print a standard Error Message if an argument does not lie in Interval
 * @param title   Title to be printed
 * @param current Current value of the argument
 * @param nmax    Maximum (inclusive) possible value
 */
void mesArg(const char* title, Id current, Id nmax)
{
  if (nmax <= 0)
    messerr("Error in %s (%d). No element of this type is recorded yet", title,
            current);
  else
    messerr("Error in %s (%d). Argument should lie within [0,%d[", title,
            current, nmax);
}

bool checkArg(const char* title, Id current, Id nmax)
{
  if (current < 0 || current >= nmax)
  {
    mesArg(title, current, nmax);
    return false;
  }
  return true;
}

/**
 * Print a message and underlines it with various formats
 * @param level  Level of the title
 * @param format Output format
 * @param ...    Additional arguments
 */
void mestitle(Id level, const char* format, ...)
{
  Id STRING_MAX = 1000;
  char STRING[1000];
  va_list ap;

  message_extern("\n");
  va_start(ap, format);
  (void)vsnprintf(STRING, sizeof(STRING), format, ap);
  va_end(ap);
  Id size = static_cast<Id>(strlen(STRING));

  (void)gslStrcat(STRING, STRING_MAX, "\n");
  message_extern(STRING);

  /* Underline the string */

  (void)gslStrcpy(STRING, STRING_MAX, "");
  for (Id i = 0; i < size; i++)
  {
    switch (level)
    {
      case 0:
        (void)gslStrcat(STRING, STRING_MAX, "=");
        break;

      case 1:
        (void)gslStrcat(STRING, STRING_MAX, "-");
        break;

      case 2:
        (void)gslStrcat(STRING, STRING_MAX, ".");
        break;
    }
  }
  (void)gslStrcat(STRING, STRING_MAX, "\n");
  message_extern(STRING);
}

/**
 * Conditionally print the progress of a procedure
 * @param string String to be printed
 * @param ntot   Total number of samples
 * @param iech   Rank of the current sample
 *
 * @remarks The value 'nproc' designates the quantile such that,
 * @remarks when changed, the printout is provoked.
 */
void mes_process(const char* string, Id ntot, Id iech)
{
  static Id quant_memo = 0;
  Id nproc             = static_cast<Id>(OptCst::query(ECst::NPROC));
  if (nproc <= 0) return;
  Id jech = iech + 1;

  /* Calculate the current quantile */

  double ratio = nproc * static_cast<double>(jech) / static_cast<double>(ntot);
  Id quant     = static_cast<Id>(ratio);

  /* Conditional printout */

  if (quant != quant_memo) message("%s - Rank : %d (Quantile : %d / %d)\n", string, iech, quant, nproc);
  quant_memo = quant;
}

/**
 * Print a message and underlines it with various formats
 * @param level  Level of the title
 * @param format Output format
 * @param ...    Additional arguments
 */
String toTitle(Id level, const char* format, ...)
{
  std::stringstream sstr;
  Id STRING_MAX = 1000;
  char STRING[1000];
  va_list ap;

  sstr << std::endl;
  va_start(ap, format);
  (void)vsnprintf(STRING, sizeof(STRING), format, ap);
  va_end(ap);
  sstr << STRING << std::endl;

  /* Underline the string */

  Id size = static_cast<Id>(strlen(STRING));
  (void)gslStrcpy(STRING, STRING_MAX, "");
  for (Id i = 0; i < size; i++)
  {
    switch (level)
    {
      case 0:
        (void)gslStrcat(STRING, STRING_MAX, "=");
        break;

      case 1:
        (void)gslStrcat(STRING, STRING_MAX, "-");
        break;

      case 2:
        (void)gslStrcat(STRING, STRING_MAX, ".");
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
void messageAbort(const char* format, ...)
{
  char STRING[1000];
  va_list ap;

  va_start(ap, format);
  (void)vsnprintf(STRING, sizeof(STRING), format, ap);
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
  if (strfmt != nullptr)
  {
    if (strfmt->hasTitle())
    {
      message_extern(strfmt->getTitle().c_str());
      message_extern("\n");
    }
  }
  message_extern(toString(strfmt).c_str());
}

void AStringable::display(Id level) const
{
  AStringFormat sf(level);
  display(&sf);
}

/**
 * @overload
 * Print the contents of a VectorDouble in a Matrix Form
 * @param title        Title of the printout
 * @param mat          Contents of a AMatrix
 * @param flagOverride true to override printout limitations
 * @param flagSkipZero when true, skip the zero values (represented by a '.' as for sparse matrix)
 *                     always true for sparse matrix
 */
String toMatrix(const String& title,
                const AMatrix& mat,
                bool flagOverride,
                bool flagSkipZero)
{
  flagSkipZero = mat.isSparse();
  return toMatrix(title, VectorString(), VectorString(), true,
                  mat.getNRows(), mat.getNCols(), mat.getValues(),
                  flagOverride, flagSkipZero);
}

/**
 * Print the contents of a VectorDouble in a Matrix Form
 * @fn String gstlrn::toMatrix(const String& title, const VectorString& colnames, const VectorString& rownames, bool bycol, Id nrows, Id ncols, const VectorDouble& tab, bool flagOverride, bool flagSkipZero)
 * @param title        Title of the printout
 * @param colnames     Names of the columns (optional)
 * @param rownames     Names of the rows (optional)
 * @param bycol        True if values as sorted by column; false otherwise
 * @param nrows        Number of rows
 * @param ncols        Number of columns
 * @param tab          VectorDouble containing the values
 * @param flagOverride True to override printout limitations
 * @param flagSkipZero when True, skip the zero values (represented by a '.' as for sparse matrix)
 */
String toMatrix(const String& title,
                const VectorString& colnames,
                const VectorString& rownames,
                bool bycol,
                Id nrows,
                Id ncols,
                const VectorDouble& tab,
                bool flagOverride,
                bool flagSkipZero)
{
  std::stringstream sstr;
  if (tab.empty() || ncols <= 0 || nrows <= 0) return sstr.str();

  return toMatrix(title, colnames, rownames, bycol, nrows, ncols, tab.data(),
                  flagOverride, flagSkipZero);
}

/** * Print the contents of a VectorDouble in a Matrix Form
 * @overload
 * @param title        Title of the printout
 * @param colnames     Names of the columns (optional)
 * @param rownames     Names of the rows (optional)
 * @param bycol        true if values as sorted by column; false otherwise
 * @param nrows        Number of rows
 * @param ncols        Number of columns
 * @param tab          VectorDouble containing the values
 * @param flagOverride true to override printout limitations
 * @param flagSkipZero when true, skip the zero values (represented by a '.' as for sparse matrix)
 */
String toMatrix(const String& title,
                const VectorString& colnames,
                const VectorString& rownames,
                bool bycol,
                Id nrows,
                Id ncols,
                const double* tab,
                bool flagOverride,
                bool flagSkipZero)
{
  std::stringstream sstr;

  /* Initializations */

  Id ncutil = ncols;
  Id nrutil = nrows;
  if (_getMaxNCols() > 0 && ncutil > _getMaxNCols() && !flagOverride) ncutil = _getMaxNCols();
  if (_getMaxNRows() > 0 && nrutil > _getMaxNRows() && !flagOverride) nrutil = _getMaxNRows();
  Id npass       = static_cast<Id>(ceil(static_cast<double>(ncutil) / static_cast<double>(_getNBatch())));
  bool multi_row = nrutil > 1 || npass > 1;

  Id colSize = 0;
  if (colnames.empty())
    colSize = _getColumnSize();
  else
  {
    colSize = MIN(_getColumnName(), getMaxStringSize(colnames) + 1);
    colSize = MAX(colSize, _getColumnSize());
  }
  Id rowSize = 0;
  if (rownames.empty())
    rowSize = _getColumnSize();
  else
    rowSize = MAX(getMaxStringSize(rownames) + 1, _getColumnSize());

  /* Print the title (optional) */

  if (!title.empty())
  {
    sstr << title;
    if (multi_row) sstr << std::endl;
  }

  // Loop on the batches

  for (Id ipass = 0; ipass < npass; ipass++)
  {
    Id jdeb = ipass * _getNBatch();
    Id jfin = MIN(jdeb + _getNBatch(), ncutil);

    /* Print the names of the columns and the column numbers */

    if (multi_row)
      sstr << _printColumnHeader(colnames, jdeb, jfin, colSize);

    /* Loop on the rows */

    for (Id iy = 0; iy < nrutil; iy++)
    {
      if (multi_row) sstr << _printRowHeader(rownames, iy, rowSize);

      /* Loop on the columns */
      for (Id ix = jdeb; ix < jfin; ix++)
      {
        Id iad = (bycol) ? iy + nrows * ix : ix + ncols * iy;
        if (flagSkipZero && ABS(tab[iad]) < EPSILON20)
          sstr << _tabPrintString(".", EJustify::RIGHT, colSize);
        else
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
 * @overload
 * Print the contents of a VectorDouble in a Matrix Form
 * @param title        Title of the printout
 * @param colnames     Names of the columns (optional)
 * @param rownames     Names of the rows (optional)
 * @param bycol        true if values as sorted by column; false otherwise
 * @param nrows        Number of rows
 * @param ncols        Number of columns
 * @param tab          VectorInt containing the values
 * @param flagOverride true to override printout limitations
 * @param flagSkipZero when true, skip the zero values (represented by a '.' as for sparse matrix)
 *
 */
String toMatrix(const String& title,
                const VectorString& colnames,
                const VectorString& rownames,
                bool bycol,
                Id nrows,
                Id ncols,
                const VectorInt& tab,
                bool flagOverride,
                bool flagSkipZero)
{
  std::stringstream sstr;
  if (tab.empty() || ncols <= 0 || nrows <= 0) return sstr.str();

  /* Initializations */

  Id ncutil = ncols;
  Id nrutil = nrows;
  if (_getMaxNCols() > 0 && ncutil > _getMaxNCols() && !flagOverride) ncutil = _getMaxNCols();
  if (_getMaxNRows() > 0 && nrutil > _getMaxNRows() && !flagOverride) nrutil = _getMaxNRows();
  Id npass       = static_cast<Id>(ceil(static_cast<double>(ncutil) / static_cast<double>(_getNBatch())));
  bool multi_row = nrutil > 1 || npass > 1;

  Id colSize = 0;
  if (colnames.empty())
    colSize = _getColumnSize();
  else
  {
    colSize = MIN(_getColumnName(), getMaxStringSize(colnames) + 1);
    colSize = MAX(colSize, _getColumnSize());
  }
  Id rowSize = 0;
  if (rownames.empty())
    rowSize = _getColumnSize();
  else
    rowSize = MAX(getMaxStringSize(rownames) + 1, _getColumnSize());

  /* Print the title (optional) */

  if (!title.empty())
  {
    sstr << title;
    if (multi_row) sstr << std::endl;
  }

  // Loop on the batches

  for (Id ipass = 0; ipass < npass; ipass++)
  {
    Id jdeb = ipass * _getNBatch();
    Id jfin = MIN(jdeb + _getNBatch(), ncutil);

    /* Print the names of the columns and the column numbers */

    if (multi_row)
      sstr << _printColumnHeader(colnames, jdeb, jfin, colSize);

    /* Loop on the rows */

    for (Id iy = 0; iy < nrutil; iy++)
    {
      if (multi_row) sstr << _printRowHeader(rownames, iy, rowSize);

      /* Loop on the columns */

      for (Id ix = jdeb; ix < jfin; ix++)
      {
        Id iad = (bycol) ? iy + nrows * ix : ix + ncols * iy;
        if (flagSkipZero && tab[iad] == 0)
          sstr << _tabPrintString(".", EJustify::RIGHT, colSize);
        else
          sstr << _tabPrintInt(tab[iad], EJustify::RIGHT, colSize);
      }
      sstr << std::endl;
    }
  }

  /* Print the trailer */

  sstr << _printTrailer(ncols, nrows, ncutil, nrutil);
  return sstr.str();
}

/**
 * Printout a vector in a formatted manner
 * @fn String gstlrn::toVector(const String& title, const VectorDouble& tab, bool flagOverride)
 * @param title Title of the printout (or empty string)
 * @param tab   Vector (real values) to be printed
 * @param flagOverride true to override printout limitations
 * @return The string (terminated with a newline)
 */
String toVector(const String& title, const VectorDouble& tab, bool flagOverride)
{
  std::stringstream sstr;
  if (tab.empty()) return sstr.str();

  Id ncols  = static_cast<Id>(tab.size());
  Id ncutil = ncols;
  if (_getMaxNCols() > 0 && ncutil > _getMaxNCols() && !flagOverride) ncutil = _getMaxNCols();
  bool multi_row = ncutil > _getNBatch();

  /* Print the title (optional) */

  if (!title.empty())
  {
    sstr << title;
    if (multi_row) sstr << std::endl;
  }

  Id lec = 0;
  if (multi_row) sstr << _printColumnHeader(VectorString(), 0, _getNBatch());

  for (Id i = 0; i < ncutil; i += _getNBatch())
  {
    if (multi_row) sstr << _printRowHeader(VectorString(), i);

    for (Id j = 0; j < _getNBatch(); j++)
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
 * Printout a vector in a formatted manner
 * @fn String gstlrn::toVector(const String& title, constvect tab, bool flagOverride)
 * @param title Title of the printout (or empty string)
 * @param tab   Vector (real values) to be printed
 * @param flagOverride true to override printout limitations
 * @return The string (terminated with a newline)
 */
String toVector(const String& title, constvect tab, bool flagOverride)
{
  std::stringstream sstr;
  if (tab.empty()) return sstr.str();

  Id ncols  = static_cast<Id>(tab.size());
  Id ncutil = ncols;
  if (_getMaxNCols() > 0 && ncutil > _getMaxNCols() && !flagOverride) ncutil = _getMaxNCols();
  bool multi_row = ncutil > _getNBatch();

  /* Print the title (optional) */

  if (!title.empty())
  {
    sstr << title;
    if (multi_row) sstr << std::endl;
  }

  Id lec = 0;
  if (multi_row) sstr << _printColumnHeader(VectorString(), 0, _getNBatch());

  for (Id i = 0; i < ncutil; i += _getNBatch())
  {
    if (multi_row) sstr << _printRowHeader(VectorString(), i);

    for (Id j = 0; j < _getNBatch(); j++)
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
 * @overload
 * @fn String gstlrn::toVector(const String& title, const VectorVectorDouble& tab, bool flagOverride)
 * @param title Title of the printout (or empty string)
 * @param tab   Vector of vectors (real values) to be printed
 * @param flagOverride true to override printout limitations
 * @return The string (terminated with a newline)
 */
String toVector(const String& title, const VectorVectorDouble& tab, bool flagOverride)
{
  std::stringstream sstr;
  if (tab.empty()) return sstr.str();

  if (!title.empty())
    sstr << title << std::endl;

  Id nrows  = static_cast<Id>(tab.size());
  Id nrutil = nrows;
  if (_getMaxNRows() > 0 && nrutil > _getMaxNRows() && !flagOverride) nrutil = _getMaxNRows();

  for (Id i = 0; i < nrutil; i++)
    sstr << toVector(String(), tab[i], flagOverride);

  // Print the trailer
  sstr << _printTrailer(0, nrows, 0, nrutil);

  return sstr.str();
}

/**
 * Printout a list of vectors in a formatted manner
 * @fn String gstlrn::toVector(const String& title, const VectorVectorInt& tab, bool flagOverride)
 * @overload
 * @param title Title of the printout (or empty string)
 * @param tab   Vector of vectors (integer values) to be printed
 * @param flagOverride true to override printout limitations
 * @return The string (terminated with a newline)
 */
String toVector(const String& title, const VectorVectorInt& tab, bool flagOverride)
{
  std::stringstream sstr;
  if (tab.empty()) return sstr.str();

  if (!title.empty()) sstr << title << std::endl;

  Id nrows  = static_cast<Id>(tab.size());
  Id nrutil = nrows;
  if (_getMaxNRows() > 0 && nrutil > _getMaxNRows() && !flagOverride)
    nrutil = _getMaxNRows();

  for (Id i = 0; i < nrutil; i++) sstr << toVector(String(), tab[i], flagOverride);

  // Print the trailer
  sstr << _printTrailer(0, nrows, 0, nrutil);

  return sstr.str();
}

/**
 * Printout a vector in a formatted manner
 * @fn String gstlrn::toVector(const String& title, const VectorString& tab, bool flagOverride)
 * @overload
 * @param title Title of the printout (or empty string)
 * @param tab   Vector (string values) to be printed
 * @param flagOverride true to override printout limitations
 * @return The string (terminated with a newline)
 */
String toVector(const String& title, const VectorString& tab, bool flagOverride)
{
  std::stringstream sstr;
  if (tab.empty()) return sstr.str();

  Id ncols  = static_cast<Id>(tab.size());
  Id ncutil = ncols;
  if (_getMaxNCols() > 0 && ncutil > _getMaxNCols() && !flagOverride) ncutil = _getMaxNCols();
  bool multi_row = ncutil > _getNBatch();

  /* Print the title (optional) */

  if (!title.empty())
  {
    sstr << title;
    if (multi_row) sstr << std::endl;
  }

  Id lec = 0;
  if (multi_row) sstr << _printColumnHeader(VectorString(), 0, _getNBatch());

  for (Id i = 0; i < ncutil; i += _getNBatch())
  {
    if (multi_row) sstr << _printRowHeader(VectorString(), i);

    for (Id j = 0; j < _getNBatch(); j++)
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
 * @fn String gstlrn::toVector(const String& title, const VectorInt& tab, bool flagOverride)
 * @overload
 * @param title Title of the printout (or empty string)
 * @param tab   Vector (integer values) to be printed
 * @param flagOverride true to override printout limitations
 * @return The string (terminated with a newline)
 */
String toVector(const String& title, const VectorInt& tab, bool flagOverride)
{
  std::stringstream sstr;
  if (tab.empty()) return sstr.str();

  Id ncols  = static_cast<Id>(tab.size());
  Id ncutil = ncols;
  if (_getMaxNCols() > 0 && ncutil > _getMaxNCols() && !flagOverride) ncutil = _getMaxNCols();
  bool multi_row = ncutil > _getNBatch();

  /* Print the title (optional) */

  if (!title.empty())
  {
    sstr << title;
    if (multi_row) sstr << std::endl;
  }

  Id lec = 0;
  if (multi_row) sstr << _printColumnHeader(VectorString(), 0, _getNBatch());

  for (Id i = 0; i < ncutil; i += _getNBatch())
  {
    if (multi_row) sstr << _printRowHeader(VectorString(), i);

    for (Id j = 0; j < _getNBatch(); j++)
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

String toStr(const String& string, const EJustify& justify, Id localSize)
{
  std::stringstream sstr;
  sstr << _tabPrintString(string, justify, localSize);
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
  for (Id i = 0; i < static_cast<Id>(values.size()); i++)
    strings.push_back(toDouble(values[i], justify));
  return strings;
}

String toInt(Id value, const EJustify& justify)
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
void tab_prints(const char* title,
                const char* string,
                Id ncol,
                const EJustify& justify)
{
  Id taille = (1 + static_cast<Id>(OptCst::query(ECst::NTCAR))) * ncol;
  Id size   = static_cast<Id>(strlen(string));
  Id neff   = MIN(taille, size);
  Id nrst   = taille - neff;
  Id n1     = nrst / 2;
  Id n2     = taille - size - n1;

  /* Encode the title (if defined) */

  if (title != nullptr) message("%s", title);

  /* Blank the string out */

  (void)gslStrcpy(TABSTR, "");

  /* Switch according to the justification */

  switch (justify.toEnum())
  {
    case EJustify::E_LEFT:
      (void)gslStrcat(TABSTR, string);
      TABSTR[neff] = '\0';
      for (Id i = 0; i < nrst; i++)
        (void)gslStrcat(TABSTR, " ");
      break;

    case EJustify::E_CENTER:
      for (Id i = 0; i < n1; i++)
        (void)gslStrcat(TABSTR, " ");
      (void)gslStrcat(TABSTR, string);
      TABSTR[n1 + neff] = '\0';
      for (Id i = 0; i < n2; i++)
        (void)gslStrcat(TABSTR, " ");
      break;

    case EJustify::E_RIGHT:
      for (Id i = 0; i < nrst; i++)
        (void)gslStrcat(TABSTR, " ");
      (void)gslStrcat(TABSTR, string);
      break;
  }
  message(TABSTR.data());
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
void tab_printg(const char* title,
                double value,
                Id ncol,
                const EJustify& justify)
{
  _buildFormat(CASE_REAL);

  if (FFFF(value))
    (void)gslStrcpy(DECODE, "N/A");
  else
  {
    // Prevent -0.00 : https://stackoverflow.com/a/12536500/3952924
    value = (ABS(value) < _getThresh()) ? 0. : value;
    (void)gslSPrintf(DECODE, FORMAT.data(), value);
  }
  tab_prints(title, DECODE.data(), ncol, justify);
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
void tab_printd(const char* title,
                double value,
                Id ncol,
                const EJustify& justify)
{
  _buildFormat(CASE_DOUBLE);

  if (FFFF(value))
    (void)gslStrcpy(DECODE, "N/A");
  else
    (void)gslSPrintf(DECODE, FORMAT.data(), value);

  tab_prints(title, DECODE.data(), ncol, justify);
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
void tab_printi(const char* title, Id value, Id ncol, const EJustify& justify)
{
  _buildFormat(CASE_INT);

  if (IFFFF(value))
    (void)gslStrcpy(DECODE, "N/A");
  else
    (void)gslSPrintf(DECODE, FORMAT.data(), value);

  tab_prints(title, DECODE.data(), ncol, justify);
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
void tab_print_rc(const char* title,
                  Id mode,
                  Id value,
                  Id ncol,
                  const EJustify& justify)
{
  _buildFormat(mode);

  (void)gslSPrintf(DECODE, FORMAT.data(), value);
  string_strip_blanks(DECODE.data(), 0);
  tab_prints(title, DECODE.data(), ncol, justify);
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
void tab_print_rowname(const char* string, Id taille)
{
  Id size = static_cast<Id>(strlen(string));
  Id neff = MIN(taille, size);
  Id nrst = taille - neff;

  /* Blank the string out */

  (void)gslStrcpy(TABSTR, "");
  (void)gslStrcat(TABSTR, string);
  TABSTR[neff] = '\0';
  for (Id i = 0; i < nrst; i++)
    (void)gslStrcat(TABSTR, " ");
  message(TABSTR.data());
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
void print_matrix(const char* title,
                  Id flag_limit,
                  Id bycol,
                  Id nx,
                  Id ny,
                  const double* sel,
                  const double* tab)
{
  if (tab == nullptr || nx <= 0 || ny <= 0) return;
  Id nx_util   = (flag_limit && static_cast<Id>(OptCst::query(ECst::NTCOL)) > 0) ? MIN(static_cast<Id>(OptCst::query(ECst::NTCOL)), nx) : nx;
  Id ny_util   = (flag_limit && static_cast<Id>(OptCst::query(ECst::NTROW)) > 0) ? MIN(static_cast<Id>(OptCst::query(ECst::NTROW)), ny) : ny;
  Id multi_row = (ny > 1 || title == nullptr);

  /* Print the title (optional) */

  if (title != nullptr)
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
    for (Id ix = 0; ix < nx_util; ix++)
      tab_print_rc(NULL, CASE_COL, ix + 1);
    message("\n");
  }

  /* Print the contents of the array */

  Id ny_done = 0;
  for (Id iy = 0; iy < ny; iy++)
  {
    if (sel != nullptr && !sel[iy]) continue;
    ny_done++;
    if (ny_done > ny_util) break;
    if (multi_row) tab_print_rc(NULL, CASE_ROW, iy + 1);
    for (Id ix = 0; ix < nx_util; ix++)
    {
      Id iad = (bycol) ? iy + ny * ix : ix + nx * iy;
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

void print_matrix(const char* title,
                  Id flag_limit,
                  const AMatrix& mat)
{
  print_matrix(title, flag_limit, true, mat.getNCols(), mat.getNRows(), nullptr, mat.getValues().data());
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
void print_trimat(const char* title, Id mode, Id neq, const double* tl)
{
#define TRI(i)    (((i) * ((i) + 1)) / 2)
#define TL1(i, j) (tl[(j) * neq + (i) - TRI(j)]) /* only for i >= j */
#define TL2(i, j) (tl[TRI(i) + (j)])             /* only for i >= j */

  /* Initializations */

  if (tl == nullptr || neq <= 0) return;

  /* Print the title (optional) */

  if (title != nullptr) message("%s\n", title);

  /* Print the header */

  tab_prints(NULL, " ");
  for (Id ix = 0; ix < neq; ix++)
    tab_print_rc(NULL, CASE_COL, ix + 1);
  message("\n");

  /* Print the contents of the array */

  for (Id iy = 0; iy < neq; iy++)
  {
    tab_print_rc(NULL, CASE_ROW, iy + 1);
    for (Id ix = 0; ix < neq; ix++)
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
void print_imatrix(const char* title,
                   Id flag_limit,
                   Id bycol,
                   Id nx,
                   Id ny,
                   const double* sel,
                   const Id* tab)
{
  if (tab == nullptr || nx <= 0 || ny <= 0) return;
  Id nx_util   = (flag_limit && static_cast<Id>(OptCst::query(ECst::NTCOL)) > 0) ? MIN(static_cast<Id>(OptCst::query(ECst::NTCOL)), nx) : nx;
  Id ny_util   = (flag_limit && static_cast<Id>(OptCst::query(ECst::NTROW)) > 0) ? MIN(static_cast<Id>(OptCst::query(ECst::NTROW)), ny) : ny;
  Id multi_row = (ny > 1 || title == nullptr);

  /* Print the title (optional) */

  if (title != nullptr)
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
    for (Id ix = 0; ix < nx_util; ix++)
      tab_print_rc(NULL, CASE_COL, ix + 1);
    message("\n");
  }

  /* Print the contents of the array */

  Id ny_done = 0;
  for (Id iy = 0; iy < ny; iy++)
  {
    if (sel != nullptr && !sel[iy]) continue;
    ny_done++;
    if (ny_done > ny_util) break;
    if (multi_row) tab_print_rc(NULL, CASE_ROW, iy + 1);
    for (Id ix = 0; ix < nx_util; ix++)
    {
      Id iad = (bycol) ? iy + ny * ix : ix + nx * iy;
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
void print_vector(const char* title,
                  Id flag_limit,
                  Id ntab,
                  const double* tab)
{
  static Id nby_def = 5;

  /* Initializations */

  if (ntab <= 0) return;
  Id nby         = (flag_limit && static_cast<Id>(OptCst::query(ECst::NTCOL)) >= 0) ? static_cast<Id>(OptCst::query(ECst::NTCOL)) : nby_def;
  bool flag_many = (ntab > nby);

  if (title != nullptr)
  {
    message("%s", title);
    if (flag_many) message("\n");
  }
  Id lec = 0;
  for (Id i = 0; i < ntab; i += nby)
  {
    if (flag_many) message(" %2d+  ", i);
    for (Id j = 0; j < nby; j++)
    {
      if (lec >= ntab) continue;
      message(" %10f", tab[lec]);
      lec++;
    }
    message("\n");
  }
}

void print_vector(const char* title,
                  Id flag_limit,
                  Id ntab,
                  const VectorDouble& tab)
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
void print_ivector(const char* title, Id flag_limit, Id ntab, const Id* itab)
{
  static Id nby_def = 5;

  /* Initializations */

  if (ntab <= 0) return;
  Id nby         = (flag_limit && static_cast<Id>(OptCst::query(ECst::NTCOL)) >= 0) ? static_cast<Id>(OptCst::query(ECst::NTCOL)) : nby_def;
  bool flag_many = (ntab > nby);

  if (title != nullptr)
  {
    message("%s", title);
    if (flag_many) message("\n");
  }
  Id lec = 0;
  for (Id i = 0; i < ntab; i += nby)
  {
    if (flag_many) message(" %2d+  ", i);
    for (Id j = 0; j < nby; j++)
    {
      if (lec >= ntab) continue;
      message(" %10d", itab[lec]);
      lec++;
    }
    message("\n");
  }
}

void print_ivector(const char* title,
                   Id flag_limit,
                   Id ntab,
                   const VectorInt& itab)
{
  print_ivector(title, flag_limit, ntab, itab.data());
}

/**
 * Print Error message
 * @param format Output format
 * @param ...    Additional arguments
 */
void messerr(const char* format, ...)
{
  char str[1000];
  va_list ap;

  va_start(ap, format);
  (void)vsnprintf(str, sizeof(str), format, ap);
  va_end(ap);

  message_extern(str);
  message_extern("\n");
}
} // namespace gstlrn
