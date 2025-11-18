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
#include "Basic/VectorNumT.hpp"
#include "Matrix/AMatrix.hpp"

#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <cstring>
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
static String FORMAT;
static String DECODE;
static String TABSTR;

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

/****************************************************************************/
/*!
 **  Tabulated printout of a string

 **
 ** \param[in]  title    optional title (NULL if not defined)
 ** \param[in]  string   String to be written
 ** \param[in]  ncol     number of columns for the printout
 ** \param[in]  justification  justification flag
 **
 *****************************************************************************/
void tab_prints(const String& title,
                const String& string,
                Id ncol,
                Id justification)
{
  Id taille = (1 + static_cast<Id>(OptCst::query(ECst::NTCAR))) * ncol;
  Id size   = string.length();
  Id neff   = MIN(taille, size);
  Id nrst   = taille - neff;
  Id n1     = nrst / 2;
  Id n2     = taille - size - n1;

  /* Encode the title (if defined) */

  if (!title.empty()) message("%s", title.c_str());

  /* Blank the string out */

  (void)gslStrcpy(TABSTR, "");

  /* Switch according to the justification */

  switch (justification)
  {
    case -1:
      (void)gslStrcat(TABSTR, string);
      TABSTR[neff] = '\0';
      for (Id i = 0; i < nrst; i++)
        (void)gslStrcat(TABSTR, " ");
      break;

    case 0:
      for (Id i = 0; i < n1; i++)
        (void)gslStrcat(TABSTR, " ");
      (void)gslStrcat(TABSTR, string);
      TABSTR[n1 + neff] = '\0';
      for (Id i = 0; i < n2; i++)
        (void)gslStrcat(TABSTR, " ");
      break;

    case 1:
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
 ** \param[in]  justification  justification flag
 **
 *****************************************************************************/
void tab_printg(const String& title,
                double value,
                Id ncol,
                Id justification)
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
  tab_prints(title, DECODE, ncol, justification);
}

/****************************************************************************/
/*!
 **  Tabulated printout of a double value
 **
 ** \param[in]  title    optional title (NULL if not defined)
 ** \param[in]  value    Value to be written
 ** \param[in]  ncol     number of columns for the printout
 ** \param[in]  justification  justification flag
 **
 *****************************************************************************/
void tab_printd(const String& title,
                double value,
                Id ncol,
                Id justification)
{
  _buildFormat(CASE_DOUBLE);

  if (FFFF(value))
    (void)gslStrcpy(DECODE, "N/A");
  else
    (void)gslSPrintf(DECODE, FORMAT.data(), value);

  tab_prints(title, DECODE, ncol, justification);
}

/****************************************************************************/
/*!
 **  Tabulated printout of an integer value
 **
 ** \param[in]  title    optional title (NULL if not defined)
 ** \param[in]  value    Value to be written
 ** \param[in]  ncol     number of columns for the printout
 ** \param[in]  justification  justification flag
 **
 *****************************************************************************/
void tab_printi(const String& title, Id value, Id ncol, Id justification)
{
  _buildFormat(CASE_INT);

  if (IFFFF(value))
    (void)gslStrcpy(DECODE, "N/A");
  else
    (void)gslSPrintf(DECODE, FORMAT.data(), value);

  tab_prints(title, DECODE, ncol, justification);
}

/****************************************************************************/
/*!
 **  Tabulated printout of a row or column value
 **
 ** \param[in]  title    optional title (NULL if not defined)
 ** \param[in]  mode     CASE_ROW or CASE_COL
 ** \param[in]  value    Value to be written
 ** \param[in]  ncol     number of columns for the printout
 ** \param[in]  justification  justification flag
 **
 *****************************************************************************/
void tab_print_rc(const String& title,
                  Id mode,
                  Id value,
                  Id ncol,
                  Id justification)
{
  _buildFormat(mode);

  (void)gslSPrintf(DECODE, FORMAT.data(), value);
  string_strip_blanks(DECODE.data(), 0);
  tab_prints(title, DECODE, ncol, justification);
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
 ** \param[in]  tab    array containing the matrix
 **
 ** \remarks The order of the dimension (nx,ny) is opposite
 ** \remarks of the one used in R-packages where dim[1]=nrow and dim[2]=ncol
 **
 *****************************************************************************/
void print_matrix(const String& title,
                  Id flag_limit,
                  Id bycol,
                  Id nx,
                  Id ny,
                  const VectorDouble& tab)
{
  if (tab.empty() || nx <= 0 || ny <= 0) return;
  Id nx_util   = (flag_limit && static_cast<Id>(OptCst::query(ECst::NTCOL)) > 0) ? MIN(static_cast<Id>(OptCst::query(ECst::NTCOL)), nx) : nx;
  Id ny_util   = (flag_limit && static_cast<Id>(OptCst::query(ECst::NTROW)) > 0) ? MIN(static_cast<Id>(OptCst::query(ECst::NTROW)), ny) : ny;
  Id multi_row = (ny > 1 || title.empty());

  /* Print the title (optional) */

  if (!title.empty())
  {
    if (multi_row)
      message("%s\n", title.c_str());
    else
      message("%s ", title.c_str());
  }

  /* Print the header */

  if (multi_row)
  {
    tab_prints(String(), " ");
    for (Id ix = 0; ix < nx_util; ix++)
      tab_print_rc(String(), CASE_COL, ix + 1);
    message("\n");
  }

  /* Print the contents of the array */

  Id ny_done = 0;
  for (Id iy = 0; iy < ny; iy++)
  {
    ny_done++;
    if (ny_done > ny_util) break;
    if (multi_row) tab_print_rc(String(), CASE_ROW, iy + 1);
    for (Id ix = 0; ix < nx_util; ix++)
    {
      Id iad = (bycol) ? iy + ny * ix : ix + nx * iy;
      tab_printg(String(), tab[iad]);
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

void print_matrix(const String& title,
                  Id flag_limit,
                  const AMatrix& mat)
{
  print_matrix(title, flag_limit, true, mat.getNCols(), mat.getNRows(), mat.getValues());
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
void print_trimat(const String& title, Id mode, Id neq, const VectorDouble& tl)
{
#define TRI(i)    (((i) * ((i) + 1)) / 2)
#define TL1(i, j) (tl[(j) * neq + (i) - TRI(j)]) /* only for i >= j */
#define TL2(i, j) (tl[TRI(i) + (j)])             /* only for i >= j */

  /* Initializations */

  if (tl.empty() || neq <= 0) return;

  /* Print the title (optional) */

  if (!title.empty()) message("%s\n", title.c_str());

  /* Print the header */

  tab_prints(String(), " ");
  for (Id ix = 0; ix < neq; ix++)
    tab_print_rc(String(), CASE_COL, ix + 1);
  message("\n");

  /* Print the contents of the array */

  for (Id iy = 0; iy < neq; iy++)
  {
    tab_print_rc(String(), CASE_ROW, iy + 1);
    for (Id ix = 0; ix < neq; ix++)
    {
      if (ix >= iy)
      {
        if (mode == 1)
          tab_printg(String(), TL1(ix, iy));
        else
          tab_printg(String(), TL2(ix, iy));
      }
      else
        tab_prints(String(), " ");
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
 ** \param[in]  tab    array containing the matrix
 **
 *****************************************************************************/
void print_imatrix(const String& title,
                   Id flag_limit,
                   Id bycol,
                   Id nx,
                   Id ny,
                   const VectorInt& tab)
{
  if (tab.empty() || nx <= 0 || ny <= 0) return;
  Id nx_util   = (flag_limit && static_cast<Id>(OptCst::query(ECst::NTCOL)) > 0) ? MIN(static_cast<Id>(OptCst::query(ECst::NTCOL)), nx) : nx;
  Id ny_util   = (flag_limit && static_cast<Id>(OptCst::query(ECst::NTROW)) > 0) ? MIN(static_cast<Id>(OptCst::query(ECst::NTROW)), ny) : ny;
  Id multi_row = (ny > 1 || title.empty());

  /* Print the title (optional) */

  if (!title.empty())
  {
    if (multi_row)
      message("%s\n", title.c_str());
    else
      message("%s ", title.c_str());
  }

  /* Print the header */

  if (multi_row)
  {
    tab_prints(String(), " ");
    for (Id ix = 0; ix < nx_util; ix++)
      tab_print_rc(String(), CASE_COL, ix + 1);
    message("\n");
  }

  /* Print the contents of the array */

  Id ny_done = 0;
  for (Id iy = 0; iy < ny; iy++)
  {
    ny_done++;
    if (ny_done > ny_util) break;
    if (multi_row) tab_print_rc(String(), CASE_ROW, iy + 1);
    for (Id ix = 0; ix < nx_util; ix++)
    {
      Id iad = (bycol) ? iy + ny * ix : ix + nx * iy;
      tab_printi(String(), tab[iad]);
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
 ** \param[in]  tab        Array to be printed
 **
 *****************************************************************************/
void print_vector(const String& title, const VectorDouble& tab)
{
  static Id nby_def = 5;

  /* Initializations */

  Id ntab = tab.size();
  if (ntab <= 0) return;
  Id nby         = nby_def;
  bool flag_many = (ntab > nby);

  if (!title.empty())
  {
    message("%s", title.c_str());
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

void print_vector(const String& title, const VectorInt& tab)
{
  static Id nby_def = 5;

  /* Initializations */

  Id ntab = tab.size();
  if (ntab <= 0) return;
  Id nby         = nby_def;
  bool flag_many = (ntab > nby);

  if (!title.empty())
  {
    message("%s", title.c_str());
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
