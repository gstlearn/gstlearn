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
#include "Basic/Message.hpp"

#include "Basic/AStringable.hpp"
#include "Basic/OptCst.hpp"
#include "Basic/String.hpp"
#include "Basic/VectorNumT.hpp"
#include "Enum/ECst.hpp"
#include "Matrix/AMatrix.hpp"

#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <cstring>
#include <sstream>

namespace gstlrn
{
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
 * @brief Check if 'current' is a valid argument
 * @brief - it must be non negative
 * @brief - it must be lower than 'nmax' (when 'nmax' is specified, i.e. positive)
 *
 * @param title Description of the argument
 * @param current Current value of the argument
 * @param nmax Maximum (exclusive) possible value (optional)
 * @return true If the argument is valid
 * @return false Otherwise
 */
bool checkArg(const char* title, Id current, Id nmax)
{
  if (nmax >= 0)
  {
    if (current < 0 || current >= nmax)
    {
      messerr("Error in %s. The index (%d) should be non-negative and lower than %d",
              title, current, nmax);
      return false;
    }
  }
  else
  {
    if (current < 0)
    {
      messerr("Error in %s. The index (%d) should be non-negative",
              title, current);
      return false;
    }
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
void printElement(const String& string,
                  const String& title,
                  Id ncol,
                  Id justification)
{
  if (!title.empty()) message("%s", title.c_str());
  String str = toStr(string, justification, 0, true, ncol);
  message(str.c_str());
}

/****************************************************************************/
/*!
 **  Tabulated printout of a real value
 **
 ** \param[in]  title    optional title (NULL if not defined)
 ** \param[in]  value    Value to be written
 ** \param[in]  ncol     number of columns for the printout
 ** \param[in]  justification  justification flag
 ** \param[in]  roundZero true to round very small double values to zero
 ** \param[in]  flagScientific true to use scientific notation
 **
 *****************************************************************************/
void printElement(double value,
                  const String& title,
                  Id ncol,
                  Id justification,
                  bool roundZero,
                  bool flagScientific)
{
  if (!title.empty()) message("%s", title.c_str());
  String string = toStr(value, justification, 0, roundZero, ncol, flagScientific);
  message(string.c_str());
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
void printElement(Id value,
                  const String& title,
                  Id ncol,
                  Id justification)
{
  if (!title.empty()) message("%s", title.c_str());
  String string = toStr(value, justification, 0, true, ncol);
  message(string.c_str());
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
void printMatrix(const VectorDouble& tab,
                 Id nx,
                 Id ny,
                 const String& title,
                 Id flag_limit,
                 Id bycol)
{
  String string = toStrMatrix(title, VectorString(), VectorString(), bycol, ny, nx, tab,
                              flag_limit == 0, false);
  message(string.c_str());
}

void printMatrix(const VectorInt& tab,
                 Id nx,
                 Id ny,
                 const String& title,
                 Id flag_limit,
                 Id bycol)
{
  String string = toStrMatrix(title, VectorString(), VectorString(), bycol, ny, nx, tab,
                              flag_limit == 0, false);
  message(string.c_str());
}

void printMatrix(const AMatrix& mat,
                 const String& title,
                 Id flag_limit)
{
  printMatrix(mat.getValues(), mat.getNCols(), mat.getNRows(), title, flag_limit, true);
}

/****************************************************************************/
/*!
 **  Print a vector of real values in a matrix form
 **
 ** \param[in]  title      Title (Optional)
 ** \param[in]  tab        Array to be printed
 ** \param[in]  flagIgnoreMaxNCols  true to ignore the maximum number of columns
 ** \param[in]  newLineAfterTitle   true to put a new line after the title (if any)
 **
 *****************************************************************************/
void printVector(const VectorDouble& tab,
                 const String& title,
                 bool flagIgnoreMaxNCols,
                 bool newLineAfterTitle)
{
  String string = toStrVector(title, tab, flagIgnoreMaxNCols, newLineAfterTitle);
  message(string.c_str());
}

void printVector(const VectorInt& tab,
                 const String& title,
                 bool flagIgnoreMaxNCols,
                 bool newLineAfterTitle)
{
  String string = toStrVector(title, tab, flagIgnoreMaxNCols, newLineAfterTitle);
  message(string.c_str());
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
void printTriMat(const VectorDouble& tl, Id neq, Id mode, const String& title)
{
  String string = toStrTrimMat(tl, neq, mode, title);
  message(string.c_str());
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