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
#include "Basic/ASerializable.hpp"
#include "Basic/AStringable.hpp"
#include <string.h>
#include <stdio.h>
#include <stdarg.h>

static char LINE[LONG_SIZE], LINE_MEM[LONG_SIZE], *LCUR;
static char *cur = NULL;

ASerializable::ASerializable()
  : _fileName()
  , _fileType()
  , _file(nullptr)
  , _currentRecord()
{
}

/****************************************************************************/
/*!
 **   Open an ASCII file
 **
 ** \return  FILE returned pointer
 **
 ** \param[in]  filename Local file name
 ** \param[in]  filetype Type of the file (optional [NULL] when 'w')
 ** \param[in]  mode     "r" or "w"
 **
 *****************************************************************************/
int ASerializable::_fileOpen(const String& filename,
                             const String& filetype,
                             const String& mode)
{
  // Preliminary check
  if (filename.empty())
  {
    messerr("The Neutral File Name cannot be left empty");
    return 1;
  }

  // Check that the file is not already opened (in 'w' mode)
  if (mode == "w")
  {
    if (_file != nullptr)
    {
      messerr("The storage ('w') in Neutral File is not impossible");
      messerr("The previous file has not been closed beforehand");
      return 1;
    }
  }

  // Open the File
  _fileName = filename;
  _fileType = filetype;
  _file = fopen(filename.c_str(), mode.c_str());
  if (_file == (FILE *) NULL)
  {
    messerr("Error when opening the Neutral File %s", filename.c_str());
    return 1;
  }

  // Preliminary action

  if (mode == "r")
  {
    char idtype[LONG_SIZE];
    if (_recordRead("File Type", "%s", idtype))
    {
      _fileClose();
      return 1;
    }
    if (strcmp(idtype,filetype.c_str()))
    {
      messerr(
          "Error: in the File (%s), its Type (%s) does not match the requested one (%s)",
          _fileName.c_str(), idtype, filetype.c_str());
      _fileClose();
      return 1;
    }
  }
  else
  {
    // Write the File ID
    if (!filetype.empty())
    {
      _recordWrite("%s", filetype.c_str());
      _recordWrite("\n");
    }
  }

  return 0;
}

int ASerializable::_fileClose()
{
  if (_file == nullptr)
  {
    messerr("Error: Attempt to close a Neutral File which is not opened");
    return 1;
  }
  fclose(_file);
  _fileName = String();
  _fileType = String();
  _currentRecord = String();
  _file = nullptr;
  return 0;
}

/**
 *
 * @param title
 * @param format
 * @return
 */
int ASerializable::_recordRead(const String& title, const String& format, ...)
{
  va_list ap;
  int error;

  error = 0;
  va_start(ap, format);
  error = _fileRead(format, ap);

  if (error > 0)
  {
    messerr("Error when reading '%s' from %s", title, _fileName);
    messerr("Current Line: %s", _currentRecord.c_str());
  }
  va_end(ap);
  return (error);
}

void ASerializable::_recordWrite(const String& format, ...)
{
  va_list ap;
  va_start(ap, format);
  _fileWrite(format, ap);
  va_end(ap);
}

/****************************************************************************/
/*!
 **  Read the next token from the file
 **
 ** \return  -1 if the end-of-file has been found
 ** \return   1 for a decoding error
 ** \return   0 otherwise
 **
 ** \param[in]  format     format
 ** \param[in]  ap         Value to be read
 **
 *****************************************************************************/
int ASerializable::_fileRead(const String& format, va_list ap)
{
  char DEL_COM = '#';
  char DEL_SEP = ' ';
  char DEL_BLK = ' ';

  const char *fmt;
  int    *ret_i;
  float  *ret_f;
  double *ret_d;
  char   *ret_s;

  /* Loop on the elements to read (from the format) */

  unsigned int ideb = 0;
  while (ideb < format.size())
  {
    /* Eliminate the blanks */

    if (format[ideb] == DEL_BLK)
    {
      ideb++;
      continue;
    }

    label_start: fmt = &format[ideb];
    if (LCUR == NULL)
    {

      /* Read the next line */

      if (fgets(LINE, LONG_SIZE, _file) == NULL) return (-1);
      LINE[strlen(LINE) - 1] = '\0';
      (void) strcpy(LINE_MEM, LINE);

      /* Eliminate the comments */

      int flag_com = 0;
      for (unsigned int i = 0; i < strlen(LINE); i++)
      {
        if (LINE[i] == DEL_COM)
        {
          flag_com = 1 - flag_com;
          LINE[i] = '\0';
        }
        else
        {
          if (flag_com) LINE[i] = '\0';
        }
      }
      cur = LINE;
    }

    /* Decode the line looking for the next token */

    LCUR = strtok(cur, &DEL_SEP);
    cur = NULL;
    if (LCUR == NULL) goto label_start;

    /* Reading */

    if (!strcmp(fmt, "%s"))
    {
      ret_s = va_arg(ap, char *);
      if (!_onlyBlanks(LCUR))
      {
        if (sscanf(LCUR, "%s", ret_s) <= 0) return (1);
      }
      ideb += 2;
    }
    else if (!strcmp(fmt, "%d"))
    {
      ret_i = va_arg(ap, int *);
      if (sscanf(LCUR, "%d", ret_i) <= 0) return (1);
      ideb += 2;
      if (*ret_i == (int) ASCII_TEST) *ret_i = ITEST;
    }
    else if (!strcmp(fmt, "%f"))
    {
      ret_f = va_arg(ap, float *);
      if (sscanf(LCUR, "%f", ret_f) <= 0) return (1);
      ideb += 2;
      if (*ret_f == ASCII_TEST) *ret_f = TEST;
    }
    else if (!strcmp(fmt, "%lf"))
    {
      ret_d = va_arg(ap, double *);
      if (sscanf(LCUR, "%lf", ret_d) <= 0) return (1);
      ideb += 3;
      if (*ret_d == ASCII_TEST) *ret_d = TEST;
    }
    else if (!strcmp(fmt, "%lg"))
    {
      ret_d = va_arg(ap, double *);
      if (sscanf(LCUR, "%lg", ret_d) <= 0) return (1);
      ideb += 3;
      if (*ret_d == ASCII_TEST) *ret_d = TEST;
    }
    else
    {
      messerr("Wrong format %s", fmt);
      va_end(ap);
      return (2);
    }
  }
  return (0);
}

/****************************************************************************/
/*!
 **  Write the next token from the file
 **
 ** \param[in]  format     Encoding format
 ** \param[in]  ap         Value to be written
 **
 *****************************************************************************/
void ASerializable::_fileWrite(const String& format, va_list ap)
{
  double ret_d;
  char *ret_s;
  bool no_blank = false;

  /* Writing */

  if (format == "%s")
  {
    ret_s = va_arg(ap, char *);
    fprintf(_file, "%s", ret_s);
  }
  else if (format == "%d")
  {
    int ret_i = va_arg(ap, int);
    if (ret_i == TEST)
      fprintf(_file, "%5.1lf", ASCII_TEST);
    else
      fprintf(_file, "%d", ret_i);
  }
  else if (format == "%f")
  {
    ret_d = va_arg(ap, double);
    if (ret_d == TEST)
      fprintf(_file, "%5.1lf", ASCII_TEST);
    else
      fprintf(_file, "%f", ret_d);
  }
  else if (format == "%lf")
  {
    ret_d = va_arg(ap, double);
    if (ret_d == TEST)
      fprintf(_file, "%5.1lf", ASCII_TEST);
    else
      fprintf(_file, "%lf", ret_d);
  }
  else if (format == "%lg")
  {
    ret_d = va_arg(ap, double);
    if (ret_d == TEST)
      fprintf(_file, "%5.1lf", ASCII_TEST);
    else
      fprintf(_file, "%lg", ret_d);
  }
  else if (format == "\n")
  {
    fprintf(_file, "\n");
    no_blank = true;
  }
  else if (format == "#")
  {
    ret_s = va_arg(ap, char *);
    fprintf(_file, "# %s\n", ret_s);
    no_blank = true;
  }
  else
  {
    messerr("Wrong format %s", format);
    return;
  }
  if (!no_blank) fprintf(_file, " ");
  return;
}

bool ASerializable::_onlyBlanks(char *string)
{
  int number = strlen(string);
  for (int i = 0; i < number; i++)
  {
    if (string[i] != ' ') return false;
  }
  return true;
}

