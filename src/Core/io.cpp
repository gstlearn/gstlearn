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
#include "geoslib_io.h"
#include "geoslib_old_f.h"

#include "Basic/Utilities.hpp"
#include "Basic/File.hpp"
#include "Basic/String.hpp"
#include "Basic/OptDbg.hpp"
#include "Core/io.hpp"

#include <string.h>
#include <stdarg.h>
#include <math.h>

/*! \cond */
#define OLD 0
#define NEW 1

/*! \endcond */

static char BUFFER[STRING_LENGTH];
static char DEL_COM = '#';
static char DEL_BLK = ' ';
const char* DEL_SEP = " ";

// TODO : No more char* and printf ! Use std::string and iostream
static void st_print(const char *string);
static void st_read(const char* prompt, char* buffer);
static void st_exit(void);
static void (*WRITE_FUNC)(const char*) = (void (*)(const char*)) st_print;
static void (*WARN_FUNC)(const char*) = (void (*)(const char*)) st_print;
static void (*READ_FUNC)(const char*, char*) = st_read;
static void (*EXIT_FUNC)(void) = st_exit;

static char LINE[LONG_SIZE], LINE_MEM[LONG_SIZE], *LCUR, *LINEB;
static char *cur = NULL;

// https://stackoverflow.com/a/26359433/3952924
#ifdef _MSC_VER
#define strncasecmp _strnicmp
#define strcasecmp _stricmp
#endif

/****************************************************************************/
/*!
 **  Exit from the gstlearn library
 **  (not killing the encapsulation if any)
 **
 *****************************************************************************/
static void st_exit(void)
{
  exit(1);
}

/****************************************************************************/
/*!
 **  Internal print from the library
 **
 **  \param[in]  string Message to be printed
 **
 *****************************************************************************/
static void st_print(const char *string)
{
  //(void) printf("%s",string); // Default printf statement
  std::cout << string << std::flush;
}

/****************************************************************************/
/*!
 **  Read a string from the Standard Input
 **
 ** \param[in]  prompt String to be prompted to ask the question
 ** \param[in]  buffer Array where the Input string is stored
 **
 *****************************************************************************/
static void st_read(const char *prompt, char *buffer)
{
  message("%s :", prompt);

  while (fgets(LINE, LONG_SIZE, stdin) == NULL);

  (void) gslStrcpy(buffer, LINE);
  buffer[strlen(buffer) - 1] = '\0';
}

/****************************************************************************/
/*!
 **  Redefine the IO routine for printing message
 **
 ** \param[in]  write_func Writing function
 **
 *****************************************************************************/
void redefine_message(void (*write_func)(const char*))
{
  if (write_func != NULL) WRITE_FUNC = write_func;
}

/****************************************************************************/
/*!
 **  Redefine the IO routine for printing error message
 **
 ** \param[in]  warn_func  Warning function
 **
 *****************************************************************************/
void redefine_error(void (*warn_func)(const char*))
{
  if (warn_func != NULL) WARN_FUNC = warn_func;
}

/****************************************************************************/
/*!
 **  Redefine the IO routine for Reading
 **
 ** \param[in]  read_func  Reading function
 **
 *****************************************************************************/
void redefine_read(void (*read_func)(const char*, char*))
{
  if (read_func != NULL) READ_FUNC = read_func;
}

/****************************************************************************/
/*!
 **  Redefine the exiting routine
 **
 ** \param[in]  exit_func  Exiting function
 **
 *****************************************************************************/
void redefine_exit(void (*exit_func)(void))
{
  if (exit_func != NULL) EXIT_FUNC = exit_func;
}

/*****************************************************************************/
/*!
 **  Strip the blanks from a string
 **
 ** \param[in,out] string  String to be cleaned
 ** \param[in]  flag_lead 1 to strip only the leading blanks
 **
 *****************************************************************************/
void string_strip_blanks(char *string, int flag_lead)

{
  int i, ecr, length, flag_test;

  flag_test = 0;
  length = static_cast<int>(strlen(string));
  for (i = ecr = 0; i < length; i++)
  {
    if (string[i] == ' ' && !flag_test) continue;
    string[ecr++] = string[i];
    if (flag_lead) flag_test = 1;
  }
  string[ecr] = '\0';
}

/*****************************************************************************/
/*!
 **  Strip the leading and trailing quotes from a string
 **
 ** \param[in,out]  string    String to be cleaned
 **
 ** \remarks The quote is searched in first position. If not found, nothing done
 ** \remarks If found, the trailing quote is stripped, if similar to the first
 ** \remarks character
 **
 *****************************************************************************/
void string_strip_quotes(char *string)

{
  int ecr, length;

  length = static_cast<int>(strlen(string));

  if (string[0] != '"') return;
  ecr = 0;
  for (int i = 1; i < length; i++)
  {
    if (string[i] == '"')
    {
      string[ecr] = '\0';
      return;
    }
    string[ecr++] = string[i];
  }
}

#if defined(_WIN32) || defined(_WIN64)
/****************************************************************************/
/*!
 **  Duplicates the strsep function (not available on Windows)
 **  Split the buffer per sentence (delimited by \n)
 **
 ** \return Pointer to the next sentence
 **
 ** \param[in,out]  stringp    Pointer to the buffer to decoded
 ** \param[in]      delim      Delimeter ca
 **
 ** \remark  In output, the buffer he input buffer
 **
 *****************************************************************************/
char * strsep(char **stringp, const char* delim)
{
  char* start = *stringp;
  char* p;

  p = (start != NULL) ? strpbrk(start, delim) : NULL;

  if (p == NULL)
  {
    *stringp = NULL;
  }
  else
  {
    *stringp = p + 1;
  }

  return start;
}
#endif

/****************************************************************************/
/*!
 **  Print a message
 **  This call comes from AStringable where initial message() has been moved
 **
 ** \param[in]  string   String to be displayed
 **
 ****************************************************************************/
void message_extern(const char *string)

{
  WRITE_FUNC(string);
}

/****************************************************************************/
/*!
 **  External function to provoke an exit of API
 **  This call comes from AStringable where initial mes_abort() has been moved
 **
 ****************************************************************************/
void exit_extern()

{
  EXIT_FUNC();
}

/****************************************************************************/
/*!
 **  Problem in memory allocation
 **
 ** \param[in]  nbyte  number of bytes to be allocated
 **
 ****************************************************************************/
void mem_error(int nbyte)

{
  message("Error: Core allocation problem.\n");
  message("       Number of bytes to be allocated = %d\n", nbyte);
}

/****************************************************************************/
/*! 
 **  Open an ASCII file
 **
 ** \return  FILE returned pointer
 **
 ** \param[in]  filename Local file name
 ** \param[in]  mode     type of file (OLD or NEW)
 **
 ** This method is not documented on purpose. It should remain private
 **
 *****************************************************************************/
FILE* _file_open(const char *filename, int mode)
{
  FILE *file;

  /* Dispatch */

  if (mode == OLD)
    file = gslFopen(filename, "r");
  else
    file = gslFopen(filename, "w");

  _erase_current_string();
  return (file);
}

/****************************************************************************/
/*!
 **  Define the file delimitors
 **
 ** \param[in]  del_com  Delimitor for comments
 ** \param[in]  del_sep  Delimitor for separator
 ** \param[in]  del_blk  Delimitor for blank
 **
 ** This method is not documented on purpose. It should remain private
 **
 *****************************************************************************/
void _file_delimitors(char del_com, const char* del_sep, char del_blk)
{
  DEL_COM = del_com;
  DEL_SEP = del_sep;
  DEL_BLK = del_blk;
}

/****************************************************************************/
/*! 
 **  Print the current line read from an ASCII file
 **
 *****************************************************************************/
void print_current_line(void)
{
  messerr("Current Line: %s", LINE_MEM);
}

/****************************************************************************/
/*! 
 **  Check if a string is composed of blanks only
 **
 ** \return  1 if it is only blanks
 **
 ** \param[in]  string     String to be checked
 **
 *****************************************************************************/
static int st_only_blanks(char *string)
{
  int number;

  number = static_cast<int>(strlen(string));
  for (int i = 0; i < number; i++)
  {
    if (string[i] != ' ') return (0);
  }
  return (1);
}

/****************************************************************************/
/*! 
 **  Read the next token from the file
 **
 ** \return  -1 if the end-of-file has been found
 ** \return   1 for a decoding error
 ** \return   0 otherwise
 **
 ** \param[in]  file       FILE structure
 ** \param[in]  format     format
 ** \param[in]  ap         Value to be read
 **
 ** This method is not documented on purpose. It should remain private
 **
 *****************************************************************************/
int _file_read(FILE *file, const char *format, va_list ap)
{
  int flag_com;
  unsigned int ideb, i;
  const char *fmt;
  int *ret_i;
  float *ret_f;
  double *ret_d;
  char *ret_s;

  /* Loop on the elements to read (from the format) */

  ideb = 0;
  while (ideb < strlen(format))
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

      if (fgets(LINE, LONG_SIZE, file) == NULL) return (-1);
      LINE[strlen(LINE) - 1] = '\0';
      (void) gslStrcpy(LINE_MEM, LINE);
      if (OptDbg::query(EDbg::INTERFACE)) message("Lecture ASCII = %s\n", LINE);

      /* Eliminate the comments */

      flag_com = 0;
      for (i = 0; i < strlen(LINE); i++)
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

    LCUR = gslStrtok(cur, DEL_SEP);
    cur = NULL;
    if (LCUR == NULL) goto label_start;
    if (OptDbg::query(EDbg::INTERFACE)) message("String to be decoded = '%s'\n", LCUR);

    /* Reading */

    if (!strcmp(fmt, "%s"))
    {
      ret_s = va_arg(ap, char*);
      if (!st_only_blanks(LCUR))
      {
        if (gslSScanf(LCUR, "%s", ret_s) <= 0) return (1);
      }
      ideb += 2;
      if (OptDbg::query(EDbg::INTERFACE)) message("Decoded String = %s\n", ret_s);
    }
    else if (!strcmp(fmt, "%d"))
    {
      ret_i = va_arg(ap, int*);
      if (gslSScanf(LCUR, "%d", ret_i) <= 0) return (1);
      ideb += 2;
      if (*ret_i == (int) ASCII_TEST) *ret_i = ITEST;
      if (OptDbg::query(EDbg::INTERFACE)) message("Decoded Integer = %i\n", *ret_i);
    }
    else if (!strcmp(fmt, "%f"))
    {
      ret_f = va_arg(ap, float*);
      if (gslSScanf(LCUR, "%f", ret_f) <= 0) return (1);
      ideb += 2;
      if (*ret_f == ASCII_TEST) *ret_f = (float) TEST;
      if (OptDbg::query(EDbg::INTERFACE)) message("Decoded Float = %s\n", *ret_f);
    }
    else if (!strcmp(fmt, "%lf"))
    {
      ret_d = va_arg(ap, double*);
      if (gslSScanf(LCUR, "%lf", ret_d) <= 0) return (1);
      ideb += 3;
      if (*ret_d == ASCII_TEST) *ret_d = TEST;
      if (OptDbg::query(EDbg::INTERFACE)) message("Decoded Double = %lf\n", *ret_d);
    }
    else if (!strcmp(fmt, "%lg"))
    {
      ret_d = va_arg(ap, double*);
      if (gslSScanf(LCUR, "%lg", ret_d) <= 0) return (1);
      ideb += 3;
      if (*ret_d == ASCII_TEST) *ret_d = TEST;
      if (OptDbg::query(EDbg::INTERFACE)) message("Decoded Double = %lg\n", *ret_d);
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
 **  Get the number of tokens in the line
 **
 ** \return   Number of tokens
 **
 ** \param[in]  file       FILE structure
 **
 ** This method is not documented on purpose. It should remain private
 **
 *****************************************************************************/
int _file_get_ncol(FILE *file)

{
  int ncol, flag_com, i;

  /* Initializations */

  ncol = 0;

  /* Read the next line */

  if (fgets(LINE, LONG_SIZE, file) == NULL) return (ncol);
  LINE[strlen(LINE) - 1] = '\0';
  if (OptDbg::query(EDbg::INTERFACE)) message("Lecture ASCII = %s\n", LINE);

  /* Eliminate the comments */

  flag_com = 0;
  for (i = 0; i < (int) strlen(LINE); i++)
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

  /* Get the number of tokens */

  if (gslStrtok(LINE, DEL_SEP) != NULL)
  {
    ncol++;
    while (gslStrtok(NULL, DEL_SEP) != NULL)
      ncol++;
  }

  if (OptDbg::query(EDbg::INTERFACE)) message("Number of columns = %d\n", ncol);
  return (ncol);
}

/****************************************************************************/
/*! 
 **  Erase the current decoding string
 **
 ** This method is not documented on purpose. It should remain private
 **
 *****************************************************************************/
void _erase_current_string(void)
{
  LCUR = NULL;
}

/****************************************************************************/
/*!
 **  Read the next token from the buffer
 **
 ** \return  -1 if the end-of-record has been found
 ** \return   1 for a decoding error
 ** \return   0 otherwise
 **
 ** \param[in]  buffer     Buffer to be read
 ** \param[in]  format     format
 **
 ** \param[out] ap         va_list containing the read variables
 **
 ** This method is not documented on purpose. It should remain private
 **
 *****************************************************************************/
int _buffer_read(char **buffer, const char *format, va_list ap)
{
  int flag_com;
  unsigned int ideb, i;
  const char *fmt;
  int *ret_i;
  float *ret_f;
  double *ret_d;
  char *ret_s;

  /* Loop on the elements to read (from the format) */

  ideb = 0;
  while (ideb < strlen(format))
  {
    /* Eliminate the blanks */

    if (format[ideb] == DEL_BLK)
    {
      ideb++;
      continue;
    }

    /* Loop on the buffer to be decode */

    label_start: fmt = &format[ideb];
    if (LCUR == NULL)
    {

      /* Read the next line */

      LINEB = strsep(buffer, "\n");
      if (LINEB == NULL) return (-1);
      (void) gslStrcpy(LINE_MEM, LINEB);
      if (OptDbg::query(EDbg::INTERFACE)) message("Lecture ASCII = %s\n", LINEB);

      /* Eliminate the comments */

      flag_com = 0;
      for (i = 0; i < strlen(LINEB); i++)
      {
        if (LINEB[i] == DEL_COM)
        {
          flag_com = 1 - flag_com;
          LINEB[i] = '\0';
        }
        else
        {
          if (flag_com) LINEB[i] = '\0';
        }
      }
      cur = LINEB;
    }

    /* Decode the line looking for the next token */

    LCUR = gslStrtok(cur, DEL_SEP);
    cur = NULL;
    if (LCUR == NULL) goto label_start;
    if (OptDbg::query(EDbg::INTERFACE))
      message("String to be decoded = '%s'\n", LCUR);

    /* Reading */

    if (!strcmp(fmt, "%s"))
    {
      ret_s = va_arg(ap, char*);
      if (gslSScanf(LCUR, "%s", ret_s) <= 0) return (1);
      ideb += 2;
      if (OptDbg::query(EDbg::INTERFACE)) message("Decoded String = %s\n", ret_s);
    }
    else if (!strcmp(fmt, "%d"))
    {
      ret_i = va_arg(ap, int*);
      if (gslSScanf(LCUR, "%d", ret_i) <= 0) return (1);
      ideb += 2;
      if (*ret_i == (int) ASCII_TEST) *ret_i = ITEST;
      if (OptDbg::query(EDbg::INTERFACE)) message("Decoded Integer = %i\n", *ret_i);
    }
    else if (!strcmp(fmt, "%f"))
    {
      ret_f = va_arg(ap, float*);
      if (gslSScanf(LCUR, "%f", ret_f) <= 0) return (1);
      ideb += 2;
      if (*ret_f == ASCII_TEST) *ret_f = (float) TEST;
      if (OptDbg::query(EDbg::INTERFACE)) message("Decoded Float = %s\n", *ret_f);
    }
    else if (!strcmp(fmt, "%lf"))
    {
      ret_d = va_arg(ap, double*);
      if (gslSScanf(LCUR, "%lf", ret_d) <= 0) return (1);
      ideb += 3;
      if (*ret_d == ASCII_TEST) *ret_d = TEST;
      if (OptDbg::query(EDbg::INTERFACE)) message("Decoded Double = %lf\n", *ret_d);
    }
    else if (!strcmp(fmt, "%lg"))
    {
      ret_d = va_arg(ap, double*);
      if (gslSScanf(LCUR, "%lg", ret_d) <= 0) return (1);
      ideb += 3;
      if (*ret_d == ASCII_TEST) *ret_d = TEST;
      if (OptDbg::query(EDbg::INTERFACE)) message("Decoded Double = %lg\n", *ret_d);
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
 ** \param[in]  file       FILE structure
 ** \param[in]  format     Encoding format
 ** \param[in]  ap         Value to be written
 **
 ** This method is not documented on purpose. It should remain private
 **
 *****************************************************************************/
void _file_write(FILE *file, const char *format, va_list ap)
{
  int ret_i, no_blank;
  double ret_d;
  char *ret_s;

  /* Initializations */

  no_blank = 0;

  /* Writing */

  if (!strcmp(format, "%s"))
  {
    ret_s = va_arg(ap, char*);
    fprintf(file, "%s", ret_s);
    if (OptDbg::query(EDbg::INTERFACE)) message("Encoded String = %s\n", ret_s);
  }
  else if (!strcmp(format, "%d"))
  {
    ret_i = va_arg(ap, int);
    if (ret_i == TEST)
      fprintf(file, "%5.1lf", ASCII_TEST);
    else
      fprintf(file, "%d", ret_i);
    if (OptDbg::query(EDbg::INTERFACE)) message("Encoded Integer = %i\n", ret_i);
  }
  else if (!strcmp(format, "%f"))
  {
    ret_d = va_arg(ap, double);
    if (ret_d == TEST)
      fprintf(file, "%5.1lf", ASCII_TEST);
    else
      fprintf(file, "%f", ret_d);
    if (OptDbg::query(EDbg::INTERFACE)) message("Encoded Float = %s\n", ret_d);
  }
  else if (!strcmp(format, "%lf"))
  {
    ret_d = va_arg(ap, double);
    if (ret_d == TEST)
      fprintf(file, "%5.1lf", ASCII_TEST);
    else
      fprintf(file, "%lf", ret_d);
    if (OptDbg::query(EDbg::INTERFACE)) message("Encoded Double = %lf\n", ret_d);
  }
  else if (!strcmp(format, "%lg"))
  {
    ret_d = va_arg(ap, double);
    if (ret_d == TEST)
      fprintf(file, "%5.1lf", ASCII_TEST);
    else
      fprintf(file, "%lg", ret_d);
    if (OptDbg::query(EDbg::INTERFACE)) message("Encoded Double = %lg\n", ret_d);
  }
  else if (!strcmp(format, "\n"))
  {
    fprintf(file, "\n");
    no_blank = 1;
  }
  else if (!strcmp(format, "#"))
  {
    ret_s = va_arg(ap, char*);
    fprintf(file, "# %s\n", ret_s);
    no_blank = 1;
    if (OptDbg::query(EDbg::INTERFACE)) message("Encoded Comment = %s\n", ret_s);
  }
  else
  {
    messerr("Wrong format %s", format);
    return;
  }
  if (!no_blank) fprintf(file, " ");
}

/****************************************************************************/
/*! 
 **  Write the next token into the buffer
 **
 ** \param[in]  buffer     Writing buffer
 ** \param[in]  format     Encoding format
 ** \param[in]  ap         va_list to be written
 **
 ** This method is not documented on purpose. It should remain private
 **
 *****************************************************************************/
void _buffer_write(char *buffer, const char *format, va_list ap)
{
  int ret_i, no_blank;
  double ret_d;
  char *ret_s;

  /* Initializations */

  no_blank = 0;

  /* Writing */

  if (!strcmp(format, "%s"))
  {
    ret_s = va_arg(ap, char*);
    (void) gslSPrintf(buffer, "%s", ret_s);
    if (OptDbg::query(EDbg::INTERFACE)) message("Encoded String = %s\n", ret_s);
  }
  else if (!strcmp(format, "%d"))
  {
    ret_i = va_arg(ap, int);
    if (ret_i == TEST)
      (void) gslSPrintf(buffer, "%5.1lf", ASCII_TEST);
    else
      (void) gslSPrintf(buffer, "%d", ret_i);
    if (OptDbg::query(EDbg::INTERFACE)) message("Encoded Integer = %i\n", ret_i);
  }
  else if (!strcmp(format, "%f"))
  {
    ret_d = va_arg(ap, double);
    if (ret_d == TEST)
      (void) gslSPrintf(buffer, "%5.1lf", ASCII_TEST);
    else
      (void) gslSPrintf(buffer, "%f", ret_d);
    if (OptDbg::query(EDbg::INTERFACE)) message("Encoded Float = %s\n", ret_d);
  }
  else if (!strcmp(format, "%lf"))
  {
    ret_d = va_arg(ap, double);
    if (ret_d == TEST)
      (void) gslSPrintf(buffer, "%5.1lf", ASCII_TEST);
    else
      (void) gslSPrintf(buffer, "%lf", ret_d);
    if (OptDbg::query(EDbg::INTERFACE)) message("Encoded Double = %lf\n", ret_d);
  }
  else if (!strcmp(format, "%lg"))
  {
    ret_d = va_arg(ap, double);
    if (ret_d == TEST)
      (void) gslSPrintf(buffer, "%5.1lf", ASCII_TEST);
    else
      (void) gslSPrintf(buffer, "%lg", ret_d);
    if (OptDbg::query(EDbg::INTERFACE)) message("Encoded Double = %lg\n", ret_d);
  }
  else if (!strcmp(format, "\n"))
  {
    (void) gslSPrintf(buffer, "\n");
    no_blank = 1;
  }
  else if (!strcmp(format, "#"))
  {
    ret_s = va_arg(ap, char*);
    (void) gslSPrintf(buffer, "# %s\n", ret_s);
    no_blank = 1;
    if (OptDbg::query(EDbg::INTERFACE)) message("Encoded Comment = %s\n", ret_s);
  }
  else
  {
    messerr("Wrong format %s", format);
    return;
  }
  if (!no_blank) (void) gslStrcat(buffer, " ");
}


/****************************************************************************/
/*! 
 **  Read astring
 **
 ** \param[in]  question  Question to be asked
 ** \param[in]  flag_def  1 if the default is authorized; 0 otherwise
 ** \param[in]  valdef    Default string
 **
 ** \param[out] answer    Answering string
 **
 ** This method is not documented on purpose. It should remain private
 **
 *****************************************************************************/
void _lire_string(const char *question,
                  int flag_def,
                  const char *valdef,
                  char *answer)
{

  loop:

  /* Compose the question */

  (void) gslSPrintf(LINE, "%s ", question);
  if (flag_def) (void) gslSPrintf(&LINE[strlen(LINE)], "(Def=%s) ", valdef);
  (void) gslStrcat(LINE, ": ");

  /* Read the answer */

  READ_FUNC(LINE, BUFFER);

  /* Handle the default value */

  if (strlen(BUFFER) <= 0)
  {
    if (flag_def)
    {
      (void) gslStrcpy(answer, valdef);
    }
    else
    {
      messerr("No default value provided");
      goto loop;
    }
  }
  else
  {
    (void) gslStrcpy(answer, BUFFER);
  }
}

/****************************************************************************/
/*! 
 **  Read an integer value
 **
 ** \return  Integer value
 **
 ** \param[in]  question  Question to be asked
 ** \param[in]  flag_def  1 if the default is authorized; 0 otherwise
 ** \param[in]  valdef    Default value or ITEST
 ** \param[in]  valmin    Minimum authorized value or ITEST
 ** \param[in]  valmax    Maximum authorized value or ITEST
 **
 ** This method is not documented on purpose. It should remain private
 **
 *****************************************************************************/
int _lire_int(const char *question,
              int flag_def,
              int valdef,
              int valmin,
              int valmax)
{
  int rep;

  loop:

  /* Compose the question */

  (void) gslSPrintf(LINE, "%s ", question);
  if (!IFFFF(valmin) && !IFFFF(valmax) && valmin > valmax)
    valmin = valmax = ITEST;
  if (!IFFFF(valmin) && !IFFFF(valdef) && valdef < valmin) valdef = valmin;
  if (!IFFFF(valmax) && !IFFFF(valdef) && valdef > valmax) valdef = valmax;
  if (flag_def && !IFFFF(valdef))
    (void) gslSPrintf(&LINE[strlen(LINE)], "(Def=%d) ", valdef);
  if (IFFFF(valmin))
    (void) gslStrcat(LINE, "[NA,");
  else
    (void) gslSPrintf(&LINE[strlen(LINE)], "[%d,", valmin);
  if (IFFFF(valmax))
    (void) gslStrcat(LINE, "NA] ");
  else
    (void) gslSPrintf(&LINE[strlen(LINE)], "%d] ", valmax);
  (void) gslStrcat(LINE, ": ");

  /* Read the answer */

  READ_FUNC(LINE, BUFFER);

  /* Handle the default value */

  if (strlen(BUFFER) <= 0)
  {
    if (flag_def && !IFFFF(valdef))
    {
      rep = valdef;
    }
    else
    {
      messerr("No default value provided");
      goto loop;
    }
  }
  else
  {
    if (!strcmp(BUFFER, STRING_NA)) return (ITEST);
    rep = atoi(BUFFER);
  }

  /* Check the bounds */

  if (!IFFFF(valmin) && rep < valmin)
  {
    messerr("Answer (%d) must be larger than Minimum (%d)", rep, valmin);
    goto loop;
  }
  if (!IFFFF(valmax) && rep > valmax)
  {
    messerr("Answer (%d) must be smaller than Maximum (%d)", rep, valmax);
    goto loop;
  }
  return (rep);
}

/****************************************************************************/
/*!
 **  Read a double value
 **
 ** \return  Double value
 **
 ** \param[in]  question  Question to be asked
 ** \param[in]  flag_def  1 if the default is authorized; 0 otherwise
 ** \param[in]  valdef    Default value or TEST
 ** \param[in]  valmin    Minimum authorized value or TEST
 ** \param[in]  valmax    Maximum authorized value or TEST
 **
 ** This method is not documented on purpose. It should remain private
 **
 *****************************************************************************/
double _lire_double(const char *question,
                    int flag_def,
                    double valdef,
                    double valmin,
                    double valmax)
{
  double rep;

  loop:

  /* Compose the question */

  (void) gslSPrintf(LINE, "%s ", question);
  if (!FFFF(valmin) && !FFFF(valmax) && valmin > valmax) valmin = valmax = TEST;
  if (!FFFF(valmin) && !FFFF(valdef) && valdef < valmin) valdef = valmin;
  if (!FFFF(valmax) && !FFFF(valdef) && valdef > valmax) valdef = valmax;
  if (flag_def && !FFFF(valdef))
    (void) gslSPrintf(&LINE[strlen(LINE)], "(Def=%lf) ", valdef);
  if (FFFF(valmin))
    (void) gslStrcat(LINE, "[NA,");
  else
    (void) gslSPrintf(&LINE[strlen(LINE)], "[%lf,", valmin);
  if (FFFF(valmax))
    (void) gslStrcat(LINE, "NA] ");
  else
    (void) gslSPrintf(&LINE[strlen(LINE)], "%lf] ", valmax);
  (void) gslStrcat(LINE, ": ");

  /* Read the answer */

  READ_FUNC(LINE, BUFFER);

  /* Handle the default value */

  if (strlen(BUFFER) <= 0)
  {
    if (flag_def)
    {
      rep = valdef;
    }
    else
    {
      messerr("No default value provided");
      goto loop;
    }
  }
  else
  {
    if (!strcmp(BUFFER, STRING_NA)) return (TEST);
    rep = atof(BUFFER);
  }

  /* Check the bounds */

  if (!FFFF(valmin) && rep < valmin)
  {
    messerr("Answer (%lf) must be larger than Minimum (%lf)", rep, valmin);
    goto loop;
  }
  if (!FFFF(valmax) && rep > valmax)
  {
    messerr("Answer (%lf) must be smaller than Maximum (%lf)", rep, valmax);
    goto loop;
  }
  return (rep);
}

/****************************************************************************/
/*!
 **  Read a boolean answer
 **
 ** \return  Integer value: 1 for 'yes' and 0 for 'no'
 **
 ** \param[in]  question  Question to be asked
 ** \param[in]  flag_def  1 if the default is authorized; 0 otherwise
 ** \param[in]  valdef    Default value (0 for NO and 1 for YES)
 **
 ** This method is not documented on purpose. It should remain private
 **
 *****************************************************************************/
int _lire_logical(const char *question, int flag_def, int valdef)
{
  loop:

  /* Compose the question */

  (void) gslSPrintf(LINE, "%s ", question);
  if (flag_def && !IFFFF(valdef))
  {
    if (valdef == 0)
      (void) gslStrcat(LINE, "(Def=n)");
    else
      (void) gslStrcat(LINE, "(Def=y)");
  }
  (void) gslStrcat(LINE, " [y,n] : ");

  /* Read the answer */

  READ_FUNC(LINE, BUFFER);

  /* Handle the default value */

  if (strlen(BUFFER) <= 0)
  {
    if (flag_def && !IFFFF(valdef))
    {
      return (valdef);
    }
    messerr("No default value provided");
    goto loop;
  }
  else
  {

    /* Check the authorized values */

    if (!strcasecmp(BUFFER, "Y")) return (1);
    if (!strcasecmp(BUFFER, "N")) return (0);
    message("The only authorized answers are 'y' or 'n'\n");
    goto loop;
  }
}

/****************************************************************************/
/*! 
 **  Read the next record
 **
 ** This method is not documented on purpose. It should remain private
 **
 *****************************************************************************/
void record_close(void)
{
  cur = NULL;
  LCUR = NULL;
  LINE[0] = '\0';
  LINE_MEM[0] = '\0';
}

/****************************************************************************/
/*! 
 **  Read the next record
 **
 ** \return Error return code
 **
 ** \param[in]  file       Pointer to the file to be read
 ** \param[in]  format     Encoding format
 ** \param[in]  ...        Value to be written
 **
 ** This method is not documented on purpose. It should remain private
 **
 *****************************************************************************/
int _record_read(FILE *file, const char *format, ...)
{
  va_list ap;
  int error;

  error = 0;
  va_start(ap, format);
  error = _file_read(file, format, ap);
  va_end(ap);
  return (error);
}

/****************************************************************************/
/*!
 **  Print the range of values in an array
 **
 ** \param[in]  title    optional title (NULL if not defined)
 ** \param[in]  ntab     number of values
 ** \param[in]  tab      array of values
 ** \param[in]  sel      (optional) selection
 **
 *****************************************************************************/
void print_range(const char *title, int ntab, const double *tab, const double *sel)
{
  if (tab == nullptr || ntab <= 0) return;
  StatResults stats = ut_statistics(ntab, tab, sel);

  /* Encode the title (if defined) */

  if (title != NULL)
    message("%s : ", title);
  else
    message("Range : ");
  message("  ");

  if (FFFF(stats.mini))
    message(STRING_NA);
  else
    message("%lf", stats.mini);
  message(" ; ");
  if (FFFF(stats.maxi))
    message(STRING_NA);
  else
    message("%lf", stats.maxi);
  message(" (%d/%d)\n", stats.nvalid, ntab);
}
