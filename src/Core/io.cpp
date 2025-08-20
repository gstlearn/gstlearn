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
#include "Core/io.hpp"
#include "Basic/File.hpp"
#include "Basic/OptDbg.hpp"
#include "Basic/String.hpp"
#include "Basic/Utilities.hpp"
#include "geoslib_define.h"
#include "geoslib_io.h"
#include <cmath>
#include <cstdarg>
#include <cstring>

/*! \cond */
#define OLD 0
#define NEW 1

/*! \endcond */

namespace gstlrn
{
static char BUFFER[STRING_LENGTH];
static char DEL_COM = '#';
static char DEL_BLK = ' ';
static char DEL_SEP = ' ';

// TODO : No more char* and printf ! Use std::string and iostream
static void st_print(const char* string);
static void st_read(const char* prompt, char* buffer);
static void st_exit(void);
static void (*WRITE_FUNC)(const char*)       = static_cast<void (*)(const char*)>(st_print);
static void (*WARN_FUNC)(const char*)        = static_cast<void (*)(const char*)>(st_print);
static void (*READ_FUNC)(const char*, char*) = st_read;
static void (*EXIT_FUNC)(void)               = st_exit;

static char* LCUR;
static String cur;
static String LINE;
static String LINE_MEM;

// https://stackoverflow.com/a/26359433/3952924
#ifdef _MSC_VER
#  define strncasecmp _strnicmp
#  define strcasecmp  _stricmp
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
static void st_print(const char* string)
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
static void st_read(const char* prompt, char* buffer)
{
  message("%s :", prompt);

  String ligne;
  if (std::getline(std::cin, ligne))
  {
    // ligne contient le texte lu (sans le '\n')
    (void)gslStrcpy(buffer, ligne.data());
    buffer[strlen(buffer) - 1] = '\0';
  }
  else
  {
    ligne.clear(); // rien lu (EOF ou erreur)
  }
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
  if (write_func != nullptr) WRITE_FUNC = write_func;
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
  if (warn_func != nullptr) WARN_FUNC = warn_func;
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
  if (read_func != nullptr) READ_FUNC = read_func;
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
  if (exit_func != nullptr) EXIT_FUNC = exit_func;
}

/*****************************************************************************/
/*!
 **  Strip the blanks from a string
 **
 ** \param[in,out] string  String to be cleaned
 ** \param[in]  flag_lead 1 to strip only the leading blanks
 **
 *****************************************************************************/
void string_strip_blanks(char* string, Id flag_lead)

{
  Id i, ecr, length, flag_test;

  flag_test = 0;
  length    = static_cast<Id>(strlen(string));
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
void string_strip_quotes(char* string)

{
  Id ecr, length;

  length = static_cast<Id>(strlen(string));

  if (string[0] != '"') return;
  ecr = 0;
  for (Id i = 1; i < length; i++)
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
char* strsep(char** stringp, const char* delim)
{
  char* start = *stringp;
  char* p;

  p = (start != nullptr) ? strpbrk(start, delim) : NULL;

  if (p == nullptr)
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
void message_extern(const char* string)

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
void mem_error(Id nbyte)

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
FILE* _file_open(const char* filename, Id mode)
{
  FILE* file;

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
void _token_delimitors(const char del_com, const char del_sep, const char del_blk)
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
  messerr("Current Line: %s", LINE_MEM.data());
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

static Id _stripToken(std::string& token, const char* format)
{
  if (strcmp(format, "%s") != 0)
  {
    // pour %d/%f/%lf, ignorer token vide
    size_t tStart = token.find_first_not_of(DEL_BLK);
    if (tStart == std::string::npos) return 1;
    size_t tEnd = token.find_last_not_of(DEL_BLK);
    token       = token.substr(tStart, tEnd - tStart + 1);
  }
  else
  {
    // pour %s, garder token même s'il est vide
    size_t tStart = token.find_first_not_of(DEL_BLK);
    size_t tEnd   = token.find_last_not_of(DEL_BLK);
    if (tStart != std::string::npos)
      token = token.substr(tStart, tEnd - tStart + 1);
    else
      token = ""; // token vide
  }
  return 0;
}

static Id _decodeToken(const std::string& token, const char* format, void* out)
{
  if (strcmp(format, "%ld") == 0)
  {
    char* endptr;
    int val = std::strtol(token.c_str(), &endptr, 10);
    if (*endptr != '\0') return 1; // conversion échouée
    *static_cast<Id*>(out) = val;
    return 0;
  }
  if (strcmp(format, "%f") == 0)
  {
    char* endptr;
    float val = std::strtof(token.c_str(), &endptr);
    if (*endptr != '\0') return 1;
    *static_cast<float*>(out) = val;
    return 0;
  }
  if (strcmp(format, "%lf") == 0 || strcmp(format, "%lg") == 0)
  {
    char* endptr;
    double val = std::strtod(token.c_str(), &endptr);
    if (*endptr != '\0') return 1;
    *static_cast<double*>(out) = val;
    return 0;
  }
  if (strcmp(format, "%s") == 0)
  {
    std::string& val = *static_cast<std::string*>(out);
    val              = token;
    return 0;
  }
  return 1;
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
Id _buffer_read(const String& line, const char* format, void* out)
{
  static std::string currentLine;
  static size_t pos = 0;

  // initialisation au premier appel
  if (currentLine.empty())
  {
    currentLine = line;
    pos         = 0;

    // supprimer les commentaires
    size_t cmt = currentLine.find(DEL_COM);
    if (cmt != std::string::npos)
      currentLine.erase(cmt);

    // supprimer CR/LF fin de ligne
    while (!currentLine.empty() &&
           (currentLine.back() == '\n' || currentLine.back() == '\r'))
      currentLine.pop_back();
  }

  if (pos >= currentLine.size()) return 1; // plus de tokens

  // extraire le prochain token jusqu'au séparateur
  size_t sepPos = currentLine.find(DEL_SEP, pos);
  std::string token;
  if (sepPos != std::string::npos)
  {
    token = currentLine.substr(pos, sepPos - pos);
    pos   = sepPos + 1;
  }
  else
  {
    token = currentLine.substr(pos);
    pos   = currentLine.size();
  }

  // supprimer blancs autour du token
  if (_stripToken(token, format)) return 1;

  // conversion selon format
  return _decodeToken(token, format, out);
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
void _file_write(FILE* file, const char* format, va_list ap)
{
  Id ret_i, no_blank;
  double ret_d;
  char* ret_s;

  /* Initializations */

  no_blank = 0;

  /* Writing */

  if (!strcmp(format, "%s"))
  {
    ret_s = va_arg(ap, char*);
    fprintf(file, "%s", ret_s);
    if (OptDbg::query(EDbg::INTERFACE)) message("Encoded String = %s\n", ret_s);
  }
  else if (!strcmp(format, "%ld"))
  {
    ret_i = va_arg(ap, Id);
    if (ret_i == TEST)
      fprintf(file, "%5.1lf", ASCII_TEST);
    else
      fprintf(file, "%ld", ret_i);
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
void _buffer_write(String& buffer, const char* format, va_list ap)
{
  Id ret_i, no_blank;
  double ret_d;
  char* ret_s;

  /* Initializations */

  no_blank = 0;

  /* Writing */

  if (!strcmp(format, "%s"))
  {
    ret_s = va_arg(ap, char*);
    (void)gslSPrintf2(buffer, "%s", ret_s);
    if (OptDbg::query(EDbg::INTERFACE)) message("Encoded String = %s\n", ret_s);
  }
  else if (!strcmp(format, "%ld"))
  {
    ret_i = va_arg(ap, Id);
    if (ret_i == TEST)
      (void)gslSPrintf2(buffer, "%5.1lf", ASCII_TEST);
    else
      (void)gslSPrintf2(buffer, "%ld", ret_i);
    if (OptDbg::query(EDbg::INTERFACE)) message("Encoded Integer = %i\n", ret_i);
  }
  else if (!strcmp(format, "%f"))
  {
    ret_d = va_arg(ap, double);
    if (ret_d == TEST)
      (void)gslSPrintf2(buffer, "%5.1lf", ASCII_TEST);
    else
      (void)gslSPrintf2(buffer, "%f", ret_d);
    if (OptDbg::query(EDbg::INTERFACE)) message("Encoded Float = %s\n", ret_d);
  }
  else if (!strcmp(format, "%lf"))
  {
    ret_d = va_arg(ap, double);
    if (ret_d == TEST)
      (void)gslSPrintf2(buffer, "%5.1lf", ASCII_TEST);
    else
      (void)gslSPrintf2(buffer, "%lf", ret_d);
    if (OptDbg::query(EDbg::INTERFACE)) message("Encoded Double = %lf\n", ret_d);
  }
  else if (!strcmp(format, "%lg"))
  {
    ret_d = va_arg(ap, double);
    if (ret_d == TEST)
      (void)gslSPrintf2(buffer, "%5.1lf", ASCII_TEST);
    else
      (void)gslSPrintf2(buffer, "%lg", ret_d);
    if (OptDbg::query(EDbg::INTERFACE)) message("Encoded Double = %lg\n", ret_d);
  }
  else if (!strcmp(format, "\n"))
  {
    (void)gslSPrintf2(buffer, "\n");
    no_blank = 1;
  }
  else if (!strcmp(format, "#"))
  {
    ret_s = va_arg(ap, char*);
    (void)gslSPrintf2(buffer, "# %s\n", ret_s);
    no_blank = 1;
    if (OptDbg::query(EDbg::INTERFACE)) message("Encoded Comment = %s\n", ret_s);
  }
  else
  {
    messerr("Wrong format %s", format);
    return;
  }
  if (!no_blank) (void)gslStrcat2(buffer, " ");
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
void _lire_string(const char* question,
                  Id flag_def,
                  const char* valdef,
                  char* answer)
{
  LINE.resize(LONG_SIZE);
loop:

  /* Compose the question */

  (void)gslSPrintf2(LINE, "%s ", question);
  if (flag_def) (void)gslSPrintfCat2(LINE, "(Def=%s) ", valdef);
  (void)gslStrcat2(LINE, ": ");

  /* Read the answer */

  READ_FUNC(LINE.data(), BUFFER);

  /* Handle the default value */

  if (strlen(BUFFER) <= 0)
  {
    if (flag_def)
    {
      (void)gslStrcpy(answer, valdef);
    }
    else
    {
      messerr("No default value provided");
      goto loop;
    }
  }
  else
  {
    (void)gslStrcpy(answer, BUFFER);
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
Id _lire_int(const char* question,
             Id flag_def,
             Id valdef,
             Id valmin,
             Id valmax)
{
  Id rep;
  LINE.resize(LONG_SIZE);
loop:

  /* Compose the question */

  (void)gslSPrintf2(LINE, "%s ", question);
  if (!IFFFF(valmin) && !IFFFF(valmax) && valmin > valmax)
    valmin = valmax = ITEST;
  if (!IFFFF(valmin) && !IFFFF(valdef) && valdef < valmin) valdef = valmin;
  if (!IFFFF(valmax) && !IFFFF(valdef) && valdef > valmax) valdef = valmax;
  if (flag_def && !IFFFF(valdef))
    (void)gslSPrintfCat2(LINE, "(Def=%d) ", valdef);
  if (IFFFF(valmin))
    (void)gslStrcat2(LINE, "[NA,");
  else
    (void)gslSPrintfCat2(LINE, "[%d,", valmin);
  if (IFFFF(valmax))
    (void)gslStrcat2(LINE, "NA] ");
  else
    (void)gslSPrintfCat2(LINE, "%d] ", valmax);
  (void)gslStrcat2(LINE, ": ");

  /* Read the answer */

  READ_FUNC(LINE.data(), BUFFER);

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
double _lire_double(const char* question,
                    Id flag_def,
                    double valdef,
                    double valmin,
                    double valmax)
{
  double rep;
  LINE.resize(LONG_SIZE);
loop:

  /* Compose the question */

  (void)gslSPrintf2(LINE, "%s ", question);
  if (!FFFF(valmin) && !FFFF(valmax) && valmin > valmax) valmin = valmax = TEST;
  if (!FFFF(valmin) && !FFFF(valdef) && valdef < valmin) valdef = valmin;
  if (!FFFF(valmax) && !FFFF(valdef) && valdef > valmax) valdef = valmax;
  if (flag_def && !FFFF(valdef))
    (void)gslSPrintfCat2(LINE, "(Def=%lf) ", valdef);
  if (FFFF(valmin))
    (void)gslStrcat2(LINE, "[NA,");
  else
    (void)gslSPrintfCat2(LINE, "[%lf,", valmin);
  if (FFFF(valmax))
    (void)gslStrcat2(LINE, "NA] ");
  else
    (void)gslSPrintfCat2(LINE, "%lf] ", valmax);
  (void)gslStrcat2(LINE, ": ");

  /* Read the answer */

  READ_FUNC(LINE.data(), BUFFER);

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
Id _lire_logical(const char* question, Id flag_def, Id valdef)
{
  LINE.resize(LONG_SIZE);
loop:

  /* Compose the question */

  (void)gslSPrintf2(LINE, "%s ", question);
  if (flag_def && !IFFFF(valdef))
  {
    if (valdef == 0)
      (void)gslStrcat2(LINE, "(Def=n)");
    else
      (void)gslStrcat2(LINE, "(Def=y)");
  }
  (void)gslStrcat2(LINE, " [y,n] : ");

  /* Read the answer */

  READ_FUNC(LINE.data(), BUFFER);

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
  cur.clear();
  LINE.clear();
  LINE_MEM.clear();
  LCUR = NULL;
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
Id _record_read(FILE* file, const char* format, void* out)
{
  static std::string currentLine;
  static size_t pos = 0; // position du prochain token dans currentLine

  while (true)
  {
    // lire une nouvelle ligne si nécessaire
    if (currentLine.empty() || pos >= currentLine.size())
    {
      char temp[1024];
      if (!fgets(temp, sizeof(temp), file))
        return 1; // EOF ou erreur

      if (OptDbg::query(EDbg::INTERFACE)) message("Lecture ASCII = %s", temp);
      currentLine = temp;
      pos         = 0;

      // supprimer les commentaires
      size_t cmt = currentLine.find(DEL_COM);
      if (cmt != std::string::npos)
        currentLine.erase(cmt);

      // supprimer CR/LF en fin de ligne
      while (!currentLine.empty() &&
             (currentLine.back() == '\n' || currentLine.back() == '\r'))
        currentLine.pop_back();

      // supprimer blancs en début et fin
      size_t start = currentLine.find_first_not_of(DEL_BLK);
      size_t end   = currentLine.find_last_not_of(DEL_BLK);
      if (start == std::string::npos)
      {
        currentLine.clear(); // ligne vide
        continue;
      }
      currentLine = currentLine.substr(start, end - start + 1);
    }

    // extraire le prochain token jusqu'au séparateur
    size_t sepPos = currentLine.find(DEL_SEP, pos);
    std::string token;
    if (sepPos != std::string::npos)
    {
      token = currentLine.substr(pos, sepPos - pos);
      pos   = sepPos + 1;
    }
    else
    {
      token = currentLine.substr(pos);
      pos   = currentLine.size();
    }

    // Check empty line
    if (_stripToken(token, format)) continue;

    // conversion selon format
    return _decodeToken(token, format, out);
  }
}
} // namespace gstlrn
