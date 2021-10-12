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
#include "geoslib_e.h"
#include "Basic/Utilities.hpp"
#include "Basic/EJustify.hpp"

/*! \cond */
#define OLD 0
#define NEW 1

#define CASE_DOUBLE 0
#define CASE_REAL   1
#define CASE_INT    2
#define CASE_COL    3
#define CASE_ROW    4

/*! \endcond */

typedef struct {
  int    mode;                  // 1 for integer; 2 for real
  int    ival;
  double rval;
  char   keyword[20];
  char   comment[STRING_LENGTH];
} Constant;

static Constant CST[CST_NUMBER] = {
  { 1, 10,      0., "NTCAR" , "Number of characters in printout" },
  { 1,  3,      0., "NTDEC" , "Number of decimal digits in printout" },
  { 1,  7,      0., "NTROW" , "Maximum number of rows in table printout" },
  { 1,  7,      0., "NTCOL" , "Maximum number of columns in table printout" },
  { 1,  0,      0., "NPROC" , "Display the Progress Bar"},
  { 1,  0,      0., "LOCMOD", "Option for updating locator of new variable"},
  { 1,  1,      0., "LOCNEW", "When defining new locator, option for old ones"},
  { 2,  0,      0., "RGL"   , "Use 'rgl' for graphic rendition"},
  { 2,  0,      0., "ASP"   , "Defaulted y/x aspect ratio for graphics"},
  { 2,  0, 1.0e-20, "TOLINV", "Tolerance for matrix inversion"},
  { 2,  0, 1.0e-20, "TOLGEN", "Tolerance for matrix generalized inversion"},
  { 2,  0, 2.3e-16, "EPSMAT", "Tolerance value for Matrix calculations"},
  { 2,  0, 1.0e-15, "EPSSVD", "Tolerance value for SVD Matrix calculations"}
};

static char TABSTR[BUFFER_LENGTH];
static char FORMAT[STRING_LENGTH];
static char DECODE[STRING_LENGTH];
static char BUFFER[STRING_LENGTH];
static char DEL_COM = '#';
static char DEL_SEP = ' ';
static char DEL_BLK = ' ';

static void st_print(const char *string);
static void st_read(const char *,char *);
static void st_exit(void);
static void (*WRITE_FUNC)(const char *) =
  (void (*)(const char*)) st_print;
static void (*WARN_FUNC)(const char *)  =
  (void (*)(const char*)) st_print;
static void (*READ_FUNC)(const char *,char *) = 
  st_read;
static void (*EXIT_FUNC)(void)   = 
  st_exit;

static char LINE[LONG_SIZE],LINE_MEM[LONG_SIZE],*LCUR,*LINEB;
static char *cur = NULL;

/****************************************************************************/
/*!
**  Exit from the Geoslib library 
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
  (void) printf("%s",string); // Default printf statement
}

/****************************************************************************/
/*!
**  Read a string from the Standard Input
**
** \param[in]  prompt String to be prompted to ask the question
** \param[in]  buffer Array where the Input string is stored
**
*****************************************************************************/
static void st_read(const char *prompt,
                    char *buffer)

{
  message("%s :",prompt);

  while (fgets(LINE,LONG_SIZE,stdin) == NULL);

  (void) strcpy(buffer,LINE);
  buffer[strlen(buffer)-1] = '\0';
}

/****************************************************************************/
/*!
**  Construct the FORMAT string
**
** \param[in]  mode  CASE_DOUBLE or CASE_REAL or CASE_INT
**
*****************************************************************************/
static void st_format(int mode)
{

  /* Dispatch */

  switch (mode)
  {
    case CASE_INT:
      (void) sprintf(FORMAT,"%%%dd",
                     CST[CST_NTCAR].ival);
      break;

    case CASE_REAL:
      (void) sprintf(FORMAT,"%%%d.%dlf",
                     CST[CST_NTCAR].ival,
                     CST[CST_NTDEC].ival);
      break;

    case CASE_DOUBLE:
      (void) sprintf(FORMAT,"%%%d.%dlg",
                     CST[CST_NTCAR].ival,
                     CST[CST_NTDEC].ival);
      break;

    case CASE_COL:
      (void) sprintf(FORMAT,"[,%%%dd]",
                     CST[CST_NTCAR].ival-3);
      break;

    case CASE_ROW:
      (void) sprintf(FORMAT,"[%%%dd,]",
                     CST[CST_NTCAR].ival-3);
      break;
  }

  return;
}

/****************************************************************************/
/*!
**  Redefine the IO routine for printing message
**
** \param[in]  write_func Writing function
**
*****************************************************************************/
GEOSLIB_API void redefine_message(void (*write_func)(const char *))
{
  if (write_func != NULL)
    WRITE_FUNC = write_func;
  return;
}

/****************************************************************************/
/*!
**  Redefine the IO routine for printing error message
**
** \param[in]  warn_func  Warning function
**
*****************************************************************************/
GEOSLIB_API void redefine_error(void (*warn_func) (const char *))
{
  if (warn_func != NULL)
    WARN_FUNC = warn_func;
  return;
}

/****************************************************************************/
/*!
**  Redefine the IO routine for Reading
**
** \param[in]  read_func  Reading function
**
*****************************************************************************/
GEOSLIB_API void redefine_read(void (*read_func) (const char *,char *))
{
  if (read_func != NULL)
    READ_FUNC  = read_func;
  return;
}

/****************************************************************************/
/*!
**  Redefine the exiting routine
**
** \param[in]  exit_func  Exiting function
**
*****************************************************************************/
GEOSLIB_API void redefine_exit(void (*exit_func) (void))
{
  if (exit_func != NULL)
    EXIT_FUNC = exit_func;
  return;
}

/*****************************************************************************/
/*!
**  Strip the blanks from a string
**
** \param[in]  string    String to be cleaned
** \param[in]  flag_lead 1 to strip only the leading blanks
**
** \param[out] string    Cleaned string
**
*****************************************************************************/
GEOSLIB_API void string_strip_blanks(char *string,
                                     int   flag_lead)

{
  int  i, ecr, length, flag_test;

  flag_test = 0;
  length = static_cast<int> (strlen(string));
  for (i = ecr = 0; i < length; i++)
  {
    if (string[i] == ' ' && ! flag_test) continue;
    string[ecr++] = string[i];
    if (flag_lead) flag_test = 1;
  }
  string[ecr] = '\0';

  return;
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
GEOSLIB_API void string_strip_quotes(char *string)

{
  int ecr,length;

  length = static_cast<int> (strlen(string));

  if (string[0] != '"') return;
  ecr = 0;
  for (int i=1; i < length; i++)
  {
    if (string[i] == '"') 
    {
      string[ecr] = '\0';
      return;
    }
    string[ecr++] = string[i];
  }
  return;
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
GEOSLIB_API char * strsep(char **stringp, const char* delim)
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
GEOSLIB_API void message_extern(const char *string)

{
  WRITE_FUNC(string);
}

/****************************************************************************/
/*!
**  External function to provoke an exit of API
**  This call comes from AStringable where initial mes_abort() has been moved
**
****************************************************************************/
GEOSLIB_API void exit_extern()

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
GEOSLIB_API void mem_error(int nbyte)

{
  message("Error: Core allocation problem.\n");
  message("       Number of bytes to be allocated = %d\n",nbyte);
  return;
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
GEOSLIB_API FILE *_file_open(const char *filename,
                             int   mode)
{
  FILE *file;

  /* Dispatch */

  if (mode == OLD)
    file = fopen(filename,"r");
  else
    file = fopen(filename,"w");

  _erase_current_string();
  return(file);
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
GEOSLIB_API void _file_delimitors(char del_com,
                                  char del_sep,
                                  char del_blk)
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
GEOSLIB_API void print_current_line(void)
{
  messerr("Current Line: %s",LINE_MEM);
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

  number = static_cast<int> (strlen(string));
  for (int i=0; i<number; i++)
  {
    if (string[i] != ' ') return(0);
  }
  return(1);
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
GEOSLIB_API int _file_read(FILE *file,
                           const char *format,
                           va_list ap)
{
  int flag_com;
  unsigned int ideb,i;
  const char   *fmt;
  int    *ret_i;
  float  *ret_f;
  double *ret_d;
  char   *ret_s;

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

  label_start:
    fmt = &format[ideb];
    if (LCUR == NULL)
    {

      /* Read the next line */

      if (fgets(LINE,LONG_SIZE,file) == NULL) return(-1);
      LINE[strlen(LINE)-1] = '\0';
      (void) strcpy(LINE_MEM,LINE);
      if (debug_query("interface"))
        message("Lecture ASCII = %s\n",LINE);

      /* Eliminate the comments */

      flag_com = 0;
      for (i=0; i<strlen(LINE); i++)
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

    LCUR = strtok(cur,&DEL_SEP);
    cur = NULL;
    if (LCUR == NULL) goto label_start;
    if (debug_query("interface"))
      message("String to be decoded = '%s'\n",LCUR);

    /* Reading */

    if (! strcmp(fmt,"%s"))
    {
      ret_s = va_arg(ap,char *);
      if (! st_only_blanks(LCUR))
      {
        if (sscanf(LCUR,"%s",ret_s) <= 0) return(1);
      }
      ideb += 2;
      if (debug_query("interface"))
        message("Decoded String = %s\n",ret_s);
    }
    else if (! strcmp(fmt,"%d"))
    {
      ret_i = va_arg(ap,int *);
      if (sscanf(LCUR,"%d",ret_i) <= 0) return(1);
      ideb += 2;
      if (*ret_i == (int) ASCII_TEST) *ret_i = ITEST;
      if (debug_query("interface"))
        message("Decoded Integer = %i\n",*ret_i);
    }
    else if (! strcmp(fmt,"%f"))
    {
      ret_f = va_arg(ap,float *);
      if (sscanf(LCUR,"%f",ret_f) <= 0) return(1);
      ideb += 2;
      if (*ret_f == ASCII_TEST) *ret_f = (float) TEST;
      if (debug_query("interface"))
        message("Decoded Float = %s\n",*ret_f);
    }
    else if (! strcmp(fmt,"%lf"))
    {
      ret_d = va_arg(ap,double *);
      if (sscanf(LCUR,"%lf",ret_d) <= 0) return(1);
      ideb += 3;
      if (*ret_d == ASCII_TEST) *ret_d = TEST;
      if (debug_query("interface"))
        message("Decoded Double = %lf\n",*ret_d);
    }
    else if (! strcmp(fmt,"%lg"))
    {
      ret_d = va_arg(ap,double *);
      if (sscanf(LCUR,"%lg",ret_d) <= 0) return(1);
      ideb += 3;
      if (*ret_d == ASCII_TEST) *ret_d = TEST;
      if (debug_query("interface"))
        message("Decoded Double = %lg\n",*ret_d);
    }
    else
    {
      messerr("Wrong format %s",fmt);
      va_end(ap);
      return(2);
    }
  }
  return(0);
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
GEOSLIB_API int _file_get_ncol(FILE *file)

{
  int ncol,flag_com,i;

  /* Initializations */

  ncol = 0;

  /* Read the next line */

  if (fgets(LINE,LONG_SIZE,file) == NULL) return(ncol);
  LINE[strlen(LINE)-1] = '\0';
  if (debug_query("interface"))
    message("Lecture ASCII = %s\n",LINE);

  /* Eliminate the comments */
  
  flag_com = 0;
  for (i=0; i<(int) strlen(LINE); i++)
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

  if (strtok(LINE,&DEL_SEP) != NULL)
  {
    ncol++;
    while (strtok(NULL,&DEL_SEP) != NULL) ncol++;
  }

  if (debug_query("interface"))
    message("Number of columns = %d\n",ncol);
  return(ncol);
}

/****************************************************************************/
/*! 
**  Erase the current decoding string
**
** This method is not documented on purpose. It should remain private
**
*****************************************************************************/
GEOSLIB_API void _erase_current_string(void)
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
GEOSLIB_API int _buffer_read(char       **buffer,
                             const char  *format,
                             va_list      ap)
{
  int flag_com;
  unsigned int ideb,i;
  const char *fmt;
  int    *ret_i;
  float  *ret_f;
  double *ret_d;
  char   *ret_s;

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

  label_start:
    fmt = &format[ideb];
    if (LCUR == NULL)
    {

      /* Read the next line */

      LINEB = strsep(buffer,"\n");
      if (LINEB == NULL) return(-1);
      (void) strcpy(LINE_MEM,LINEB);
      if (debug_query("interface"))
        message("Lecture ASCII = %s\n",LINEB);
      
      /* Eliminate the comments */
      
      flag_com = 0;
      for (i=0; i<strlen(LINEB); i++)
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

    LCUR = strtok(cur,&DEL_SEP);
    cur = NULL;
    if (LCUR == NULL) goto label_start;
    if (debug_query("interface"))
      message("String to be decoded = '%s'\n",LCUR);

    /* Reading */

    if (! strcmp(fmt,"%s"))
    {
      ret_s = va_arg(ap,char *);
      if (sscanf(LCUR,"%s",ret_s) <= 0) return(1);
      ideb += 2;
      if (debug_query("interface"))
        message("Decoded String = %s\n",ret_s);
    }
    else if (! strcmp(fmt,"%d"))
    {
      ret_i = va_arg(ap,int *);
      if (sscanf(LCUR,"%d",ret_i) <= 0) return(1);
      ideb += 2;
      if (*ret_i == (int) ASCII_TEST) *ret_i = ITEST;
      if (debug_query("interface"))
        message("Decoded Integer = %i\n",*ret_i);
    }
    else if (! strcmp(fmt,"%f"))
    {
      ret_f = va_arg(ap,float *);
      if (sscanf(LCUR,"%f",ret_f) <= 0) return(1);
      ideb += 2;
      if (*ret_f == ASCII_TEST) *ret_f = (float) TEST;
      if (debug_query("interface"))
        message("Decoded Float = %s\n",*ret_f);
    }
    else if (! strcmp(fmt,"%lf"))
    {
      ret_d = va_arg(ap,double *);
      if (sscanf(LCUR,"%lf",ret_d) <= 0) return(1);
      ideb += 3;
      if (*ret_d == ASCII_TEST) *ret_d = TEST;
      if (debug_query("interface"))
        message("Decoded Double = %lf\n",*ret_d);
    }
    else if (! strcmp(fmt,"%lg"))
    {
      ret_d = va_arg(ap,double *);
      if (sscanf(LCUR,"%lg",ret_d) <= 0) return(1);
      ideb += 3;
      if (*ret_d == ASCII_TEST) *ret_d = TEST;
      if (debug_query("interface"))
        message("Decoded Double = %lg\n",*ret_d);
    }
    else
    {
      messerr("Wrong format %s",fmt);
      va_end(ap);
      return(2);
    }
  }
  return(0);
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
GEOSLIB_API void _file_write(FILE *file,
                             const char *format,
                             va_list ap)
{
  int     ret_i,no_blank;
  double  ret_d;
  char   *ret_s;

  /* Initializations */

  no_blank = 0;

  /* Writing */

  if (! strcmp(format,"%s"))
  {
    ret_s = va_arg(ap,char *);
    fprintf(file,"%s",ret_s);
    if (debug_query("interface"))
      message("Encoded String = %s\n",ret_s);
  }
  else if (! strcmp(format,"%d"))
  {
    ret_i = va_arg(ap,int);
    if (ret_i == TEST) 
      fprintf(file,"%5.1lf",ASCII_TEST);
    else
      fprintf(file,"%d",ret_i);
    if (debug_query("interface"))
      message("Encoded Integer = %i\n",ret_i);
  }
  else if (! strcmp(format,"%f"))
  {
    ret_d = va_arg(ap,double);
    if (ret_d == TEST) 
      fprintf(file,"%5.1lf",ASCII_TEST);
    else
      fprintf(file,"%f",ret_d);
    if (debug_query("interface"))
      message("Encoded Float = %s\n",ret_d);
  }
  else if (! strcmp(format,"%lf"))
  {
    ret_d = va_arg(ap,double);
    if (ret_d == TEST) 
      fprintf(file,"%5.1lf",ASCII_TEST);
    else
      fprintf(file,"%lf",ret_d);
    if (debug_query("interface"))
      message("Encoded Double = %lf\n",ret_d);
  }
  else if (! strcmp(format,"%lg"))
  {
    ret_d = va_arg(ap,double);
    if (ret_d == TEST) 
      fprintf(file,"%5.1lf",ASCII_TEST);
    else
      fprintf(file,"%lg",ret_d);
    if (debug_query("interface"))
      message("Encoded Double = %lg\n",ret_d);
  }
  else if (! strcmp(format,"\n"))
  {
    fprintf(file,"\n");
    no_blank = 1;
  }
  else if (! strcmp(format,"#"))
  {
    ret_s = va_arg(ap,char *);
    fprintf(file,"# %s\n",ret_s);
    no_blank = 1;
    if (debug_query("interface"))
      message("Encoded Comment = %s\n",ret_s);
  }
  else
  {
    messerr("Wrong format %s",format);
    return;
  }
  if (! no_blank) fprintf(file," ");
  return;
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
GEOSLIB_API void _buffer_write(char       *buffer,
                               const char *format,
                               va_list     ap)
{
  int     ret_i,no_blank;
  double  ret_d;
  char   *ret_s;

  /* Initializations */

  no_blank = 0;

  /* Writing */

  if (! strcmp(format,"%s"))
  {
    ret_s = va_arg(ap,char *);
    (void) sprintf(buffer,"%s",ret_s);
    if (debug_query("interface"))
      message("Encoded String = %s\n",ret_s);
  }
  else if (! strcmp(format,"%d"))
  {
    ret_i = va_arg(ap,int);
    if (ret_i == TEST) 
      (void) sprintf(buffer,"%5.1lf",ASCII_TEST);
    else
      (void) sprintf(buffer,"%d",ret_i);
    if (debug_query("interface"))
      message("Encoded Integer = %i\n",ret_i);
  }
  else if (! strcmp(format,"%f"))
  {
    ret_d = va_arg(ap,double);
    if (ret_d == TEST) 
      (void) sprintf(buffer,"%5.1lf",ASCII_TEST);
    else
      (void) sprintf(buffer,"%f",ret_d);
    if (debug_query("interface"))
      message("Encoded Float = %s\n",ret_d);
  }
  else if (! strcmp(format,"%lf"))
  {
    ret_d = va_arg(ap,double);
    if (ret_d == TEST) 
      (void) sprintf(buffer,"%5.1lf",ASCII_TEST);
    else
      (void) sprintf(buffer,"%lf",ret_d);
    if (debug_query("interface"))
      message("Encoded Double = %lf\n",ret_d);
  }
  else if (! strcmp(format,"%lg"))
  {
    ret_d = va_arg(ap,double);
    if (ret_d == TEST) 
      (void) sprintf(buffer,"%5.1lf",ASCII_TEST);
    else
      (void) sprintf(buffer,"%lg",ret_d);
    if (debug_query("interface"))
      message("Encoded Double = %lg\n",ret_d);
  }
  else if (! strcmp(format,"\n"))
  {
    (void) sprintf(buffer,"\n");
    no_blank = 1;
  }
  else if (! strcmp(format,"#"))
  {
    ret_s = va_arg(ap,char *);
    (void) sprintf(buffer,"# %s\n",ret_s);
    no_blank = 1;
    if (debug_query("interface"))
      message("Encoded Comment = %s\n",ret_s);
  }
  else
  {
    messerr("Wrong format %s",format);
    return;
  }
  if (! no_blank) (void) strcat(buffer," ");
  return;
}

/****************************************************************************/
/*!
**  Reset the IO parameters
**
*****************************************************************************/
GEOSLIB_API void constant_reset(void)

{
  CST[CST_NTCAR].ival  =      10;
  CST[CST_NTDEC].ival  =       3;
  CST[CST_NTROW].ival  =       7;
  CST[CST_NTCOL].ival  =       7;
  CST[CST_NPROC].ival  =       0;
  CST[CST_LOCMOD].ival =       1;
  CST[CST_LOCNEW].ival =       0;
  CST[CST_RGL].ival    =       0;
  CST[CST_ASP].ival    =       0;
  CST[CST_TOLINV].rval = matrix_constant_query(CST_TOLINV);
  CST[CST_TOLGEN].rval = matrix_constant_query(CST_TOLGEN);
  CST[CST_EPSMAT].rval = matrix_constant_query(CST_EPSMAT);
  CST[CST_EPSSVD].rval = matrix_constant_query(CST_EPSSVD);

  return;
}

/****************************************************************************/
/*!
**  Print the authorized keywords
**
*****************************************************************************/
static void st_constant_list(void)
{
  int i;

  message("The keywords for IO parameter definition are:\n");
  for (i=0; i<CST_NUMBER; i++)
    message("%6s : %s\n",CST[i].keyword,CST[i].comment);

  return;
}

/****************************************************************************/
/*!
**  Define one IO parameter
**
** \param[in]  name   Name of the parameter to be defined
** \li                NTCAR  : Number of characters in printout
** \li                NTDEC  : Number of decimal digits in printout
** \li                NTROW  : Maximum number of rows in table printout
** \li                NTCOL  : Maximum number of columns in table printout
** \li                NPROC  : Display the Progress Bar
** \li                LOCMOD : Option for updating locator of new variable
** \li                LOCNEW : When defining new locator, option for old ones
** \li                RGL    : Using 'rgl' for graphic rendition
** \li                ASP    : Default y/x aspcet ratio for graphics
** \li                TOLINV : Tolerance for matrix inversion
** \li                TOLGEN : Tolerance for matrix generalized inversion
** \li                EPSMAT : Tolerance value for Matrix calculations
** \li                EPSSVD : Tolerance value for SVD Matrix calculations
** \param[in]  value  New value for the IO parameter
**
*****************************************************************************/
GEOSLIB_API void constant_define(const  char *name,
                                 double value)
{
  int i,found,flag_defined;

  /* Look for an authorized keyword */

  for (i=0, found= -1; i<CST_NUMBER; i++)
    if (! strcasecmp(name,CST[i].keyword)) found = i;

  if (found < 0)
  {
    st_constant_list();
    message("The keyword '%s' is unknown\n",name);
  }
  else
  {
    flag_defined = ! FFFF(value);

    if (CST[found].mode == 1)
    {
      CST[found].ival = (flag_defined) ? (int) value : -1;
      if (found == CST_NTCAR) setFormatColumnSize(static_cast<int> (value));
      if (found == CST_NTDEC) setFormatDecimalNumber(static_cast<int> (value));
      if (found == CST_NTCOL) setFormatMaxNCols(static_cast<int> (value));
      if (found == CST_NTROW) setFormatMaxNRows(static_cast<int> (value));
    }
    else
    {
      CST[found].rval = (flag_defined) ? value : -1;
    }

    if (found == CST_TOLINV ||
        found == CST_TOLGEN ||
        found == CST_EPSMAT ||
        found == CST_EPSSVD)
      matrix_constant_define(found,value);
  }

  return;
}

/****************************************************************************/
/*!
**  Query one IO parameter
**
** \return  Value of the IO parameter
**
** \param[in]  name  Name of the IO parameter to be asked
**
*****************************************************************************/
GEOSLIB_API double constant_query(const char *name)

{
  int    i,found;
  double value;

  /* Look for an authorized keyword */

  for (i=0, found= -1; i<CST_NUMBER; i++)
    if (! strcasecmp(name,CST[i].keyword)) found = i;

  if (found < 0)
  {
    st_constant_list();
    message("The keyword '%s' is unknown\n",name);
    value = 0;
  }
  else
  {
    if (CST[found].mode == 1)
      value = (double) CST[found].ival;
    else
      value = CST[found].rval;
  }

  return(value);
}

/****************************************************************************/
/*!
**  Print the constants for IO
**
*****************************************************************************/
GEOSLIB_API void constant_print(void)

{
  int    i,ival;
  double rval;

  mestitle(1,"Parameters for printout");
  for (i=0; i<CST_NUMBER; i++)
  {
    message (". %-50s [%6s] = ",CST[i].comment,CST[i].keyword);
    if (CST[i].mode == 1)
    {
      ival = CST[i].ival;
      if (ival > 0)
        message ("%d\n",ival);
      else
        message ("NA\n");
    }
    else
    {
      rval = CST[i].rval;
      if (rval > 0)
        message ("%lg\n",rval);
      else
        message ("NA\n");
    }
  }
  message("Use 'constant.define' to modify previous values\n");
}

/****************************************************************************/
/*!
**  Tabulated printout of a string
**
** \param[in]  title    optional title (NULL if not defined)
** \param[in]  ncol     number of columns for the printout
** \param[in]  justify  justification flag
**                      (EJustify::LEFT, EJustify::CENTER or EJustify::RIGHT)
** \param[in]  string   String to be written
**
*****************************************************************************/
GEOSLIB_API void tab_prints(const char*     title,
                            int             ncol,
                            const EJustify& justify,
                            const char*     string)
{
  int i,size,neff,nrst,n1,n2,taille;

  taille = CST[CST_NTCAR].ival * ncol;
  size   = static_cast<int> (strlen(string));
  neff   = MIN(taille, size);
  nrst   = taille - neff;
  n1     = nrst / 2;
  n2     = taille - size - n1;

  /* Encode the title (if defined) */

  if (title != NULL) message("%s",title);

  /* Blank the string out */

  (void) strcpy(TABSTR,"");

  /* Switch according to the justification */

  switch (justify.toEnum())
  {
    case EJustify::E_LEFT:
      (void) strncpy(TABSTR,string,neff);
      TABSTR[neff] = '\0';
      for (i=0; i<nrst; i++) (void) strcat(TABSTR," ");
      break;

    case EJustify::E_CENTER:
      for (i=0; i<n1; i++) (void) strcat(TABSTR," ");
      (void) strncpy(&TABSTR[n1],string,neff);
      TABSTR[n1+neff] = '\0';
      for (i=0; i<n2; i++) (void) strcat(TABSTR," ");
      break;

    case EJustify::E_RIGHT:
      for (i=0; i<nrst; i++) (void) strcat(TABSTR," ");
      (void) strncpy(&TABSTR[nrst],string,neff);
      TABSTR[nrst+neff] = '\0';
      break;
  }
  message(TABSTR);
  return;
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
GEOSLIB_API void tab_print_rowname(const char *string,
                                   int   taille)
{
  int i,size,neff,nrst;

  size   = static_cast<int> (strlen(string));
  neff   = MIN(taille, size);
  nrst   = taille - neff;

  /* Blank the string out */

  (void) strcpy(TABSTR,"");
  (void) strncpy(TABSTR,string,neff);
  TABSTR[neff] = '\0';
  for (i=0; i<nrst; i++) (void) strcat(TABSTR," ");
  message(TABSTR);
  return;
}

/****************************************************************************/
/*!
**  Tabulated printout of a real value
**
** \param[in]  title    optional title (NULL if not defined)
** \param[in]  ncol     number of columns for the printout
** \param[in]  justify  justification flag
**                      (EJustify::LEFT, EJustify::CENTER or EJustify::RIGHT)
** \param[in]  value    Value to be written
**
*****************************************************************************/
GEOSLIB_API void tab_printg(const char*     title,
                            int             ncol,
                            const EJustify& justify,
                            double          value)
{
  st_format(CASE_REAL);

  if (FFFF(value))
    (void) strcpy(DECODE,"N/A");
  else
    (void) sprintf(DECODE,FORMAT,value);

  tab_prints(title,ncol,justify,DECODE);

  return;
}

/****************************************************************************/
/*!
**  Tabulated printout of a double value
**
** \param[in]  title    optional title (NULL if not defined)
** \param[in]  ncol     number of columns for the printout
** \param[in]  justify  justification flag
**                      (EJustify::LEFT, EJustify::CENTER or EJustify::RIGHT)
** \param[in]  value    Value to be written
**
*****************************************************************************/
GEOSLIB_API void tab_printd(const char*     title,
                            int             ncol,
                            const EJustify& justify,
                            double          value)
{
  st_format(CASE_DOUBLE);

  if (FFFF(value))
    (void) strcpy(DECODE,"N/A");
  else
    (void) sprintf(DECODE,FORMAT,value);

  tab_prints(title,ncol,justify,DECODE);

  return;
}

/****************************************************************************/
/*!
**  Tabulated printout of an integer value
**
** \param[in]  title    optional title (NULL if not defined)
** \param[in]  ncol     number of columns for the printout
** \param[in]  justify  justification flag
**                      (EJustify::LEFT, EJustify::CENTER or EJustify::RIGHT)
** \param[in]  value    Value to be written
**
*****************************************************************************/
GEOSLIB_API void tab_printi(const char*     title,
                            int             ncol,
                            const EJustify& justify,
                            int             value)
{
  st_format(CASE_INT);

  if (IFFFF(value))
    (void) strcpy(DECODE,"N/A");
  else
    (void) sprintf(DECODE,FORMAT,value);

  tab_prints(title,ncol,justify,DECODE);

  return;
}

/****************************************************************************/
/*!
**  Tabulated printout of a row or column value
**
** \param[in]  title    optional title (NULL if not defined)
** \param[in]  ncol     number of columns for the printout
** \param[in]  justify  justification flag
**                      (EJustify::LEFT, EJustify::CENTER or EJustify::RIGHT)
** \param[in]  mode     CASE_ROW or CASE_COL
** \param[in]  value    Value to be written
**
*****************************************************************************/
GEOSLIB_API void tab_print_rc(const char*     title,
                              int             ncol,
                              const EJustify& justify,
                              int             mode,
                              int             value)
{
  st_format(mode);

  (void) sprintf(DECODE,FORMAT,value);
  string_strip_blanks(DECODE,0);

  tab_prints(title,ncol,justify,DECODE);

  return;
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
GEOSLIB_API void print_matrix(const char   *title,
                              int     flag_limit,
                              int     bycol,
                              int     nx,
                              int     ny,
                              const double *sel,
                              const double *tab)
{
  int ix,iy,nx_util,ny_util,ny_done,multi_row,iad;

  /* Initializations */

  if (tab == (double *) NULL || nx <= 0 || ny <= 0) return;
  nx_util = (flag_limit && CST[CST_NTCOL].ival > 0) ?
    MIN(CST[CST_NTCOL].ival,nx) : nx;
  ny_util = (flag_limit && CST[CST_NTROW].ival > 0) ?
    MIN(CST[CST_NTROW].ival,ny) : ny;
  multi_row = (ny > 1 || title == NULL);

  /* Print the title (optional) */

  if (title != NULL) 
  {
    if (multi_row)
      message("%s\n",title);
    else
      message("%s ",title);
  }

  /* Print the header */

  if (multi_row)
  {
    tab_prints(NULL,1,EJustify::RIGHT," ");
    for (ix=0; ix<nx_util; ix++)
      tab_print_rc(NULL,1,EJustify::RIGHT,CASE_COL,ix+1);
    message("\n");
  }

  /* Print the contents of the array */

  ny_done = 0;
  for (iy=0; iy<ny; iy++)
  {
    if (sel != (double *) NULL && ! sel[iy]) continue;
    ny_done++;
    if (ny_done > ny_util) break;
    if (multi_row) tab_print_rc(NULL,1,EJustify::RIGHT,CASE_ROW,iy+1);
    for (ix=0; ix<nx_util; ix++)
    {
      iad = (bycol) ? iy + ny * ix : ix + nx * iy;
      tab_printg(NULL,1,EJustify::RIGHT,tab[iad]);
    }
    message("\n");
  }

  /* Print the trailor */

  if (nx != nx_util || ny != ny_util)
  {
    if (nx == nx_util)
      message("(Ncol=%d",nx);
    else
      message("(Ncol=%d[from %d]",nx_util,nx);

    if (ny == ny_util)
      message(",Nrow=%d)",ny);
    else
      message(",Nrow=%d[from %d])",ny_util,ny);
    message("\n");
  }

  return;
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
GEOSLIB_API void print_trimat(const char   *title,
                              int     mode,
                              int     neq,
                              const double *tl)
{
#define TRI(i)        (((i) * ((i) + 1)) / 2)
#define TL1(i,j)      (tl[(j)*neq+(i)-TRI(j)])  /* only for i >= j */
#define TL2(i,j)      (tl[TRI(i)+(j)])          /* only for i >= j */
  int ix,iy;

  /* Initializations */

  if (tl == (double *) NULL || neq <= 0) return;

  /* Print the title (optional) */

  if (title != NULL) message("%s\n",title);

  /* Print the header */

  tab_prints(NULL,1,EJustify::RIGHT," ");
  for (ix=0; ix<neq; ix++)
    tab_print_rc(NULL,1,EJustify::RIGHT,CASE_COL,ix+1);
  message("\n");

  /* Print the contents of the array */

  for (iy=0; iy<neq; iy++)
  {
    tab_print_rc(NULL,1,EJustify::RIGHT,CASE_ROW,iy+1);
    for (ix=0; ix<neq; ix++)
    {
      if (ix >= iy)
      {
        if (mode == 1)
          tab_printg(NULL,1,EJustify::RIGHT,TL1(ix,iy));
        else
          tab_printg(NULL,1,EJustify::RIGHT,TL2(ix,iy));
      }
      else
        tab_prints(NULL,1,EJustify::RIGHT," ");
    }
    message("\n");
  }

  return;
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
GEOSLIB_API void print_imatrix(const char   *title,
                               int     flag_limit,
                               int     bycol,
                               int     nx,
                               int     ny,
                               const double *sel,
                               const int    *tab)
{
  int ix,iy,nx_util,ny_util,ny_done,multi_row,iad;

  /* Initializations */

  if (tab == (int *) NULL || nx <= 0 || ny <= 0) return;
  nx_util = (flag_limit && CST[CST_NTCOL].ival > 0) ?
    MIN(CST[CST_NTCOL].ival,nx) : nx;
  ny_util = (flag_limit && CST[CST_NTROW].ival > 0) ?
    MIN(CST[CST_NTROW].ival,ny) : ny;
  multi_row = (ny > 1 || title == NULL);

  /* Print the title (optional) */

  if (title != NULL) 
  {
    if (multi_row)
      message("%s\n",title);
    else
      message("%s ",title);
  }

  /* Print the header */

  if (multi_row)
  {
    tab_prints(NULL,1,EJustify::RIGHT," ");
    for (ix=0; ix<nx_util; ix++)
      tab_print_rc(NULL,1,EJustify::RIGHT,CASE_COL,ix+1);
    message("\n");
  }

  /* Print the contents of the array */

  ny_done = 0;
  for (iy=0; iy<ny; iy++)
  {
    if (sel != (double *) NULL && ! sel[iy]) continue;
    ny_done++;
    if (ny_done > ny_util) break;
    if (multi_row) tab_print_rc(NULL,1,EJustify::RIGHT,CASE_ROW,iy+1);
    for (ix=0; ix<nx_util; ix++)
    {
      iad = (bycol) ? iy + ny * ix : ix + nx * iy;
      tab_printi(NULL,1,EJustify::RIGHT,tab[iad]);
    }
    message("\n");
  }

  /* Print the trailor */

  if (nx != nx_util || ny != ny_util)
  {
    if (nx == nx_util)
      message("(Ncol=%d",nx);
    else
      message("(Ncol=%d[from %d]",nx_util,nx);

    if (ny == ny_util)
      message(",Nrow=%d)",ny);
    else
      message(",Nrow=%d[from %d])",ny_util,ny);
    message("\n");
  }

  return;
}

/****************************************************************************/
/*!
**  Print a vector of real values in a matrix form
**
** \param[in]  title      Title (Optional)
** \param[in]  flag_limit 1 if CST[CST_NTCOL] is used; 0 otherwise
** \param[in]  ntab       Number of elements in the array
** \param[in]  tab        Array to be printed
**
*****************************************************************************/
GEOSLIB_API void print_vector(const char *title,
                              int     flag_limit,
                              int     ntab,
                              const double *tab)
{
  int i,j,lec,nby,flag_many;
  static int nby_def = 5;

  /* Initializations */

  if (ntab <= 0) return;
  nby = (flag_limit && CST[CST_NTCOL].ival >= 0) ?
    CST[CST_NTCOL].ival : nby_def;
  flag_many = (ntab > nby);

  if (title != NULL) 
  {
    message("%s",title);
    if (flag_many) message("\n");
  }
  for (i=lec=0; i<ntab; i+=nby)
  {
    if (flag_many) 
      message(" %2d+  ",i);
    for (j=0; j<nby; j++)
    {
      if (lec >= ntab) continue;
      message(" %10f",tab[lec]);
      lec++;
    }
    message("\n");
  }
  return;
}

GEOSLIB_API void print_vector(const char *title,
                              int flag_limit,
                              int ntab,
                              const VectorDouble& tab)
{
  print_vector(title, flag_limit, ntab, tab.data());
}

/****************************************************************************/
/*!
**  Print a vector of integer values in a matrix form
**
** \param[in]  title      Title (Optional)
** \param[in]  flag_limit 1 if CST[CST_NTCOL] is used; 0 otherwise
** \param[in]  ntab       Number of elements in the array
** \param[in]  itab       Array to be printed
**
*****************************************************************************/
GEOSLIB_API void print_ivector(const char *title,
                               int     flag_limit,
                               int     ntab,
                               const int *itab)
{
  int i,j,lec,nby,flag_many;
  static int nby_def = 5;

  /* Initializations */

  if (ntab <= 0) return;
  nby = (flag_limit && CST[CST_NTCOL].ival >= 0) ?
    CST[CST_NTCOL].ival : nby_def;
  flag_many = (ntab > nby);

  if (title != NULL) 
  {
    message("%s",title);
    if (flag_many) message("\n");
  }
  for (i=lec=0; i<ntab; i+=nby)
  {
    if (flag_many) 
      message(" %2d+  ",i);
    for (j=0; j<nby; j++)
    {
      if (lec >= ntab) continue;
      message(" %10d",itab[lec]);
      lec++;
    }
    message("\n");
  }
  return;
}

GEOSLIB_API void print_ivector(const char *title,
                               int     flag_limit,
                               int     ntab,
                               const VectorInt& itab)
{
  print_ivector(title,flag_limit,ntab,itab.data());
}

/****************************************************************************/
/*!
**  Print the names of the columns
**
** \param[in]  nx     number of columns in the matrix
** \param[in]  ranks  Indirection array (optional)
** \param[in]  names  Array of locator names
**
*****************************************************************************/
GEOSLIB_API void print_names(int    nx,
                             int   *ranks,
                             VectorString names)
{
  int ix,iix,nx_util;

  /* Initializations */

  nx_util = (CST[CST_NTCOL].ival < 0) ? nx : MIN(CST[CST_NTCOL].ival,nx);

  /* Loop on the columns */

  tab_prints(NULL,1,EJustify::RIGHT," ");
  for (iix=0; iix<nx_util; iix++)
  {
    ix = (ranks == (int *) NULL) ? iix : ranks[iix];
    tab_prints(NULL,1,EJustify::RIGHT,names[ix].c_str());
  }
  message("\n");
  return;
}

/****************************************************************************/
/*! 
**  Read a keyword
**
** \return  Rank of the keyword (starting from 0) or -1 if not recognized
**
** \param[in]  question Question to be asked
** \param[in]  nkeys    Number of authorized keys
** \param[in]  keys     List of the keywords
**
** This method is not documented on purpose. It should remain private
**
*****************************************************************************/
GEOSLIB_API int _lire_key(const char  *question,
                          int    nkeys,
                          const char **keys)
{
  int i;

label_ques:

  /* Compose the question */

  (void) sprintf(LINE,"%s ",question);
  (void) strcat(LINE,"[");
  for (i=0; i<nkeys; i++)
  {
    if (i > 0) strcat(LINE,",");
    (void) sprintf(&LINE[strlen(LINE)],"%s",keys[i]);
  }
  (void) strcat(LINE,"] : ");

  /* Read the answer */

  READ_FUNC(LINE,BUFFER);

  /* Interruption */

  if (! strcasecmp(BUFFER,"STOP")) return(-1);

  /* Check that the answer if authorized */

  string_strip_blanks(BUFFER,0);
  for (i=0; i<nkeys; i++)
    if (! strcasecmp(BUFFER,keys[i])) return(i);
  message("Error: the keyword '%s' is unknown\n",BUFFER);
  message("The only keywords authorized are : ");
  for (i=0; i<nkeys; i++) message(" %s",keys[i]);
  message("\n");
  goto label_ques;
  return -1; // Just to prevent from an eclispe warning
}

/****************************************************************************/
/*! 
**  Read astring
**
** \return  String read
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
GEOSLIB_API void _lire_string(const char *question,
                              int   flag_def,
                              const char *valdef,
                              char *answer)
{

loop:

  /* Compose the question */

  (void) sprintf(LINE,"%s ",question);
  if (flag_def)
    (void) sprintf(&LINE[strlen(LINE)],"(Def=%s) ",valdef);
  (void) strcat(LINE,": ");

  /* Read the answer */

  READ_FUNC(LINE,BUFFER);

  /* Handle the default value */

  if (strlen(BUFFER) <= 0)
  {
    if (flag_def)
    {
      (void) strcpy(answer,valdef);
    }
    else
    {
      messerr("No default value provided");
      goto loop;
    }
  }
  else
  {
    (void) strcpy(answer,BUFFER);
  }

  return;
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
GEOSLIB_API int _lire_int(const char *question,
                          int   flag_def,
                          int   valdef,
                          int   valmin,
                          int   valmax)
{
  int rep;

loop:

  /* Compose the question */

  (void) sprintf(LINE,"%s ",question);
  if (! IFFFF(valmin) && ! IFFFF(valmax) && valmin > valmax)
    valmin = valmax = ITEST;
  if (! IFFFF(valmin) && ! IFFFF(valdef) && valdef < valmin)
    valdef = valmin;
  if (! IFFFF(valmax) && ! IFFFF(valdef) && valdef > valmax)
    valdef = valmax;
  if (flag_def && ! IFFFF(valdef))
    (void) sprintf(&LINE[strlen(LINE)],"(Def=%d) ",valdef);
  if (IFFFF(valmin))
    (void) strcat(LINE,"[NA,");
  else
    (void) sprintf(&LINE[strlen(LINE)],"[%d,",valmin);
  if (IFFFF(valmax))
    (void) strcat(LINE,"NA] ");
  else
    (void) sprintf(&LINE[strlen(LINE)],"%d] ",valmax);
  (void) strcat(LINE,": ");

  /* Read the answer */

  READ_FUNC(LINE,BUFFER);

  /* Handle the default value */

  if (strlen(BUFFER) <= 0)
  {
    if (flag_def && ! IFFFF(valdef))
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
    if (! strcmp(BUFFER,"NA")) return(ITEST);
    rep = atoi(BUFFER);
  }

  /* Check the bounds */

  if (! IFFFF(valmin) && rep < valmin)
  {
    messerr("Answer (%d) must be larger than Minimum (%d)",rep,valmin);
    goto loop;
  }
  if (! IFFFF(valmax) && rep > valmax)
  {
    messerr("Answer (%d) must be smaller than Maximum (%d)",rep,valmax);
    goto loop;
  }
  return(rep);
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
GEOSLIB_API double _lire_double(const char  *question,
                                int    flag_def,
                                double valdef,
                                double valmin,
                                double valmax)
{
  double rep;

loop:

  /* Compose the question */

  (void) sprintf(LINE,"%s ",question);
  if (! FFFF(valmin) && ! FFFF(valmax) && valmin > valmax)
    valmin = valmax = TEST;
  if (! FFFF(valmin) && ! FFFF(valdef) && valdef < valmin)
    valdef = valmin;
  if (! FFFF(valmax) && ! FFFF(valdef) && valdef > valmax)
    valdef = valmax;
  if (flag_def && ! FFFF(valdef))
    (void) sprintf(&LINE[strlen(LINE)],"(Def=%lf) ",valdef); 
  if (FFFF(valmin))
    (void) strcat(LINE,"[NA,");
  else
    (void) sprintf(&LINE[strlen(LINE)],"[%lf,",valmin);
  if (FFFF(valmax))
    (void) strcat(LINE,"NA] ");
  else
    (void) sprintf(&LINE[strlen(LINE)],"%lf] ",valmax);
  (void) strcat(LINE,": ");

  /* Read the answer */

  READ_FUNC(LINE,BUFFER);

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
    if (! strcmp(BUFFER,"NA")) return(TEST);
    rep = atof(BUFFER);
  }

  /* Check the bounds */

  if (! FFFF(valmin) && rep < valmin)
  {
    messerr("Answer (%lf) must be larger than Minimum (%lf)",rep,valmin);
    goto loop;
  }
  if (! FFFF(valmax) && rep > valmax)
  {
    messerr("Answer (%lf) must be smaller than Maximum (%lf)",rep,valmax);
    goto loop;
  }
  return(rep);
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
GEOSLIB_API int _lire_logical(const char *question,
                              int   flag_def,
                              int   valdef)
{
loop:

  /* Compose the question */

  (void) sprintf(LINE,"%s ",question);
  if (flag_def && ! IFFFF(valdef))
  {
    if (valdef == 0)
      (void) strcat(LINE,"(Def=n)");
    else
      (void) strcat(LINE,"(Def=y)");
  }
  (void) strcat(LINE," [y,n] : ");

  /* Read the answer */

  READ_FUNC(LINE,BUFFER);

  /* Handle the default value */

  if (strlen(BUFFER) <= 0)
  {
    if (flag_def && ! IFFFF(valdef))
    {
      return(valdef);
    }
    else
    {
      messerr("No default value provided");
      goto loop;
    }
  }
  else
  {

    /* Check the aurhotized values */

    if (! strcasecmp(BUFFER,"Y")) return(1);
    if (! strcasecmp(BUFFER,"N")) return(0);
    message("The only authorized answers are 'y' or 'n'\n");
    goto loop;
  }
}

/****************************************************************************/
/*!
**  Returns the name of the next file in the directory
**
** \return  Name of the next file or NULL
**
** \param[in]  dirname    Name of the directory
** \param[in]  in_string  Compulsory substring
** \param[in]  ex_string  Forbidden substring
**
** This method is not documented on purpose. It should remain private
**
*****************************************************************************/
GEOSLIB_API char *_next_file(char *dirname,
                             char *in_string,
                             char *ex_string)
{
  static DIR *dir_ptr = NULL;
  struct dirent *dirent = NULL;

  if (dir_ptr == NULL)
  {
    dir_ptr = opendir(dirname);
    if (dir_ptr == NULL) return(NULL);
  }

label_loop:
  dirent = readdir(dir_ptr);
  if (dirent == NULL)
  {
    (void) closedir(dir_ptr);
    dir_ptr = NULL;
    return(NULL);
  }

  /* Discard the file whose name do not contain the compulsory string */

  if (in_string != NULL && strlen(in_string) > 0 &&
      strstr(dirent->d_name,in_string) == NULL) goto label_loop;

  /* Discard the file whose name contains the forbidden string */

  if (ex_string != NULL && strlen(ex_string) > 0 &&
      strstr(dirent->d_name,ex_string) != NULL) goto label_loop;
  return(dirent->d_name);
}

/****************************************************************************/
/*!
**  Conditionally print the progress of a procedure
**
** \param[in]  string   String to be printed
** \param[in]  ntot     Total number of samples
** \param[in]  iech     Rank of the current sample
**
*****************************************************************************/
GEOSLIB_API void mes_process(const char *string,
                             int   ntot,
                             int   iech)
{
  static int memo = 0;
  double ratio;
  int nproc,jech,percent;

  nproc = CST[CST_NPROC].ival;
  if (nproc <= 0) return;
  jech = iech + 1;

  /* Calculate the current percentage */

  ratio   = 100. * (double) jech / (double) ntot;
  percent = (int) (ratio / (double) nproc) * nproc;

  /* Conditional printout */

  if (percent != memo)
    message("%s : %d (percent)\n",string,percent);
  memo = percent;

  return;
}

/****************************************************************************/
/*! 
**  Read the next record
**
** \return Error return code
**
** This method is not documented on purpose. It should remain private
**
*****************************************************************************/
GEOSLIB_API void record_close(void)
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
GEOSLIB_API int _record_read(FILE *file,
                             const char *format,
                             ...)
{
  va_list ap;
  int error;

  error = 0;
  va_start(ap,format);
  error = _file_read(file,format,ap);

  va_end(ap);
  return(error);
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
GEOSLIB_API void print_range(const char *title,
                             int     ntab,
                             double *tab,
                             double *sel)
{
  double mini,maxi;
  int nvalid;

  if (tab == (double *) NULL || ntab <= 0) return;
  mini = maxi = TEST;
  nvalid = 0;
  ut_stats_mima(ntab,tab,sel,&nvalid,&mini,&maxi);

  /* Encode the title (if defined) */

  if (title != NULL) 
    message("%s : ",title);
  else
    message("Range : ");
  message("  ");

  if (FFFF(mini)) 
    message("NA");
  else
    message("%lf",mini);
  message(" ; ");
  if (FFFF(maxi))
    message("NA");
  else
    message("%lf",maxi);
  message(" (%d/%d)\n",nvalid,ntab);
  return;
}

/****************************************************************************/
/*!
**  Encode a real value
**
** \param[in]  string   Char array where the value is encoded
** \param[in]  ntcar    Number of characters (for real values)
** \param[in]  ntdec    Number of decimals (for real values)
** \param[in]  value    Value to be written
**
*****************************************************************************/
GEOSLIB_API void encode_printg(char  *string,
                               int    ntcar,
                               int    ntdec,
                               double value)
{
  (void) sprintf(FORMAT,"%%%d.%dlg",ntcar,ntdec);

  if (FFFF(value))
    (void) strcpy(string,"N/A");
  else
    (void) sprintf(string,FORMAT,value);
  string_strip_blanks(string,0);

  return;
}

/****************************************************************************/
/*!
**  Dump a vector of real values in a file
**  (used for debugging)
**
** \param[in]  ntab     Number of values to be dumped
** \param[in]  tab      Vector of values to be dumped out
**
*****************************************************************************/
GEOSLIB_API void file_dump(int ntab,
                           double *tab)
{
  FILE *file;
  char Local[] = "/home/drenard/Bureau/Dump_trunk";

  file = _file_open(Local,NEW);
  if (file == nullptr) return;
  for (int i = 0; i < ntab; i++)
    fprintf(file,"%d %30.20lf\n",i+1,tab[i]);
  fclose(file);
}
