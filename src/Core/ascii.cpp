/******************************************************************************/
/* COPYRIGHT ARMINES, ALL RIGHTS RESERVED                                     */
/*                                                                            */
/* THE CONTENT OF THIS WORK CONTAINS CONFIDENTIAL AND PROPRIETARY             */
/* INFORMATION OF ARMINES. ANY DUPLICATION, MODIFICATION,                     */
/* DISTRIBUTION, OR DISCLOSURE IN ANY FORM, IN WHOLE, OR IN PART, IS STRICTLY */
/*                                                                            */
/* TAG_SOURCE_CG                                                              */
/******************************************************************************/
#include "geoslib_f.h"
#include "geoslib_old_f.h"
#include "geoslib_f_private.h"
#include "geoslib_enum.h"
#include "Variogram/Vario.hpp"
#include "Anamorphosis/AAnam.hpp"
#include "Anamorphosis/AnamDiscreteDD.hpp"
#include "Anamorphosis/AnamDiscreteIR.hpp"
#include "Anamorphosis/AnamEmpirical.hpp"
#include "Anamorphosis/AnamHermite.hpp"
#include "Basic/Utilities.hpp"
#include "Basic/File.hpp"
#include "Basic/String.hpp"
#include "Basic/OptDbg.hpp"
#include "Covariances/CovAniso.hpp"
#include "Db/Db.hpp"
#include "LithoRule/Rule.hpp"
#include "Neigh/ANeighParam.hpp"
#include "Model/Model.hpp"

#include <string.h>
#include <algorithm>

/*! \cond */
#define OLD 0
#define NEW 1

#define NODES(inode,i)   (nodes[6 * (inode) + (i)])
#define FROM_TYPE(inode) (nodes[6 * (inode) + 0])
#define FROM_RANK(inode) (nodes[6 * (inode) + 1])
#define FROM_VERS(inode) (nodes[6 * (inode) + 2])
#define NODE_TYPE(inode) (nodes[6 * (inode) + 3])
#define NODE_RANK(inode) (nodes[6 * (inode) + 4])
#define FACIES(inode)    (nodes[6 * (inode) + 5])

static int ASCII_BUFFER_LENGTH = 0;
static int ASCII_BUFFER_QUANT = 1000;
static char *ASCII_BUFFER = NULL;
static FILE *FILE_MEM = NULL;
static char FILE_NAME_MEM[BUFFER_LENGTH];

/*! \endcond */

static char STUDY[BUFFER_LENGTH] = "./";
static char EXT_DAT[]         = "dat";
static char EXT_OUT[]         = "out";
static char Fichier_environ[] = "Environ";
static char Fichier_donnees[] = "Data";
static char Fichier_grid[]    = "Grid";
static char Fichier_vario[]   = "Vario";
static char Fichier_model[]   = "Model";
static char Fichier_neigh[]   = "Neigh";
static char Fichier_polygon[] = "Polygon";
static char Fichier_option[]  = "Option";
static char Fichier_rule[]    = "Rule";
static char Fichier_simu[]    = "Simu";
static char Fichier_frac[]    = "Frac";

/****************************************************************************/
/*!
 **  Read the next record
 **
 ** \return Error return code
 **
 ** \param[in]  title      Name of the quantity to be read
 ** \param[in]  format     Encoding format
 ** \param[in]  ...        Value to be written
 **
 *****************************************************************************/
static int st_record_read(const char *title, const char *format, ...)
{
  va_list ap;
  int error;

  error = 0;
  va_start(ap, format);

  if (FILE_MEM != nullptr)
  {
    error = _file_read(FILE_MEM, format, ap);
  }
  else
  {
    error = _buffer_read(&ASCII_BUFFER, format, ap);
  }

  if (error > 0)
  {
    messerr("Error when reading '%s' from %s", title, FILE_NAME_MEM);
    print_current_line();
  }

  va_end(ap);
  return (error);
}

/****************************************************************************/
/*!
 **  Write the next record
 **
 ** \param[in]  format     Encoding format
 ** \param[in]  ...        Value to be written
 **
 *****************************************************************************/
static void st_record_write(const char *format, ...)
{
  va_list ap;
  char buf[1000];
  int long1, long2;

  va_start(ap, format);
  if (FILE_MEM != nullptr)
  {
    _file_write(FILE_MEM, format, ap);
  }
  else
  {
    _buffer_write(buf, format, ap);
    long1 = static_cast<int>(strlen(buf));
    long2 = (ASCII_BUFFER != NULL) ? static_cast<int>(strlen(ASCII_BUFFER)) :
                                     0;
    while (long1 + long2 > ASCII_BUFFER_LENGTH)
    {
      ASCII_BUFFER_LENGTH += ASCII_BUFFER_QUANT;
      ASCII_BUFFER = mem_realloc(ASCII_BUFFER, ASCII_BUFFER_LENGTH, 1);
    }
    (void) gslStrcat(ASCII_BUFFER, buf);
  }

  va_end(ap);
}

/****************************************************************************/
/*!
 **   Create the File Name by patching the generic name
 **
 ** \param[in]  ref_name Reference Name
 ** \param[in]  rank     Rank of the name
 ** \param[in]  mode     Choice of the added extension
 ** \li                  0 for reading - extension ".dat"
 ** \li                  1 for writing - extension ".out"
 ** \li                 -1 no extension
 **
 ** \param[out] file_name Output filename
 **
 ** \remark  When the rank is 0, the generic name is returned
 ** \remark  Otherwise the rank is combined in the name
 **
 *****************************************************************************/
static void st_filename_patch(const char *ref_name,
                              int rank,
                              int mode,
                              char *file_name)
{
  if (rank == 0)
  {
    switch (mode)
    {
      case 0:
        (void) gslSPrintf(file_name, "%s/%s.%s", STUDY, ref_name, EXT_DAT);
        break;

      case 1:
        (void) gslSPrintf(file_name, "%s.%s", ref_name, EXT_OUT);
        break;

      case -1:
        (void) gslSPrintf(file_name, "%s", ref_name);
        break;
    }
  }
  else
  {
    switch (mode)
    {
      case 0:
        (void) gslSPrintf(file_name, "%s/%s%1d.%s",  STUDY, ref_name, rank,
                          EXT_DAT);
        break;

      case 1:
        (void) gslSPrintf(file_name, "%s%1d.%s",  ref_name, rank,
                          EXT_OUT);
        break;

      case -1:
        (void) gslSPrintf(file_name, "%s%1d",  ref_name, rank);
        break;
    }
  }
//  if (rank == 0)
//  {
//    switch (mode)
//    {
//      case 0:
//        (void) gslSPrintf(file_name, "%s/%s.%s", STUDY, ref_name, EXT_DAT);
//        break;
//
//      case 1:
//        (void) gslSPrintf(file_name, "%s/%s.%s", STUDY, ref_name, EXT_OUT);
//        break;
//
//      case -1:
//        (void) gslSPrintf(file_name, "%s/%s", STUDY, ref_name);
//        break;
//    }
//  }
//  else
//  {
//    switch (mode)
//    {
//      case 0:
//        (void) gslSPrintf(file_name, "%s/%s%1d.%s", STUDY, ref_name, rank,
//                          EXT_DAT);
//        break;
//
//      case 1:
//        (void) gslSPrintf(file_name, "%s/%s%1d.%s", STUDY, ref_name, rank,
//                          EXT_OUT);
//        break;
//
//      case -1:
//        (void) gslSPrintf(file_name, "%s/%s%1d", STUDY, ref_name, rank);
//        break;
//    }
//  }
//
  return;
}

/****************************************************************************/
/*!
 **   Returns the name of the file
 **
 ** \param[in]  type      Type of the file to be named
 ** \param[in]  rank      Rank of the file (optional)
 ** \param[in]  mode      0 for read; 1 for write
 **
 ** \param[out] filename  Output filename
 **
 *****************************************************************************/
void ascii_filename(const char *type, int rank, int mode, char *filename)
{
  if (!strcmp(type, "Environ"))
    st_filename_patch(Fichier_environ, rank, mode, filename);
  else if (!strcmp(type, "Data"))
    st_filename_patch(Fichier_donnees, rank, mode, filename);
  else if (!strcmp(type, "Grid"))
    st_filename_patch(Fichier_grid, rank, mode, filename);
  else if (!strcmp(type, "Vario"))
    st_filename_patch(Fichier_vario, rank, mode, filename);
  else if (!strcmp(type, "Model"))
    st_filename_patch(Fichier_model, rank, mode, filename);
  else if (!strcmp(type, "Neigh"))
    st_filename_patch(Fichier_neigh, rank, mode, filename);
  else if (!strcmp(type, "Rule"))
    st_filename_patch(Fichier_rule, rank, mode, filename);
  else if (!strcmp(type, "Simu"))
    st_filename_patch(Fichier_simu, rank, mode, filename);
  else if (!strcmp(type, "Polygon"))
    st_filename_patch(Fichier_polygon, rank, mode, filename);
  else if (!strcmp(type, "Option"))
    st_filename_patch(Fichier_option, rank, mode, filename);
  else if (!strcmp(type, "Frac"))
    st_filename_patch(Fichier_frac, rank, mode, filename);
  else
  {
    messageAbort("The file type %s is not referenced", type);
  }
}

/****************************************************************************/
/*!
 **   Set the Study for the test data
 **
 ** \param[in]  study Local name of the study
 **
 *****************************************************************************/
void ascii_study_define(const char *study)

{
  (void) gslStrcpy(STUDY, study);
  return;
}

/****************************************************************************/
/*!
 **   Close the ASCII file
 **
 ** \param[in]  file       FILE structure to be close
 **
 *****************************************************************************/
static void st_file_close(FILE *file)
{
  FILE_MEM = NULL;
  fclose(file);
}

/****************************************************************************/
/*!
 **   Open an ASCII file
 **
 ** \return  FILE returned pointer
 **
 ** \param[in]  filename Local file name
 ** \param[in]  filetype Type of the file (optional [NULL] when NEW)
 ** \param[in]  mode     type of file (OLD or NEW)
 ** \param[in]  verbose  Verbose option if the file cannot be opened
 **
 *****************************************************************************/
static FILE* st_file_open(const char *filename,
                          const char *filetype,
                          int mode,
                          int verbose)
{
  FILE *file;
  char idtype[LONG_SIZE];

  /* Open the file */

  file = FILE_MEM = _file_open(filename, mode);
  (void) gslStrcpy(FILE_NAME_MEM, filename);

  if (file == nullptr)
  {
    if (verbose) messerr("Error when opening the file %s", filename);
    FILE_MEM = NULL;
    return (file);
  }

  if (OptDbg::query(EDbg::INTERFACE)) message("Opening the File = %s\n", filename);

  /* Check against the file type */

  if (mode == OLD)
  {
    if (st_record_read("File Type", "%s", idtype))
    {
      FILE_MEM = NULL;
      return (NULL);
    }
    if (strcmp(idtype, filetype))
    {
      messerr(
          "Error: in the File (%s), its Type (%s) does not match the requested one (%s)",
          filename, idtype, filetype);
      FILE_MEM = NULL;
      return (NULL);
    }
  }
  else
  {
    if (filetype != NULL)
    {
      st_record_write("%s", filetype);
      st_record_write("\n");
    }
  }

  return (file);
}

/****************************************************************************/
/*!
 **   Deliver a message for successfull creation
 **
 ** \param[in]  filename Local file name
 ** \param[in]  filetype Type of the file
 ** \param[in]  verbose  Verbose option if the file cannot be opened
 **
 *****************************************************************************/
static void st_create_message(const char *filename,
                              const char *filetype,
                              int verbose)
{
  if (!verbose) return;

  if (filename != NULL)
    message("File %s (%s) created successfully\n", filename, filetype);
  else
    message("Buffer created successfully\n", filetype);

  return;
}

/****************************************************************************/
/*!
 **   Read the Environment definition file
 **
 ** \param[in] file_name  Name of the ASCII file
 ** \param[in] verbose    Verbose option if the file cannot be opened
 **
 *****************************************************************************/
void ascii_environ_read(char *file_name, int verbose)

{
  FILE *file;
  char name[10];
  int debug;

  /* Opening the Data file */

  file = st_file_open(file_name, "Environ", OLD, verbose);
  if (file == nullptr) return;

  /* Reading the environment */

  while (1)
  {
    if (st_record_read("Debug Keyword", "%s", name)) goto label_end;
    if (st_record_read("Debug Value", "%d", &debug)) goto label_end;
    String s = toUpper(String(name));
    if (debug == 1)
      OptDbg::defineByKey(s);
    else
      OptDbg::undefineByKey(s);
  }

  label_end: st_file_close(file);
  return;
}

/****************************************************************************/
/*!
 **   Read the Simulation Characteristics
 **
 ** \param[in]  file_name  Name of the ASCII file
 ** \param[in]  verbose    Verbose option if the file cannot be opened
 **
 ** \param[out]  nbsimu    Number of simulations
 ** \param[out]  nbtuba    Number of turning bands
 ** \param[out]  seed      Seed for the random number generator
 **
 *****************************************************************************/
void ascii_simu_read(char *file_name,
                     int verbose,
                     int *nbsimu,
                     int *nbtuba,
                     int *seed)
{
  FILE *file;

  /* Initializations */

  (*nbsimu) = 0;
  (*nbtuba) = 100;
  (*seed) = 0;

  /* Opening the Simulation Definition file */

  file = st_file_open(file_name, "Simu", OLD, verbose);
  if (file == nullptr) return;

  /* Read the parameters */

  if (st_record_read("Number of simulations", "%d", nbsimu)) return;
  if (st_record_read("Number of Turning Bands", "%d", nbtuba)) return;
  if (st_record_read("Random Seed", "%d", seed)) return;

  st_file_close(file);
  return;
}

/****************************************************************************/
/*!
 **   Check if an option is defined in the Options ASCII file
 **
 ** \return  1 if the option is defined and 0 otherwise
 **
 ** \param[in]  file_name    Name of the ASCII file
 ** \param[in]  verbose      Verbose option if the file cannot be opened
 ** \param[in]  option_name  Keyword for the requested option
 ** \param[in]  type         Answer type
 ** \li                      0 : Logical (returned as 0 or 1)
 ** \li                      1 : integer
 ** \li                      2 : real (returned as a double)
 **
 ** \param[out]  answer      Answer
 **
 *****************************************************************************/
int ascii_option_defined(const char *file_name,
                                         int verbose,
                                         const char *option_name,
                                         int type,
                                         void *answer)
{
  FILE *file;
  char keyword[100], keyval[100];
  double rval;
  int lrep, ival;

  /* Initializations */

  lrep = 0;

  /* Opening the Data file */

  file = st_file_open(file_name, "Option", OLD, verbose);
  if (file == nullptr) return (lrep);

  /* Implicit loop on the lines of the file */

  while (1)
  {
    if (st_record_read("Option Keyword", "%s", keyword)) goto label_end;
    if (st_record_read("Option Key-value", "%s", keyval)) goto label_end;
    if (strcmp(keyword, option_name)) continue;

    /* The keyword matches the option name */
    switch (type)
    {
      case 0:
        ival = 0;
        if (!strcmp(keyval, "Y") || !strcmp(keyval, "YES")
            || !strcmp(keyval, "y") || !strcmp(keyval, "yes")
            || atoi(keyval) == 1) ival = 1;
        *((int*) answer) = ival;
        break;

      case 1:
        ival = atoi(keyval);
        *((int*) answer) = ival;
        break;

      case 2:
        rval = atof(keyval);
        *((double*) answer) = rval;
        break;
    }
    lrep = 1;
    goto label_end;
  }

  label_end: st_file_close(file);
  return (lrep);
}

/****************************************************************************/
/*!
 **   Write a set of Fracture
 **
 ** \return  Error returned code
 **
 ** \param[in]  file_name  Name of the ASCII file
 ** \param[in]  frac       Pointer to the Frac_Environ structure to be written
 ** \param[in]  verbose    Verbose option if the file cannot be opened
 **
 *****************************************************************************/
int ascii_frac_write(const char *file_name,
                                     Frac_Environ *frac,
                                     int verbose)
{
  FILE *file;
  int family, ifault, error;

  /* Opening the Data file */

  file = st_file_open(file_name, "Frac", NEW, verbose);
  if (file == nullptr) return (1);

  /* Create the Frac_Environ structure */

  st_record_write("%d", frac->nfamilies);
  st_record_write("#", "Number of families");
  st_record_write("%d", frac->nfaults);
  st_record_write("#", "Number of main faults");
  st_record_write("%lf", frac->xmax);
  st_record_write("#", "Maximum horizontal distance");
  st_record_write("%lf", frac->ymax);
  st_record_write("#", "Maximum vertical distance");
  st_record_write("%lf", frac->deltax);
  st_record_write("#", "Dilation along the horizontal axis");
  st_record_write("%lf", frac->deltay);
  st_record_write("#", "Dilation along the vertical axis");
  st_record_write("%lf", frac->mean);
  st_record_write("#", "Mean of thickness distribution");
  st_record_write("%lf", frac->stdev);
  st_record_write("#", "Stdev of thickness distribution");

  /* Loop on the families */

  for (family = 0; family < frac->nfamilies; family++)
  {
    st_record_write("#", "Characteristics of family");
    const Frac_Fam &frac_fam = frac->frac_fams[family];
    st_record_write("%lf", frac_fam.orient);
    st_record_write("#", "Mean orientation");
    st_record_write("%lf", frac_fam.dorient);
    st_record_write("#", "Tolerance for orientation");
    st_record_write("%lf", frac_fam.theta0);
    st_record_write("#", "Reference Poisson intensity");
    st_record_write("%lf", frac_fam.alpha);
    st_record_write("#", "Power dependency between layer and intensity");
    st_record_write("%lf", frac_fam.ratcst);
    st_record_write("#", "Ratio of cste vs. shaped intensity");
    st_record_write("%lf", frac_fam.prop1);
    st_record_write("#", "Survival probability (constant term)");
    st_record_write("%lf", frac_fam.prop2);
    st_record_write("#", "Survival probability (length dependent term)");
    st_record_write("%lf", frac_fam.aterm);
    st_record_write("#", "Survival probability (cumulative length exponent)");
    st_record_write("%lf", frac_fam.bterm);
    st_record_write("#", "Survival probability (layer thickness exponent)");
    st_record_write("%lf", frac_fam.range);
    st_record_write("#", "Fracture repulsion area Range");
  }

  /* Loop on the main faults */

  for (ifault = 0; ifault < frac->nfaults; ifault++)
  {
    st_record_write("#", "Characteristics of main fault");
    const Frac_Fault &frac_fault = frac->frac_faults[ifault];
    st_record_write("%lf", frac_fault.coord);
    st_record_write("#", "Abscissae of the first Fault point");
    st_record_write("%lf", frac_fault.orient);
    st_record_write("#", "Fault orientation");

    /* Loop on the families */

    for (family = 0; family < frac->nfamilies; family++)
    {
      st_record_write("%lf", frac_fault.thetal[family]);
      st_record_write("%lf", frac_fault.rangel[family]);
      st_record_write("%lf", frac_fault.thetar[family]);
      st_record_write("%lf", frac_fault.ranger[family]);
      st_record_write("#", "Max-left, Range-Left, Max-Right, Range-Right");
    }
  }

  /* Set the error return code */

  st_create_message(file_name, "Frac", verbose);
  error = 0;

  st_file_close(file);
  return (error);
}

/****************************************************************************/
/*!
 **   Read the fracture definition file
 **
 ** \return  Pointer to the Frac_Environ structure
 **
 ** \param[in]  file_name  Name of the ASCII file
 ** \param[in]  verbose    Verbose option if the file cannot be opened
 **
 *****************************************************************************/
Frac_Environ* ascii_frac_read(const char *file_name,
                                              int verbose)
{
  Frac_Environ *frac;
  FILE *file;
  int family, ifault, nfamilies, nfaults;
  double thetal, rangel, thetar, ranger, xmax, ymax, deltax, deltay, mean,
      stdev;
  double orient, dorient, coord, theta0, alpha, prop1, prop2, aterm, bterm,
      ratcst, range;

  /* Initializations */

  frac = nullptr;

  /* Opening the Data file */

  file = st_file_open(file_name, "Frac", OLD, verbose);
  if (file == nullptr) return (frac);

  /* Create the structure */

  if (st_record_read("Number of Families", "%d", &nfamilies)) goto label_end;
  if (st_record_read("Number of Faults", "%d", &nfaults)) goto label_end;
  if (st_record_read("Maximum X-value", "%lf", &xmax)) goto label_end;
  if (st_record_read("Maximum Y-value", "%lf", &ymax)) goto label_end;
  if (st_record_read("Increment along X", "%lf", &deltax)) goto label_end;
  if (st_record_read("Increment along Y", "%lf", &deltay)) goto label_end;
  if (st_record_read("Mean value", "%lf", &mean)) goto label_end;
  if (st_record_read("Standard Deviation value", "%lf", &stdev)) goto label_end;

  frac = fracture_alloc_environ(nfamilies, xmax, ymax, deltax, deltay, mean,
                                stdev);

  /* Loop on the families */

  for (family = 0; family < nfamilies; family++)
  {
    if (st_record_read("Fracture Orientation", "%lf", &orient)) goto label_end;
    if (st_record_read("Tolerance on the Orientation", "%lf", &dorient))
      goto label_end;
    if (st_record_read("Reference Angle", "%lf", &theta0)) goto label_end;
    if (st_record_read("Fracture 'alpha' coefficient", "%lf", &alpha))
      goto label_end;
    if (st_record_read("Constant Ratio", "%lf", &ratcst)) goto label_end;
    if (st_record_read("First Proportion", "%lf", &prop1)) goto label_end;
    if (st_record_read("Second Proportion", "%lf", &prop2)) goto label_end;
    if (st_record_read("Coefficient 'aterm'", "%lf", &aterm)) goto label_end;
    if (st_record_read("Coefficient 'bterm'", "%lf", &bterm)) goto label_end;
    if (st_record_read("Range value", "%lf", &range)) goto label_end;

    fracture_update_family(frac, family, orient, dorient, theta0, alpha, ratcst,
                           prop1, prop2, aterm, bterm, range);
  }

  /* Loop on the faults */

  for (ifault = 0; ifault < nfaults; ifault++)
  {
    if (st_record_read("Fault coordinate", "%lf", &coord)) goto label_end;
    if (st_record_read("Fault orientation", "%lf", &orient)) goto label_end;
    (void) fracture_add_fault(frac, coord, orient);

    for (family = 0; family < nfamilies; family++)
    {
      if (st_record_read("Coefficient 'theta' on the left", "%lf", &thetal))
        goto label_end;
      if (st_record_read("Coefficient 'range' on the left", "%lf", &rangel))
        goto label_end;
      if (st_record_read("Coefficient 'theta' on the right", "%lf", &thetar))
        goto label_end;
      if (st_record_read("Coefficient 'range' on the right", "%lf", &ranger))
        goto label_end;

      fracture_update_fault(frac, ifault, family, thetal, thetar, rangel,
                            ranger);
    }
  }

  label_end: if (OptDbg::query(EDbg::INTERFACE)) fracture_print(frac);
  st_file_close(file);
  return (frac);
}

/****************************************************************************/
/*!
 **   Read a CSV file and load the results into a Db
 **
 ** \return  Pointer to the Db descriptor
 **
 ** \param[in]  file_name     Name of the ASCII file
 ** \param[in]  verbose       Verbose option if the file cannot be opened
 ** \param[in]  csvfmt        CSVformat structure
 ** \param[in]  ncol_max      Maximum number of columns (or -1)
 ** \param[in]  nrow_max      Maximum number of rows (or -1)
 ** \param[in]  flag_add_rank 1 To add the rank number
 **
 *****************************************************************************/
Db* db_read_csv(const char *file_name,
                const CSVformat& csvfmt,
                int verbose,
                int ncol_max,
                int nrow_max,
                int flag_add_rank)
{
  Db *db;
  int ncol, nrow;
  VectorString names;
  VectorDouble tab;

  /* Initializations */

  db = nullptr;

  /* Reading the CSV file */

  if (csv_table_read(file_name, csvfmt, verbose, ncol_max, nrow_max, &ncol, &nrow, names, tab))
    goto label_end;

  /* Creating the Db */

  db = db_create_point(nrow, ncol, ELoadBy::SAMPLE, flag_add_rank, tab);
  if (db == nullptr) goto label_end;

  /* Loading the names */

  for (int i = 0; i < ncol; i++)
  {
    int j = (flag_add_rank) ? i + 1 :
                              i;
    if (db_name_set(db, j, names[i])) messerr("Error in db_name_set");
  }

  /* Core deallocation */

  label_end: return (db);
}
