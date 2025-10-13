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
#include "Core/Ascii.hpp"
#include "Anamorphosis/AAnam.hpp"
#include "Anamorphosis/AnamDiscreteIR.hpp"
#include "Anamorphosis/AnamHermite.hpp"
#include "Basic/OptDbg.hpp"
#include "Basic/String.hpp"
#include "Core/CSV.hpp"
#include "Core/io.hpp"
#include "Db/Db.hpp"
#include "LithoRule/Rule.hpp"
#include "Model/Model.hpp"

/*! \cond */
#define OLD 0
#define NEW 1

#define NODES(inode, i)  (nodes[6 * (inode) + (i)])
#define FROM_TYPE(inode) (nodes[6 * (inode) + 0])
#define FROM_RANK(inode) (nodes[6 * (inode) + 1])
#define FROM_VERS(inode) (nodes[6 * (inode) + 2])
#define NODE_TYPE(inode) (nodes[6 * (inode) + 3])
#define NODE_RANK(inode) (nodes[6 * (inode) + 4])
#define FACIES(inode)    (nodes[6 * (inode) + 5])

namespace gstlrn
{

static Id ASCII_BUFFER_LENGTH = 0;
static Id ASCII_BUFFER_QUANT  = 1000;
static String ASCII_BUFFER;
static FILE* FILE_MEM = NULL;
static String FILE_NAME_MEM;

/*! \endcond */

static String STUDY;
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
 ** \param[in]  vout       Returned argument
 **
 *****************************************************************************/
static Id st_record_read(const char* title, const char* format, void* vout)
{
  Id error = 1;

  if (FILE_MEM != nullptr)
  {
    error = _record_read(FILE_MEM, format, vout);
  }

  if (error > 0)
  {
    messerr("Error when reading '%s' from %s", title, FILE_NAME_MEM.data());
    print_current_line();
  }

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
static void st_record_write(const char* format, ...)
{
  va_list ap;
  String buf;
  Id long1, long2;

  va_start(ap, format);
  if (FILE_MEM != nullptr)
  {
    _file_write(FILE_MEM, format, ap);
  }
  else
  {
    _buffer_write(buf, format, ap);
    long1 = static_cast<Id>(buf.size());
    long2 = (!ASCII_BUFFER.empty()) ? static_cast<Id>(ASCII_BUFFER.size()) : 0;
    while (long1 + long2 > ASCII_BUFFER_LENGTH)
    {
      ASCII_BUFFER_LENGTH += ASCII_BUFFER_QUANT;
      ASCII_BUFFER.resize(ASCII_BUFFER_LENGTH);
    }
    (void)gslStrcat(ASCII_BUFFER, buf.data());
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
 ** \param[out] filename Output filename
 **
 ** \remark  When the rank is 0, the generic name is returned
 ** \remark  Otherwise the rank is combined in the name
 **
 *****************************************************************************/
static void st_filename_patch(const char* ref_name,
                              Id rank,
                              Id mode,
                              String& filename)
{
  if (rank == 0)
  {
    switch (mode)
    {
      case 0:
        (void)gslSPrintf(filename, "%s/%s.%s",
                         STUDY.data(), ref_name, EXT_DAT);
        break;

      case 1:
        (void)gslSPrintf(filename, "%s.%s", ref_name, EXT_OUT);
        break;

      case -1:
        (void)gslSPrintf(filename, "%s", ref_name);
        break;
    }
  }
  else
  {
    switch (mode)
    {
      case 0:
        (void)gslSPrintf(filename, "%s/%s%1d.%s",
                         STUDY.data(), ref_name, rank,
                         EXT_DAT);
        break;

      case 1:
        (void)gslSPrintf(filename, "%s%1d.%s", ref_name, rank,
                         EXT_OUT);
        break;

      case -1:
        (void)gslSPrintf(filename, "%s%1d", ref_name, rank);
        break;
    }
  }
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
void ascii_filename(const char* type, Id rank, Id mode, String& filename)
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
void ascii_study_define(const char* study)

{
  (void)gslStrcpy(STUDY, study);
}

/****************************************************************************/
/*!
 **   Close the ASCII file
 **
 ** \param[in]  file       FILE structure to be close
 **
 *****************************************************************************/
static void st_file_close(FILE* file)
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
static FILE* st_file_open(const String& filename,
                          const char* filetype,
                          Id mode,
                          bool verbose)
{
  FILE* file;
  String idtype;

  /* Open the file */

  file = FILE_MEM = _file_open(filename.data(), mode);
  FILE_NAME_MEM   = filename;

  if (file == nullptr)
  {
    if (verbose) messerr("Error when opening the file %s", filename.data());
    FILE_MEM = NULL;
    return (file);
  }

  if (OptDbg::query(EDbg::INTERFACE))
    message("Opening the File = %s\n", filename.data());

  /* Check against the file type */

  if (mode == OLD)
  {
    if (st_record_read("File Type", "%s", &idtype))
    {
      FILE_MEM = NULL;
      return (NULL);
    }
    if (idtype != filetype)
    {
      messerr("Error: in the File (%s), its Type (%s) does not match the requested one (%s)",
              filename.data(), idtype.data(), filetype);
      FILE_MEM = NULL;
      return (NULL);
    }
  }
  else
  {
    if (filetype != nullptr)
    {
      st_record_write("%s", filetype);
      st_record_write("\n");
    }
  }

  return (file);
}

/****************************************************************************/
/*!
 **   Read the Environment definition file
 **
 ** \param[in] filename  Name of the ASCII file
 ** \param[in] verbose    Verbose option if the file cannot be opened
 **
 *****************************************************************************/
void ascii_environ_read(String& filename, bool verbose)

{
  FILE* file;
  String name;
  Id debug;

  /* Opening the Data file */

  file = st_file_open(filename, "Environ", OLD, verbose);
  if (file == nullptr) return;

  /* Reading the environment */

  while (1)
  {
    if (st_record_read("Debug Keyword", "%s", &name)) goto label_end;
    if (st_record_read("Debug Value", "%ld", &debug)) goto label_end;
    String s = toUpper(String(name));
    if (debug == 1)
      OptDbg::defineByKey(s);
    else
      OptDbg::undefineByKey(s);
  }

label_end:
  st_file_close(file);
}

/****************************************************************************/
/*!
 **   Read the Simulation Characteristics
 **
 ** \param[in]  filename  Name of the ASCII file
 ** \param[in]  verbose    Verbose option if the file cannot be opened
 **
 ** \param[out]  nbsimu    Number of simulations
 ** \param[out]  nbtuba    Number of turning bands
 ** \param[out]  seed      Seed for the random number generator
 **
 *****************************************************************************/
void ascii_simu_read(String& filename,
                     bool verbose,
                     Id* nbsimu,
                     Id* nbtuba,
                     Id* seed)
{
  FILE* file;

  /* Initializations */

  (*nbsimu) = 0;
  (*nbtuba) = 100;
  (*seed)   = 0;

  /* Opening the Simulation Definition file */

  file = st_file_open(filename, "Simu", OLD, verbose);
  if (file == nullptr) return;

  /* Read the parameters */

  if (st_record_read("Number of simulations", "%ld", nbsimu)) return;
  if (st_record_read("Number of Turning Bands", "%ld", nbtuba)) return;
  if (st_record_read("Random Seed", "%ld", seed)) return;

  st_file_close(file);
}

/****************************************************************************/
/*!
 **   Check if an option is defined in the Options ASCII file
 **
 ** \return  True if the option is defined and False otherwise
 **
 ** \param[in]  filename    Name of the ASCII file
 ** \param[in]  option_name Keyword for the requested option
 ** \param[in]  verbose     Verbose option
 **
 ** \param[out]  answer      Answer
 **
 *****************************************************************************/
bool ascii_option_defined(const String& filename,
                          const char* option_name,
                          Id* answer,
                          bool verbose)
{
  FILE* file;
  String keyword;

  /* Opening the Data file */

  file = st_file_open(filename, "Option", OLD, verbose);
  if (file == nullptr) return false;

  while (1)
  {
    if (st_record_read("Debug Keyword", "%s", &keyword)) goto label_end;
    if (st_record_read("Debug Value", "%ld", answer)) goto label_end;
    if (keyword == option_name)
    {
      st_file_close(file);
      return true;
    }
  }

label_end:
  st_file_close(file);
  return false;
}

/****************************************************************************/
/*!
 **   Read a CSV file and load the results into a Db
 **
 ** \return  Pointer to the Db descriptor
 **
 ** \param[in]  filename     Name of the ASCII file
 ** \param[in]  verbose       Verbose option if the file cannot be opened
 ** \param[in]  csvfmt        CSVformat structure
 ** \param[in]  ncol_max      Maximum number of columns (or -1)
 ** \param[in]  nrow_max      Maximum number of rows (or -1)
 ** \param[in]  flagAddSampleRank True To add the rank number
 **
 *****************************************************************************/
Db* db_read_csv(const String& filename,
                const CSVformat& csvfmt,
                bool verbose,
                Id ncol_max,
                Id nrow_max,
                bool flagAddSampleRank)
{
  Db* db;
  Id ncol, nrow;
  VectorString names;
  VectorDouble tab;

  /* Initializations */

  db = nullptr;

  /* Reading the CSV file */

  if (csv_table_read(filename,
                     csvfmt, verbose, ncol_max, nrow_max, &ncol, &nrow, names, tab))
    goto label_end;

  /* Creating the Db */

  db = Db::createFromSamples(nrow, ELoadBy::SAMPLE, tab, VectorString(),
                             VectorString(), flagAddSampleRank);
  if (db == nullptr) goto label_end;

  /* Loading the names */

  for (Id i = 0; i < ncol; i++)
  {
    Id j = (flagAddSampleRank) ? i + 1 : i;
    db->setNameByUID(j, names[i]);
  }

  /* Core deallocation */

label_end:
  return (db);
}
} // namespace gstlrn
