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
#include "geoslib_f.h"
#include "geoslib_old_f.h"
#include "geoslib_f_private.h"
#include "geoslib_define.h"

#include "Basic/AException.hpp"
#include "Basic/Utilities.hpp"
#include "Basic/File.hpp"
#include "Basic/String.hpp"
#include "Basic/ASerializable.hpp"
#include "Basic/VectorHelper.hpp"
#include "Db/Db.hpp"
#include "OutputFormat/GridEclipse.hpp"
#include "OutputFormat/GridIfpEn.hpp"
#include "OutputFormat/GridXYZ.hpp"
#include "OutputFormat/GridZycor.hpp"
#include "OutputFormat/GridIrap.hpp"
#include "OutputFormat/GridBmp.hpp"
#include "OutputFormat/GridArcGis.hpp"
#include "OutputFormat/GridF2G.hpp"
#include "OutputFormat/FileVTK.hpp"
#include "OutputFormat/FileLAS.hpp"
#include "vtk.h"

#include <string.h>
#include <sstream>
#include <string>
#include <stdio.h>
#include <algorithm>

/*! \cond */

struct CSV_Encoding
{
  FILE *file;                // Stream used for writing into CSV file
  int nitem;                 // Number of items per line
  int current;               // Rank of the current item
  int nlines;                // Number of lines printed
  bool flag_integer;         // true for Integer encoding
  char char_sep;             // Separator between consecutive fields
  String na_string;          // Substitute for NA
};

static CSV_Encoding *CSV_ENCODE = NULL;

/*! \endcond */

/****************************************************************************/
/*!
 **   Read / Write a File (Grid or Not) according to different format
 **
 ** \return  Error return code
 **
 *****************************************************************************/
int db_grid_write_XYZ(const char *filename, DbGrid *db, int icol)
{
  GridXYZ aof(filename, db);
  aof.setCol(icol);
  if (! aof.isAuthorized()) return 1;
  if (aof.writeInFile()) return 1;
  return 0;
}
int db_grid_write_zycor(const char *filename, DbGrid *db, int icol)
{
  GridZycor aof(filename, db);
  aof.setCol(icol);
  if (! aof.isAuthorized()) return 1;
  if (aof.writeInFile()) return 1;
  return 0;
}
DbGrid* db_grid_read_zycor(const char* filename, int /* verbose*/)
{
  GridZycor aof(filename);
  DbGrid* dbgrid = aof.readGridFromFile();
  return dbgrid;
}
int db_grid_write_arcgis(const char *filename, DbGrid *db, int icol)
{
  GridArcGis aof(filename, db);
  aof.setCol(icol);
  if (! aof.isAuthorized()) return 1;
  if (aof.writeInFile()) return 1;
  return 0;
}
int db_grid_write_bmp(const char *filename,
                      DbGrid *db,
                      int icol,
                      int nsamplex,
                      int nsampley,
                      int nmult,
                      int ncolor,
                      int flag_low,
                      int flag_high,
                      double valmin,
                      double valmax,
                      int *red,
                      int *green,
                      int *blue,
                      int mask_red,
                      int mask_green,
                      int mask_blue,
                      int ffff_red,
                      int ffff_green,
                      int ffff_blue,
                      int low_red,
                      int low_green,
                      int low_blue,
                      int high_red,
                      int high_green,
                      int high_blue)
{
  VectorInt reds = VH::initVInt(red, ncolor);
  VectorInt greens = VH::initVInt(green, ncolor);
  VectorInt blues = VH::initVInt(blue, ncolor);

  GridBmp aof(filename, db);
  aof.setCol(icol);
  aof.setNsamplex(nsamplex);
  aof.setNsampley(nsampley);
  aof.setNmult(nmult);
  aof.setNcolor(ncolor);
  aof.setFlagLow(flag_low);
  aof.setFlagHigh(flag_high);
  aof.setValmin(valmin);
  aof.setValmax(valmax);
  aof.setMask(mask_red, mask_green, mask_blue);
  aof.setFFFF(ffff_red, ffff_green, ffff_blue);
  aof.setLow (low_red,  low_green,  low_blue);
  aof.setHigh(high_red, high_green, high_blue);
  aof.setColors(reds, greens, blues);
  if (! aof.isAuthorized()) return 1;
  if (aof.writeInFile()) return 1;
  return 0;
}

DbGrid* db_grid_read_bmp(const char* filename, int /*verbose*/)
{
  GridBmp aof(filename);
  DbGrid* dbgrid = aof.readGridFromFile();
  return dbgrid;
}

int db_grid_write_irap(const char *filename,
                       DbGrid *db,
                       int icol,
                       int nsamplex,
                       int nsampley)
{
  GridIrap aof(filename, db);
  aof.setCol(icol);
  aof.setNsamplex(nsamplex);
  aof.setNsampley(nsampley);
  if (! aof.isAuthorized()) return 1;
  if (aof.writeInFile()) return 1;
  return 0;
}
int db_grid_write_ifpen(const char *filename, DbGrid *db, int ncol, int *icols)
{
  GridIfpEn aof(filename, db);
  aof.setCols(ncol,icols);
  if (! aof.isAuthorized()) return 1;
  if (aof.writeInFile()) return 1;
  return 0;
}
DbGrid* db_grid_read_ifpen(const char* filename, int /*verbose*/)
{
  GridIfpEn aof(filename);
  DbGrid* dbgrid = aof.readGridFromFile();
  return dbgrid;
}
int db_grid_write_eclipse(const char *filename, DbGrid *db, int icol)
{
  GridEclipse aof(filename, db);
  aof.setCol(icol);
  if (! aof.isAuthorized()) return 1;
  if (aof.writeInFile()) return 1;
  return 0;
}
int db_write_vtk(const char *filename,
                 DbGrid *db,
                 const VectorInt &cols)
{
  FileVTK aof(filename, db);
  aof.setCols(cols);
  if (! aof.isAuthorized()) return 1;
  if (aof.writeInFile()) return 1;
  return 0;
}
Db* db_well_read_las(const char *filename,
                     double xwell,
                     double ywell,
                     double cwell,
                     int /*verbose*/)
{
  FileLAS aof(filename);
  Db* db = aof.readGridFromFile();
  aof.setXwell(xwell);
  aof.setYwell(ywell);
  aof.setCwell(cwell);
  return db;
}
DbGrid* db_grid_read_f2g(const char* filename, int /* verbose*/)
{
  GridF2G aof(filename);
  DbGrid* dbgrid = aof.readGridFromFile();
  return dbgrid;
}

/****************************************************************************/
/*!
 **   Write a STRING element into the (opened) CSV file
 **
 ** \param[in]  string       String to be written
 **
 ** \remark: This function uses CSV_ENCODING static structure
 ** \remark: which must have been initiated beforehand
 **
 *****************************************************************************/
static void st_csv_print_string(const char *string)
{
  if (CSV_ENCODE == NULL)
  my_throw("You must initiate CSV_ENCODING first");

  (void) fprintf(CSV_ENCODE->file, "%s", string);
  if (CSV_ENCODE->current < CSV_ENCODE->nitem - 1)
  {
    (void) fprintf(CSV_ENCODE->file, "%c", CSV_ENCODE->char_sep);
    CSV_ENCODE->current++;
  }
  else
  {
    (void) fprintf(CSV_ENCODE->file, "\n");
    CSV_ENCODE->nlines++;
    CSV_ENCODE->current = 0;
  }
}

/****************************************************************************/
/*!
 **   Force the printing of End-Of-Line into the (opened) CSV file
 **
 ** \remark: This function uses CSV_ENCODING static structure
 ** \remark: which must have been initiated beforehand
 **
 *****************************************************************************/
static void st_csv_print_eol(void)
{
  if (CSV_ENCODE->current <= 0) return;

  (void) fprintf(CSV_ENCODE->file, "\n");
  CSV_ENCODE->current = 0;
  CSV_ENCODE->nlines++;
}

/****************************************************************************/
/*!
 **   Write a DOUBLE element into the (opened) CSV file
 **
 ** \param[in]  value        Real value to be written
 **
 ** \remark: This function uses CSV_ENCODING static structure
 ** \remark: which must have been initiated beforehand
 **
 *****************************************************************************/
void csv_print_double(double value)
{
  if (CSV_ENCODE == NULL)
  my_throw("You must initiate CSV_ENCODING first");

  if (FFFF(value))
    (void) fprintf(CSV_ENCODE->file, "%s", CSV_ENCODE->na_string.c_str());
  else
  {
    if (CSV_ENCODE->flag_integer)
      (void) fprintf(CSV_ENCODE->file, "%d", (int) value);
    else
      (void) fprintf(CSV_ENCODE->file, "%lf", value);
  }
  if (CSV_ENCODE->current < CSV_ENCODE->nitem - 1)
  {
    (void) fprintf(CSV_ENCODE->file, "%c", CSV_ENCODE->char_sep);
    CSV_ENCODE->current++;
  }
  else
  {
    (void) fprintf(CSV_ENCODE->file, "\n");
    CSV_ENCODE->nlines++;
    CSV_ENCODE->current = 0;
  }
}

/****************************************************************************/
/*!
 **   Manage the Utility to write into a CSV file
 **
 ** \return  Error return code
 **
 ** \param[in]  filename     Name of the CSV file
 ** \param[in]  csv          CSVFormat description
 ** \param[in]  mode         1 for opening File; -1 for closing File
 ** \param[in]  nitem        Number of items per line
 ** \param[in]  flag_integer true if the numerical values must be printed as integer
 ** \param[in]  verbose      Verbose flag
 **
 ** \remark: This procedure manages an internal structure (declared as static)
 ** \remark: When opened, you can use csv_print_string() or csv_print_double()
 ** \remark: in order to store items in the file
 ** \remark: Do not forget to use csv_manage(-1,...) to close the file
 **
 *****************************************************************************/
int csv_manage(const char *filename,
               const CSVformat& csv,
               int mode,
               int nitem,
               bool flag_integer,
               bool verbose)
{
  char char_sep = csv.getCharSep();
  String na_string = csv.getNaString();

  // Dispatch

  if (mode > 0)
  {
    // Initiate the CSV_ENCODE structure

    if (CSV_ENCODE != NULL)
      CSV_ENCODE = (CSV_Encoding*) mem_free((char* ) CSV_ENCODE);
    CSV_ENCODE = (CSV_Encoding*) mem_alloc(sizeof(CSV_Encoding), 1);
    CSV_ENCODE->file = gslFopen(filename, "w");
    if (CSV_ENCODE->file == nullptr)
    {
      messerr("Error when opening the CSV file %s for writing", filename);
      (void) csv_manage(filename, csv, -1, nitem, flag_integer);
      return 1;
    }
    CSV_ENCODE->nitem = nitem;
    CSV_ENCODE->current = 0;
    CSV_ENCODE->nlines = 0;
    CSV_ENCODE->flag_integer = flag_integer;
    CSV_ENCODE->char_sep = char_sep;
    CSV_ENCODE->na_string = na_string;

    // Optional printout

    if (verbose)
    {
      if (CSV_ENCODE->flag_integer)
        mestitle(1, "CSV Integer Encoding");
      else
        mestitle(1, "CSV Float Encoding\n");
      message("File Name                      = %s\n", filename);
      message("Number of items per line       = %d\n", CSV_ENCODE->nitem);
      message("Separator between items        = %s\n", CSV_ENCODE->char_sep);
      message("String for missing information = %s\n", CSV_ENCODE->na_string.c_str());
    }
  }
  else
  {
    // Write the last record (if necessary)
    st_csv_print_eol();

    if (CSV_ENCODE->file != NULL) fclose(CSV_ENCODE->file);

    // Option printout
    if (verbose)
    {
      if (CSV_ENCODE->flag_integer)
        message("CSV Integer Encoding : Summary\n");
      else
        message("CSV Float Encoding : Summary\n");
      message("Number of lines successfully written = %d\n",
              CSV_ENCODE->nlines);
    }

    if (CSV_ENCODE != NULL)
      CSV_ENCODE = (CSV_Encoding*) mem_free((char* ) CSV_ENCODE);
  }
  return 0;
}

/****************************************************************************/
/*!
 **   Write the Data frame into a CSV file. Reserved for numerical data frame.
 **
 ** \return  Error return code
 **
 ** \param[in]  db           Name of the Db
 ** \param[in]  filename     Name of the CSV file
 ** \param[in]  csvfmt       CSVformat structure
 ** \param[in]  flag_allcol  1 if all the columns available must be dumped out
 ** \param[in]  flag_coor    1 if the coordinates must be dumped out
 ** \param[in]  flag_integer true if the numerical values must be printed as integer
 **
 ** \remarks: This procedure dumps the Z-variables and optionally the X-variables
 **
 *****************************************************************************/
int db_write_csv(Db *db,
                 const char *filename,
                 const CSVformat& csvfmt,
                 int flag_allcol,
                 int flag_coor,
                 bool flag_integer)
{
  if (db == nullptr) return 1;
  int ncol = db->getColumnNumber();
  int ndim = db->getNDim();
  int nech = db->getSampleNumber();
  int nvar = db->getLocNumber(ELoc::Z);
  bool flag_header = csvfmt.getFlagHeader();

  // Count the number of items per line

  int nitem = 0;
  if (flag_allcol)
    nitem = ncol;
  else
  {
    nitem = nvar;
    if (flag_coor) nitem += ndim;
  }

  // Initiate the CSV_Encoding structure

  if (csv_manage(filename, csvfmt, 1, nitem, flag_integer)) return 1;

  /* Dump the header */

  if (flag_header)
  {
    // Case where all columns are dumped out

    if (flag_allcol)
    {
      for (int rank = 0; rank < ncol; rank++)
      {
        st_csv_print_string(db_name_get_by_att(db, rank).c_str());
      }
    }
    else
    {
      int rank = 0;
      if (flag_coor) for (int idim = 0; idim < ndim; idim++)
      {
        int iatt = db_attribute_identify(db, ELoc::X, idim);
        st_csv_print_string(db_name_get_by_att(db, iatt).c_str());
        rank++;
      }
      for (int ivar = 0; ivar < nvar; ivar++)
      {
        int iatt = db_attribute_identify(db, ELoc::Z, ivar);
        st_csv_print_string(db_name_get_by_att(db, iatt).c_str());
        rank++;
      }
    }
  }

  // Dump the samples (one sample per line)

  for (int iech = 0; iech < nech; iech++)
  {
    if (!db->isActive(iech)) continue;

    if (flag_allcol)
    {
      for (int rank = 0; rank < ncol; rank++)
        csv_print_double(db->getValueByColIdx(iech, rank));
    }
    else
    {
      int rank = 0;
      if (flag_coor) for (int idim = 0; idim < ndim; idim++)
      {
        int iatt = db_attribute_identify(db, ELoc::X, idim);
        csv_print_double(db->getCoordinate(iech, iatt));
        rank++;
      }
      for (int ivar = 0; ivar < nvar; ivar++)
      {
        int iatt = db_attribute_identify(db, ELoc::Z, ivar);
        csv_print_double(db->getLocVariable(ELoc::Z,iech, iatt));
        rank++;
      }
    }
  }

  // Close the file

  (void) csv_manage(filename, csvfmt, -1, nitem, flag_integer);

  return 0;
}

/****************************************************************************/
/*!
 **   Read the Data frame from a CSV file. Reserved for numerical data frame.
 **
 ** \return  Error return code
 **
 ** \param[in]  filename    Name of the CSV file
 ** \param[in]  csvfmt      CSVformat structure
 ** \param[in]  verbose     1 for a verbose output; 0 otherwise
 ** \param[in]  ncol_max    Maximum number of columns (or -1)
 ** \param[in]  nrow_max    Maximum number of rows (or -1)
 **
 ** \param[out]  ncol_arg   Number of columns
 ** \param[out]  nrow_arg   Number of rows
 ** \param[out]  names      Array containing the variable names
 ** \param[out]  tab        Array of values
 **
 ** \remarks The returned array 'tab' is organized by sample
 **
 *****************************************************************************/
int csv_table_read(const String &filename,
                   const CSVformat& csvfmt,
                   int verbose,
                   int ncol_max,
                   int nrow_max,
                   int *ncol_arg,
                   int *nrow_arg,
                   VectorString &names,
                   VectorDouble &tab)
{
  bool flag_header = csvfmt.getFlagHeader();
  int nskip = csvfmt.getNSkip();
  char char_sep = csvfmt.getCharSep();
  char char_dec = csvfmt.getCharDec();
  String na_string = csvfmt.getNaString();

  String line;
  String filepath = ASerializable::buildFileName(filename, true);
  // Open new stream
  std::ifstream file;
  file.open(filepath, std::ios::in);

//  std::ifstream file(filename.c_str());
  if (!file.is_open())
  {
    messerr("Error when opening the CSV file %s for reading", filename.c_str());
    return 1;
  }

  // Remove windows stuff at the file beginning
  skipBOM(file);

  // Initialization
  names.clear();
  tab.clear();
  int ncol = 0;

  // Define the variable names
  if (flag_header)
  {
    std::getline(file, line);
    if (!line.empty())
    {
      line = trimRight(line);
      std::istringstream iss(line);
      std::string word;
      while (std::getline(iss, word, char_sep))
      {
        word = trim(word, "\"\'");
        word = trim(word);
        names.push_back(word);
        if (verbose) message("Column Name (%d): %s\n", ncol + 1, word.c_str());
        ncol++;
        if (ncol_max > 0 && ncol >= ncol_max) break;
      }
    }

    if (verbose) message("Number of columns = %d\n", ncol);
  }

  // Skip some lines (optional)
  if (nskip > 0)
  {
    int iskip = 0;
    while (iskip < nskip && !file.eof())
    {
      std::getline(file, line);
      iskip++;
    }
  }

  // Read the values:
  int ncol2 = 0;
  int nrow = 0;
  while (!file.eof())
  {
    std::getline(file, line);
    if (!line.empty())
    {
      ncol2 = 0;
      std::istringstream iss(line);
      std::string word;
      while (std::getline(iss, word, char_sep))
      {
        if (word == na_string)
          tab.push_back(TEST);
        else
          tab.push_back(toDouble(word, char_dec));
        ncol2++;
        if (ncol_max > 0 && ncol2 >= ncol_max) break;
        if (ncol > 0 && ncol2 >= ncol) break;
      }
      if (ncol <= 0) ncol = ncol2;
      nrow++;
    }
    if (nrow_max > 0 && nrow >= nrow_max) break;
  }

  // Optional printout
  if (verbose)
  {
    message("Data table read (%s) successfully\n", gslBaseName(filename, true).c_str());
    message("- Number of columns = %d\n", ncol);
    message("- Number of rows    = %d\n", nrow);
  }

  *ncol_arg = ncol;
  *nrow_arg = nrow;

  return 0;
}

