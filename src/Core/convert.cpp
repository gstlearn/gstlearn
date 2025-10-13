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
#include "Basic/File.hpp"
#include "Basic/String.hpp"
#include "Basic/Utilities.hpp"
#include "Basic/VectorHelper.hpp"
#include "Core/CSV.hpp"
#include "Db/Db.hpp"
#include "Db/DbGrid.hpp"
#include "OutputFormat/AOF.hpp"
#include "OutputFormat/FileLAS.hpp"
#include "OutputFormat/FileVTK.hpp"
#include "OutputFormat/GridArcGis.hpp"
#include "OutputFormat/GridBmp.hpp"
#include "OutputFormat/GridEclipse.hpp"
#include "OutputFormat/GridF2G.hpp"
#include "OutputFormat/GridIfpEn.hpp"
#include "OutputFormat/GridIrap.hpp"
#include "OutputFormat/GridXYZ.hpp"
#include "OutputFormat/GridZycor.hpp"
#include "geoslib_define.h"

#include <cstdio>
#include <cstring>
#include <sstream>
#include <string>

/*! \cond */

namespace gstlrn
{
struct CSV_Encoding
{
  FILE* file;       // Stream used for writing into CSV file
  Id nitem;         // Number of items per line
  Id current;       // Rank of the current item
  Id nlines;        // Number of lines printed
  bool flagInteger; // true for Integer encoding
  char char_sep;    // Separator between consecutive fields
  String na_string; // Substitute for NA
};

static CSV_Encoding CSV_ENCODE;

/*! \endcond */

/****************************************************************************/
/*!
 **   Read / Write a File (Grid or Not) according to different format
 **
 ** \return  Error return code
 **
 *****************************************************************************/
Id db_grid_write_XYZ(const char* filename, DbGrid* db, Id icol)
{
  GridXYZ aof(filename, db);
  aof.setCol(icol);
  if (!aof.isAuthorized()) return 1;
  if (aof.writeInFile()) return 1;
  return 0;
}
Id db_grid_write_zycor(const char* filename, DbGrid* db, Id icol)
{
  GridZycor aof(filename, db);
  aof.setCol(icol);
  if (!aof.isAuthorized()) return 1;
  if (aof.writeInFile()) return 1;
  return 0;
}
DbGrid* db_grid_read_zycor(const char* filename, Id /* verbose*/)
{
  GridZycor aof(filename);
  DbGrid* dbgrid = aof.readGridFromFile();
  return dbgrid;
}
Id db_grid_write_arcgis(const char* filename, DbGrid* db, Id icol)
{
  GridArcGis aof(filename, db);
  aof.setCol(icol);
  if (!aof.isAuthorized()) return 1;
  if (aof.writeInFile()) return 1;
  return 0;
}
Id db_grid_write_bmp(const char* filename,
                     DbGrid* db,
                     Id icol,
                     Id nsamplex,
                     Id nsampley,
                     Id nmult,
                     Id ncolor,
                     Id flag_low,
                     Id flag_high,
                     double valmin,
                     double valmax,
                     Id* red,
                     Id* green,
                     Id* blue,
                     Id mask_red,
                     Id mask_green,
                     Id mask_blue,
                     Id ffff_red,
                     Id ffff_green,
                     Id ffff_blue,
                     Id low_red,
                     Id low_green,
                     Id low_blue,
                     Id high_red,
                     Id high_green,
                     Id high_blue)
{
  VectorInt reds   = VH::initVInt(red, ncolor);
  VectorInt greens = VH::initVInt(green, ncolor);
  VectorInt blues  = VH::initVInt(blue, ncolor);

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
  aof.setLow(low_red, low_green, low_blue);
  aof.setHigh(high_red, high_green, high_blue);
  aof.setColors(reds, greens, blues);
  if (!aof.isAuthorized()) return 1;
  if (aof.writeInFile()) return 1;
  return 0;
}

DbGrid* db_grid_read_bmp(const char* filename, Id /*verbose*/)
{
  GridBmp aof(filename);
  DbGrid* dbgrid = aof.readGridFromFile();
  return dbgrid;
}

Id db_grid_write_irap(const char* filename,
                      DbGrid* db,
                      Id icol,
                      Id nsamplex,
                      Id nsampley)
{
  GridIrap aof(filename, db);
  aof.setCol(icol);
  aof.setNsamplex(nsamplex);
  aof.setNsampley(nsampley);
  if (!aof.isAuthorized()) return 1;
  if (aof.writeInFile()) return 1;
  return 0;
}
Id db_grid_write_ifpen(const char* filename, DbGrid* db, Id ncol, Id* icols)
{
  GridIfpEn aof(filename, db);
  aof.setCols(ncol, icols);
  if (!aof.isAuthorized()) return 1;
  if (aof.writeInFile()) return 1;
  return 0;
}
DbGrid* db_grid_read_ifpen(const char* filename, Id /*verbose*/)
{
  GridIfpEn aof(filename);
  DbGrid* dbgrid = aof.readGridFromFile();
  return dbgrid;
}
Id db_grid_write_eclipse(const char* filename, DbGrid* db, Id icol)
{
  GridEclipse aof(filename, db);
  aof.setCol(icol);
  if (!aof.isAuthorized()) return 1;
  if (aof.writeInFile()) return 1;
  return 0;
}
Id db_write_vtk(const char* filename,
                DbGrid* db,
                const VectorInt& cols)
{
  FileVTK aof(filename, db);
  aof.setCols(cols);
  if (!aof.isAuthorized()) return 1;
  if (aof.writeInFile()) return 1;
  return 0;
}
Db* db_well_read_las(const char* filename,
                     double xwell,
                     double ywell,
                     double cwell,
                     Id /*verbose*/)
{
  FileLAS aof(filename);
  Db* db = aof.readGridFromFile();
  aof.setXwell(xwell);
  aof.setYwell(ywell);
  aof.setCwell(cwell);
  return db;
}
DbGrid* db_grid_read_f2g(const char* filename, Id /* verbose*/)
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
static void st_csv_print_string(const char* string)
{
  (void)fprintf(CSV_ENCODE.file, "%s", string);
  if (CSV_ENCODE.current < CSV_ENCODE.nitem - 1)
  {
    (void)fprintf(CSV_ENCODE.file, "%c", CSV_ENCODE.char_sep);
    CSV_ENCODE.current++;
  }
  else
  {
    (void)fprintf(CSV_ENCODE.file, "\n");
    CSV_ENCODE.nlines++;
    CSV_ENCODE.current = 0;
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
  if (CSV_ENCODE.current <= 0) return;

  (void)fprintf(CSV_ENCODE.file, "\n");
  CSV_ENCODE.current = 0;
  CSV_ENCODE.nlines++;
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
  if (FFFF(value))
    (void)fprintf(CSV_ENCODE.file, "%s", CSV_ENCODE.na_string.c_str());
  else
  {
    if (CSV_ENCODE.flagInteger)
      (void)fprintf(CSV_ENCODE.file, "%ld", static_cast<Id>(value));
    else
      (void)fprintf(CSV_ENCODE.file, "%lf", value);
  }
  if (CSV_ENCODE.current < CSV_ENCODE.nitem - 1)
  {
    (void)fprintf(CSV_ENCODE.file, "%c", CSV_ENCODE.char_sep);
    CSV_ENCODE.current++;
  }
  else
  {
    (void)fprintf(CSV_ENCODE.file, "\n");
    CSV_ENCODE.nlines++;
    CSV_ENCODE.current = 0;
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
 ** \param[in]  flagInteger  true if the numerical values must be printed as integer
 ** \param[in]  verbose      Verbose flag
 **
 ** \remark: This procedure manages an internal structure (declared as static)
 ** \remark: When opened, you can use csv_print_string() or csv_print_double()
 ** \remark: in order to store items in the file
 ** \remark: Do not forget to use csv_manage(-1,...) to close the file
 **
 *****************************************************************************/
Id csv_manage(const char* filename,
              const CSVformat& csv,
              Id mode,
              Id nitem,
              bool flagInteger,
              bool verbose)
{
  char char_sep    = csv.getCharSep();
  String na_string = csv.getNaString();

  // Dispatch

  if (mode > 0)
  {
    // Initiate the CSV_ENCODE structure

    CSV_ENCODE.file = gslFopen(filename, "w");
    if (CSV_ENCODE.file == nullptr)
    {
      messerr("Error when opening the CSV file %s for writing", filename);
      (void)csv_manage(filename, csv, -1, nitem, flagInteger);
      return 1;
    }
    CSV_ENCODE.nitem       = nitem;
    CSV_ENCODE.current     = 0;
    CSV_ENCODE.nlines      = 0;
    CSV_ENCODE.flagInteger = flagInteger;
    CSV_ENCODE.char_sep    = char_sep;
    CSV_ENCODE.na_string   = na_string;

    // Optional printout

    if (verbose)
    {
      if (CSV_ENCODE.flagInteger)
        mestitle(1, "CSV Integer Encoding");
      else
        mestitle(1, "CSV Float Encoding\n");
      message("File Name                      = %s\n", filename);
      message("Number of items per line       = %d\n", CSV_ENCODE.nitem);
      message("Separator between items        = %s\n", CSV_ENCODE.char_sep);
      message("String for missing information = %s\n", CSV_ENCODE.na_string.c_str());
    }
  }
  else
  {
    // Write the last record (if necessary)
    st_csv_print_eol();

    if (CSV_ENCODE.file != nullptr) fclose(CSV_ENCODE.file);

    // Option printout
    if (verbose)
    {
      if (CSV_ENCODE.flagInteger)
        message("CSV Integer Encoding : Summary\n");
      else
        message("CSV Float Encoding : Summary\n");
      message("Number of lines successfully written = %d\n",
              CSV_ENCODE.nlines);
    }
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
 ** \param[in]  flagInteger true if the numerical values must be printed as integer
 **
 ** \remarks: This procedure dumps the Z-variables and optionally the X-variables
 **
 *****************************************************************************/
Id db_write_csv(Db* db,
                const char* filename,
                const CSVformat& csvfmt,
                Id flag_allcol,
                Id flag_coor,
                bool flagInteger)
{
  if (db == nullptr) return 1;
  Id ncol          = db->getNColumn();
  Id ndim          = db->getNDim();
  Id nech          = db->getNSample();
  Id nvar          = db->getNLoc(ELoc::Z);
  bool flag_header = csvfmt.getFlagHeader();

  // Count the number of items per line

  Id nitem = 0;
  if (flag_allcol)
    nitem = ncol;
  else
  {
    nitem = nvar;
    if (flag_coor) nitem += ndim;
  }

  // Initiate the CSV_Encoding structure

  if (csv_manage(filename, csvfmt, 1, nitem, flagInteger)) return 1;

  /* Dump the header */

  if (flag_header)
  {
    // Case where all columns are dumped out

    if (flag_allcol)
    {
      for (Id rank = 0; rank < ncol; rank++)
      {
        st_csv_print_string(db->getNameByUID(rank).c_str());
      }
    }
    else
    {
      if (flag_coor)
        for (Id idim = 0; idim < ndim; idim++)
        {
          Id iatt = db->getUIDByLocator(ELoc::X, idim);
          st_csv_print_string(db->getNameByUID(iatt).c_str());
        }
      for (Id ivar = 0; ivar < nvar; ivar++)
      {
        Id iatt = db->getUIDByLocator(ELoc::Z, ivar);
        st_csv_print_string(db->getNameByUID(iatt).c_str());
      }
    }
  }

  // Dump the samples (one sample per line)

  for (Id iech = 0; iech < nech; iech++)
  {
    if (!db->isActive(iech)) continue;

    if (flag_allcol)
    {
      for (Id rank = 0; rank < ncol; rank++)
        csv_print_double(db->getValueByColIdx(iech, rank));
    }
    else
    {
      if (flag_coor)
        for (Id idim = 0; idim < ndim; idim++)
        {
          Id iatt = db->getUIDByLocator(ELoc::X, idim);
          csv_print_double(db->getCoordinate(iech, iatt));
        }
      for (Id ivar = 0; ivar < nvar; ivar++)
      {
        Id iatt = db->getUIDByLocator(ELoc::Z, ivar);
        csv_print_double(db->getZVariable(iech, iatt));
      }
    }
  }

  // Close the file

  (void)csv_manage(filename, csvfmt, -1, nitem, flagInteger);

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
Id csv_table_read(const String& filename,
                  const CSVformat& csvfmt,
                  bool verbose,
                  Id ncol_max,
                  Id nrow_max,
                  Id* ncol_arg,
                  Id* nrow_arg,
                  VectorString& names,
                  VectorDouble& tab)
{
  bool flag_header = csvfmt.getFlagHeader();
  auto nskip       = csvfmt.getNSkip();
  char char_sep    = csvfmt.getCharSep();
  char char_dec    = csvfmt.getCharDec();
  String na_string = csvfmt.getNaString();

  String line;
  // Open new stream
  std::ifstream file;
  file.open(filename, std::ios::in);

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
  Id ncol = 0;

  // Define the variable names
  if (flag_header)
  {
    // std::getline(file, line);
    gslSafeGetline(file, line);
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
    Id iskip = 0;
    while (iskip < nskip && !file.eof())
    {
      // std::getline(file, line);
      gslSafeGetline(file, line);
      iskip++;
    }
  }

  // Read the values:
  Id ncol2 = 0;
  Id nrow  = 0;
  while (!file.eof())
  {
    // std::getline(file, line);
    gslSafeGetline(file, line);
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
} // namespace gstlrn
