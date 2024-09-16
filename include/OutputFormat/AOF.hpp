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
#pragma once

#include "gstlearn_export.hpp"

#include "Basic/VectorNumT.hpp"
#include <stdio.h>

class Db;
class DbGrid;

class GSTLEARN_EXPORT AOF
{
public:
  AOF(const String& filename, const Db* db = nullptr);
  AOF(const AOF& r);
  AOF& operator=(const AOF& r);
  virtual ~AOF();

  virtual bool mustBeGrid() const { return false; }
  virtual bool mustBeOneVariable() const { return false; }
  virtual bool mustBeForNDim(int /*ndim*/) const { return true; }
  virtual bool mustBeForRotation(int /*mode*/) const { return true; }
  virtual bool isAuthorized() const;
  virtual int  writeInFile()  { return 1; }
  virtual Db* readFromFile() { return nullptr; }
  virtual DbGrid* readGridFromFile() { return nullptr; }

  bool isValidForGrid() const;
  bool isValidForVariable() const;
  bool isValidForNDim() const;
  bool isValidForRotation() const;

  void setCols(const VectorInt& cols) { _cols = cols; }
  void setCols(int ncol, const int* icols);
  void setCol(int icol);

  const String& getFilename() const { return _filename; }

protected:
  int  _fileWriteOpen();
  int  _fileReadOpen();
  void _fileClose();

protected:
  String _filename;
  const Db*     _db;
  const DbGrid* _dbgrid;
  VectorInt     _cols;
  FILE* _file;
};

GSTLEARN_EXPORT DbGrid* db_grid_read_f2g(const char* filename, int verbose = 0);
GSTLEARN_EXPORT int db_grid_write_zycor(const char* filename, DbGrid* db, int icol);
GSTLEARN_EXPORT DbGrid* db_grid_read_zycor(const char* filename,
                                           int verbose = 0);
GSTLEARN_EXPORT int db_grid_write_arcgis(const char* filename, DbGrid* db, int icol);
GSTLEARN_EXPORT int db_grid_write_XYZ(const char* filename, DbGrid* db, int icol);
GSTLEARN_EXPORT int db_write_vtk(const char* filename, DbGrid* db, const VectorInt& cols);
GSTLEARN_EXPORT int db_grid_write_bmp(const char* filename,
                                      DbGrid* db,
                                      int icol,
                                      int nsamplex   = 1,
                                      int nsampley   = 1,
                                      int nmult      = 1,
                                      int ncolor     = 1,
                                      int flag_low   = 1,
                                      int flag_high  = 1,
                                      double valmin  = TEST,
                                      double valmax  = TEST,
                                      int* red       = nullptr,
                                      int* green     = nullptr,
                                      int* blue      = nullptr,
                                      int mask_red   = 0,
                                      int mask_green = 0,
                                      int mask_blue  = 0,
                                      int ffff_red   = 232,
                                      int ffff_green = 232,
                                      int ffff_blue  = 0,
                                      int low_red    = 255,
                                      int low_green  = 255,
                                      int low_blue   = 255,
                                      int high_red   = 255,
                                      int high_green = 0,
                                      int high_blue  = 0);
GSTLEARN_EXPORT DbGrid* db_grid_read_bmp(const char* filename, int verbose = 0);
GSTLEARN_EXPORT int db_grid_write_irap(const char* filename,
                                       DbGrid* db,
                                       int icol,
                                       int nsamplex = 1,
                                       int nsampley = 1);
GSTLEARN_EXPORT int db_grid_write_ifpen(const char* filename, DbGrid* db, int ncol, int* icols);
GSTLEARN_EXPORT DbGrid* db_grid_read_ifpen(const char* filename,
                                           int verbose = 0);
GSTLEARN_EXPORT int db_grid_write_eclipse(const char* filename, DbGrid* db, int icol);
GSTLEARN_EXPORT Db* db_well_read_las(const char* filename,
                                     double xwell,
                                     double ywell,
                                     double cwell,
                                     int verbose = 0);
