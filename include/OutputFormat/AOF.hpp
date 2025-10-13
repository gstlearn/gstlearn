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
#include <cstdio>

namespace gstlrn
{
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
  virtual bool mustBeForNDim(Id /*ndim*/) const { return true; }
  virtual bool mustBeForRotation(Id /*mode*/) const { return true; }
  virtual bool isAuthorized() const;
  virtual Id writeInFile() { return 1; }
  virtual Db* readFromFile() { return nullptr; }
  virtual DbGrid* readGridFromFile() { return nullptr; }

  bool isValidForGrid() const;
  bool isValidForVariable() const;
  bool isValidForNDim() const;
  bool isValidForRotation() const;

  void setCols(const VectorInt& cols) { _cols = cols; }
  void setCols(Id ncol, const Id* icols);
  void setCol(Id icol);

  const String& getFilename() const { return _filename; }

protected:
  Id _fileWriteOpen();
  Id _fileReadOpen();
  void _fileClose();

protected:
  String _filename;
  const Db* _db;
  const DbGrid* _dbgrid;
  VectorInt _cols;
  FILE* _file;
};

GSTLEARN_EXPORT DbGrid* db_grid_read_f2g(const char* filename, Id verbose = 0);
GSTLEARN_EXPORT Id db_grid_write_zycor(const char* filename, DbGrid* db, Id icol);
GSTLEARN_EXPORT DbGrid* db_grid_read_zycor(const char* filename,
                                           Id verbose = 0);
GSTLEARN_EXPORT Id db_grid_write_arcgis(const char* filename, DbGrid* db, Id icol);
GSTLEARN_EXPORT Id db_grid_write_XYZ(const char* filename, DbGrid* db, Id icol);
GSTLEARN_EXPORT Id db_write_vtk(const char* filename, DbGrid* db, const VectorInt& cols);
GSTLEARN_EXPORT Id db_grid_write_bmp(const char* filename,
                                      DbGrid* db,
                                      Id icol,
                                      Id nsamplex   = 1,
                                      Id nsampley   = 1,
                                      Id nmult      = 1,
                                      Id ncolor     = 1,
                                      Id flag_low   = 1,
                                      Id flag_high  = 1,
                                      double valmin  = TEST,
                                      double valmax  = TEST,
                                      Id* red       = nullptr,
                                      Id* green     = nullptr,
                                      Id* blue      = nullptr,
                                      Id mask_red   = 0,
                                      Id mask_green = 0,
                                      Id mask_blue  = 0,
                                      Id ffff_red   = 232,
                                      Id ffff_green = 232,
                                      Id ffff_blue  = 0,
                                      Id low_red    = 255,
                                      Id low_green  = 255,
                                      Id low_blue   = 255,
                                      Id high_red   = 255,
                                      Id high_green = 0,
                                      Id high_blue  = 0);
GSTLEARN_EXPORT DbGrid* db_grid_read_bmp(const char* filename, Id verbose = 0);
GSTLEARN_EXPORT Id db_grid_write_irap(const char* filename,
                                       DbGrid* db,
                                       Id icol,
                                       Id nsamplex = 1,
                                       Id nsampley = 1);
GSTLEARN_EXPORT Id db_grid_write_ifpen(const char* filename, DbGrid* db, Id ncol, Id* icols);
GSTLEARN_EXPORT DbGrid* db_grid_read_ifpen(const char* filename,
                                           Id verbose = 0);
GSTLEARN_EXPORT Id db_grid_write_eclipse(const char* filename, DbGrid* db, Id icol);
GSTLEARN_EXPORT Db* db_well_read_las(const char* filename,
                                     double xwell,
                                     double ywell,
                                     double cwell,
                                     Id verbose = 0);
} // namespace gstlrn