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
#include "OutputFormat/AOF.hpp"

namespace gstlrn
{ 
class Db;

class GSTLEARN_EXPORT FileLAS: public AOF
{
public:
  FileLAS(const char* filename, const Db* db = nullptr);
  FileLAS(const FileLAS& r);
  FileLAS& operator=(const FileLAS& r);
  virtual ~FileLAS();

  bool mustBeGrid() const override { return false; }
  bool mustBeOneVariable() const override { return true; }
  bool mustBeForNDim(Id ndim) const override { return ndim <= 3; }
  bool mustBeForRotation(Id mode) const override { return mode == 0; }
  Db* readFromFile() override;

  void setCwell(double cwell) { _cwell = cwell; }
  void setXwell(double xwell) { _xwell = xwell; }
  void setYwell(double ywell) { _ywell = ywell; }

private:
  Id _readFind(I32 s_length, const char* target, Id* numline, char* string);
  Id _readNext(I32 s_length, Id flag_up, Id* numline, char* string);
  static void _stringToUppercase(char *string);

private:
  double _xwell;
  double _ywell;
  double _cwell;
};
} // namespace gstlrn
