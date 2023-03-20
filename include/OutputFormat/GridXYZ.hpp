/******************************************************************************/
/* gstlearn C++ Library                                                       */
/*                                                                            */
/* Authors: <authors>                                                         */
/*                                                                            */
/* License: BSD 3 Clause                                                      */
/******************************************************************************/
#pragma once

#include "gstlearn_export.hpp"
#include "OutputFormat/AOF.hpp"

class Db;

class GSTLEARN_EXPORT GridXYZ: public AOF
{
public:
  GridXYZ(const char* filename, const Db* db = nullptr);
  GridXYZ(const GridXYZ& r);
  GridXYZ& operator=(const GridXYZ& r);
  virtual ~GridXYZ();

  bool mustBeGrid() const override { return true; }
  bool mustBeOneVariable() const override { return true; }
  bool mustBeForNDim(int ndim) const override { return ndim == 2; }
  int  writeInFile() override;
};
