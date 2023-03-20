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

class GSTLEARN_EXPORT GridF2G: public AOF
{
public:
  GridF2G(const char* filename, const Db* db = nullptr);
  GridF2G(const GridF2G& r);
  GridF2G& operator=(const GridF2G& r);
  virtual ~GridF2G();

  bool mustBeGrid() const override { return true; }
  bool mustBeOneVariable() const override { return false; }
  bool mustBeForNDim(int ndim) const override { return ndim <= 3; }
  bool mustBeForRotation(int mode) const override { return mode <= 1; }
  DbGrid* readGridFromFile() override;
};
