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

class Db;

class GSTLEARN_EXPORT GridZycor: public AOF
{
public:
  GridZycor(const char* filename, const Db* db = nullptr);
  GridZycor(const GridZycor& r);
  GridZycor& operator=(const GridZycor& r);
  virtual ~GridZycor();

  bool mustBeGrid() const override { return true; }
  bool mustBeOneVariable() const override { return false; }
  bool mustBeForNDim(int ndim) const override { return ndim == 2; }
  bool mustBeForRotation(int mode) const override { return mode == 0; }
  int  writeInFile() override;
  DbGrid* readGridFromFile() override;
};
