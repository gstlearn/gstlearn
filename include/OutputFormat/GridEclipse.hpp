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

class GSTLEARN_EXPORT GridEclipse: public AOF
{
public:
  GridEclipse(const char* filename, const Db* db = nullptr);
  GridEclipse(const GridEclipse& r);
  GridEclipse& operator=(const GridEclipse& r);
  virtual ~GridEclipse();

  bool mustBeGrid() const override { return true; }
  bool mustBeOneVariable() const override { return true; }
  bool mustBeForNDim(int /*ndim*/) const override { return true; }
  bool mustBeForRotation(int mode) const override { return mode <= 1; }
  int  writeInFile() override;
};
