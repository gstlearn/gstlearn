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

class GSTLEARN_EXPORT GridIfpEn: public AOF
{
public:
  GridIfpEn(const char* filename, const Db* db = nullptr);
  GridIfpEn(const GridIfpEn& r);
  GridIfpEn& operator=(const GridIfpEn& r);
  virtual ~GridIfpEn();

  bool mustBeGrid() const override { return true; }
  bool mustBeOneVariable() const override { return false; }
  bool mustBeForNDim(Id /*ndim*/) const override { return true; }
  bool mustBeForRotation(Id mode) const override { return mode <= 1; }
  Id  writeInFile() override;
  DbGrid* readGridFromFile() override;

private:
  void _writeLine(Id mode,
                  const char *comment,
                  Id valint,
                  double valrel,
                  const char *combis);
  Id _readLine(Id mode, const char *comment, Id *valint, double *valrel);
};
}