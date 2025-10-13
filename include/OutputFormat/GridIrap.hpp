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

class GSTLEARN_EXPORT GridIrap: public AOF
{
public:
  GridIrap(const char* filename, const Db* db = nullptr);
  GridIrap(const GridIrap& r);
  GridIrap& operator=(const GridIrap& r);
  virtual ~GridIrap();

  bool mustBeGrid() const override { return true; }
  bool mustBeOneVariable() const override { return true; }
  bool mustBeForNDim(Id ndim) const override { return ndim == 2; }
  bool mustBeForRotation(Id mode) const override { return mode == 0; }
  Id  writeInFile() override;

  Id getNsamplex() const { return _nsamplex; }
  void setNsamplex(Id nsamplex) { _nsamplex = nsamplex; }
  Id getNsampley() const { return _nsampley; }
  void setNsampley(Id nsampley) { _nsampley = nsampley; }

private:
  Id _nsamplex;
  Id _nsampley;
};
}