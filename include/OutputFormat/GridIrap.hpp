/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://gstlearn.org                                              */
/* License: BSD 3 clauses                                                     */
/*                                                                            */
/******************************************************************************/
#pragma once

#include "gstlearn_export.hpp"
#include "OutputFormat/AOF.hpp"

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
  bool mustBeForNDim(int ndim) const override { return ndim == 2; }
  bool mustBeForRotation(int mode) const override { return mode == 0; }
  int  writeInFile() override;

  int getNsamplex() const { return _nsamplex; }
  void setNsamplex(int nsamplex) { _nsamplex = nsamplex; }
  int getNsampley() const { return _nsampley; }
  void setNsampley(int nsampley) { _nsampley = nsampley; }

private:
  int _nsamplex;
  int _nsampley;
};
