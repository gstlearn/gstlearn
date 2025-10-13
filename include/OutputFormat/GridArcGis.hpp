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

class GSTLEARN_EXPORT GridArcGis: public AOF
{
public:
  GridArcGis(const char* filename, const Db* db = nullptr);
  GridArcGis(const GridArcGis& r);
  GridArcGis& operator=(const GridArcGis& r);
  virtual ~GridArcGis();

  bool mustBeGrid() const override { return true; }
  bool mustBeOneVariable() const override { return true; }
  bool mustBeForNDim(Id ndim) const override { return ndim == 2; }
  bool mustBeForRotation(Id mode) const override { return mode == 0; }
  bool isAuthorized() const override;
  Id  writeInFile() override;
};
}