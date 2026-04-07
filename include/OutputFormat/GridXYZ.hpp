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

#include "Db/Db.hpp"
#include "OutputFormat/AOF.hpp"
#include "gstlearn_export.hpp"

namespace gstlrn
{
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

    bool mustBeForNDim(Id ndim) const override { return ndim == 2; }

    Id writeInFile() override;
  };
} // namespace gstlrn
