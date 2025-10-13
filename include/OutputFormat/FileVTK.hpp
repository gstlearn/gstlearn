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

class GSTLEARN_EXPORT FileVTK: public AOF
{
public:
  FileVTK(const char* filename, const Db* db = nullptr);
  FileVTK(const FileVTK& r);
  FileVTK& operator=(const FileVTK& r);
  virtual ~FileVTK();

  bool mustBeGrid() const override { return false; }
  bool mustBeOneVariable() const override { return true; }
  bool mustBeForNDim(Id ndim) const override { return ndim <= 3; }
  bool mustBeForRotation(Id mode) const override { return mode == 0; }
  Id  writeInFile() override;

  void setFactvar(float factvar) { _factvar = factvar; }
  void setFactx(Id factx) { _factx = factx; }
  void setFacty(Id facty) { _facty = facty; }
  void setFactz(Id factz) { _factz = factz; }
  void setFlagBinary(bool flagBinary) { _flagBinary = flagBinary; }

private:
  bool _flagBinary;
  Id _factx;
  Id _facty;
  Id _factz;
  float _factvar;
};
}