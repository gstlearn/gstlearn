/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://github.com/gstlearn                                       */
/* License: BSD 3 clauses                                                     */
/*                                                                            */
/******************************************************************************/
#pragma once

#include "gstlearn_export.hpp"
#include "OutputFormat/AOF.hpp"

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
  bool mustBeForNDim(int ndim) const override { return ndim <= 3; }
  bool mustBeForRotation(int mode) const override { return mode == 0; }
  int  writeInFile() override;

  void setFactvar(float factvar) { _factvar = factvar; }
  void setFactx(int factx) { _factx = factx; }
  void setFacty(int facty) { _facty = facty; }
  void setFactz(int factz) { _factz = factz; }
  void setFlagBinary(bool flagBinary) { _flagBinary = flagBinary; }

private:
  bool _flagBinary;
  int _factx;
  int _facty;
  int _factz;
  float _factvar;
};
