/******************************************************************************/
/* COPYRIGHT ARMINES, ALL RIGHTS RESERVED                                     */
/*                                                                            */
/* THE CONTENT OF THIS WORK CONTAINS CONFIDENTIAL AND PROPRIETARY             */
/* INFORMATION OF ARMINES. ANY DUPLICATION, MODIFICATION,                     */
/* DISTRIBUTION, OR DISCLOSURE IN ANY FORM, IN WHOLE, OR IN PART, IS STRICTLY */
/* PROHIBITED WITHOUT THE PRIOR EXPRESS WRITTEN PERMISSION OF ARMINES         */
/*                                                                            */
/* TAG_SOURCE_CG                                                              */
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
