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

class GSTLEARN_EXPORT FileLAS: public AOF
{
public:
  FileLAS(const char* filename, const Db* db = nullptr);
  FileLAS(const FileLAS& r);
  FileLAS& operator=(const FileLAS& r);
  virtual ~FileLAS();

  bool mustBeGrid() const override { return false; }
  bool mustBeOneVariable() const override { return true; }
  bool mustBeForNDim(int ndim) const override { return ndim <= 3; }
  bool mustBeForRotation(int mode) const override { return mode == 0; }
  Db* readFromFile() override;

  void setCwell(double cwell) { _cwell = cwell; }
  void setXwell(double xwell) { _xwell = xwell; }
  void setYwell(double ywell) { _ywell = ywell; }

private:
  int _readFind(int s_length, const char *target, int *numline, char *string);
  int _readNext(int s_length, int flag_up, int *numline, char *string);
  void _stringToUppercase(char *string) const;

private:
  double _xwell;
  double _ywell;
  double _cwell;
};
