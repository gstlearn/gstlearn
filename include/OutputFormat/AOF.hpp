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
#include "Basic/Vector.hpp"
#include <stdio.h>

class Db;
class DbGrid;

class GSTLEARN_EXPORT AOF
{
public:
  AOF(const String& filename, const Db* db = nullptr);
  AOF(const AOF& r);
  AOF& operator=(const AOF& r);
  virtual ~AOF();

  virtual bool mustBeGrid() const { return false; }
  virtual bool mustBeOneVariable() const { return false; }
  virtual bool mustBeForNDim(int /*ndim*/) const { return true; }
  virtual bool mustBeForRotation(int /*mode*/) const { return true; }
  virtual bool isAuthorized() const;
  virtual int  writeInFile()  { return 1; }
  virtual Db* readFromFile() { return nullptr; }
  virtual DbGrid* readGridFromFile() { return nullptr; }

  bool isValidForGrid() const;
  bool isValidForVariable() const;
  bool isValidForNDim() const;
  bool isValidForRotation() const;

  void setCols(const VectorInt& cols) { _cols = cols; }
  void setCols(int ncol, int* icols);
  void setCol(int icol);

  const String& getFilename() const { return _filename; }

protected:
  int  _fileWriteOpen();
  int  _fileReadOpen();
  void _fileClose();

protected:
  String _filename;
  const Db*     _db;
  const DbGrid* _dbgrid;
  VectorInt     _cols;
  FILE* _file;
};
