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

class GSTLEARN_EXPORT GridIfpEn: public AOF
{
public:
  GridIfpEn(const char* filename, const Db* db = nullptr);
  GridIfpEn(const GridIfpEn& r);
  GridIfpEn& operator=(const GridIfpEn& r);
  virtual ~GridIfpEn();

  bool mustBeGrid() const override { return true; }
  bool mustBeOneVariable() const override { return false; }
  bool mustBeForNDim(int /*ndim*/) const override { return true; }
  bool mustBeForRotation(int mode) const override { return mode <= 1; }
  int  writeInFile() override;
  DbGrid* readGridFromFile() override;

private:
  void _writeLine(int mode,
                  const char *comment,
                  int valint,
                  double valrel,
                  const char *combis);
  int _readLine(int mode, const char *comment, int *valint, double *valrel);
};
