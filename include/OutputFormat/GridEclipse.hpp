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

class GSTLEARN_EXPORT GridEclipse: public AOF
{
public:
  GridEclipse(const char* filename, const Db* db = nullptr);
  GridEclipse(const GridEclipse& r);
  GridEclipse& operator=(const GridEclipse& r);
  virtual ~GridEclipse();

  bool mustBeGrid() const override { return true; }
  bool mustBeOneVariable() const override { return true; }
  bool mustBeForNDim(int /*ndim*/) const override { return true; }
  bool mustBeForRotation(int mode) const override { return mode <= 1; }
  int  writeInFile() override;
};
