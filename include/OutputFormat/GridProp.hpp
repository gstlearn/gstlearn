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

class GSTLEARN_EXPORT GridProp: public AOF
{
public:
  GridProp(const char* filename, const Db* db);
  GridProp(const GridProp& r);
  GridProp& operator=(const GridProp& r);
  virtual ~GridProp();

  bool mustBeGrid() const override { return true; }
  bool mustBeOneVariable() const override { return false; }
  bool isAuthorized() const override;
  int  dumpFile() override;

private:
  void _writeLine(int mode,
                  const char *comment,
                  int valint,
                  double valrel,
                  const char *combis);
};
