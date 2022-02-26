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

class GSTLEARN_EXPORT GridArcGis: public AOF
{
public:
  GridArcGis(const char* filename, const Db* db = nullptr);
  GridArcGis(const GridArcGis& r);
  GridArcGis& operator=(const GridArcGis& r);
  virtual ~GridArcGis();

  bool mustBeGrid() const override { return true; }
  bool mustBeOneVariable() const override { return true; }
  bool mustBeForNDim(int ndim) const override { return ndim == 2; }
  bool mustBeForRotation(int mode) const override { return mode == 0; }
  bool isAuthorized() const override;
  int  writeInFile() override;
};
