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
#include "Basic/AStringFormat.hpp"

#include "geoslib_define.h"

class GSTLEARN_EXPORT PCAStringFormat: public AStringFormat
{
public:
  PCAStringFormat(int level = 1);
  PCAStringFormat(const PCAStringFormat& r);
  PCAStringFormat& operator=(const PCAStringFormat& r);
  virtual ~PCAStringFormat();

  bool getflagCenter() const { return _flagCenter; }
  bool getflagStats() const { return _flagStats; }
  void setflagCenter(bool flagCenter) { _flagCenter = flagCenter; }
  void setflagStats(bool flagStats) { _flagStats = flagStats; }

private:
  bool _flagCenter;
  bool _flagStats;
};
