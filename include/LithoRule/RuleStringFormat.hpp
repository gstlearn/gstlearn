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

class GSTLEARN_EXPORT RuleStringFormat: public AStringFormat
{
public:
  RuleStringFormat(int level = 1);
  RuleStringFormat(const RuleStringFormat& r);
  RuleStringFormat& operator=(const RuleStringFormat& r);
  virtual ~RuleStringFormat();

  bool getFlagProp() const { return _flagProp; }
  bool getFlagThresh() const { return _flagThresh; }
  void setFlagProp(bool flagProp) { _flagProp = flagProp; }
  void setFlagThresh(bool flagThresh) { _flagThresh = flagThresh; }

private:
  bool _flagProp;
  bool _flagThresh;
};
