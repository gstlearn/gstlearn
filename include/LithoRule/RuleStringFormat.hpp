/******************************************************************************/
/* gstlearn C++ Library                                                       */
/*                                                                            */
/* Authors: <authors>                                                         */
/*                                                                            */
/* License: BSD 3 Clause                                                      */
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
