/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES Paris / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://gstlearn.org                                              */
/* License: BSD 3-clause                                                      */
/*                                                                            */
/******************************************************************************/
#pragma once

#include "gstlearn_export.hpp"

#include "geoslib_define.h"

namespace gstlrn
{ 
class GSTLEARN_EXPORT AStringFormat
{
public:
  AStringFormat(Id level = 1);
  AStringFormat(const String& title);
  AStringFormat(const AStringFormat& r);
  AStringFormat& operator=(const AStringFormat& r);
  virtual ~AStringFormat();

  Id getLevel() const { return _level; }
  bool hasTitle() const { return !_title.empty(); }
  String getTitle() const { return _title; }

  void setTitle(const String& title) { _title = title; }

private:
  Id _level;
  String _title;
};
}