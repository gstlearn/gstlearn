/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://gstlearn.org                                              */
/* License: BSD 3 clauses                                                     */
/*                                                                            */
/******************************************************************************/
#pragma once

#include "gstlearn_export.hpp"

class GSTLEARN_EXPORT AStringFormat
{
public:
  AStringFormat(int level = 1);
  AStringFormat(const AStringFormat& r);
  AStringFormat& operator=(const AStringFormat& r);
  virtual ~AStringFormat();

  int getLevel() const { return _level; }

private:
  int _level;
};
