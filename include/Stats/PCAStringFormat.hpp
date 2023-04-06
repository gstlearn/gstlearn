/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://github.com/gstlearn                                       */
/* License: BSD 3 clauses                                                     */
/*                                                                            */
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
