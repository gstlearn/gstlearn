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
#include "Basic/AStringFormat.hpp"

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
