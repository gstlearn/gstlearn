/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://gstlearn.org                                              */
/* License: BSD 3-clause                                                      */
/*                                                                            */
/******************************************************************************/
#pragma once

#include "gstlearn_export.hpp"

class GSTLEARN_EXPORT ACalculator
{
public:
  ACalculator();
  ACalculator(const ACalculator &r) = delete;
  ACalculator& operator=(const ACalculator &r) = delete;
  virtual ~ACalculator();

  bool run();

protected:
  virtual bool _run() = 0;

  virtual bool _check() { return true; }
  virtual bool _preprocess()  { return true; }
  virtual bool _postprocess() { return true; }
  virtual void _rollback() { }
};
