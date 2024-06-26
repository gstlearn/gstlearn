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

#include "Basic/String.hpp"
#include <map>

/**
 * Operate the list of Constant options (defined by the user within an open list)
 */
class GSTLEARN_EXPORT OptCustom
{
public:
  static double query(const String& name, double valdef = 0.);
  static void   define(const String& name, double value);
  static void   undefine(const String& name);
  static void   display(void);

private:
  static std::map<const String, double> _cst;
};
