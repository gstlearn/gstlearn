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

#include "Enum/ECst.hpp"

#include <map>

/**
 * Operate the list of Constant options
 */
class GSTLEARN_EXPORT OptCst
{
public:
  static double query(const ECst& option);
  static double queryByKey(const String& name);
  static void   define(const ECst& option, double value);
  static void   defineByKey(const String& name, double value);
  static void   display(void);

private:
//  static std::map<const ECst, double> _cst;
  static std::map<int, double> _cst;
};
