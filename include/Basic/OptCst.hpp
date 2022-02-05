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

#include <map>
#include "ECst.hpp"

/**
 * Operate the list of Constant options
 */
class GSTLEARN_EXPORT OptCst
{
public:
  static double query(const ECst& option);
  static double queryByKey(const String& name);
  static void define(const ECst& option, double value);
  static void defineByKey(const String& name, double value);
  static void display(void);

private:
//  static std::map<const ECst, double> _cst;
  static std::map<int, double> _cst;
};
