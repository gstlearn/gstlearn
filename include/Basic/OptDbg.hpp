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

#include "EOptDbg.hpp"

/**
 * Operate the list of active Debug options
 */
class GSTLEARN_EXPORT OptDbg
{
public:
  static void reset();
  static bool query(const EDbg& option);
  static bool queryByKey(const String& name);
  static void define(const EDbg& option, bool status);
  static void defineByKey(const String& name, bool status);
  static void defineAll(bool status);
  static void display(void);

  static void setIndex(int cur_index) { _currentIndex = cur_index; }
  static bool isReferenceDefined() { return _reference >= 0; }
  static void setReference(int index) { _reference = index; }
  static bool force();

private:
  static std::vector<EDbg> _dbg;
  static int _currentIndex;
  static int _reference;
};
