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
#include "geoslib_define.h"

#include "Enum/EDbg.hpp"

/**
 * Operate the list of active Debug options
 */
class GSTLEARN_EXPORT OptDbg
{
public:
  static void reset();

  static bool queryByKey(const String& name);
  static void defineByKey(const String& name);
  static void undefineByKey(const String& name);

  static bool query(const EDbg& option, bool discardForce = false);
  static void define(const EDbg& option);
  static void undefine(const EDbg& option);

  static void defineAll();
  static void undefineAll();
  static void display();

  static void setCurrentIndex(int cur_index) { _currentIndex = cur_index; }
  static bool isReferenceDefined() { return _reference >= 0; }
  static void setReference(int index) { _reference = index; }
  static int  getReference() { return _reference; }
  static bool force();

  static int getCurrentIndex() { return _currentIndex; }

private:
  static std::vector<EDbg> _dbg;
  static int _currentIndex;
  static int _reference;
};
