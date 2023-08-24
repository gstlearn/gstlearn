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
