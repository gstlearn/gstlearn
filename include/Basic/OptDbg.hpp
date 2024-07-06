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
#include "geoslib_define.h"

#include "Enum/EDbg.hpp"

/**
 * Operate the list of active Debug options
 * These options correspond to various keywords chosen from a close list (see EDbg.hpp).
 * When one of these options is switched ON, some statements are printed each time
 * this particular option is concerned
 *
 * As example, to produce specific output during all Kriging steps, you should use;
 *
 *     OptCst::define(EDbg::KRIGING)
 *
 * To cancel this option, it suffices to use:
 *
 *     OptCst::undefine(EDbg::KRIGING)
 *
 * To know the current status of all these environmental parameters,
 * use the display() function.
 *
 * A complementary option is to use the Reference: this is an index such that, when the
 * rank of the target matches this number, all the flags are turned ON automatically.
 * This Reference index is provided as a 1-based number.
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
