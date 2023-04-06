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

#ifdef SWIG
// Prevent SWIG wrapper generation from crying about REPEAT macro
// Look at diff for more explanation
#include "Basic/RepeatMacroSwig.hpp"
#else
// This one is compatible with MSVC, Mingw and g++
#include "Basic/RepeatMacroOrig.hpp"
#endif
