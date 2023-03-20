/******************************************************************************/
/* gstlearn C++ Library                                                       */
/*                                                                            */
/* Authors: <authors>                                                         */
/*                                                                            */
/* License: BSD 3 Clause                                                      */
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
