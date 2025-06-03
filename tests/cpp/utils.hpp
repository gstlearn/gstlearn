/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2025) MINES Paris / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://gstlearn.org                                              */
/* License: BSD 3-clause                                                      */
/*                                                                            */
/******************************************************************************/

#pragma once

#include <geoslib_define.h>
#include <Basic/File.hpp>

#include <filesystem>

#ifndef CMAKE_SOURCE_DIR
#  warning "CMAKE_SOURCE_DIR not defined, using GSTLEARN_DIR env var instead"
#  define CMAKE_SOURCE_DIR gslGetEnv("GSTLEAR_DIR")
#endif // CMAKE_SOURCE_DIR

namespace gstlearn
{
  /**
   * This method returns the absolute path to a Test Data file
   * This can only be used in non-regression test (NOT in any Python or R stand-alone script)
   *
   * @param subdir Sub directory (in doc/data folder) containing the required file
   * @param filename Name of the required data file
   *
   * @return
   */
  inline String getTestData(const String& subdir, const String& filename)
  {
    std::filesystem::path p {CMAKE_SOURCE_DIR};
    return (p / "doc" / "data" / subdir / filename).string();
  }

} // namespace gstlearn
