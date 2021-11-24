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

#include "Basic/String.hpp"
#include "Basic/ASerializable.hpp"

  // TODO : Cross-platform way to build file path (use boost ?)
String getTestData(const String& filename)
{
  String exec_dir = ASerializable::getExecDirectory();
  // This path is compatible with CMake generation
  String filepath(exec_dir + "../../doc/data/" + filename);

  return filepath;
}
