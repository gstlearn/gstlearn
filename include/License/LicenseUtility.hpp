#pragma once

#include "gstlearn_export.hpp"
#include "geoslib_define.h"

#include <vector>
#include <string>

class GSTLEARN_EXPORT LicenseUtility {

public:
  static VectorString get_activation_code();
  static std::string to_alpha_num(const unsigned char* src, int length);
} ;

