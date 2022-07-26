#pragma once

#include "gstlearn_export.hpp"
#include "geoslib_define.h"

#include "Basic/VectorT.hpp"

class GSTLEARN_EXPORT LicenseUtility {

public:
  static VectorString get_activation_code();
  static String to_alpha_num(const unsigned char* src, int length);
} ;

