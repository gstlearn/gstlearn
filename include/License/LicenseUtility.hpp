#pragma once

#include "gstlearn_export.hpp"

#include <vector>
#include <string>

class GSTLEARN_EXPORT LicenseUtility {

public:
  static std::vector<std::string> get_activation_code();
  static std::string to_alpha_num(const unsigned char* src, int length);
} ;

