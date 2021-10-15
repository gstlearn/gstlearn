#ifndef LICENSE_UTILITY_H
#define LICENSE_UTILITY_H

#include <vector>
#include <string>

class LicenseUtility {

public:
  static std::vector<std::string> get_activation_code();
  static std::string to_alpha_num(const unsigned char* src, int length);
} ;

#endif // LICENSE_UTILITY_H
