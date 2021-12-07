#pragma once

#include "gstlearn_export.hpp"

#include <string>

class GSTLEARN_EXPORT RegistryUtility {
public :
  static std::string get_environ(const std::string& varname);
  static std::string get_value(const std::string& prog,
                               const std::string& key);
  static int         set_value(const std::string& prog,
                               const std::string& key,
                               const std::string& value);
};
