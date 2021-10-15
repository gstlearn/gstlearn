#ifndef REGISTRY_UTILITY_H
#define REGISTRY_UTILITY_H

#include <string>

class RegistryUtility {
public :
  static std::string get_environ(const std::string& varname);
  static std::string get_value(const std::string& prog,
                               const std::string& key);
  static int         set_value(const std::string& prog,
                               const std::string& key,
                               const std::string& value);
};

#endif // REGISTRY_UTILITY_H
