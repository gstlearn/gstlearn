#include "License/RegistryUtility.hpp"
#include "INIParser.hpp"
#include "Basic/File.hpp"
#include <iostream>

#if defined(_WIN32) || defined(_WIN64)
#include <windows.h>
#endif
//
// Read the value from an Environment variable
//
std::string RegistryUtility::get_environ(const std::string& varname)
{
  return gslGetEnv(varname);
}

//
// Read the Key from the Environment
//
std::string RegistryUtility::get_value(const std::string& prog,
                                       const std::string& key)
{
  std::string val;
#if defined(_WIN32) || defined(_WIN64)
  std::string filename = gslGetEnv("APPDATA");
  filename += "\\";
  filename += prog;
  filename += "\\";
  filename += "License.ini";
  /*
    std::string rkey = "SOFTWARE\\Armines\\";
    rkey += prog;
    val = hkey_get(HKEY_CURRENT_USER, rkey, key);
  */
#elif defined(__linux__) || defined(__APPLE__)
  std::string filename = gslGetEnv("HOME");
  filename += "/.";
  filename += prog;
#endif
  INIParser registry_file;
  if(registry_file.load(filename)) {
    val = registry_file.GetValue<std::string>("General", key, "");
  }
  return val;
}

//
// Write the Key into the Environment
//
int RegistryUtility::set_value(const std::string& prog,
                               const std::string& key,
                               const std::string& value)
{
  int ret = 1;
#if defined(_WIN32) || defined(_WIN64)
  std::string filename = gslGetEnv("APPDATA");
  filename += "\\";
  filename += prog;
  if (CreateDirectory(filename.c_str(), NULL) ||
      ERROR_ALREADY_EXISTS == GetLastError()) {
    filename += "\\";
    filename += "License.ini";
  }
  else {
    std::cout << "Error creating directory:" << filename << std::endl;
  }
#elif defined(__linux__) || defined(__APPLE__)
  std::string filename = gslGetEnv("HOME");
  filename += "/.";
  filename += prog;
#endif
  INIParser registry_file;
  if(registry_file.load(filename)) {
    registry_file.SetValue<std::string>("General", key, value);
    if (registry_file.save(filename))
      ret = 0;
    else
      std::cout << "Error writing to file:" << filename << std::endl;
  }
  else {
    std::cout << "Error opening file:" << filename << std::endl;
  }
  return ret;
}

