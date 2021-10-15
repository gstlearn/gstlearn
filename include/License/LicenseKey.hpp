#ifndef LICENSE_KEY_H
#define LICENSE_KEY_H

#include "License/MD5Utility.hpp"
#include <map>
#include <string>
#include <vector>

typedef std::map<std::string, std::vector<std::string> > FeatureList;
 
class LicenseKey {

public:
  static bool encodeLicenseFile(const std::string &features_file,
                                const std::string &license_file,
                                const std::string &activation_code);
  static bool isAuthorized(const std::string& feature);
  static bool checkLicense(std::string & lic_key);
  static bool checkLicenseFromKey(const std::string & lic_key);
  static bool checkLicenseFromFile(const std::string & lic_file,
                                   std::string & lic_key);
  static bool registerLicense(const std::string & target);
  static bool registerLicenseFromKey(const std::string & target,
                                     const std::string & lic_key);
  static bool registerLicenseFromFile(const std::string & target,
                                      const std::string & lic_file);
  
private:
  LicenseKey () { }
  ~LicenseKey () { }
  static std::string firstLine(const std::string & file);
  static void writeRegistry(const std::string &key,
                            const std::string &string);
  static std::string readRegistry(const std::string &key);
  static std::string readEnvironVariable(const std::string &varname);
  static void messageLicenseInvalid(void);
  static std::string stringClean(const std::string str);
  static bool isValidVersion(const std::string &version);
  static std::string appendVersion(const std::string & version,
                                   md5_state_t & state);
  static std::string appendFeature(const std::string feature,
                                   const std::string target,
                                   md5_state_t & state);
  static std::string  appendActivationCode(const std::string activation_code,
                                           md5_state_t & state);
  static std::string encrypt(const md5_state_t state);

private: 
  static std::string _targetName;
  static FeatureList _featureTable;
} ;

#endif // LICENSE_KEY_H
