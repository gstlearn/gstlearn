#pragma once

#include "gstlearn_export.hpp"
#include "geoslib_define.h"
#include "Basic/VectorT.hpp"
#include "License/MD5Utility.hpp"

#include <map>

typedef std::map<String, VectorString > FeatureList;
 
class GSTLEARN_EXPORT LicenseKey {

public:
  static bool encodeLicenseFile(const String &features_file,
                                const String &license_file,
                                const String &activation_code);
  static bool isAuthorized(const String& feature);
  static bool checkLicense(String & lic_key);
  static bool checkLicenseFromKey(const String & lic_key);
  static bool checkLicenseFromFile(const String & lic_file,
                                   String & lic_key);
  static bool registerLicense(const String & target);
  static bool registerLicenseFromKey(const String & target,
                                     const String & lic_key);
  static bool registerLicenseFromFile(const String & target,
                                      const String & lic_file);
  
private:
  LicenseKey () { }
  ~LicenseKey () { }
  static String firstLine(const String & file);
  static void writeRegistry(const String &key,
                            const String &string);
  static String readRegistry(const String &key);
  static String readEnvironVariable(const String &varname);
  static void messageLicenseInvalid(void);
  static String stringClean(const String str);
  static bool isValidVersion(const String &version);
  static String appendVersion(const String & version,
                              md5_state_t & state);
  static String appendFeature(const String feature,
                              const String target,
                              md5_state_t & state);
  static String  appendActivationCode(const String activation_code,
                                      md5_state_t & state);
  static String encrypt(const md5_state_t state);

private: 
  static String _targetName;
  static FeatureList _featureTable;
} ;

