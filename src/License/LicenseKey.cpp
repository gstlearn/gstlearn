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
#include "License/LicenseKey.hpp"
#include "License/LicenseUtility.hpp"
#include "License/TimeUtility.hpp"
#include "License/RegistryUtility.hpp"
#include "version.h"
#include "geoslib_old_f.h"

#include <fstream>

#define DEBUG 0

/****************************************************************************/
/*!
 * Define the default value of the target
 *
 ****************************************************************************/
std::string LicenseKey::_targetName = "";

/****************************************************************************/
/*!
 * Define the default contents of the keyTable
 *
 ****************************************************************************/
FeatureList LicenseKey::_featureTable = FeatureList();

/****************************************************************************/
/*!
 * Read the first line of a file
 *
 * \param[in]  file  Path to the file
 *
 ****************************************************************************/
std::string LicenseKey::firstLine(const std::string & file)
{
  // Check if the file path is empty
  if (file.empty())
    return "";

  // Read the contents of the file
  std::ifstream fstream(file.c_str(), std::ifstream::in);
  if (! fstream.good()) 
  {
    fstream.close();
    return "";
  }
  
  // Load the first line
  std::string str;
  if (! std::getline(fstream, str)) 
  {
    fstream.close();
    return "";
  }
  
  fstream.close();
  return str;
}

/****************************************************************************/
/*!
 * Cleaning the end-of-line symbols
 *
 * \param[in]  str String to be cleaned
 *
 ****************************************************************************/
std::string LicenseKey::stringClean(const std::string str)
{
  std::string key = str;
  if (key.empty()) return(key);

  // Suppress <CR>
  if (key[key.length()-1] == '\n' || key[key.length()-1] == '\r')
    key.erase(key.length()-1);
  if (key[key.length()-1] == '\n' || key[key.length()-1] == '\r')
    key.erase(key.length()-1);

  // Suppress trailing '\0'
  if (key[key.length()-1] == '\0')
    key.erase(key.length()-1);

  // Suppress trailing ';'
  if (key[key.length()-1] == ';')
    key.erase(key.length()-1);

  // TODO : use trim global function ?
  // trim trailing spaces
  size_t endpos = key.find_last_not_of(" \t");
  key = key.substr( 0, endpos+1 );

  // trim leading spaces
  size_t startpos = key.find_first_not_of(" \t");
  key = key.substr( startpos );

  return(key);
}

/****************************************************************************/
/*!
 * Print the message when the License is invalid
 *
 ****************************************************************************/
void LicenseKey::messageLicenseInvalid(void)
{
  messerr(" ");
  messerr("gstlearn %s - %s", GSTLEARN_RELEASE, GSTLEARN_DATE);
  messerr(" ");
  messerr("Problem with the gstlearn License");
  messerr("================================");

  // Get the activation codes

  std::vector<std::string> all_codes = LicenseUtility::get_activation_code();
  if (all_codes.size() <= 0)
  {
    // Print network interface error
    messerr("Your computer does not have any network interface activated!");
    messerr("This makes the usage of gstlearn license impossible :-(");
    messerr("Please, try to activate one network interface (Wifi, Ethernet Card, other...).");
    messerr(" ");
  }
  else
  {
    // Print the default message
    messerr("You must register a valid License for gstlearn!");
    messerr("To identify yourself, please carefully note the activation code:");
    messerr(" ");
    messerr("%s", all_codes[0].c_str());
    messerr(" ");
    messerr("Provide this code to get your relevant license file produced");
    messerr(" ");

    messerr("When you receive the encoded License File, place it in a specific location...");
    messerr(" ");
    messerr("For LINUX:");
    messerr("- Store its complete path into the environment variable:");
    messerr("          setenv GSTLEARN_LICENSE complete_path");
    messerr("- Its contents will be copied in the local (hidden) file:");
    messerr("          .gstlearn (in the HOME Directory");
    messerr("You can then remove the License File");

    messerr(" ");
    messerr("In case of problem, please contact: gstlearn@groupes.mines-paristech.fr");
  }
  messerr(" ");
}

/****************************************************************************/
/*!
 * Check the registered License
 *
 * \return  Returned code (TRUE is correct)
 * 
 * \param[out] lic_key Registered license key 
 *
 ****************************************************************************/
bool LicenseKey::checkLicense(std::string & lic_key)
{
  // Check if the License Key is already defined (in the Registry)
  lic_key = readRegistry("LICKEY");
  if (checkLicenseFromKey(lic_key)) 
    return true;

  // Check if the filename is passed through the environment variable
  // named GSTLEARN_LICENSE
  std::string lic_file1 = readEnvironVariable("GSTLEARN_LICENSE");
  if (checkLicenseFromFile(lic_file1, lic_key))
    return true;

  // Wrong License
  return false;
}


/****************************************************************************/
/*!
 *  Check the given License Key and reset available features
 *
 * \return  Returned code (TRUE if License is correct)
 *
 * \param[in]  lic_key      License Key
 *
 ****************************************************************************/
bool LicenseKey::checkLicenseFromKey(const std::string & lic_key)
{
  md5_state_t state;
  md5_init(&state);
  std::string left, right, comment, serial_number;
  std::string version = "";
  bool flag_site = false;
  _featureTable.clear();

  if (lic_key.empty())
    return false;

  // Skipping the first statement
  std::istringstream line(lic_key);
  if (! std::getline(line, comment, ';')) 
    return false;
  
  // Reading the records of the License Key
  while (std::getline(line, left, ','))
  {
    if (! std::getline(line, right, ';'))
    {
      // No more pair, this is the serial number
      serial_number = stringClean(left);
      break;
    }

    // Particular case of the keyword "Site"
    if (left == "Site")
    {
      flag_site = true;
      appendFeature(left,right,state);
    }

    // Particular case of the keyword "Version"
    else if (left == "Version")
    {
      // Check the version number
      if (! isValidVersion(right))
        return false;
      appendVersion(right,state);
    }

    // Standard Feature
    else
    {
      appendFeature(left,right,state);
      _featureTable[right].push_back(left);
    }
  }
  
  // Empty serial number
  if (serial_number.empty())
    return false;

  if (flag_site)
  {
    // Encrypt the License contents (no activation code)
    std::string serial = encrypt(state);
    
    if (DEBUG) message("Serial tested: %s\n",serial.c_str());
    
    // Check serial number
    if (serial_number == serial)
      return true;
  }
  else
  {
    // Read the list of possible activation codes 
    std::vector<std::string> all_codes = LicenseUtility::get_activation_code();
    if (all_codes.size() <= 0) 
      return false;

    // Test each activation code
    md5_state_t state_ref = state;
    std::vector<std::string>::const_iterator it = all_codes.begin();
    while (it != all_codes.end()) 
    {
      state = state_ref;
      std::string activation_code(*it);
      (void)appendActivationCode(activation_code,state);

      // Encrypt the License contents
      std::string serial = encrypt(state);

      // Check serial number
      if (serial_number == serial) 
        return true;
      
      it++;
    }
    
    // In case of error: print the different possible activation codes
    messerr("The License Key does not correspond to the current machine");
    messerr("List of current activation codes available (debug) :");
    it = all_codes.begin();
    while (it != all_codes.end()) 
    {
      messerr("...%s...",it->c_str());
      it++;
    }
  }
  
  // Wrong License
  return false;
}

/****************************************************************************/
/*!
 * Check the given License File
 *
 * \param[in]  lic_file  Path to the License File
 * \param[out] lic_key Registered license key 
 *
 ****************************************************************************/
bool LicenseKey::checkLicenseFromFile(const std::string &lic_file,
                                      std::string & lic_key)
{
  // Check if the License File path is empty
  if (lic_file.empty())
    return false;

  // Read the contents of the License File
  lic_key = firstLine(lic_file);
  if (lic_key.empty())
  {
    messerr("This License File doesn't exist: %s\n", lic_file.c_str());
    return false;
  }
  
  // Check the License Key
  if (checkLicenseFromKey(lic_key))
    return true;
  
  // Wrong License
  return false;
}


/****************************************************************************/
/*!
 *  Check and register the License
 *
 * \return  Returned code (TRUE if License is correct)
 *
 * \param[in]  target       Target name
 *
 ****************************************************************************/
bool LicenseKey::registerLicense(const std::string & target)
{
  std::string lic_key;

  /* Let "RGeostats" or "Demonstration" connect without check */

  // Magic target keyword
  if (target == "Demonstration" || target == "RGeostats")
  {
    // Store the target name
    _targetName = target;
    
    return true;
  }

  if (checkLicense(lic_key))
  {
    // Register the license from the given key and target
    if (registerLicenseFromKey(target, lic_key))
      return true;

    // Wrong License (message already dumped)
    return false;
  }

  // Wrong License
  messageLicenseInvalid();
  return false;
}

/****************************************************************************/
/*!
 *  Check and register the given License Key
 *
 * \return  Returned code (TRUE if License is correct)
 *
 * \param[in]  target       Target name
 * \param[in]  lic_key      License Key
 *
 ****************************************************************************/
bool LicenseKey::registerLicenseFromKey(const std::string & target,
                                        const std::string & lic_key)
{
  // Check the License
  if (checkLicenseFromKey(lic_key))
  {
    // Check if the given target is registered
    if (_featureTable.find(target) != _featureTable.end())
    {
      // The License is correct: store it in the registry
      writeRegistry("LICKEY",lic_key);
      
      // Store the target name
      _targetName = target;
      
      return true;
    }
    
    // Magic target keyword
    if (target == "Demonstration" || target == "RGeostats")
    {
      // Store the target name
      _targetName = target;
      
      return true;
    }
    
    // Non registered target
    messerr("Non registered target: %s", target.c_str());
  }
  
  // Wrong License
  messageLicenseInvalid();
  return false;
}

/****************************************************************************/
/*!
 * Check and register the given License File
 *
 * \param[in]  target       Target name
 * \param[in]  lic_file     Path to the License File
 *
 ****************************************************************************/
bool LicenseKey::registerLicenseFromFile(const std::string & target,
                                         const std::string & lic_file)
{
  // Check if the License File path is empty
  if (lic_file.empty())
    return false;

  // Read the contents of the License File
  std::string lic_key = firstLine(lic_file);
  if (lic_key.empty())
  {
    messerr("This License File doesn't exist: %s\n", lic_file.c_str());
    return false;
  }
  
  // Check the License Key
  if (registerLicenseFromKey(target,lic_key))
    return true;

  // Wrong License (message already dumped)
  return false;
}

/****************************************************************************/
/*!
 *   Encode the License File
 *
 * \param[in]  features_file    Path to the input Features File
 * \param[in]  license_file     Path to the output License File
 * \param[in]  activation_code  Activation Code
 *
 ****************************************************************************/
bool LicenseKey::encodeLicenseFile(const std::string &features_file,
                                   const std::string &license_file,
                                   const std::string &activation_code)
{
  bool flag_site = false;
  std::ifstream feat_file(features_file.c_str(), std::ios::in);
  if (! feat_file.good())
  {
    messerr("Unable to open the input Features File : %s",features_file.c_str());
    return false;
  }
  std::ofstream lic_file(license_file.c_str(), std::ios::trunc);
  if (! lic_file.good())
  {
    messerr("Unable to open the output License File : %s",license_file.c_str());
    return false;
  }

  md5_state_t state;
  md5_init(&state);

  // Append the header statement (not digested by MD5)
  lic_file << "Never edit this License File;";

  // Append the Features File contents
  std::string left, right;
  while (std::getline(feat_file, left, ','))
  {
    if (left == "Site") flag_site = "true";
    if (! std::getline(feat_file, right))
    {
      messerr("Error in the Features File");
      return false;
    }
    lic_file << appendFeature(left,right,state) << ";";
  }
  
  // Append the version number
  std::stringstream ssversion;
  ssversion << GSTLEARN_VERSION;
  lic_file << appendVersion(ssversion.str(), state) << ";";

  // Append the activation code
  if (! flag_site)
    (void)appendActivationCode(activation_code,state);
  
  // Append serial number
  std::string serial = encrypt(state);
  lic_file << serial << ";";
  if (DEBUG) message("Serial number: %s\n",serial.c_str());

  feat_file.close();
  lic_file.close();

  return true;
}

/****************************************************************************/
/*!
 *  Write in a Registry according to a given Key
 *
 * \param[in]  key           Name of the Key
 * \param[in]  string        String to be encoded
 *
 ****************************************************************************/
void LicenseKey::writeRegistry(const std::string &key,
                               const std::string &string)
{
  RegistryUtility::set_value("gstlearn", key, string);
}

/****************************************************************************/
/*!
 * FUNCTION: readRegistry
 *
 *  Read from a Registry according to a given Key
 *
 * \param[in]  key          Name of the Key
 *
 ****************************************************************************/
std::string LicenseKey::readRegistry(const std::string &key)
{
  std::string string = RegistryUtility::get_value("gstlearn", key);
  return string;
}

/****************************************************************************/
/*!
 * FUNCTION: readEnvironVariable
 *
 *  Read from the Environment Variable GSTLEARN_LICENSE
 *
 ****************************************************************************/
std::string LicenseKey::readEnvironVariable(const std::string& varname)
{
  std::string string = RegistryUtility::get_environ(varname);
  return string;
}

/****************************************************************************/
/*!
 *  Check if feature is authorized
 *
 * \param[in]  feature  Name of the feature
 *
 ****************************************************************************/
bool LicenseKey::isAuthorized(const std::string& feature)
{
  if (DEBUG) message("%s: Testing %s...", _targetName.c_str(), feature.c_str());
  if (_targetName == "Demonstration") return true;
  if (_targetName == "RGeostats") return true;

  if (_featureTable.empty() || _targetName.empty()) {
    messageLicenseInvalid();
    return false;
  }

  std::vector<std::string> features = _featureTable[_targetName];
  if (features.empty())
  {
    messerr("\n---> Target \'%s\' is not registered!\n", _targetName.c_str());
    return false;
  }
  
  std::vector<std::string>::const_iterator iter = features.begin();
  while (iter < features.end())
  {
    // If the feature names are equal
    if (*iter == feature) return true;
  
    // If the registered feature is "any"
    if (*iter == "any") return true;

    iter++;
  }

  messerr("gstlearn %s - %s", GSTLEARN_RELEASE, GSTLEARN_DATE);
  messerr("\n---> Feature \'%s\' for Target \'%s\'",
          feature.c_str(),_targetName.c_str());
  messerr("     is not authorized with the current License!\n");
  return false;
}

/****************************************************************************/
/*!
 *  Append the version in the coding string
 *
 * \param[in]      version  Version to add
 *
 * \param[in,out]  state    MD5 string in construction
 *
 ****************************************************************************/
std::string LicenseKey::appendVersion(const std::string & version,
                                      md5_state_t & state)
{
  std::string pair  = "Version," + stringClean(version);
  if (DEBUG) message("Append Version: %s\n",pair.c_str());
  md5_append(&state, (const md5_byte_t*) pair.c_str(),
             static_cast<int> (pair.size()));
  return pair;
}

/****************************************************************************/
/*!
 *  Append the feature name for a given target in the coding string
 *
 * \param[in]  feature       Feature Name
 * \param[in]  target        Target Name
 *
 * \param[in,out]  state     MD5 string in construction
 *
 ****************************************************************************/
std::string LicenseKey::appendFeature(const std::string feature,
                                      const std::string target,
                                      md5_state_t & state)
{
  std::string pair = stringClean(feature) + "," + stringClean(target); 
  if (DEBUG) message("Append Feature: %s\n",pair.c_str());
  md5_append(&state, (const md5_byte_t*) pair.c_str(),
             static_cast<int> (pair.size()));
  return pair;
}

/****************************************************************************/
/*!
 *  Append the activation code in the coding string
 *
 * \param[in]  activation_code  Activation Code
 *
 * \param[in,out]  state        MD5 string in construction
 *
 ****************************************************************************/
std::string LicenseKey::appendActivationCode(const std::string activation_code,
                                             md5_state_t & state)
{
  std::string code = stringClean(activation_code);
  if (DEBUG) message("Append Code: %s\n",code.c_str());
  md5_append(&state, (const md5_byte_t*) code.c_str(),
             static_cast<int> (code.size()));
  return code;
}

/****************************************************************************/
/*!
 *  Check if the Version is valid
 *
 * \param[in]  version  Version
 *
 ****************************************************************************/
bool LicenseKey::isValidVersion(const std::string &version)
{
  if (version.empty()) return true;
  if (atoi(version.c_str()) < GSTLEARN_VERSION)
  {
    messerr("The License does not correspond to the correct gstlearn version: %s", version.c_str());
    return false;
  }
  return true;
}

/****************************************************************************/
/*
 *  Encrypt the MD5 string
 *
 * \param[in,out]  state        MD5 string in construction
 *
 ****************************************************************************/
std::string LicenseKey::encrypt(md5_state_t state)

{
  std::string code;
  md5_state_t state_ref = state;

  md5_byte_t digest[16];
  md5_finish(&state_ref, digest);
  code = LicenseUtility::to_alpha_num(digest, 16);
  code = stringClean(code);
  return code;
}

