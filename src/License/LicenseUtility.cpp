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
#include "License/LicenseUtility.hpp"
#include "License/MD5Utility.hpp"
#include "License/MACAddressUtility.hpp"
#include "License/RegistryUtility.hpp"
#include "License/TimeUtility.hpp"
#include "geoslib_old_f.h"
//
// Conversion from numeric to alphanumeric
//
std::string LicenseUtility::to_alpha_num(const unsigned char* src, 
                                         int length) 
{
  // Convert src into length characters in [1-9] and [A-N,P-Z] (excluding 'O')
  // 49 <=> '1'
  // 57 <=> '9'
  // 65 <=> 'A'
  // 79 <=> 'O'
  // 90 <=> 'Z'
  char one = '1';
  char nine = '9';
  char a = 'A';
  char o = 'O';
  char z = 'Z';
  std::string dst;
  char nb_char = z-a+1 - 1; // Number of letters (excluding 'O')
  nb_char += nine-one+1;    // Number of digits  (excluding '0')
  for( int i=0;i<length;i++) {
    // src[i] is in [0, 255] (255 = MAX(unsigned char))
    char ratio = (char)round(src[i]*nb_char/256);

    // Now, ratio is in [0, nb_char-1]
    ratio += one;
    // Now, ratio is in [49, 83]
    if(ratio > nine)
      ratio += a-nine-1;
    if(ratio >= o)
      ratio += 1;

    // Now, ratio is in [49-57] and [65-78,80-90]
    dst += ratio;
  }
  return dst;
}

//
// Return the activation code corresponding to the current computer
//
std::vector<std::string> LicenseUtility::get_activation_code() 
{
  std::vector<std::string> activation_codes;
  // Retrieving MAC address  
  std::vector<MacAddress> mac_addresses;
  if(MACAddressUtility::GetMACAddress(mac_addresses) == 0)
  {
    // Reverse iteration to use older/static network interfaces first
    std::vector<MacAddress>::const_reverse_iterator it = mac_addresses.rbegin();
    while(it != mac_addresses.rend()) 
    {
      // Calculating and adding to the activation code
      activation_codes.push_back(to_alpha_num(it->data(), it->size()));
      it++;
    }
  }
  else {
    messerr("Cannot retrieve MAC address");
  }

  if (activation_codes.size() <= 0)
  {
    messerr("Your computer does not have any network card connected");
    messerr("This makes the usage of Geoslib license impossible");
  }
  return activation_codes;
}
