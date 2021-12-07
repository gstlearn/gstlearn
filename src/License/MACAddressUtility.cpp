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
#include "License/MACAddressUtility.hpp"
#include "Basic/String.hpp"
#include "geoslib_old_f.h"

#include <fstream>
#include <sstream>

/****************************************************************************
 **
 ** Returns the vector MacAddresses
 **
 *****************************************************************************/
long MACAddressUtility::GetMACAddress(std::vector<MacAddress>& result) 
{
  result.clear();
  // Call appropriate function for each platform
#if defined(_WIN32) || defined(_WIN64)
  return GetMACAddressMSW(result);
#elif defined(__linux__)
  return GetMACAddressLinux(result);
#elif defined(__APPLE__)
  return GetMACAddressMAC(result);
#endif
  // If platform is not supported then return error code
  return -1;
}

/****************************************************************************
 ** Version for WINDOWS
 ** Using GetAdaptatersAddresses from :
 ** http://msdn.microsoft.com/en-us/library/aa365915%28v=vs.85%29.aspx
 *****************************************************************************/
#if defined(_WIN32) || defined(_WIN64)

long MACAddressUtility::GetMACAddressMSW(std::vector<MacAddress>& result) 
{
  /* Declare and initialize variables */

  DWORD dwRetVal = 0;

  // Set the flags to pass to GetAdaptersAddresses
  ULONG flags = GAA_FLAG_INCLUDE_PREFIX;

  // default to unspecified address family (both)
  ULONG family = AF_UNSPEC;

  PIP_ADAPTER_ADDRESSES pAddresses = NULL;
  ULONG outBufLen = 0;
  ULONG Iterations = 0;

  PIP_ADAPTER_ADDRESSES pCurrAddresses = NULL;

  // Only IPv4 support
  family = AF_INET;
  //family = AF_INET6; // for IPv6

  // Allocate a 15 KB buffer to start with.
  outBufLen = WORKING_BUFFER_SIZE;

  do
  {
    pAddresses = (IP_ADAPTER_ADDRESSES *) MALLOC(outBufLen);
    if (pAddresses == NULL)
    {
      messerr("Not enough memory retrieving network interfaces!");
      return -1;
    }
    dwRetVal = GetAdaptersAddresses(family, flags, NULL, pAddresses, &outBufLen);

    if (dwRetVal == ERROR_BUFFER_OVERFLOW)
    {
      FREE(pAddresses);
      pAddresses = NULL;
    }
    else
    {
      break;
    }

    Iterations++;
  }
  while ((dwRetVal == ERROR_BUFFER_OVERFLOW) && (Iterations < MAX_TRIES));

  if (dwRetVal == NO_ERROR)
  {
    // If successful, output some information from the data we received
    pCurrAddresses = pAddresses;
    while (pCurrAddresses)
    {
      // Filter loopback and point to point (VPN)
      if(pCurrAddresses->IfType != IF_TYPE_SOFTWARE_LOOPBACK &&
         pCurrAddresses->IfType != IF_TYPE_PPP) {
        if (pCurrAddresses->PhysicalAddressLength >= 6)
        {
          result.push_back(MacAddress(pCurrAddresses->PhysicalAddress,
                                      pCurrAddresses->PhysicalAddressLength));
        }
        else
        {
          messerr("Cannot retrieve network interface properties for:%s",
                  pCurrAddresses->AdapterName);
        }
      }
      pCurrAddresses = pCurrAddresses->Next;
    }
  }
  else
  {
    messerr("Cannot retrieve network interfaces");
    if (pAddresses) FREE(pAddresses);
    return -1;
  }

  if (pAddresses) FREE(pAddresses);

  return 0;
}

/****************************************************************************
 ** Version for LINUX
 *****************************************************************************/
#elif defined(__linux__)

long MACAddressUtility::GetMACAddressLinux(std::vector<MacAddress>& result)
{
  int sd = socket(AF_INET, SOCK_DGRAM, 0);
  if (sd == -1)
  {
    messerr("Cannot access to network interfaces");
    return -1;
  }

  std::string line;
  std::ifstream proc_net_dev("/proc/net/dev");
  if (!proc_net_dev.is_open())
  {
    messerr("Cannot parse /proc/net/dev file");
    return -1;
  }
  int nskip_line = 2;
  while (nskip_line > 0)
  {
    getline(proc_net_dev, line);
    nskip_line--;
  }
  while (getline(proc_net_dev, line))
  {
    std::stringstream sstr(line);
    std::string name;
    std::string mac_address;
    getline(sstr, name, ':');

    // Trim
    name = name.substr (name.find_first_not_of (" " ));
    name = name.substr (0, name.find_last_not_of (" " ) + 1);

    if (GetMACAddressLinuxByName(sd, name, mac_address) == 0)
    {
      result.push_back(MacAddress(mac_address));
    }
  }

  if (result.size() <= 0)
  {
    messerr("No network interface found");
    return -1;
  }

  return 0;
}

long MACAddressUtility::GetMACAddressLinuxByName(int sd, 
                                                 const std::string& name, 
                                                 std::string& mac_address) 
{
  struct ifreq ifreq;
  gslStrcpy(ifreq.ifr_name, name.c_str());
  if (ioctl (sd, SIOCGIFHWADDR, &ifreq) < 0 ||
      ioctl (sd, SIOCGIFFLAGS , &ifreq) < 0)
  {
    messerr("Cannot retrieve network interface properties for %s", ifreq.ifr_name);
    return -1;
  }

  // keep only broadcast interfaces
  if (ifreq.ifr_flags & IFF_BROADCAST)
  {
    mac_address = std::string((char*)ifreq.ifr_hwaddr.sa_data, 6);
    return 0;
  }
  return -1;
}

const char * MACAddressUtility::flags(int sd, const char * name)
{
  static char buf[1024];

  static struct ifreq ifreq;
  gslStrcpy(ifreq.ifr_name, name);

  int r = ioctl(sd, SIOCGIFFLAGS, (char *)&ifreq);
  if (r != 0) messageAbort("Error in MACAdressUtility");

  int l = 0;
#define FLAG(b) if(ifreq.ifr_flags & b) l += snprintf(buf + l, sizeof(buf) - l, #b " ")
  FLAG(IFF_UP);
  FLAG(IFF_BROADCAST);
  FLAG(IFF_DEBUG);
  FLAG(IFF_LOOPBACK);
  FLAG(IFF_POINTOPOINT);
  FLAG(IFF_RUNNING);
  FLAG(IFF_NOARP);
  FLAG(IFF_PROMISC);
  FLAG(IFF_NOTRAILERS);
  FLAG(IFF_ALLMULTI);
  FLAG(IFF_MASTER);
  FLAG(IFF_SLAVE);
  FLAG(IFF_MULTICAST);
  FLAG(IFF_PORTSEL);
  FLAG(IFF_AUTOMEDIA);
  FLAG(IFF_DYNAMIC);
#undef FLAG
  return buf;
}

/****************************************************************************
 ** Version for MAC
 *****************************************************************************/
#elif defined(__APPLE__)
 
static kern_return_t FindEthernetInterfaces(io_iterator_t *matchingServices)
{
  kern_return_t           kernResult;
  CFMutableDictionaryRef  matchingDict;
  CFMutableDictionaryRef  propertyMatchDict;
  
  matchingDict = IOServiceMatching(kIOEthernetInterfaceClass);
  
  if (NULL != matchingDict)
  {
    propertyMatchDict = CFDictionaryCreateMutable(kCFAllocatorDefault, 0,
                                                  &kCFTypeDictionaryKeyCallBacks,
                                                  &kCFTypeDictionaryValueCallBacks);
    
    if (NULL != propertyMatchDict)
    {
      CFDictionarySetValue(propertyMatchDict, CFSTR(kIOPrimaryInterface), kCFBooleanTrue);
      CFDictionarySetValue(matchingDict, CFSTR(kIOPropertyMatchKey), propertyMatchDict);
      CFRelease(propertyMatchDict);
    }
  }
  kernResult = IOServiceGetMatchingServices(kIOMasterPortDefault, matchingDict, matchingServices);
  return kernResult;
}

static kern_return_t GetMACAddress(io_iterator_t intfIterator,
                                   UInt8 *MACAddress,
                                   UInt8 bufferSize)
{
  io_object_t     intfService;
  io_object_t     controllerService;
  kern_return_t   kernResult = KERN_FAILURE;
  
  if (bufferSize < kIOEthernetAddressSize) 
  {
    return kernResult;
  }
 
  bzero(MACAddress, bufferSize);
  
  while ((intfService = IOIteratorNext(intfIterator)))
  {
    CFTypeRef   MACAddressAsCFData;       
    
    // IONetworkControllers can't be found directly by the IOServiceGetMatchingServices call,
    // since they are hardware nubs and do not participate in driver matching. In other words,
    // registerService() is never called on them. So we've found the IONetworkInterface and will
    // get its parent controller by asking for it specifically.
    
    // IORegistryEntryGetParentEntry retains the returned object,
    // so release it when we're done with it.
    kernResult = IORegistryEntryGetParentEntry(intfService,
                                               kIOServicePlane,
                                               &controllerService);
    
    if (KERN_SUCCESS != kernResult)
    {
      message("IORegistryEntryGetParentEntry returned 0x%08x\n", kernResult);
    }
    else
    {
      // Retrieve the MAC address property from the I/O Registry in the form of a CFData
      MACAddressAsCFData = IORegistryEntryCreateCFProperty(controllerService,
                                                           CFSTR(kIOMACAddress),
                                                           kCFAllocatorDefault,
                                                           0);
      if (MACAddressAsCFData)
      {
        //CFShow(MACAddressAsCFData); // for display purposes only; output goes to stderr

        // Get the raw bytes of the MAC address from the CFData
        CFDataGetBytes((CFDataRef)MACAddressAsCFData, CFRangeMake(0, kIOEthernetAddressSize), MACAddress);
        CFRelease(MACAddressAsCFData);
      }
      
      // Done with the parent Ethernet controller object so we release it.
      (void) IOObjectRelease(controllerService);
    }
    
    // Done with the Ethernet interface object so we release it.
    (void) IOObjectRelease(intfService);
  }
  return kernResult;
}

long MACAddressUtility::GetMACAddressMAC(std::vector<MacAddress>& result)
{
  unsigned char   one_result[10];
  io_iterator_t   intfIterator;
  kern_return_t   kernResult = KERN_FAILURE;
  do
  {
    kernResult = ::FindEthernetInterfaces(&intfIterator);
    if (KERN_SUCCESS != kernResult) break;
    kernResult = ::GetMACAddress(intfIterator, (UInt8*)one_result, 6);
    MacAddress mac_address((unsigned char*)one_result, 6);
    result.push_back(mac_address);
  }
  while(false);
  
  (void) IOObjectRelease(intfIterator);

  return 0;
}

#endif

