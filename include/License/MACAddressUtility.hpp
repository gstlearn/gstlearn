#ifndef _MACADDRESS_UTILITY_H
#define _MACADDRESS_UTILITY_H

#include <vector>
#include <string>
#include <string.h>

#if defined(_WIN32) || defined(_WIN64)
  #include <winsock2.h>
  // Link with Iphlpapi.lib
  #include <iphlpapi.h>
  #if defined(_MSC_VER)
    #pragma comment(lib, "IPHLPAPI.lib")
  #endif
  #define WORKING_BUFFER_SIZE 15000
  #define MAX_TRIES 3
  #define MALLOC(x) HeapAlloc(GetProcessHeap(), 0, (x))
  #define FREE(x) HeapFree(GetProcessHeap(), 0, (x))
  #include <windows.h>
#elif defined(__linux__)
  #include <sys/ioctl.h>
  #include <sys/socket.h>
  #include <netinet/in.h>
  #include <linux/if.h>
  #include <arpa/inet.h>
  #include <unistd.h>
  #define MAX_IFS 64
#elif defined(__APPLE__)
  #include <CoreFoundation/CoreFoundation.h>
  #include <IOKit/IOKitLib.h>
  #include <IOKit/network/IOEthernetInterface.h>
  #include <IOKit/network/IONetworkInterface.h>
  #include <IOKit/network/IOEthernetController.h>
#endif


class MacAddress {
public:
  MacAddress(unsigned char* data, int size) : _data(NULL), _size(size) {
    _data = new unsigned char[_size]();
    memcpy(_data, data, _size);
  }
  MacAddress(std::string addr) : _data(NULL), _size(addr.size()) {
    _data = new unsigned char[_size]();
    memcpy(_data, addr.c_str(), _size);
  }
  MacAddress(const MacAddress& _from) {
    _size = _from._size;
    _data = new unsigned char[_size]();
    memcpy(_data, _from._data, _size);
  }
  ~MacAddress() {
    delete _data;
  }
  const unsigned char* data() const { return _data; }
  int                  size() const { return _size; }

private:
  MacAddress();
  unsigned char* _data;
  int _size;
};

class MACAddressUtility
{
public:
  static long GetMACAddress(std::vector<MacAddress>& result);
private:
#if defined(_WIN32) || defined(_WIN64)
  static long GetMACAddressMSW(std::vector<MacAddress>& result);
#elif defined(__linux__)
  static long GetMACAddressLinux(std::vector<MacAddress>& result);
  static const char * flags(int sd, const char * name);
  static long GetMACAddressLinuxByName(int sd, const std::string& name, 
                                       std::string& mac_address);
#elif defined(__APPLE__)
  static kern_return_t FindEthernetInterfaces(io_iterator_t *matchingServices);
  static kern_return_t GetMACAddress(io_iterator_t intfIterator,
                                     UInt8 *MACAddress, UInt8 bufferSize);
  static long GetMACAddressMAC(std::vector<MacAddress>& result);
#endif

};

#endif

