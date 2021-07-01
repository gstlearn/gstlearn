#ifndef _MACADDRESS_UTILITY_H
#define _MACADDRESS_UTILITY_H

#include "geoslib_e.h"

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

