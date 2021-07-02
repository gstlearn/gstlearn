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
#ifndef GEOSLIB_E_H
#define GEOSLIB_E_H

#ifndef WINVER
  #define WINVER 0x0600
#endif

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

#ifndef NULL
#define NULL 0      ///< NULL macro
#endif /* NULL */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h> 
#include <stdarg.h>
#include <ctype.h>
#include <sys/types.h>

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <string>
#include <iostream>
#include <list>
#include <map>
#include <set>
#include <sstream>
#include <fstream>
#include <algorithm>

#if defined(_WIN32) || defined(_WIN64)
#include <cstdlib>
#endif

#if !defined(__APPLE__)
#include <malloc.h>
#endif

#include "geoslib_d_private.h"
#include "geoslib_f_private.h"
#include "geoslib_d.h"
#include "geoslib_f.h"
#include "version.h"
#include "ctpl.h"

#include "License/MD5Utility.hpp"
#include "License/LicenseKey.hpp"

#endif
