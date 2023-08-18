/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://gstlearn.org                                              */
/* License: BSD 3-clause                                                      */
/*                                                                            */
/******************************************************************************/
#include "Basic/Memory.hpp"

#if defined(WIN32) || defined(_WIN32) || defined(__WIN32) && !defined(__CYGWIN__)
#include <windows.h>

unsigned long long getTotalSystemMemory()
{
  MEMORYSTATUSEX status;
  status.dwLength = sizeof(status);
  if (GlobalMemoryStatusEx(&status) == 0) return 0;
  return status.ullTotalPhys;
}
#elif __APPLE__
#include <sys/types.h>
#include <sys/sysctl.h>
#include <sys/vmmeter.h>

unsigned long long getTotalSystemMemory()
{
  int rc;
  u_int page_size;
  struct vmtotal vmt;
  size_t vmt_size, uint_size;

  vmt_size = sizeof(vmt);
  uint_size = sizeof(page_size);

  rc = sysctlbyname("vm.vmtotal", &vmt, &vmt_size, NULL, 0);
  if (rc < 0) return 0;

  rc = sysctlbyname("vm.stats.vm.v_page_size", &page_size, &uint_size, NULL, 0);
  if (rc < 0) return 0;

  return vmt.t_avm * (u_int64_t)page_size;
}

#else // assume Linux
#include <unistd.h>

unsigned long long getTotalSystemMemory()
{
  long avail_pages = sysconf(_SC_AVPHYS_PAGES);
  long page_size = sysconf(_SC_PAGE_SIZE);
  if (avail_pages < 0 || page_size < 0) return 0;
  return avail_pages * page_size;
}
#endif

