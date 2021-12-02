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
#pragma once

#ifndef WINVER
  #define WINVER 0x0600
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

#include "geoslib_f_private.h"
#include "geoslib_d.h"
#include "geoslib_f.h"
#include "version.h"

