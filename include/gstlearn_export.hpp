#pragma once

// Do not use cmake GenerateExportHeader which is not adapted to double targets (static/shared)
#if defined(_WIN32) || defined (__CYGWIN__) // Windows
#  ifdef GSTLEARN_BUILD_SHARED // Shared export
     // We are compiling the shared library
#    define GSTLEARN_EXPORT __declspec(dllexport)
#  else
#    ifdef GSTLEARN_IMPORT_SHARED // Shared import
       // We are importing the shared library")
#      define GSTLEARN_EXPORT __declspec(dllimport)
#    else 
#      ifdef GSTLEARN_BUILD_STATIC // Static export
         // We are compiling the static library
#      else // Static import
         // We are importing the static library
#      endif
#      define GSTLEARN_EXPORT
#    endif
#  endif
#else // Linux
#  if defined(GSTLEARN_BUILD_SHARED) || defined (GSTLEARN_IMPORT_SHARED)
#    define GSTLEARN_EXPORT __attribute__((visibility("default")))
#  else
#    define GSTLEARN_EXPORT
#  endif
// else TODO MacOS
#endif

#ifdef SWIG
#  undef GSTLEARN_EXPORT
#  define GSTLEARN_EXPORT
#endif
