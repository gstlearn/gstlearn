/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES Paris / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://gstlearn.org                                              */
/* License: BSD 3-clause                                                      */
/*                                                                            */
/******************************************************************************/
#pragma once

// Thanks to https://www.fluentcpp.com/2019/08/30/how-to-disable-a-warning-in-cpp/

#if defined(_MSC_VER)
  #define DISABLE_WARNING_PUSH           __pragma(warning( push ))
  #define DISABLE_WARNING_POP            __pragma(warning( pop ))
  #define DISABLE_WARNING(warningNumber) __pragma(warning( disable : warningNumber ))

  #define DISABLE_WARNING_DECLARATION_MASKED               DISABLE_WARNING(4456)
  #define DISABLE_WARNING_EXPR_COND_ASSIGNMENT             DISABLE_WARNING(4706)
  #define DISABLE_WARNING_COND_EXPR_CONSTANT               DISABLE_WARNING(4127)
  #define DISABLE_WARNING_NOT_EXPORTED_FROM_DLL            DISABLE_WARNING(4251)
  #define DISABLE_WARNING_BASE_NOT_EXPORTED_FROM_DLL       DISABLE_WARNING(4275)
  #define DISABLE_WARNING_UNREFERENCED_FORMAL_PARAMETER    DISABLE_WARNING(4100)
  #define DISABLE_WARNING_UNUSED_BUT_SET_VARIABLE          // TODO : Warning equivalence ?
  #define DISABLE_WARNING_FREE_NON_HEAP_OBJECT             //
  #define DISABLE_WARNING_UNREFERENCED_FUNCTION            DISABLE_WARNING(4505)
  #define DISABLE_WARNING_ARRAY_BOUNDS                     // No equivalence? or maybe https://stackoverflow.com/a/38732408
  #define DISABLE_WARNING_STRIC_OVERFLOW                   // No equivalence?
  #define DISABLE_WARNING_RESTRIC                          // No equivalence?
  // other warnings you want to deactivate...

#elif defined(__GNUC__) || defined(__clang__)
  #if defined(__GNUC__)
    #define DO_PRAGMA(X) _Pragma(#X)
    #define DISABLE_WARNING_PUSH           DO_PRAGMA(GCC diagnostic push)
    #define DISABLE_WARNING_POP            DO_PRAGMA(GCC diagnostic pop)
    #define DISABLE_WARNING(warningName)   DO_PRAGMA(GCC diagnostic ignored #warningName)
  #else // clang
    #define DO_PRAGMA(X) _Pragma(#X)
    #define DISABLE_WARNING_PUSH           DO_PRAGMA(clang diagnostic push)
    #define DISABLE_WARNING_POP            DO_PRAGMA(clang diagnostic pop)
    #define DISABLE_WARNING(warningName)   DO_PRAGMA(clang diagnostic ignored #warningName)
  #endif

  #define DISABLE_WARNING_DECLARATION_MASKED               // TODO : Warning equivalence ?
  #define DISABLE_WARNING_EXPR_COND_ASSIGNMENT             //
  #define DISABLE_WARNING_COND_EXPR_CONSTANT               //
  #define DISABLE_WARNING_NOT_EXPORTED_FROM_DLL            //
  #define DISABLE_WARNING_BASE_NOT_EXPORTED_FROM_DLL       //
  #define DISABLE_WARNING_UNREFERENCED_FORMAL_PARAMETER    DISABLE_WARNING(-Wunused-parameter)
  #define DISABLE_WARNING_UNUSED_BUT_SET_VARIABLE          DISABLE_WARNING(-Wunused-but-set-variable)
  #define DISABLE_WARNING_FREE_NON_HEAP_OBJECT             DISABLE_WARNING(-Wfree-nonheap-object)
  #define DISABLE_WARNING_UNREFERENCED_FUNCTION            DISABLE_WARNING(-Wunused-function)
  #define DISABLE_WARNING_ARRAY_BOUNDS                     DISABLE_WARNING(-Warray-bounds)
  #define DISABLE_WARNING_STRIC_OVERFLOW                   DISABLE_WARNING(-Wstrict-overflow)
  #define DISABLE_WARNING_RESTRIC                          DISABLE_WARNING(-Wrestrict)
  // other warnings you want to deactivate...

#else
  #define DISABLE_WARNING_PUSH
  #define DISABLE_WARNING_POP

  // TODO : Warning equivalence ?
  #define DISABLE_WARNING_DECLARATION_MASKED
  #define DISABLE_WARNING_EXPR_COND_ASSIGNMENT
  #define DISABLE_WARNING_COND_EXPR_CONSTANT
  #define DISABLE_WARNING_NOT_EXPORTED_FROM_DLL
  #define DISABLE_WARNING_BASE_NOT_EXPORTED_FROM_DLL
  #define DISABLE_WARNING_UNREFERENCED_FORMAL_PARAMETER
  #define DISABLE_WARNING_UNUSED_BUT_SET_VARIABLE
  #define DISABLE_WARNING_FREE_NON_HEAP_OBJECT
  #define DISABLE_WARNING_UNREFERENCED_FUNCTION
  #define DISABLE_WARNING_ARRAY_BOUNDS
  #define DISABLE_WARNING_STRIC_OVERFLOW
  #define DISABLE_WARNING_RESTRIC
  // other warnings you want to deactivate...

#endif
