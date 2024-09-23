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

#include "gstlearn_export.hpp"

#define mem_free(tab)          mem_free_(__FILE__, __LINE__, tab)
#define mem_alloc(a, b)        mem_alloc_(__FILE__, __LINE__, a, b)
#define mem_calloc(a, b, c)    mem_calloc_(__FILE__, __LINE__, a, b, c)
#define mem_realloc(tab, a, b) mem_realloc_(__FILE__, __LINE__, tab, a, b)
#define mem_copy(tab, a, b)    mem_copy_(__FILE__, __LINE__, tab, a, b)

GSTLEARN_EXPORT void memory_leak_set(int flag);
GSTLEARN_EXPORT void memory_leak_reset(void);
GSTLEARN_EXPORT void memory_leak_report(void);
GSTLEARN_EXPORT char* mem_alloc_(const char* call_file,
                                 unsigned int call_line,
                                 int size,
                                 int flag_fatal);
GSTLEARN_EXPORT char* mem_calloc_(const char* call_file,
                                  unsigned int call_line,
                                  int size_t,
                                  int size,
                                  int flag_fatal);
GSTLEARN_EXPORT char* mem_realloc_(const char* call_file,
                                   unsigned int call_line,
                                   char* tab,
                                   int size,
                                   int flag_fatal);
GSTLEARN_EXPORT char* mem_copy_(const char* call_file,
                                unsigned int call_line,
                                char* tabin,
                                int size,
                                int flag_fatal);
GSTLEARN_EXPORT char*
mem_free_(const char* call_file, unsigned int call_line, char* tab);
GSTLEARN_EXPORT void mem_debug_set(int flag);
GSTLEARN_EXPORT void memory_status(const char* title);
GSTLEARN_EXPORT double** mem_tab_free(double** tab, int nvar);
GSTLEARN_EXPORT double** mem_tab_alloc(int nvar, int size, int flag_fatal);
GSTLEARN_EXPORT unsigned long long getTotalSystemMemory();
