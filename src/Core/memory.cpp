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
#include "Basic/Memory.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/String.hpp"
#include "Core/Keypair.hpp"

#include <cstring>

/*! \cond */
#define SHIFT() ((MEMORY_DEBUG) ? (Id)sizeof(Id) : 0)

namespace gstlrn
{
typedef struct
{
  String call_file;
  Id call_line;
  size_t size;
  void* ptr;
} MemChunk;

/*! \endcond */

static Id MEMORY_LEAK   = 0;
static Id MEMORY_DEBUG  = 0;
static Id MEMORY_TOTAL  = 0;
static Id MEMORY_MAX    = 0;
static Id MEMORY_MIN_PT = 1000000;
static std::vector<MemChunk> MemLeak;

/****************************************************************************/
/*!
 ** Update the Memory management
 **
 ** \param[in]  size    Size of the allocation/deallocation
 **
 *****************************************************************************/
static void st_mem_update(Id size)
{
  MEMORY_TOTAL += size;
  if (MEMORY_TOTAL > MEMORY_MAX) MEMORY_MAX = MEMORY_TOTAL;
}

/****************************************************************************/
/*!
 ** Add a Memory Chunk
 **
 ** \param[in]  call_file     Name of the calling file
 ** \param[in]  call_line     Line of the calling function
 ** \param[in]  size          Dimension of the Chunk
 ** \param[in]  ptr           Address of the Chunk
 **
 *****************************************************************************/
static void st_memory_leak_add(const char* call_file,
                               Id call_line,
                               size_t size,
                               void* ptr)
{
  if (!MEMORY_LEAK) return;

  Id nb_mem_chunk = MemLeak.size();
  MemLeak.resize(nb_mem_chunk + 1);
  MemChunk& chunk = MemLeak[nb_mem_chunk];
  gslStrcpy2(chunk.call_file, call_file);
  chunk.call_line = call_line;
  chunk.size      = size;
  chunk.ptr       = ptr;
}

/****************************************************************************/
/*!
 ** Delete a Memory Chunk
 **
 ** \param[in]  call_file Name of the calling file
 ** \param[in]  call_line Line in the calling file
 ** \param[in]  ptr       Address of the Chunk to be freed
 **
 *****************************************************************************/
static void st_memory_leak_delete(const char* call_file,
                                  size_t call_line,
                                  void* ptr)
{
  if (!MEMORY_LEAK) return;

  // Look for the chunk to be freed

  Id found = -1;
  for (Id i = 0, n = MemLeak.size(); i < n && found < 0; i++)
  {
    MemChunk& chunk = MemLeak[i];
    if (chunk.ptr == ptr) found = i;
  }

  // The Chunk to be freed does not seem to be allocated

  if (found < 0)
  {
    messerr("Attempt to free a Chunk which seems not to be allocated (called from %s : %d)",
            call_file, call_line);
    return;
  }

  // Compress the array of Memory chunks

  Id nb_mem_chunk = MemLeak.size();
  MemLeak[found]  = MemLeak[nb_mem_chunk - 1];
  MemLeak.resize(nb_mem_chunk - 1);
}

/****************************************************************************/
/*!
 ** Print a memory debugging message
 **
 ** \param[in]  call_file Name of the calling file
 ** \param[in]  call_line Line in the calling file
 ** \param[in]  format    Output format
 ** \param[in]  oper      Sign of the operation (1 for allocation; -1 for free)
 ** \param[in]  size      Number of bytes treated
 **
 ** \remarks The printout is performed only for memory chunks whose size is
 ** \remarks larger than MEMORY_MIN_PT
 ** \remarks It can be modified using keypair with keyword "Minimum_Debug_Size"
 **
 *****************************************************************************/
static void st_mem_message(const char* call_file,
                           size_t call_line,
                           const char* format,
                           Id oper,
                           Id size)
{
  Id minsize;

  minsize = static_cast<Id>(get_keypone("Minimum_Debug_Size", MEMORY_MIN_PT));

  if (MEMORY_DEBUG > 1 && size > minsize)
  {
    if (oper > 0)
      message("%s (%15s : %5d): +%5d Nbytes - Still allocated (%6d)\n", format,
              call_file, call_line, size, MEMORY_TOTAL);
    else
      message("%s (%15s : %5d): -%5d Nbytes - Still allocated (%6d)\n", format,
              call_file, call_line, size, MEMORY_TOTAL);
  }
}

/****************************************************************************/
/*!
 ** Core deallocation
 **
 ** \return  Pointer to the freed array
 **
 ** \param[in]  call_file Name of the calling file
 ** \param[in]  call_line Line in the calling file
 ** \param[in]  tab       Array to be freed
 **
 *****************************************************************************/
char* mem_free_(const char* call_file,
                size_t call_line,
                char* tab)
{
  Id size_eff;
  char* tab_aux;

  if (tab == nullptr) return (tab);
  tab_aux = &tab[-SHIFT()];

  if (MEMORY_DEBUG)
  {
    (void)memcpy(reinterpret_cast<char*>(&size_eff), tab_aux, sizeof(Id));
    st_mem_update(-size_eff);
    st_mem_message(call_file, call_line, "De-allocation", -1, size_eff);
  }
  if (MEMORY_LEAK)
  {
    st_memory_leak_delete(call_file, call_line, tab_aux);
  }

  free(tab_aux);
  tab = NULL;
  return (tab);
}

/****************************************************************************/
/*!
 ** Core allocation routine
 **
 ** \return  Pointer to the array to be allocated
 **
 ** \param[in]  call_file  Name of the calling file
 ** \param[in]  call_line  Line in the calling file
 ** \param[in]  size       Number of bytes
 ** \param[in]  flag_fatal Error status (1 = the program stops)
 **
 *****************************************************************************/
char* mem_alloc_(const char* call_file,
                 size_t call_line,
                 Id size,
                 Id flag_fatal)
{
  Id size_eff;
  char *tab, *tab_aux;

  tab = tab_aux = nullptr;
  if (size <= 0) return (NULL);
  size_eff = size;
  size     = size_eff + SHIFT();

  tab_aux = (char*)malloc(size);
  if (tab_aux == nullptr)
  {
    mem_error(size_eff);
    if (flag_fatal) messageAbort("Fatal error");
    return (NULL);
  }

  if (MEMORY_DEBUG)
  {
    (void)memcpy((char*)tab_aux, (char*)&size_eff, sizeof(Id));
    st_mem_update(size_eff);
    st_mem_message(call_file, call_line, "Allocation   ", +1, size_eff);
  }
  if (MEMORY_LEAK)
  {
    st_memory_leak_add(call_file, call_line, size, tab_aux);
  }

  tab = &tab_aux[SHIFT()];
  return (tab);
}

/****************************************************************************/
/*!
 ** Core routine for allocating and copying
 **
 ** \return  Pointer to the newly allocated (and copied) array
 **
 ** \param[in]  call_file  Name of the calling file
 ** \param[in]  call_line  Line in the calling file
 ** \param[in]  tabin      Array to be copied
 ** \param[in]  size       Number of bytes
 ** \param[in]  flag_fatal Error status (1 = the program stops)
 **
 *****************************************************************************/
char* mem_copy_(const char* call_file,
                size_t call_line,
                char* tabin,
                Id size,
                Id flag_fatal)
{
  Id size_eff;
  char *tab, *tab_aux;

  tab = tab_aux = nullptr;
  if (size <= 0) return (NULL);
  size_eff = size;
  size     = size_eff + SHIFT();

  tab_aux = (char*)malloc(size);
  if (tab_aux == nullptr)
  {
    mem_error(size_eff);
    if (flag_fatal) messageAbort("Fatal error");
    return (NULL);
  }

  if (MEMORY_DEBUG)
  {
    (void)memcpy((char*)tab_aux, (char*)&size_eff, sizeof(Id));
    st_mem_update(size_eff);
    st_mem_message(call_file, call_line, "Allocation   ", +1, size_eff);
  }
  if (MEMORY_LEAK)
  {
    st_memory_leak_add(call_file, call_line, size, tab_aux);
  }

  tab = &tab_aux[SHIFT()];

  /* Copy the input array */

  (void)memcpy(tab, tabin, size);

  return (tab);
}

/****************************************************************************/
/*!
 ** Core allocation routine
 **
 ** \return  Pointer to the array to be allocated
 **
 ** \param[in]  call_file  Name of the calling file
 ** \param[in]  call_line  Line in the calling file
 ** \param[in]  size       Number of elements
 ** \param[in]  size_elem  Number of bytes per element
 ** \param[in]  flag_fatal Error status (1 = the program stops)
 **
 *****************************************************************************/
char* mem_calloc_(const char* call_file,
                  size_t call_line,
                  Id size,
                  Id size_elem,
                  Id flag_fatal)
{
  Id size_eff;
  char *tab, *tab_aux;

  tab = tab_aux = nullptr;
  if (size <= 0) return (NULL);
  size_eff = size * size_elem;
  size     = size_eff + SHIFT();

  tab_aux = (char*)calloc(size_elem, size);
  if (tab_aux == nullptr)
  {
    mem_error(size_eff);
    if (flag_fatal) messageAbort("Fatal error");
    return (NULL);
  }

  if (MEMORY_DEBUG)
  {
    (void)memcpy((char*)tab_aux, (char*)&size_eff, sizeof(Id));
    st_mem_update(size_eff);
    st_mem_message(call_file, call_line, "Allocation   ", +1, size_eff);
  }
  if (MEMORY_LEAK)
  {
    st_memory_leak_add(call_file, call_line, size, tab_aux);
  }

  tab = &tab_aux[SHIFT()];
  return (tab);
}

/****************************************************************************/
/*!
 * Core re-allocation routine
 *
 * \return  Pointer to the array to be re_allocated
 *
 * \param[in]  call_file  Name of the calling file
 * \param[in]  call_line  Line in the calling file
 * \param[in]  tab        Array to be reallocated
 * \param[in]  size       New number of bytes
 * \param[in]  flag_fatal Error status (1 = the program stops)
 *
 *****************************************************************************/
char* mem_realloc_(const char* call_file,
                   size_t call_line,
                   char* tab,
                   Id size,
                   Id flag_fatal)
{
  Id size_eff, size_old;
  char* tab_aux;

  size_eff = size;
  size     = size_eff + SHIFT();

  if (size_eff > 0)
  {
    // The new dimension is positive

    if (tab == nullptr)
    {

      // The memory chunk does not already exist

      tab_aux = (char*)malloc(size);
      if (MEMORY_DEBUG)
      {
        (void)memcpy((char*)tab_aux, (char*)&size_eff, sizeof(Id));
        st_mem_update(size_eff);
        st_mem_message(call_file, call_line, "Allocation   ", +1, size_eff);
      }
      if (MEMORY_LEAK)
      {
        st_memory_leak_add(call_file, call_line, size, tab_aux);
      }
    }
    else
    {

      // The memory chunk already exists

      tab_aux = &tab[-SHIFT()];
      if (MEMORY_DEBUG)
      {
        (void)memcpy((char*)&size_old, (char*)tab_aux, sizeof(Id));
        st_mem_update(-size_old);
        st_mem_message(call_file, call_line, "Re_allocation", -1, size_old);
      }
      if (MEMORY_LEAK)
      {
        st_memory_leak_delete(call_file, call_line, tab_aux);
      }
      auto* placeholder = realloc(tab_aux, size);
      tab_aux           = (char*)placeholder;
      if (MEMORY_DEBUG)
      {
        (void)memcpy((char*)tab_aux, (char*)&size_eff, sizeof(Id));
        st_mem_update(size_eff);
        st_mem_message(call_file, call_line, "Re-allocation", +1, size_eff);
      }
      if (MEMORY_LEAK)
      {
        st_memory_leak_add(call_file, call_line, size, tab_aux);
      }
    }

    if (tab_aux == nullptr)
    {
      mem_error(size_eff);
      if (flag_fatal) messageAbort("Fatal error");
    }
    tab = &tab_aux[SHIFT()];
  }
  else
  {

    // The new dimension is zero

    if (tab != nullptr)
    {
      tab_aux = &tab[-SHIFT()];
      if (MEMORY_DEBUG)
      {
        (void)memcpy((char*)&size_old, (char*)tab_aux, sizeof(Id));
        st_mem_update(-size_old);
        st_mem_message(call_file, call_line, "Re-allocation", -1, size_old);
      }
      if (MEMORY_LEAK)
      {
        st_memory_leak_delete(call_file, call_line, tab_aux);
      }
      free(tab_aux);
      tab = NULL;
    }
  }

  return (tab);
}

/****************************************************************************/
/*!
 ** Core deallocation of an array of pointers
 **
 ** \return  Pointer to the freed array
 **
 ** \param[in]  tab   array of pointers to be freed
 ** \param[in]  nvar  Number of elements in the array
 **
 *****************************************************************************/
double** mem_tab_free(double** tab, Id nvar)
{
  Id ivar;

  if (tab == nullptr) return (tab);
  for (ivar = 0; ivar < nvar; ivar++)
    tab[ivar] = (double*)mem_free((char*)tab[ivar]);
  tab = (double**)mem_free((char*)tab);
  return (tab);
}

/****************************************************************************/
/*!
 ** Core allocation of an array of double
 **
 ** \return  Pointer to the array of pointers to be allocated
 **
 ** \param[in]  nvar        number of elements in the array
 ** \param[in]  size        number of double values
 ** \param[in]  flag_fatal  error status (1 = the program stops)
 **
 *****************************************************************************/
double** mem_tab_alloc(Id nvar, Id size, Id flag_fatal)
{
  double** tab;
  Id ivar, i;

  /* Allocate the array */

  tab = (double**)mem_alloc(sizeof(double*) * nvar, flag_fatal);
  if (tab == nullptr) return (tab);
  for (ivar = 0; ivar < nvar; ivar++)
    tab[ivar] = nullptr;

  for (ivar = 0; ivar < nvar; ivar++)
  {
    tab[ivar] = (double*)mem_alloc(sizeof(double) * size, flag_fatal);
    if (tab[ivar] == nullptr)
    {
      tab = mem_tab_free(tab, nvar);
      return (tab);
    }
    for (i = 0; i < size; i++)
      tab[ivar][i] = 0.;
  }
  return (tab);
}
} // namespace gstlrn