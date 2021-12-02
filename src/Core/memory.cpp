//****************************************************************************//
// COPYRIGHT ARMINES, ALL RIGHTS RESERVED                                     //
//                                                                            //
// THE CONTENT OF THIS WORK CONTAINS CONFIDENTIAL AND PROPRIETARY             //
// INFORMATION OF ARMINES. ANY DUPLICATION, MODIFICATION,                     //
// DISTRIBUTION, OR DISCLOSURE IN ANY FORM, IN WHOLE, OR IN PART, IS STRICTLY //
// PROHIBITED WITHOUT THE PRIOR EXPRESS WRITTEN PERMISSION OF ARMINES         //
//                                                                            //
// TAG_SOURCE_CG                                                              //
//****************************************************************************//
//#include "geoslib_e.h"
#include "geoslib_old_f.h"
#include "Basic/String.hpp"

#include <string.h>

/*! \cond */
#define STORE_NAME_LENGTH 10
#define SHIFT()  ((MEMORY_DEBUG) ? (int) sizeof(int) : 0)

typedef struct
{
  char call_name[STORE_NAME_LENGTH];
  int ncalls;
  int msec;
} TimeChunk;

typedef struct
{
  char call_file[STORE_NAME_LENGTH];
  int call_line;
  size_t size;
  void *ptr;
} MemChunk;

/*! \endcond */

static int MEMORY_LEAK = 0;
static int MEMORY_DEBUG = 0;
static int MEMORY_TOTAL = 0;
static int MEMORY_MAX = 0;
static int MEMORY_MIN_PT = 1000000;
static int NB_MEM_CHUNK = 0;
static int NB_TIME_CHUNK = 0;
static int TIME_FOCUS = -1;
static clock_t TIME_CURRENT = 0;
static MemChunk **MemLeak = NULL;
static TimeChunk **TimeStat = NULL;

/****************************************************************************/
/*!
 ** Initialize the Timer
 **
 *****************************************************************************/
void time_start(void)
{
  TIME_CURRENT = clock();
}

/****************************************************************************/
/*!
 ** Close the current Time Chunk
 **
 *****************************************************************************/
static void st_time_chunk_close(void)
{
  TimeChunk *tchunk;
  clock_t newtime, difftime;

  if (TIME_FOCUS < 0) return;

  // Update the current Time Chunk

  tchunk = TimeStat[TIME_FOCUS];
  newtime = clock();
  difftime = newtime - TIME_CURRENT;
  tchunk->msec += difftime / CLOCKS_PER_SEC;

  TIME_CURRENT = newtime;
  TIME_FOCUS = -1;
}

/****************************************************************************/
/*!
 ** Reset the Time
 **
 *****************************************************************************/
void time_reset(void)
{
  // Close the current Time Chunk (if any)

  st_time_chunk_close();

  for (int i = 0; i < NB_TIME_CHUNK; i++)
    free((char*) TimeStat[i]);
  free((char*) TimeStat);
  TimeStat = NULL;
  NB_TIME_CHUNK = 0;
  TIME_FOCUS = -1;
}

/****************************************************************************/
/*!
 ** Initialize a Time Chunk
 **
 ** \param[in]  call_name       Name of the Chunk
 **
 *****************************************************************************/
void time_chunk_add(const char *call_name)
{
  TimeChunk *tchunk;
  int found;

  // Test if the Time Chunk already exists

  found = -1;
  for (int i = 0; i < NB_TIME_CHUNK; i++)
  {
    tchunk = TimeStat[i];
    if (!strncmp(tchunk->call_name, call_name, STORE_NAME_LENGTH - 1))
      found = i;
  }

  if (found < 0)
  {

    // Initialize the new chunk

    tchunk = (TimeChunk*) malloc(sizeof(TimeChunk));
    if (tchunk == NULL)
    {
      messerr("Memory problem: Timer procedure is interrupter");
      time_reset();
      return;
    }
    (void) gslStrncpy(tchunk->call_name, call_name, STORE_NAME_LENGTH);
    tchunk->call_name[STORE_NAME_LENGTH - 1] = '\0';
    tchunk->ncalls = 0;
    tchunk->msec = 0;

    // Glue the new Time Chunk to the Global array

    TimeStat = (TimeChunk**) realloc((char*) TimeStat,
                                     (NB_TIME_CHUNK + 1) * sizeof(TimeChunk));
    if (TimeStat == NULL)
    {
      messerr("Memory problem: Timer procedure is interrupted");
      time_reset();
      return;
    }
    TimeStat[NB_TIME_CHUNK] = tchunk;
    found = NB_TIME_CHUNK;
    NB_TIME_CHUNK++;
  }
  TimeStat[found]->ncalls++;
  TIME_FOCUS = found;
}

/****************************************************************************/
/*!
 ** Report the Time Stats
 **
 *****************************************************************************/
void time_report(void)
{
  TimeChunk *tchunk;

  if (NB_TIME_CHUNK <= 0) return;

  st_time_chunk_close();

  for (int i = 0; i < NB_TIME_CHUNK; i++)
  {
    tchunk = TimeStat[i];
    message("Time %s : %d s (%d calls)\n", tchunk->call_name, tchunk->msec,
            tchunk->ncalls);
  }
}

/****************************************************************************/
/*!
 ** Update the Memory management
 **
 ** \param[in]  size    Size of the allocation/deallocation
 **
 *****************************************************************************/
static void st_mem_update(int size)
{
  MEMORY_TOTAL += size;
  if (MEMORY_TOTAL > MEMORY_MAX) MEMORY_MAX = MEMORY_TOTAL;
}

/****************************************************************************/
/*!
 ** Reset the Memory Leak processing structure
 **
 *****************************************************************************/
void memory_leak_reset(void)
{
  if (!MEMORY_LEAK) return;

  for (int i = 0; i < NB_MEM_CHUNK; i++)
    free((char*) MemLeak[i]);
  free((char*) MemLeak);
  MemLeak = NULL;
  NB_MEM_CHUNK = 0;
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
static void st_memory_leak_add(const char *call_file,
                               int call_line,
                               size_t size,
                               void *ptr)
{
  MemChunk *chunk;

  if (!MEMORY_LEAK) return;

  // Allocate the new Chunk

  chunk = (MemChunk*) malloc(sizeof(MemChunk));
  if (chunk == NULL)
  {
    messerr("Memory problem: Memory Leak procedure is interrupted");
    memory_leak_reset();
    return;
  }
  (void) gslStrncpy(chunk->call_file, call_file, STORE_NAME_LENGTH);
  chunk->call_file[STORE_NAME_LENGTH - 1] = '\0';
  chunk->call_line = call_line;
  chunk->size = size;
  chunk->ptr = ptr;

  // Glue the new chunk to the Global array

  MemLeak = (MemChunk**) realloc((char*) MemLeak,
                                 (NB_MEM_CHUNK + 1) * sizeof(MemChunk));
  if (MemLeak == NULL)
  {
    messerr("Memory problem: Memory Leak procedure is interrupted");
    memory_leak_reset();
    return;
  }
  MemLeak[NB_MEM_CHUNK] = chunk;
  NB_MEM_CHUNK++;
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
static void st_memory_leak_delete(const char *call_file,
                                  unsigned int call_line,
                                  void *ptr)
{
  MemChunk *chunk;
  int found;

  if (!MEMORY_LEAK) return;

  // Look for the chunk to be freed

  found = -1;
  for (int i = 0; i < NB_MEM_CHUNK && found < 0; i++)
  {
    chunk = MemLeak[i];
    if (chunk->ptr == ptr) found = i;
  }

  // The Chunk to be freed does not seem to be allocated

  if (found < 0)
  {
    messerr("A Chunk seems not to be allocated (called from %s : %d)",
            call_file, call_line);
    return;
  }

  // Free the Chunk

  free((char*) MemLeak[found]);

  // Compress the array of Memory chunks

  MemLeak[found] = MemLeak[NB_MEM_CHUNK - 1];
  MemLeak = (MemChunk**) realloc((char*) MemLeak,
                                 (NB_MEM_CHUNK - 1) * sizeof(MemChunk));
  NB_MEM_CHUNK--;
}

/****************************************************************************/
/*!
 ** Report Memory Leak
 **
 *****************************************************************************/
void memory_leak_report(void)
{
  MemChunk *chunk;
  int total;

  if (!MEMORY_LEAK) return;

  if (NB_MEM_CHUNK <= 0)
  {
    message("No Memory Leak\n");
  }
  else
  {
    total = 0;
    for (int i = 0; i < NB_MEM_CHUNK; i++)
    {
      chunk = MemLeak[i];
      message("Leak %s (line:%d) : %d words\n", chunk->call_file,
              chunk->call_line, chunk->size);
      total += static_cast<int>(chunk->size);
    }
    message("Total leak = %d\n", total);
  }
}

/****************************************************************************/
/*!
 ** Set the status of the memory
 **
 ** \param[in]  flag      Activiation flag
 **
 *****************************************************************************/
void mem_debug_set(int flag)
{
  MEMORY_DEBUG = flag;
}

/****************************************************************************/
/*!
 ** Set the memory leak mechanism
 **
 ** \param[in]  flag      Activation flag
 **
 *****************************************************************************/
void memory_leak_set(int flag)
{
  MEMORY_LEAK = flag;
  if (flag == 1)
    memory_leak_reset();
  else
    memory_leak_report();
}

/****************************************************************************/
/*!
 ** Print the status of the memory
 **
 ** \param[in] title   Title printed when checking memory
 **
 *****************************************************************************/
void memory_status(const char *title)

{
  if (!MEMORY_DEBUG) return;
  if (title == NULL)
    message("Memory currently allocated: %d (Max: %d)\n", MEMORY_TOTAL,
            MEMORY_MAX);
  else
    message("Memory currently allocated in %s: %d (Max: %d)\n", title,
            MEMORY_TOTAL, MEMORY_MAX);
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
static void st_mem_message(const char *call_file,
                           unsigned int call_line,
                           const char *format,
                           int oper,
                           int size)
{
  int minsize;

  minsize = (int) get_keypone("Minimum_Debug_Size", MEMORY_MIN_PT);

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
char* mem_free_(const char *call_file, unsigned int call_line, char *tab)
{
  int size_eff;
  char *tab_aux;

  if (tab == nullptr) return (tab);
  tab_aux = &tab[-SHIFT()];

  if (MEMORY_DEBUG)
  {
    (void) memcpy((char*) &size_eff, (char*) tab_aux, sizeof(int));
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
char* mem_alloc_(const char *call_file,
                 unsigned int call_line,
                 int size,
                 int flag_fatal)
{
  int size_eff;
  char *tab, *tab_aux;

  tab = tab_aux = nullptr;
  if (size <= 0) return (NULL);
  size_eff = size;
  size = size_eff + SHIFT();

  tab_aux = (char*) malloc(size);
  if (tab_aux == NULL)
  {
    mem_error(size_eff);
    if (flag_fatal) messageAbort("Fatal error");
    return (NULL);
  }

  if (MEMORY_DEBUG)
  {
    (void) memcpy((char*) tab_aux, (char*) &size_eff, sizeof(int));
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
char* mem_copy_(const char *call_file,
                unsigned int call_line,
                char *tabin,
                int size,
                int flag_fatal)
{
  int size_eff;
  char *tab, *tab_aux;

  tab = tab_aux = nullptr;
  if (size <= 0) return (NULL);
  size_eff = size;
  size = size_eff + SHIFT();

  tab_aux = (char*) malloc(size);
  if (tab_aux == NULL)
  {
    mem_error(size_eff);
    if (flag_fatal) messageAbort("Fatal error");
    return (NULL);
  }

  if (MEMORY_DEBUG)
  {
    (void) memcpy((char*) tab_aux, (char*) &size_eff, sizeof(int));
    st_mem_update(size_eff);
    st_mem_message(call_file, call_line, "Allocation   ", +1, size_eff);
  }
  if (MEMORY_LEAK)
  {
    st_memory_leak_add(call_file, call_line, size, tab_aux);
  }

  tab = &tab_aux[SHIFT()];

  /* Copy the input array */

  (void) memcpy(tab, tabin, size);

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
char* mem_calloc_(const char *call_file,
                  unsigned int call_line,
                  int size,
                  int size_elem,
                  int flag_fatal)
{
  int size_eff;
  char *tab, *tab_aux;

  tab = tab_aux = nullptr;
  if (size <= 0) return (NULL);
  size_eff = size * size_elem;
  size = size_eff + SHIFT();

  tab_aux = (char*) calloc(size_elem, size);
  if (tab_aux == NULL)
  {
    mem_error(size_eff);
    if (flag_fatal) messageAbort("Fatal error");
    return (NULL);
  }

  if (MEMORY_DEBUG)
  {
    (void) memcpy((char*) tab_aux, (char*) &size_eff, sizeof(int));
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
char* mem_realloc_(const char *call_file,
                   unsigned int call_line,
                   char *tab,
                   int size,
                   int flag_fatal)
{
  int size_eff, size_old;
  char *tab_aux;

  size_eff = size;
  size = size_eff + SHIFT();

  if (size_eff > 0)
  {
    // The new dimension is positive

    if (tab == NULL)
    {

      // The memory chunk does not already exist

      tab_aux = (char*) malloc(size);
      if (MEMORY_DEBUG)
      {
        (void) memcpy((char*) tab_aux, (char*) &size_eff, sizeof(int));
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
        (void) memcpy((char*) &size_old, (char*) tab_aux, sizeof(int));
        st_mem_update(-size_old);
        st_mem_message(call_file, call_line, "Re_allocation", -1, size_old);
      }
      if (MEMORY_LEAK)
      {
        st_memory_leak_delete(call_file, call_line, tab_aux);
      }
      tab_aux = (char*) realloc(tab_aux, size);
      if (MEMORY_DEBUG)
      {
        (void) memcpy((char*) tab_aux, (char*) &size_eff, sizeof(int));
        st_mem_update(size_eff);
        st_mem_message(call_file, call_line, "Re-allocation", +1, size_eff);
      }
      if (MEMORY_LEAK)
      {
        st_memory_leak_add(call_file, call_line, size, tab_aux);
      }
    }

    if (tab_aux == NULL)
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
        (void) memcpy((char*) &size_old, (char*) tab_aux, sizeof(int));
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
double** mem_tab_free(double **tab, int nvar)
{
  int ivar;

  if (tab == nullptr) return (tab);
  for (ivar = 0; ivar < nvar; ivar++)
    tab[ivar] = (double*) mem_free((char* ) tab[ivar]);
  tab = (double**) mem_free((char* ) tab);
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
double** mem_tab_alloc(int nvar, int size, int flag_fatal)
{
  double **tab;
  int ivar, i;

  /* Allocate the array */

  tab = (double**) mem_alloc(sizeof(double*) * nvar, flag_fatal);
  if (tab == nullptr) return (tab);
  for (ivar = 0; ivar < nvar; ivar++)
    tab[ivar] = nullptr;

  for (ivar = 0; ivar < nvar; ivar++)
  {
    tab[ivar] = (double*) mem_alloc(sizeof(double) * size, flag_fatal);
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
