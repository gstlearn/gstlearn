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
#include "Basic/Law.hpp"
#include "geoslib_e.h"
#include <unistd.h>

#if defined (_WIN32) || defined(_WIN64)
#  define psleep(sec) Sleep ((sec) * 1000)
#else
#  define psleep(sec) sleep ((sec))
#endif

#define INITIAL_STOCK     20
#define RELOAD_STOCK      20
#define TOTAL_STOCK      510
#define NB_TOTAL_THREADS  21	// Store + 20 clients
#define VERBOSE            0

typedef struct{
  int avail;			// Available number of items in stock
  int total;			// Total number of items in stock
  int reload;			// Amount by which the stock is reloaded
  int finish;     // Stop when this number of items is reached
  int verbose;		// Verbose flag
  Th_Mutex mutex_stock;
  Th_Cond  cond_stock;
  Th_Cond  cond_clients;
} Shared_Area;

typedef struct{
  int rank;			// Rank: 0 (Store); Starting from 1 for Clients
  int cumul;			// Total number of items retrieved
  Shared_Area *shared;		// Pointer to the shared area
} Thread_Area;
  
/* Initialization of the Shared area */
static void st_init_shared(Shared_Area *shared)
{
  shared->avail        = 0;
  shared->total        = 0;
  shared->reload       = RELOAD_STOCK;
  shared->finish       = TOTAL_STOCK;
  shared->verbose      = VERBOSE;
  th_mutex_init(&shared->mutex_stock);
  th_cond_init (&shared->cond_stock);
  th_cond_init (&shared->cond_clients);
}

/* Initialization of the Threaded area */
static void st_init_thread(Thread_Area *area,
                           Shared_Area *shared,
                           int rank)
{
  area->rank   = rank;
  area->cumul  = 0;
  area->shared = shared;
}

/* Function for the thread of the store */

static void *fn_store(void *p_data)
{
  Thread_Area *area   = (Thread_Area *) p_data;
  Shared_Area *shared = area->shared;
  int newadd;
  
  while (1)
  {
    newadd = MIN(shared->finish, shared->total+shared->reload) - shared->total;
    
    th_mutex_lock(&shared->mutex_stock);
    th_cond_wait (&shared->cond_stock, &shared->mutex_stock);
    
    if (newadd > 0)
    {
      area->cumul   += newadd;
      shared->avail += newadd;
      shared->total += newadd;
      if (shared->verbose)
        message("Load the stock by %d reaching %d (Global = %d/%d)\n", 
                newadd,shared->avail,shared->total,shared->finish);
    }
    th_cond_release_all(&shared->cond_clients);
    th_mutex_unlock(&shared->mutex_stock);
    if (shared->total >= shared->finish) break;
  }
  th_delete();
  return NULL;
}

/* Function for thread of clients */

static void *fn_clients(void *p_data)
{
  Thread_Area *area   = (Thread_Area *) p_data;
  Shared_Area *shared = area->shared;
  law_set_random_seed(13241);

  while (1)
  {
    psleep(law_uniform(0.,3.));
    int retrait = (int) law_uniform(0.,6.);
    if (retrait <= 0) continue;

    th_mutex_lock(&shared->mutex_stock);
    
    if (retrait > shared->avail)
    {
      if (shared->total >= shared->finish)
      {
        retrait = 0;
      }
      else
      {
        th_cond_signal(&shared->cond_stock);
        th_cond_wait  (&shared->cond_clients, &shared->mutex_stock);
        retrait = MIN(retrait, shared->avail);
      }
    }
    
    if (retrait > 0)
    {
      area->cumul   += retrait;
      shared->avail -= retrait;
      if (shared->verbose)
        message("Client %2d withdraws %d (His total = %2d). Stock = %2d\n",
                area->rank, retrait, area->cumul, shared->avail);
    }
    th_mutex_unlock(&shared->mutex_stock);
    if (retrait <= 0) break;
  }
  th_delete();
  return NULL;
}
 
/*********************/
/* Program principal */
/*********************/

int main(int argc, char *argv[])

{
  Shared_Area shared;
  Thread_Area areas[NB_TOTAL_THREADS];
  Th_Rank     th_ranks[NB_TOTAL_THREADS];
  int total;

  /* Connect the Geoslib Library */

  if (setup_license("Demonstration")) goto label_end;

  /* Setup constants */

  debug_reset();
  constant_reset();

  /* Initialize the shared area */

  st_init_shared(&shared);

  /* Initialize the threads */

  for (int i=0; i<NB_TOTAL_THREADS; i++)
    st_init_thread(&areas[i],&shared,i);

  /* Creating the threads */

  for (int i=0; i<NB_TOTAL_THREADS; i++)
    if (i == 0)
      th_ranks[i] = th_create(fn_store,&areas[i]);
    else
      th_ranks[i] = th_create(fn_clients,&areas[i]);
  
  /* Wait for the end of the threads */

  (void) th_wait_for_all(NB_TOTAL_THREADS,th_ranks);

  /* Printout of Statistics */

  message("\nFinal statistics (%d Clients)\n",NB_TOTAL_THREADS-1);
  total = 0;
  for (int i=0; i<NB_TOTAL_THREADS; i++)
  {
    if (i == 0) continue;
    total += areas[i].cumul;
// The next message has been commented out
// in order to get no difference in the regression test
//    message("Client %2d cumulates %2d items\n",i,areas[i].cumul);
  }
  message("Total number of items = %d\n",total);

label_end:
  return(0);
}
