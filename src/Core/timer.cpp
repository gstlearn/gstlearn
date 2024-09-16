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
#include "Basic/AStringable.hpp"
#include "Basic/String.hpp"
#include "Basic/Timer.hpp"

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

/*! \endcond */

static int NB_TIME_CHUNK = 0;
static int TIME_FOCUS = -1;
static hrc::time_point TIME_CURRENT = hrc::time_point();
static TimeChunk **TimeStat = NULL;

/****************************************************************************/
/*!
 ** Initialize the Timer
 **
 *****************************************************************************/
void time_start(void)
{
  TIME_CURRENT = hrc::now();
}

/****************************************************************************/
/*!
 ** Close the current Time Chunk
 **
 *****************************************************************************/
static void st_time_chunk_close(void)
{
  TimeChunk *tchunk;

  if (TIME_FOCUS < 0) return;

  // Update the current Time Chunk

  tchunk = TimeStat[TIME_FOCUS];
  hrc::time_point newtime = hrc::now();
  ms difftime = std::chrono::duration_cast<ms>(newtime - TIME_CURRENT);
  tchunk->msec += (int) difftime.count();

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

    auto* placeholder =
      realloc((char*)TimeStat, (NB_TIME_CHUNK + 1) * sizeof(TimeChunk*));
    TimeStat = (TimeChunk**)placeholder;
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
    message("Time %s : %d ms (%d calls)\n", tchunk->call_name, tchunk->msec,
            tchunk->ncalls);
  }
}

