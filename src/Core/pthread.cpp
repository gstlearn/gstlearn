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
#include "geoslib_e.h"

/*! \cond */

/*! \endcond */

/*****************************************************************************/
/*!
**  Join a thread
**
*****************************************************************************/
static void st_th_join(Th_Rank th_rank)
{
  pthread_join(th_rank,NULL);
  return;
}

/*****************************************************************************/
/*!
**  Delete the current thread
**
*****************************************************************************/
GEOSLIB_API Th_Rank th_delete()
{
  pthread_exit(NULL);
  return(0);
}

/*****************************************************************************/
/*!
**  Create a thread
**
** \return  Th_Rank created
**
** \param[in]  func    Function to be launched
** \param[in]  arg     Argument passed to function
**
*****************************************************************************/
GEOSLIB_API Th_Rank th_create(void *(*func)(void *),
                              void *arg)
{
  Th_Rank th_rank;
  int error;

  error = pthread_create(&th_rank,NULL,func,arg);
  if (error) messageAbort("%s",strerror(error));
  return(th_rank);
}

/*****************************************************************************/
/*!
**  Wait for all threads to be performed
**
** \return Error return code
**
** \param[in]  number   Number of threads to be freed
** \param[in]  th_ranks Array of threads to be freed (Dimension: number)
**
*****************************************************************************/
GEOSLIB_API int th_wait_for_all(int      number,
                                Th_Rank *th_ranks)
{
  for (int i=0; i<number; i++)
    st_th_join(th_ranks[i]);
  return(0);
}

/*****************************************************************************/
/*!
**  Initialize a mutex
**
** \param[in]  th_mutex Th_Mutex structure
**
*****************************************************************************/
GEOSLIB_API void th_mutex_init(Th_Mutex *th_mutex)
{
  int error;
  error = pthread_mutex_init(th_mutex,NULL);
  if (error) messageAbort("%s",strerror(error));
}

/*****************************************************************************/
/*!
**  Switch ON the mutex
**
** \param[in]  th_mutex Th_Mutex structure
**
*****************************************************************************/
GEOSLIB_API void th_mutex_lock(Th_Mutex *th_mutex)
{
  int error;

  error = pthread_mutex_lock(th_mutex);
  if (error) messageAbort("%s",strerror(error));
}

/*****************************************************************************/
/*!
**  Switch OFF the mutex
**
** \param[in]  th_mutex Th_Mutex structure
**
*****************************************************************************/
GEOSLIB_API void th_mutex_unlock(Th_Mutex *th_mutex)
{
  int error;

  error = pthread_mutex_unlock(th_mutex);
  if (error) messageAbort("%s",strerror(error));
}

/*****************************************************************************/
/*!
**  Send a Condition Signal
**
** \param[in]  th_cond  Th_Cond structure
**
*****************************************************************************/
GEOSLIB_API void th_cond_signal(Th_Cond *th_cond)
{
  int error;

  error = pthread_cond_signal(th_cond);
  if (error) messageAbort("%s",strerror(error));
}

/*****************************************************************************/
/*!
**  Initialize a Condition Signal
**
** \param[in]  th_cond  Th_Cond structure
**
*****************************************************************************/
GEOSLIB_API void th_cond_init(Th_Cond *th_cond)
{
  int error;

  error = pthread_cond_init(th_cond,NULL);
  if (error) messageAbort("%s",strerror(error));
}

/*****************************************************************************/
/*!
**  Set the current thread in a waiting status for a Condition Signal
**
** \param[in]  th_cond  Th_Cond structure
** \param[in]  th_mutex Th_Mutex structure
**
*****************************************************************************/
GEOSLIB_API void th_cond_wait(Th_Cond  *th_cond,
                              Th_Mutex *th_mutex)
{
  int error;

  error = pthread_cond_wait(th_cond,th_mutex);
  if (error) messageAbort("%s",strerror(error));
}

/*****************************************************************************/
/*!
**  Release all threads waiting for this Condition
**
** \param[in]  th_cond  Th_Cond structure
**
*****************************************************************************/
GEOSLIB_API void th_cond_release_all(Th_Cond *th_cond)
{
  int error;

  error = pthread_cond_broadcast(th_cond);
  if (error) messageAbort("%s",strerror(error));
}

/*****************************************************************************/
/*!
**  Get the number of Cores available
**
** \return  Number of cores available (or 0)
**
*****************************************************************************/
GEOSLIB_API int ThC_get_ncores(void)
{
  int ncores;

  unsigned concurentThreadsSupported = std::thread::hardware_concurrency();
  ncores = concurentThreadsSupported;
  return(ncores);
}

/*****************************************************************************/
/*!
**  Create a Pool of threads
**
** \param[in] tp          Pointer to the Thread_pool container
** \param[in] futures     std::vector<std::future<void>> structure
** \param[in] nTasks      Number of tasks to be launched
** \param[in] func        Function to be launched per task
**                        (Its prototype should be honored)
** \param[in] arg         Set of arguments passed to the function 'func'
**
*****************************************************************************/
GEOSLIB_API void ThC_create(ctpl::thread_pool &tp,
                            std::vector<std::future <void>> &futures,
                            int   nTasks,
                            void *(*func)(int, int, void *),
                            void *arg)
{
  for (int i=0; i<nTasks; i++)
    futures.push_back(tp.push([arg,i,func](int id){ func(id, i, arg); }));
}

/*****************************************************************************/
/*!
**  Wait for all threads of a Thread_pool to be finished
**
** \param[in]  futures  std::vector<std::future<void>> structure
**
*****************************************************************************/
GEOSLIB_API void ThC_wait_for_all(std::vector<std::future<void>> & futures)

{
  int nTasks = futures.size();

  for (int i=0; i<nTasks; i++)
    futures[i].get();

  return;
}
