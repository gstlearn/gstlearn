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
/* The C clustering library.                                                  */
/* Copyright (C) 2002 Michiel Jan Laurens de Hoon.                            */
/*                                                                            */
/* This library was written at the Laboratory of DNA Information Analysis,    */
/* Human Genome Center, Institute of Medical Science, University of Tokyo,    */
/* 4-6-1 Shirokanedai, Minato-ku, Tokyo 108-8639, Japan.                      */
/******************************************************************************/
#include "Basic/Law.hpp"
#include "geoslib_e.h"
#include "geoslib_old_f.h"

#define DBL_MAXIMUM 1.e30
#define DATA(iech,ivar)         (data[(iech)    * nvar + (ivar)])
#define DATA1(iech,ivar)        (data1[(iech)   * nvar + (ivar)])
#define DATA2(iech,ivar)        (data2[(iech)   * nvar + (ivar)])
#define CDATA(icl,ivar)         (cdata[(icl)    * nvar + (ivar)])
#define DISTMATRIX(i,j)         (distmatrix[(i) * nech + (j)])

/*****************************************************************************/
/*!
 **  Calculate the distance between two samples
 **
 ** \param[in]  nvar      Number of variables
 ** \param[in]  data1     Array of values for first part
 ** \param[in]  data2     Array of values for second part
 ** \param[in]  index1    Rank of the first sample
 ** \param[in]  index2    Rank of the second sample
 **
 *****************************************************************************/
static double st_distance(int nvar,
                          double *data1,
                          double *data2,
                          int index1,
                          int index2)
{
  double result, term, weight;

  result = weight = 0.;
  for (int ivar = 0; ivar < nvar; ivar++)
  {
    term = DATA1(index1,ivar) - DATA2(index2, ivar);
    result += term * term;
    weight += 1.;
  }
  if (weight > 0) result /= weight;
  return (result);

}

/*****************************************************************************/
/*!
 **  Initial random assignment of samples to clusters
 **
 ** \param[in]  nclusters Number if clusters
 ** \param[in]  nech      Number of samples
 **
 ** \param[out] clusterid Array of cluster number for each sample
 **
 *****************************************************************************/
static void st_randomassign(int nclusters, int nech, int *clusterid)
{
  int i, j;
  int k = 0;
  double p;

  int n = nech - nclusters;
  /* Draw the number of elements in each cluster from a multinomial
   * distribution, reserving ncluster elements to set independently
   * in order to guarantee that none of the clusters are empty.
   */
  for (i = 0; i < nclusters - 1; i++)
  {
    p = 1.0 / (nclusters - i);
    j = law_binomial(n, p);
    n -= j;
    j += k + 1; /* Assign at least one element to cluster i */
    for (; k < j; k++)
      clusterid[k] = i;
  }

  /* Assign the remaining elements to the last cluster */
  for (; k < nech; k++)
    clusterid[k] = i;

  /* Create a random permutation of the cluster assignments */
  for (i = 0; i < nech; i++)
  {
    j = (int) (i + (nech - i) * law_uniform(0., 1.));
    k = clusterid[j];
    clusterid[j] = clusterid[i];
    clusterid[i] = k;
  }
}

/*****************************************************************************/
/*!
 **  Print the list of cluster per sample
 **
 ** \param[in]  nclusters Number if clusters
 ** \param[in]  nech      Number of samples
 ** \param[in]  clusterid Array of cluster number for each sample
 **
 *****************************************************************************/
static void st_printclusterlist(int nclusters, int nech, int *clusterid)
{
  message("Population of %d clusters\n", nclusters);
  for (int i = 0; i < nech; i++)
    message("Sample %3d: cluster %d\n", i + 1, clusterid[i]);
}

/*****************************************************************************/
/*!
 **  Print the number of samples per cluster
 **
 ** \param[in]  nclusters Number if clusters
 ** \param[in]  nech      Number of samples
 ** \param[in]  clusterid Array of cluster number for each sample
 **
 *****************************************************************************/
static void st_printclustercount(int nclusters, int nech, int *clusterid)
{
  int j, count;

  message("Population of %d clusters\n", nclusters);
  for (int i = 0; i < nclusters; i++)
  {
    count = 0;
    for (int k = 0; k < nech; k++)
    {
      j = clusterid[k];
      if (i == j) count++;
    }
    message("Cluster %2d : %d\n", i + 1, count);
  }
}

/*****************************************************************************/
/*!
 **  Find the center of a cluster
 **
 ** \param[in]  nclusters Number if clusters
 ** \param[in]  nvar      Number of samples
 ** \param[in]  nech      Number of variables
 ** \param[in]  data      Array of values
 ** \param[in]  clusterid Array of cluster number for each sample
 **
 ** \param[out] cdata     Array containing the centroids
 ** \param[out] cmask     Array containing the number of sample per cluster
 **
 *****************************************************************************/
static void st_getclustermeans(double *data,
                               int nvar,
                               int nech,
                               int nclusters,
                               int *clusterid,
                               double *cdata,
                               int *cmask)
{
  int i, j, k;

  for (i = 0; i < nclusters; i++)
  {
    cmask[i] = 0;
    for (j = 0; j < nvar; j++)
      CDATA(i,j) = 0.;
  }

  for (k = 0; k < nech; k++)
  {
    i = clusterid[k];
    cmask[i]++;
    for (j = 0; j < nvar; j++)
      CDATA(i,j) += DATA(k, j);
  }

  for (i = 0; i < nclusters; i++)
    for (j = 0; j < nvar; j++)
      if (cmask[i] > 0) CDATA(i,j) /= cmask[i];
}

/*****************************************************************************/
/*!
 **  Find the median of a cluster
 **
 ** \param[in]  nclusters Number if clusters
 ** \param[in]  nvar      Number of samples
 ** \param[in]  nech      Number of variables
 ** \param[in]  data      Array of values
 ** \param[in]  clusterid Array of cluster number for each sample
 **
 ** \param[out] cdata     Array containing the centroids
 ** \param[out] cache     Array for storing data of a cluster
 **
 *****************************************************************************/
static void st_getclustermedian(double *data,
                                int nvar,
                                int nech,
                                int nclusters,
                                int *clusterid,
                                double *cdata,
                                double *cache)
{
  int i, j, k, count;

  for (i = 0; i < nclusters; i++)
  {
    for (j = 0; j < nvar; j++)
    {
      count = 0;
      for (k = 0; k < nech; k++)
      {
        if (i == clusterid[k])
        {
          cache[count] = DATA(k, j);
          count++;
        }
      }
      if (count > 0)
        CDATA(i,j) = ut_median(cache, count);
      else
        CDATA(i,j) = 0.;
    }
  }
}

/*****************************************************************************/
/*!
 **  Evaluate the distance matrix
 **
 ** \param[in]  data      Array of values
 ** \param[in]  nvar      Number of samples
 ** \param[in]  nech      Number of variables
 **
 *****************************************************************************/
static double* st_get_distmatrix(double *data, int nvar, int nech)
{
  double *distmatrix;

  /* Initializations */

  distmatrix = nullptr;

  /* Core allocation */

  distmatrix = (double*) mem_alloc(sizeof(double) * nech * nech, 0);
  if (distmatrix == nullptr) return (distmatrix);

  /* Calculate the distance matrix */

  for (int i = 0; i < nech; i++)
    for (int j = 0; j < i; j++)
      DISTMATRIX(i,j) = st_distance(nvar, data, data, i, j);

  return (distmatrix);
}

/*****************************************************************************/
/*!
 **  Find the center of a cluster
 **
 ** \param[in]  distmatrix Matrix of distances
 ** \param[in]  nech       Number of variables
 ** \param[in]  nclusters  Number if clusters
 ** \param[in]  clusterid  Array of cluster number for each sample
 **
 ** \param[out] centroids  Array of centroids
 ** \param[out] errors     Array containing within-cluster distances
 **
 *****************************************************************************/
static void st_getclustermedoids(int nech,
                                 int nclusters,
                                 double *distmatrix,
                                 int *clusterid,
                                 int *centroids,
                                 double *errors)
{
  double d;
  int j;

  for (int i = 0; i < nclusters; i++)
    errors[i] = DBL_MAXIMUM;

  for (int i = 0; i < nech; i++)
  {
    d = 0.0;
    j = clusterid[i];
    for (int k = 0; k < nech; k++)
    {
      if (i == k || clusterid[k] != j) continue;
      d += (i < k ? DISTMATRIX(k, i) :
                    DISTMATRIX(i, k));
      if (d > errors[j]) break;
    }
    if (d < errors[j])
    {
      errors[j] = d;
      centroids[j] = i;
    }
  }
}

/*****************************************************************************/
/*!
 **  Perform k-means clustering on a given set of variables.
 **  The number of clusters is given in input
 **
 ** \returns Array of cluster coordinates
 **
 ** \param[in]  data      Array of values
 ** \param[in]  nvar      Number of samples
 ** \param[in]  nech      Number of variables
 ** \param[in]  nclusters Number if clusters
 ** \param[in]  npass     Number of times clustering is performed
 ** \param[in]  mode      0 for k-means and 1 for k-medians
 ** \param[in]  verbose   Verbose option
 **
 ****************************************************************************/
GSTLEARN_EXPORT double* kclusters(double *data,
                                  int nvar,
                                  int nech,
                                  int nclusters,
                                  int npass,
                                  int mode,
                                  int verbose)
{
  int i, j, k, ifound, error, niter, period, flag_same;
  int *clusterid, *tclusterid, *counts, *cmask, *mapping, *saved;
  double *cdata, *cache, total, previous, distance, tdistance, distot;

  /* Initializations */

  error = 1;
  ifound = 0;
  clusterid = tclusterid = saved = counts = cmask = mapping = nullptr;
  cdata = cache = nullptr;
  if (nclusters >= nech)
  {
    messerr(
        "The number of clusters (%d) cannot be larger than the number of points (%d)",
        nclusters, nech);
    goto label_end;
  }

  /* Core allocation */

  clusterid = (int*) mem_alloc(sizeof(int) * nech, 0);
  if (clusterid == nullptr) goto label_end;
  tclusterid = (int*) mem_alloc(sizeof(int) * nech, 0);
  if (tclusterid == nullptr) goto label_end;
  saved = (int*) mem_alloc(sizeof(int) * nech, 0);
  if (saved == nullptr) goto label_end;
  cache = (double*) mem_alloc(sizeof(double) * nech, 0);
  if (cache == nullptr) goto label_end;
  counts = (int*) mem_alloc(sizeof(int) * nclusters, 0);
  if (counts == nullptr) goto label_end;
  cmask = (int*) mem_alloc(sizeof(int) * nclusters, 0);
  if (cmask == nullptr) goto label_end;
  mapping = (int*) mem_alloc(sizeof(int) * nclusters, 0);
  if (mapping == nullptr) goto label_end;
  cdata = (double*) mem_alloc(sizeof(double) * nclusters * nvar, 0);
  if (cdata == nullptr) goto label_end;
  for (i = 0; i < nech; i++)
    clusterid[i] = 0;

  ifound = 1;
  distot = DBL_MAXIMUM;

  for (int ipass = 0; ipass < npass; ipass++)
  {
    if (verbose) message("Pass #%d\n", ipass + 1);
    total = DBL_MAXIMUM;
    niter = 0;
    period = 10;

    /* First, randomly assign elements to clusters. */
    if (verbose) message("Assign the clusters randomly\n");
    st_randomassign(nclusters, nech, tclusterid);

    for (i = 0; i < nclusters; i++)
      counts[i] = 0;
    for (i = 0; i < nech; i++)
      counts[tclusterid[i]]++;

    /* Start the loop */
    while (1)
    {
      previous = total;
      total = 0.0;

      /* Save the current cluster */
      if (niter % period == 0)
      {
        for (i = 0; i < nech; i++)
          saved[i] = tclusterid[i];
        if (period < INT_MAX / 2) period *= 2;
      }
      niter++;

      /* Find the center */
      if (mode == 0)
        st_getclustermeans(data, nvar, nech, nclusters, tclusterid, cdata,
                           cmask);
      else
        st_getclustermedian(data, nvar, nech, nclusters, tclusterid, cdata,
                            cache);

      for (i = 0; i < nech; i++)
      {
        /* Calculate the distances */
        k = tclusterid[i];
        if (counts[k] == 1) continue;

        distance = st_distance(nvar, data, cdata, i, k);
        for (j = 0; j < nclusters; j++)
        {
          if (j == k) continue;
          tdistance = st_distance(nvar, data, cdata, i, j);
          if (tdistance < distance)
          {
            distance = tdistance;
            counts[tclusterid[i]]--;
            tclusterid[i] = j;
            counts[j]++;
          }
        }
        total += distance;
      }
      if (verbose)
        message("Iteration #%d - Total Distance = %lf\n", niter, total);
      if (total >= previous) break;

      /* Check for identical clustering */
      flag_same = 1;
      for (i = 0; i < nech && flag_same; i++)
        if (saved[i] != tclusterid[i]) flag_same = 0;
      if (flag_same) break;
    }

    if (npass <= 1)
    {
      distot = total;
      break;
    }

    for (i = 0; i < nclusters; i++)
      mapping[i] = -1;
    for (i = 0; i < nech; i++)
    {
      j = tclusterid[i];
      k = clusterid[i];
      if (mapping[k] == -1)
        mapping[k] = j;
      else if (mapping[k] != j)
      {
        if (total < distot)
        {
          ifound = 1;
          distot = total;
          for (j = 0; j < nech; j++)
            clusterid[j] = tclusterid[j];
        }
        break;
      }
    }
    if (i == nech) ifound++;
  }

  /* Final printout (optional) */

  if (verbose)
  {
    message("Minimum distance = %lf\n", distot);
    st_printclusterlist(nclusters, nech, tclusterid);
    st_printclustercount(nclusters, nech, tclusterid);
  }

  /* Set the error return code */

  error = 0;

  label_end: if (error || ifound <= 0)
    cdata = (double*) mem_free((char* ) cdata);
  clusterid = (int*) mem_free((char* ) clusterid);
  tclusterid = (int*) mem_free((char* ) tclusterid);
  saved = (int*) mem_free((char* ) saved);
  counts = (int*) mem_free((char* ) counts);
  cmask = (int*) mem_free((char* ) cmask);
  mapping = (int*) mem_free((char* ) mapping);
  cache = (double*) mem_free((char* ) cache);
  cache = (double*) mem_free((char* ) cache);
  return (cdata);
}

/*****************************************************************************/
/*!
 **  Perform k-medoids clustering on a given set of variables.
 **  The number of clusters is given in input
 **
 ** \returns Array of cluster indices
 **
 ** \param[in]  data      Array of values
 ** \param[in]  nvar      Number of samples
 ** \param[in]  nech      Number of variables
 ** \param[in]  nclusters Number if clusters
 ** \param[in]  npass     Number of times clustering is performed
 ** \param[in]  verbose   Verbose option
 **
 *****************************************************************************/
GSTLEARN_EXPORT int* kmedoids(double *data,
                              int nvar,
                              int nech,
                              int nclusters,
                              int npass,
                              int verbose)
{
  int i, j, error, niter, ifound, period, flag_same;
  int *clusterid, *tclusterid, *saved, *centroids;
  double *errors, *distmatrix, distot, total, previous, distance;

  /* Initializations */

  error = 1;
  ifound = 0;
  saved = centroids = tclusterid = clusterid = nullptr;
  errors = distmatrix = nullptr;
  if (nclusters >= nech)
  {
    messerr(
        "The number of clusters (%d) cannot be larger than the number of points (%d)",
        nclusters, nech);
    goto label_end;
  }

  /* Core allocation */

  clusterid = (int*) mem_alloc(sizeof(int) * nech, 0);
  if (clusterid == nullptr) goto label_end;
  tclusterid = (int*) mem_alloc(sizeof(int) * nech, 0);
  if (tclusterid == nullptr) goto label_end;
  saved = (int*) mem_alloc(sizeof(int) * nech, 0);
  if (saved == nullptr) goto label_end;
  centroids = (int*) mem_alloc(sizeof(int) * nclusters, 0);
  if (centroids == nullptr) goto label_end;
  errors = (double*) mem_alloc(sizeof(double) * nclusters, 0);
  if (errors == nullptr) goto label_end;
  for (i = 0; i < nech; i++)
    clusterid[i] = 0;

  /* Calculate the distance matrix */

  distmatrix = st_get_distmatrix(data, nvar, nech);
  if (distmatrix == nullptr) goto label_end;

  ifound = 1;
  distot = DBL_MAXIMUM;

  for (int ipass = 0; ipass < npass; ipass++)
  {
    if (verbose) message("Pass #%d\n", ipass + 1);
    total = DBL_MAXIMUM;
    niter = 0;
    period = 10;

    /* First, randomly assign elements to clusters. */
    if (verbose) message("Assign the clusters randomly\n");
    st_randomassign(nclusters, nech, tclusterid);

    while (1)
    {
      previous = total;
      total = 0.0;

      /* Save the current cluster */
      if (niter % period == 0)
      {
        for (i = 0; i < nech; i++)
          saved[i] = tclusterid[i];
        if (period < INT_MAX / 2) period *= 2;
      }
      niter++;

      /* Find the center */
      st_getclustermedoids(nech, nclusters, distmatrix, tclusterid, centroids,
                           errors);

      for (i = 0; i < nech; i++)
      {
        /* Calculate the distance */

        distance = DBL_MAXIMUM;
        for (int icluster = 0; icluster < nclusters; icluster++)
        {
          double tdistance;
          j = centroids[icluster];
          if (i == j)
          {
            distance = 0.0;
            tclusterid[i] = icluster;
            break;
          }
          tdistance = (i > j) ? DISTMATRIX(i, j) :
                                DISTMATRIX(j, i);
          if (tdistance < distance)
          {
            distance = tdistance;
            tclusterid[i] = icluster;
          }
        }
        total += distance;
      }
      if (verbose)
        message("Iteration #%d - Total Distance = %lf\n", niter, total);
      if (total >= previous) break;

      /* Check for identical clustering */
      flag_same = 1;
      for (i = 0; i < nech && flag_same; i++)
        if (saved[i] != tclusterid[i]) flag_same = 0;
      if (flag_same) break;
    }

    if (npass <= 1)
    {
      ifound = 1;
      distot = total;
      for (j = 0; j < nech; j++)
      {
        clusterid[j] = centroids[tclusterid[j]];
      }
      break;
    }

    for (i = 0; i < nech; i++)
    {
      if (clusterid[i] != centroids[tclusterid[i]])
      {
        if (total < distot)
        {
          ifound = 1;
          distot = total;
          for (j = 0; j < nech; j++)
          {
            clusterid[j] = centroids[tclusterid[j]];
          }
        }
        break;
      }
    }
    if (i == nech) ifound++;
  }

  /* Set the error return code */

  error = 0;

  label_end: if (error || ifound <= 0) free(tclusterid);
  tclusterid = (int*) mem_free((char* ) tclusterid);
  centroids = (int*) mem_free((char* ) centroids);
  saved = (int*) mem_free((char* ) saved);
  errors = (double*) mem_free((char* ) errors);
  distmatrix = (double*) mem_free((char* ) distmatrix);
  return (clusterid);
}
