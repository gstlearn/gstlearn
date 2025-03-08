/*
                                      ball

Original Author: Eung Bum Lee
Website: https://42.fr
License: MIT
*/

/* ************************************************************************** */
/*                                                                            */
/*                                                        :::      ::::::::   */
/*   neighbors_heap.c                                   :+:      :+:    :+:   */
/*                                                    +:+ +:+         +:+     */
/*   By: elee <elee@student.42.us.org>              +#+  +:+       +#+        */
/*                                                +#+#+#+#+#+   +#+           */
/*   Created: 2017/06/28 10:45:06 by elee              #+#    #+#             */
/*   Updated: 2017/06/28 20:17:58 by elee             ###   ########.fr       */
/*                                                                            */
/* ************************************************************************** */

/*
Modified by MINES Paris / ARMINES (2023)
Authors: gstlearn Team
Website: https://gstlearn.org
License: BSD 3-clause
*/

#include "Tree/ball_algorithm.h"

void dual_swap(double* darr, int* iarr, int i1, int i2)
{
  double dtmp = darr[i1];
  darr[i1]    = darr[i2];
  darr[i2]    = dtmp;
  int itmp    = iarr[i1];
  iarr[i1]    = iarr[i2];
  iarr[i2]    = itmp;
}

t_nheap* nheap_init(int n_pts, int n_nbrs)
{
  t_nheap* h   = (t_nheap*)malloc(sizeof(t_nheap));
  h->n_pts     = n_pts;
  h->n_nbrs    = n_nbrs;
  h->distances = (double**)malloc(sizeof(double*) * n_pts);
  for (int i = 0; i < n_pts; i++)
  {
    h->distances[i] = (double*)malloc(sizeof(double) * n_nbrs);
    for (int j = 0; j < n_nbrs; j++) h->distances[i][j] = INFINITY;
  }
  h->indices = (int**)malloc(sizeof(int*) * n_pts);
  for (int i = 0; i < n_pts; i++)
  {
    h->indices[i] = (int*)calloc(n_nbrs, sizeof(int));
    for (int j = 0; j < n_nbrs; j++) h->indices[i][j] = ITEST;
  }
  return (h);
}

t_nheap* nheap_free(t_nheap* heap)
{
  if (heap == nullptr) return heap;
  int n_pts = heap->n_pts;
  for (int i = 0; i < n_pts; i++) free(heap->distances[i]);
  free(heap->distances);
  for (int i = 0; i < n_pts; i++) free(heap->indices[i]);
  free(heap->indices);
  free(heap);
  heap = nullptr;
  return heap;
}

void nheap_load(t_nheap* heap, t_btree* b, const double** x)
{
  double dist;
  for (int i = 0; i < heap->n_pts; i++)
  {
    dist = min_dist(b, 0, x[i]);
    query_depth_first(b, 0, x[i], i, heap, dist);
  }
}

double nheap_largest(t_nheap* h, int row)
{
  return h->distances[row][0];
}

int nheap_push(t_nheap* h, int row, double val, int i_val)
{
  int ic1, ic2, i_swap;

  int size         = h->n_nbrs;
  double* dist_arr = h->distances[row];
  int* ind_arr     = h->indices[row];

  // if distance is already greater than the furthest element, don't push
  if (val > dist_arr[0]) return (0);

  // insert the values at position 0
  dist_arr[0] = val;
  ind_arr[0]  = i_val;

  // descend the heap, swapping values until the max heap criterion is met
  int i = 0;
  while (TRUE)
  {
    ic1 = 2 * i + 1;
    ic2 = ic1 + 1;

    if (ic1 >= size)
      break;
    if (ic2 >= size)
    {
      if (dist_arr[ic1] > val)
        i_swap = ic1;
      else
        break;
    }
    else if (dist_arr[ic1] >= dist_arr[ic2])
    {
      if (val < dist_arr[ic1])
        i_swap = ic1;
      else
        break;
    }
    else
    {
      if (val < dist_arr[ic2])
        i_swap = ic2;
      else
        break;
    }
    dist_arr[i] = dist_arr[i_swap];
    ind_arr[i]  = ind_arr[i_swap];
    i           = i_swap;
  }

  dist_arr[i] = val;
  ind_arr[i]  = i_val;

  return (0);
}

void simultaneous_sort(double* dist, int* idx, int size)
{
  int pivot_idx, store_idx;
  double pivot_val;

  if (size <= 1)
    ;
  else if (size == 2)
  {
    if (dist[0] > dist[1]) dual_swap(dist, idx, 0, 1);
  }
  else if (size == 3)
  {
    if (dist[0] > dist[1]) dual_swap(dist, idx, 0, 1);
    if (dist[1] > dist[2])
    {
      dual_swap(dist, idx, 1, 2);
      if (dist[0] > dist[1]) dual_swap(dist, idx, 0, 1);
    }
  }
  else
  {
    pivot_idx = size / 2;
    if (dist[0] > dist[size - 1]) dual_swap(dist, idx, 0, size - 1);
    if (dist[size - 1] > dist[pivot_idx])
    {
      dual_swap(dist, idx, size - 1, pivot_idx);
      if (dist[0] > dist[size - 1]) dual_swap(dist, idx, 0, size - 1);
    }
    pivot_val = dist[size - 1];

    store_idx = 0;
    for (int i = 0; i < size - 1; i++)
    {
      if (dist[i] < pivot_val)
      {
        dual_swap(dist, idx, i, store_idx);
        store_idx++;
      }
    }
    dual_swap(dist, idx, store_idx, size - 1);
    pivot_idx = store_idx;
    if (pivot_idx > 1) simultaneous_sort(dist, idx, pivot_idx);
    if (pivot_idx * 2 < size)
      simultaneous_sort(dist + pivot_idx + 1, idx + pivot_idx + 1,
                        size - pivot_idx - 1);
  }
}

void nheap_sort(t_nheap* h)
{
  for (int row = 0; row < h->n_pts; row++)
    simultaneous_sort(h->distances[row], h->indices[row], h->n_nbrs);
}
