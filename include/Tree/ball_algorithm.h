/* ************************************************************************** */
/*                                                                            */
/*                                                        :::      ::::::::   */
/*   ball.h                                             :+:      :+:    :+:   */
/*                                                    +:+ +:+         +:+     */
/*   By: elee <elee@student.42.us.org>              +#+  +:+       +#+        */
/*                                                +#+#+#+#+#+   +#+           */
/*   Created: 2017/06/28 10:55:27 by elee              #+#    #+#             */
/*   Updated: 2017/06/28 20:56:18 by elee             ###   ########.fr       */
/*                                                                            */
/* ************************************************************************** */
#pragma once

#include "gstlearn_export.hpp"
#include "Basic/VectorNumT.hpp"

# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include <limits.h>
# include <math.h>
# include <float.h>
# include <sys/stat.h>
# include <time.h>

# define TRUE 1
# define FALSE 0

typedef struct s_nheap
{
  double **distances;
  int **indices;
  int n_pts;
  int n_nbrs;
} t_nheap;

typedef struct s_knn
{
  double **distances;
  int **indices;
  int n_samples;
  int n_neighbors;
} t_knn;

typedef struct s_nodedata
{
  int idx_start;
  int idx_end;
  int is_leaf;
  double radius;
} t_nodedata;

typedef struct s_btree
{
  double **data;
  int *idx_array;
  t_nodedata *node_data;
  double ***node_bounds;

  int n_samples;
  int n_features;

  int leaf_size;
  int n_levels;
  int dist_type;
  int n_nodes;
} t_btree;

/*
** ball.c
*/

GSTLEARN_EXPORT double **copy_double_arrAsVVD(VectorVectorDouble& arr);
GSTLEARN_EXPORT double **copy_double_arr(double **arr, int row, int col);
GSTLEARN_EXPORT t_btree* btree_init(double **data,
                                    int n_samples,
                                    int n_features,
                                    int leaf_size,
                                    int dist_type);
GSTLEARN_EXPORT t_knn btree_query(t_btree *b,
                                  double **x,
                                  int n_samples,
                                  int n_features,
                                  int k);
GSTLEARN_EXPORT void  free_2d_double(double **arr, int row);
GSTLEARN_EXPORT void  free_2d_int(int **arr, int row);
GSTLEARN_EXPORT void  free_tree(t_btree *tree);
GSTLEARN_EXPORT void  free_knn(t_knn knn, int row);


/*
** metrics.c
*/

double manhattan_dist(const double *x1, const double *x2, int size);
double euclidean_dist(const double *x1, const double *x2, int size);

/*
** neighbors_heap.c
*/

t_nheap	*nheap_init(int n_pts, int n_nbrs);
double	nheap_largest(t_nheap *h, int row);
int		nheap_push(t_nheap *h, int row, double val, int i_val);
t_knn	nheap_get_arrays(t_nheap *h);

