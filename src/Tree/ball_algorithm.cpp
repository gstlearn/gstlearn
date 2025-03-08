/*
                                      ball

Original Author: Eung Bum Lee
Website: https://42.fr
License: MIT
*/

/* ************************************************************************** */
/*                                                                            */
/*                                                        :::      ::::::::   */
/*   ball.c                                             :+:      :+:    :+:   */
/*                                                    +:+ +:+         +:+     */
/*   By: elee <elee@student.42.us.org>              +#+  +:+       +#+        */
/*                                                +#+#+#+#+#+   +#+           */
/*   Created: 2017/06/28 10:45:02 by elee              #+#    #+#             */
/*   Updated: 2017/06/28 20:56:01 by elee             ###   ########.fr       */
/*                                                                            */
/* ************************************************************************** */

/*
Modified by MINES Paris / ARMINES (2023)
Authors: gstlearn Team
Website: https://gstlearn.org
License: BSD 3-clause
*/

#include "Tree/ball_algorithm.h"
#include "Tree/KNN.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/VectorHelper.hpp"
#include "Space/SpacePoint.hpp"

static double (*st_distance_function)(const double*, const double*, int) = euclidean_distance;

double **copy_double_arrAsVVD(const VectorVectorDouble& arr)
{
  int col = (int) arr.size();
  int row = (int) arr[0].size();

  double** copy = (double**)malloc(sizeof(double*) * row);
  for (int i = 0; i < row; i++)
  {
    copy[i] = (double*)malloc(sizeof(double) * col);
    for (int j = 0; j < col; j++)
      copy[i][j] = arr[j][i];
  }
  return (copy);
}

double **copy_double_arr(const double **arr, int row, int col)
{
	double** copy = (double**)malloc(sizeof(double*) * row);
	for (int i = 0; i < row; i++)
	{
		copy[i] = (double*)malloc(sizeof(double) * col);
		for (int j = 0; j < col; j++)
			copy[i][j] = arr[i][j];
	}
	return (copy);
}

bool* init_bool_arr(int col, bool status)
{
  bool* copy = (bool*)malloc(sizeof(bool) * col);
  for (int j = 0; j < col; j++)
    copy[j] = status;
  return copy;
}

VectorVectorDouble copy_double_toVVD(const double** arr, int row, int col)
{
  VectorVectorDouble copy(row);
  for (int i = 0; i < row; i++)
  {
    copy[i].resize(col);
    for (int j = 0; j < col; j++) copy[i][j] = arr[i][j];
  }
  return (copy);
}

int **copy_int_arr(const int **arr, int row, int col)
{
  int** copy = (int**)malloc(sizeof(int*) * row);
  for (int i = 0; i < row; i++)
  {
    copy[i] = (int*)malloc(sizeof(int) * col);
    for (int j = 0; j < col; j++)
      copy[i][j] = arr[i][j];
  }
  return (copy);
}

VectorVectorInt copy_int_toVVI(const int** arr, int row, int col)
{
  VectorVectorInt copy(row);
  for (int i = 0; i < row; i++)
  {
    copy[i].resize(col);
    for (int j = 0; j < col; j++) copy[i][j] = arr[i][j];
  }
  return (copy);
}

void swap(int *arr, int i1, int i2)
{
	if (i1 == i2) return;
	int tmp = arr[i1];
	arr[i1] = arr[i2];
	arr[i2] = tmp;
}

void btree_zero(t_btree *b)
{
	b->data = NULL;
  b->accept = NULL;
	b->idx_array = NULL;
	b->node_data = NULL;
	b->node_bounds = NULL;

	b->leaf_size = 40;
	b->n_levels = 0;
	b->n_nodes = 0;
}

int init_node(t_btree *b, int i_node, int idx_start, int idx_end)
{
  int n_features = b->n_features;
  int n_points = idx_end - idx_start;
  double* centroid = b->node_bounds[0][i_node];

  for (int j = 0; j < n_features; j++)
    centroid[j] = 0.0;

  for (int i = idx_start; i < idx_end; i++)
    for (int j = 0; j < n_features; j++)
      centroid[j] += b->data[b->idx_array[i]][j];

  for (int j = 0; j < n_features; j++)
    centroid[j] /= n_points;

  double radius = 0.0;
  for (int i = idx_start; i < idx_end; i++)
    radius = fmax(radius, st_distance_function(centroid, b->data[b->idx_array[i]], n_features));

  b->node_data[i_node].radius = radius;
  b->node_data[i_node].idx_start = idx_start;
  b->node_data[i_node].idx_end = idx_end;
  return (0);
}

int find_node_split_dim(double **data, const int *node_indices, int n_features, int n_points)
{
	double	min_val, max_val, val, spread;

	int j_max = 0;
	double max_spread = 0;
	for (int j = 0; j < n_features; j++)
	{
		max_val = data[node_indices[0]][j];
		min_val = max_val;
		for (int i = 1; i < n_points; i++)
		{
			val = data[node_indices[i]][j];
			max_val = fmax(max_val, val);
			min_val = fmin(min_val, val);
		}
		spread = max_val - min_val;
		if (spread > max_spread)
		{
			max_spread = spread;
			j_max = j;
		}
	}
	return (j_max);
}

int partition_node_indices(double **data, int *node_indices, int split_dim, int n_points, int split_index)
{
  int   midindex;
  double  d1, d2;

  int left = 0;
  int right = n_points - 1;

  while (TRUE)
  {
    midindex = left;
    for (int i = left; i < right; i++)
    {
      d1 = data[node_indices[i]][split_dim];
      d2 = data[node_indices[right]][split_dim];
      if (d1 < d2)
      {
        swap(node_indices, i, midindex);
        midindex++;
      }
    }
    swap(node_indices, midindex, right);
    if (midindex == split_index)
      break ;
    if (midindex < split_index)
      left = midindex + 1;
    else
      right = midindex - 1;
  }

  return (0);
}

void recursive_build(t_btree *b, int i_node, int idx_start, int idx_end)
{
	int	imax;
	int n_features = b->n_features;
	int n_points = idx_end - idx_start;
	int n_mid = n_points / 2;

	//initialize the node data
	init_node(b, i_node, idx_start, idx_end);

	if (2 * i_node + 1 >= b->n_nodes)
	{
		b->node_data[i_node].is_leaf = TRUE;
		if (idx_end - idx_start > 2 * b->leaf_size)
			messerr("Memory layout is flawed: not enough nodes allocated");
	}
	else if (idx_end - idx_start < 2)
	{
		messerr("Memory layout is flawed: too many nodes allocated");
		b->node_data[i_node].is_leaf = TRUE;
	}
	else
	{
		b->node_data[i_node].is_leaf = FALSE;
		imax = find_node_split_dim(b->data, b->idx_array, n_features, n_points);
		partition_node_indices(b->data, &b->idx_array[idx_start], imax, n_points, n_mid);
		recursive_build(b, 2 * i_node + 1, idx_start, idx_start + n_mid);
		recursive_build(b, 2 * i_node + 2, idx_start + n_mid, idx_end);
	}
}

void define_dist_function(double (*dist_function)(const double* x1,
                                                  const double* x2,
                                                  int size),
                          int default_distance_function)
{
  if (dist_function != nullptr)
  {
    st_distance_function = dist_function;
  }
  else
  {
    if (default_distance_function == 1)
      st_distance_function = euclidean_distance;
    if (default_distance_function == 2)
      st_distance_function = manhattan_distance;
  }
}

t_btree* btree_init(const double** data,
                    int n_samples,
                    int n_features,
                    bool has_constraints,
                    double (*dist_function)(const double* x1,
                                            const double* x2,
                                            int size),
                    int leaf_size,
                    int default_distance_function)
{
  t_btree* b = (t_btree*)malloc(sizeof(t_btree));
  btree_zero(b);

  b->data = copy_double_arr(data, n_samples, n_features);
  if (has_constraints)
    b->accept = init_bool_arr(n_samples, false);
  b->leaf_size = leaf_size;

  if (leaf_size < 1)
  {
    messerr("leaf_size must be greater than or equal to 1\n");
    return nullptr;
  }

  // Define the relevant distance function
  define_dist_function(dist_function, default_distance_function);

	b->n_samples = n_samples;
	b->n_features = n_features;

	b->n_levels = log2(fmax(1, (b->n_samples - 1) / b->leaf_size)) + 1;
	b->n_nodes = pow(2.0, b->n_levels) - 1;

	b->idx_array = (int*)malloc(sizeof(int) * b->n_samples);
	for (int i = 0; i < b->n_samples; i++)
		b->idx_array[i] = i;
	b->node_data = (t_nodedata*)calloc(b->n_nodes, sizeof(t_nodedata));
	b->node_bounds = (double***)malloc(sizeof(double**));
	b->node_bounds[0] = (double**)malloc(sizeof(double*) * b->n_nodes);
	for (int i = 0; i < b->n_nodes; i++)
	{
		b->node_bounds[0][i] = (double*)malloc(sizeof(double) * b->n_features);
		for (int j = 0; j < b->n_features; j++)
			b->node_bounds[0][i][j] = 0.0;
	}
	recursive_build(b, 0, 0, b->n_samples);

	return (b);
}

/**
 * @brief Calculate the distance between the current 'pt' and the centroid of node 'i_node'
 * Returns 0 if 'pt' belongs to the node
 * 
 * @param tree t_btree structure
 * @param i_node Rank of the target node
 * @param pt     Characteristics of the target SpacePoint
 * @return double Minimum distance or 0
 */
double min_dist(t_btree *tree, int i_node, const double *pt)
{
  double dist_pt = st_distance_function(pt, tree->node_bounds[0][i_node], tree->n_features);
  return (fmax(0.0, dist_pt - tree->node_data[i_node].radius));
}

int query_depth_first(t_btree *b, int i_node, const double *pt, int i_pt, t_nheap *heap, double dist)
{
  t_nodedata node_info = b->node_data[i_node];
  double dist_pt, dist1, dist2;
  int i1, i2;

  // case 1: query point is outside node radius: trim it from the query
  if (dist > nheap_largest(heap, i_pt))
  {
    ;
  }
  // case 2: this is a leaf node. Update set of nearby points
  else if (node_info.is_leaf)
  {
    for (int i = node_info.idx_start; i < node_info.idx_end; i++)
    {
      int j = b->idx_array[i];
      if (b->accept != nullptr && !b->accept[j]) 
        continue;
      dist_pt = st_distance_function(pt, b->data[j], b->n_features);
      if (dist_pt < nheap_largest(heap, i_pt))
        nheap_push(heap, i_pt, dist_pt, j);
    }
  }
  // case 3: Node is not a leaf, Recursively query sub-nodes starting with the closest
  else
  {
    i1    = 2 * i_node + 1;
    i2    = i1 + 1;
    dist1 = min_dist(b, i1, pt); // implement min_rdist
    dist2 = min_dist(b, i2, pt);
    if (dist1 <= dist2)
    {
      query_depth_first(b, i1, pt, i_pt, heap, dist1);
      query_depth_first(b, i2, pt, i_pt, heap, dist2);
    }
    else
    {
      query_depth_first(b, i2, pt, i_pt, heap, dist2);
      query_depth_first(b, i1, pt, i_pt, heap, dist1);
    }
  }
  return (0);
}

void free_2d_double(double **arr, int row)
{
	for (int i = 0; i < row; i++)
		free(arr[i]);
	free(arr);
}

void free_2d_int(int **arr, int row)
{
	for (int i = 0; i < row; i++)
		free(arr[i]);
	free(arr);
}

void free_tree(t_btree *tree)
{
  if (tree == nullptr) return;
	free_2d_double(tree->data, tree->n_samples);
	free(tree->idx_array);
	free(tree->node_data);
	free_2d_double(tree->node_bounds[0], tree->n_nodes);
	free(tree->node_bounds);
	free(tree);
}

void btree_display(const t_btree *tree, int level)
{
  if (tree == nullptr) return;

  mestitle(0, "Ball Tree");
  message("- Number of samples = %d\n", tree->n_samples);
  message("- Number of Features = %d\n", tree->n_features);
  message("- Number of levels = %d\n", tree->n_levels);
  message("- Number of nodes = %d\n", tree->n_nodes);
  message("- Size of leaf = %d\n", tree->leaf_size);
  if (level < 0) return;

  // Loop on the nodes

  mestitle(1, "List of nodes");
  for (int i_node = 0; i_node < tree->n_nodes; i_node++)
  {
    t_nodedata* info = &tree->node_data[i_node];

    message("Node #%3d/%3d - Indices [%5d; %5d[ - Radius = %lf - Centroid = ",
            i_node, tree->n_nodes, info->idx_start, info->idx_end, info->radius);
    for (int j = 0; j < tree->n_features; j++)
      message("%lf ", tree->node_bounds[0][i_node][j]);
    message("\n");

    if (level > 0 && info->is_leaf)
    {
      message(" Sample indices = ");
      for (int is = info->idx_start; is < info->idx_end; is++)
        message(" %d", tree->idx_array[is]);
      message("\n");
    }
  }
}

/**
 * Returns the Manhattan distance between two points
 * @param x1 Vector of coordinates for the first point
 * @param x2 Vector of coordinates for the second point
 * @param size Number of coordinates
 * @return
 */
double manhattan_distance(const double* x1, const double* x2, int size)
{
  double delta;
  double d1 = 0.;
  for (int i = 0; i < size; i++)
  {
    delta = fabs(x1[i] - x2[i]);
    d1 += delta;
  }
  return (d1);
}

/**
 * Returns the Standard distance between two points
 * @param x1 Vector of coordinates for the first point
 * @param x2 Vector of coordinates for the second point
 * @param size Number of coordinates
 * @return
 */
double euclidean_distance(const double* x1, const double* x2, int size)
{
  SpacePoint p1;
  SpacePoint p2;
  p1.setCoords(x1, size);
  p2.setCoords(x2, size);
  return p1.getDistance(p2);
}
