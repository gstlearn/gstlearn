/* MIT License
 *
 * Copyright (C) 2012, 2013, 2014 James McLaughlin et al.  All rights reserved.
 * https://github.com/udp/json-parser
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *   notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *   notice, this list of conditions and the following disclaimer in the
 *   documentation and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE AUTHOR AND CONTRIBUTORS ``AS IS'' AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED.  IN NO EVENT SHALL THE AUTHOR OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
 * OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
 * OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
 * SUCH DAMAGE.
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

// Local modifications by D. RENARD
// Insert 'Basic/AStringable.hpp'
// Replace 'printf' by 'messerr'
// Replace 'exit(-1)' by 'return nullptr' (which must be checked by calling function)
// Return 'exit(-1)' by returning a t_knn structure where n-samples is set to -1.

#include "Basic/AStringable.hpp"
#include "Tree/ball.hpp"

double	**copy_double_arr(double **arr, int row, int col)
{
	double	**copy;
	int		i, j;
	
	copy = (double**)malloc(sizeof(double*) * row);
	for (i = 0; i < row; i++)
	{
		copy[i] = (double*)malloc(sizeof(double) * col);
		for (j = 0; j < col; j++)
			copy[i][j] = arr[i][j];
	}
	return (copy);
}

void	swap(int *arr, int i1, int i2)
{
	int	tmp;
	
	tmp = arr[i1];
	arr[i1] = arr[i2];
	arr[i2] = tmp;
}

void	btree_zero(t_btree *b)
{
	b->data = NULL;
	b->idx_array = NULL;
	b->node_data = NULL;
	b->node_bounds = NULL;

	b->leaf_size = 40;
	b->n_levels = 0;
	b->n_nodes = 0;
}

int		init_node(t_btree *b, int i_node, int idx_start, int idx_end)
{
	int		n_features;
	int		n_points;
	int		i, j;
	double	radius;
	double	*centroid;

	n_features = b->n_features;
	n_points = idx_end - idx_start;
	centroid = b->node_bounds[0][i_node];
	
	for (j = 0; j < n_features; j++)
		centroid[j] = 0.0;
	
	for (i = idx_start; i < idx_end; i++)
		for (j = 0; j < n_features; j++)
			centroid[j] += b->data[b->idx_array[i]][j];

	for (j = 0; j < n_features; j++)
		centroid[j] /= n_points;

	radius = 0.0;
	for (i = idx_start; i < idx_end; i++)
		radius = fmax(radius, manhattan_dist(centroid, b->data[b->idx_array[i]], n_features));

	b->node_data[i_node].radius = radius;
	b->node_data[i_node].idx_start = idx_start;
	b->node_data[i_node].idx_end = idx_end;
	return (0);
}

int		find_node_split_dim(double **data, int *node_indices, int n_features, int n_points)
{
	double	min_val, max_val, val, spread, max_spread;
	int		i, j, j_max;

	j_max = 0;
	max_spread = 0;
	for (j = 0; j < n_features; j++)
	{
		max_val = data[node_indices[0]][j];
		min_val = max_val;
		for (i = 1; i < n_points; i++)
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

int		partition_node_indices(double **data, int *node_indices, int split_dim, int split_index,
								int n_features, int n_points)
{
	(void)n_features;
	int		left, right, midindex, i;
	double	d1, d2;

	left = 0;
	right = n_points - 1;

	while (TRUE)
	{
		midindex = left;
		for (i = left; i < right; i++)
		{
			d1 = data[node_indices[i]][split_dim];
			d2 = data[node_indices[right]][split_dim];
			if (d1 < d2)
			{
				swap(node_indices, i, midindex);
				midindex += 1;
			}
		}
		swap(node_indices, midindex, right);
		if (midindex == split_index)
			break ;
		else if (midindex < split_index)
			left = midindex + 1;
		else
			right = midindex - 1;
	}
	return (0);
}

void	recursive_build(t_btree *b, int i_node, int idx_start, int idx_end)
{
	int	imax;
	int	n_features;
	int	n_points;
	int	n_mid;

	n_features = b->n_features;
	n_points = idx_end - idx_start;
	n_mid = n_points / 2;

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
		partition_node_indices(b->data, b->idx_array, imax, n_mid, n_features, n_points);
		recursive_build(b, 2 * i_node + 1, idx_start, idx_start + n_mid);
		recursive_build(b, 2 * i_node + 2, idx_start + n_mid, idx_end);
	}
}

t_btree	*btree_init(double **data, int n_samples, int n_features, int leaf_size)
{
	t_btree	*b;
	int		i, j;

	b = (t_btree*)malloc(sizeof(t_btree));
	btree_zero(b);

	b->data = copy_double_arr(data, n_samples, n_features);
	b->leaf_size = leaf_size;
	
	if (leaf_size < 1)
	{
		messerr("leaf_size must be greater than or equal to 1\n");
		return nullptr;
	}

	b->n_samples = n_samples;
	b->n_features = n_features;

	b->n_levels = log2(fmax(1, (b->n_samples - 1) / b->leaf_size)) + 1;
	b->n_nodes = pow(2.0, b->n_levels) - 1;

	b->idx_array = (int*)malloc(sizeof(int) * b->n_samples);
	for (i = 0; i < b->n_samples; i++)
		b->idx_array[i] = i;
	b->node_data = (t_nodedata*)calloc(b->n_nodes, sizeof(t_nodedata));
	b->node_bounds = (double***)malloc(sizeof(double**));
	b->node_bounds[0] = (double**)malloc(sizeof(double*) * b->n_nodes);
	for (i = 0; i < b->n_nodes; i++)
	{
		b->node_bounds[0][i] = (double*)malloc(sizeof(double) * b->n_features);
		for (j = 0; j < b->n_features; j++)
			b->node_bounds[0][i][j] = 0.0;
	}
	recursive_build(b, 0, 0, b->n_samples);
	return (b);
}

int		query_depth_first(t_btree *b, int i_node, double *pt, int i_pt, t_nheap *heap, double dist)
{
	t_nodedata	node_info = b->node_data[i_node];
	double		dist_pt, dist1, dist2;
	int			i, i1, i2;

	//case 1: query point is outside node radius: trim it from the query
	if (dist > nheap_largest(heap, i_pt))
	{
		;
	}
	//case 2: this is a leaf node. Update set of nearby points
	else if (node_info.is_leaf)
	{
		for (i = node_info.idx_start; i < node_info.idx_end; i++)
		{
			dist_pt = manhattan_dist(pt, b->data[b->idx_array[i]], b->n_features);
			if (dist_pt < nheap_largest(heap, i_pt))
				nheap_push(heap, i_pt, dist_pt, b->idx_array[i]);
		}
	}
	//case 3: Node is not a leaf, Recursively query sub-nodes starting with the closest
	else
	{
		i1 = 2 * i_node +1;
		i2 = i1 +1;
		dist1 = min_dist(b, i1, pt); //implement min_rdist
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

t_knn	btree_query(t_btree *b, double **x, int n_samples, int n_features, int k)
{
	t_nheap	*heap;
	double	dist;
	int		i;
	t_knn	output;

	if (n_features != b->n_features)
	{
		messerr("query data dimension must match training data dimension.");
		output.n_samples = -1;
		return output;
	}
	if (b->n_samples < k)
	{
		messerr("k must be less than or equal to the number of training points.");
    output.n_samples = -1;
    return output;
	}
	heap = nheap_init(n_samples, k);
	for (i = 0; i < n_samples; i++)
	{
		dist = min_dist(b, 0, x[i]);
		query_depth_first(b, 0, x[i], i, heap, dist);
	}
	output = nheap_get_arrays(heap);
	return (output);
}

void	free_2d_double(double **arr, int row)
{
	int	i;

	for (i = 0; i < row; i++)
		free(arr[i]);
	free(arr);
}

void	free_2d_int(int **arr, int row)
{
	int	i;

	for (i = 0; i < row; i++)
		free(arr[i]);
	free(arr);
}

void	free_tree(t_btree *tree)
{
	free_2d_double(tree->data, tree->n_samples);
	free(tree->idx_array);
	free(tree->node_data);
	free_2d_double(tree->node_bounds[0], tree->n_nodes);
	free(tree->node_bounds);
	free(tree);
}

void	free_knn(t_knn knn, int row)
{
	free_2d_double(knn.distances, row);
	free_2d_int(knn.indices, row);
}
