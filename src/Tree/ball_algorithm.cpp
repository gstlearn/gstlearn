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
#include "Space/SpacePoint.hpp"

namespace gstlrn
{

  void swap(Id* arr, Id i1, Id i2)
  {
    std::swap(arr[i1], arr[i2]);
  }

  Id t_btree::init_node(Id i_node, Id idx_start, Id idx_end)
  {
    Id n_featuresLocal = this->n_features;
    Id n_points = idx_end - idx_start;
    auto centroid = this->node_bounds.getRow(i_node);

    for (Id j = 0; j < n_featuresLocal; j++) centroid[j] = 0.0;

    for (Id i = idx_start; i < idx_end; i++)
      for (Id j = 0; j < n_featuresLocal; j++)
        centroid[j] += this->data(this->idx_array[i], j);

    for (Id j = 0; j < n_featuresLocal; j++) centroid[j] /= n_points;

    double radius = 0.0;
    const auto dist_func = this->default_distance_function == 1
                           ? euclidean_distance
                           : manhattan_distance;
    for (Id i = idx_start; i < idx_end; i++)
      radius = fmax(radius,
                    dist_func(centroid.data(),
                              this->data.getRow(this->idx_array[i]).data(),
                              n_featuresLocal));

    this->node_data[i_node].radius = radius;
    this->node_data[i_node].idx_start = idx_start;
    this->node_data[i_node].idx_end = idx_end;
    return (0);
  }

  Id find_node_split_dim(const MatrixT<double>& data,
                         const VectorInt& node_indices,
                         Id n_features,
                         Id n_points)
  {
    double min_val, max_val, val, spread;

    Id j_max = 0;
    double max_spread = 0;
    for (Id j = 0; j < n_features; j++)
    {
      max_val = data(node_indices[0], j);
      min_val = max_val;
      for (Id i = 1; i < n_points; i++)
      {
        val = data(node_indices[i], j);
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

  Id partition_node_indices(const MatrixT<double>& data,
                            Id* node_indices,
                            Id split_dim,
                            Id n_points,
                            Id split_index)
  {
    Id midindex;
    double d1, d2;

    Id left = 0;
    Id right = n_points - 1;

    while (true)
    {
      midindex = left;
      for (Id i = left; i < right; i++)
      {
        d1 = data(node_indices[i], split_dim);
        d2 = data(node_indices[right], split_dim);
        if (d1 < d2)
        {
          swap(node_indices, i, midindex);
          midindex++;
        }
      }
      swap(node_indices, midindex, right);
      if (midindex == split_index) break;
      if (midindex < split_index)
        left = midindex + 1;
      else
        right = midindex - 1;
    }

    return (0);
  }

  void t_btree::recursive_build(Id i_node, Id idx_start, Id idx_end)
  {
    Id imax;
    Id n_featuresLocal = this->n_features;
    Id n_points = idx_end - idx_start;
    Id n_mid = n_points / 2;

    // initialize the node data
    init_node(i_node, idx_start, idx_end);

    if (2 * i_node + 1 >= this->n_nodes)
    {
      this->node_data[i_node].is_leaf = true;
      if (idx_end - idx_start > 2 * this->leaf_size)
        messerr("Memory layout is flawed: not enough nodes allocated");
    }
    else if (idx_end - idx_start < 2)
    {
      messerr("Memory layout is flawed: too many nodes allocated");
      this->node_data[i_node].is_leaf = true;
    }
    else
    {
      this->node_data[i_node].is_leaf = false;
      imax = find_node_split_dim(this->data,
                                 this->idx_array,
                                 n_featuresLocal,
                                 n_points);
      partition_node_indices(this->data,
                             &this->idx_array[idx_start],
                             imax,
                             n_points,
                             n_mid);
      recursive_build(2 * i_node + 1, idx_start, idx_start + n_mid);
      recursive_build(2 * i_node + 2, idx_start + n_mid, idx_end);
    }
  }

  t_btree::t_btree(MatrixT<double>&& data,
                   Id n_samples,
                   Id n_features,
                   bool all_available,
                   Id leaf_size,
                   Id default_distance_function)
    : data(std::move(data))
    , leaf_size(40)
    , n_levels(0)
    , n_nodes(0)
    , default_distance_function(default_distance_function)
  {
    this->available.resize(n_samples, all_available);
    this->leaf_size = leaf_size;

    if (leaf_size < 1)
    {
      messerr("leaf_size must be greater than or equal to 1\n");
      return;
    }

    this->n_samples = n_samples;
    this->n_features = n_features;

    this->n_levels =
      log2(fmax(1, static_cast<double>(this->n_samples - 1) / this->leaf_size))
      + 1;
    this->n_nodes = pow(2.0, this->n_levels) - 1;

    this->idx_array.resize(this->n_samples);
    for (Id i = 0; i < this->n_samples; i++) this->idx_array[i] = i;
    this->node_data.resize(this->n_nodes);
    this->node_bounds.resize(this->n_nodes, this->n_features, 0.0);
    recursive_build(0, 0, this->n_samples);
  }

  /**
   * @brief Calculate the distance between the current 'pt' and the centroid of node 'i_node'
   * Returns 0 if 'pt' belongs to the node
   *
   * @param i_node Rank of the target node
   * @param pt     Characteristics of the target SpacePoint
   * @return double Minimum distance or 0
   */
  double t_btree::min_dist(Id i_node, const constvect pt) const
  {
    const auto dist_func = this->default_distance_function == 1
                           ? euclidean_distance
                           : manhattan_distance;
    double dist_pt = dist_func(pt.data(),
                               this->node_bounds.getRow(i_node).data(),
                               this->n_features);
    return (fmax(0.0, dist_pt - this->node_data[i_node].radius));
  }

  Id t_btree::query_depth_first(Id i_node,
                                const constvect pt,
                                Id i_pt,
                                t_nheap& heap,
                                double dist) const
  {
    t_nodedata node_info = this->node_data[i_node];
    double dist_pt, dist1, dist2;
    Id i1, i2;

    // case 1: query point is outside node radius: trim it from the query
    if (dist > heap.largest(i_pt))
    {
      ;
    }
    // case 2: this is a leaf node. Update set of nearby points
    else if (node_info.is_leaf)
    {
      const auto dist_func = this->default_distance_function == 1
                             ? euclidean_distance
                             : manhattan_distance;
      for (Id i = node_info.idx_start; i < node_info.idx_end; i++)
      {
        Id j = this->idx_array[i];
        if (!this->available[j]) continue;
        dist_pt =
          dist_func(pt.data(), this->data.getRow(j).data(), this->n_features);
        if (dist_pt < heap.largest(i_pt)) heap.push(i_pt, dist_pt, j);
      }
    }
    // case 3: Node is not a leaf, Recursively query sub-nodes starting with the closest
    else
    {
      i1 = 2 * i_node + 1;
      i2 = i1 + 1;
      dist1 = min_dist(i1, pt); // implement min_rdist
      dist2 = min_dist(i2, pt);
      if (dist1 <= dist2)
      {
        query_depth_first(i1, pt, i_pt, heap, dist1);
        query_depth_first(i2, pt, i_pt, heap, dist2);
      }
      else
      {
        query_depth_first(i2, pt, i_pt, heap, dist2);
        query_depth_first(i1, pt, i_pt, heap, dist1);
      }
    }
    return (0);
  }

  void t_btree::display(Id level) const
  {
    mestitle(0, "Ball Tree");
    message("- Number of samples = %d\n", this->n_samples);
    message("- Number of Features = %d\n", this->n_features);
    message("- Number of levels = %d\n", this->n_levels);
    message("- Number of nodes = %d\n", this->n_nodes);
    message("- Size of leaf = %d\n", this->leaf_size);
    if (level < 0) return;

    // Loop on the nodes

    mestitle(1, "List of nodes");
    for (Id i_node = 0; i_node < this->n_nodes; i_node++)
    {
      const auto& info = this->node_data[i_node];

      message("Node #%3d/%3d - Indices [%5d; %5d[ - Radius = %lf - Centroid = ",
              i_node,
              this->n_nodes,
              info.idx_start,
              info.idx_end,
              info.radius);
      for (Id j = 0; j < this->n_features; j++)
        message("%lf ", this->node_bounds(i_node, j));
      message("\n");

      if (level > 0 && info.is_leaf)
      {
        message(" Sample indices = ");
        for (Id is = info.idx_start; is < info.idx_end; is++)
          message(" %d", this->idx_array[is]);
        message("\n");
      }
    }
  }

  /**
   * Returns the Manhattan distance between two points
   * @param x1 Vector of coordinates for the first point
   * @param x2 Vector of coordinates for the second point
   * @param n_features Number of coordinates
   * @return
   */
  double manhattan_distance(const double* x1, const double* x2, Id n_features)
  {
    double delta;
    double d1 = 0.;
    for (Id i = 0; i < n_features; i++)
    {
      delta = fabs(x1[i] - x2[i]);
      d1 += delta;
    }
    return (d1);
  }

  /**
   * Returns the Standard Euclidean distance between two points
   * @param x1 Vector of coordinates for the first point
   * @param x2 Vector of coordinates for the second point
   * @param n_features Number of coordinates
   * @return
   */
  double euclidean_distance(const double* x1, const double* x2, Id n_features)
  {
    // TODO[space]: Default space is used here! Hopefully, it's RN (euclidean)... but n_features must be 2!
    thread_local SpacePoint p1;
    thread_local SpacePoint p2;
    if (p1.getSpace() != getDefaultSpaceSh())
    {
      p1.setSpace(getDefaultSpaceSh());
      p2.setSpace(getDefaultSpaceSh());
    }
    p1.setCoords(x1, n_features);
    p2.setCoords(x2, n_features);
    return p1.getDistance(p2);
  }
} // namespace gstlrn
