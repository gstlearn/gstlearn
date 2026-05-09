/*
                                      ball

Original Author: Eung Bum Lee
Website: https://42.fr
License: MIT
*/

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

/*
Modified by MINES Paris / ARMINES (2023)
Authors: gstlearn Team
Website: https://gstlearn.org
License: BSD 3-clause
*/

#pragma once

#include "Basic/VectorNumT.hpp"
#include "Matrix/MatrixT.hpp"

namespace gstlrn
{

  struct t_btree;

  struct t_nheap
  {
    MatrixT<double> distances;
    MatrixT<Id> indices;
    Id n_pts;
    Id n_nbrs;

    t_nheap() = default;
    void resize(Id n_pts, Id n_nbrs);

    /*
    ** neighbors_heap.c
    */
    double largest(Id row) const;
    Id push(Id row, double val, Id i_val);
    void sort();
    void load(const t_btree& b, const MatrixT<double>& x);
  };

  struct t_nodedata
  {
    Id idx_start;
    Id idx_end;
    double radius;
    bool is_leaf;
  };

  struct t_btree
  {
    MatrixT<double> data;
    VectorBool available;
    VectorInt idx_array;
    std::vector<t_nodedata> node_data;
    MatrixT<double> node_bounds;

    Id n_samples;
    Id n_features;

    Id leaf_size;
    Id n_levels;
    Id n_nodes;
    Id default_distance_function{1};

    t_btree() = default;
    t_btree(
      MatrixT<double>&& data,
      Id n_samples,
      Id n_features,
      bool all_available,
      Id leaf_size,
      Id default_distance_function);

    void display(Id level = -1) const;
    Id init_node(Id i_node, Id idx_start, Id idx_end);
    void recursive_build(Id i_node, Id idx_start, Id idx_end);
    double min_dist(Id i_node, const constvect pt) const;
    Id query_depth_first(
      Id i_node,
      const constvect pt,
      Id i_pt,
      t_nheap& heap,
      double dist) const;
    void query_radius_depth_first(
      Id i_node,
      const constvect pt,
      double radius,
      VectorInt& results) const;
  };

  /*
  ** metrics.c
  */
  double manhattan_distance(const double* x1, const double* x2, Id n_features);
  double euclidean_distance(const double* x1, const double* x2, Id n_features);

} // namespace gstlrn
