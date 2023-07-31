/* ************************************************************************** */
/*                                                                            */
/*                                                        :::      ::::::::   */
/*   metrics.c                                          :+:      :+:    :+:   */
/*                                                    +:+ +:+         +:+     */
/*   By: elee <elee@student.42.us.org>              +#+  +:+       +#+        */
/*                                                +#+#+#+#+#+   +#+           */
/*   Created: 2017/06/28 10:45:32 by elee              #+#    #+#             */
/*   Updated: 2017/06/28 16:52:59 by elee             ###   ########.fr       */
/*                                                                            */
/* ************************************************************************** */

#include "Tree/ball.h"

double manhattan_dist(double *x1, double *x2, int size)
{
	double	d = 0;

	for (int i = 0; i < size; i++)
		d += fabs(x1[i] - x2[i]);
	return (d);
}

double euclidean_dist(double *x1, double *x2, int size)
{
  double  d2 = 0;
  for (int i = 0; i < size; i++)
  {
    double delta = (x1[i] - x2[i]);
    d2 += delta * delta;
  }
  return (sqrt(d2));
}

double min_dist(t_btree *tree, int i_node, double *pt)
{
	double	dist_pt;

	dist_pt = manhattan_dist(pt, tree->node_bounds[0][i_node], tree->n_features);
	return (fmax(0.0, dist_pt - tree->node_data[i_node].radius));
}
