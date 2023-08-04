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
#include "Tree/ball_algorithm.h"

/**
 * Returns the Manhattan distance between two points
 * @param x1 Vector of coordinates for the first point
 * @param x2 Vector of coordinates for the second point
 * @param size Number of coordinates
 * @return
 */
double manhattan_dist(const double *x1, const double *x2, int size)
{
  double delta;
	double d1 = 0;
	for (int i = 0; i < size; i++)
	{
	  delta = fabs(x1[i] - x2[i]);
		d1 += delta;
	}
	return (d1);
}

/**
 * Returns the Euclidean distance between two points
 * @param x1 Vector of coordinates for the first point
 * @param x2 Vector of coordinates for the second point
 * @param size Number of coordinates
 * @return
 */
double euclidean_dist(const double *x1, const double *x2, int size)
{
  double delta;
  double d2 = 0;
  for (int i = 0; i < size; i++)
  {
    delta = (x1[i] - x2[i]);
    d2 += delta * delta;
  }
  return sqrt(d2);
}
