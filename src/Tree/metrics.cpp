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
#include "Space/SpacePoint.hpp"

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
double euclidean_dist(const double *x1, const double *x2, int size)
{
  SpacePoint p1;
  SpacePoint p2;
  for (int i = 0; i < size; i++)
  {
    p1.setCoord(i, x1[i]);
    p2.setCoord(i, x2[i]);
  }
  return p1.getDistance(p2);
}
