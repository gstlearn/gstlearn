/* ************************************************************************** */
/*                                                                            */
/*                                                        :::      ::::::::   */
/*   main.c                                             :+:      :+:    :+:   */
/*                                                    +:+ +:+         +:+     */
/*   By: elee <elee@student.42.us.org>              +#+  +:+       +#+        */
/*                                                +#+#+#+#+#+   +#+           */
/*   Created: 2017/06/28 17:45:08 by elee              #+#    #+#             */
/*   Updated: 2017/06/28 21:53:50 by elee             ###   ########.fr       */
/*                                                                            */
/* ************************************************************************** */
#include "geoslib_define.h"

#include "Db/Db.hpp"
#include "Basic/File.hpp"
#include "Basic/Timer.hpp"
#include "Basic/AStringable.hpp"
#include "Tree/Ball.hpp"
#include "Tree/ball.h"

int main(int argc, char *argv[])
{
	t_knn	knn;
  Timer timer;

  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str(), argc, argv);

	// Constructing a Data Db
	int ndim = 2;
	int ndata = 1000;
	Db* data = Db::createFillRandom(ndata, ndim, 0, 0, 0., 0.,
                                  VectorDouble(), VectorDouble(), VectorDouble(),
                                  2451);

	// Constructing a Target Db
	int ntarget = 1000;
	Db* target = Db::createFillRandom(ntarget, ndim, 0, 0, 0., 0.,
	                                  VectorDouble(), VectorDouble(), VectorDouble(),
	                                  131323);

	// Parameters
  int leaf_size = 10;
  int dist_type = 0;
  int neigh_size = 10;
  mestitle(1, "Testing BallTree performance");
  message(" - Data Db contains %d samples\n", ndata);
  message(" - Testing Db contains %d targets\n", ntarget);
  message(" - Space Dimension = %d\n", ndim);
  message(" - Leaf size is %d\n", leaf_size);
  message(" - Neighborhood size = %d\n", neigh_size);
  message(" - Distance type = %d\n", dist_type);

  // Constructing the Ball
  timer.reset();
  Ball ball(data, leaf_size, dist_type);
  timer.displayIntervalMilliseconds("Instantiating BallTree class", 0);

  timer.reset();
  (void) ball.build();
  timer.displayIntervalMilliseconds("Building BallTree",0);

  timer.reset();
  VectorDouble coor(ndim);
	for (int is = 0, ns = target->getSampleNumber(); is < ns; is++)
	{
	  target->getCoordinatesPerSampleInPlace(is, coor);
	  knn = ball.queryOne(coor, neigh_size);
	  free_knn(knn, 1);
	}
  timer.displayIntervalMilliseconds("Querying the Ball Tree Per Target",0);

  timer.reset();
  VectorVectorDouble coors = target->getAllCoordinates();
  knn = ball.queryAsVVD(coors, neigh_size);
  timer.displayIntervalMilliseconds("Querying the Ball Tree for all Targets at once",0);
  free_knn(knn, ntarget);

	//free stuff
  delete data;
  delete target;
	return (0);
}
