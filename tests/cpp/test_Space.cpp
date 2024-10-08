/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES Paris / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://gstlearn.org                                              */
/* License: BSD 3-clause                                                      */
/*                                                                            */
/******************************************************************************/
#include "Basic/File.hpp"
#include "Space/ASpaceObject.hpp"
#include "Space/SpaceComposite.hpp"
#include "Space/SpaceRN.hpp"
#include "Space/SpacePoint.hpp"
#include <iostream>

/**
 * Demonstration of space composite features
 */
int main(int argc, char *argv[])
{
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str(), argc, argv);

  // 2D Space + Time (for example)
  SpaceComposite sp({new SpaceRN(2), new SpaceRN(1)});
  defineDefaultSpace(&sp);

  const ASpace* psp = getDefaultSpace();
  psp->display();  // Long description (level 1)
  psp->display(0); // Short description (level 0)

  SpacePoint pt1({4.5, 6.5, 10});
  SpacePoint pt2({3.0, 2.5, 15});

  std::cout << "Global dimension: " << pt1.getNDim() << std::endl;
  std::cout << "Dimension space #0: " << pt1.getNDim(0) << std::endl;
  std::cout << "Dimension space #1: " << pt1.getNDim(1) << std::endl;

  // Why this works only under linux ?
  // std::cout << "Global distances: " << pt1.getDistances(pt2) << std::endl;
  std::cout << "Global distances: " << pt1.getDistances(pt2).toString() << std::endl;

  std::cout << "Distance space #0: " << pt1.getDistance(pt2, 0) << std::endl;
  std::cout << "Distance space #1: " << pt1.getDistance(pt2, 1) << std::endl;

  std::cout << "Increments space #0: " << pt1.getIncrement(pt2, 0).toString() << std::endl;
  std::cout << "Increments space #1: " << pt1.getIncrement(pt2, 1).toString() << std::endl;

  return 0;
}
