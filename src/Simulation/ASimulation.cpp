/******************************************************************************/
/* COPYRIGHT ARMINES, ALL RIGHTS RESERVED                                     */
/*                                                                            */
/* THE CONTENT OF THIS WORK CONTAINS CONFIDENTIAL AND PROPRIETARY             */
/* INFORMATION OF ARMINES. ANY DUPLICATION, MODIFICATION,                     */
/* DISTRIBUTION, OR DISCLOSURE IN ANY FORM, IN WHOLE, OR IN PART, IS STRICTLY */
/* PROHIBITED WITHOUT THE PRIOR EXPRESS WRITTEN PERMISSION OF ARMINES         */
/*                                                                            */
/* TAG_SOURCE_CG                                                              */
/******************************************************************************/
#include "Simulation/ASimulation.hpp"


ASimulation::ASimulation(int nbsimu, int seed)
    : _nbsimu(nbsimu),
      _seed(seed)
{
}

ASimulation::ASimulation(const ASimulation &r)
  : _nbsimu(r._nbsimu),
    _seed(r._seed)
{
}

ASimulation& ASimulation::operator=(const ASimulation &r)
{
  if (this != &r)
  {
    _nbsimu = r._nbsimu;
    _seed = r._seed;
  }
  return *this;
}

ASimulation::~ASimulation()
{
}
