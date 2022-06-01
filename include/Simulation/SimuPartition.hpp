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
#pragma once

#include "gstlearn_export.hpp"
#include "geoslib_define.h"

#include "Simulation/ASimulation.hpp"
#include "Basic/Plane.hpp"

class SimuPartitionParam;
class Db;
class DbGrid;

typedef struct
{
  double valref;
  double valsim;
} Stack;

class GSTLEARN_EXPORT SimuPartition: public ASimulation
{
public:
  SimuPartition(int nbsimu = 0, int seed = 4324324);
  SimuPartition(const SimuPartition &r) = delete;
  SimuPartition& operator=(const SimuPartition &r) = delete;
  virtual ~SimuPartition();

  int voronoi(DbGrid *dbgrid,
              Model *model,
              const SimuPartitionParam& parparam,
              int iptr,
              bool verbose = false);
  int poisson(DbGrid *dbgrid,
              Model *model,
              const SimuPartitionParam& parparam,
              int iptr,
              bool verbose = false);

private:
  double _stackSearch(const std::vector<Stack>& stacks, double valref);
};
