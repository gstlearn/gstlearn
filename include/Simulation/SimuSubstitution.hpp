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

class SimuSubstitutionParam;
class Db;
class DbGrid;

class GSTLEARN_EXPORT SimuSubstitution: public ASimulation
{
public:
  SimuSubstitution(int nbsimu = 0, int seed = 4324324);
  SimuSubstitution(const SimuSubstitution &r);
  SimuSubstitution& operator=(const SimuSubstitution &r);
  virtual ~SimuSubstitution();

  int simulate(DbGrid *dbgrid,
               const SimuSubstitutionParam& subparam,
               int iptr,
               int verbose);

private:
  void _calculValue(int ip, double factor, const VectorDouble& vector);
  VectorDouble _transToProp(const SimuSubstitutionParam& subparam,
                            int verbose,
                            double eps = EPSILON5);

private:
  std::vector<Plane> _planes;
};
