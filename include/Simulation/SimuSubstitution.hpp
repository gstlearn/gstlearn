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
#include "Basic/AStringable.hpp"

class SimuSubstitutionParam;
class Db;
class DbGrid;

class GSTLEARN_EXPORT SimuSubstitution: public ASimulation, public AStringable
{
public:
  SimuSubstitution(int nbsimu = 0, int seed = 4324324);
  SimuSubstitution(const SimuSubstitution &r);
  SimuSubstitution& operator=(const SimuSubstitution &r);
  virtual ~SimuSubstitution();

  /// Interface to AStringable
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  int simulate(Db *dbin,
               DbGrid *dbout,
               const SimuSubstitutionParam& subparam,
               int iptr_simu,
               int iptr_rank,
               bool verbose = false);
};
