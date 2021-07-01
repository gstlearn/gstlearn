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

#include "Basic/Vector.hpp"
#include "Space/ASpaceObject.hpp"
#include "Space/SpacePoint.hpp"

class Db;

class ADrift : public ASpaceObject
{
public:
  ADrift(const ASpace* space = nullptr);
  ADrift(const ADrift &r);
  ADrift& operator=(const ADrift &r);
  virtual ~ADrift();

  virtual int getNVariables() const = 0;
  virtual double eval(const Db* db, int iech) const = 0;
};
