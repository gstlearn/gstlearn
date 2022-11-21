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

#include "Space/ASpaceObject.hpp"
#include "Space/SpacePoint.hpp"

class Db;

class GSTLEARN_EXPORT ADrift : public ASpaceObject
{
public:
  ADrift(const ASpace* space = nullptr);
  ADrift(const ADrift &r);
  ADrift& operator=(const ADrift &r);
  virtual ~ADrift();

  virtual int getNVariables() const = 0;
  /// TODO : Change ADrift::eval args from {Db,iech} to {SpacePoint p1, SpacePoint p2}
  virtual double eval(const Db* db, int iech) const = 0;
};
