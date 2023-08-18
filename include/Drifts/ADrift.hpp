/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://gstlearn.org                                              */
/* License: BSD 3-clause                                                      */
/*                                                                            */
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
