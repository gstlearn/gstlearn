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
#pragma once

#include "gstlearn_export.hpp"

#include "Mesh/AMesh.hpp"
#include "Basic/AFunctional.hpp"
#include "Model/ANoStat.hpp"

/**
 * This class concerns the non-stationarity defined as a function (hence its name).
 * It can be considered as an example of a 2-D implementation of a spiral with
 * a single non-stationary parameter (Example: the angle)
 */
class GSTLEARN_EXPORT NoStatFunctional : public ANoStat
{
public:
	NoStatFunctional(const VectorString& code = {"A"});
	NoStatFunctional(const AFunctional* func, const VectorString& code = {"A"});
  NoStatFunctional(const NoStatFunctional &m);
  NoStatFunctional& operator=(const NoStatFunctional &m);
  virtual ~NoStatFunctional();

  /// ICloneable interface
  IMPLEMENT_CLONING(NoStatFunctional)

  /// AStringable Interface
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  int  attachToMesh(const AMesh* mesh, bool center = true, bool verbose = false) const override;
  int  attachToDb(Db* db, int icas, bool verbose = false) const override;

  double getValue(const EConsElem &type,
                  int icas,
                  int rank,
                  int icov = 0,
                  int iv1 = -1,
                  int iv2 = -1,
                  int igrf = -1) const override;
  double getValueByParam(int ipar, int icas, int rank) const override;

private:
  const AFunctional* _func;
};
