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

#include "Basic/AFunctional.hpp"
#include "Covariances/ANoStat.hpp"

/**
 * This class concerns the non-stationarity defined as a function (hence its name).
 */
class GSTLEARN_EXPORT NoStatFunctional : public ANoStat
{
public:
	NoStatFunctional(const AFunctional* func);
  NoStatFunctional(const NoStatFunctional &m) = delete;
  NoStatFunctional& operator=(const NoStatFunctional &m) = delete;
  virtual ~NoStatFunctional();
  String toString(const AStringFormat* strfmt) const override;

private :
  void _informField(const VectorVectorDouble& coords,
                    VectorDouble& tab,
                    bool verbose = false) override; 

private:
  const AFunctional* _func;
};
