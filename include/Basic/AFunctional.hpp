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
#include "Basic/VectorNumT.hpp"

namespace gstlrn
{ 
class Db;

class GSTLEARN_EXPORT AFunctional
{
public:
  AFunctional(Id ndim);
  AFunctional(const AFunctional &m);
  AFunctional& operator=(const AFunctional &m);
  virtual ~AFunctional();

  virtual double getFunctionValue(const VectorDouble& pos) const = 0;

  Id  getNdim() const { return _ndim; }
  void setNdim(Id ndim) { _ndim = ndim; }

  VectorDouble getFunctionValues(const Db *db, bool useSel = true) const;

private:
  Id _ndim;
};
}