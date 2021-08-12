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

class Db;

class AFunctional
{
public:
  AFunctional(int ndim);
  AFunctional(const AFunctional &m);
  AFunctional& operator=(const AFunctional &m);
  virtual ~AFunctional();

  virtual double getFunctionValue(const VectorDouble& pos) const = 0;

  int  getNdim() const { return _ndim; }
  void setNdim(int ndim) { _ndim = ndim; }

  VectorDouble getFunctionValues(const Db *db, bool useSel = true) const;

private:
  int _ndim;
};
