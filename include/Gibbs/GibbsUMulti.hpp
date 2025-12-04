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

#include "Gibbs/GibbsMulti.hpp"
#include "Matrix/MatrixSquare.hpp"

namespace gstlrn
{
class Db;
class Model;

class GSTLEARN_EXPORT GibbsUMulti: public GibbsMulti
{
public:
  GibbsUMulti();
  GibbsUMulti(Db* db, Model* model);
  GibbsUMulti(const GibbsUMulti& r);
  GibbsUMulti& operator=(const GibbsUMulti& r);
  virtual ~GibbsUMulti();

  void update(VectorVectorDouble& y, Id isimu, Id ipgs, Id iter) override;
  Id covmatAlloc(bool verbose, bool verboseTimer = false) override;

private:
  Id _getSize() const;
  double _getVariance(Id iecr) const;
  double _getEstimate(Id ipgs, Id iecr, VectorVectorDouble& y);

private:
  MatrixSquare _covmat;
};
} // namespace gstlrn