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

#include "Gibbs/AGibbs.hpp"
#include "gstlearn_export.hpp"

namespace gstlrn
{
class Db;
class Model;

class GSTLEARN_EXPORT GibbsMulti: public AGibbs
{
public:
  GibbsMulti();
  GibbsMulti(Db* db, Model* model);
  GibbsMulti(const GibbsMulti& r);
  GibbsMulti& operator=(const GibbsMulti& r);
  virtual ~GibbsMulti();

  /// Interface for AGibbs
  Id calculInitialize(VectorVectorDouble& y, Id isimu, Id ipgs) override;
  double getSimulate(VectorVectorDouble& y,
                     double yk,
                     double sk,
                     Id icase,
                     Id ipgs,
                     Id ivar,
                     Id iact,
                     Id iter) override;
  Id checkGibbs(const VectorVectorDouble& y, Id isimu, Id ipgs) override;

  Model* getModel() const { return _model; } // protect using const asap

private:
  Model* _model;
};
} // namespace gstlrn