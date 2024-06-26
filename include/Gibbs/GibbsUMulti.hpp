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

class Db;
class Model;

class GSTLEARN_EXPORT GibbsUMulti : public GibbsMulti
{
public:
  GibbsUMulti();
  GibbsUMulti(Db* db, Model* model);
  GibbsUMulti(const GibbsUMulti &r);
  GibbsUMulti& operator=(const GibbsUMulti &r);
  virtual ~GibbsUMulti();

  void update(VectorVectorDouble &y, int isimu, int ipgs, int iter) override;
  int covmatAlloc(bool verbose, bool verboseTimer = false) override;

private:
  int    _getSize() const;
  double _getVariance(int iecr) const;
  double _getEstimate(int ipgs, int iecr, VectorVectorDouble& y);

private:
  VectorDouble _covmat;
};
