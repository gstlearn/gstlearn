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

  class GSTLEARN_EXPORT GibbsMultiMono: public AGibbs
  {
  public:
    GibbsMultiMono();
    GibbsMultiMono(Db* db, const std::vector<Model*>& models, double rho);
    GibbsMultiMono(const GibbsMultiMono& r);
    GibbsMultiMono& operator=(const GibbsMultiMono& r);
    virtual ~GibbsMultiMono();

    Model* getModels(Id ivar) const
    {
      return _models[ivar];
    } // TODO: protect by const asap

    double getRho() const { return _rho; }

    Id getNVar() const { return static_cast<Id>(_models.size()); }

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

  private:
    std::vector<Model*> _models;
    double _rho;
  };
} // namespace gstlrn
