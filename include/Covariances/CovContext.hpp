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

#include "Basic/AStringable.hpp"
#include "Space/ASpace.hpp"

namespace gstlrn
{
  class Db;
  class Vario;

  class GSTLEARN_EXPORT CovContext: public AStringable
  {
  public:
    CovContext(Id nvar = 1, const ASpaceSharedPtr& space = ASpaceSharedPtr());
    CovContext(Id nvar, Id ndim, const VectorDouble& covar0 = VectorDouble());
    CovContext(const Vario* vario);
    // CovContext(const Vario* vario, const ASpaceSharedPtr& space = ASpaceSharedPtr());
    CovContext(const CovContext& r);
    CovContext& operator=(const CovContext& r);
    virtual ~CovContext();

    /// AStringable interface
    String toString(const AStringFormat* strfmt = nullptr) const override;

    /// Indicate if I am consistent with the provided space
    // bool isConsistent(const ASpace* space) const override;

    static CovContext* create(Id nvar, Id ndim);

    bool isEqual(const CovContext& r) const;

    Id getNVar() const { return _nVar; }

    double getField() const { return _field; }

    const VectorDouble& getCovar0() const { return _covar0; }

    double getCovar0(Id ivar, Id jvar) const;

    void setNVar(Id nvar)
    {
      _nVar = nvar;
      _update();
    }

    void setField(double field) { _field = field; }

    void setCovar0s(const VectorDouble& covar0);
    void setCovar0(Id ivar, Id jvar, double covar0);

    void copyCovContext(const CovContext& ctxt, bool severe = false);

    const CovContext* createReduce(const VectorInt& validVars) const;

    ASpaceSharedPtr getSpace() const { return _space; }

    size_t getNDim(Id ispace = -1) const { return _space->getNDim(ispace); }

    void setSpace(const ASpaceSharedPtr& space) { _space = space; }

  private:
    Id _nVar; /*! Number of variables */
    double _field; /*! Field maximum size */
    VectorDouble
      _covar0; /*! Variance-Covariance matrix (used for covariances) */
    mutable ASpaceSharedPtr _space;

  private:
    Id _getIndex(Id ivar, Id jvar) const;
    void _update();
  };
} // namespace gstlrn
