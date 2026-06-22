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
#include "Matrix/MatrixSymmetric.hpp"

namespace gstlrn
{
  class GSTLEARN_EXPORT CorPGS
  {
  public:
    CorPGS(Id option = 0, bool flag_rho = false, double rho = 0.);
    CorPGS(const CorPGS& r) = default;
    CorPGS& operator=(const CorPGS& r) = default;
    virtual ~CorPGS() = default;

    void define(Id option, bool flag_rho, double rho);
    void
      buildCorrel(const VectorDouble& params_in, MatrixSymmetric& correl) const;
    void updateConstraintsWithJJ(
      VectorDouble& Grad,
      MatrixSymmetric& Hess,
      MatrixSymmetric& JJ);
    void derivativeEigen(
      double eigval,
      const MatrixSquare* ev,
      VectorDouble& d1,
      MatrixSymmetric& d2);
    double paramExpand(Id igrf, Id jgrf, Id idir);
    void defineModifMatrix();
    void initializeParams();
    void setOptCorrel(Id opt);

    Id getNpar() const { return _npar; }

    double getRho() const { return _rho; }

    Id getOptCorrel() const { return _optCorrel; }

    VectorDouble getParams() const { return _params; }

    double getParam(Id i) const { return _params[i]; }

    void setRho(double rho) { _rho = rho; }

    void setParams(const VectorDouble& params) { _params = params; }

    void setModif(Id i, Id j, double value) { _modif.setValue(i, j, value); }

  private:
    static Id _F(Id i, Id j);
    static Id _invgen(MatrixSymmetric& a, MatrixSymmetric& tabout);
    VectorDouble _computeParams(const VectorDouble& params_in) const;
    void _updateConstraints(VectorDouble& Grad, MatrixSymmetric& Hess);

  private:
    Id _optCorrel;
    Id _npar;
    bool _flagRho;
    double _rho;
    VectorDouble _params;
    MatrixSymmetric _modif;
  };
} // namespace gstlrn
