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
#include "PluriGaussian/CorPGS.hpp"
#include "Enum/EOperator.hpp"
#include "Matrix/EigenVectors.hpp"
#include "geoslib_define.h"

namespace gstlrn
{
  CorPGS::CorPGS(Id option, bool flag_rho, double rho)
    : _optCorrel(option)
    , _npar(0)
    , _flagRho(flag_rho)
    , _rho(rho)
    , _params()
    , _modif()
  {
    _params.resize(4, 0.);
    _modif.resetFromValue(4, 4, 0.);
  }

  void CorPGS::define(Id option, bool flag_rho, double rho)
  {
    _optCorrel = option;
    _flagRho = flag_rho;
    _rho = rho;
    defineModifMatrix();
  }

  /****************************************************************************/
  /*!
   **  Establish the total vector C1(h), C12(h), C21(h), C2(h)
   ** \param[in]  params_in   Parameters (Dimension corpgs.npar)
   **
   ** \return params      Parameters (Dimension = 4)
   **
   *****************************************************************************/
  VectorDouble CorPGS::_computeParams(const VectorDouble& params_in) const
  {
    VectorDouble params(4);

    switch (_optCorrel)
    {
      case 0:
        params[0] = params_in[0];
        params[1] = params_in[1];
        params[2] = params_in[2];
        params[3] = params_in[3];
        break;

      case 1: /* Symmetrical case */
        params[0] = params_in[0];
        params[1] = params[2] = params_in[1];
        params[3] = params_in[2];
        break;

      case 2: /* Residual case */
        double rho2 = _rho * _rho;
        params[0] = params_in[0];
        params[1] = params[2] = _rho * params_in[0];
        params[3] = rho2 * params_in[0] + (1. - rho2) * params_in[1];
        break;
    }
    return params;
  }

  /****************************************************************************/
  /*!
   **  Establish the correlation for C1(h), C12(h), C21(h), C2(h)
   **
   ** \param[in]  params_in   Array of parameters
   **
   ** \param[out] correl      Correlation matrix (Dimension = 4*4)
   **
   *****************************************************************************/
  void CorPGS::buildCorrel(
    const VectorDouble& params_in,
    MatrixSymmetric& correl) const
  {
    VectorDouble params = _computeParams(params_in);

    correl.fill(0.);

    correl.setValue(0, 0, 1.);
    correl.setValue(1, 1, 1.);
    correl.setValue(2, 2, 1.);
    correl.setValue(3, 3, 1.);

    correl.setValue(2, 0, _rho);
    correl.setValue(3, 1, _rho);

    correl.setValue(1, 0, params[0]);
    correl.setValue(2, 1, params[2]);
    correl.setValue(3, 0, params[1]);
    correl.setValue(3, 2, params[3]);
  }

  /****************************************************************************/
  /*!
   **  Update the following matrices according to constraints on model
   **
   ** \param[in,out]  Grad        Vector of gradients (Dimension = npar)
   ** \param[in,out]  Hess        Matrix of Hessian (Dimension = npar * npar)
   **
   *****************************************************************************/
  void CorPGS::_updateConstraints(VectorDouble& Grad, MatrixSymmetric& Hess)
  {
    /* Update the Grad */

    VectorDouble v = Grad;
    for (Id i = 0; i < _npar; i++)
    {
      double value = 0.;
      for (Id j = 0; j < 4; j++) value += v[j] * _modif.getValue(i, j);
      Grad[i] = value;
    }

    /* Update the Hessian */

    MatrixSymmetric m = Hess;
    Hess.fill(0.);
    for (Id i = 0; i < _npar; i++)
      for (Id j = 0; j < _npar; j++)
        for (Id k = 0; k < 4; k++)
          for (Id l = 0; l < 4; l++)
            Hess.updValue(
              i, j, EOperator::ADD,
              _modif.getValue(i, k) * m.getValue(k, l) * _modif.getValue(j, l));
  }

  /****************************************************************************/
  /*!
   **  Update the following matrices according to constraints on model
   **
   ** \param[in,out]  Grad        Vector of gradients (Dimension = npar)
   ** \param[in,out]  Hess        Matrix of Hessian (Dimension = npar * npar)
   ** \param[in,out]  JJ          Matrix of t(JJ) * JJ (Dimension = npar * npar)
   **
   *****************************************************************************/
  void CorPGS::updateConstraintsWithJJ(
    VectorDouble& Grad,
    MatrixSymmetric& Hess,
    MatrixSymmetric& JJ)
  {

    /* Update the Grad and Hessian */

    _updateConstraints(Grad, Hess);

    /* Update JJ */

    MatrixSymmetric m = JJ;
    JJ.fill(0.);
    for (Id i = 0; i < _npar; i++)
      for (Id j = 0; j < _npar; j++)
        for (Id k = 0; k < 4; k++)
          for (Id l = 0; l < 4; l++)
            JJ.updValue(
              i, j, EOperator::ADD,
              _modif.getValue(i, k) * m.getValue(k, l) * _modif.getValue(j, l));
  }

  /****************************************************************************/
  /*!
   **  Calculate the generalized inverse of a square symmetric matrix
   **
   ** \return  Error returned code
   **
   ** \param[in]  a         Matrix to be inverted
   **
   ** \param[out] tabout    Inverted matrix
   **
   *****************************************************************************/
  Id CorPGS::_invgen(MatrixSymmetric& a, MatrixSymmetric& tabout)
  {
    auto neq = a.getNRows();
    tabout.fill(0.);

    /* Calculate the eigen vectors */

    auto eigenvectors = EigenVectors(a);
    const auto& eigval = eigenvectors.getEigenValues();
    const MatrixSquare* eigvec = &eigenvectors.getEigenVectors();

    /* Calculate the generalized inverse */

    for (Id i = 0; i < neq; i++)
      for (Id j = 0; j < neq; j++)
      {
        double value = 0.;
        for (Id k = 0; k < neq; k++)
        {
          if (ABS(eigval[k]) > 1e-10)
            value +=
              eigvec->getValue(k, i) * eigvec->getValue(k, j) / eigval[k];
        }
        tabout.setValue(i, j, value);
      }
    return 0;
  }

  /****************************************************************************/
  /*!
   ** Calculate the indexes of each parameter
   **
   ** \param[in]  i Index
   ** \param[in]  j Index
   **
   *****************************************************************************/
  Id CorPGS::_F(Id i, Id j)
  {
    Id value = 0;
    switch (i)
    {
      case 0: value = (j == 0) ? 1 : 0; break;

      case 1: value = (j == 0) ? 3 : 0; break;

      case 2: value = (j == 0) ? 2 : 1; break;

      case 3: value = (j == 0) ? 3 : 2; break;
    }
    return (value);
  }

  /****************************************************************************/
  /*!
   **  Compute the derivatives (first and second) of the smallest eigenvalue
   **
   ** \param[in] eigval  Current eigen value
   **
   ** \param[out] ev     Output array
   ** \param[out] d1     First order derivative
   ** \param[out] d2     Second order derivative
   **
   *****************************************************************************/
  void CorPGS::derivativeEigen(
    double eigval,
    const MatrixSquare* ev,
    VectorDouble& d1,
    MatrixSymmetric& d2)
  {
    MatrixSymmetric temp(4);
    temp.fill(0.);
    MatrixSymmetric invGn(4);
    d2.fill(0.);
    buildCorrel(_params, temp);
    AMatrix::linearCombinationInPlace(temp, 0., -1., temp);

    for (Id i = 0; i < 4; i++) temp.updValue(i, i, EOperator::ADD, eigval);

    _invgen(temp, invGn);

    for (Id i = 0; i < 4; i++)
      d1[i] = 2 * ev->getValue(3, _F(i, 0)) * ev->getValue(3, _F(i, 1));

    for (Id i = 0; i < 4; i++)
      for (Id j = 0; j < i; j++)
      {
        double value = 0;
        value = ev->getValue(3, _F(i, 0)) * ev->getValue(3, _F(j, 0))
              * invGn.getValue(_F(i, 1), _F(j, 1));
        value += ev->getValue(3, _F(i, 1)) * ev->getValue(3, _F(j, 0))
               * invGn.getValue(_F(i, 0), _F(j, 1));
        value += ev->getValue(3, _F(i, 0)) * ev->getValue(3, _F(j, 1))
               * invGn.getValue(_F(i, 1), _F(j, 0));
        value += ev->getValue(3, _F(i, 1)) * ev->getValue(3, _F(j, 1))
               * invGn.getValue(_F(i, 0), _F(j, 0));
        d2.setValue(i, j, 2 * value);
      }

    _updateConstraints(d1, d2);
  }

  /****************************************************************************/
  /*!
   **  Expand the vector of parameters into C1, C12, C21 and C2
   **  according to the constraints
   **
   ** \return  Returned parameter
   **
   ** \param[in]  igrf       Rank of the first variable
   ** \param[in]  jgrf       Rank of the second variable
   ** \param[in]  idir       positive (1) or negative (-1) distance
   **
   *****************************************************************************/
  double CorPGS::paramExpand(Id igrf, Id jgrf, Id idir)
  {
    double rho2;

    switch (_optCorrel)
    {
      case 0:
        if (igrf == 0 && jgrf == 0)
          return (_params[0]);
        else if (igrf == 1 && jgrf == 1)
          return (_params[3]);
        else
        {
          if (idir > 0) return (_params[1]);
          return (_params[2]);
        }
        break;

      case 1:
        if (igrf == 0 && jgrf == 0) return (_params[0]);
        if (igrf == 1 && jgrf == 1) return (_params[2]);
        return (_params[1]);
        break;

      case 2:
        rho2 = _rho * _rho;
        if (igrf == 0 && jgrf == 0) return (_params[0]);
        if (igrf == 1 && jgrf == 1)
          return (_params[0] * rho2 + _params[1] * (1. - rho2));
        return (_params[0] * _rho);
        break;
    }
    return (0.);
  }

  /****************************************************************************/
  /*!
   **  Compute the modif matrix
   **
   *****************************************************************************/
  void CorPGS::defineModifMatrix()
  {
    _modif.fill(0.);

    switch (_optCorrel)
    {
      case 0: /* Full parameters */
        _npar = 4;
        _modif.setValue(0, 0, 1);
        _modif.setValue(1, 1, 1);
        _modif.setValue(2, 2, 1);
        _modif.setValue(3, 3, 1);
        break;

      case 1: /* Symmetrical case */
        _npar = 3;
        _modif.setValue(0, 0, 1);
        _modif.setValue(1, 1, 1);
        _modif.setValue(1, 2, 1);
        _modif.setValue(2, 3, 1);
        break;

      case 2: /* Residual case */
        double rho2 = _rho * _rho;
        _npar = 2;
        _modif.setValue(0, 0, 1);
        _modif.setValue(0, 1, _rho);
        _modif.setValue(0, 2, _rho);
        _modif.setValue(0, 3, rho2);
        _modif.setValue(1, 3, 1. - rho2);
        break;
    }
  }

  /****************************************************************************/
  /*!
   **  Initialize the parameters
   **
   *****************************************************************************/
  void CorPGS::initializeParams()
  {
    switch (_optCorrel)
    {
      case 0:
        _params[0] = _params[3] = fabs(_rho);
        _params[1] = _params[2] = fabs(_rho) * _rho;
        break;

      case 1:
        _params[0] = _params[2] = fabs(_rho);
        _params[1] = fabs(_rho) * _rho;
        break;

      case 2: _params[0] = _params[1] = fabs(_rho); break;
    }
  }

  /****************************************************************************/
  /*!
   **  Set the model-type (opt_correl)
   **
   ** \param[in]  opt     The model-type to set
   **
   ****************************************************************************/
  void CorPGS::setOptCorrel(Id opt)
  {
    VectorDouble params = _computeParams(_params);

    switch (opt)
    {
      case 0:
        _params[0] = params[0];
        _params[1] = params[1];
        _params[2] = params[2];
        _params[3] = params[3];
        break;

      case 1:
        _params[0] = params[0];
        _params[1] = _params[2] = (params[1] + params[2]) / 2.;
        _params[2] = params[3];
        break;

      case 2:
        _params[0] = params[0];
        _params[1] = params[3];
        break;
    }
    _optCorrel = opt;
    defineModifMatrix();
  }

} // namespace gstlrn
