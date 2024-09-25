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

#include "LinearOp/ALinearOp.hpp"
#include "Polynomials/APolynomial.hpp"
#include "Basic/ICloneable.hpp"
#include "geoslib_define.h"

#include <functional>


class AFunction;
class ALinearOpMulti;
class cs; /// TODO : Dependency to csparse to be removed

class GSTLEARN_EXPORT Chebychev: public APolynomial
{
public:
  Chebychev();
  virtual ~Chebychev();

  /// ICloneable interface
  IMPLEMENT_CLONING(Chebychev)

  /// Interface for Apolynomial
#ifndef SWIG
  void evalOp(MatrixSparse* S, const constvect x, vect y) const override;
  void addEvalOp(ALinearOp* Op, const constvect inv, vect outv) const override;

  /* void evalOp(const ALinearOpMulti *Op,
              const std::vector<Eigen::VectorXd> &inv,
              std::vector<Eigen::VectorXd> &outv) const override; */
 #endif
  double eval(double x) const override;
  int fit(const std::function<double(double)>& f,
          double a = 0.,
          double b = 1.,
          double tol = EPSILON5) override;

  void init(int ncMax=10001,int nDisc=100,double a = 0.,double b=1.,bool verbose=false);
  static Chebychev* createFromCoeffs(const VectorDouble& coeffs);
  void setCoeffs(const VectorDouble& coeffs){_coeffs = coeffs;}
  int getNcMax() const {return _ncMax;}
  int getNDisc() const {return _nDisc;}
  double getA() const {return _a;}
  double getB() const {return _b;}
  bool getVerbose() const {return _verbose;}
  void setA(double a){_a=a;}
  void setB(double b){_b=b;}
  void setNcMax(int ncMax){_ncMax=ncMax;}
  void setNDisc(int nDisc){_nDisc=nDisc;}
  void setVerbose(bool verbose){_verbose = verbose;}

  int fit2(AFunction *f, double a = 0., double b = 1., double tol = EPSILON5);

private:
  bool _isReady() const { return !_coeffs.empty(); }
  void _fillCoeffs(const std::function<double(double)>& f, double a, double b);
  int _countCoeffs(const std::function<double(double)>& f,
                   double x,
                   double a,
                   double b,
                   double tol = EPSILON5) const;

  int _ncMax ;
  int _nDisc ;
  double _a;
  double _b;
  bool _verbose;
};
