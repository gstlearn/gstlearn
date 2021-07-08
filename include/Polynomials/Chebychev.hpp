/*
 * Chebychev.hpp
 *
 *  Created on: 16 juin 2021
 *      Author: drenard
 */
#pragma once
#include "Polynomials/APolynomial.hpp"
#include "Basic/IClonable.hpp"
#include "geoslib_define.h"
#include "csparse_d.h"

#include <functional>

class Chebychev: public APolynomial
{
public:
  Chebychev();
  virtual ~Chebychev();
  IClonable* clone() const override {return new Chebychev(*this);}
  void init(int ncMax=10001,int nDisc=100,double a = 0.,double b=1.,bool verbose=false);
  int getNcMax() const {return _ncMax;}
  int getNDisc() const {return _nDisc;}
  bool getVerbose() const {return _verbose;}
  void setA(double a){_a=a;}
  void setB(double b){_b=b;}
  void setNcMax(int ncMax){_ncMax=ncMax;}
  void setNDisc(int nDisc){_nDisc=nDisc;}
  void setVerbose(bool verbose){_verbose = verbose;}

  void evalOp(cs* Op,const VectorDouble& in,VectorDouble& out,bool cumul=false) const override;
  double eval(double x) const override;
  int fit(std::function<double(double)> f,
          double a = 0.,
          double b = 1.,
          double tol = EPSILON5) override;

private:
  bool _isReady() const { return !_coeffs.empty(); }
  void _fillCoeffs(std::function<double(double)>, double a, double b);
  int _countCoeffs(std::function<double(double)> f,
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
