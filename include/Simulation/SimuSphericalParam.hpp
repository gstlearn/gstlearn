/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://github.com/gstlearn                                       */
/* License: BSD 3 clauses                                                     */
/*                                                                            */
/******************************************************************************/
#pragma once

#include "gstlearn_export.hpp"

#include "Basic/AStringable.hpp"
#include "Basic/VectorNumT.hpp"

class GSTLEARN_EXPORT SimuSphericalParam: public AStringable
{
public:
  SimuSphericalParam(int special = 0,
                     int nbf = 1,
                     int nfmax = -1,
                     int degmax = -1,
                     int ndisc = 360,
                     double tol = 1.e-5);
  SimuSphericalParam(const SimuSphericalParam &r);
  SimuSphericalParam& operator=(const SimuSphericalParam &r);
  virtual ~SimuSphericalParam();

  /// Interface to AStringable
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  int getNbf() const { return _nbf; }
  void setNbf(int nbf) { _nbf = nbf; }
  int getNdisc() const { return _ndisc; }
  void setNdisc(int ndisc) { _ndisc = ndisc; }
  int getNfmax() const { return _nfmax; }
  void setNfmax(int nfmax) { _nfmax = nfmax; }
  int getSpecial() const { return _special; }
  void setSpecial(int special) { _special = special; }
  double getTol() const { return _tol; }
  void setTol(double tol) { _tol = tol; }
  int getDegmax() const { return _degmax; }
  void setDegmax(int degmax) { _degmax = degmax; }

private:
  int _special; // 0: standard; 1 : Chentsov; 2 : Exponential
  int _nbf;     // Number of basic functions
  int _nfmax;   // Maximum number of frequencies (or <0)
  int _degmax;  // Maximum Degree
  int _ndisc;   // Number of discretization
  double _tol;  // Spectrum tolerance
};
