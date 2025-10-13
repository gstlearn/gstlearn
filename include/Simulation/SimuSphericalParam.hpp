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

#include "Basic/AStringable.hpp"

namespace gstlrn
{
class GSTLEARN_EXPORT SimuSphericalParam: public AStringable
{
public:
  SimuSphericalParam(Id special = 0,
                     Id nbf = 1,
                     Id nfmax = -1,
                     Id degmax = -1,
                     Id ndisc = 360,
                     double tol = 1.e-5);
  SimuSphericalParam(const SimuSphericalParam &r);
  SimuSphericalParam& operator=(const SimuSphericalParam &r);
  virtual ~SimuSphericalParam();

  /// Interface to AStringable
  String toString(const AStringFormat* strfmt = nullptr) const override;

  Id getNbf() const { return _nbf; }
  void setNbf(Id nbf) { _nbf = nbf; }
  Id getNdisc() const { return _ndisc; }
  void setNdisc(Id ndisc) { _ndisc = ndisc; }
  Id getNfmax() const { return _nfmax; }
  void setNfmax(Id nfmax) { _nfmax = nfmax; }
  Id getSpecial() const { return _special; }
  void setSpecial(Id special) { _special = special; }
  double getTol() const { return _tol; }
  void setTol(double tol) { _tol = tol; }
  Id getDegmax() const { return _degmax; }
  void setDegmax(Id degmax) { _degmax = degmax; }

private:
  Id _special; // 0: standard; 1 : Chentsov; 2 : Exponential
  Id _nbf;     // Number of basic functions
  Id _nfmax;   // Maximum number of frequencies (or <0)
  Id _degmax;  // Maximum Degree
  Id _ndisc;   // Number of discretization
  double _tol;  // Spectrum tolerance
};
}