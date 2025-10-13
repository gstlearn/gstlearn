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
#include "LinearOp/CGParam.hpp"

namespace gstlrn
{
class ALinearOp;

/**
 * @brief Definition of the parameters used within SPDE
 *
 * refineK Discretization factor used for Kriging
 * refineS Discretization factor used for Simulation
 * border  Border size
 * nxmax   Maximum number of vertices in the internal mesh (0 : no limit)
 * epsNugget Nugget effect
 * useStencil Default option for no Cholesky (can only be used in stationary case for Turbo Meshing)
 * nMC Number of Monte-Carlo simulations (used for Variance and logdet)
 * seedMC Seed for the random number generator (used for Variance and logdet)
 *
 * cgparams Parameters for the Conjugate Gradient method
 */
class GSTLEARN_EXPORT SPDEParam: public AStringable
{

public:
  SPDEParam(Id refineK             = 11,
            Id refineS             = 18,
            Id border              = 8,
            bool flag_polarized     = true,
            Id nxmax               = 300,
            double epsNugget        = EPSILON2,
            bool useStencil         = true,
            Id nMC                 = 10,
            Id seedMC              = 134341,
            const CGParam& cgparams = CGParam());
  SPDEParam(const SPDEParam& m);
  SPDEParam& operator=(const SPDEParam& m);
  virtual ~SPDEParam();

  /// AStringable Interface
  String toString(const AStringFormat* strfmt = nullptr) const override;

  static SPDEParam* create(Id refineK             = 11,
                           Id refineS             = 18,
                           Id border              = 8,
                           bool flag_polarized     = true,
                           Id nxmax               = 300,
                           double epsNugget        = EPSILON2,
                           bool useStencil         = true,
                           Id nMC                 = 10,
                           Id seedMC              = 134341,
                           const CGParam& cgparams = CGParam());

  Id getBorder() const { return _border; }
  CGParam getCGparams() const { return _CGparams; }
  double getEpsNugget() const { return _epsNugget; }
  Id getRefineK() const { return _refineK; }
  Id getRefineS() const { return _refineS; }
  bool isPolarized() const { return _flagPolarized; }
  void setPolarized(bool flagPolarized) { _flagPolarized = flagPolarized; }
  Id getNxMax() const { return _nxmax; }
  bool getUseStencil() const { return _useStencil; }
  Id getNMC() const { return _nMC; }
  Id getSeedMC() const { return _seedMC; }

  void setBorder(Id border) { _border = border; }
  void setCGparams(const CGParam& CGparams) { _CGparams = CGparams; }
  void setEpsNugget(double epsNugget) { _epsNugget = epsNugget; }
  void setRefineK(Id refineK) { _refineK = refineK; }
  void setRefineS(Id refineS) { _refineS = refineS; }
  void setNxMax(Id nxmax) { _nxmax = nxmax; }
  void setUseStencil(bool useStencil) { _useStencil = useStencil; }
  void setNMC(Id nMC) { _nMC = nMC; }
  void setSeedMC(Id seedMC) { _seedMC = seedMC; }

private:
  Id _refineK;
  Id _refineS;
  Id _border;
  bool _flagPolarized;
  Id _nxmax;
  double _epsNugget;
  bool _useStencil;
  Id _nMC;    // Number of Monte-Carlo simulations
  Id _seedMC; // Seed for the random number generator
  CGParam _CGparams;
};
} // namespace gstlrn