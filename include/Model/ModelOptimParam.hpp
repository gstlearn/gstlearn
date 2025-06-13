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

// WARNING: Make this include list as small as possible!
#include "Basic/AStringable.hpp"

/**
 * This class defines the options and parameters used during the Variogram Fitting.
 * All the parameters described hereafter are either available in the construction,
 * or can be set using a specific get() function.
 * - auth_aniso: When True, the inference looks for an anisotropic fit.
 * - auth_rotation: When True, the inference looks for a possible rotation
 * - lock_samerot: When True, the inference locks the same anisotropy for all basic structures
 * - lock_rot2d: When True, the anisotropy is restricted to a rotation around Z-axis only
 * - lock_no3D: When True, the inference parameters are limited to the 2-D space.
 * - lock_iso2d: When True, the inference looks for a 2-D isotropic model
 * - flag_goulard:  This very efficient algorithm can be used for inferring
 *                  the sill matrix for each basic structure. This option can be
 *                  switched OFF on purpose (and replaced by the FOXLEG algorithm).
 * - flag_instrinsic: When True, fit a Model which includes at least one Intrinsic basic Structure
 * - wmode: Weighting option (see comments on the setWMode() function)
 * - maxiter: Maximum number of iterations
 * - tolred: Define the relative criterion used for stopping the iterations
 */

class GSTLEARN_EXPORT ModelOptimParam: public AStringable
{
 public:
  ModelOptimParam();
  ModelOptimParam(const ModelOptimParam &m);
  ModelOptimParam& operator= (const ModelOptimParam &m);
  virtual ~ModelOptimParam();

  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  static ModelOptimParam* create(bool auth_aniso     = true,
                                 bool auth_rotation  = true,
                                 bool lock_samerot   = false,
                                 bool lock_rot2d     = false,
                                 bool lock_no3d      = false,
                                 bool lock_iso2d     = false,
                                 bool flag_goulard   = true,
                                 bool flag_intrinsic = false,
                                 int wmode           = 2,
                                 int maxiter         = 1000,
                                 double tolred       = EPSILON6);

  bool getAuthAniso() const { return _auth_aniso; }
  void setAuthAniso(bool authAniso) { _auth_aniso = authAniso; }
  bool getAuthRotation() const { return _auth_rotation; }
  void setAuthRotation(bool authRotation) { _auth_rotation = authRotation; }
  bool getFlagGoulard() const { return _flag_goulard; }
  void setFlagGoulard(bool flagGoulard) { _flag_goulard = flagGoulard; }
  bool getLockIso2d() const { return _lock_iso2d; }
  void setLockIso2d(bool lockIso2d) { _lock_iso2d = lockIso2d; }
  bool getLockNo3d() const { return _lock_no3d; }
  void setLockNo3d(bool lockNo3d) { _lock_no3d = lockNo3d; }
  bool getLockRot2d() const { return _lock_rot2d; }
  void setLockRot2d(bool lockRot2d) { _lock_rot2d = lockRot2d; }
  bool getLockSamerot() const { return _lock_samerot; }
  void setLockSamerot(bool lockSamerot) { _lock_samerot = lockSamerot; }
  bool getFlagIntrinsic() const { return _flag_intrinsic; }
  void setFlagIntrinsic(bool flagIntrinsic) { _flag_intrinsic = flagIntrinsic; }
  /**
   * Set the type of the weighting function used in the fitting procedure.
   * This function is defined in the case of several directional experimental variograms,
   * calculated in a multivariate case:
   * 0: The weight is constant
   * 1: The weight is proportional to the number of pairs
   * 2: The weight is proportional to the number of pairs and inverse proportional to the distance
   * 3: The weight is inverse proportional to the number of lags for each direction
   * @param wmode       type of weighting function (0, 1, 2 or 3, see above)
   * @note The default value for wmode is 2
   */
  int    getWmode() const { return _wmode; }
  void   setWmode(int wmode) { _wmode = wmode; }
  int    getMaxiter() const { return _maxiter; }
  void   setMaxiter(int maxiter) { _maxiter = maxiter; }
  double getTolred() const { return _tolred; }
  void   setTolred(double tolred) { _tolred = tolred; }

private:
  bool _auth_aniso;        /* Authorize the anisotropy */
  bool _auth_rotation;     /* Authorize the rotation of the anisotropy */
  bool _lock_samerot;      /* Lock the anisotropy rotation for all str */
  bool _lock_rot2d;        /* Lock the anisotropy rotation around Z only */
  bool _lock_no3d;         /* Lock the parameters in 2-D */
  bool _lock_iso2d;        /* Lock isotropy for 2-D */
  bool _flag_goulard;      /* True if Goulard must be used (for sills) */
  bool _flag_intrinsic;    /* Ask for an intrinsic model */
  int  _wmode;             /* Weighting option (used in Goulard) */
  int  _maxiter;           /* Maximum number of iterations */
  double _tolred;          /* Scaled tolerance (used in calculations) */
};
