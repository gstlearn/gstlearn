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
 * - flag_noreduce: in the different iterations, some structures can be discarded
 *                  if their importance is considered as too small.
 *                  The current option forbids this simplification, which ensures
 *                  that all the basic structures are kept and their number remains unchanged.
 * - flag_goulard_used: This very efficient algorithm can be used for inferring
 *                  the sill matrix for each basic structure. This option can be
 *                  switched OFF on purpose (and replaced by the FOXLEG algorithm).
 * - auth_aniso: When True, the inference looks for an anisotropic fit.
 * - auth_rotation: When True, the inference looks for a possible rotation
 * - lock_samerot: When True, the inference locks the same anisotropy for all basic structures
 * - lock_rot2d: When True, the anisotropy is restricted to a rotation around Z-axis only
 * - lock_no3D: When True, the inference parameters are limited to the 2-D space.
 * - lock_iso2d: When True, the inference looks for a 2-D isotropic model
 * - keep_instr: When True, at least ONE basic structure must be kept in the Model
 * - flag_instrinsic: When True, fit a Model which includes at least one Intrinsic basic Structure
 */

class GSTLEARN_EXPORT Option_VarioFit: public AStringable
{
 public:
  Option_VarioFit(bool flag_noreduce = false,
                  bool auth_aniso = true,
                  bool auth_rotation = true,
                  bool lock_samerot = false,
                  bool lock_rot2d = false,
                  bool lock_no3d = false,
                  bool lock_iso2d = false);
  Option_VarioFit(const Option_VarioFit &m);
  Option_VarioFit& operator= (const Option_VarioFit &m);
  virtual ~Option_VarioFit();

  virtual String toString(const AStringFormat* strfmt = nullptr) const override;
  bool getAuthAniso() const { return _auth_aniso; }
  void setAuthAniso(bool authAniso) { _auth_aniso = authAniso; }
  bool getAuthRotation() const { return _auth_rotation; }
  void setAuthRotation(bool authRotation) { _auth_rotation = authRotation; }
  bool getFlagGoulardUsed() const { return _flag_goulard_used; }
  void setFlagGoulardUsed(bool flagGoulardUsed) { _flag_goulard_used = flagGoulardUsed; }
  bool getFlagNoreduce() const { return _flag_noreduce; }
  void setFlagNoreduce(bool flagNoreduce) { _flag_noreduce = flagNoreduce; }
  bool getKeepIntstr() const { return _keep_intstr; }
  void setKeepIntstr(bool keepIntstr) { _keep_intstr = keepIntstr; }
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

private:
  bool _flag_noreduce;     /* Forbid discarding useless basic structures */
  bool _flag_goulard_used; /* True if Goulard must be used (for sills) */
  bool _auth_aniso;        /* Authorize the anisotropy */
  bool _auth_rotation;     /* Authorize the rotation of the anisotropy */
  bool _lock_samerot;      /* Lock the anisotropy rotation for all str */
  bool _lock_rot2d;        /* Lock the anisotropy rotation around Z only */
  bool _lock_no3d;         /* Lock the parameters in 2-D */
  bool _lock_iso2d;        /* Lock isotropy for 2-D */
  bool _keep_intstr;       /* Keep at least one intrinsic structure */
  bool _flag_intrinsic;    /* Ask for an intrinsic model */
};
