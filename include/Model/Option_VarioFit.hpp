/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://gstlearn.org                                              */
/* License: BSD 3-clause                                                      */
/*                                                                            */
/******************************************************************************/
#pragma once

#include "gstlearn_export.hpp"

// WARNING: Make this include list as small as possible!
#include "Basic/AStringable.hpp"

class GSTLEARN_EXPORT Option_VarioFit : public AStringable
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
  bool getFlagCheckBounds() const { return _flag_check_bounds; }
  void setFlagCheckBounds(bool flagCheckBounds) { _flag_check_bounds = flagCheckBounds; }
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

 private:
   bool _flag_noreduce; /* Forbid discarding useless basic structures */
   bool _flag_check_bounds; /* Do not infer parameter when bounds are equal */
   bool _flag_goulard_used; /* 1 if Goulard must be used (for sills) */
   /* This is switch OFF when ANAM properties are defined */
   bool _auth_aniso; /* Authorize the anisotropy */
   bool _auth_rotation; /* Authorize the rotation of the anisotropy */
   bool _lock_samerot; /* Lock the anisotropy rotation for all str */
   bool _lock_rot2d; /* Lock the anisotropy rotation around Z only */
   bool _lock_no3d; /* Lock the parameters in 2-D */
   bool _lock_iso2d; /* Lock isotropy for 2-D */
   bool _keep_intstr; /* Keep at least one intrinsic structure */
};
