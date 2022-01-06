/******************************************************************************/
/* COPYRIGHT ARMINES, ALL RIGHTS RESERVED                                     */
/*                                                                            */
/* THE CONTENT OF THIS WORK CONTAINS CONFIDENTIAL AND PROPRIETARY             */
/* INFORMATION OF ARMINES. ANY DUPLICATION, MODIFICATION,                     */
/* DISTRIBUTION, OR DISCLOSURE IN ANY FORM, IN WHOLE, OR IN PART, IS STRICTLY */
/* PROHIBITED WITHOUT THE PRIOR EXPRESS WRITTEN PERMISSION OF ARMINES         */
/*                                                                            */
/* TAG_SOURCE_CG                                                              */
/******************************************************************************/
#pragma once

#include "gstlearn_export.hpp"

// WARNING: Make this include list as small as possible!
#include "Basic/AStringable.hpp"

class GSTLEARN_EXPORT Option_VarioFit : public AStringable
{
private:
  int _flag_noreduce; /* Forbid discarding useless basic structures */
  int _flag_check_bounds; /* Do not infer parameter when bounds are equal */
  int _flag_goulard_used; /* 1 if Goulard must be used (for sills) */
  /* This is switch OFF when ANAM properties are defined */
  int _auth_aniso; /* Authorize the anisotropy */
  int _auth_rotation; /* Authorize the rotation of the anisotropy */
  int _lock_samerot; /* Lock the anisotropy rotation for all str */
  int _lock_rot2d; /* Lock the anisotropy rotation around Z only */
  int _lock_no3d; /* Lock the parameters in 2-D */
  int _lock_iso2d; /* Lock isotropy for 2-D */
  int _keep_intstr; /* Keep at least one intrinsic structure */

 public:
  Option_VarioFit();
  Option_VarioFit(const Option_VarioFit &m);
  Option_VarioFit& operator= (const Option_VarioFit &m);
  virtual ~Option_VarioFit();

  virtual String toString(const AStringFormat* strfmt = nullptr) const override;
  int getAuthAniso() const { return _auth_aniso; }
  void setAuthAniso(int authAniso) { _auth_aniso = authAniso; }
  int getAuthRotation() const { return _auth_rotation; }
  void setAuthRotation(int authRotation) { _auth_rotation = authRotation; }
  int getFlagCheckBounds() const { return _flag_check_bounds; }
  void setFlagCheckBounds(int flagCheckBounds) { _flag_check_bounds = flagCheckBounds; }
  int getFlagGoulardUsed() const { return _flag_goulard_used; }
  void setFlagGoulardUsed(int flagGoulardUsed) { _flag_goulard_used = flagGoulardUsed; }
  int getFlagNoreduce() const { return _flag_noreduce; }
  void setFlagNoreduce(int flagNoreduce) { _flag_noreduce = flagNoreduce; }
  int getKeepIntstr() const { return _keep_intstr; }
  void setKeepIntstr(int keepIntstr) { _keep_intstr = keepIntstr; }
  int getLockIso2d() const { return _lock_iso2d; }
  void setLockIso2d(int lockIso2d) { _lock_iso2d = lockIso2d; }
  int getLockNo3d() const { return _lock_no3d; }
  void setLockNo3d(int lockNo3d) { _lock_no3d = lockNo3d; }
  int getLockRot2d() const { return _lock_rot2d; }
  void setLockRot2d(int lockRot2d) { _lock_rot2d = lockRot2d; }
  int getLockSamerot() const { return _lock_samerot; }
  void setLockSamerot(int lockSamerot) { _lock_samerot = lockSamerot; }
};
