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
#include "Model/Option_VarioFit.hpp"
#include "Basic/AStringable.hpp"

Option_VarioFit::Option_VarioFit()
    : AStringable(),
      _flag_noreduce(false),
      _flag_check_bounds(false),
      _flag_goulard_used(true),
      _auth_aniso(true),
      _auth_rotation(true),
      _lock_samerot(false),
      _lock_rot2d(false),
      _lock_no3d(false),
      _lock_iso2d(false),
      _keep_intstr(false)
{
}

Option_VarioFit::Option_VarioFit(const Option_VarioFit &m)
    : AStringable(m),
      _flag_noreduce(m._flag_noreduce),
      _flag_check_bounds(m._flag_check_bounds),
      _flag_goulard_used(m._flag_goulard_used),
      _auth_aniso(m._auth_aniso),
      _auth_rotation(m._auth_rotation),
      _lock_samerot(m._lock_samerot),
      _lock_rot2d(m._lock_rot2d),
      _lock_no3d(m._lock_no3d),
      _lock_iso2d(m._lock_iso2d),
      _keep_intstr(m._keep_intstr)
{

}

Option_VarioFit& Option_VarioFit::operator=(const Option_VarioFit &m)
{
  if (this != &m)
  {
    AStringable::operator=(m);
    _flag_noreduce = m._flag_noreduce;
    _flag_check_bounds = m._flag_check_bounds;
    _flag_goulard_used = m._flag_goulard_used;
    _auth_aniso = m._auth_aniso;
    _auth_rotation = m._auth_rotation;
    _lock_samerot = m._lock_samerot;
    _lock_rot2d = m._lock_rot2d;
    _lock_no3d = m._lock_no3d;
    _lock_iso2d = m._lock_iso2d;
    _keep_intstr = m._keep_intstr;
  }
  return *this;
}

Option_VarioFit::~Option_VarioFit()
{

}

String Option_VarioFit::toString(const AStringFormat* /*strfmt*/) const
{
  std::stringstream sstr;
  static const char *NOK[] = {"OFF" , "ON"};

  sstr << "- Anisotropy                " << NOK[getAuthAniso()]       << std::endl;
  sstr << "- Anisotropy Rotation       " << NOK[getAuthRotation()]    << std::endl;
  sstr << "- Global Rotation           " << NOK[getLockSamerot()]     << std::endl;
  sstr << "- Rotation around Z only    " << NOK[getLockRot2d()]       << std::endl;
  sstr << "- Lock third dimension      " << NOK[getLockNo3d()]        << std::endl;
  sstr << "- Lock 2-D Isotropy         " << NOK[getLockIso2d()]       << std::endl;
  sstr << "- Keep Intrinsic structure  " << NOK[getKeepIntstr()]      << std::endl;
  sstr << "- Use the Goulard option    " << NOK[getFlagGoulardUsed()] << std::endl;
  sstr << "- Keep all structures       " << NOK[getFlagNoreduce()]    << std::endl;

  return sstr.str();
}
