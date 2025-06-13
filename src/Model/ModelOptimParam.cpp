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
#include "Model/ModelOptimParam.hpp"
#include "Basic/AStringable.hpp"

ModelOptimParam::ModelOptimParam()
  : AStringable()
  , _auth_aniso(true)
  , _auth_rotation(true)
  , _lock_samerot(false)
  , _lock_rot2d(false)
  , _lock_no3d(false)
  , _lock_iso2d(false)
  , _flag_goulard(true)
  , _flag_intrinsic(false)
  , _wmode(2)
  , _maxiter(1000)
  , _tolred(EPSILON6)
{
}

ModelOptimParam::ModelOptimParam(const ModelOptimParam& m)
  : AStringable(m)
  , _auth_aniso(m._auth_aniso)
  , _auth_rotation(m._auth_rotation)
  , _lock_samerot(m._lock_samerot)
  , _lock_rot2d(m._lock_rot2d)
  , _lock_no3d(m._lock_no3d)
  , _lock_iso2d(m._lock_iso2d)
  , _flag_goulard(m._flag_goulard)
  , _flag_intrinsic(m._flag_intrinsic)
  , _wmode(m._wmode)
  , _maxiter(m._maxiter)
  , _tolred(m._tolred)
{
}

ModelOptimParam& ModelOptimParam::operator=(const ModelOptimParam &m)
{
  if (this != &m)
  {
    AStringable::operator=(m);

    _auth_aniso     = m._auth_aniso;
    _auth_rotation  = m._auth_rotation;
    _lock_samerot   = m._lock_samerot;
    _lock_rot2d     = m._lock_rot2d;
    _lock_no3d      = m._lock_no3d;
    _lock_iso2d     = m._lock_iso2d;
    _flag_goulard   = m._flag_goulard;
    _flag_intrinsic = m._flag_intrinsic;
    _wmode          = m._wmode;
    _maxiter        = m._maxiter;
    _tolred         = m._tolred;
  }
  return *this;
}

ModelOptimParam::~ModelOptimParam()
{

}

ModelOptimParam* ModelOptimParam::create(bool auth_aniso,
                                         bool auth_rotation,
                                         bool lock_samerot,
                                         bool lock_rot2d,
                                         bool lock_no3d,
                                         bool lock_iso2d,
                                         bool flag_goulard,
                                         bool flag_intrinsic,
                                         int wmode,
                                         int maxiter,
                                         double tolred)
{
  ModelOptimParam* mop = new ModelOptimParam();
  mop->_auth_aniso     = auth_aniso;
  mop->_auth_rotation  = auth_rotation;
  mop->_lock_samerot   = lock_samerot;
  mop->_lock_rot2d     = lock_rot2d;
  mop->_lock_no3d      = lock_no3d;
  mop->_lock_iso2d     = lock_iso2d;
  mop->_flag_goulard   = flag_goulard;
  mop->_flag_intrinsic = flag_intrinsic;
  mop->_wmode          = wmode;
  mop->_maxiter        = maxiter;
  mop->_tolred         = tolred;
  return mop;
}

String ModelOptimParam::toString(const AStringFormat* /*strfmt*/) const
{
  std::stringstream sstr;
  static const char *NOK[] = {"OFF" , "ON"};

  sstr << "- Anisotropy                   " << NOK[getAuthAniso()]       << std::endl;
  sstr << "- Anisotropy Rotation          " << NOK[getAuthRotation()]    << std::endl;
  sstr << "- Global Rotation              " << NOK[getLockSamerot()]     << std::endl;
  sstr << "- Rotation around Z only       " << NOK[getLockRot2d()]       << std::endl;
  sstr << "- Lock third dimension         " << NOK[getLockNo3d()]        << std::endl;
  sstr << "- Lock 2-D Isotropy            " << NOK[getLockIso2d()]       << std::endl;
  sstr << "- Use the Goulard option       " << NOK[getFlagGoulard()]     << std::endl;
  if (getFlagIntrinsic())
    sstr << "- Multivariate Model must be Intrinsic" << std::endl;
  sstr << "- Optimization weighting mode  " << getWmode() << std::endl;
  sstr << "- Maximum number of iterations " << getMaxiter() << std::endl;
  sstr << "- Stopping Criterion (scaled)  " << getTolred() << std::endl;

  return sstr.str();
}
