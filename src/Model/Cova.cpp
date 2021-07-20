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
#include "Model/Cova.hpp"
#include "geoslib_f.h"
#include "geoslib_f_private.h"

Cova::Cova()
    : _type(0),
      _nVar(0),
      _nDim(0),
      _flagFilter(0),
      _flagAniso(0),
      _flagRotation(0),
      _covChoice(0),
      _range(0),
      _param(0),
      _sill(),
      _anisoAngles(),
      _anisoCoeffs(),
      _anisoRotMat()
{
  _st_cov_external = nullptr;
}

Cova::Cova(const Cova &m)
: _type(m._type),
  _nVar(m._nVar),
  _nDim(m._nDim),
  _flagFilter(m._flagFilter),
  _flagAniso(m._flagAniso),
  _flagRotation(m._flagRotation),
  _covChoice(m._covChoice),
  _range(m._range),
  _param(m._param),
  _sill(m._sill),
  _anisoAngles(m._anisoAngles),
  _anisoCoeffs(m._anisoCoeffs),
  _anisoRotMat(m._anisoRotMat)
{
  _st_cov_external = m._st_cov_external;
}
Cova& Cova::operator=(const Cova &m)
{
  if (this != &m)
  {
    _type = m._type;
    _nVar = m._nVar;
    _nDim = m._nDim;
    _flagFilter = m._flagFilter;
    _flagAniso = m._flagAniso;
    _flagRotation = m._flagRotation;
    _covChoice = m._covChoice;
    _range = m._range;
    _param = m._param;
    _sill = m._sill;
    _anisoAngles = m._anisoAngles;
    _anisoCoeffs = m._anisoCoeffs;
    _anisoRotMat = m._anisoRotMat;
    _st_cov_external = m._st_cov_external;
  }
  return *this;
}
Cova::~Cova()
{

}

void Cova::init(int ndim,
                int nvar)
{
  _nVar = nvar;
  _nDim = ndim;
  _type = 0;
  _flagAniso = 0;
  _flagRotation = 0;
  _range = 0;
  _param = 0;
  _sill.resize(_nVar * _nVar);
  _anisoRotMat.resize(_nDim * _nDim);
  _flagFilter = 0;
  _covChoice = MODEL_CALCUL_NATURAL;
  _anisoAngles.resize(_nDim,0.);
  _anisoCoeffs.resize(_nDim,1.);
}

void Cova::init(int ndim,
                int nvar,
                int type,
                double range,
                double param,
                int flag_aniso,
                int flag_rotation,
                const VectorDouble& aniso_ranges,
                const VectorDouble& aniso_rotmat,
                const VectorDouble& sill)
{
  _nVar = nvar;
  _nDim = ndim;
  _type = type;
  _flagAniso = flag_aniso;
  _flagRotation = flag_rotation;
  _range = range;
  _param = param;

  // Sill Matrix
  if (sill.empty())
  {
    _sill.resize(_nVar * _nVar);
    int ecr = 0;
    for (int ivar=0; ivar<_nVar; ivar++)
      for (int jvar=0; jvar<_nVar; jvar++)
        _sill[ecr++] = ivar == jvar;
  }
  else
    _sill = sill;

  // Anisotropy Rotation Matrix
  if (aniso_rotmat.empty())
  {
    _anisoRotMat.resize(_nDim * _nDim);
    int ecr = 0;
    for (int idim=0; idim<_nDim; idim++)
      for (int jdim=0; jdim<_nDim; jdim++)
        _anisoRotMat[ecr++] = idim == jdim;
  }
  else
    _anisoRotMat = aniso_rotmat;

  _flagFilter = 0;
  _covChoice = MODEL_CALCUL_NATURAL;

  // Anisotropy Rotation Angles
  // fromRotMatToAngles();

  // Anisotropy coefficients
  // fromRangesToCoeffs(aniso_ranges);
}

double _st_cov_exp2dfact(Cova *covariance,
                         int ndim,
                         double field,
                         VectorDouble d,
                         double h)
{
  double h1, h2, cov;
  int i;

  cov = 1.;
  h1 = 0.;
  if (! d.empty()) for (i = 0; i < MIN(ndim, 2); i++)
    h1 += d[i] * d[i];
  if (h1 > MAX_EXP) return (0.);
  cov = exp(-h1);
  for (i = 2; i < ndim; i++)
  {
    h2 = (! d.empty()) ? ABS(d[i]) :
                                  0.;
    if (h2 > MAX_EXP2) return (0.);
    cov *= exp(-h2);
  }

  return (cov);
}

double _st_cov_expfact(Cova *covariance,
                       int ndim,
                       double field,
                       VectorDouble d,
                       double h)
{
  double cov;
  int i;

  cov = 1.;
  for (i = 0; i < ndim; i++)
  {
    h = 0.;
    if (! d.empty()) h = ABS(d[i]);
    if (h > MAX_EXP) return (0.);
    cov *= exp(-h);
  }

  return (cov);
}

double _st_cov_func(Cova *covariance,
                    int ndim,
                    double field,
                    VectorDouble d,
                    double h)
{
  External_Cov E_Cov;

  fill_external_cov_kriging(E_Cov);
  fill_external_cov_model(E_Cov);

  return covariance->evaluateExternalCov(E_Cov, h, d);
}

//static Def_Cova DEF_FUNC =
//{"External Cov Function", 0,0,-1,-1, 0, 0, 1.      ,  0., _st_cov_func,
//  NULL, NULL, NULL, NULL};

//  {"Exp2dfact"    , 1,0,-1,-1, 1, 1, 2.995732,  0.,_st_cov_exp2dfact,
//  NULL, NULL, NULL, NULL},
//  {"Expfact"      , 1,0,-1,-1, 1, 1, 2.995732,  0.,_st_cov_expfact,
//   NULL, NULL, NULL, NULL},


