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
#include "Simulation/TurningBandOperate.hpp"

 TurningBandOperate:: TurningBandOperate()
    : _nt0(0),
      _flagScaled(false),
      _vexp(0.),
      _tdeb(0.),
      _omega(0.),
      _phi(0.),
      _offset(0.),
      _scale(1.),
      _t(),
      _v0(),
      _v1(),
      _v2()
{
}

 TurningBandOperate:: TurningBandOperate(const  TurningBandOperate &r)
    : _nt0(r._nt0),
      _flagScaled(r._flagScaled),
      _vexp(r._vexp),
      _tdeb(r._tdeb),
      _omega(r._omega),
      _phi(r._phi),
      _offset(r._offset),
      _scale(r._scale),
      _t(r._t),
      _v0(r._v0),
      _v1(r._v1),
      _v2(r._v2)
{
}

 TurningBandOperate&  TurningBandOperate::operator=(const  TurningBandOperate &r)
{
  if (this != &r)
  {
    _nt0 = r._nt0;
    _flagScaled = r._flagScaled;
    _vexp = r._vexp;
    _tdeb = r._tdeb;
    _omega = r._omega;
    _phi = r._phi;
    _offset = r._offset;
    _scale = r._scale;
    _t = r._t;
    _v0 = r._v0;
    _v1 = r._v1;
    _v2 = r._v2;
  }
  return *this;
}

 TurningBandOperate::~ TurningBandOperate()
{
}

void TurningBandOperate::reset()
{
  _nt0 = 0;
  _vexp = 0.;
  _tdeb = 0.;
  _omega = 0.;
  _phi = 0.;
  _offset = 0.;
  _scale = 1.;
  _t.clear();
  _v0.clear();
  _v1.clear();
  _v2.clear();
}

double TurningBandOperate::shotNoiseAffineOne(double t0)
{
  double scale = getScale();
  double tdeb = getTdeb() / scale;
  if (! isFlagScaled()) t0 /= scale;

  double dt = t0 - tdeb;
  int nt0 = (int) (dt);
  double dt0 = dt - nt0;
  return _t[nt0] * (2. * dt0 - 1.);
}

double TurningBandOperate::shotNoiseCubicOne(double t0)
{
  double scale = getScale();
  double tdeb = getTdeb() / scale;
  if (! isFlagScaled()) t0 /= scale;

  double dt = t0 - tdeb;
  int nt0 = (int) (dt);
  double dt0 = dt - nt0;
  return _t[nt0] * dt0 * (dt0 - 0.5) * (dt0 - 1.);
}

double TurningBandOperate::spectralOne(double t0)
{
  int nt0 = _rankInPoisson(getNt0(), t0, _t);
  setNt0(nt0);

  double vexp = getVexp();
  return (2. * t0 > _t[nt0 + 1] + _t[nt0]) ? -vexp : vexp;
}

double TurningBandOperate::IRFProcessOne(double t0)
{
  int nt0 = _rankInPoisson(getNt0(), t0, _t);
  setNt0(nt0);

  return _irfProcessSample(nt0, t0);
}

double TurningBandOperate::cosineOne(double t0) const
{
  double offset = getOffset();
  if (isFlagScaled()) return t0 - offset;
  double omega = getOmega();
  double phi   = getPhi();
  return cos(omega * t0 + phi) - offset;
}

/*****************************************************************************/
/*!
 **  Sample the Wiener-Levy (integrated) process
 **
 ** \param[in]  nt0    Rank of the Poisson point
 ** \param[in]  t0     starting time
 **
 *****************************************************************************/
double TurningBandOperate::_irfProcessSample(int nt0, double t0)
{
  double value;

  /* Initializations */

  double delta = t0 - _t[nt0];
  if (_v0.empty()) return TEST;

  /* Wiener-Levy process */

  value = _v0[nt0];
  if (_v1.empty()) return value;

  /* First integration of the Wiener-Levy process */

  value = _v1[nt0] + _v0[nt0] * delta;
  if (_v2.empty()) return value;

  /* Second integration of the Wiener-Levy process */

  value = _v2[nt0] + _v1[nt0] * delta + _v0[nt0] * delta * delta / 2.;
  return value;
}

/*****************************************************************************/
/*!
 **  Returns the rank of the point t0 in the Poisson point process
 **
 ** \param[in]  def_rank Rank of the Poisson point
 ** \param[in]  t0 starting time
 ** \param[in]  t  Poisson point process
 **
 *****************************************************************************/
int TurningBandOperate::_rankInPoisson(int def_rank,
                                       double t0,
                                       const VectorDouble &t)
{
  int it, itp, itn;

  /* First, try with the default interval then the next one and finally
   the previous one */

  int nt = (int) t.size();
  if (t0 >= t[def_rank] && t0 < t[def_rank + 1])
    return (def_rank);
  if (def_rank < (nt - 2) && t0 >= t[def_rank + 1] && t0 < t[def_rank + 2])
    return (def_rank + 1);
  if (def_rank > 0 && t0 >= t[def_rank - 1] && t0 < t[def_rank])
    return (def_rank - 1);

  /* The default value is not good ==> dichotomy */

  itp = 0;
  itn = nt - 1;
  while (itn - itp > 1)
  {
    it = (itn + itp) / 2;
    if (t0 >= t[it])
      itp = it;
    else
      itn = it;
  }
  return (itp);
}

void TurningBandOperate::pushT(double value)
{
  _t.push_back(value);
}

void TurningBandOperate::pushV0(double value)
{
  _v0.push_back(value);
}

void TurningBandOperate::pushV1(double value)
{
  _v1.push_back(value);
}

void TurningBandOperate::pushV2(double value)
{
  _v2.push_back(value);
}
