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

#include "ACalcSimulation.hpp"
#include "Basic/VectorNumT.hpp"

class Db;

/**
 * Class for management of Directions used in Turning Band algorithm
 * Remark: The 3-D definition is compulsory (even in 2-D)
 */
class GSTLEARN_EXPORT TurningBandOperate
{
public:
  TurningBandOperate();
  TurningBandOperate(const TurningBandOperate& r);
  TurningBandOperate& operator=(const TurningBandOperate& r);
  virtual ~TurningBandOperate();

  const VectorDouble& getT()  const { return _t; }
  const VectorDouble& getV0() const { return _v0; }
  const VectorDouble& getV1() const { return _v1; }
  const VectorDouble& getV2() const { return _v2; }
  int getNt0() const { return _nt0; }
  bool isFlagScaled() const { return _flagScaled; }
  double getVexp() const { return _vexp; }
  double getTdeb() const { return _tdeb; }
  double getOmega() const { return _omega; }
  double getPhi() const { return _phi; }
  double getOffset() const { return _offset; }
  double getScale() const { return _scale; }

  void setT(const VectorDouble &t)   { _t  = t; }
  void setV0(const VectorDouble &v0) { _v0 = v0; }
  void setV1(const VectorDouble &v1) { _v1 = v1; }
  void setV2(const VectorDouble &v2) { _v2 = v2; }
  void setNt0(int nt0) { _nt0 = nt0; }
  void setFlagScaled(bool flagScaled) { _flagScaled = flagScaled; }
  void setVexp(double vexp) { _vexp = vexp; }
  void setTdeb(double tdeb) { _tdeb = tdeb; }
  void setOmega(double omega) { _omega = omega; }
  void setPhi(double phi) { _phi = phi; }
  void setOffset(double offset) { _offset = offset; }
  void setScale(double scale) { _scale = scale; }

  int getTsize() const { return _t.size(); }
  void pushT(double value);
  void pushV0(double value);
  void pushV1(double value);
  void pushV2(double value);

  void reset();
  double shotNoiseAffineOne(double t0);
  double shotNoiseCubicOne(double t0);
  double spectralOne(double t0);
  double IRFProcessOne(double t0);
  double cosineOne(double t0);

private:
  double _irfProcessSample(int nt0, double t0);
  int _rankInPoisson(int def_rank, double t0, const VectorDouble &t);

private:
  int _nt0;
  bool _flagScaled;
  double _vexp;
  double _tdeb;
  double _omega;
  double _phi;
  double _offset;
  double _scale; // Scale factor by which 't0' should be scaled here
  VectorDouble _t;
  VectorDouble _v0;
  VectorDouble _v1;
  VectorDouble _v2;
};
