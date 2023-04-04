/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://github.com/gstlearn                                       */
/* License: BSD 3 clause                                                      */
/*                                                                            */
/******************************************************************************/
#pragma once

#include "gstlearn_export.hpp"

#include "Basic/VectorNumT.hpp"
#include "Model/ElemNostat.hpp"
#include "Basic/AStringable.hpp"

class ElemNostat;
class CovAniso;

class GSTLEARN_EXPORT ModelNostat : public AStringable
{
public:
  ModelNostat();
  ModelNostat(const ModelNostat &m);
  ModelNostat& operator= (const ModelNostat &m);
  virtual ~ModelNostat();

  void init(int ndim);
  void init(const CovAniso* cova);
  void define(int icov, const CovAniso* cova);
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;
  ElemNostat* addElemNostat();

  int getFlagNoStat() const { return getNElems() > 0; }

  const VectorDouble& getAngles1() const { return _angles1; }
  double getAngles1(int i) const { return _angles1[i]; }
  const VectorDouble& getAngles2() const { return _angles2; }
  double getAngles2(int i) const { return _angles2[i]; }
  int getNDim() const { return _nDim; }
  int getNElems() const { return static_cast<int> (_elems.size()); }
  double getParam1() const { return _param1; }
  double getParam2() const { return _param2; }
  double getScadef1() const { return _scadef1; }
  double getScadef2() const { return _scadef2; }
  const VectorDouble& getScale1() const { return _scale1; }
  double getScale1(int i) const { return _scale1[i]; }
  const VectorDouble& getScale2() const { return _scale2; }
  double getScale2(int i) const { return _scale2[i]; }
  double getSill1() const { return _sill1; }
  double getSill2() const { return _sill2; }
  const std::vector<ElemNostat*>& getElems() const { return _elems; }
  ElemNostat* getElems(int i) const { return _elems[i]; }
  void setAngles1(int idim, double value) { _setAngles1(idim, value); }
  void setAngles2(int idim, double value) { _setAngles2(idim, value); }
  void setScale1(VectorDouble scales) { _setScale1(scales); }
  void setScale1(int idim, double scale) { _setScale1(idim,scale); }
  void setScale2(VectorDouble scales) { _setScale2(scales); }
  void setScale2(int idim, double scale) { _setScale2(idim,scale); }

private:
  void _setAngles1(int idim, double value) { _angles1[idim] = value; }
  void _setAngles2(int idim, double value) { _angles2[idim] = value; }
  void _setScale1(VectorDouble scales) { _scale1 = scales; }
  void _setScale1(int idim, double scale) { _scale1[idim] = scale; }
  void _setScale2(VectorDouble scales) { _scale2 = scales; }
  void _setScale2(int idim, double scale) { _scale2[idim] = scale; }

private:
  int _nDim;
  std::vector<ElemNostat*> _elems; /* List of non-stationary parameters */
  double _sill1; /* First Sill */
  double _sill2; /* Second Sill */
  double _param1; /* First version of the third parameter */
  double _param2; /* Second version of the third parameter */
  double _scadef1; /* First normalization scale */
  double _scadef2; /* Second normalization scale */
  VectorDouble _angles1; /* Array for first angles (Dimension: _nDim) */
  VectorDouble _angles2; /* Array for second angles (Dimension: _nDim) */
  VectorDouble _scale1; /* Array of first theoretical ranges (Dim: _nDim) */
  VectorDouble _scale2; /* Array of second theoretical ranges (Dim: _nDim) */
};
