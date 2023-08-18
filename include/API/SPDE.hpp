/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://gstlearn.org                                              */
/* License: BSD 3 clauses                                                     */
/*                                                                            */
/******************************************************************************/
#pragma once

#include "gstlearn_export.hpp"
#include "geoslib_define.h"

#include "Enum/ESPDECalcMode.hpp"

#include "Basic/NamingConvention.hpp"
#include "LinearOp/PrecisionOpCs.hpp"
#include "LinearOp/PrecisionOpMultiConditional.hpp"

class ShiftOpCs;
class Db;
class DbGrid;
class PrecisionOp;
class PrecisionOpCs;
class Model;
class MeshETurbo;

class GSTLEARN_EXPORT SPDE
{
public:
  SPDE();
  SPDE(Model *model,
       const DbGrid* field,
       const Db* data = nullptr,
       const ESPDECalcMode &calc = ESPDECalcMode::fromKey("SIMUCOND"));
  SPDE(const SPDE& r) = delete;
  SPDE& operator=(const SPDE& r) = delete;
  virtual ~SPDE();

  static SPDE* create(Model *model,
                      const DbGrid *field,
                      const Db *data = nullptr,
                      const ESPDECalcMode &calc = ESPDECalcMode::fromKey(
                          "SIMUCOND"));

  void init(Model* model,
            const DbGrid* field,
            const Db* data = nullptr,
            const ESPDECalcMode &calc = ESPDECalcMode::fromKey("SIMUCOND"),
            const AMesh* mesh = nullptr,
            bool verbose = false);
  void compute(int nbsimus = 1, int seed = 131323);
  void computeLk() const;
  void computeKriging() const;
  void computeSimuNonCond(int nbsimus = 1, int seed = 131323) const;
  void computeSimuCond(int nbsimus = 1, int seed = 131323) const;
  double computeLogLike(int nbsimus = 1, int seed = 131323) const;
  double computeProfiledLogLike(int nbsimus = 1, int seed = 131323) const;
  VectorDouble getCoeffs();
  void centerByDrift(const VectorDouble& dataVect,int ivar=0,bool useSel=true) const;

  void setDriftCoeffs(VectorDouble coeffs);
  double computeLogDet(int nbsimus = 1,int seed = 1234) const;
  int query(Db *db,
            const NamingConvention &namconv = NamingConvention("spde")) const;
  const PrecisionOpCs* getPrecisionOp(int i = 0) const  { return (PrecisionOpCs*)_pilePrecisions[i];}
  const ProjMatrix* getProj(int i = 0) const  { return _pileProjMatrix[i];}
  const PrecisionOpMultiConditional* getPrecisionKriging() const { return _precisionsKriging;}

  double computeQuad() const;
  const Db* getData() const {return  _data;}

  void setEps(double eps) { _eps = eps; }
  void setNIterMax(int nitermax) { _nIterMax = nitermax; }

private:
  void _computeDriftCoeffs() const;
  void _purge();
  MeshETurbo* _createMeshing(const CovAniso &cova,
                             const DbGrid &field,
                             double discr,
                             double ext = 0.);
  bool _calculSimu() const;
  bool _calculKriging() const;

private:
  const Db*_data;
  ESPDECalcMode _calcul;
  PrecisionOpMultiConditional* _precisionsKriging;
  PrecisionOpMultiConditional* _precisionsSimu;
  std::vector<PrecisionOp*>    _pilePrecisions;
  std::vector<ProjMatrix*>     _pileProjMatrix;
  std::vector<const AMesh*>    _simuMeshing;
  std::vector<const AMesh*>    _krigingMeshing;
  mutable VectorDouble         _driftCoeffs;
  Model*                       _model;
  mutable VectorVectorDouble   _workKriging;
  mutable VectorVectorDouble   _workingSimu;
  mutable VectorDouble         _workingData;
  mutable VectorDouble         _workingData2;
  std::vector<ProjMatrix*>     _projOnDbOut;
  VectorInt                    _adressesICov;
  double _nugget;
  VectorVectorDouble _driftTab;
  bool _requireCoeffs;
  mutable bool _isCoeffsComputed;
  bool _deleteMesh;
  // query sur aproj ou // TODO ??

  // Parameters specific invertion using Conjugate Gradient (used for Kriging)
  int _nIterMax;
  double _eps;
};
