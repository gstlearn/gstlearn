#pragma once

#include "gstlearn_export.hpp"
#include "geoslib_define.h"

// Enums
#include "API/ESPDECalcMode.hpp"

#include "Basic/NamingConvention.hpp"
#include "LinearOp/PrecisionOpMultiConditional.hpp"

class ShiftOpCs;
class Db;
class PrecisionOpCs;
class Model;
class MeshETurbo;

class GSTLEARN_EXPORT SPDE
{
public:
  SPDE();
  SPDE(Model *model,
       const Db* field,
       const Db* dat = nullptr,
       const ESPDECalcMode &calc = ESPDECalcMode::SIMUCOND);
  SPDE(const SPDE& r) = delete;
  SPDE& operator=(const SPDE& r) = delete;
  virtual ~SPDE();

  void init(Model* model,
            const Db* field,
            const Db* dat = nullptr,
            const ESPDECalcMode &calc = ESPDECalcMode::SIMUCOND);
  void compute(int nbsimus = 1, int seed = 131323) const; // TODO What this seed ?
  void computeKriging() const;
  void computeSimuNonCond(int nbsimus = 1, int seed = 131323) const;
  void computeSimuCond(int nbsimus = 1, int seed = 131323) const;
  VectorDouble computeCoeffs() const;
  int query(Db *db,
            const NamingConvention &namconv = NamingConvention("spde")) const;

private:
  void _purge();
  MeshETurbo* _createMeshing(const CovAniso &cova,
                             const Db &field,
                             double discr,
                             double ext = 0.);
  bool _calculSimu() const
  {
    return _calcul == ESPDECalcMode::SIMUCOND
        || _calcul == ESPDECalcMode::SIMUNONCOND;
  }
  bool _calculKriging() const
  {
    return ((_calcul == ESPDECalcMode::SIMUCOND
          || _calcul == ESPDECalcMode::KRIGING)
          && _data != nullptr);
  }

private:
  const Db*_data;
  ESPDECalcMode _calcul;
  PrecisionOpMultiConditional _precisionsKriging;
  PrecisionOpMultiConditional _precisionsSimu;
  std::vector<ShiftOpCs*> _pileShiftOp;
  std::vector<PrecisionOpCs*> _pilePrecisions;
  std::vector<ProjMatrix*> _pileProjMatrix;
  std::vector<MeshETurbo*> _simuMeshing;
  std::vector<MeshETurbo*> _krigingMeshing;
  mutable VectorDouble     _driftCoeffs;
  Model *_model;
  mutable VectorVectorDouble _workKriging;
  mutable VectorVectorDouble _workingSimu;
  mutable VectorDouble       _workingData;
  std::vector<ProjMatrix*> _projOnDbOut;
  std::vector<int>           _adressesICov;
  double _nugget;
  VectorVectorDouble _driftTab;
  bool _computeCoeffs;
  // query sur aproj ou // TODO ??
};
