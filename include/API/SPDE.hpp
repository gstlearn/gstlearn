#pragma once

#include "Basic/NamingConvention.hpp"
#include "Db/Db.hpp"
#include "LinearOp/PrecisionOpCs.hpp"
#include "LinearOp/PrecisionOpMultiConditional.hpp"
#include "Model/Model.hpp"
#include "Mesh/MeshETurbo.hpp"
#include "API/ESPDECalcMode.hpp"
#include <vector>

class ShiftOpCs;
class SPDE
{
public:
  SPDE();
  SPDE(Model& model,
       const Db& field,
       const Db* dat = nullptr,
       const ESPDECalcMode& calc = ESPDECalcMode::SIMUCOND);
  virtual ~SPDE();

  void init(Model& model,
            const Db& field,
            const Db* dat = nullptr,
            const ESPDECalcMode& calc = ESPDECalcMode::SIMUCOND);
  void compute(int nbsimus = 1, int seed = 131323) const;  // TODO What this seed ?
  void computeKriging(const VectorDouble& vect) const;
  void computeSimuNonCond(int nbsimus = 1, int seed = 131323) const;
  void computeSimuCond(int nbsimus = 1, int seed = 131323) const;
  VectorDouble computeCoeffs()const;
  int query(Db* db,NamingConvention namconv = NamingConvention("spde")) const;

private:
  void _purge();
  MeshETurbo* _createMeshing(const CovAniso& cova,
                             const Db& field,
                             double discr,
                             double ext = 0.);
  bool _calculSimu() const{return _calcul == ESPDECalcMode::SIMUCOND ||
                                  _calcul == ESPDECalcMode::SIMUNONCOND;}
  bool _calculKriging() const
  {
    return ((_calcul == ESPDECalcMode::SIMUCOND ||
             _calcul == ESPDECalcMode::KRIGING) &&
            _data != nullptr);
  }

private:
  const Db*                   _data;
  ESPDECalcMode               _calcul;
  PrecisionOpMultiConditional _precisionsKriging;
  PrecisionOpMultiConditional _precisionsSimu;
  std::vector<ShiftOpCs*>     _pileShiftOp;
  std::vector<PrecisionOpCs*> _pilePrecisions;
  std::vector<ProjMatrix*>    _pileProjMatrix;
  std::vector<MeshETurbo*>    _simuMeshing;
  std::vector<MeshETurbo*>    _krigingMeshing;
  Model*                      _model;
  mutable VectorVectorDouble  _workKriging;
  mutable VectorVectorDouble  _workingSimu;
  std::vector<ProjMatrix*>    _projOnDbOut;
  // query sur aproj ou // TODO ??
};
