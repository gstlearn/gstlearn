#pragma once

#include "Db/Db.hpp"
#include "LinearOp/PrecisionOpCs.hpp"
#include "LinearOp/PrecisionOpMultiConditional.hpp"
#include "Model/Model.hpp"
#include "Model/Cova.hpp"
#include "Mesh/MeshETurbo.hpp"

#include "geoslib_enum.h"

class ShiftOpCs;
class SPDE
{
public:
  SPDE(Model& model,
       const Db& field,
       const Db* dat=nullptr,
       ENUM_CALCUL_MODE = CALCUL_SIMUCOND);

  void init(Model& model,
            const Db& field,
            const Db* dat=nullptr,
            ENUM_CALCUL_MODE = CALCUL_SIMUCOND);
  void computeKriging() const;
  void query(Db* db);
  virtual ~SPDE();

private:
  void _purge();
  MeshETurbo* _createMeshing(const CovAniso& cova,
                              const Db& field,
                              double discr,
                              double ext = 0.);
  bool _calculSimu()const{return _calcul == CALCUL_SIMUCOND || _calcul == CALCUL_SIMUNONCOND;}
  bool _calculKriging()const
  {
    return (_calcul == CALCUL_SIMUCOND || _calcul == CALCUL_KRIGING) && _data!=nullptr;
  }

private:
  const Db*                   _data;
  ENUM_CALCUL_MODE            _calcul;
  PrecisionOpMultiConditional _precisionsKriging;
  PrecisionOpMultiConditional _precisionsSimu;
  std::vector<ShiftOpCs*>     _pileShiftOp;
  std::vector<PrecisionOpCs*> _pilePrecisions;
  std::vector<ProjMatrix*>    _pileProjMatrix;
  std::vector<MeshETurbo*>    _simuMeshing;
  std::vector<MeshETurbo*>    _krigingMeshing;
  Model* _model;
  mutable VectorVectorDouble _workKriging;
  std::vector<ProjMatrix*>   _projOnDbOut;
  // query sur aproj ou
};
