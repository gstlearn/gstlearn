#pragma once

#include "Basic/NamingConvention.hpp"
#include "Db/Db.hpp"
#include "LinearOp/PrecisionOpCs.hpp"
#include "LinearOp/PrecisionOpMultiConditional.hpp"
#include "Model/Model.hpp"
#include "Model/Cova.hpp"
#include "Mesh/MeshETurbo.hpp"

#include "geoslib_enum.h"
#include <vector>

class ShiftOpCs;
class SPDE
{
public:
  SPDE();
  SPDE(Model& model,
       const Db& field,
       const Db* dat=nullptr,
       ENUM_CALCUL_MODE = CALCUL_SIMUCOND);

  void init(Model& model,
            const Db& field,
            const Db* dat=nullptr,
            ENUM_CALCUL_MODE = CALCUL_SIMUCOND);
  void compute(int nbsimus = 1, int seed = 131323) const;
  void computeKriging(const VectorDouble& vect) const;
  void computeSimuNonCond(int nbsimus = 1, int seed=131323) const;
  void computeSimuCond(int nbsimus = 1, int seed=131323) const;
  VectorDouble computeCoeffs(const VectorVectorDouble& x)const;
  void query(Db* db,NamingConvention namconv = NamingConvention("spde"));
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
  mutable VectorVectorDouble _workingSimu;
  std::vector<ProjMatrix*>   _projOnDbOut;
  int _seed;
  // query sur aproj ou
};
