#pragma once

#include "Db/Db.hpp"
#include "LinearOp/PrecisionOpCs.hpp"
#include "LinearOp/PrecisionOpMultiConditional.hpp"
#include "Model/Model.hpp"
#include "Model/Cova.hpp"
#include "Mesh/MeshETurbo.hpp"

class ShiftOpCs;
class SPDE
{
public:
  SPDE(Model& model,
       const Db& field,
       const Db* dat=nullptr);

  void init(Model& model,
            const Db& field,
            const Db* dat=nullptr);

  MeshETurbo* createMeshing(const CovAniso& cova,
                            const Db& field,
                            double discr,
                            double ext = 0.);
  virtual ~SPDE();

private:
  std::vector<ShiftOpCs*>     _pileShiftOp;
  PrecisionOpMultiConditional _precisionsKriging;
  std::vector<PrecisionOpCs>  _precistionLists;

  Model* _model;
};
