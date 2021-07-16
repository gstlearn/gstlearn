#pragma once

#include "Db/Db.hpp"
#include "LinearOp/PrecisionOpMultiConditional.hpp"
#include "Model/Model.hpp"
#include "Model/Cova.hpp"
#include "Mesh/MeshETurbo.hpp"

class SPDE
{
public:
  SPDE(Model& model,
       const Db& field,
       ANoStat* nostat=nullptr,
       const Db* dat=nullptr);

  void init(Model& model,
            const Db& field,
            ANoStat* nostat,
            const Db* dat=nullptr);

  MeshETurbo createMeshing(const CovAniso& cova, const Db& field,double discr);
  virtual ~SPDE();
private:
  PrecisionOpMultiConditional _precisionsKriging;
  std::vector<PrecisionOp> _precistionLists;
};
