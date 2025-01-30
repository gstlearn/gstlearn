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
#include "Basic/VectorNumT.hpp"
#include "Basic/AException.hpp"
#include "Basic/String.hpp"
#include "Drifts/DriftFactory.hpp"
#include "Drifts/ADrift.hpp"
#include "Drifts/DriftList.hpp"
#include "Drifts/DriftM.hpp"
#include "Drifts/DriftF.hpp"

/**
 * This Drift identification is used for interpreting old serialized files
 * where the drift function was encoded by its rank.
 * @param rank     Rank of the drift function (in a deprecated Enum)
 * @param rank_fex Rank of the External Drift variable
 * @return
 */
ADrift* DriftFactory::createDriftByRank(int rank, int rank_fex)
{
  switch(rank)
  {
    case 0:  return new DriftM(VectorInt());
    case 1:  return new DriftM(VectorInt({1}));
    case 2:  return new DriftM({0,1});
    case 3:  return new DriftM({0,0,1});
    case 4:  return new DriftM(VectorInt({2}));
    case 5:  return new DriftM({0,2});
    case 6:  return new DriftM({1,1});
    case 7:  return new DriftM({0,0,2});
    case 8:  return new DriftM({1,0,1});
    case 9:  return new DriftM({0,1,1});
    case 10: return new DriftM(VectorInt({3}));
    case 11: return new DriftM({2,1});
    case 12: return new DriftM({1,2});
    case 13: return new DriftM({0,3});
    case 14: return new DriftM({0,0,3});
    case 15: return new DriftF(rank_fex);
    default: break;
  }
  my_throw_impossible("Drift function not yet implemented!");
  return nullptr;
}

/**
 * Create a Drift Item defined by a symbol
 * This function is left for compatibility with RGeostats code.
 * @param symbol Name of the symbol
 * @return
 */
ADrift* DriftFactory::createDriftBySymbol(const String &symbol)
{
  ADrift* drift;

  String ds = toUpper(symbol);
  if (ds == "1") return new DriftM(VectorInt());
  if (ds == "X") return new DriftM(VectorInt( { 1 }));
  if (ds == "X2") return new DriftM(VectorInt( { 2 }));
  if (ds == "X2Y") return new DriftM( { 2, 1 });
  if (ds == "X3") return new DriftM(VectorInt( { 3 }));
  if (ds == "XY") return new DriftM( { 1, 1 });
  if (ds == "XY2") return new DriftM( { 1, 2 });
  if (ds == "XZ") return new DriftM( { 1, 0, 1 });
  if (ds == "Y") return new DriftM( { 0, 1 });
  if (ds == "Y2") return new DriftM( { 0, 2 });
  if (ds == "Y3") return new DriftM( { 0, 3 });
  if (ds == "YZ") return new DriftM( { 0, 1, 1 });
  if (ds == "Z") return new DriftM( { 0, 0, 1 });
  if (ds == "Z2") return new DriftM( { 0, 0, 2 });
  if (ds == "Z3") return new DriftM( { 0, 0, 3 });
  if (ds == "F")
  {
      int rank_fex = 0;
      if (decodeInString("F", symbol, &rank_fex, false) == 0)
        rank_fex = rank_fex-1;
      drift = new DriftF(rank_fex);
      return drift;
  }

  message("Drift Symbol %s is unknown", symbol.c_str());
  return nullptr;
}

ADrift* DriftFactory::createDriftByIdentifier(const String& driftname)
{
  ADrift* drift;

  // Look for a standard monomial drift
  drift = DriftM::createByIdentifier(driftname);

  // Look for an external drift
  if (drift == nullptr)
    drift = DriftF::createByIdentifier(driftname);

  if (drift == nullptr)
    messerr("Error: Drift Name(%s) is unknown", driftname.c_str());
  return drift;
}

/**
 * Creating the list of Drift functions correspondaing to the following constraints:
 * - Rank of the IRF
 * - Number of external drift functions
 * @param order Rank of the IRF
 * @param nfex  Number of external drift functions
 * @param ctxt  Cov_context
 * @return
 *
 * @remarks: this function is limited to order<=2 and ndim<= 3
 */
DriftList* DriftFactory::createDriftListFromIRF(int order,
                                                int nfex,
                                                const CovContext &ctxt)
{
  DriftList* drifts = new DriftList(ctxt);
  int ndim = ctxt.getNDim();

  // Standard monomials
  switch (order)
  {
    case -1:
      // In the strict stationary case, no drift is defined (even external)
      return drifts;
      break;

    case 0:
      drifts->addDrift(new DriftM(VectorInt()));                        // 1
      break;

    case 1:
      drifts->addDrift(new DriftM(VectorInt()));                        // 1
      if (ndim >= 1)
        drifts->addDrift(new DriftM(VectorInt({1})));                   // X
      if (ndim >= 2)
        drifts->addDrift(new DriftM(VectorInt({0,1})));                 // Y
      if (ndim >= 3)
        drifts->addDrift(new DriftM(VectorInt({0,0,1})));               // Z
      break;

    case 2:
      drifts->addDrift(new DriftM());                                         // 1
      if (ndim >= 1)
      {
        drifts->addDrift(new DriftM(VectorInt({1})));                   // X
        drifts->addDrift(new DriftM(VectorInt({2})));                   // X^2
      }
      if (ndim >= 2)
      {
        drifts->addDrift(new DriftM(VectorInt({0,1})));                 // Y
        drifts->addDrift(new DriftM(VectorInt({1,1})));                 // YX
        drifts->addDrift(new DriftM(VectorInt({0,2})));                 // Y^2
      }
      if (ndim >= 3)
      {
        drifts->addDrift(new DriftM(VectorInt({0,0,1})));               // Z
        drifts->addDrift(new DriftM(VectorInt({1,0,1})));               // ZX
        drifts->addDrift(new DriftM(VectorInt({0,1,1})));               // ZY
        drifts->addDrift(new DriftM(VectorInt({0,0,2})));               // Z^2
      }
  }

  if (nfex > 0)
  {
    // Adding the external drift(s)
    for (int ifex = 0; ifex < nfex; ifex++)
      drifts->addDrift(new DriftF(ifex));
  }

  drifts->resetDriftList();
  return drifts;
}

/**
 * Create the list of drift functions for a Cokriging system
 * starting from a list of initial drift functions
 * It is implemented for ndim=1 or 2
 * @param olddrifts list of drift funcitons
 * @param ctxt CovContext structure
 * @return
 */
DriftList* DriftFactory::createDriftListForGradients(const DriftList* olddrifts, const CovContext &ctxt)
{
  DriftM drft;
  DriftList* newdrifts = new DriftList(ctxt);
  newdrifts->setFlagLinked(true);
  int ndim  = ctxt.getNDim();
  int order = olddrifts->getDriftMaxIRFOrder();

  if (olddrifts->hasExternalDrift())
  {
    messerr("This method is not valid when an External Drift is present");
    return newdrifts;
  }
  if (ndim != 1 && ndim != 2)
  {
    messerr("This method is limited to 2 or 3 space dimension");
    return newdrifts;
  }
  if (order > 2)
  {
    messerr("This method is limited to order <= 2");
    return newdrifts;
  }

  if (order == -1) return newdrifts;

  // Universality condition
  drft = DriftM(VectorInt());
  newdrifts->addDrift(&drft);

  if (order >= 1)
  {
    // Order-1 drift terms
    if (ndim >= 1)
    {
      drft = DriftM(VectorInt({1}));
      newdrifts->addDrift(&drft);
    }
    if (ndim >= 2)
    {
      drft = DriftM(VectorInt({0,1}));
      newdrifts->addDrift(&drft);
    }
  }
  if (order >= 2)
  {
    // Order-2 drift terms
    if (ndim >= 1)
    {
      drft = DriftM(VectorInt({2}));
      newdrifts->addDrift(&drft);
    }
    if (ndim >= 2)
    {
      drft = DriftM(VectorInt({0,2}));
      newdrifts->addDrift(&drft);
      drft = DriftM(VectorInt({1,1}));
      newdrifts->addDrift(&drft);
    }
  }

  // Updating the auxiliary arrays
  newdrifts->resetDriftList();

  // Defining the coefficients

  if (order == 1)
  {
    // Order-1 drift terms
    if (ndim == 1)
    {
      newdrifts->setDriftCLByPart(1, 0, {0, 0}); // d(X) / dX
      newdrifts->setDriftCLByPart(1, 1, {1, 0}); // d(X) / dX
    }
    if (ndim == 2)
    {
      newdrifts->setDriftCLByPart(1, 0, {0, 0, 0}); // d(1) / dX
      newdrifts->setDriftCLByPart(1, 1, {1, 0, 0}); // d(X) / dX
      newdrifts->setDriftCLByPart(1, 2, {0, 0, 0}); // d(Y) / dX

      newdrifts->setDriftCLByPart(2, 0, {0, 0, 0}); // d(1) / dY
      newdrifts->setDriftCLByPart(2, 1, {0, 0, 0}); // d(X) / dY
      newdrifts->setDriftCLByPart(2, 2, {1, 0, 0}); // d(Y) / dY
    }
  }
  if (order == 2)
  {
    // Order-2 drift terms
    if (ndim == 1)
    {
      newdrifts->setDriftCLByPart(1, 0, {0, 0, 0}); // d(1) / dX
      newdrifts->setDriftCLByPart(1, 1, {1, 0, 0}); // d(X) / dX
      newdrifts->setDriftCLByPart(1, 1, {0, 2, 0}); // d(X^2) / dX
    }
    if (ndim == 2)
    {
      newdrifts->setDriftCLByPart(1, 0, {0, 0, 0, 0, 0, 0}); // d(1) / dX
      newdrifts->setDriftCLByPart(1, 1, {1, 0, 0, 0, 0, 0}); // d(X) / dX
      newdrifts->setDriftCLByPart(1, 2, {0, 0, 0, 0, 0, 0}); // d(Y) / dX
      newdrifts->setDriftCLByPart(1, 3, {0, 2, 0, 0, 0, 0}); // d(X^2) / dX
      newdrifts->setDriftCLByPart(1, 4, {0, 0, 0, 0, 0, 0}); // d(Y^2) / dX
      newdrifts->setDriftCLByPart(1, 5, {0, 0, 1, 0, 0, 0}); // d(XY) / dX

      newdrifts->setDriftCLByPart(2, 0, {0, 0, 0, 0, 0, 0}); // d(1) / dY
      newdrifts->setDriftCLByPart(2, 1, {0, 0, 0, 0, 0, 0}); // d(X) / dY
      newdrifts->setDriftCLByPart(2, 2, {1, 0, 0, 0, 0, 0}); // d(Y) / dY
      newdrifts->setDriftCLByPart(2, 3, {0, 0, 0, 0, 0, 0}); // d(X^2) / dY
      newdrifts->setDriftCLByPart(2, 4, {0, 0, 2, 0, 0, 0}); // d(Y^2) / dY
      newdrifts->setDriftCLByPart(2, 5, {0, 1, 0, 0, 0, 0}); // d(XY) / dY
    }
  }

  // Copy the 'filter' status
  for (int il = 0, ndrift = olddrifts->getNDrift(); il < ndrift; il++)
  {
    newdrifts->setFiltered(il, olddrifts->isDriftFiltered(il));
  }

  return newdrifts;
}
