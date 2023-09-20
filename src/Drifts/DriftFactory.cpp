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
#include "Basic/Utilities.hpp"
#include "Basic/VectorNumT.hpp"
#include "Basic/AException.hpp"
#include "Basic/String.hpp"
#include "Drifts/DriftFactory.hpp"
#include "Drifts/ADriftElem.hpp"
#include "Drifts/DriftList.hpp"
#include "Drifts/DriftM.hpp"
#include "Drifts/DriftF.hpp"

#include <iostream>

/**
 * This Drift identification is used for interpreting old serialized files
 * where the drift function was encoded by its rank.
 * @param rank     Rank of the drift function (in the deprecated EDrift Enum)
 * @param rank_fex Rank of the External Drift variable
 * @param ctxt     Cov_Context()
 * @return
 */
ADriftElem* DriftFactory::createDriftByRank(int rank,
                                            int rank_fex,
                                            const CovContext &ctxt)
{
  switch(rank)
  {
    case 0:  return new DriftM(VectorInt(), ctxt);
    case 1:  return new DriftM({1},         ctxt);
    case 2:  return new DriftM({0,1},       ctxt);
    case 3:  return new DriftM({0,0,1},     ctxt);
    case 4:  return new DriftM({2},         ctxt);
    case 5:  return new DriftM({0,2},       ctxt);
    case 6:  return new DriftM({1,1},       ctxt);
    case 7:  return new DriftM({0,0,2},     ctxt);
    case 8:  return new DriftM({1,0,1},     ctxt);
    case 9:  return new DriftM({0,1,1},     ctxt);
    case 10: return new DriftM({3},         ctxt);
    case 11: return new DriftM({2,1},       ctxt);
    case 12: return new DriftM({1,2},       ctxt);
    case 13: return new DriftM({0,3},       ctxt);
    case 14: return new DriftM({0,0,3},     ctxt);
    case 15: return new DriftF(rank_fex,    ctxt);
    default: break;
  }
  my_throw ("Drift function not yet implemented!");
  return nullptr;
}

/**
 * Create a Drift Item defined by a symbol
 * This function is left for compatibility with RGeostats code.
 * @param symbol Name of the symbol
 * @param ctxt   CovContext which specifies the space dimension and number of variables
 * @return
 */
ADriftElem* DriftFactory::createDriftBySymbol(const String &symbol,
                                              const CovContext &ctxt)
{
  ADriftElem* drift;

  String ds = toUpper(symbol);
  if (ds == "1")    return new DriftM(VectorInt(), ctxt);
  if (ds == "X")    return new DriftM({1},         ctxt);
  if (ds == "X2")   return new DriftM({2},         ctxt);
  if (ds == "X2Y")  return new DriftM({2,1},       ctxt);
  if (ds == "X3")   return new DriftM({3},         ctxt);
  if (ds == "XY")   return new DriftM({1,1},       ctxt);
  if (ds == "XY2")  return new DriftM({1,2},       ctxt);
  if (ds == "XZ")   return new DriftM({1,0,1},     ctxt);
  if (ds == "Y")    return new DriftM({0,1},       ctxt);
  if (ds == "Y2")   return new DriftM({0,2},       ctxt);
  if (ds == "Y3")   return new DriftM({0,3},       ctxt);
  if (ds == "YZ")   return new DriftM({0,1,1},     ctxt);
  if (ds == "Z")    return new DriftM({0,0,1},     ctxt);
  if (ds == "Z2")   return new DriftM({0,0,2},     ctxt);
  if (ds == "Z3")   return new DriftM({0,0,3},     ctxt);
  if (ds == "F")
  {
      int rank_fex = 0;
      if (decodeInString("F", symbol, &rank_fex, false) == 0)
        rank_fex = rank_fex-1;
      drift = new DriftF(rank_fex, ctxt);
      return drift;
  }

  message("Drift Symbol %s is unknown", symbol.c_str());
  return nullptr;
}

ADriftElem* DriftFactory::createDriftByIdentifier(const String& driftname,
                                                  const CovContext &ctxt)
{
  // Look for Universality Condition
  ADriftElem* drift = new DriftM(VectorInt(), ctxt);
  String locname = drift->getDriftName();
  if (driftname.find(locname) == 0) return drift;
  delete drift;

  // Look for an External Drift (testing the 5 first variable ranks)
  for (int itest = 0; itest < 5; itest++)
  {
    ADriftElem* drift = new DriftF(itest, ctxt);
    String locname = drift->getDriftName();
    if (driftname.find(locname) == 0) return drift;
    delete drift;
  }

  // Looking for a standard monomial drift
  // For simplicity, we only check the 14 first possibilities defined in the function createDriftByRank
  for (int itest = 0; itest <= 14; itest++)
  {
    ADriftElem* drift = DriftFactory::createDriftByRank(itest,0,ctxt);
    String locname = drift->getDriftName();
    if (driftname.find(locname) == 0) return drift;
    delete drift;
  }

  messerr("Error: Drift Name(%s) is unknown", driftname);
  return nullptr;
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
 * @Note: this function is limited to order<=2 and ndim<= 3
 */
DriftList* DriftFactory::createDriftListFromIRF(int order,
                                                int nfex,
                                                const CovContext &ctxt)
{
  DriftList* drifts = new DriftList();
  int ndim = ctxt.getNDim();

  // Standard monomials
  switch (order)
  {
    case -1:
      // In the strict stationary case, no drift is defined (even external)
      return drifts;
      break;

    case 0:
      drifts->addDrift(new DriftM(VectorInt(), ctxt));                        // 1
      break;

    case 1:
      drifts->addDrift(new DriftM(VectorInt(), ctxt));                        // 1
      if (ndim >= 1)
        drifts->addDrift(new DriftM(VectorInt({1}), ctxt));                   // X
      if (ndim >= 2)
        drifts->addDrift(new DriftM(VectorInt({0,1}), ctxt));                 // Y
      if (ndim >= 3)
        drifts->addDrift(new DriftM(VectorInt({0,0,1}), ctxt));               // Z
      break;

    case 2:
      drifts->addDrift(new DriftM());                                         // 1
      if (ndim >= 1)
      {
        drifts->addDrift(new DriftM(VectorInt({1}), ctxt));                   // X
        drifts->addDrift(new DriftM(VectorInt({2}), ctxt));                   // X^2
      }
      if (ndim >= 2)
      {
        drifts->addDrift(new DriftM(VectorInt({0,1}), ctxt));                 // Y
        drifts->addDrift(new DriftM(VectorInt({1,1}), ctxt));                 // YX
        drifts->addDrift(new DriftM(VectorInt({0,2}), ctxt));                 // Y^2
      }
      if (ndim >= 3)
      {
        drifts->addDrift(new DriftM(VectorInt({0,0,1}), ctxt));               // Z
        drifts->addDrift(new DriftM(VectorInt({1,0,1}), ctxt));               // ZX
        drifts->addDrift(new DriftM(VectorInt({0,1,1}), ctxt));               // ZY
        drifts->addDrift(new DriftM(VectorInt({0,0,2}), ctxt));               // Z^2
      }
  }

  if (nfex > 0)
  {
    // Adding the external drift(s)
    for (int ifex = 0; ifex < nfex; ifex++)
      drifts->addDrift(new DriftF(ifex, ctxt));
  }

  drifts->updateDriftList();
  return drifts;
}

/**
 * Create the list of drift functions for a Cokriging system
 * starting from a list of initial drift functions
 * It is implemented for ndim=1 or 2
 * @param Input list of drift funcitons
 * @param ctxt CovContext structure
 * @return
 */
DriftList* DriftFactory::createDriftListForGradients(const DriftList& inputlist, const CovContext &ctxt)
{
  DriftList* newdrifts = new DriftList();
  int ndim = ctxt.getNDim();

  int ndrift = driftlist.getDriftNumber();

  drifts.setFlagLinked(flag_linked);
  for (int il = 0; il < ndrift; il++)
  {
    ADriftElem* drft = dynamic_cast<ADriftElem*>(inputlist.getDrift(il)->clone());
    String drift_name = drft->getDriftName();

    // Creating the Z variable

    newdrifts->addDrift(drift)
  }

  for (int il = 0; il < nbfl; il++)
  {
    ADriftElem* drft = dynamic_cast<ADriftElem*>(model->getDrift(il)->clone());
    drft->setCtxt(ctxt);
    drifts.addDrift(drft);
    drifts.setFiltered(il, model->isDriftFiltered(il));
  }
  new_model->setDriftList(&drifts);
  for (int il = 0; il < nbfl; il++)
  {
    new_model->setDriftFiltered(il, model->isDriftFiltered(il));
  }

  // Update the drift for the derivatives


  if (ndim == 1)
  {

  }
}
