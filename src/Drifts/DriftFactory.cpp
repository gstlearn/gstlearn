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
#include "Drifts/DriftM.hpp"
#include "Drifts/DriftF.hpp"

#include <iostream>

ADriftElem* DriftFactory::createDriftByType(const EDrift &type,
                                            int rank_fex,
                                            const CovContext &ctxt)
{
  switch(type.toEnum())
  {
    case EDrift::E_UC:  return new DriftM(VectorInt(), 1., VectorDouble(), ctxt);
    case EDrift::E_X:   return new DriftM({1},         0., VectorDouble(), ctxt);
    case EDrift::E_X2:  return new DriftM({2},         0., VectorDouble(), ctxt);
    case EDrift::E_X2Y: return new DriftM({2,1},       0., VectorDouble(), ctxt);
    case EDrift::E_X3:  return new DriftM({3},         0., VectorDouble(), ctxt);
    case EDrift::E_XY:  return new DriftM({1,1},       0., VectorDouble(), ctxt);
    case EDrift::E_XY2: return new DriftM({1,2},       0., VectorDouble(), ctxt);
    case EDrift::E_XZ:  return new DriftM({1,0,1},     0., VectorDouble(), ctxt);
    case EDrift::E_Y:   return new DriftM({0,1},       0., VectorDouble(), ctxt);
    case EDrift::E_Y2:  return new DriftM({0,2},       0., VectorDouble(), ctxt);
    case EDrift::E_Y3:  return new DriftM({0,3},       0., VectorDouble(), ctxt);
    case EDrift::E_YZ:  return new DriftM({0,1,1},     0., VectorDouble(), ctxt);
    case EDrift::E_Z:   return new DriftM({0,0,1},     0., VectorDouble(), ctxt);
    case EDrift::E_Z2:  return new DriftM({0,0,2},     0., VectorDouble(), ctxt);
    case EDrift::E_Z3:  return new DriftM({0,0,3},     0., VectorDouble(), ctxt);
    case EDrift::E_F:   return new DriftF(rank_fex, ctxt);
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
  if (ds == "1")    return new DriftM(VectorInt(), 1., VectorDouble(), ctxt);
  if (ds == "X")    return new DriftM({1},         0., VectorDouble(), ctxt);
  if (ds == "X2")   return new DriftM({2},         0., VectorDouble(), ctxt);
  if (ds == "X2Y")  return new DriftM({2,1},       0., VectorDouble(), ctxt);
  if (ds == "X3")   return new DriftM({3},         0., VectorDouble(), ctxt);
  if (ds == "XY")   return new DriftM({1,1},       0., VectorDouble(), ctxt);
  if (ds == "XY2")  return new DriftM({1,2},       0., VectorDouble(), ctxt);
  if (ds == "XZ")   return new DriftM({1,0,1},     0., VectorDouble(), ctxt);
  if (ds == "Y")    return new DriftM({0,1},       0., VectorDouble(), ctxt);
  if (ds == "Y2")   return new DriftM({0,2},       0., VectorDouble(), ctxt);
  if (ds == "Y3")   return new DriftM({0,3},       0., VectorDouble(), ctxt);
  if (ds == "YZ")   return new DriftM({0,1,1},     0., VectorDouble(), ctxt);
  if (ds == "Z")    return new DriftM({0,0,1},     0., VectorDouble(), ctxt);
  if (ds == "Z2")   return new DriftM({0,0,2},     0., VectorDouble(), ctxt);
  if (ds == "Z3")   return new DriftM({0,0,3},     0., VectorDouble(), ctxt);
  if (ds == "F")
  {
      drift = new DriftF(0, ctxt);

      int rank_fex = 0;
      if (decodeInString("F", symbol, &rank_fex, false) == 0)
        rank_fex = rank_fex-1;
      drift->setRankFex(rank_fex);
      return drift;
  }

  message("Drift Symbol %s is unknown", symbol.c_str());
  return nullptr;
}

/**
 * Prints the list of Drift functions available
 */
void DriftFactory::displayList(const CovContext& ctxt)
{
  message("List of authorized Drift Functions:\n");
  auto it = EDrift::getIterator();
  while (it.hasNext())
  {
    if (*it != EDrift::UNKNOWN)
    {
      ADriftElem* drift = createDriftByType(*it, 0, ctxt);
      message("%2d - %s\n", it.getValue(), drift->getDriftName().c_str());
      delete drift;
    }
    it.toNext();
  }
}
