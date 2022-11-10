/******************************************************************************/
/* COPYRIGHT ARMINES, ALL RIGHTS RESERVED                                     */
/*                                                                            */
/* THE CONTENT OF THIS WORK CONTAINS CONFIDENTIAL AND PROPRIETARY             */
/* INFORMATION OF ARMINES. ANY DUPLICATION, MODIFICATION,                     */
/* DISTRIBUTION, OR DISCLOSURE IN ANY FORM, IN WHOLE, OR IN PART, IS STRICTLY */
/* PROHIBITED WITHOUT THE PRIOR EXPRESS WRITTEN PERMISSION OF ARMINES         */
/*                                                                            */
/* TAG_SOURCE_CG                                                              */
/******************************************************************************/
#include "Basic/Utilities.hpp"
#include "Basic/VectorNumT.hpp"
#include "Basic/AException.hpp"
#include "Basic/String.hpp"
#include "Drifts/DriftFactory.hpp"
#include "Drifts/ADriftElem.hpp"
#include "Drifts/Drift1.hpp"
#include "Drifts/DriftF.hpp"
#include "Drifts/DriftX.hpp"
#include "Drifts/DriftX2.hpp"
#include "Drifts/DriftX2Y.hpp"
#include "Drifts/DriftX3.hpp"
#include "Drifts/DriftXY.hpp"
#include "Drifts/DriftXY2.hpp"
#include "Drifts/DriftXZ.hpp"
#include "Drifts/DriftY.hpp"
#include "Drifts/DriftY2.hpp"
#include "Drifts/DriftY3.hpp"
#include "Drifts/DriftYZ.hpp"
#include "Drifts/DriftZ.hpp"
#include "Drifts/DriftZ2.hpp"
#include "Drifts/DriftZ3.hpp"

#include <iostream>

ADriftElem* DriftFactory::createDriftFunc(const EDrift& type, const CovContext& ctxt, int rank_fex)
{
  ADriftElem* drift;

  switch(type.toEnum())
  {
    case EDrift::E_UC:  return new Drift1(ctxt);
    case EDrift::E_X:   return new DriftX(ctxt);
    case EDrift::E_X2:  return new DriftX2(ctxt);
    case EDrift::E_X2Y: return new DriftX2Y(ctxt);
    case EDrift::E_X3:  return new DriftX3(ctxt);
    case EDrift::E_XY:  return new DriftXY(ctxt);
    case EDrift::E_XY2: return new DriftXY2(ctxt);
    case EDrift::E_XZ:  return new DriftXZ(ctxt);
    case EDrift::E_Y:   return new DriftY(ctxt);
    case EDrift::E_Y2:  return new DriftY2(ctxt);
    case EDrift::E_Y3:  return new DriftY3(ctxt);
    case EDrift::E_YZ:  return new DriftYZ(ctxt);
    case EDrift::E_Z:   return new DriftZ(ctxt);
    case EDrift::E_Z2:  return new DriftZ2(ctxt);
    case EDrift::E_Z3:  return new DriftZ3(ctxt);
    case EDrift::E_F:
      drift = new DriftF(ctxt);
      drift->setRankFex(rank_fex);
      return drift;

    default: break;
  }
  my_throw ("Drift function not yet implemented!");
  return nullptr;
}

ADriftElem* DriftFactory::duplicateDriftFunc(const ADriftElem& drift)
{
  switch(drift.getType().toEnum())
  {
    // Warning : if a crash with "bad cast" occurs, please check the type of your Drift
    case EDrift::E_UC:  return new Drift1(   dynamic_cast<const Drift1&>   (drift));
    case EDrift::E_F:   return new DriftF(   dynamic_cast<const DriftF&>   (drift));
    case EDrift::E_X:   return new DriftX(   dynamic_cast<const DriftX&>   (drift));
    case EDrift::E_X2:  return new DriftX2(  dynamic_cast<const DriftX2&>  (drift));
    case EDrift::E_X2Y: return new DriftX2Y( dynamic_cast<const DriftX2Y&> (drift));
    case EDrift::E_X3:  return new DriftX3(  dynamic_cast<const DriftX3&>  (drift));
    case EDrift::E_XY:  return new DriftXY(  dynamic_cast<const DriftXY&>  (drift));
    case EDrift::E_XY2: return new DriftXY2( dynamic_cast<const DriftXY2&> (drift));
    case EDrift::E_XZ:  return new DriftXZ(  dynamic_cast<const DriftXZ&>  (drift));
    case EDrift::E_Y:   return new DriftY(   dynamic_cast<const DriftY&>   (drift));
    case EDrift::E_Y2:  return new DriftY2(  dynamic_cast<const DriftY2&>  (drift));
    case EDrift::E_Y3:  return new DriftY3(  dynamic_cast<const DriftY3&>  (drift));
    case EDrift::E_YZ:  return new DriftYZ(  dynamic_cast<const DriftYZ&>  (drift));
    case EDrift::E_Z:   return new DriftZ(   dynamic_cast<const DriftZ&>   (drift));
    case EDrift::E_Z2:  return new DriftZ2(  dynamic_cast<const DriftZ2&>  (drift));
    case EDrift::E_Z3:  return new DriftZ3(  dynamic_cast<const DriftZ3&>  (drift));
    default: break;
  }
  my_throw ("Drift function not yet implemented!");
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
      ADriftElem* drift = createDriftFunc(*it, ctxt);
      message("%2d - %s\n", it.getValue(), drift->getDriftName().c_str());
      delete drift;
    }
    it.toNext();
  }
}

/**
 * Return the EDrift object from the given drift symbol.
 * The symbol must correspond to one of the getDriftSymbol().
 * If the symbol doesn't exists, this method returns EDrift::UNKNOWN
 * and displays available drifts functions for the given context.
 *
 * @param symbol  Symbol of the required drift function
 * @param rank    Rank of the drift for the given symbol
 * @param ctxt    Context from which we want authorized drift functions
 */
EDrift DriftFactory::identifyDrift(const String& symbol,
                                   int* rank,
                                   const CovContext& ctxt)
{
  auto it = EDrift::getIterator();
  while (it.hasNext())
  {
    // Test drift symbol using ACovFunc::getDriftSymbol (not the EDrift keys!)
    // (This permits to ensure RGeostats scripts retro compatibility)
    if (*it != EDrift::UNKNOWN)
    {
      ADriftElem* drift = createDriftFunc(*it, ctxt);
      if (*it != EDrift::F)
      {
        String ds = toUpper(symbol);
        String dds = toUpper(drift->getDriftSymbol());
        if (ds == dds)
        {
          *rank = drift->getRankFex();
          return *it;
        }
      }
      else
      {
        int rank_loc = 0;
        if (decodeInString(drift->getDriftSymbol(), symbol, &rank_loc, false) == 0)
        {
          *rank = rank_loc-1;
          return *it;
        }
      }
      delete drift;
    }
    it.toNext();
  }
  messerr("Unknown drift function symbol:%s!", symbol.c_str());
  displayList(ctxt);
  return EDrift::UNKNOWN;
}

