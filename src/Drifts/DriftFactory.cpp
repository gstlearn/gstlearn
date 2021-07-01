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

#include "geoslib_f.h"

#include <iostream>

#include "Basic/Utilities.hpp"
#include "Basic/Vector.hpp"
#include "Basic/AException.hpp"
#include "Basic/String.hpp"

ADriftElem* DriftFactory::createDriftFunc(const ENUM_DRIFTS& type, const CovContext& ctxt)
{
  switch(type)
  {
    case DRIFT_1:   return new Drift1(ctxt);
    case DRIFT_F:   return new DriftF(ctxt);
    case DRIFT_X:   return new DriftX(ctxt);
    case DRIFT_X2:  return new DriftX2(ctxt);
    case DRIFT_X2Y: return new DriftX2Y(ctxt);
    case DRIFT_X3:  return new DriftX3(ctxt);
    case DRIFT_XY:  return new DriftXY(ctxt);
    case DRIFT_XY2: return new DriftXY2(ctxt);
    case DRIFT_XZ:  return new DriftXZ(ctxt);
    case DRIFT_Y:   return new DriftY(ctxt);
    case DRIFT_Y2:  return new DriftY2(ctxt);
    case DRIFT_Y3:  return new DriftY3(ctxt);
    case DRIFT_YZ:  return new DriftYZ(ctxt);
    case DRIFT_Z:   return new DriftZ(ctxt);
    case DRIFT_Z2:  return new DriftZ2(ctxt);
    case DRIFT_Z3:  return new DriftZ3(ctxt);
    default:
      std::cout << "Error unknown drift !" << std::endl; break;
  }
  my_throw ("Covariance function not yet implemented !");
  return nullptr;
}

ADriftElem* DriftFactory::duplicateDriftFunc(const ADriftElem& drift)
{
  switch(drift.getType())
  {
    // Warning : if a crash with "bad cast" occurs, please check the type of your CovFunc
    case DRIFT_1:   return new Drift1(   dynamic_cast<const Drift1&>   (drift));
    case DRIFT_F:   return new DriftF(   dynamic_cast<const DriftF&>   (drift));
    case DRIFT_X:   return new DriftX(   dynamic_cast<const DriftX&>   (drift));
    case DRIFT_X2:  return new DriftX2(  dynamic_cast<const DriftX2&>  (drift));
    case DRIFT_X2Y: return new DriftX2Y( dynamic_cast<const DriftX2Y&> (drift));
    case DRIFT_X3:  return new DriftX3(  dynamic_cast<const DriftX3&>  (drift));
    case DRIFT_XY:  return new DriftXY(  dynamic_cast<const DriftXY&>  (drift));
    case DRIFT_XY2: return new DriftXY2( dynamic_cast<const DriftXY2&> (drift));
    case DRIFT_XZ:  return new DriftXZ(  dynamic_cast<const DriftXZ&>  (drift));
    case DRIFT_Y:   return new DriftY(   dynamic_cast<const DriftY&>   (drift));
    case DRIFT_Y2:  return new DriftY2(  dynamic_cast<const DriftY2&>  (drift));
    case DRIFT_Y3:  return new DriftY3(  dynamic_cast<const DriftY3&>  (drift));
    case DRIFT_YZ:  return new DriftYZ(  dynamic_cast<const DriftYZ&>  (drift));
    case DRIFT_Z:   return new DriftZ(   dynamic_cast<const DriftZ&>   (drift));
    case DRIFT_Z2:  return new DriftZ2(  dynamic_cast<const DriftZ2&>  (drift));
    case DRIFT_Z3:  return new DriftZ3(  dynamic_cast<const DriftZ3&>  (drift));
    default:
      break;
  }
  message("Drift type = %d\n",drift.getType());
  my_throw ("Drift function not yet implemented !");
  return nullptr;
}

/**
 * Prints the list of Drift functions available
 */
void DriftFactory::displayList(const CovContext& ctxt)
{
  message("List of authorized Drift Functions:\n");
  for (int i = 0; i < DRIFT_NUMBER; i++)
  {
    ADriftElem* drift = createDriftFunc((ENUM_DRIFTS) i, ctxt);
    message("%2d - %s\n", i + 1, drift->getDriftName().c_str());
    delete drift;
  }
}

int DriftFactory::identifyDrift(const String& symbol, ENUM_DRIFTS *type, int *rank,
                                const CovContext& ctxt)
{
  int rank_loc;

  for (int i = 0; i < DRIFT_NUMBER; i++)
  {
    ADriftElem* drift = createDriftFunc((ENUM_DRIFTS) i, ctxt);

    if ((ENUM_DRIFTS) i != DRIFT_F)
    {
      if (matchRegexp(drift->getDriftSymbol(),symbol,false))
      {
        *type = (ENUM_DRIFTS) i;
        *rank = drift->getRankFex();
        return 0;
      }
    }
    else
    {
      if (decodeInString(drift->getDriftSymbol(), symbol, &rank_loc, false) == 0)
      {
        *type = (ENUM_DRIFTS) i;
        *rank = rank_loc-1;
        return 0;
      }
    }
    delete drift;
  }
  messerr("Unknown Drift Symbol : %s",symbol.c_str());
  DriftFactory::displayList(ctxt);
  return 1;
}

int DriftFactory::getDriftNumber()
{
  return DRIFT_NUMBER;
}
