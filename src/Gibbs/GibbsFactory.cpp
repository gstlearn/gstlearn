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
#include "Gibbs/GibbsFactory.hpp"

#include "../../include/Gibbs/GibbsUMulti.hpp"
#include "../../include/Gibbs/GibbsUMultiMono.hpp"
#include "../../include/Gibbs/GibbsUPropMono.hpp"
#include "Gibbs/GibbsMoving.hpp"
#include "Basic/AStringable.hpp"
#include "geoslib_f.h"

GibbsFactory::GibbsFactory()
{
}

GibbsFactory::~GibbsFactory()
{
}

AGibbs* GibbsFactory::createGibbs(Db* db,
                                  Model* model,
                                  Neigh* neigh,
                                  bool flag_multi_mono,
                                  bool flag_propagation)
{
  if (neigh != (Neigh *) NULL && neigh->getType() == NEIGH_MOVING)
  {
    if (flag_propagation)
    {
      messerr("Option 'flag_propagation' is incompatible with 'moving' neighborhood");
      return nullptr;
    }
    if (flag_multi_mono)
    {
      messerr("Option 'flag_multi_mono' is incompatible with 'moving' neighborhood");
      return nullptr;
    }

    // Moving Neighborhood

    GibbsMoving* gibbs = new GibbsMoving(db, model, neigh);
    return ((AGibbs *) gibbs);
  }
  else
  {

    // Unique Neighborhood

    if (model->getVariableNumber() == 1)
    {

      // Monovariate

      if (flag_propagation)
      {
        if (db->getLowerBoundNumber() >= 0 ||
            db->getUpperBoundNumber() >= 0)
        {
          messerr("The option 'flag_propagation' is incompatible with presence of Bounds");
          return nullptr;
        }

        // Propagation algorithm

        GibbsUPropMono* gibbs = new GibbsUPropMono(db, model);
        return ((AGibbs *) gibbs);
      }
      else
      {

        // Standard case

        GibbsUMultiMono* gibbs = new GibbsUMultiMono(db, model);
        return ((AGibbs *) gibbs);
      }
    }
    else
    {

      if (flag_propagation)
      {
        messerr("The option 'flag_propagation' is not compatible with 'multivariate'");
        return nullptr;
      }

      if (flag_multi_mono)
      {

        // Multi-monovariate

        GibbsUMultiMono* gibbs = new GibbsUMultiMono(db, model);
          return ((AGibbs *) gibbs);
      }
      else
      {

        // Multivariate

        GibbsUMulti* gibbs = new GibbsUMulti(db, model);
          return ((AGibbs *) gibbs);
      }
    }
  }

  messerr("No relevant option found in Gibbs Factory");
  return nullptr;
}
