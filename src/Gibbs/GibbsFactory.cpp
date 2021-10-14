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
#include "Gibbs/GibbsUMultiMono.hpp"
#include "Gibbs/GibbsUMulti.hpp"
#include "Gibbs/GibbsUPropMono.hpp"
#include "Gibbs/GibbsMMulti.hpp"
#include "Model/Model.hpp"
#include "Neigh/Neigh.hpp"
#include "Neigh/ENeigh.hpp"
#include "Basic/AStringable.hpp"
#include "geoslib_f.h"

GibbsFactory::GibbsFactory()
{
}

GibbsFactory::~GibbsFactory()
{
}

/**
 * Create the relevant Gibbs with Multivariate complete model
 * @param db     Db structure
 * @param model  Multivariate structure
 * @param neigh  Neighborhood structure
 * @return
 */
AGibbs* GibbsFactory::createGibbs(Db* db,
                                  Model* model,
                                  Neigh* neigh)
{
  if (neigh != (Neigh *) NULL && neigh->getType() == ENeigh::MOVING)
  {

    // Moving Neighborhood

    GibbsMMulti* gibbs = new GibbsMMulti(db, model, neigh);
    return (static_cast<AGibbs *> (gibbs));
  }
  else
  {

    // Unique Neighborhood

    GibbsUMulti* gibbs = new GibbsUMulti(db, model);
    return (static_cast<AGibbs *> (gibbs));
  }

  messerr("No relevant option found in Gibbs Factory");
  return nullptr;
}

/**
 * Create the Gibbs instance in the case of Multi-Mono model
 * @param db     Db structure
 * @param models Vector of monovariate models
 * @param rho    Correlation coefficient (current to first model)
 * @param flag_propagation Propagation flag
 * @return
 */
AGibbs* GibbsFactory::createGibbs(Db* db,
                                  std::vector<Model *> models,
                                  double rho,
                                  bool flag_propagation)
{

  // Unique Neighborhood

  if (models.size() == 1)
  {

    // Monovariate

    if (flag_propagation)
    {
      if (db->getLowerBoundNumber() >= 0 || db->getUpperBoundNumber() >= 0)
      {
        messerr(
            "The option 'flag_propagation' is incompatible with presence of Bounds");
        return nullptr;
      }

      // Propagation algorithm

      GibbsUPropMono* gibbs = new GibbsUPropMono(db, models, 1.);
      return (static_cast<AGibbs *>(gibbs));
    }
    else
    {

      // Standard case

      GibbsUMultiMono* gibbs = new GibbsUMultiMono(db, models, rho);
      return (static_cast<AGibbs *>(gibbs));
    }
  }
  else
  {

    if (flag_propagation)
    {
      messerr(
          "The option 'flag_propagation' is not compatible with 'multivariate'");
      return nullptr;
    }

    GibbsUMultiMono* gibbs = new GibbsUMultiMono(db, models, rho);
    return (static_cast<AGibbs *>(gibbs));
  }

  messerr("No relevant option found in Gibbs Factory");
  return nullptr;
}
