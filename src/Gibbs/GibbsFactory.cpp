/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://github.com/gstlearn                                       */
/* License: BSD 3 clauses                                                     */
/*                                                                            */
/******************************************************************************/
#include "Enum/ENeigh.hpp"

#include "Gibbs/GibbsFactory.hpp"
#include "Gibbs/GibbsUMultiMono.hpp"
#include "Gibbs/GibbsUMulti.hpp"
#include "Gibbs/GibbsUPropMono.hpp"
#include "Gibbs/GibbsMMulti.hpp"
#include "Model/Model.hpp"
#include "Neigh/ANeighParam.hpp"
#include "Basic/AStringable.hpp"
#include "Db/Db.hpp"

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
 * @param flagMoving True if a Moving Neighborhood must be used
 * @return
 */
AGibbs* GibbsFactory::createGibbs(Db* db,
                                  Model* model,
                                  bool flagMoving)
{
  if (flagMoving)
  {

    // Moving Neighborhood

    GibbsMMulti* gibbs = new GibbsMMulti(db, model);
    return (static_cast<AGibbs *> (gibbs));
  }

  // Unique Neighborhood

  GibbsUMulti* gibbs = new GibbsUMulti(db, model);
  return (static_cast<AGibbs *> (gibbs));
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
      if (db->getLocNumber(ELoc::L) >= 0 || db->getLocNumber(ELoc::U) >= 0)
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
