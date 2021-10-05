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
#include "Basic/Law.hpp"
#include "Space/Space.hpp"
#include "Covariances/CovContext.hpp"
#include "Covariances/CovAniso.hpp"
#include "geoslib_d.h"
#include "geoslib_f.h"

/****************************************************************************/
/*!
** Main Program
**
*****************************************************************************/
int main(int argc, char *argv[])

{
  int seed     = 31415;
  int ndim     = 2;
  int nvar     = 1;
  int nbsimu   = 1;
  bool flag_moving      = true;
  bool flag_propagation = false;
  bool flag_multi_mono  = false;

  // Setup constants

  ASpaceObject::createGlobalSpace(SPACE_RN, ndim);
  law_set_random_seed(seed);
  
  // Data file

  int nech = 10;
  Db* db = new Db(nech, {0.,0.}, {1.,1.});
  VectorDouble sel;
  sel.resize(nech,1);
  sel[3] = 0;
  sel[5] = 0;
  db->addFields(sel,"Selection",LOC_SEL);
  
  // Model

  CovContext ctxt(nvar,2,1.);
  Model* model = new Model(ctxt);
  CovAniso cova(ECov::BESSEL_K,ctxt);
  cova.setRange(0.2);
  cova.setParam(1.);
  VectorDouble sills;
  if (nvar == 2)
    sills = {1., 0.2, 0.2, 2.};
  else
    sills = {1.5};
  cova.setSill(sills);
  model->addCova(&cova);
  model->display();

  // Neighborhood

  Neigh* neigh = nullptr;
  if (flag_moving)
  {
    int nmaxi = 20;
    double radius = 1.;
    neigh = new Neigh(ndim, nmaxi, radius);
    neigh->display();
  }

  // Gibbs

  int niter    = 1000;
  int nburn    = 10;
  bool verbose = false;
  gibbs_sampler(db, model, neigh, nbsimu, seed, nburn, niter, false,
                flag_multi_mono, flag_propagation, 2,
                5., EPSILON3, false, false, verbose);
  db->displayMore(FLAG_STATS);

  // Cleaning structures

  db    = db_delete(db);
  model = model_free(model);
  neigh = neigh_free(neigh);
  return(0);
}
