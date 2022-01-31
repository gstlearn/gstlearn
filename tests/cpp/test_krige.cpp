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
#include "geoslib_d.h"
#include "geoslib_f.h"
#include "Db/Db.hpp"
#include "Model/Model.hpp"
#include "Covariances/CovAniso.hpp"
#include "Covariances/CovLMC.hpp"
#include "Basic/Law.hpp"
#include "Neigh/ANeighParam.hpp"
#include "Neigh/NeighUnique.hpp"

/****************************************************************************/
/*!
 ** Main Program
 **
 *****************************************************************************/
int main(int /*argc*/, char */*argv*/[])

{
  // Global parameters
  int ndim = 2;
  int nvar = 1;

  setup_license("Demonstration");
  ASpaceObject::defineDefaultSpace(SPACE_RN, ndim);

  // Generate the output grid
  VectorInt nx = {100,100};
  DbGrid* grid = DbGrid::create(nx);
  grid->display();

  // Generate the data base
  int nech = 100;
  VectorDouble tab;
  // Coordinates
  for (int idim=0; idim<ndim; idim++)
    for (int iech=0; iech<nech; iech++)
    {
      tab.push_back(law_uniform(0,100));
    }
  // Variable
  for (int ivar=0; ivar<nvar; ivar++)
    for (int iech=0; iech<nech; iech++)
    {
      tab.push_back(10 * law_gaussian());
    }

  Db* data = Db::createFromSamples(nech,ELoadBy::COLUMN,tab);
  data->setNameByAttribute(1,"xcoor1");
  data->setNameByAttribute(2,"xcoor2");
  data->setNameByAttribute(3,"var");
  data->setLocatorByAttribute(1,ELoc::X,0);
  data->setLocatorByAttribute(2,ELoc::X,1);
  data->setLocatorByAttribute(3,ELoc::Z);
  data->display();

  // Create the Model
  CovContext ctxt(nvar); // use default space
  Model model(ctxt);
  CovLMC covs(ctxt.getSpace());
  CovAniso cova(ECov::SPHERICAL, 5., 0., 45, ctxt);
  covs.addCov(&cova);
  model.setCovList(&covs);
  model.display();

  // Creating a Neighborhood
  NeighUnique* neighU = NeighUnique::create(ndim,false);
  neighU->display();

  // Launch kriging
  kriging(data, grid, &model, neighU);

  delete data;
  delete grid;
  delete neighU;

  return (0);
}
