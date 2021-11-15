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
#include "Db/Db.hpp"
#include "Model/Model.hpp"
#include "Covariances/CovAniso.hpp"
#include "Basic/Law.hpp"
#include "Neigh/Neigh.hpp"
#include "geoslib_d.h"
#include "geoslib_f.h"

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
  Db* grid = new Db(nx);
  grid->display(0);

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

  Db* data = new Db(nech,ELoadBy::COLUMN,tab);
  data->setName(1,"xcoor1");
  data->setName(2,"xcoor2");
  data->setName(3,"var");
  data->setLocatorByAttribute(1,ELoc::X,0);
  data->setLocatorByAttribute(2,ELoc::X,1);
  data->setLocatorByAttribute(3,ELoc::Z);
  data->display(0);

  // Create the Model
  CovContext ctx(nvar); // use default space
  Model* model = new Model(ctx);
  CovAniso* cov = new CovAniso(ECov::SPHERICAL, 5., 0., 45, ctx);
  model->addCova(cov);
  model->display();

  // Creating a Neighborhood
  Neigh* neigh = new Neigh(ndim);
  neigh->display();

  // Launch kriging
  kriging(data, grid, model, neigh);

  delete model;
  delete neigh;
  delete data;
  delete grid;

  return (0);
}
