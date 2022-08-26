#include "geoslib_old_f.h"

#include "API/SPDE.hpp"
#include "Model/Model.hpp"
#include "Model/NoStatArray.hpp"
#include "Basic/AException.hpp"
#include "Basic/Vector.hpp"
#include "Basic/Law.hpp"
#include "Basic/File.hpp"
#include "Space/SpaceRN.hpp"
#include "Covariances/CovAniso.hpp"
#include "Covariances/CovFactory.hpp"
#include "Covariances/CovLMC.hpp"
#include "Covariances/ECov.hpp"
#include "Db/Db.hpp"
#include "Db/DbGrid.hpp"

#include <algorithm>
#include <math.h>
#include <iostream>
#include <numeric>
#include <random>
#include <string>
#include <vector>
#include <random>

#ifndef M_PI // Under Visual, this is not defined
  const double M_PI = 3.14159265358979323846;
#endif

double fa(double x,double y,double a,double b)
{
    return a*x + b*y;
}

double spirale(VectorDouble pos)
{

  auto  a=0.;
  auto  b=-1.4;
  auto  c=1.;
  auto  d=1.;

  auto x = pos[0];
  auto y = pos[1];
  auto u1=fa(x-50.,y-50.,a,b);
  auto u2=fa(x-50.,y-50.,c,d);
  auto norm = sqrt(u1*u1+u2*u2);
  if(norm >0)
  {
    return acos(u2/norm)/M_PI*180.* ((u1>=0)?1.:-1.);
  }
  else
  {
    return 0.;
  }
}

/****************************************************************************/
/*!
 ** Main Program
 **
 *****************************************************************************/
int main(int /*argc*/, char */*argv*/[])

{
  // Standard output redirection to file
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str());

  // R_model_characteristics
  int    flag_range,flag_param,min_order,max_ndim,ndim,flag_aniso;
  int    flag_int_1d,flag_int_2d,flag_rotation;
  double scale,parmax;
  char   cova_name[STRING_LENGTH];
  std::string cov_name;

  /* Initializations */

  cov_name = "Nugget Effect";
  ndim     = 2;

  /* Identify the covariance */

  SpaceRN space = SpaceRN(ndim);
  CovContext ctxt = CovContext(1, ndim);
  ECov type = CovFactory::identifyCovariance(cov_name, ctxt);
  if (type == ECov::UNKNOWN)
    return 1;

  /* Asking for complementary information */

  model_cova_characteristics(type, cova_name,
                             &flag_range,&flag_param,&min_order,&max_ndim,
                             &flag_int_1d,&flag_int_2d,
                             &flag_aniso,&flag_rotation,&scale,&parmax);

  int seed = 432432;
  law_set_random_seed(seed);
  std::mt19937 gen;
  gen.seed(seed);

  std::normal_distribution<double> d{0,1};

  return 0;
  ///////////////////////
  // Creating the Grid Db
  auto nx = { 101, 101 };
  DbGrid* workingDbc = DbGrid::create(nx);
  VectorDouble angle;
  for(auto &e : workingDbc->getAllCoordinates())
  {
    angle.push_back(spirale(e));
  }
  workingDbc->addColumns(angle,"angle",ELoc::NOSTAT);

  ///////////////////////
  // Creating the Model
  Model* model = Model::createFromDb(workingDbc);
  CovLMC covs(ctxt.getSpace());
  CovAniso cova = CovAniso(ECov::BESSEL_K,model->getContext());
  cova.setRanges({4,45});
  covs.addCov(&cova);
  model->setCovList(&covs);

  // Non-stationary part
  NoStatArray nostat({"A"}, workingDbc);
  model->addNoStat(&nostat);

  ///////////////////////////////////////////////////
  // Simulation (Chebyshev)

  // Creating the Data
  int ndata = 1000;
  Db* dat = Db::createFromBox(ndata, {0.,0.}, {100.,100.}, 4324);
  VectorDouble tab;

  for (int iech = 0; iech < dat->getSampleNumber(); iech++)
  {
    tab.push_back(law_gaussian());
  }

  dat->addColumns(tab, "Simu", ELoc::Z);

  SPDE spde(model,workingDbc,dat,ESPDECalcMode::KRIGING);
  spde.compute();
  spde.query(workingDbc);

  delete dat;
  delete workingDbc;
  delete model;
  return 0;
}
