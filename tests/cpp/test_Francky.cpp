#include "Basic/AException.hpp"
#include "Basic/Vector.hpp"
#include "Basic/Law.hpp"
#include "Covariances/CovAniso.hpp"
#include "Db/Db.hpp"
#include "geoslib_e.h"

#include "geoslib_f.h"

#include "Model/Model.hpp"
#include "Model/NoStatArray.hpp"


#include "MatrixC/MatrixCRectangular.hpp"

#include <algorithm>
#include <math.h>
#include <iostream>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#define __USE_MATH_DEFINES
#include <cmath>

double fa(double x,double y,double a,double b)
{
    return a*x + b*y;
}

double spirale(std::vector<double> pos)
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
int main(int argc, char *argv[])

{
  int seed = 10355;
  law_set_random_seed(seed);

  // Creating the 2-D Db
  auto nx = { 101, 101 };
  Db workingDbc(nx);

  // Generating a vector of angles (spirale)
  VectorDouble angle;
  for(auto &e : workingDbc.getCoordinates())
    angle.push_back(spirale(e));
  workingDbc.addFields(angle,"angle",LOC_NOSTAT);

  // Creating the Non-stationary Model
  Model model = Model(&workingDbc);
  CovAniso cova = CovAniso(COV_BESSEL_K,model.getContext());
  cova.setRanges({10,45});
  model.addCova(&cova);
  NoStatArray NoStat({"A"},&workingDbc);
  model.addNoStat(&NoStat);

  // Creating the 2-D Data Db with a Normal Variable
  int ndata    = 10;
  auto coormin = {0., 0.};
  auto coormax = {100., 100.};
  Db dat = Db(ndata, coormin, coormax);
  VectorDouble Z = ut_vector_simulate_gaussian(ndata);
  dat.addFields(Z, "Z");
  dat.setLocator("Z", LOC_Z);

  // Creating the Neighborhood (Unique)
  Neigh neigh(2);

  // Testing Kriging
  kriging(&dat,&workingDbc,&model,&neigh);
  dat.display(1);
  workingDbc.display(1);

  // Testing the storage of the non-stationary parameters within Db
  db_model_nostat(&dat,&model,0);
  dat.display(1);

  return 0;
}
