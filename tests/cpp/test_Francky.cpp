#include "Basic/AException.hpp"
#include "Basic/Vector.hpp"
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
  std::normal_distribution<double> d{0,1};
  std::uniform_real_distribution <double> u{0,1};
  std::mt19937 gen{seed};

  ///////////////////////
  // Création de la db //

  auto nx={ 101,101 };
  Db workingDbc(nx);
  VectorDouble angle;
  // Génération des angles
  for(auto &e : workingDbc.getCoordinates())
  {
    angle.push_back(spirale(e));
  }

  workingDbc.addFields(angle,"angle",LOC_NOSTAT);

  ///////////////////////
  // Création du modèle
  Model model = Model(&workingDbc);
  CovAniso cova = CovAniso(COV_BESSEL_K,model.getContext());
  cova.setRanges({10,45});
  model.addCova(&cova);


  NoStatArray NoStat({"A"},&workingDbc);
  model.addNoStat(&NoStat);


  // Création des données (peut-être qu'on pourrait utiliser un équivalent de db.grid.init)
  int ndata = 10;


  VectorDouble coordsX,coordsY,Z;

  for (int iech = 0; iech < ndata; iech++)
  {
      coordsX.push_back(99. * u(gen));
      coordsY.push_back(99. * u(gen));
      Z.push_back(d(gen));
  }

//
  Db dat;
  dat.addFields(coordsX,"X");
  dat.addFields(coordsY,"Y");
  VectorString vct={"X","Y"};
  dat.setLocator(vct,LOC_X);
  dat.addFields(Z,"Z");
  dat.setLocator("Z",LOC_Z);

  Neigh neigh(2);
  kriging(&dat,&workingDbc,&model,&neigh);
//


  return 0;
}
