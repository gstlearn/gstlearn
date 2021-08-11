#include "Basic/AException.hpp"
#include "Basic/Vector.hpp"
#include "Covariances/CovAniso.hpp"
#include "Db/Db.hpp"
#include "geoslib_e.h"

#include "LinearOp/PrecisionOpMultiConditional.hpp"
#include "LinearOp/ProjMatrix.hpp"
#include "API/SPDE.hpp"
#include "Model/Model.hpp"
#include "Model/NoStatArray.hpp"
#include "Mesh/AMesh.hpp"
#include "Mesh/MeshETurbo.hpp"
#include "Mesh/MeshFactory.hpp"

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
  auto pygst = std::string(std::getenv("PYGSTLEARN_DIR"));
  int seed = 10355;
  std::normal_distribution<double> d{0,1};
  std::uniform_real_distribution <double> u{0,1};
  std::mt19937 gen{seed};

  ///////////////////////
  // Création de la db //

  auto nx={ 101,101 };
  Db workingDbc(nx);
  VectorDouble angle;
  // Generating angles
  for(auto &e : workingDbc.getCoordinates())
  {
    angle.push_back(spirale(e));
  }

  workingDbc.addFields(angle,"angle",LOC_NOSTAT);

  ///////////////////////
  // Creating the Model
  Model model = Model(&workingDbc);
  CovAniso cova = CovAniso(COV_BESSEL_K,model.getContext());
  cova.setRanges({10,45});
  model.addCova(&cova);

  //////////////////////
  //Creating the Mesh
  MeshETurbo mesh(workingDbc);

  /////////////////////////////////////////////////////
  // Creating the Precision Operator for simulation

  NoStatArray NoStat({"A"},&workingDbc);
  model.addNoStat(&NoStat);
  SPDE spde(model,workingDbc);
  std::cout<<"end creation "<<std::endl;

  ShiftOpCs S(&mesh, &model, &workingDbc);
  PrecisionOp Qsimu(&S, &cova, POPT_MINUSHALF);

//
//  ///////////////////////////////////////////////////
//  // Simulation (Chebyshev)
//
  VectorDouble tab;
  VectorDouble resultSimu;

  for (int iech = 0; iech < mesh.getNApices(); iech++)
  {
    tab.push_back(d(gen));
  }

  resultSimu.resize(tab.size());
  Qsimu.eval(tab,resultSimu);
  workingDbc.addFields(resultSimu,"Simu",LOC_Z);


//  // Création des données (peut-être qu'on pourrait utiliser un équivalent de db.grid.init)
  int ndata = 1000;


  VectorDouble coordsX,coordsY;

  for (int iech = 0; iech < ndata; iech++)
  {
      coordsX.push_back(99. * u(gen));
      coordsY.push_back(99. * u(gen));
  }


  Db dat;
  dat.addFields(coordsX,"X");
  dat.addFields(coordsY,"Y");
  VectorString vct={"X","Y"};
  dat.setLocator(vct,LOC_X);

  // Simulating Data points

   ProjMatrix B(&dat,&mesh);
   VectorDouble datval(ndata);
   B.mesh2point(resultSimu,datval);
   dat.addFields(datval,"Simu",LOC_Z);

   // Kriging
   double nug = 0.1;
   VectorDouble rhs(S.getSize());
   B.point2mesh(dat.getField("Simu"),rhs);
   for(auto &e:rhs)
   {
     e/=nug;
   }

   PrecisionOp Qkriging(&S, &cova,POPT_ONE);
   PrecisionOpMultiConditional A;
   A.push_back(&Qkriging,&B);
   A.setVarianceData(0.01);


   VectorVectorDouble Rhs,resultvc;
   VectorDouble vc(S.getSize());

   resultvc.push_back(vc);
   Rhs.push_back(VectorDouble(rhs));

   A.evalInverse(Rhs,resultvc);
   workingDbc.addFields(resultvc[0],"Kriging");
   workingDbc.serialize(pygst + "spde.ascii");
//
  return 0;
}
