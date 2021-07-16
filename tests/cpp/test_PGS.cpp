/******************************************************************************/
/* COPYRIGHT ARMINES, ALL RIGHTS RESERVED                                     */
/*                                                                            */
/* THE CONTENT OF THIS WORK CONTAINS CONFIDENTIAL AND PROPRIETARY             */
/* INFORMATION OF ARMINES. ANY DUPLICATION, MODIFICATION,                     */
/* DISTRIBUTION, OR DISCLOSURE IN ANY FORM, IN WHOLE, OR IN PART, IS STRICTLY */
/* PROHIBITED WITHOUT THE PRIOR EXPRESS WRITTEN PERMISSION OF ARMINES         */
/*                                                                            */
/* TAG_SOURCE_CG                                                              */
/*                                                                            */
/* This file is meant to demonstrate the process of using PGS                 */
/*                                                                            */
/******************************************************************************/
#include "geoslib_f.h"
#include "Neigh/Neigh.hpp"
#include "Model/Model.hpp"
#include "Variogram/Vario.hpp"

/****************************************************************************/
/*!
** Main Program for testing the sparse matrix algebra
**
*****************************************************************************/
int main(int argc, char *argv[])

{
  auto pygst  = std::string(std::getenv("PYGSTLEARN_DIR"));
  int error = 0;
  int ndim = 2;
  CovContext ctxt(1,2,1.);

  // Creating a Point Data base in the 1x1 square with 'nech' samples
  int nech = 1000;
  Db db(nech,VectorDouble(2,0.),VectorDouble(2,1.));
  db.display(FLAG_STATS);

  // Creating the Model(s) of the Underlying GRF(s)
  Model model1(ctxt);
  double range1 = 0.2;
  CovAniso cova1(COV_BESSEL_K,range1,1.,1.,ctxt);
  model1.addCova(&cova1);
  model1.display();
  model1.serialize(pygst+ "data/truemodel1.ascii");

  Model model2(ctxt);
  double range2 = 0.3;
  CovAniso cova2(COV_EXPONENTIAL,range2,1.,1.,ctxt);
  model2.addCova(&cova2);
  model2.display();
  model2.serialize(pygst+ "data/truemodel2.ascii");

  // Creating the Neighborhood
  Neigh neigh = Neigh(ndim);
  neigh.display();

  // Pour voir
  //error = simtub(nullptr,&db,&model,&neigh,4);
  //db.display(1);
  //db.setLocator("Simu*",LOC_F);
  //db.display(1);

  // Creating the Rule
  Rule rule({"S","T","F1","F2","F3"});
  rule.display();

  // Perform a non-conditional simulation on the Db
  VectorDouble props({0.2, 0.5, 0.3});
  error = simpgs(nullptr,&db,nullptr,&rule,&model1,&model2,&neigh,props);
  db.setLocator(db.getLastName(),LOC_Z);

  // Adding constant proportions as Vectors in the Db
  int nfac = props.size();
  VectorString names = generateMultipleNames("Props",nfac);
  for (int ifac = 0; ifac < nfac; ifac++)
    db.addFields(1,props[ifac],names[ifac]);
  db.setLocator(names,LOC_P);
  db.display(FLAG_STATS);

  // Determination of the variogram of the Underlying GRF
  Vario cov = Vario();
  int nlag = 19;
  Dir dir = Dir(ndim, nlag, 0.5 / nlag);
  cov.addDirs(dir);
  error = variogram_pgs(&db,&cov,&rule,props);
  Vario vario(cov,VectorInt(),VectorInt(),true);
  vario.display(1);

  // Fitting the experimental variogram o Underlying GRF (with constraint that total sill is 1)
  Model modelPGS(ctxt);
  Option_AutoFit option = Option_AutoFit();
  option.setConstantSillValue(1.);

  std::vector<ENUM_COVS> covs {COV_BESSEL_K, COV_EXPONENTIAL};
  modelPGS.fit(&vario,covs,true,option);
  modelPGS.display();

  vario.serialize(pygst+ "data/variopgs.ascii");
  modelPGS.serialize(pygst+ "data/modelfitpgs.ascii");

  // Compute the experimental variograms of the indicators

  Dir dir2 = Dir(ndim, nlag, 0.5 / nlag);
  Vario vario2 = Vario();
  vario2.addDirs(dir2);
  error = vario2.computeIndic(&db);
  vario2.serialize(pygst+ "data/varioindic.ascii");

  error = model_pgs(&db, &vario2, &rule, &modelPGS, nullptr, props);
  vario2.serialize(pygst+ "data/modelpgs.ascii");
//
////  vario.display(1);


  return(error);
}
