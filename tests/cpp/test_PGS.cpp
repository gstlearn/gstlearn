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
  auto pygst = std::string(std::getenv("PYGSTLEARN_DIR"));
  int error = 0;
  int ndim = 2;
  CovContext ctxt(1,2,1.);

  // Prepare the Discrete process with Discretized Option
  set_test_discrete(0);

  // Creating a Point Data base in the 1x1 square with 'nech' samples
  int nech = 1000;
  Db db(nech,{0.,0.},{1.,1.});
  db.display(FLAG_STATS);

  // Creating the Model(s) of the Underlying GRF(s)
  Model model1(ctxt);
  double range1 = 0.2;
  CovAniso cova1(COV_BESSEL_K,range1,1.,1.,ctxt);
  model1.addCova(&cova1);
  model1.display();
  model1.serialize(pygst+ "truemodel1.ascii");

  Model model2(ctxt);
  double range2 = 0.3;
  CovAniso cova2(COV_EXPONENTIAL,range2,1.,1.,ctxt);
  model2.addCova(&cova2);
  model2.display();
  model2.serialize(pygst+ "truemodel2.ascii");

  // Creating the Neighborhood
  Neigh neigh = Neigh(ndim);
  neigh.display();

  // Creating the Rule
  Rule rule({"S","T","F1","F2","F3"});
  rule.display();
  rule.serialize(pygst+ "truerule.ascii");

  // Prepare proportions
  VectorDouble props({0.2, 0.5, 0.3});
  int nfac = props.size();
  VectorString names = generateMultipleNames("Props",nfac);
  for (int ifac = 0; ifac < nfac; ifac++)
    db.addFields(1,props[ifac],names[ifac]);
  db.setLocator(names,LOC_P);

  // Perform a non-conditional simulation on the Db
  error = simpgs(nullptr,&db,nullptr,&rule,&model1,&model2,&neigh,props);
  db.setLocator(db.getLastName(),LOC_Z);

  // Determination of the variogram of the Underlying GRF
  Vario cov = Vario();
  int nlag = 19;
  Dir dir = Dir(ndim, nlag, 0.5 / nlag);
  cov.addDirs(dir);

  RuleProp ruleprop = RuleProp(&rule, props);
  error = variogram_pgs(&db,&cov,&ruleprop);
  Vario vario1(cov,VectorInt(1,0),VectorInt(),true);
  Vario vario2(cov,VectorInt(1,1),VectorInt(),true);
  vario1.display(1);
  vario2.display(1);

  // Fitting the experimental variogram o Underlying GRF (with constraint that total sill is 1)
  Model modelPGS1(ctxt);
  Model modelPGS2(ctxt);
  Option_AutoFit option = Option_AutoFit();
  option.setConstantSillValue(1.);

  std::vector<ENUM_COVS> covs {COV_BESSEL_K, COV_EXPONENTIAL};
  modelPGS1.fit(&vario1,covs,true,option);
  modelPGS1.display();

  vario1.serialize(pygst+ "variopgs1.ascii");
  modelPGS1.serialize(pygst+ "modelfitpgs1.ascii");

  modelPGS2.fit(&vario2,covs,true,option);
  modelPGS2.display();

  vario2.serialize(pygst+ "variopgs2.ascii");
  modelPGS2.serialize(pygst+ "modelfitpgs2.ascii");

  // Compute the experimental variograms of the indicators

  Vario varioParam = Vario();
  Dir dir3 = Dir(ndim, nlag, 0.5 / nlag);
  varioParam.addDirs(dir3);
  varioParam.setCalculName("vg");

  RuleProp ruleprop2 = RuleProp((Rule*) NULL, props);
  Rule* ruleFit = rule_auto(&db,&varioParam,&ruleprop2,1);
  ruleprop2.setRule(ruleFit);
  ruleFit->display(1);
  ruleFit->serialize(pygst + "ruleFit.ascii");

  Dir dir2 = Dir(ndim, nlag, 0.5 / nlag);
  Vario varioIndic = Vario();
  varioIndic.addDirs(dir2);
  error = varioIndic.computeIndic(&db);
  varioIndic.serialize(pygst+ "varioindic.ascii");

  error = model_pgs(&db, &varioIndic, &ruleprop2, &modelPGS1, &modelPGS2);
  varioIndic.serialize(pygst+ "modelpgs.ascii");
  varioIndic.display(1);


  return(error);
}
