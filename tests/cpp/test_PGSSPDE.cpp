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
#include "Model/Model.hpp"
#include "API/SPDE.hpp"
#include "API/PGSSPDE.hpp"

/****************************************************************************/
/*!
** Main Program for testing the sparse matrix algebra
**
*****************************************************************************/
int main(int argc, char *argv[])

{
  ASerializable::setContainerName(true);
  ASerializable::setPrefixName("PGSSPDE-");
  int error = 0;
  int ndim = 2;
  ASpaceObject::defineDefaultSpace(SPACE_RN, ndim);
  CovContext ctxt(1,2,1.);

  // Prepare the Discrete process with Discretized Option
  set_test_discrete(false);

  // Creating a Point Data base in the 1x1 square with 'nech' samples
  int nech = 10;
  Db db(nech,{0.,0.},{1.,1.});
  db.display(FLAG_STATS);

  auto nx={ 101,101 };
  Db workingDbc(nx);

  Db dbprop= Db({100,100},{0.01,0.01},{0.,0.});

  VectorDouble props({0.2, 0.5, 0.3});
  int nfac = props.size();
  VectorString names = generateMultipleNames("Props",nfac);
  for (int ifac = 0; ifac < nfac; ifac++)
    dbprop.addFields(1,props[ifac],names[ifac]);
  dbprop.setLocator(names,LOC_P);

  // Creating the Model(s) of the Underlying GRF(s)
  Model model1(ctxt);
  double range1 = 20;
  CovAniso cova1(ECov::BESSEL_K,range1,1.,1.,ctxt);
  model1.addCova(&cova1);
  model1.display();
  model1.serialize("truemodel1.ascii");

  Model model2(ctxt);
  double range2 = 40;
  CovAniso cova2(ECov::BESSEL_K,range2,2.,1.,ctxt);
  model2.addCova(&cova2);
  model2.display();
  model2.serialize("truemodel2.ascii");

  std::vector<Model*> models;
  models.push_back(&model1);
  models.push_back(&model2);

  // Creating the Rule
  Rule rule({"S","T","F1","F2","F3"});
  RuleProp ruleProp = RuleProp(&rule, props);

  auto ndata = 100;
  Db dat = Db(ndata, { 0., 0. }, { 100., 100. });
  VectorDouble z = ut_vector_simulate_gaussian(ndata);
  dat.addFields(z,"variable",LOC_Z);

  PGSSPDE sCond(models,workingDbc,ruleProp,&dat);
  PGSSPDE sNonCond(models,workingDbc,ruleProp);

  PGSSPDE spgs = sNonCond;
  spgs.simulate();
  spgs.query(&workingDbc);
  workingDbc.display();
  workingDbc.serialize("pgs.ascii");

  return(error);
}
