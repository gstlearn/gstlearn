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
  int error = 0;
  int ndim = 2;
  CovContext ctxt(1,2,1.);

  // Creating a Point Data base in the 1x1 square with 'nech' samples
  int nech = 100;
  Db db(nech,VectorDouble(2,0.),VectorDouble(2,1.));
  db.display(FLAG_STATS);

  // Creating the Model of the Underlying GRF
  Model model(ctxt);
  double range = 0.2;
  CovAniso cova(COV_CUBIC,range,1.,1.,ctxt);
  model.addCova(&cova);
  model.display();

  // Creating the Neighborhood
  Neigh neigh = Neigh(ndim);
  neigh.display();

  // Pour voir
  //error = simtub(nullptr,&db,&model,&neigh,4);
  //db.display(1);
  //db.setLocator("Simu*",LOC_F);
  //db.display(1);

  // Creating the Rule
  Rule rule({"S","S","F1","F2","F3"});
  rule.display();

  // Perform a non-conditional simulation on the Db
  VectorDouble props({0.2, 0.5, 0.3});
  error = simpgs(nullptr,&db,nullptr,&rule,&model,nullptr,&neigh,props);
  db.setLocator(db.getLastName(),LOC_Z);

  // Adding constant proportions as Vectors in the Db
  int nfac = props.size();
  VectorString names = generateMultipleNames("Props",nfac);
  for (int ifac = 0; ifac < nfac; ifac++)
  {
    db.addFields(1,props[ifac]);
    db.setName(db.getLastName(),names[ifac]);
  }
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

  std::vector<ENUM_COVS> covs {COV_CUBIC};
  modelPGS.fit(&vario,covs,true,option);
  modelPGS.display();

  // Deriving the experimental variograms of the indicators

  error = model_pgs(&db, &vario, &rule, &modelPGS, nullptr, props);
  vario.display(1);

  return(error);
}
