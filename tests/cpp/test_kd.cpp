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
#include "Space/Space.hpp"
#include "Space/ASpaceObject.hpp"
#include "Db/Db.hpp"
#include "Db/DbStringFormat.hpp"
#include "Model/Model.hpp"
#include "Model/Constraints.hpp"
#include "Model/Option_VarioFit.hpp"
#include "Covariances/CovAniso.hpp"
#include "Covariances/CovLMC.hpp"
#include "Covariances/ECov.hpp"
#include "Drifts/Drift1.hpp"
#include "Drifts/DriftX.hpp"
#include "Drifts/DriftY.hpp"
#include "Basic/Law.hpp"
#include "Basic/File.hpp"
#include "Basic/OptDbg.hpp"
#include "Basic/OptCustom.hpp"
#include "Neigh/ANeighParam.hpp"
#include "Neigh/NeighMoving.hpp"
#include "Anamorphosis/AnamHermite.hpp"
#include "Anamorphosis/AnamContinuous.hpp"
#include "Variogram/VarioParam.hpp"
#include "Variogram/Vario.hpp"
#include "Simulation/CalcSimuTurningBands.hpp"

#include <math.h>

/****************************************************************************/
/*!
 ** Main Program
 **
 *****************************************************************************/
int main(int /*argc*/, char */*argv*/[])
{
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str());

  // Global parameters
  int ndim = 2;
  law_set_random_seed(32131);

  ASpaceObject::defineDefaultSpace(SPACE_RN, ndim);
  DbStringFormat dbfmt(FLAG_STATS);

  // Generate initial grid
  DbGrid* grid = DbGrid::create({100,100}, {0.01,0.01});

  // Create grid of (single) Panel
  double dx_P = 0.250;
  double x0_P = 0.375;
  DbGrid* panel = DbGrid::create({1,1}, {dx_P, dx_P}, {x0_P, x0_P});
  panel->display();

  // Discretization with blocs
  int nx_B = 5;
  double x0_B = x0_P + dx_P / nx_B / 2.;
  double dx_B = dx_P / nx_B;
  DbGrid* blocs = DbGrid::create({nx_B,nx_B}, {dx_B,dx_B}, {x0_B,x0_B});
  blocs->display();

  // Simulation of the Data
  Model* model_init = Model::createFromParam(ECov::EXPONENTIAL, 0.1, 1.);
  (void) simtub(nullptr, grid, model_init);
  grid->setName("Simu", "Y");
  grid->display();

  // Non-linear transform
  double m_Z = 1.5;
  double s_Z = 0.5;
  VectorDouble Zval = grid->getColumn("Y");
  for (int i = 0; i < (int) Zval.size(); i++)
    Zval[i] = m_Z * exp(s_Z * Zval[i] - 0.5 * s_Z * s_Z);
  grid->addColumns(Zval, "Z");
  grid->display(&dbfmt);

  // Data extraction
  int np = 500;
  Db* data = Db::createSamplingDb(grid, 0., np, {"x1","x2","Y","Z"});
  data->setLocator("Z", ELoc::Z);
  data->display(&dbfmt);

  // Gaussian Anamorphosis with 10 coefficients
  AnamHermite* anam = AnamHermite::create(20);
  anam->fit(data);
  anam->display();

  // Selectivity
  Selectivity* selectivity =
      Selectivity::createByCodes( { ESelectivity::Q, ESelectivity::T },
                                  { 0., 0.5 }, true, true);

  // Global experimental selectivity
  data->display();
  selectivity->eval(data).display();

  // Selectivity in the model
  selectivity->eval(anam).display();

  // Transform Data into Gaussian
  (void) anam->RawToGaussian(data);
  data->setName("Y.Z","Gauss.Z");
  data->display();

  // Calculate the variogram
  VarioParam* varioparam = VarioParam::createOmniDirection(ndim, 10, 0.025);
  Vario* vario = Vario::computeFromDb(varioparam, data);
  vario->display();

  // Fitting the Model on the Gaussian transformed variable
  Model* model = new Model(1, ndim);
  Constraints constraints = Constraints();
  constraints.setConstantSillValue(1.);
  (void) model->fit(vario, {ECov::EXPONENTIAL, ECov::EXPONENTIAL}, false,
                     Option_AutoFit(), constraints);
  model->display();

  // Attach the Anamorphosis
  model->setAnam(anam);

  // Computing the Point factors
  int nfactor = 3;
  (void) calculateHermiteFactors(data, nfactor);
  data->display();

  // Creating a Moving Neighborhood
  int nmini = 5;
  int nmaxi = 5;
  double radius = 1.;
  NeighMoving* neigh = NeighMoving::create(ndim, false, nmaxi, radius, nmini);
  neigh->display();

  // Setting the trace
  OptDbg::setReference(1);

  // ====================== Point Disjunctive Kriging =====================

  // Simple Point Kriging over the blocks
  (void) dk(data, blocs, model, neigh, EKrigOpt::PONCTUAL, VectorInt(),
            true, true, NamingConvention("DK_Pts"));
  blocs->display();

  // Simple Block Kriging over the blocks

  VectorInt ndisc_B = {5,5};
  (void) dk(data, blocs, model, neigh, EKrigOpt::BLOCK, ndisc_B,
            true, true, NamingConvention("DK_Blk"));
  blocs->display();

  // Simple Block Kriging over the panel(s)
  VectorInt ndisc_P = { 10,10};
  (void) dk(data, panel, model, neigh, EKrigOpt::BLOCK, ndisc_P,
            true, true, NamingConvention("DK_Blk"));
  panel->display();

  // ====================== Block Disjunctive Kriging (DGM-1) =====================

  // Calculate the change of support coefficient
  model->setAnamIClass(0); // Z variable
  double cvv = model->evalCvv(blocs->getDXs(), ndisc_B, blocs->getAngles());
  double r1 = sqrt(anam->invertVariance(cvv));
  message("Change of Support coefficient (DGM-1)= %6.3lf\n", r1);

  // Update the Model with Block anamorphosis
  AnamHermite* anam_b1 = anam->clone();
  anam_b1->setRCoef(r1);

  // Regularization of the point model by the block support
  Vario* vario_b1_Z = Vario::createRegularizeFromModel(model, varioparam, blocs->getDXs(),
                                                 ndisc_B, blocs->getAngles());
  Vario* vario_b1_Y = Vario::createTransformZToY(vario_b1_Z, anam, cvv);

  // Fitting the regularized model on the point Gaussian variable
  Model* model_b1_Y = new Model(1, ndim);
  constraints.setConstantSillValue(1);
  (void) model_b1_Y->fit(vario_b1_Y, { ECov::CUBIC, ECov::EXPONENTIAL }, false,
                      Option_AutoFit(), constraints);
  model_b1_Y->setAnam(anam_b1);
  model_b1_Y->display();

  // Simple Point Kriging over the blocs(s) with Model with Change of Support
  (void) dk(data, blocs, model_b1_Y, neigh, EKrigOpt::PONCTUAL, VectorInt(), true,
            true, NamingConvention("DK_DGM1"));
  blocs->display();

  // Simple Point Kriging over the panel(s) with Model with Change of Support
  (void) dk(data, panel, model_b1_Y, neigh, EKrigOpt::BLOCK, { nx_B, nx_B }, true,
            true, NamingConvention("DK_DGM1"));
  panel->display();

  // ====================== Block Disjunctive Kriging (DGM-2) =====================

  // Calculate the change of support coefficient
  model->setAnamIClass(1); // Y Variable
  double r2 = sqrt(model->evalCvv(blocs->getDXs(), ndisc_B, blocs->getAngles()));
  message("Change of Support coefficient (DGM2)= %6.3lf\n",r2);

  // Regularization of the point model by the block support
  Vario* vario_b2_Y = Vario::createRegularizeFromModel(model, varioparam, blocs->getDXs(),
                                                 ndisc_B, blocs->getAngles());

  // Fitting the regularized model on the point Gaussian variable
  Model* model_b2_Y = new Model(1, ndim);
  constraints.setConstantSillValue(r2 * r2);
  (void) model_b2_Y->fit(vario_b2_Y, { ECov::CUBIC, ECov::EXPONENTIAL }, false,
                     Option_AutoFit(), constraints);
  model_b2_Y->display();

  // Normalization of the block model to a total sill equal to 1.0
  model_b2_Y->normalize(1.0);

  // Update the Model with Block anamorphosis
  AnamHermite* anam_b2 = anam->clone();
  anam_b2->setRCoef(r2);
  model_b2_Y->setAnam(anam_b2);

  // Simple Point Kriging over the blocs(s) with Model with Change of Support
  (void) dk(data, blocs, model_b2_Y, neigh, EKrigOpt::PONCTUAL, VectorInt(), true,
            true, NamingConvention("DK_DGM2"));
  blocs->display();

  // Simple Point Kriging over the panel(s) with Model with Change of Support
  (void) dk(data, panel, model_b2_Y, neigh, EKrigOpt::BLOCK, {nx_B, nx_B},
            true, true, NamingConvention("DK_DGM2"));
  panel->display();

  // ====================== Selectivity Function ==================================

  anamFactor2Selectivity(blocs, anam, selectivity,
                         blocs->getNames("DK_Pts*estim"),blocs->getNames("DK_Pts*stdev"));
  blocs->display();

  // ====================== Free pointers ==================================

  if (data       != nullptr) delete data;
  if (grid       != nullptr) delete grid;
  if (model_init != nullptr) delete model_init;
  if (model      != nullptr) delete model;
  if (model_b1_Y != nullptr) delete model_b1_Y;
  if (model_b2_Y != nullptr) delete model_b2_Y;
  if (panel      != nullptr) delete panel;
  if (blocs      != nullptr) delete blocs;
  if (anam       != nullptr) delete anam;
  if (anam_b1    != nullptr) delete anam_b1;
  if (anam_b2    != nullptr) delete anam_b2;
  if (varioparam != nullptr) delete varioparam;
  if (vario      != nullptr) delete vario;
  if (vario_b1_Z != nullptr) delete vario_b1_Z;
  if (vario_b1_Y != nullptr) delete vario_b1_Y;
  if (vario_b2_Y != nullptr) delete vario_b2_Y;

  if (neigh      != nullptr) delete neigh;

  return (0);
}
