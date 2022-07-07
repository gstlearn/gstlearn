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
//  StdoutRedirect sr(sfn.str());

  // Global parameters
  int ndim = 2;
  law_set_random_seed(32131);

  ASpaceObject::defineDefaultSpace(SPACE_RN, ndim);
  DbStringFormat dbfmt(FLAG_STATS);

  // Generate initial grid
  DbGrid* grd = DbGrid::create({100,100}, {0.01,0.01});
  Model* model_init = Model::createFromParam(ECov::EXPONENTIAL, 0.1, 1.);
  (void) simtub(nullptr, grd, model_init);
  grd->display();
  grd->setName("Simu", "Y");

  // Non-linear transform
  double m_Z = 1.5;
  double s_Z = 0.5;
  VectorDouble Zval = grd->getColumn("Y");
  for (int i = 0; i < (int) Zval.size(); i++)
    Zval[i] = m_Z * exp(s_Z * Zval[i] - 0.5 * s_Z * s_Z);
  grd->addColumns(Zval, "Z");
  grd->display(&dbfmt);

  // Data extraction
  int np = 500;
  Db* data = Db::createSamplingDb(grd, 0., np, {"x1","x2","Y","Z"});
  data->setLocator("Z", ELoc::Z);
  data->display(&dbfmt);

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

  // ====================== Point Disjunctive Kriging =====================

  // Gaussian Anamorphosis with 10 coefficients
  AnamHermite* anam = AnamHermite::create(20);
  anam->fit(data);
  anam->display();

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
  Option_VarioFit optvar = Option_VarioFit();
  optvar.setAuthAniso(0);
  (void) model->fit(vario, {ECov::EXPONENTIAL, ECov::EXPONENTIAL}, true,
                     Option_AutoFit(), constraints, optvar);
  model->display();
  // Attach the Anamorphosis
  model->addAnam(anam);

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

  // Simple Point kriging over the blocks
  (void) dk(data, blocs, model, neigh, EKrigOpt::PONCTUAL, VectorInt(),
            true, true, NamingConvention("DK_Pts"));
  blocs->display();

  // Simple Block Kriging over the blocks
  (void) dk(data, blocs, model, neigh, EKrigOpt::BLOCK, {5,5},
            true, true, NamingConvention("DK_Blk"));
  blocs->display();

  // Simple Block Kriging over the panel(s)
  (void) dk(data, panel, model, neigh, EKrigOpt::BLOCK, {10,10},
            true, true, NamingConvention("DB_BLK"));

  // ====================== Block Disjunctive Kriging =====================


//  // Modeling the block model (anamorphosis and variogram) using DGM2
//    r <- sqrt(model.cvv(v.mesh = c(1,1) * dx_P / nx_B, model = mod_m, ndisc=10))
//    anam.Blc <- anam.point2block(anam = anam.Pts, coeff = r)
//    // il faut aussi definir le model gaussien des Gaussiennes de blocs \rho(v, v') !
//  // ====================== Free pointers ==================================
  if (data       != nullptr) delete data;
  if (grd        != nullptr) delete grd;
  if (model_init != nullptr) delete model_init;
  if (panel      != nullptr) delete panel;
  if (blocs      != nullptr) delete blocs;
  if (anam       != nullptr) delete anam;
  if (varioparam != nullptr) delete varioparam;
  if (vario      != nullptr) delete vario;
  if (model      != nullptr) delete model;
  if (neigh      != nullptr) delete neigh;

  return (0);
}
