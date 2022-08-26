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
#include "Db/Db.hpp"
#include "Basic/File.hpp"
#include "Db/DbStringFormat.hpp"
#include "Model/Model.hpp"
#include "Covariances/CovAniso.hpp"
#include "Covariances/CovLMC.hpp"
#include "Simulation/CalcSimuTurningBands.hpp"
#include "Matrix/MatrixSquareGeneral.hpp"
#include "Geometry/Geometry.hpp"

/**
 * This file is meant to perform any test that needs to be coded for a quick trial
 */
int main(int /*argc*/, char */*argv*/[])

{
  // Standard output redirection to file
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str());

  double dzoverdx = 0.5;
  double dzoverdy = 0.2;
  message("Gradients = %lf %lf\n", dzoverdx, dzoverdy);

  MatrixSquareGeneral rotmat = util_gradXYToRotmat(dzoverdx, dzoverdy);
  rotmat.display();

  VectorDouble angles = util_rotmatToEuler(rotmat);
  ut_vector_display("angles",angles);

  MatrixSquareGeneral rotmat2 = util_EulerToRotmat(angles);
  rotmat2.display();

  double diff = rotmat.compare(rotmat2);
  message("\nDifference between two rotation matrices = %lg\n",diff);

  return (0);
}
