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
#include "geoslib_old_f.h"

#include "Anamorphosis/PPMT.hpp"
#include "Matrix/MatrixRectangular.hpp"
#include "Matrix/MatrixSquareSymmetric.hpp"
#include "Db/Db.hpp"
#include "Db/DbStringFormat.hpp"
#include "Basic/VectorNumT.hpp"
#include "Basic/MathFunc.hpp"
#include "Basic/Law.hpp"
#include "Basic/File.hpp"
#include "Basic/OptCst.hpp"
#include "csparse_f.h"

/****************************************************************************/
/*!
** Main Program for testing the sparse matrix algebra
**
*****************************************************************************/
int main(int /*argc*/, char */*argv*/[])
{
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str());

  ASerializable::setContainerName(true);
  ASerializable::setPrefixName("PPMT-");

  defineDefaultSpace(ESpaceType::RN, 2);
  OptCst::define(ECst::NTROW, 15);
  
  Db* data = Db::createFromNF("Data.ascii");

  // Creating PPMT model
  int ndir = 10;
  int niter = 100;
  PPMT ppmt(ndir, niter);

  // Fitting the PPMT model
  ppmt.rawToGaussian(data, {"Y1","Y2"}, false, true, NamingConvention("U"));
  data->display();

  return 0;
}
