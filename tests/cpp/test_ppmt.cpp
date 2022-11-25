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
#include "Basic/MathFunc.hpp"
#include "Basic/Law.hpp"
#include "Basic/File.hpp"
#include "csparse_f.h"

/****************************************************************************/
/*!
** Main Program for testing the sparse matrix algebra
**
*****************************************************************************/
int main(int /*argc*/, char */*argv*/[])
{
  // Standard output redirection to file
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
    StdoutRedirect sr(sfn.str());

  ASerializable::setContainerName(true);
  ASerializable::setPrefixName("PPMT-");

  defineDefaultSpace(ESpaceType::RN, 2);
  
  Db* data = Db::createFromNF("Data.ascii");
  int nech = data->getSampleNumber();

  VectorDouble XX = data->getColumns({"Y1","Y2"});
  MatrixRectangular X(nech, 2);
  X.setValues(XX.data(),true);

  PPMT ppmt;
  MatrixSquareSymmetric E = ppmt.sphering(X);

  return(0);
}
