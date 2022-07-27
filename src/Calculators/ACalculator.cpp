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
#include "Calculators/ACalculator.hpp"
#include "Db/Db.hpp"
#include "Model/Model.hpp"
#include "Neigh/ANeighParam.hpp"
#include "Basic/AException.hpp"
#include "Basic/AStringable.hpp"

ACalculator::ACalculator()
{
}

ACalculator::~ACalculator()
{
}

bool ACalculator::run()
{
  try
  {
    if (! _check())
      my_throw("Check has failed. Calculation aborted");
    if (! _preprocess())
      my_throw("Pre-processing has failed. Calculation aborted");
    if (! _run())
      my_throw("Run has failed. Calculation aborted");
    if (! _postprocess())
      my_throw("Post-processing has failed.");
  }
  catch(const AException& e)
  {
    messerr("Calculator has failed: %s",e.what());
    _rollback();
    return false;
  }
  catch(const std::exception& e)
  {
    messerr("Calculator has failed: %s",e.what());
    _rollback();
    return false;
  }
  return true;
}
