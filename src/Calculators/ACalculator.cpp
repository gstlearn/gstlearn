/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://gstlearn.org                                              */
/* License: BSD 3 clauses                                                     */
/*                                                                            */
/******************************************************************************/
#include "Calculators/ACalculator.hpp"
#include "Db/Db.hpp"
#include "Model/Model.hpp"
#include "Basic/AException.hpp"
#include "Basic/AStringable.hpp"

ACalculator::ACalculator()
{
}

ACalculator::~ACalculator()
{
}

/**
 * Run the calculator
 *
 * \return false if error, true otherwise
 */
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
