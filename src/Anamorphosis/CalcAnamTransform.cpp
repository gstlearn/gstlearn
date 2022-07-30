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
#include "geoslib_f.h"
#include "geoslib_old_f.h"
#include "Db/DbGrid.hpp"
#include "Db/Db.hpp"
#include "Anamorphosis/CalcAnamTransform.hpp"
#include "Anamorphosis/AnamContinuous.hpp"

#include <math.h>

CalcAnamTransform::CalcAnamTransform(AAnam* anam)
    : ACalcDbVarCreator(),
      _iatt(-1),
      _anam(anam)
{
}

CalcAnamTransform::~CalcAnamTransform()
{
}

bool CalcAnamTransform::_check()
{
  if (! hasDb()) return false;

  if (_anam == nullptr)
  {
    messerr("The argument 'anam'  must be defined");
    return false;
  }
  int number = getDb()->getLocatorNumber(ELoc::Z);
  if (number <= 0)
  {
    messerr("The argument 'db'  must have some variable(s) defined");
    return false;
  }
  AnamContinuous* anamC = dynamic_cast<AnamContinuous*>(_anam);
  if (anamC == nullptr)
  {
    messerr("The argument 'anam'  must be of type AnamContinuous");
    return false;
  }
  return true;
}

bool CalcAnamTransform::_preprocess()
{
  int nvar = _getNVar();
  if (nvar < 0) return false;
  _iatt = getDb()->addColumnsByConstant(nvar);
  if (_iatt < 0) return false;
  return true;
}

bool CalcAnamTransform::_postprocess()
{
  int nvar = _getNVar();
  _renameVariable(ELoc::Z, nvar, _iatt, String(), 1);
  return true;
}

void CalcAnamTransform::_rollback()
{
  _cleanVariableDb(1);
}

bool CalcAnamTransform::_run()
{
  int number = _getNVar();
  const AnamContinuous* anamC = dynamic_cast<const AnamContinuous*>(getAnam());
  for (int item = 0; item < number; item++)
  {
    VectorDouble z = getDb()->getColumnByLocator(ELoc::Z, item, true);
    if (z.size() <= 0) continue;
    VectorDouble y = anamC->RawToGaussianVector(z);
    getDb()->setColumnByUID(y, _iatt + item, true);
  }
  return true;
}

int RawToGaussian(Db *db,
                  AAnam* anam,
                  const ELoc &locatorType,
                  const NamingConvention &namconv)
{
  CalcAnamTransform transfo(anam);
  transfo.setDb(db);
  transfo.setNamingConvention(namconv);

  // Run the calculator
  int error = (transfo.run()) ? 0 : 1;
  return error;
}
