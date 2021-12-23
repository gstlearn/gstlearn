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
#include "Basic/Utilities.hpp"
#include "Model/ModTrans.hpp"
#include "Model/EModelProperty.hpp"
#include "Anamorphosis/AnamDiscreteDD.hpp"
#include "Anamorphosis/AnamDiscreteIR.hpp"
#include "Anamorphosis/AnamHermite.hpp"

ModTrans::ModTrans()
  : AStringable(),
    _modTransMode(EModelProperty::NONE)
{
}

ModTrans::ModTrans(const ModTrans &m)
    : AStringable(m),
      _modTransMode(m._modTransMode)
{
}

ModTrans& ModTrans::operator=(const ModTrans &m)
{
  if (this != &m)
  {
    AStringable::operator=(m);
    _modTransMode = m._modTransMode;
  }
  return (*this);
}

ModTrans::~ModTrans()
{
}

void ModTrans::cancelProperty()
{
  _modTransMode = EModelProperty::NONE;
}

String ModTrans::toString(int level) const
{
  std::stringstream sstr;

  if (_modTransMode == EModelProperty::NONE) return sstr.str();
  mestitle(1,"Additional Properties");

  return sstr.str();
}
