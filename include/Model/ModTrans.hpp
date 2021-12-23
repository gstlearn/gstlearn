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
#pragma once

#include "gstlearn_export.hpp"
#include "Model/EModelProperty.hpp"
#include "Anamorphosis/Anam.hpp"
#include "Basic/Vector.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/IClonable.hpp"

class GSTLEARN_EXPORT ModTrans : public AStringable, public IClonable
{
public:
  ModTrans();
  ModTrans(const ModTrans &m);
  ModTrans& operator= (const ModTrans &m);
  virtual ~ModTrans();

  virtual String toString(int level = 0) const override;
  virtual IClonable* clone() const override { return new ModTrans(*this); };

  void cancelProperty();

  const EModelProperty& getModTransMode() const { return _modTransMode; }

private:
  EModelProperty _modTransMode;
};
