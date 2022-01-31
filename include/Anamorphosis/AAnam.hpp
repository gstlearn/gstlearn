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

#include "Anamorphosis/EAnam.hpp"

#include "Basic/AStringable.hpp"
#include "Basic/ASerializable.hpp"

class GSTLEARN_EXPORT AAnam : public AStringable, public ASerializable
{
public:
  AAnam();
  AAnam(const AAnam &m);
  AAnam& operator= (const AAnam &m);
  virtual ~AAnam();

  /// Interface for AAnam
  virtual const EAnam&  getType() const = 0;
};
