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

#include "Basic/Vector.hpp"
#include "Covariances/CovContext.hpp"
#include "geoslib_enum.h"

class ADriftElem;

class DriftFactory
{
public:
  static int          getDriftNumber();
  static ADriftElem*  createDriftFunc(const ENUM_DRIFTS& type, const CovContext& ctxt);
  static ADriftElem*  duplicateDriftFunc(const ADriftElem& cov);
  static void         displayList(const CovContext& ctxt);
  static int          identifyDrift(const String& symbol, ENUM_DRIFTS *type, int *rank,
                                    const CovContext& ctxt);
 };
