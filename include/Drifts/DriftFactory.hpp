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

#include "Enum/EDrift.hpp"

#include "Covariances/CovContext.hpp"

class ADriftElem;

class GSTLEARN_EXPORT DriftFactory
{
public:
  static ADriftElem* createDriftFunc(const EDrift &type,
                                     const CovContext &ctxt = CovContext(),
                                     int rank_fex = 0);
  static ADriftElem*  duplicateDriftFunc(const ADriftElem& cov);
  static void         displayList(const CovContext& ctxt);
  static EDrift       identifyDrift(const String& symbol, int* rank, const CovContext& ctxt);
 };
