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
#include "geoslib_enum.h"

class CovAniso;
class ACovFunc;
class CovContext;

class CovFactory
{
public:
  static int          getCovarianceNumber();
  static ACovFunc*    createCovFunc(const ENUM_COVS& type, const CovContext& ctxt);
  static ACovFunc*    duplicateCovFunc(const ACovFunc& cov);
  static void         displayList(const CovContext& ctxt);
  static VectorString getCovList(const CovContext& ctxt);
  static int identifyCovariance(const String& cov_name,
                                ENUM_COVS *rank,
                                const CovContext& ctxt);
};

