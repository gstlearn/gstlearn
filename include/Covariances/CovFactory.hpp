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

#include "Enum/ECov.hpp"

#include "Basic/Vector.hpp"

class CovAniso;
class ACovFunc;
class CovContext;

class GSTLEARN_EXPORT CovFactory
{
public:
  static ACovFunc*    createCovFunc(const ECov& type, const CovContext& ctxt);
  static ACovFunc*    duplicateCovFunc(const ACovFunc& cov);
  static void         displayList(const CovContext& ctxt);
  static VectorString getCovList(const CovContext& ctxt);
  static ECov         identifyCovariance(const String& cov_name,
                                         const CovContext& ctxt);
};

