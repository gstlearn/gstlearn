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
#include "Basic/Vector.hpp"
#include "Basic/NamingConvention.hpp"
#include "Stats/EStatOption.hpp"

class Db;

GSTLEARN_EXPORT VectorString statsNames(const std::vector<EStatOption>& opers);

GSTLEARN_EXPORT void dbStatisticsVariables(Db *db,
                                           const VectorInt &iatts,
                                           const std::vector<EStatOption>& opers,
                                           int iattn,
                                           double vmin = TEST,
                                           double vmax = TEST,
                                           double proba = TEST);

GSTLEARN_EXPORT VectorDouble dbStatisticsMono(Db *db,
                                              const VectorInt &iatts,
                                              const std::vector<EStatOption>& opers = {EStatOption::MEAN},
                                              bool flagIso = true,
                                              double proba = TEST,
                                              double vmin = TEST,
                                              double vmax = TEST);
GSTLEARN_EXPORT VectorDouble dbStatisticsMulti(Db *db,
                                               const VectorInt &iatts,
                                               bool flagIso = true);
GSTLEARN_EXPORT VectorDouble dbStatisticsFacies(Db *db);
GSTLEARN_EXPORT double dbStatisticsIndicator(Db *db);

GSTLEARN_EXPORT String statisticsMonoPrint(const VectorDouble &tab,
                                           const std::vector<EStatOption>& opers = {EStatOption::MEAN},
                                           const VectorString &varnames = VectorString(),
                                           const String &title = "");
GSTLEARN_EXPORT String statisticsMultiPrint(const VectorDouble &cov,
                                            const VectorString &varnames = VectorString(),
                                            const String &title = "");

