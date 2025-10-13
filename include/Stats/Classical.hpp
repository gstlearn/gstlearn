/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES Paris / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://gstlearn.org                                              */
/* License: BSD 3-clause                                                      */
/*                                                                            */
/******************************************************************************/
#pragma once

#include "gstlearn_export.hpp"

#include "Enum/EStatOption.hpp"

#include "Matrix/MatrixSymmetric.hpp"
#include "Basic/VectorNumT.hpp"
#include "Basic/NamingConvention.hpp"

namespace gstlrn
{

class Polygons;
class Table;
class VarioParam;
class Db;
class DbGrid;

GSTLEARN_EXPORT VectorString statOptionToName(const std::vector<EStatOption>& opers);
GSTLEARN_EXPORT std::vector<EStatOption> KeysToStatOptions(const VectorString& opers);

/**
 * \defgroup STATS Statistics calculations
 *
 **/

/** @addtogroup STATS_0 Statistics on variables of a Db
 * \ingroup STATS
 *
 * @param db         Db structure
 * @param names      Vector of n describing the target variables
 * @param title      Title given to the output (if defined)
 *  @{
 */
GSTLEARN_EXPORT Table dbStatisticsMono(Db *db,
                                       const VectorString &names,
                                       const std::vector<EStatOption> &opers = EStatOption::fromKeys({ "MEAN" }),
                                       bool flagIso = true,
                                       double proba = TEST,
                                       double vmin = TEST,
                                       double vmax = TEST,
                                       const String& title = "");
GSTLEARN_EXPORT Table dbStatisticsCorrel(Db *db,
                                         const VectorString &names,
                                         bool flagIso = true,
                                         const String& title = "");
GSTLEARN_EXPORT void dbStatisticsPrint(const Db *db,
                                       const VectorString &names,
                                       const std::vector<EStatOption> &opers = EStatOption::fromKeys({ "MEAN" }),
                                       bool flagIso = false,
                                       bool flagCorrel = false,
                                       const String &title = "",
                                       const String &radix = "");
GSTLEARN_EXPORT Table dbStatisticsMulti(Db *db,
                                        const VectorString &names,
                                        const EStatOption &oper = EStatOption::fromKey("MEAN"),
                                        bool flagMono = true,
                                        const String& title = "");

/**@}*/

/** @addtogroup STATS_1 Statistics from Db to DbGrid
 * \ingroup STATS
 *
 * @param  db         Input Db
 * @param  dbgrid     Output Grid Db
 * @param  name1      Name of the primary variable
 * @param  name2      Name of the secondary variable
 * @param  oper       Statistical operator
 * @param  cuts       Array of cutoffs (when needed)
 *  @{
 */
GSTLEARN_EXPORT VectorDouble dbStatisticsPerCell(Db *db,
                                                 DbGrid *dbgrid,
                                                 const EStatOption &oper,
                                                 const String& name1,
                                                 const String& name2 = "",
                                                 const VectorDouble &cuts = VectorDouble());
/**@}*/

GSTLEARN_EXPORT Id statisticsProportion(DbGrid *dbin,
                                         DbGrid *dbout,
                                         Id pos,
                                         Id nfacies,
                                         Id radius);
GSTLEARN_EXPORT Id statisticsTransition(DbGrid *dbin,
                                         DbGrid *dbout,
                                         Id pos,
                                         Id nfacies,
                                         Id radius,
                                         Id orient);

GSTLEARN_EXPORT VectorDouble dbStatisticsFacies(Db *db);
GSTLEARN_EXPORT double dbStatisticsIndicator(Db *db);

GSTLEARN_EXPORT MatrixSquare* sphering(const AMatrix* X);

GSTLEARN_EXPORT VectorVectorInt correlationPairs(Db *db1,
                                                 Db *db2,
                                                 const String &name1,
                                                 const String &name2,
                                                 bool flagFrom1 = false,
                                                 bool verbose = false);
GSTLEARN_EXPORT VectorVectorInt hscatterPairs(Db *db,
                                              const String &name1,
                                              const String &name2,
                                              VarioParam *varioparam,
                                              Id ilag = 0,
                                              Id idir = 0,
                                              bool verbose = false);
GSTLEARN_EXPORT Id correlationIdentify(Db *db1,
                                        Db *db2,
                                        Id icol1,
                                        Id icol2,
                                        Polygons *polygon);
GSTLEARN_EXPORT VectorVectorDouble condexp(Db *db1,
                                           Db *db2,
                                           Id icol1,
                                           Id icol2,
                                           double mini,
                                           double maxi,
                                           Id nclass,
                                           bool verbose = false);

GSTLEARN_EXPORT std::map<Id, Id> contingencyTable(const VectorInt& values);
GSTLEARN_EXPORT std::map<Id, std::map<Id, Id>>contingencyTable2(const VectorInt& values, const VectorInt& bins);
GSTLEARN_EXPORT MatrixSymmetric dbVarianceMatrix(const Db* db);

#ifndef SWIG
// All the following functions assume that the variables in the output Db used
// to store the results are already created. This is the reason why they are
// not supposed to be presented to the Target Language.

/** @addtogroup STATS_2 Statistics stored in already created variables
 * \ingroup STATS
 *
 * @param db Input Data Base
 * @param names List of target variables
 * @param opers List of statistical operators
 * @param iptr0 Starting address for storage
 * @param proba Probability (used for calculations)
 * @param vmin Minimum threshold (or TEST)
 * @param vmax Maximum threshold (or TEST)
 *
 * @{
 *
 */
GSTLEARN_EXPORT void dbStatisticsVariables(Db* db,
                                           const VectorString& names,
                                           const std::vector<EStatOption>& opers,
                                           Id iptr0,
                                           double proba = TEST,
                                           double vmin  = TEST,
                                           double vmax  = TEST);
/**@}*/

GSTLEARN_EXPORT Id dbStatisticsInGridTool(Db* db,
                                           DbGrid* dbgrid,
                                           const VectorString& names,
                                           const EStatOption& oper,
                                           Id radius,
                                           Id iptr0);

#endif // SWIG
}
