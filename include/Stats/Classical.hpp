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

#include "Matrix/MatrixSquareSymmetric.hpp"
#include "Basic/VectorNumT.hpp"
#include "Basic/NamingConvention.hpp"

class Db;
class DbGrid;
class Polygons;
class Table;
class VarioParam;

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

GSTLEARN_EXPORT int statisticsProportion(DbGrid *dbin,
                                         DbGrid *dbout,
                                         int pos,
                                         int nfacies,
                                         int radius);
GSTLEARN_EXPORT int statisticsTransition(DbGrid *dbin,
                                         DbGrid *dbout,
                                         int pos,
                                         int nfacies,
                                         int radius,
                                         int orient);

GSTLEARN_EXPORT VectorDouble dbStatisticsFacies(Db *db);
GSTLEARN_EXPORT double dbStatisticsIndicator(Db *db);

GSTLEARN_EXPORT MatrixSquareGeneral* sphering(const AMatrix* X);

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
                                              int ipas = 0,
                                              int idir = 0,
                                              bool verbose = false);
GSTLEARN_EXPORT int correlationIdentify(Db *db1,
                                        Db *db2,
                                        int icol1,
                                        int icol2,
                                        Polygons *polygon);
GSTLEARN_EXPORT VectorVectorDouble condexp(Db *db1,
                                           Db *db2,
                                           int icol1,
                                           int icol2,
                                           double mini,
                                           double maxi,
                                           int nclass,
                                           bool verbose = false);

#ifndef SWIG
// All the following functions assume that the variables in the output Db used to store
// the results are already created.
// This is the reason why they are not supposed to be presented to the Target Language.

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
GSTLEARN_EXPORT void dbStatisticsVariables(Db *db,
                                           const VectorString &names,
                                           const std::vector<EStatOption> &opers,
                                           int iptr0,
                                           double proba = TEST,
                                           double vmin = TEST,
                                           double vmax = TEST);
/**@}*/

GSTLEARN_EXPORT int dbStatisticsInGridTool(Db *db,
                                           DbGrid *dbgrid,
                                           const VectorString &names,
                                           const EStatOption &oper,
                                           int radius,
                                           int iptr0);

#endif // SWIG
