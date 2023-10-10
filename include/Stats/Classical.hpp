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
class Table;

GSTLEARN_EXPORT VectorString statOptionToName(const std::vector<EStatOption>& opers);
GSTLEARN_EXPORT std::vector<EStatOption> KeysToStatOptions(const VectorString& opers);

/**
 * \defgroup Stats1 Statistics on variables of a Db
 *
 * \brief Compute one or several statistics (oper/opers) on a set of variables (names)
 * contained in a Db and produce the results as a Table or as variables added to the input Db
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
                                         const String& title = String());
GSTLEARN_EXPORT void dbStatisticsPrint(const Db *db,
                                       const VectorString &names,
                                       const std::vector<EStatOption> &opers = EStatOption::fromKeys({ "MEAN" }),
                                       bool flagIso = false,
                                       bool flagCorrel = false,
                                       const String &title = String(),
                                       const String &radix = String());
GSTLEARN_EXPORT Table dbStatisticsMulti(Db *db,
                                        const VectorString &names,
                                        const EStatOption &oper = EStatOption::MEAN,
                                        bool flagMono = true,
                                        const String& title = String());

/**@}*/

/**
 * \defgroup Stats3 Statistics from Db to DbGrid
 *
 * \brief Compute Statistics of points per cell
 *
 * @param  db         Input Db
 * @param  dbgrid     Output Grid Db
 * @param  oper       Statistical operator
 * @param  cuts       Array of cutoffs (when needed)
 * @return Vector of results
 *  @{
 */
GSTLEARN_EXPORT VectorDouble dbStatisticsPerCell(Db *db,
                                                 DbGrid *dbgrid,
                                                 const EStatOption &oper,
                                                 const String& name1,
                                                 const String& name2 = String(),
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

GSTLEARN_EXPORT MatrixRectangular* sphering(const AMatrix* X);

#ifndef SWIG
// All the following functions assume that the variables in the output Db used to store
// the results are already created.
// This is the reason why they are not supposed to be presented to the Target Language.

/**
 * \defgroup Stats4 Statistics stored in already created variables
 *
 * \brief Store several statistics calculated on a set of variables of a Db and store them
 * in this same Db in variables already created.
 * These functions should not be used in Target Language.
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
