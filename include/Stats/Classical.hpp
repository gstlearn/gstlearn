/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Authors: <authors>                                                         */
/* Website: <website>                                                         */
/* License: BSD 3 clause                                                      */
/*                                                                            */
/******************************************************************************/
#pragma once

#include "gstlearn_export.hpp"

#include "Enum/EStatOption.hpp"

#include "Matrix/MatrixSquareSymmetric.hpp"
#include "Basic/VectorNumT.hpp"
#include "Basic/NamingConvention.hpp"

typedef struct
{
  int count;
  int nvar;
  bool flagCste;
  VectorDouble coeffs;
  double variance;
  double varres;
} ResRegr;

class Db;
class Table;

GSTLEARN_EXPORT VectorString statOptionToName(const std::vector<EStatOption>& opers);
GSTLEARN_EXPORT std::vector<EStatOption> KeysToStatOptions(const VectorString& opers);
GSTLEARN_EXPORT VectorDouble regrDeming(const VectorDouble &x,
                                        const VectorDouble &y,
                                        double delta = 1);


/**
 * \defgroup Stats2 Monovariate Statistics on variables
 *
 * \brief Compute several monovariate statistics (opers) on variables contained in a Db and produce the results
 * in a convenient output format (Vector of values or Table).
 *
 * @param  db         Db structure
 * @param  opers      List of the operator ranks
 * @param  flagIso    Restrain statistics to isotopic samples
 * @param  proba      Probability value (between 0 and 1)
 * @param  vmin       Minimum threshold
 * @param  vmax       Maximum threshold
 *
 *  @{
 */
GSTLEARN_EXPORT VectorDouble dbStatisticsMonoByUID(Db *db,
                                                   const VectorInt &iuids,
                                                   const std::vector<EStatOption> &opers = EStatOption::fromKeys(
                                                       { "MEAN" }),
                                                   bool flagIso = true,
                                                   double proba = TEST,
                                                   double vmin = TEST,
                                                   double vmax = TEST);
GSTLEARN_EXPORT VectorDouble dbStatisticsMono(Db *db,
                                              const VectorString& names,
                                              const std::vector<EStatOption>& opers = EStatOption::fromKeys({"MEAN"}),
                                              bool flagIso = true,
                                              double proba = TEST,
                                              double vmin = TEST,
                                              double vmax = TEST);
GSTLEARN_EXPORT Table dbStatisticsMonoT(Db *db,
                                        const VectorString &names,
                                        const std::vector<EStatOption> &opers = EStatOption::fromKeys(
                                            { "MEAN" }),
                                        bool flagIso = true,
                                        double proba = TEST,
                                        double vmin = TEST,
                                        double vmax = TEST);
/**@}*/

/**
 * \defgroup Stats3 Correlation Matrix on variables
 *
 * \brief Compute correlation matrix on variables contained in a Db and produce the results
 * in a convenient output format (Vector or a Matrix).
 *
 * @param  db         Db structure
 * @param  flagIso    Restrain statistics to isotopic samples
 *

 *  @{
 */
GSTLEARN_EXPORT VectorDouble dbStatisticsCorrelByUID(Db *db,
                                                     const VectorInt &iuids,
                                                     bool flagIso = true);
GSTLEARN_EXPORT VectorDouble dbStatisticsCorrel(Db *db,
                                                const VectorString &names,
                                                bool flagIso = true);

GSTLEARN_EXPORT MatrixSquareSymmetric dbStatisticsCorrelT(Db *db,
                                                          const VectorString &names,
                                                          bool flagIso = true);
/**@}*/

/**
 * \defgroup Stats4 Statistics on variables (printed)
 * \brief Compute statistics on variables contained in a Db and print out the results
 *
 * @param  db         Db structure
 * @param  opers      List of the operator ranks
 * @param  flagIso    Restrain statistics to isotopic samples
 * @param  flagCorrel Print the correlation matrix
 * @param  title      Title given to the printout
 * @param  radix      Radix given to the printout
 *  @{
 */
GSTLEARN_EXPORT void dbStatisticsPrintByUID(const Db *db,
                                            const VectorInt &iuids = VectorInt(),
                                            const std::vector<EStatOption> &opers = EStatOption::fromKeys(
                                                { "MEAN" }),
                                            bool flagIso = false,
                                            bool flagCorrel = false,
                                            const String &title = String(),
                                            const String &radix = String());
GSTLEARN_EXPORT void dbStatisticsPrint(const Db *db,
                                       const VectorString &names,
                                       const std::vector<EStatOption> &opers = EStatOption::fromKeys(
                                           { "MEAN" }),
                                       bool flagIso = false,
                                       bool flagCorrel = false,
                                       const String &title = String(),
                                       const String &radix = String());
/**@}*/

/**
 * \defgroup Stats5 Multivariate Statistic on variables
 * \brief Compute a Multivariate statistic on variables contained in a Db and return the results
 *
 * @param  db         Db structure
 * @param  oper       Statistical operator
 * @param  flagMono   When True, statistics by variable; otherwise, statistics by pair of variables
 * @param  verbose    Verbosity flag
 *
 * @return the vector of results
 *  @{
 */
GSTLEARN_EXPORT VectorDouble dbStatisticsMultiByColIdx(Db *db,
                                                       const VectorInt &cols,
                                                       const EStatOption &oper = EStatOption::MEAN,
                                                       bool flagMono = true,
                                                       bool verbose = false);
GSTLEARN_EXPORT VectorDouble dbStatisticsMulti(Db *db,
                                               const VectorString &names,
                                               const EStatOption &oper = EStatOption::MEAN,
                                               bool flagMono = true,
                                               bool verbose = false);
GSTLEARN_EXPORT Table dbStatisticsMultiT(Db *db,
                                         const VectorString &names,
                                         const EStatOption &oper = EStatOption::MEAN,
                                         bool flagMono = true,
                                         bool verbose = false);
/**@}*/

/**
 * \defgroup Stats6 Statistics of points per cell
 * \brief Statistics of points per cell
 *
 * @param  db         Input Db
 * @param  dbgrid     Output Grid Db
 * @param  oper       Statistical operator
 * @param  cuts       Array of cutoffs (when needed)
 * @return Vector of results
 *  @{
 */
GSTLEARN_EXPORT VectorDouble dbStatisticsPerCellByUID(Db *db,
                                                      DbGrid *dbgrid,
                                                      const EStatOption &oper,
                                                      int iuid,
                                                      int juid = -1,
                                                      const VectorDouble &cuts = VectorDouble());
GSTLEARN_EXPORT VectorDouble dbStatisticsPerCell(Db *db,
                                                 DbGrid *dbgrid,
                                                 const EStatOption &oper,
                                                 const String& name1,
                                                 const String& name2 = String(),
                                                 const VectorDouble &cuts = VectorDouble());
/**@}*/


GSTLEARN_EXPORT String statisticsMonoPrint(const VectorDouble &tab,
                                           const std::vector<EStatOption>& opers = EStatOption::fromKeys({"MEAN"}),
                                           const VectorString& names = VectorString(),
                                           const String &title = "");
GSTLEARN_EXPORT String statisticsMultiPrint(const VectorDouble &cov,
                                            const VectorString& names = VectorString(),
                                            const String &title = "");

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

GSTLEARN_EXPORT ResRegr regressionByUID(Db *db1,
                                        int icol0,
                                        const VectorInt &icols = VectorInt(),
                                        int mode = 0,
                                        bool flagCste = false,
                                        Db *db2 = nullptr,
                                        const Model *model = nullptr,
                                        bool verbose = false);
GSTLEARN_EXPORT ResRegr regression(Db *db1,
                                   const String& name0,
                                   const VectorString& names = VectorString(),
                                   int mode = 0,
                                   bool flagCste = false,
                                   Db *db2 = nullptr,
                                   const Model *model = nullptr,
                                   bool verbose = false);
GSTLEARN_EXPORT int regressionApply(Db *db1,
                                    int iptr0,
                                    const String &name0,
                                    const VectorString &names,
                                    int mode = 0,
                                    bool flagCste = false,
                                    Db *db2 = nullptr,
                                    const Model *model = nullptr);

GSTLEARN_EXPORT VectorDouble dbStatisticsFacies(Db *db);
GSTLEARN_EXPORT double dbStatisticsIndicator(Db *db);

GSTLEARN_EXPORT MatrixRectangular* sphering(const AMatrix* X);

#ifndef SWIG
// All the following functions assume that the variables in the output Db used to store
// the results are already created.
// This is the reason why they are not supposed to be presented to the Target Language.

/**
 * \defgroup Stats1 Statistics stored in already created variables
 *
 * \brief Store several statistics calculated on a set of variables of a Db and store them
 * in this same Db in variables already created.
 * These functions should not be used in Target Language.
 *
 * @param db Input Data Base
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
GSTLEARN_EXPORT void dbStatisticsVariablesByUID(Db *db,
                                                const VectorInt &iuids,
                                                const std::vector<EStatOption> &opers,
                                                int iptr0,
                                                double probas = TEST,
                                                double vmin = TEST,
                                                double vmax = TEST);
/**@}*/

GSTLEARN_EXPORT int dbStatisticsInGridToolByUID(Db *db,
                                                DbGrid *dbgrid,
                                                const VectorInt &iuids,
                                                const EStatOption &oper,
                                                int radius,
                                                int iptr0);
GSTLEARN_EXPORT int dbStatisticsInGridTool(Db *db,
                                           DbGrid *dbgrid,
                                           const VectorString &names,
                                           const EStatOption &oper,
                                           int radius,
                                           int iptr0);

#endif // SWIG
