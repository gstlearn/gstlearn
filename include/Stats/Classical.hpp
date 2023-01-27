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

GSTLEARN_EXPORT void dbStatisticsVariables(Db *db,
                                           const VectorInt &iatts,
                                           const std::vector<EStatOption>& opers,
                                           int iptr0,
                                           double proba = TEST,
                                           double vmin = TEST,
                                           double vmax = TEST);
GSTLEARN_EXPORT VectorDouble dbStatisticsMonoByUID(Db *db,
                                                    const VectorInt &iatts,
                                                    const std::vector<
                                                        EStatOption> &opers = EStatOption::fromKeys({ "MEAN" }),
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
GSTLEARN_EXPORT VectorDouble dbStatisticsMultiByUID(Db *db,
                                                    const VectorInt &iatts,
                                                    bool flagIso = true);
GSTLEARN_EXPORT VectorDouble dbStatisticsMultiByUID(Db *db,
                                               const VectorString& names,
                                               bool flagIso = true);
GSTLEARN_EXPORT void dbStatisticsPrintByUID(const Db *db,
                                            const VectorInt &iatts = VectorInt(),
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
GSTLEARN_EXPORT VectorDouble dbStatisticsFacies(Db *db);
GSTLEARN_EXPORT double dbStatisticsIndicator(Db *db);

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
GSTLEARN_EXPORT MatrixRectangular* sphering(const AMatrix* X);

