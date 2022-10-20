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

#include "Basic/Vector.hpp"
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

GSTLEARN_EXPORT VectorString statsNames(const std::vector<EStatOption>& opers);

GSTLEARN_EXPORT void dbStatisticsVariables(Db *db,
                                           const VectorInt &iatts,
                                           const std::vector<EStatOption>& opers,
                                           int iptr0,
                                           double proba = TEST,
                                           double vmin = TEST,
                                           double vmax = TEST);

GSTLEARN_EXPORT VectorDouble dbStatisticsMono(Db *db,
                                              const VectorInt &iatts,
                                              const std::vector<EStatOption>& opers = {EStatOption::MEAN},
                                              bool flagIso = true,
                                              double proba = TEST,
                                              double vmin = TEST,
                                              double vmax = TEST);
GSTLEARN_EXPORT VectorDouble dbStatisticsMono(Db *db,
                                              const VectorString& names,
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

GSTLEARN_EXPORT bool regressionCheck(Db *db1,
                                     int icol0,
                                     const VectorInt &icols,
                                     int mode,
                                     Db *db2,
                                     const Model* model = nullptr);
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
GSTLEARN_EXPORT bool regressionLoad(Db *db1,
                                    Db *db2,
                                    int iech,
                                    int icol0,
                                    const VectorInt &icols,
                                    int mode,
                                    int flagCste,
                                    const Model* model,
                                    double *value,
                                    VectorDouble &x);
GSTLEARN_EXPORT void regrprint(const ResRegr& regr);
