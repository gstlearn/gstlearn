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
#include "Basic/NamingConvention.hpp"

class Db;

typedef enum
{
  STAT_UNKNOWN,
  STAT_NUM,
  STAT_MEAN,
  STAT_VAR,
  STAT_STDV,
  STAT_MINI,
  STAT_MAXI,
  STAT_SUM,
  STAT_PROP,
  STAT_QUANT,
  STAT_T,
  STAT_Q,
  STAT_M,
  STAT_B
} ENUM_STATS;

String statsName(int ioper);
VectorInt statsList(const VectorString& opers);
VectorString statsNames(const VectorInt& iopers);

void dbStatisticsVariables(Db *db,
                           const VectorInt& iatts,
                           const VectorInt& iopers,
                           int iattn,
                           double vmin = TEST,
                           double vmax = TEST,
                           double proba = TEST);

VectorDouble dbStatisticsMono(Db *db,
                              const VectorInt& iatts,
                              const VectorInt& iopers = VectorInt(STAT_MEAN),
                              bool flagIso = true,
                              double proba = TEST,
                              double vmin = TEST,
                              double vmax = TEST);
VectorDouble dbStatisticsMulti(Db *db,
                               const VectorInt& iatts,
                               bool flagIso = true);

String statisticsMonoPrint(const VectorDouble& tab,
                           const VectorInt& iopers = VectorInt(STAT_MEAN),
                           const VectorString& varnames = VectorString(),
                           const String& title = "");
String statisticsMultiPrint(const VectorDouble& cov,
                            const VectorString& varnames = VectorString(),
                            const String& title = "");

