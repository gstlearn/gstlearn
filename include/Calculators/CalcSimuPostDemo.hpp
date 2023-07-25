/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://github.com/gstlearn                                       */
/* License: BSD 3 clauses                                                     */
/*                                                                            */
/******************************************************************************/
#pragma once

#include "gstlearn_export.hpp"

#include "Calculators/CalcSimuPost.hpp"

#include "Enum/EPostUpscale.hpp"
#include "Enum/EPostStat.hpp"

#include "Db/DbGrid.hpp"
#include "Basic/NamingConvention.hpp"
#include "Basic/VectorNumT.hpp"

class GSTLEARN_EXPORT CalcSimuPostDemo: public CalcSimuPost
{
public:
  CalcSimuPostDemo();
  CalcSimuPostDemo(const CalcSimuPostDemo &r) = delete;
  CalcSimuPostDemo& operator=(const CalcSimuPostDemo &r) = delete;
  virtual ~CalcSimuPostDemo();

  int getTransfoNvar() const override;
  VectorVectorDouble transformFunction(const VectorVectorDouble& tab) const override;
};

GSTLEARN_EXPORT int simuPostDemo(Db *dbin,
                                 DbGrid *dbout,
                                 const VectorString &names,
                                 bool flag_match = false,
                                 const EPostUpscale &upscale = EPostUpscale::fromKey("MEAN"),
                                 const std::vector<EPostStat> &stats = EPostStat::fromKeys({"MEAN"}),
                                 bool verbose = false,
                                 int rank_check = 0,
                                 const NamingConvention &namconv = NamingConvention("Post"));
