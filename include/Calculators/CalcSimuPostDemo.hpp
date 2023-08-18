/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://gstlearn.org                                              */
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

protected:
  int _getTransfoNvar() const override;
  void _transformFunction(const VectorDouble& tabin, VectorDouble& tabout) const override;
};

GSTLEARN_EXPORT int simuPostDemo(Db *dbin,
                                 DbGrid *dbout,
                                 const VectorString &names,
                                 bool flag_match = false,
                                 const EPostUpscale &upscale = EPostUpscale::fromKey("MEAN"),
                                 const std::vector<EPostStat> &stats = EPostStat::fromKeys({"MEAN"}),
                                 bool verbose = false,
                                 const VectorInt& check_targets = VectorInt(),
                                 int check_level = 0,
                                 const NamingConvention &namconv = NamingConvention("Post"));
