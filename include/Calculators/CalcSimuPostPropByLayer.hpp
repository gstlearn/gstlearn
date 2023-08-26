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

#include "Calculators/CalcSimuPost.hpp"

#include "Enum/EPostUpscale.hpp"
#include "Enum/EPostStat.hpp"

#include "Db/DbGrid.hpp"
#include "Basic/NamingConvention.hpp"
#include "Basic/VectorNumT.hpp"

/**
 * This particular Multivariate Simulation post_processing considers each simulated variable as a thickness
 * of ordered layers (define in R_N). For each cell of the output grid (defined in R_{N+1}, we calculate the
 * proportion of each layer within the cell
 */
class GSTLEARN_EXPORT CalcSimuPostPropByLayer: public CalcSimuPost
{
public:
  CalcSimuPostPropByLayer();
  CalcSimuPostPropByLayer(const CalcSimuPostPropByLayer &r) = delete;
  CalcSimuPostPropByLayer& operator=(const CalcSimuPostPropByLayer &r) = delete;
  virtual ~CalcSimuPostPropByLayer();

  void setFlagTopToBase(bool topToBase)        { _flagTopToBase = topToBase; }

protected:
  bool _check() override;

  int _getTransfoNvar() const override;
  void _transformFunction(const VectorDouble& tabin, VectorDouble& tabout) const override;

private:
  const DbGrid* _dbgrid; //!< Pointer to 'dbout'. Not to be freed
  bool _flagTopToBase;   // define the direction to compute the proportion

};

GSTLEARN_EXPORT int simuPostPropByLayer(Db *dbin,
                                        DbGrid *dbout,
                                        const VectorString &names,
                                        bool flag_match = false,
										bool flag_topToBase = false,
                                        const EPostUpscale &upscale = EPostUpscale::fromKey("MEAN"),
                                        const std::vector<EPostStat> &stats = EPostStat::fromKeys({"MEAN"}),
                                        bool verbose = false,
                                        const VectorInt& check_targets = VectorInt(),
                                        int check_level = 0,
                                        const NamingConvention &namconv = NamingConvention("Prop"));
