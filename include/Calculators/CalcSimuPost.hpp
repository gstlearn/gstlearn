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

#include "ACalcDbToDb.hpp"

#include <Enum/EPostStat.hpp>
#include <Enum/EPostUpscale.hpp>

#include "Basic/NamingConvention.hpp"
#include "Basic/VectorNumT.hpp"
#include "Db/DbGrid.hpp"

namespace gstlrn
{
class GSTLEARN_EXPORT CalcSimuPost: public ACalcDbToDb
{
public:
  CalcSimuPost();
  CalcSimuPost(const CalcSimuPost& r)            = delete;
  CalcSimuPost& operator=(const CalcSimuPost& r) = delete;
  virtual ~CalcSimuPost();

  void setNames(const VectorString& names) { _names = names; }
  void setNfact(const VectorInt& nfact) { _nfact = nfact; }
  void setUpscale(const EPostUpscale& upscale) { _upscale = upscale; }
  void setVerbose(bool verbose) { _verbose = verbose; }
  void setFlagMatch(bool match) { _flagMatch = match; }
  void setFlagUpscale(bool flagUpscale) { _flagUpscale = flagUpscale; }
  void setStats(const std::vector<EPostStat>& stats) { _stats = stats; }
  void setCheckTargets(const VectorInt& ranks) { _checkTargets = ranks; }
  void setCheckLevel(Id level) { _checkLevel = level; }

protected:
  bool _check() override;
  bool _preprocess() override;
  bool _run() override;
  bool _postprocess() override;
  void _rollback() override;

protected:
  virtual Id _getTransfoNvar() const { return 0; }
  virtual void _transformFunction(const VectorDouble& tabin, VectorDouble& tabout) const { DECLARE_UNUSED(tabin, tabout); }

  Id _getIechout() const { return _iechout; }
  bool _getFlagUpscale() const { return _flagUpscale; }

private:
  Id _defineNames();
  void _defineIterations();
  Id _process();
  Id _getNiter() const { return _niter; }

  Id _getNVarout() const { return _nvarOut; }
  Id _getNStats() const { return static_cast<Id>(_stats.size()); }
  Id _getNEff() const;

  VectorVectorInt _getIndices() const;
  VectorInt _samplesInCellIdenticalSpaceDimension(const VectorInt& indblock) const;
  VectorInt _samplesInCellDifferentSpaceDimension() const;
  void _upscaleFunction(const VectorVectorDouble& Y_p_k_s, VectorDouble& tabout) const;
  void _readIn(Id iech, const VectorInt& indices, VectorDouble& tabin) const;
  void _statisticsFunction(const VectorVectorDouble& Y_p, VectorDouble& tabout) const;
  void _printIndices(const VectorVectorInt& indices) const;
  Id _defineVaroutNumber();
  void _writeOut(Id iech, const VectorDouble& tabout) const;
  void _environPrint() const;
  bool _mustBeChecked(Id level = 0) const;
  Id _getSortingCase() const;

private:
  bool _verbose;
  bool _flagMatch;
  bool _flagUpscale;
  Id _checkLevel;
  VectorInt _checkTargets;
  EPostUpscale _upscale;
  std::vector<EPostStat> _stats;
  VectorString _names;

  mutable Id _iechout;
  mutable Id _iter;
  mutable Id _iattOut;
  mutable Id _niter;
  mutable Id _nvarOut;
  mutable VectorInt _nfact;
  mutable VectorVectorInt _iuids;
};

GSTLEARN_EXPORT Id simuPost(Db* dbin,
                             DbGrid* dbout,
                             const VectorString& names,
                             bool flag_match                     = false,
                             const EPostUpscale& upscale         = EPostUpscale::fromKey("MEAN"),
                             const std::vector<EPostStat>& stats = EPostStat::fromKeys({"MEAN"}),
                             bool verbose                        = false,
                             const VectorInt& check_targets      = VectorInt(),
                             Id check_level                     = 0,
                             const NamingConvention& namconv     = NamingConvention("Post"));
} // namespace gstlrn