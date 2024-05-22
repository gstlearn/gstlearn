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

#include <Enum/EPostUpscale.hpp>
#include <Enum/EPostStat.hpp>

#include "Db/DbGrid.hpp"
#include "Basic/NamingConvention.hpp"
#include "Basic/VectorNumT.hpp"

class GSTLEARN_EXPORT CalcSimuPost: public ACalcDbToDb
{
public:
  CalcSimuPost();
  CalcSimuPost(const CalcSimuPost &r) = delete;
  CalcSimuPost& operator=(const CalcSimuPost &r) = delete;
  virtual ~CalcSimuPost();


  void setNames(VectorString names)            { _names = names; }
  void setNfact(VectorInt nfact)               { _nfact = nfact; }
  void setUpscale(const EPostUpscale &upscale) { _upscale = upscale; }
  void setVerbose(bool verbose)                { _verbose = verbose; }
  void setFlagMatch(bool match)                { _flagMatch = match; }
  void setFlagUpscale(bool flagUpscale)        { _flagUpscale = flagUpscale; }
  void setStats(std::vector<EPostStat> stats)  { _stats = stats; }
  void setCheckTargets(const VectorInt& ranks) { _checkTargets = ranks; }
  void setCheckLevel(int level)                { _checkLevel = level; }

protected:
  /// Interface for ACalcDbToDb
  virtual bool _check() override;
  virtual bool _preprocess() override;
  virtual bool _run() override;
  virtual bool _postprocess() override;
  virtual void _rollback() override;

protected:
  virtual int _getTransfoNvar() const { return 0; }
  virtual void _transformFunction(const VectorDouble& tabin, VectorDouble& tabout) const { DECLARE_UNUSED(tabin, tabout); return; }

  int  _getNVar() const override { return (int) _names.size(); }
  int  _getIechout() const { return _iechout; }
  bool _getFlagUpscale() const { return _flagUpscale; }

private:
  int  _defineNames();
  void _defineIterations();
  int  _process();
  int  _getNiter() const { return _niter; }

  int  _getNVarout() const { return _nvarOut; }
  int  _getNStats() const { return (int) _stats.size(); }
  int  _getNEff() const;

  VectorVectorInt _getIndices() const;
  VectorInt _samplesInCellIdenticalSpaceDimension(const VectorInt& indblock) const;
  VectorInt _samplesInCellDifferentSpaceDimension() const;
  void _upscaleFunction(const VectorVectorDouble& Y_p_k_s, VectorDouble& tabout) const;
  void _readIn(int iech, const VectorInt& indices, VectorDouble& tabin) const;
  void _statisticsFunction(const VectorVectorDouble& Y_p, VectorDouble& tabout) const;
  void _printIndices(const VectorVectorInt &indices) const;
  int  _defineVaroutNumber();
  void _writeOut(int iech, const VectorDouble& tabout) const;
  void _environPrint() const;
  bool _mustBeChecked(int level = 0) const;
  int  _getSortingCase() const;

private:
  bool _verbose;
  bool _flagMatch;
  bool _flagUpscale;
  int  _checkLevel;
  VectorInt _checkTargets;
  EPostUpscale _upscale;
  std::vector<EPostStat> _stats;
  VectorString _names;

  mutable int  _iechout;
  mutable int  _iter;
  mutable int  _iattOut;
  mutable int  _niter;
  mutable int  _nvarOut;
  mutable VectorInt _nfact;
  mutable VectorVectorInt _iuids;
};

GSTLEARN_EXPORT int simuPost(Db *dbin,
                             DbGrid *dbout,
                             const VectorString& names,
                             bool flag_match = false,
                             const EPostUpscale& upscale = EPostUpscale::fromKey("MEAN"),
                             const std::vector<EPostStat>& stats = EPostStat::fromKeys({"MEAN"}),
                             bool verbose = false,
                             const VectorInt& check_targets = VectorInt(),
                             int check_level = 0,
                             const NamingConvention &namconv = NamingConvention("Post"));
