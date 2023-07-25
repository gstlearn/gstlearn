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

  virtual int getTransfoNvar() const;
  virtual VectorVectorDouble transformFunction(const VectorVectorDouble& tab) const;

  void setNames(VectorString names)            { _names = names; }
  void setNfact(VectorInt nfact)               { _nfact = nfact; }
  void setUpscale(const EPostUpscale &upscale) { _upscale = upscale; }
  void setVerbose(bool verbose)                { _verbose = verbose; }
  void setFlagMatch(bool match)                { _flagMatch = match; }
  void setStats(std::vector<EPostStat> stats)  { _stats = stats; }
  void setRankCheck(int rankCheck)             { _rankCheck = rankCheck; }

  int  getIechout() const { return _iechout; }

protected:
  virtual bool _check() override;
  virtual bool _preprocess() override;
  virtual bool _run() override;
  virtual bool _postprocess() override;
  virtual void _rollback() override;

private:
  int  _defineNames();
  void _defineIterations();
  VectorInt _getIndices(int rank) const;
  int _process();
  int _getNiter() const { return _niter; }
  int _getNVar() const { return (int) _names.size(); }
  int _getNVarout() const { return _nvarOut; }
  int _getNStats() const { return (int) _stats.size(); }

  VectorDouble _readIn(int iech, const VectorInt& indices) const;
  VectorDouble _upscaleFunction(const VectorVectorDouble& Y_p_k_s) const;
  void _statisticsFunction(const VectorVectorDouble& Y_p, VectorDouble& tabout) const;
  void _printIndices(int rank, const VectorInt &indices) const;
  int  _defineVaroutNumber();
  void _writeOut(int iech, const VectorDouble& tabout) const;
  void _environPrint() const;
  bool _mustBeChecked() const;

private:
  bool _verbose;
  bool _flagMatch;
  int  _rankCheck;
  EPostUpscale _upscale;
  std::vector<EPostStat> _stats;
  VectorString _names;

  mutable int _iechout;
  mutable int _iattOut;
  mutable int _niter;
  mutable int _nvarOut;
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
                             int rank_check = 0,
                             const NamingConvention &namconv = NamingConvention("Post"));
