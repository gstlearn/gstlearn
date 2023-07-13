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

  void setNames(VectorString names) { _names = names; }
  void setNfact(VectorInt nfact)    { _nfact = nfact; }
  void setVerbose(bool verbose)     { _verbose = verbose; }
  void setFlagMatch(bool match)     { _flagMatch = match; }

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
  int _getNVar() const { return (int) _names.size(); }

private:
  bool _verbose;
  bool _flagMatch;
  int  _niter;
  VectorString _names;
  VectorInt    _nfact;
  VectorVectorInt _iuids;
};

GSTLEARN_EXPORT int simuPost(Db *dbin,
                             DbGrid *dbout,
                             const VectorString& names,
                             bool flag_match = false,
                             bool verbose = false,
                             const NamingConvention &namconv = NamingConvention(
                                 "Post"));
