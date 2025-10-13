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

#include "Matrix/MatrixT.hpp"
#include "gstlearn_export.hpp"

#include "Basic/VectorNumT.hpp"
#include "memory"

namespace gstlrn
{
class Db;

/**
 * \brief
 * Class returning the list of sample IDs for a quick search within a Db
 *
 * The main functionality of this class is to return the list of samples
 * per variable, within a given list of elligible sample ranks (neighborhood)
 */
class GSTLEARN_EXPORT RankHandler
{
public:
  RankHandler(const Db* db = nullptr,
              bool useSel = true,
              bool useZ   = true,
              bool useVerr = false,
              bool useExtD = true);
  RankHandler(const RankHandler& r);
  RankHandler& operator=(const RankHandler& r);
  virtual ~RankHandler();

  void defineSampleRanks(const VectorInt& nbgh = VectorInt());

  const VectorInt& getSampleRanks(Id ivar) const { return _index[ivar]; }
  const VectorVectorInt& getSampleRanks() const { return _index; }
  VectorInt& getSampleRanksByVariable(Id ivar)  { return _index[ivar]; }
  std::shared_ptr<VectorDouble>& getZflatten()  { return _Zflatten; }
  Id getNumber() const;
  Id getCount(Id ivar) const;
  Id getTotalCount() const;
  Id identifyVariableRank(Id ipos) const;
  Id identifySampleRank(Id ipos) const;

  void dump(bool flagFull = false) const;

private:
  void _initElligible();

private:
  bool _useSel;
  bool _useZ;
  bool _useVerr;
  bool _useExtD;

  Id  _nvar;
  Id  _nExtD;
  Id  _iptrSel;
  VectorInt _iptrZ;
  VectorInt _iptrVerr;
  VectorInt _iptrExtD;
  MatrixT<bool> _elligible; 
  constvectint _nbgh; // Span of internal buffer

  VectorVectorInt _index; // Vector of sample ranks per variable
  std::shared_ptr<VectorDouble> _Zflatten; // Vector of Z values (fpr active samples of target variables)

  const Db* _db;       // Pointer to Db
  VectorInt _workNbgh; // Vector of ellible sample absolute ranks
};
}