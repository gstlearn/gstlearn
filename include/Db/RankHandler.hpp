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

  const VectorInt& getSampleRanks(int ivar) const { return _index[ivar]; }
  const VectorVectorInt& getSampleRanks() const { return _index; }
  VectorInt& getSampleRanksByVariable(int ivar)  { return _index[ivar]; }
  std::shared_ptr<VectorDouble>& getZflatten()  { return _Zflatten; }
  int getNumber() const;
  int getCount(int ivar) const;
  int getTotalCount() const;
  int identifyVariableRank(int ipos) const;
  int identifySampleRank(int ipos) const;

  void dump(bool flagFull = false) const;

private:
  void _initElligible();

private:
  bool _useSel;
  bool _useZ;
  bool _useVerr;
  bool _useExtD;

  int  _nvar;
  int  _nExtD;
  int  _iptrSel;
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
