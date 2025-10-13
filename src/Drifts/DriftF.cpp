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
#include "Drifts/DriftF.hpp"
#include "Db/Db.hpp"
#include "Drifts/ADrift.hpp"
#include <iostream>
#include <sstream>

#include "Drifts/DriftF.hpp"
#include "Drifts/ADrift.hpp"
#include "Db/Db.hpp"

namespace gstlrn
{
DriftF::DriftF(Id rank_fex)
  : ADrift()
  , _rankFex(rank_fex)
{
}

DriftF::DriftF(const DriftF& r)
  : ADrift(r)
  , _rankFex(r._rankFex)
{
}

DriftF& DriftF::operator=(const DriftF& r)
{
  if (this != &r)
  {
    ADrift::operator=(r);
    _rankFex = r._rankFex;
  }
  return *this;
}

DriftF::~DriftF()
{
}

double DriftF::eval(const Db* db, Id iech) const
{
  return db->getLocVariable(ELoc::F, iech, _rankFex);
}

String DriftF::getDriftName() const
{
  std::stringstream sstr;
  sstr << "External_Drift:" << _rankFex;
  return sstr.str();
}

DriftF* DriftF::createByIdentifier(const String& driftname)
{
  String substring = {"External_Drift:"};

  std::size_t found = driftname.find(substring);
  if (found != 0) return nullptr;
  String string_rank = driftname.substr(substring.size(), driftname.size() - 1);
  Id rank_fex       = atoi(string_rank.c_str());
  return new DriftF(rank_fex);
}
}