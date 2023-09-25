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
#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <iterator>

#include "Drifts/DriftF.hpp"
#include "Drifts/ADriftElem.hpp"
#include "Db/Db.hpp"

DriftF::DriftF(int rank_fex, const CovContext& ctxt)
    : ADriftElem(ctxt),
      _rankFex(rank_fex)
{
}

DriftF::DriftF(const DriftF &r)
    : ADriftElem(r),
      _rankFex(r._rankFex)
{
}

DriftF& DriftF::operator=(const DriftF &r)
{
  if (this != &r)
  {
    ADriftElem::operator =(r);
    _rankFex = r._rankFex;
  }
  return *this;
}

DriftF::~DriftF()
{
}

double DriftF::eval(const Db* db, int iech) const
{
  return db->getLocVariable(ELoc::F,iech,_rankFex);
}

String DriftF::getDriftName() const
{
  std::stringstream sstr;
  sstr << "External_Drift:" << _rankFex;
  return sstr.str();
}

DriftF* DriftF::createByIdentifier(const String &driftname,
                                   const CovContext &ctxt)
{
  String substring = {"External_Drift:"};

  std::size_t found = driftname.find(substring);
  if (found != 0) return nullptr;
  String string_rank = driftname.substr(substring.size(), driftname.size()-1);
  int rank_fex = atoi(string_rank.c_str());
  return new DriftF(rank_fex, ctxt);
}
