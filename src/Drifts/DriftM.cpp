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
#include "Drifts/DriftM.hpp"
#include "Drifts/ADrift.hpp"
#include "Db/Db.hpp"

DriftM::DriftM(const VectorInt &powers)
    : ADrift(),
      _monomialPower(powers)
{
}

DriftM::DriftM(const DriftM &r)
    : ADrift(r),
      _monomialPower(r._monomialPower)
{
}

DriftM& DriftM::operator=(const DriftM &r)
{
  if (this != &r)
  {
    ADrift::operator =(r);
    _monomialPower = r._monomialPower;
  }
  return *this;
}

DriftM::~DriftM()
{
}

double DriftM::eval(const Db* db, int iech) const
{
  double value = 1.;
  for (int idim = 0, ndim = (int) _monomialPower.size(); idim < ndim; idim++)
  {
    double locoor = db->getCoordinate(iech,idim);
    double locpow = _monomialPower[idim];
    value *= pow(locoor, locpow);
  }
  return value;
}

int DriftM::getOrderIRF() const
{
  int irf = -1;
  for (int idim = 0, ndim = (int) _monomialPower.size(); idim < ndim; idim++)
  {
    double locpow = _monomialPower[idim];
    if (locpow > irf) irf = locpow;
  }
  return irf;
}

int DriftM::getOrderIRFIdim(int idim) const
{
  if (idim < getDriftNDimMax()) return -1;
  return _monomialPower[idim];
}

int DriftM::getDriftNDimMax() const
{
  return (int) _monomialPower.size();
}

String DriftM::getDriftName() const
{
  std::stringstream sstr;
  if (_monomialPower.empty())
    sstr << "Universality_Condition";
  else
  {
    sstr << "Drift:";
    bool flag_first = true;
    for (int idim = 0, ndim = (int) _monomialPower.size(); idim < ndim; idim++)
    {
      double locpow = _monomialPower[idim];
      if (locpow > 0)
      {
        if (!flag_first) sstr << "*";
        sstr << "x" << idim+1;
        if (locpow > 1) sstr << "^" << locpow;
        flag_first = false;
      }
    }
  }
  return sstr.str();
}

DriftM* DriftM::createByIdentifier(const String &driftname)
{
  String input = driftname;
  String substring;
  std::size_t found;

  // Looking for Universality Condition
  substring = "Universality_Condition";
  found = input.find(substring);
  if (found == 0) return new DriftM();

  // Looking for other drift conditions
  substring = "Drift:";
  found = input.find(substring);
  if (found != 0) return nullptr;

  // Decode the rest of the string
  input = input.substr(substring.size(), input.size()-1);

  // Initiate a vector of powers of the monomials to an extreme dimension: it will be resized at the end
  VectorInt powers(10, 0);
  int rank_max = 0;
  while (input.size() > 0)
  {
    // Decode the character "x"
    substring = "x";
    found = input.find(substring);
    if (found != 0) return nullptr;
    input = input.substr(substring.size(), input.size()-1);

    // Decode the power
    int rank = atoi(input.c_str());
    input = input.substr(1, input.size()-1);
    if (rank > rank_max) rank_max = rank;

    // Attempt to read the exponentiation
    int power = 1;
    substring = "^";
    found = input.find(substring);
    if (found == 0)
    {
      // Attempt to read the exponent
      input = input.substr(substring.size(), input.size()-1);
      power = atoi(input.c_str());
      input = input.substr(1, input.size()-1);
    }

    // Attempt to read the character "*"
    substring = "*";
    found = input.find(substring);

    // Concatenate the results
    powers[rank-1] = power;

    if (found != 0) break;
    input = input.substr(substring.size(), input.size()-1);
  }

  // Final Resizing
  powers.resize(rank_max);
  return new DriftM(powers);
}
