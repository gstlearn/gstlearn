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
#include "LinearOp/ALinearOp.hpp"
#include "Basic/VectorNumT.hpp"

ALinearOp::ALinearOp()
: _usefactor(false)
, _idfactor(0.)
, _factor(1.)
{
}

VectorDouble ALinearOp::evalDirect(const VectorDouble& in) const
{
  VectorDouble res;
  evalDirect(in,res);
  return res;
}

int ALinearOp::addToDest(const Eigen::VectorXd& inv,
                         Eigen::VectorXd& outv) const
{
  constvect ins(inv.data(), inv.size());
  vect outs(outv.data(), outv.size());
  return addToDest(ins, outs);
}

int ALinearOp::addToDest(const constvect inv, vect outv) const
{
  
  if (!_usefactor)
    return _addToDest(inv, outv);
  
  _temp.resize(outv.size());
  vect ctemp(_temp.data(),_temp.size());
  std::fill(ctemp.begin(), ctemp.end(), 0.);
  int err = _addToDest(inv, ctemp);
  for (int i = 0; i < (int)outv.size(); i++)
  {
    outv[i] = _idfactor * inv[i] + _factor * ctemp[i];
  }
  return err;
}

int ALinearOp::evalDirect(constvect inv, vect outv) const
{
  std::fill(outv.begin(), outv.end(), 0.);
  return addToDest(inv, outv);
}

int ALinearOp::evalDirect(const VectorDouble& inv, VectorDouble& outv) const
{ 
  outv.resize(inv.size());
  constvect in(inv);
  vect out(outv);
  return evalDirect(in,out);
}

void ALinearOp::multiplyByValueAndAddDiagonal(double v1,double v2)
{
  _usefactor = true;
  _idfactor = v2;
  _factor = v1;
}

void ALinearOp::resetModif()
{
  _usefactor = false;
  _idfactor = 0.;
  _factor = 1.;
}
