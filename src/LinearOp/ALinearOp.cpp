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

namespace gstlrn
{
ALinearOp::ALinearOp()
  : _usefactor(false)
  , _idfactor(0.)
  , _factor(1.)
  , _temp()
{
}

// ALinearOp::ALinearOp(const ALinearOp &m)
//   : _usefactor(m._usefactor)
//   , _idfactor(m._idfactor)
//   , _factor(m._factor)
//   , _temp(m._temp)
// {
// }

// ALinearOp& ALinearOp::operator=(const ALinearOp& m)
// {
//   if (this != &m)
//   {
//     _usefactor = m._usefactor;
//     _idfactor  = m._idfactor;
//     _factor    = m._factor;
//     _temp      = m._temp;
//   }
//   return *this;
// }

VectorDouble ALinearOp::evalDirect(const VectorDouble& in) const
{
  VectorDouble res;
  evalDirect(in, res);
  return res;
}

Id ALinearOp::addToDest(const ::Eigen::VectorXd& inv,
                        ::Eigen::VectorXd& outv) const
{
  constvect ins(inv.data(), inv.size());
  vect outs(outv.data(), outv.size());
  return addToDest(ins, outs);
}

Id ALinearOp::addToDest(const constvect inv, vect outv) const
{

  if (!_usefactor)
    return _addToDest(inv, outv);

  _temp.resize(outv.size());
  vect ctemp(_temp.data(), _temp.size());
  std::fill(ctemp.begin(), ctemp.end(), 0.);
  Id err = _addToDest(inv, ctemp);
  for (Id i = 0; i < static_cast<Id>(outv.size()); i++)
    outv[i] = _idfactor * inv[i] + _factor * ctemp[i];
  return err;
}

Id ALinearOp::evalDirect(constvect inv, vect outv) const
{
  std::fill(outv.begin(), outv.end(), 0.);
  return addToDest(inv, outv);
}

Id ALinearOp::evalDirect(const VectorDouble& inv, VectorDouble& outv) const
{
  outv.resize(inv.size());
  constvect in(inv);
  vect out(outv);
  return evalDirect(in, out);
}

void ALinearOp::multiplyByValueAndAddDiagonal(double v1, double v2) const
{
  _usefactor = true;
  _idfactor  = v2;
  _factor    = v1;
}

void ALinearOp::resetModif() const
{
  _usefactor = false;
  _idfactor  = 0.;
  _factor    = 1.;
}
} // namespace gstlrn