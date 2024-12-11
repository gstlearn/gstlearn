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
#include "LinearOp/AShiftOp.hpp"  
#include "Covariances/CovAniso.hpp"
#include "LinearOp/ALinearOp.hpp"
#include <math.h>

AShiftOp::AShiftOp(CovAniso* cova, int napices)
: _Lambda()
, _napices(napices)
, _cova(cova)
{
}


AShiftOp::AShiftOp(const AShiftOp& shift)
: _Lambda(shift._Lambda)
, _napices(shift._napices)
, _cova(shift._cova)
{
}


AShiftOp& AShiftOp::operator=(const AShiftOp &shift)
{
  if (this != &shift)
  {
    _napices = shift._napices;
    _cova = shift._cova;
    _Lambda = shift._Lambda;
  }
  return *this;
}

AShiftOp::~AShiftOp()
{
}


void AShiftOp::prodLambda(const VectorDouble& x,
                          vect y,
                           const EPowerPT& power) const
{
  constvect xv(x.data(),x.size());
  prodLambda(xv,y,power);
  
}

void AShiftOp::prodLambda(const constvect x,
                           VectorDouble& y,
                           const EPowerPT& power) const
{
    vect yv(y.data(),y.size());
    prodLambda(x,yv,power);
}

void AShiftOp::prodLambda(const constvect x,
                           vect y,
                           const EPowerPT& power) const
{
  std::fill(y.begin(),y.end(),0.);
  addProdLambda(x,y,power);
}

void AShiftOp::prodLambda(const VectorDouble& x,
                           VectorDouble& y,
                           const EPowerPT& power) const
{
  constvect xv(x.data(),x.size());
  vect yv(y.data(),y.size());
  prodLambda(xv,yv,power);
}

void AShiftOp::addProdLambda(const constvect x,
                             vect y,
                             const EPowerPT& power) const
{
  if (power == EPowerPT::ONE)
  {
    for (int i = 0, n = getSize(); i < n; i++)
      y[i] += x[i] * getLambda(i);
  }
  else if (power == EPowerPT::MINUSONE)
  {
    for (int i = 0, n = getSize(); i < n; i++)
      y[i] += x[i] / getLambda(i);
  }
  else if (power == EPowerPT::HALF)
  {
    for (int i = 0, n = getSize(); i < n; i++)
      y[i] += x[i] * sqrt(getLambda(i));
  }
  else if (power == EPowerPT::MINUSHALF)
  {
    for (int i = 0, n = getSize(); i < n; i++)
      y[i] += x[i] / sqrt(getLambda(i));
  }
  else
  {
    my_throw("Unexpected value for argument 'power'");
  }
}

std::shared_ptr<CovAniso> AShiftOp::cloneAndCast(const std::shared_ptr<CovAniso> &cova)
{
    return std::shared_ptr<CovAniso>((CovAniso*)cova->clone());

}

std::shared_ptr<CovAniso> AShiftOp::cloneAndCast(const CovAniso* cova)
{
    return std::shared_ptr<CovAniso>((CovAniso*)cova->clone());

}

void AShiftOp::normalizeLambdaBySills(const AMesh* mesh)
{
  VectorDouble tab;

  if (_cova->isNoStatForVariance())
  {
    _cova->informMeshByApexForSills(mesh);
    int number = (int) _Lambda.size();
                       
    
    for (int imesh = 0; imesh < number; imesh++)
    {
      _cova->updateCovByMesh(imesh,false);
      double sill = _cova->getSill(0,0);
      double invsillsq = 1. / sqrt(sill);
      _Lambda[imesh] *= invsillsq;
    }
  }
  else 
  {
    double invsillsq = 1. / sqrt(_cova->getSill(0,0));
    for (auto &e:_Lambda)
    {
      e *= invsillsq;
    }
  }
}

bool AShiftOp::_isNoStat()
{
  return _getCovAniso()->isNoStat();
}

bool AShiftOp::_isGlobalHH()
{
  return !_cova->isNoStatForAnisotropy();
}


void AShiftOp::_setCovAniso(const CovAniso* cova)
{
  _cova = cloneAndCast(cova);
}

std::shared_ptr<CovAniso> &AShiftOp::_getCovAniso()
{
  return _cova;
}
