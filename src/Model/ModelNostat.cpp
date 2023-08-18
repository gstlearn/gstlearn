/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://gstlearn.org                                              */
/* License: BSD 3-clause                                                      */
/*                                                                            */
/******************************************************************************/
#include "geoslib_old_f.h"

#include "Model/ModelNostat.hpp"
#include "Model/ElemNostat.hpp"
#include "Basic/Utilities.hpp"
#include "Covariances/CovAniso.hpp"

ModelNostat::ModelNostat()
    : AStringable(),
      _nDim(0),
      _elems(),
      _sill1(0),
      _sill2(0),
      _param1(0),
      _param2(0),
      _scadef1(0),
      _scadef2(0),
      _angles1(),
      _angles2(),
      _scale1(),
      _scale2()
{

}

ModelNostat::ModelNostat(const ModelNostat &m)
    : AStringable(m),
      _nDim(m._nDim),
      _elems(),
      _sill1(m._sill1),
      _sill2(m._sill2),
      _param1(m._param1),
      _param2(m._param2),
      _scadef1(m._scadef1),
      _scadef2(m._scadef2),
      _angles1(m._angles1),
      _angles2(m._angles2),
      _scale1(m._scale1),
      _scale2(m._scale2)
{
  for (int i = 0; i < (int) m._elems.size(); i++)
    _elems.push_back(new ElemNostat(*m._elems[i]));
}

ModelNostat& ModelNostat::operator=(const ModelNostat &m)
{
  if (this != &m)
  {
    AStringable::operator=(m);
    _nDim = m._nDim;
    _sill1 = m._sill1;
    _sill2 = m._sill2;
    _param1 = m._param1;
    _param2 = m._param2;
    _scadef1 = m._scadef1;
    _scadef2 = m._scadef2;
    _angles1 = m._angles1;
    _angles2 = m._angles2;
    _scale1 = m._scale1;
    _scale2 = m._scale2;
    for (int i = 0; i < (int) m._elems.size(); i++)
      _elems.push_back(new ElemNostat(*m._elems[i]));
  }
  return *this;
}

ModelNostat::~ModelNostat()
{
  for (int i = 0; i < (int) _elems.size(); i++)
  {
    delete _elems[i];
  }
}

void ModelNostat::init(int ndim)
{
  _nDim    = ndim;
  _elems.resize(0);
  _angles1.resize(ndim);
  _angles2.resize(ndim);
  _scale1.resize(ndim);
  _scale2.resize(ndim);
}

void ModelNostat::init(const CovAniso* cova)
{
  if (cova->getFlagAniso() && cova->getFlagRotation())
    _angles1 = cova->getAnisoAngles();
  else
    for (int idim = 0; idim < _nDim; idim++)
      _angles1[idim] = 0.;

  _angles2 = _angles1;

  _sill1 = cova->getSill(0, 0);
  _sill2 = cova->getSill(0, 0);
  _param1 = cova->getParam();
  _param2 = cova->getParam();
  _scadef1 = cova->getScadef();
  _scadef2 = cova->getScadef();
  _scale1 = cova->getRanges();
  _scale2 = cova->getRanges();

  for (int idim = 0; idim < _nDim; idim++)
  {
    _scale1[idim] /= _scadef1;
    _scale2[idim] /= _scadef2;
  }
}

ElemNostat* ModelNostat::addElemNostat()
{
  int nparam = static_cast<int> (_elems.size());
  _elems.resize(nparam + 1);
  ElemNostat* elem = new(ElemNostat);
  _elems[nparam] = elem;
  return elem;
}

void ModelNostat::define(int icov, const CovAniso* cova)
{
  int flag_range = 0;
  for (int ipar=0; ipar<(int) _elems.size(); ipar++)
  {
    ElemNostat* elem = _elems[ipar];

    if (elem->getRankStr() != icov) continue;

    if (elem->getLocType() == EConsElem::ANGLE)
    {
      _angles1[elem->getRankV1()] = elem->getVal1();
      _angles2[elem->getRankV1()] = elem->getVal2();
    }
    else if (elem->getLocType() == EConsElem::SILL)
    {
      _sill1 = elem->getVal1();
      _sill2 = elem->getVal2();
    }
    else if (elem->getLocType() == EConsElem::PARAM)
    {
      // Convert from scale to range
      if (flag_range)
        for (int idim=0; idim<_nDim; idim++)
        {
          _scale1[idim] *= _scadef1;
          _scale2[idim] *= _scadef2;
        }
      _param1  = elem->getVal1();
      _param2  = elem->getVal2();
      _scadef1 = cova_get_scale_factor(cova->getType(),_param1);
      _scadef2 = cova_get_scale_factor(cova->getType(),_param2);

      // Convert back from range to scale
      if (flag_range)
        for (int idim=0; idim<_nDim; idim++)
        {
          _scale1[idim] /= _scadef1;
          _scale2[idim] /= _scadef2;
        }
    }
    else if (elem->getLocType() == EConsElem::RANGE)
    {
      flag_range = 1;
      _scale1[elem->getRankV1()] =
          FFFF(elem->getVal1()) ? TEST : elem->getVal1() / _scadef1;
      _scale2[elem->getRankV1()] =
          FFFF(elem->getVal2()) ? TEST : elem->getVal2() / _scadef2;
    }
    else if (elem->getLocType() == EConsElem::SCALE)
    {
      _scale1[elem->getRankV1()] = elem->getVal1();
      _scale2[elem->getRankV1()] = elem->getVal2();
    }
    else
    {
      messerr("Error in the non-stationary parameters");
    }
  }

  if (! cova->getFlagAniso())
  {
    for (int idim=1; idim<_nDim; idim++)
    {
      _scale1[idim] = _scale1[0];
      _scale2[idim] = _scale2[0];
    }
  }
}

String ModelNostat::toString(const AStringFormat* strfmt) const
{
  std::stringstream sstr;
  if (_elems.size() <= 0) return sstr.str();
  sstr << "Non_stationary Model is switched ON" << std::endl;

  sstr << "List of non-stationary parameters:" << std::endl;
  for (int ipar=0; ipar<(int) _elems.size(); ipar++)
  {
    ElemNostat* elem = _elems[ipar];
    sstr << "Parameter #" << ipar+1 << std::endl;
    sstr << elem->toString(strfmt);
    sstr << std::endl;
  }
  return sstr.str();
}
