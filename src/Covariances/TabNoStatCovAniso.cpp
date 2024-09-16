#include "Covariances/TabNoStatCovAniso.hpp"
#include "Basic/AStringable.hpp"
#include "Covariances/TabNoStat.hpp"
#include "Enum/EConsElem.hpp"

TabNoStatCovAniso::TabNoStatCovAniso()
: _nAngles(0)
, _nRanges(0)
, _nScales(0)
, _nTensor(0)
, _param(false)
, _definedForAnisotropy(false)
, _definedByAnglesAndScales(false)
, _definedForRotation(false)
, _definedForTensor(false)
{}

TabNoStatCovAniso::TabNoStatCovAniso(const TabNoStatCovAniso &m):
    TabNoStat(m)
{
    this->_definedByAnglesAndScales = m._definedByAnglesAndScales;
    this->_definedForAnisotropy = m._definedForAnisotropy;
    this->_definedForRotation = m._definedForRotation;
    this->_definedForTensor = m._definedForTensor;
    this->_nTensor = m._nTensor;
    this->_nAngles = m._nAngles;
    this->_nRanges = m._nRanges;
    this->_nScales = m._nScales;
    this->_param = m._param;
}


TabNoStatCovAniso& TabNoStatCovAniso::operator= (const TabNoStatCovAniso &m)
{
    TabNoStat::operator =(m);
    if (this != &m)
    {
        _definedByAnglesAndScales = m._definedByAnglesAndScales;
        _definedForAnisotropy = m._definedForAnisotropy;
        _definedForRotation = m._definedForRotation;
        _definedForTensor = m._definedForTensor;
        _nTensor = m._nTensor;
        _nAngles = m._nAngles;
        _nRanges = m._nRanges;
        _nScales = m._nScales;
        _param = m._param;
    }
    return *this;
}

/**
 * Look if a Non-stationary parameter for Anisotropy is defined
 * either by Tensor or by (angle/range/scale)
 * @return
 */
bool TabNoStatCovAniso::isDefinedForAnisotropy() const
{
  return _definedForAnisotropy;
}

/*
 * Look if a Non-stationary parameter for Anisotropy is defined
 * only by angle / range / scale
 * @return
 */
 
bool TabNoStatCovAniso::isDefinedForRotation() const
{
    return _definedForRotation;
}

void TabNoStatCovAniso::_updateDescription()
{
  _definedForRotation = (_nAngles > 0) || (_nRanges > 0) || (_nScales > 0);
  _definedForTensor   = _nTensor > 0;
  _definedForAnisotropy = _definedForRotation || _definedForTensor;
}


int TabNoStatCovAniso::addElem(std::shared_ptr<ANoStat> &nostat,const EConsElem &econs, int iv1, int iv2) 
{
    if (econs == EConsElem::RANGE)
    {
        if (isElemDefined(EConsElem::SCALE, iv1, iv2) && _nScales == 1)
        {
            removeElem(EConsElem::SCALE,iv1,iv2);
            messerr("Warning, you gave a non-stationary specification for the range");
            messerr("but it was already given for the scale.");
            messerr("The new specification has replaced the previous one.");
        }
        else if(_nScales > 0)
        {
            messerr("You try to specify non stationarities for range whereas");
            messerr("you had already specified one for the scale in another dimension.");
            messerr("It is invalid");
            return 0;
        }          
    }
    if (econs == EConsElem::SCALE)
    {
        if (isElemDefined(EConsElem::RANGE, iv1, iv2) && _nRanges == 1)
        {
            removeElem(EConsElem::RANGE,iv1,iv2);
            messerr("Warning, you gave a non-stationary specification for the scale");
            messerr("but it was already given for the range.");
            messerr("The new specification has replaced the previous one.");
        }
        else if (_nRanges > 0) 
        {
            messerr("You try to specify non stationarities for scale whereas");
            messerr("you had already specified one for the range in another dimension.");
            messerr("It is invalid");
            return 0;
        }
    }
    int res = TabNoStat::addElem(nostat, econs, iv1,iv2);
    if (res == 0) return res;
    if (econs == EConsElem::PARAM)
        _param = true;
    if (econs == EConsElem::TENSOR)
        _nTensor += res;
    if (econs == EConsElem::RANGE) 
        _nRanges += res;
    if (econs == EConsElem::SCALE) 
        _nScales += res;
    if (econs == EConsElem::ANGLE)
        _nAngles += res;
    updateDescription();
    return res;
}

int TabNoStatCovAniso::removeElem(const EConsElem &econs, int iv1, int iv2) 
{
    int res = TabNoStat::removeElem(econs, iv1,iv2);
    if (res == 0) return res;
    if (econs == EConsElem::PARAM)
        _param = false;
    if (econs == EConsElem::TENSOR)
        _nTensor -= res;
    if (econs == EConsElem::RANGE)
        _nRanges -= res;
    if (econs == EConsElem::SCALE)
        _nScales -= res;
    if (econs == EConsElem::ANGLE)
        _nAngles -= res;
    updateDescription();
    return res;
}

bool TabNoStatCovAniso::_isValid(const EConsElem& econs) const
{
    return (econs == EConsElem::RANGE || econs == EConsElem::ANGLE ||
            econs == EConsElem::SCALE || econs == EConsElem::TENSOR|| 
            econs == EConsElem::PARAM);
}

TabNoStatCovAniso::~TabNoStatCovAniso()
{

}

