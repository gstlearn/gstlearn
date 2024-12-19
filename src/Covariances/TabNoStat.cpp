#include "Covariances/TabNoStat.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/String.hpp"
#include "Basic/VectorNumT.hpp"
#include "Covariances/ParamId.hpp"
#include "Enum/EConsElem.hpp"
#include "geoslib_define.h"
#include <memory>
#include "Db/Db.hpp"

TabNoStat::TabNoStat()
  : _items()
  , _dbNoStatRef(nullptr)
  , _definedForVariance(false)
  , _nSills(0)
{
}

TabNoStat::TabNoStat(const TabNoStat& m)
  : AStringable(m)
{
  this->_nSills             = m._nSills;
  this->_dbNoStatRef        = m._dbNoStatRef;
  this->_definedForVariance = m._definedForVariance;
  this->_items              = m._items;
}

TabNoStat& TabNoStat::operator=(const TabNoStat& m)
{
  if (this != &m)
  {
    _nSills             = m._nSills;
    _dbNoStatRef        = m._dbNoStatRef;
    _definedForVariance = m._definedForVariance;
    _items              = m._items;
  }
  return *this;
}

int TabNoStat::removeElem(const EConsElem &econs, int iv1, int iv2)
{
    ParamId param(econs,iv1,iv2);
    int res = _items.erase(param);
    if (econs == EConsElem::SILL)
        _nSills -= res;
    updateDescription();
    return res;

}

bool TabNoStat::isValid(const EConsElem& econs) const
{
    bool res = _isValid(econs);
    if (!res)
    {
        messerr("Invalid type of parameters for this covariance structure");
    }
    return res ;
}
bool TabNoStat::_isValid(const EConsElem& econs) const
{
     return (econs == EConsElem::SILL);
}

void TabNoStat::updateDescription()
{
    _definedForVariance = _nSills > 0;
    _updateDescription();
}

bool TabNoStat::isElemDefined(const EConsElem &econs, int iv1, int iv2) const
{
    ParamId conselem(econs,iv1,iv2);
    return _items.count(conselem) > 0;
}

std::shared_ptr<ANoStat> TabNoStat::getElem(const EConsElem &econs, int iv1, int iv2)
{
    ParamId conselem(econs,iv1,iv2);
    return _items[conselem];
}

String TabNoStat::toString(const AStringFormat* strfmt) const
{
    return toStringInside(strfmt,0);
}
String TabNoStat::toStringInside(const AStringFormat* strfmt,int i) const
{
    std::stringstream sstr;
    if (_items.empty()) return sstr.str();
    
    for (const auto &e: getTable())
    {
        sstr << std::to_string(i+1) << " - ";
        sstr << e.first.toString(strfmt);
        sstr << e.second->toString(strfmt);
        i++;
    }
    return sstr.str();
}

void TabNoStat::informCoords(const VectorVectorDouble &coords,
                             const EConsElem &econs,
                             int iv1, 
                             int iv2, 
                             VectorDouble& result) const
{
    ParamId conselem(econs,iv1,iv2);
    if (isElemDefined(econs, iv1,iv2))
        _items.at(conselem)->informField(coords,result);
}

int TabNoStat::addElem(std::shared_ptr<ANoStat> &nostat, const EConsElem &econs, int iv1, int iv2)
{
    if (!isValid(econs))
        return 0;
    ParamId param(econs,iv1,iv2);
    int res = _items.count(param);
    _items[param] = nostat;
    if (res == 1)
    {
        messerr("Warning, this non stationarity was already specified. It has been replaced");
        messerr("with the new specifications.");
    }
    res = 1 - res;
    if (econs == EConsElem::SILL)
        _nSills += res;
    updateDescription();
    return res;
}

void TabNoStat::informMeshByMesh(const AMesh* amesh) const
{
    for (const auto &e : _items)
    {
        e.second->informMeshByMesh(amesh);
    }
}
void TabNoStat::informMeshByApex(const AMesh* amesh) const
{
    for (const auto &e : _items)
    {
        e.second->informMeshByApex(amesh);
    }

}
void TabNoStat::informDbIn(const Db* dbin) const
{
    for (const auto &e : _items)
    {
        e.second->informDbIn(dbin);
    }
}
void TabNoStat::informDbOut(const Db* dbout) const 
{
    for (const auto &e : _items)
    {
        e.second->informDbOout(dbout);
    }
}

void TabNoStat::informMeshByMesh(const AMesh* amesh, const EConsElem & econs) const
{
    for (const auto &e : _items)
    {
        if (e.first.getType() == econs)
            e.second->informMeshByMesh(amesh);
    }
}

void TabNoStat::informMeshByApex(const AMesh* amesh, const EConsElem & econs) const
{
    for (const auto &e : _items)
    {
        if (e.first.getType() == econs)
            e.second->informMeshByApex(amesh);
    }

}
void TabNoStat::informDbIn(const Db* dbin, const EConsElem & econs) const
{
    for (const auto &e : _items)
    {
        if (e.first.getType() == econs)
            e.second->informDbIn(dbin);
    }
}
void TabNoStat::informDbOut(const Db* dbout, const EConsElem & econs) const
{
    for (const auto &e : _items)
    {
        if (e.first.getType() == econs)
            e.second->informDbOout(dbout);
    }
}

TabNoStat::~TabNoStat()
{

}