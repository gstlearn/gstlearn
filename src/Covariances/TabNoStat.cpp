#include "Covariances/TabNoStat.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/VectorNumT.hpp"
#include "Covariances/ParamId.hpp"
#include "Enum/EConsElem.hpp"
#include "geoslib_define.h"
#include <memory>
#include "Db/Db.hpp"

TabNoStat::TabNoStat()
  : _items()
  , _dbNoStatRef(nullptr)
{
}

TabNoStat::TabNoStat(const TabNoStat& m)
  : AStringable(m)
{
  _dbNoStatRef = m._dbNoStatRef;
  _items       = m._items;
}

TabNoStat& TabNoStat::operator=(const TabNoStat& m)
{
  if (this != &m)
  {
    _dbNoStatRef = m._dbNoStatRef;
    _items       = m._items;
  }
  return *this;
}

int TabNoStat::removeElem(const EConsElem& econs, int iv1, int iv2)
{
  ParamId param(econs, iv1, iv2);
  int res = _items.erase(param);
  updateDescription();
  return res;
}

void TabNoStat::clear()
{
  _items.clear();
  _clear();
}
bool TabNoStat::isValid(const EConsElem& econs) const
{
  bool res = _isValid(econs);
  if (!res)
  {
    messerr("%s is an invalid parameter for this covariance structure",
            std::string(econs.getKey()).c_str());
  }
  return res;
}

void TabNoStat::updateDescription()
{
  _updateDescription();
}

bool TabNoStat::isElemDefined(const EConsElem& econs, int iv1, int iv2) const
{
  ParamId conselem(econs, iv1, iv2);
  return _items.count(conselem) > 0; // Warning : use count for C++17 compatibility
}

std::shared_ptr<ANoStat> TabNoStat::getElem(const EConsElem& econs, int iv1, int iv2)
{
  ParamId conselem(econs, iv1, iv2);
  return _items[conselem];
}

String TabNoStat::toString(const AStringFormat* strfmt) const
{
  return toStringInside(strfmt, 0);
}
String TabNoStat::toStringInside(const AStringFormat* strfmt, int i) const
{
  std::stringstream sstr;
  if (_items.empty()) return sstr.str();

  for (const auto& e: getTable())
  {
    sstr << std::to_string(i + 1) << " - ";
    sstr << e.first.toString(strfmt);
    sstr << e.second->toString(strfmt);
    i++;
  }
  return sstr.str();
}

void TabNoStat::informCoords(const VectorVectorDouble& coords,
                             const EConsElem& econs,
                             int iv1,
                             int iv2,
                             VectorDouble& result) const
{
  ParamId conselem(econs, iv1, iv2);
  if (isElemDefined(econs, iv1, iv2))
    _items.at(conselem)->informField(coords, result);
}

int TabNoStat::addElem(std::shared_ptr<ANoStat>& nostat, const EConsElem& econs, int iv1, int iv2)
{
  if (!isValid(econs))
    return 0;
  ParamId param(econs, iv1, iv2);
  int res       = _items.count(param);
  _items[param] = nostat;
  if (res == 1)
  {
    messerr("Warning, this non stationarity was already specified. It has been replaced");
    messerr("with the new specifications.");
  }
  res = 1 - res;
  updateDescription();
  return res;
}

void TabNoStat::setDbNoStatRef(const Db* dbref)
{
  if (dbref != nullptr)
    _dbNoStatRef = std::shared_ptr<const Db>((Db*)dbref->clone());
  else
    _dbNoStatRef = nullptr;
}

void TabNoStat::setDbNoStatRef(std::shared_ptr<const Db>& dbref)
{
  _dbNoStatRef = dbref;
}

void TabNoStat::informMeshByMesh(const AMesh* amesh) const
{
  for (const auto& e: _items)
  {
    e.second->informMeshByMesh(amesh);
  }
}
void TabNoStat::informMeshByApex(const AMesh* amesh) const
{
  for (const auto& e: _items)
  {
    e.second->informMeshByApex(amesh);
  }
}
void TabNoStat::informDbIn(const Db* dbin) const
{
  for (const auto& e: _items)
  {
    e.second->informDbIn(dbin);
  }
}
void TabNoStat::informDbOut(const Db* dbout) const
{
  for (const auto& e: _items)
  {
    e.second->informDbOout(dbout);
  }
}

std::shared_ptr<const Db> TabNoStat::getDbNoStatRef() const
{
  return _dbNoStatRef;
}

const Db* TabNoStat::getDbNoStatRefRaw() const
{
  return _dbNoStatRef.get();
}

void TabNoStat::informMeshByMesh(const AMesh* amesh, const EConsElem& econs) const
{
  for (const auto& e: _items)
  {
    if (e.first.getType() == econs)
      e.second->informMeshByMesh(amesh);
  }
}

void TabNoStat::informMeshByApex(const AMesh* amesh, const EConsElem& econs) const
{
  for (const auto& e: _items)
  {
    if (e.first.getType() == econs)
      e.second->informMeshByApex(amesh);
  }
}
void TabNoStat::informDbIn(const Db* dbin, const EConsElem& econs) const
{
  for (const auto& e: _items)
  {
    if (e.first.getType() == econs)
      e.second->informDbIn(dbin);
  }
}
void TabNoStat::informDbOut(const Db* dbout, const EConsElem& econs) const
{
  for (const auto& e: _items)
  {
    if (e.first.getType() == econs)
      e.second->informDbOout(dbout);
  }
}

TabNoStat::~TabNoStat()
{
}