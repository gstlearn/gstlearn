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
#include "Covariances/TabNoStat.hpp"
#include "Basic/ASerializable.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/SerializeHDF5.hpp"
#include "Basic/VectorNumT.hpp"
#include "Covariances/NoStatArray.hpp"
#include "Covariances/ParamId.hpp"
#include "Db/Db.hpp"
#include "Db/DbGrid.hpp"
#include "Enum/EConsElem.hpp"
#include "geoslib_define.h"
#include <memory>

namespace gstlrn
{

TabNoStat::TabNoStat()
  : AStringable()
  , ASerializable()
  , _items()
  , _dbNoStatRef(nullptr)
{
}

TabNoStat::TabNoStat(const TabNoStat& m)
  : AStringable(m)
  , ASerializable(m)
{
  _dbNoStatRef = m._dbNoStatRef;
  _items       = m._items;
}

TabNoStat& TabNoStat::operator=(const TabNoStat& m)
{
  if (this != &m)
  {
    AStringable::operator=(m);
    ASerializable::operator=(m);
    _dbNoStatRef = m._dbNoStatRef;
    _items       = m._items;
  }
  return *this;
}

Id TabNoStat::removeElem(const EConsElem& econs, Id iv1, Id iv2)
{
  ParamId param(econs, iv1, iv2);
  Id res = static_cast<Id>(_items.erase(param));
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

bool TabNoStat::isElemDefined(const EConsElem& econs, Id iv1, Id iv2) const
{
  ParamId conselem(econs, iv1, iv2);
#ifdef USE_BOOST_SPAN
  return _items.count(conselem) > 0;
#else
  return _items.contains(conselem);
#endif
}

std::shared_ptr<ANoStat> TabNoStat::getElem(const EConsElem& econs, Id iv1, Id iv2)
{
  ParamId conselem(econs, iv1, iv2);
  return _items[conselem];
}

String TabNoStat::toString(const AStringFormat* strfmt) const
{
  return toStringInside(strfmt, 0);
}

String TabNoStat::toStringInside(const AStringFormat* strfmt, Id i) const
{
  std::stringstream sstr;
  if (_items.empty()) return sstr.str();

  // sort content of std::unordered_map before serialization
  std::vector<std::pair<std::string, std::string>> sortedTable(getTable().size());
  size_t j {};
  for (const auto& e: getTable())
  {
    sortedTable[j++] = {e.first.toString(strfmt), e.second->toString(strfmt)};
  }
  std::sort(sortedTable.begin(), sortedTable.end());
  for (const auto& e: sortedTable)
  {
    sstr << std::to_string(i + 1) << " - ";
    sstr << e.first;
    sstr << e.second;
    i++;
  }
  return sstr.str();
}

void TabNoStat::informCoords(const VectorVectorDouble& coords,
                             const EConsElem& econs,
                             Id iv1,
                             Id iv2,
                             VectorDouble& result) const
{
  ParamId conselem(econs, iv1, iv2);
  if (isElemDefined(econs, iv1, iv2))
    _items.at(conselem)->informField(coords, result);
}

Id TabNoStat::addElem(std::shared_ptr<ANoStat>& nostat, const EConsElem& econs, Id iv1, Id iv2)
{
  if (!isValid(econs))
    return 0;
  ParamId param(econs, iv1, iv2);
  Id res        = static_cast<Id>(_items.count(param));
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
  {
    _dbNoStatRef = std::shared_ptr<const Db>(dbref->clone());
  }
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

bool TabNoStat::variableExistsInDb(const String& namecol) const
{
  if (_dbNoStatRef->getUID(namecol) < 0)
  {
    messerr("The Name of the Non-stationary variable (%s) does not exist in the Reference Db",
            namecol.c_str());
    return false;
  }
  return true;
}

#ifdef HDF5
bool TabNoStat::deserializeH5(H5::Group& grp)
{
  bool ret = true;
  Id nrank = 0;
  ret      = ret && SerializeHDF5::readValue(grp, "Number of items", nrank);

  // Explicit loop to read the non-stationary parameters
  Id type = 0;
  Id iv1  = 0;
  Id iv2  = 0;
  String colName;
  bool isGrid = false;
  for (Id rank = 0; rank < nrank; rank++)
  {
    String locName = "Param_" + std::to_string(rank + 1);
    auto itemG     = SerializeHDF5::getGroup(grp, locName, false);
    if (!itemG) return true;

    ret = ret && SerializeHDF5::readValue(*itemG, "Type", type);
    ret = ret && SerializeHDF5::readValue(*itemG, "IV1", iv1);
    ret = ret && SerializeHDF5::readValue(*itemG, "IV2", iv2);
    ret = ret && SerializeHDF5::readValue(*itemG, "ColName", colName);
    ret = ret && SerializeHDF5::readValue(*itemG, "dbIsGrid", isGrid);

    const EConsElem& econs = EConsElem::fromValue(type);
    if (isGrid)
    {
      auto* dbgrid = new DbGrid();
      ret          = ret && dbgrid->deserializeH5(*itemG);
      setDbNoStatRef(dbgrid);
      if (!variableExistsInDb(colName)) return false;
      delete dbgrid;
    }
    else
    {
      auto* db = new Db();
      ret      = ret && db->deserializeH5(*itemG);
      setDbNoStatRef(db);
      if (!variableExistsInDb(colName)) return false;
      delete db;
    }

    std::shared_ptr<ANoStat> ns;
    ns = std::shared_ptr<ANoStat>(new NoStatArray(getDbNoStatRef(), colName));
    addElem(ns, econs, iv1, iv2);
  }

  return ret;
}

bool TabNoStat::serializeH5(H5::Group& grp) const
{
  bool ret = true;
  ret      = ret && SerializeHDF5::writeValue(grp, "Number of items", size());

  // Implicit loop to write the non-stationary parameters (in NoStatArray style)
  Id rank = 0;
  for (const auto& [paramId, noStatPtr]: _items)
  {
    rank++;
    auto& sna = dynamic_cast<NoStatArray&>(*noStatPtr);

    const auto& dbnostat = sna.getDbNoStat();
    bool isGrid          = dbnostat->isGrid();

    auto itemG = grp.createGroup("Param_" + std::to_string(rank));
    ret        = ret && SerializeHDF5::writeValue(itemG, "Type", paramId.getType().getValue());
    ret        = ret && SerializeHDF5::writeValue(itemG, "IV1", paramId.getIV1());
    ret        = ret && SerializeHDF5::writeValue(itemG, "IV2", paramId.getIV2());
    ret        = ret && SerializeHDF5::writeValue(itemG, "ColName", sna.getColName());
    ret        = ret && SerializeHDF5::writeValue(itemG, "dbIsGrid", isGrid);
    ret        = ret && dbnostat->serializeH5(itemG);
  }
  return ret;
}
#endif

} // namespace gstlrn
