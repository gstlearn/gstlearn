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
#include "DataBase/DbData.hpp"
#include "Basic/SerializeHDF5.hpp"
#include <boost/math/special_functions/detail/fp_traits.hpp>

namespace gstlrn
{
  DbData* DbData::createFromNF(const String& NFFilename, bool verbose)
  {
    auto* dbdata = new DbData;
    if (dbdata->_fileOpenAndDeserialize(NFFilename, verbose)) return dbdata;
    delete dbdata;
    return nullptr;
  }

#ifdef HDF5
  bool DbData::serializeH5(H5::Group& grp) const
  {
    auto dbG = grp.createGroup("DbData");
    SerializeHDF5::writeValue(dbG, "NColumn", getNColumns());

    for (Id i = 0; i < getNColumns(); i++)
    {
      auto colG = dbG.createGroup("Column_" + std::to_string(i));

      // Information contained in DbData
      SerializeHDF5::writeValue(colG, "Role", _roleIDs[i].getRole().getValue());
      SerializeHDF5::writeValue(colG, "Rank", _roleIDs[i].getIndex());

      // Information contained in DbCol
      if (!_cols[i].serializeH5(colG)) return false;
    }

    return true;
  }

  bool DbData::deserializeH5(H5::Group& grp)
  {
    auto dbG = grp.openGroup("DbData");

    Id ncols;
    Id roleValue;
    Id rank;

    // Read number of columns
    SerializeHDF5::readValue(dbG, "NColumn", ncols);

    // Clear previous contents of the Data Base
    deleteAllColumns();

    for (Id i = 0; i < ncols; i++)
    {
      auto colG = dbG.openGroup("Column_" + std::to_string(i));

      // Read DbData information
      SerializeHDF5::readValue(colG, "Role", roleValue);
      SerializeHDF5::readValue(colG, "Rank", rank);

      // Read DbCol information
      _cols.push_back(DbCol::createEmpty());

      if (!_cols.back().deserializeH5(colG)) return false;

      _roleIDs.emplace_back(ERole::fromValue(roleValue), rank);
    }

    return true;
  }
#endif

  bool DbData::renameColumn(ColID&& colid, const String& newName)
  {
    const auto icol = _getColumnIndex(colid);
    if (!icol) return false;
    this->_cols[*icol].setName(newName);
    return true;
  }

  void DbData::deleteColumn(ColID&& colid)
  {
    const auto icol = _getColumnIndex(colid);
    if (!icol) return;

    auto roleID = _roleIDs[*icol];
    if (static_cast<Id>(this->_cols.size()) > (*icol))
    {
      this->_cols.erase(this->_cols.begin() + (*icol));
      this->_roleIDs.erase(this->_roleIDs.begin() + (*icol));
    }

    // Update the RoleID list after deletion
    _updateRoleIDDeletion(roleID);
  }

  void DbData::deleteAllColumns()
  {
    this->_cols.clear();
    this->_roleIDs.clear();
  }

  /**
   * @brief Identify and returns a reference on the Column of interest
   *
   * @param colid Column Indentifier
   * @return std::optional<std::reference_wrapper<DbCol>>
   */
  std::optional<std::reference_wrapper<DbCol>>
    DbData::_identifyColumn(ColID&& colid)
  {
    const auto icol = _getColumnIndex(colid);
    if (!icol) return std::nullopt;
    return {this->_cols[*icol]};
  }

  std::optional<std::reference_wrapper<const DbCol>>
    DbData::_identifyColumn(ColID&& colid) const
  {
    const auto icol = _getColumnIndex(colid);
    if (!icol) return std::nullopt;
    return {this->_cols[*icol]};
  }

  bool DbData::hasColumn(ColID&& colid) const
  {
    const auto icol = _getColumnIndex(colid);
    return static_cast<bool>(icol);
  }

  /**
   * @brief Get the Name object
   *
   * @param colid Column Indentifier
   * @return String
   */
  String DbData::getName(ColID&& colid) const
  {
    const auto icol = _getColumnIndex(colid);
    return icol ? _cols[*icol].getName() : String();
  }

  Id DbData::getUniqueIndex(ColID&& colid) const
  {
    const auto icol = _getColumnIndex(colid);
    return icol ? _cols[*icol].getUniqueIndex() : -1;
  }

  /**
   * @brief Get the Index of the Column
   *
   * @param colid Column indentifier
   * @return Id
   */
  Id DbData::getICol(ColID&& colid) const
  {
    const auto icol = _getColumnIndex(colid, false);
    return icol ? *icol : -1;
  }

  /**
   * @brief Get the RoleID of the Column
   *
   * @param colid Column indentifier
   * @return RoleID
   */
  RoleID DbData::getRoleID(ColID&& colid) const
  {
    const auto icol = _getColumnIndex(colid);
    return icol ? _roleIDs[*icol] : RoleID();
  }

  const ERole& DbData::getRole(ColID&& colid) const
  {
    const auto icol = _getColumnIndex(colid);
    return icol ? _roleIDs[*icol].getRole() : ERole::UNDEFINED;
  }

  ColID DbData::getColID(const ColID& colid) const
  {
    const auto icol = _getColumnIndex(colid, false);
    if (!icol) return ColID(-1);

    ColID colIDout(*icol);
    colIDout.setName(_cols[*icol].getName());
    colIDout.setRoleID(_roleIDs[*icol]);
    return colIDout;
  }

  void DbData::setRoleID(ColID&& colid, const RoleID& roleID)
  {
    const auto icol = _getColumnIndex(colid);
    if (!icol) return;
    if (roleID.isEqual(_roleIDs[*icol], true)) return;

    auto roleIDLocal = roleID;
    _updateRoleIDAddition(*icol, roleIDLocal);
    _roleIDs[*icol] = roleIDLocal;
  }

  void DbData::setName(ColID&& colid, const String& newName)
  {
    const auto icol = _getColumnIndex(colid);
    if (!icol) return;
    if (newName == _cols[*icol].getName()) return;

    auto proposedName = newName;
    _updateName(proposedName);
    _cols[*icol].setName(proposedName);
  }

  /**
   * @brief Return a set of Column Identifiers starting from a Column Name
   *
   * @param name Column Name (with regexp facility)
   * @return std::vector<ColID>
   */
  std::vector<ColID> DbData::getColIDs(const String& name) const
  {
    // Looking for matching names
    VectorString matchNames = expandList(getNames(), name);

    // Loop to create the list of Column Identifiers
    std::vector<ColID> colIDs;
    for (const auto& matchName: matchNames)
    {
      const auto colID = getColID(ColID(matchName));
      if (colID.getICol() >= 0) colIDs.push_back(colID);
    }
    return colIDs;
  }

  /**
   * @brief Return a set of Column Identifiers starting from a set of Column Names
   *
   * @param names Column Names (with regexp facility)
   * @return std::vector<ColID>
   */
  std::vector<ColID> DbData::getColIDs(const VectorString& names) const
  {
    VectorString matchNames = expandList(getNames(), names);
    std::vector<ColID> colIDs;
    for (const auto& matchName: matchNames)
    {
      const auto colID = getColID(ColID(matchName));
      if (colID.getICol() >= 0) colIDs.push_back(colID);
    }
    return colIDs;
  }

  /**
   * @brief Returns a set of Column Identifiers matching a given Role
   *
   * @param role Target Role to be searched
   * @return std::vector<ColID>
   */
  std::vector<ColID> DbData::getColIDs(const ERole& role) const
  {
    std::vector<ColID> colIDs;
    for (const auto& id: this->_roleIDs)
    {
      if (id.getRole().isEqual(role))
      {
        const auto colID = getColID(ColID(id));
        if (colID.getICol() >= 0) colIDs.push_back(colID);
      }
    }
    return colIDs;
  }

  Id DbData::getNSamples() const
  {
    if (getNColumns() <= 0) return 0;
    return _cols[0].getNSamples();
  }

  void DbData::removeRole(ColID&& colid)
  {
    const auto icol = _getColumnIndex(colid);
    if (!icol) return;
    _roleIDs[*icol].removeRole();
  }

  void DbData::removeAllRoles()
  {
    for (auto& roleID: _roleIDs)
    {
      roleID.removeRole();
    }
  }

  Id DbData::getNVersions(ColID&& colid) const
  {
    const auto icol = _getColumnIndex(colid);
    if (!icol) return 0;
    return _cols[*icol].getNVersions();
  }

  /**
   * @brief Count the number of entries with a given Role
   *
   * @param role Target Role to be matched
   */
  Id DbData::getNRoles(const ERole& role) const
  {
    Id count = 0;
    for (const auto& roleID: _roleIDs)
      if (roleID.getRole().isEqual(role)) ++count;
    return count;
  }

  Id DbData::getNRoles(ColID&& colid) const
  {
    const auto icol = _getColumnIndex(colid, false);
    if (!icol) return 0;

    const auto& role = _roleIDs[*icol].getRole();
    return getNRoles(role);
  }

  /**
   * @brief Erase the Role of all Columns with the given Role
   *
   * @param role Target Role to be erased
   */
  void DbData::clearRole(const ERole& role)
  {
    std::vector<ColID> colIDs = getColIDs(role);
    for (const auto& colID: colIDs) _roleIDs[colID.getICol()].removeRole();
  }

  void DbData::clearAllRoles()
  {
    for (auto& roleID: _roleIDs) roleID.removeRole();
  }

  /**
   * @brief Produces the contents of the DbData content
   */
  void DbData::printContents(const String& title) const
  {
    Id ncols = _cols.size();

    if (!title.empty()) std::cout << title << '\n';
    std::cout << "The Data Base contains " << ncols << " column(s) of "
              << getNSamples() << " samples\n";
    for (Id icol = 0; icol < ncols; ++icol)
    {
      const auto& c = _cols[icol];
      const auto& id = _roleIDs[icol];

      std::cout << "Column " << icol << "/" << ncols;
      std::cout << " : " << c.getDescr();
      std::cout << " {Role: " << id.getName() << "}";
      if (c.forbidNA()) std::cout << " [NA forbidden]";
      std::cout << std::endl;
    }
  }

  /**
   * @brief Returns the Column index corresponding to the given Column Indentifier (ColID)
   *
   * @param colid Column Indentifier
   * @param verbose Whether to print error messages if the Column is not found (default = true)
   * @return std::optional<Id>
   *
   * @remark: The Column Indentifier (ColID) is searched in the following order:
   * - by Column Name
   * - by Column Role (and Rank)
   * - by Column index
   */
  std::optional<Id>
    DbData::_getColumnIndex(const ColID& colid, bool verbose) const
  {
    auto ncol = getNColumns();

    // Try to identify by Column Name
    const String& localName = colid.getName();
    if (!localName.empty())
    {
      for (Id icol = 0; icol < ncol; ++icol)
      {
        // 'localName' may contain a regular expression, so we use matchRegexp to check for a match
        if (matchRegexp(_cols[icol].getName(), localName)) return icol;
        //  if (_cols[icol].getName() == localName) return icol;
      }
      if (verbose) _unknownName(localName);
      return std::nullopt;
    }

    // Try to identify by Column Role
    if (colid.getRole().isDifferent(ERole::UNDEFINED))
    {
      const RoleID& roleID = colid.getRoleID();
      for (Id icol = 0; icol < ncol; ++icol)
      {
        if (_roleIDs[icol].isEqual(roleID, true)) return icol;
      }
      if (verbose) _unknownRoleID(roleID);
      return std::nullopt;
    }

    // Try to identify by Column rank
    if (colid.getICol() >= 0) return colid.getICol();
    if (verbose) messerr("Column does not exist.");
    return std::nullopt;
  }

  /**
   * @brief Check if the new name is compatible with existing ones
   *
   * @param name Name of the new column to be added (possibly modified)
   *
   * @remark If the new name matches an already existing one, it is modified to be unique
   */
  void DbData::_updateName(String& name) const
  {
    // Establish the list of already existing names
    VectorString proposedNames = getNames();
    auto ncol = static_cast<Id>(proposedNames.size());

    // Add the new proposal to the list of already existing names
    proposedNames.push_back(name);

    // Modify the 'ncol' proposal and retrive the modified value
    correctNamesForDuplicates(proposedNames, ncol);
    name = proposedNames[ncol];
  }

  /**
   * @brief Check the addition of a new column with the roleID provided as argument
   *
   * @param icol0 Index of the column to be added (<0 for a non existing column)
   * @param roleID RoleID of the new column to be added (possibly modified)
   *
   * @remark If the Role of the new Column is already present in the already defined ones:
   * - if the Rank of the new Column matches the one of the old matching Column:
   *   the Old matching Column is moved to an UNDEFINED Role (and a Rank set to 0).
   *   the New Column keeps its Role and Rank unchanged
   * - if the Rank of the new Column does not match the one of the old matching Column,
   *   this rank is calculated as the largest Rank found in matching Columns incremented by 1.
   */
  void DbData::_updateRoleIDAddition(Id icol0, RoleID& roleID)
  {
    auto ncol = getNColumns();
    auto wasDefined = icol0 >= 0 && _roleIDs[icol0].isDefined();
    const auto newRole = roleID.getRole();
    const auto newIndex = roleID.getIndex();

    // Check if the target Column does not already have the same roleID (ERole and Index)
    if (icol0 >= 0 && roleID.isEqual(_roleIDs[icol0], true)) return;

    // The ERole of the new Column is defined.
    if (roleID.isDefined())
    {
      Id rankMin = -1;

      // Look for a Column with the same roleID (same ERole and same Index)
      for (Id icol = 0; icol < ncol; icol++)
      {
        if (icol == icol0) continue;

        // Look for a column (different from the new one) with the same RoleID
        if (roleID.isEqual(_roleIDs[icol], true))
        {
          // Set it to UNDEFINED
          _roleIDs[icol] = RoleID(ERole::UNDEFINED, 0);
        }

        // Calculate the highest index of the same RoleID already present in the DataBase
        if (roleID.isEqual(_roleIDs[icol], false))
        {
          auto curIndex = _roleIDs[icol].getIndex();
          if (curIndex > rankMin) rankMin = curIndex;
        }
      }

      // Modify the Index of the new Column (if necessary)
      if (rankMin >= 0 && newIndex > rankMin) roleID.setIndex(rankMin + 1);
    }
    else
    {
      // The new Column exists: it is undefined, but was defined previously
      if (wasDefined)
      {
        // Look for all Columns with same ERole (but possibly different Rank)
        for (Id icol = 0; icol < ncol; icol++)
        {
          if (icol == icol0) continue;
          if (!_roleIDs[icol0].isEqual(_roleIDs[icol], false)) continue;

          auto oldIndex = _roleIDs[icol0].getIndex();
          auto curIndex = _roleIDs[icol].getIndex();
          if (curIndex > oldIndex)
          {
            // Same role and a larger rank already exists: decrease the index of the old Column by one
            _roleIDs[icol].setIndex(curIndex - 1);
          }
        }
      }
    }
  }

  void DbData::_updateRoleIDDeletion(RoleID& roleID)
  {
    const auto& newRole = roleID.getRole();
    const auto newIndex = roleID.getIndex();

    // Look for already existing Columns with the same RoleID
    for (auto& id: this->_roleIDs)
    {
      if (id.getRole().isEqual(newRole) && id.getIndex() > newIndex)
      {
        // Same role and a larger rank already exists:
        // Decrease the rank of the old one by 1
        id.setIndex(id.getIndex() - 1);
      }
    }
  }

  /**
   * @brief Check if the new Category Column is Valid when forbidNA is true
   *
   * @param tab VectorCategory of the new column to be added
   */
  bool DbData::_checkForbidNA(const VectorCategory& tab)
  {
    for (size_t i = 0; i < tab.size(); i++)
    {
      if (!tab[i].has_value())
      {
        messerr("Column forbids NA values, but the input tab contains some.");
        return false;
      }
    }
    return true;
  }

  /**
   * @brief Get the Names of all the Columns
   *
   * @return VectorString
   */
  VectorString DbData::getNames() const
  {
    VectorString names;
    for (const auto& col: _cols)
    {
      names.push_back(col.getName());
    }
    return names;
  }

  void DbData::addSamples(Id nadd, const double valinit)
  {
    if (nadd <= 0) return;
    for (Id icol = 0, ncol = getNColumns(); icol < ncol; icol++)
    {
      _cols[icol].addSamples(nadd, valinit);
    }
  }

  void DbData::deleteSample(Id idel)
  {
    for (Id icol = 0, ncol = getNColumns(); icol < ncol; icol++)
    {
      _cols[icol].deleteSample(idel);
    }
  }

  Id DbData::getColMatchUniqueIndex(Id uniqueIndex) const
  {
    // This trick is meant to save time as this function systematically searches
    // for a match over the 'ncol' columns.
    // If two consecutive calls are made with the same uniqueIndex, the second call
    // will start the loop at the position of the previous search
    static Id lastUniqueIndex = 0;

    const Id ncol = getNColumns();

    for (Id icol = 0; icol < ncol; icol++)
    {
      const Id jcol = (icol + lastUniqueIndex) % ncol;

      if (_cols[jcol].getUniqueIndex() == uniqueIndex)
      {
        lastUniqueIndex = jcol; // Save this for the next search
        return jcol;
      }
    }

    return -1;
  }

  void DbData::_checkVersion(Id& nversion)
  {
    if (nversion <= 0)
    {
      messerr(
        "The number of versions (%d) must be strictly positive. It has "
        "been set to 1.",
        nversion);
      nversion = 1;
    }
  }

  String DbData::summaryRoles(void) const
  {
    std::stringstream sstr;

    auto it = ERole::getIterator();
    while (it.hasNext())
    {
      const auto& role = *it;

      auto colIDs = getColIDs(role);
      if (!colIDs.empty())
      {
        sstr << "Role: " << role.getKey() << " (" << getNRoles(role) << ")";
        sstr << " - Columns = ";
        for (const auto& colID: colIDs) sstr << colID.getICol() << " ";
        sstr << std::endl;
      }
      it.toNext();
    }
    return sstr.str();
  }

  void DbData::_unknownName(const String& name)
  {
    messerr("Column '%s' does not exist.", name.c_str());
  }

  void DbData::_unknownRoleID(const RoleID& roleID)
  {
    messerr("Role '%s' does not exist.", roleID.getName().c_str());
  }

  /************************************************************************/
  /* Internal helper classes for column access syntax                     */
  /* Only used for aliasing and slicing operations, not for data storage. */
  /************************************************************************/
  DbData::ColProxy DbData::X(Id rank)
  {
    return ColProxy(*this, ColID(ERole::X, rank));
  }

  DbData::ColProxy DbData::Z(Id rank)
  {
    return ColProxy(*this, ColID(ERole::Z, rank));
  }

  DbData::ColProxy DbData::W(Id rank)
  {
    return ColProxy(*this, ColID(ERole::W, rank));
  }

  DbData::ColProxy DbData::F(Id rank)
  {
    return ColProxy(*this, ColID(ERole::F, rank));
  }

  DbData::ColProxy DbData::SEL(Id rank)
  {
    return ColProxy(*this, ColID(ERole::SEL, rank));
  }

  DbData::ColProxy DbData::role(const ERole& role, Id rank)
  {
    return ColProxy(*this, ColID(role, rank));
  }

} // namespace gstlrn
