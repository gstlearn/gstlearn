/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2025) MINES Paris / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://gstlearn.org                                              */
/* License: BSD 3-clause                                                      */
/*                                                                            */
/******************************************************************************/
#pragma once

#include "Basic/ASerializable.hpp"
#include "DataBase/ColID.hpp"
#include "DataBase/DbCol.hpp"
#include "DataBase/Dictionary.hpp"
#include "DataBase/MultID.hpp"
#include "DataBase/RoleID.hpp"
#include "DataBase/VectorCategory.hpp"
#include "gstlearn_export.hpp"

#include <functional>
#include <map>
#include <optional>
#include <type_traits>

namespace gstlrn
{
  /**
   * @brief Tag structure and constant used to represent full dimension selection '*'
   *        in slicing operations.
   *        Only used for aliasing and slicing operations, not for data storage.
   */
  struct AllType
  {
  };

#ifndef SWIG
  constexpr AllType ALL{};
  constexpr AllType _{};
#endif

  class GSTLEARN_EXPORT DbData: public ASerializable
  {
  public:
#ifndef SWIG
    class ColProxy;
    class ValueProxy;
    class SliceProxy;
#endif

    static DbData* createFromNF(const String& NFFilename, bool verbose);

    /// ASerializable interface
    String getNFName() const override { return "DbData"; }
#ifdef HDF5
    bool deserializeH5(H5::Group& grp) override;
    bool serializeH5(H5::Group& grp) const override;
#endif

    template<typename VectorType>
    void addColumnEmpty(
      const String& name,
      Id nsamples = 0,
      Id nversion = 1,
      const RoleID& roleID = RoleID(),
      std::optional<typename VectorType::value_type> valinit = std::nullopt,
      bool forbidNA = false,
      const Dictionary* dict = nullptr);

    template<typename VectorType>
    void addColumn(
      String&& name,
      VectorType&& array,
      const RoleID& roleID = RoleID(),
      Id nversion = 1,
      bool forbidNA = false);

    template<typename VectorType>
    void addColumn(
      const String& name,
      VectorType&& array,
      const RoleID& roleID = RoleID(),
      Id nversion = 1,
      bool forbidNA = false);

#ifndef SWIG
    template<typename T>
    std::optional<T>
      getValue(ColID&& colid, Id isample, bool verbose = true) const;

    template<typename T>
    bool setValue(ColID&& colid, const Id isample, const T& value);

    template<typename VectorType>
    VectorType& getColumn(ColID&& colid);

    template<typename VectorType>
    const VectorType& getColumn(ColID&& colid) const;

    template<typename VectorType>
    span<typename VectorType::value_type>
      getVersion(ColID&& colid, Id version = 0);

    template<typename VectorType>
    span<const typename VectorType::value_type>
      getVersion(ColID&& colid, Id version = 0) const;

    template<typename VectorType>
    bool setColumn(ColID&& colid, const VectorType& values);

    template<typename VectorType>
    bool setVersion(ColID&& colid, const VectorType& values, Id version = 0);
#endif

#ifndef SWIG
    ColProxy X(Id rank = 0);
    ColProxy Z(Id rank = 0);
    ColProxy W(Id rank = 0);
    ColProxy F(Id rank = 0);
    ColProxy SEL(Id rank = 0);
    ColProxy role(const ERole& role, Id rank = 0);
#endif

    bool renameColumn(ColID&& colid, const String& newName);
    void deleteColumn(ColID&& colid);
    void deleteAllColumns();
    bool hasColumn(ColID&& colid) const;
    bool hasValidVersion(ColID&& colid) const;
    String getName(ColID&& colid, bool withVersion = false) const;
    VectorString getNames() const;
    Id getColumnIndex(ColID&& colid) const;
    Id getUniqueIndex(ColID&& colid) const;
    RoleID getRoleID(ColID&& colid) const;
    const ERole& getRole(ColID&& colid) const;
    ColID getColID(const ColID& colid) const;
    void removeRole(ColID&& colid);
    void removeAllRoles();
    ColID fromUID(Id iuid, Id version = 0) const;
    MultID fromUIDs(const VectorInt& iuids) const;
    ColIDs getColIDs(const String& name) const;
    ColIDs getColIDs(const VectorString& name) const;
    ColIDs getColIDs(const ERole& role) const;

    Id getNColumns() const { return static_cast<Id>(_cols.size()); }

    Id getNVersions(ColID&& colid) const;
    Id getNSamples() const;
    Id getNRoles(const ERole& role) const;
    Id getNRoles(ColID&& colid) const;
    void setName(ColID&& colid, const String& newName);
    void setRoleID(ColID&& colid, const RoleID& roleID);
    void printContents(const String& title = "") const;
    void clearRole(const ERole& role);
    void clearAllRoles();
    void addSamples(Id nadd, const double valinit);
    void deleteSample(Id idel);
    void deleteVersion(ColID&& colid, Id version);

    Id getUniqueIndexCounter() const { return _uniqueIndexCounter; }

    String summaryRoles(void) const;

    bool isValidColumn(Id columnIndex) const;
    ColIDs getColumnIndices(const MultID& multid) const;

  private:
    std::optional<std::reference_wrapper<DbCol>> _identifyColumn(ColID&& colid);
    std::optional<std::reference_wrapper<const DbCol>>
      _identifyColumn(ColID&& colid, bool verbose = false) const;
    std::optional<Id>
      _getColumnIndex(const ColID& colid, bool verbose = true) const;

    template<class VectorType>
    VectorType _createEmptyVector(
      Id n,
      std::optional<typename VectorType::value_type> value,
      const Dictionary* dict = nullptr);

    template<typename VectorType>
    static bool _checkForbidNA(const VectorType& tab);

    template<typename VectorType>
    bool _checkNSample(const VectorType& array, Id nversion) const;

    void _updateName(String& name) const;
    void _updateRoleIDAddition(Id columnIndex0, RoleID& roleID);
    void _updateRoleIDDeletion(RoleID& roleID);
    void _rebuildRoleMap();
    void _rebuildUIDMap();

    Id _searchColumnByUID(Id iuid) const;

    static bool _checkForbidNA(const VectorCategory& tab);
    static void _checkNVersion(Id& nversion);
    static void _unknownName(const String& name);
    static void _unknownRoleID(const ColID& colID);
    static void _unknownColumnIndex(Id columnIndex, Id ncol);

    ColID _makeColIDFromIndex(Id icol) const;

    bool _checkVersion(const ColID& colid, Id columnIndex0) const;
    Id _checkUID(const ColID& colid, Id columnIndex0) const;
    Id _searchByRole(const RoleID& roleid, bool verbose) const;
    Id _searchByName(const String& name, bool verbose) const;
    Id _searchByUID(Id uniqueIndex) const;

  private:
    // Private section of DbData
    std::vector<DbCol> _cols;
    std::vector<RoleID> _roleIDs;
    std::map<RoleID, Id> _roleToColMap;
    std::map<Id, Id> _uidToColMap;
    Id _uniqueIndexCounter = 0;
  };

  /***************************************************************************/
  /* Internal Proxy Implementation                                           */
  /* Only used for aliasing and slicing operations, not for data storage.    */
  /***************************************************************************/
#ifndef SWIG
  /**
   * @brief Type trait to detect if a type is an iterable vector container.
   */
  template<typename T, typename = void>
  struct is_vector_like: std::false_type
  {
  };

  template<typename T>
  struct is_vector_like<
    T,
    std::void_t<typename T::value_type, decltype(std::declval<T>().size())>>
    : std::true_type
  {
  };

  template<typename T>
  inline constexpr bool is_vector_like_v = is_vector_like<T>::value;

  /**
   * @brief Proxy helper class for scalar read/write operations targeting an individual sample.
   */
  class DbData::ValueProxy
  {
  public:
    ValueProxy(DbData& db, const ColID& colid, Id isample)
      : _db(db)
      , _colid(colid)
      , _isample(isample)
    {
    }

    template<typename T>
    ValueProxy& operator=(const T& value)
    {
      _db.setValue<T>(ColID(_colid), _isample, value);
      return *this;
    }

    operator double() const { return _get<double>(); }

    operator int() const { return _get<int>(); }

    operator Id() const { return _get<Id>(); }

    operator String() const { return _get<String>(); }

  private:
    // Private section of ValueProxy
    template<typename T>
    T _get() const
    {
      auto val = _db.getValue<T>(ColID(_colid), _isample);
      return val ? *val : getNA<T>();
    }

    DbData& _db;
    ColID _colid;
    Id _isample;
  };

  /**
   * @brief Proxy class representing a multi-element slice of a DbData column.
   *        Supports reading into VectorDouble/VectorType and direct assignment/writing back into DbData.
   */
  class GSTLEARN_EXPORT DbData::SliceProxy
  {
  public:
    enum class SliceMode
    {
      ALL_SAMPLES_SINGLE_VERSION,
      SINGLE_SAMPLE_ALL_VERSIONS,
      FULL_MATRIX
    };

    SliceProxy(
      DbData& db,
      const ColID& colid,
      SliceMode mode,
      Id sampleOrVersion)
      : _db(db)
      , _colid(colid)
      , _mode(mode)
      , _targetIdx(sampleOrVersion)
    {
    }

    /**
     * @brief Build and extract the actual data slice as a vector.
     */
    template<typename VectorType = VectorDouble>
    VectorType toVector() const
    {
      if (_mode == SliceMode::ALL_SAMPLES_SINGLE_VERSION)
      {
        span<const double> span =
          _db.getVersion<VectorType>(ColID(_colid), _targetIdx);
        return VectorType(span.begin(), span.end());
      }
      if (_mode == SliceMode::SINGLE_SAMPLE_ALL_VERSIONS)
      {
        Id nversions = _db.getNVersions(ColID(_colid));
        VectorType result(nversions);
        ColID colid(_colid);
        for (Id v = 0; v < nversions; ++v)
        {
          auto val =
            _db.getValue<double>(std::move(colid.setVersion(v)), _targetIdx);
          result[v] = val ? *val : getNA<double>();
        }
        return result;
      }
      // FULL_MATRIX
      Id nversions = _db.getNVersions(ColID(_colid));
      Id nsamples = _db.getNSamples();
      VectorType result(nsamples * nversions);

      Id offset = 0;
      const double naVal = getNA<double>();
      ColID colIDref(_colid);

      for (Id v = 0; v < nversions; ++v)
      {
        colIDref.setVersion(v);
        for (Id i = 0; i < nsamples; ++i)
        {
          auto val = _db.getValue<double>(ColID(colIDref), i);
          result[offset++] = val ? *val : naVal;
        }
      }

      return result;
    }

    /**
     * @brief Implicit conversions to VectorType (e.g. VectorDouble).
     */
    template<typename VectorType = VectorDouble>
    operator VectorType() const
    {
      return toVector<VectorType>();
    }

    /**
     * @brief Direct pass-through method for .dump()
     */
    void dump(const String& title = "") const
    {
      toVector<VectorDouble>().dump(title);
    }

    // ------------------------------------------------------------------------
    // Write Access Operators (Modifies the DbData)
    // ------------------------------------------------------------------------

    /**
     * @brief Assign a scalar value to every element in the slice.
     *        Enables: data.F()(_, 3) = 13;
     *                 data.F()(2, _) = 14;
     *                 data.F()() = 1234.;
     */
    template<typename T>
    std::enable_if_t<!is_vector_like_v<T>, SliceProxy&>
      operator=(const T& value)
    {
      if (_mode == SliceMode::ALL_SAMPLES_SINGLE_VERSION)
      {
        Id nsamples = _db.getNSamples();
        ColID colIDtarget = ColID(_colid).setVersion(_targetIdx);
        for (Id i = 0; i < nsamples; ++i)
        {
          _db.setValue<T>(ColID(colIDtarget), i, value);
        }
      }
      else if (_mode == SliceMode::SINGLE_SAMPLE_ALL_VERSIONS)
      {
        Id nversions = _db.getNVersions(ColID(_colid));
        ColID colid(_colid);
        for (Id v = 0; v < nversions; ++v)
        {
          _db.setValue<T>(std::move(colid.setVersion(v)), _targetIdx, value);
        }
      }
      else // FULL_MATRIX
      {
        Id nversions = _db.getNVersions(ColID(_colid));
        Id nsamples = _db.getNSamples();
        ColID colIDver(_colid);
        for (Id v = 0; v < nversions; ++v)
        {
          colIDver.setVersion(v);
          for (Id i = 0; i < nsamples; ++i)
          {
            _db.setValue<T>(ColID(colIDver), i, value);
          }
        }
      }
      return *this;
    }

    /**
     * @brief Assign a vector of values to the slice.
     *        Enables: data.F()(_, 1) = myVector;
     */
    template<typename VectorType>
    std::enable_if_t<is_vector_like_v<VectorType>, SliceProxy&>
      operator=(const VectorType& values)
    {
      if (_mode == SliceMode::ALL_SAMPLES_SINGLE_VERSION)
      {
        _db.setVersion<VectorType>(ColID(_colid), values, _targetIdx);
      }
      else if (_mode == SliceMode::SINGLE_SAMPLE_ALL_VERSIONS)
      {
        Id nversions = _db.getNVersions(ColID(_colid));
        const Id n = std::min(nversions, static_cast<Id>(values.size()));
        for (Id v = 0; v < n; ++v)
        {
          _db.setValue(ColID(_colid).setVersion(v), _targetIdx, values[v]);
        }
      }
      else // FULL_MATRIX
      {
        _db.setColumn<VectorType>(ColID(_colid), values);
      }
      return *this;
    }

  private:
    DbData& _db;
    ColID _colid;
    SliceMode _mode;
    Id _targetIdx;
  };

  /**
   * @brief Column and Data Slicing Proxy (DbData).
   */
  class GSTLEARN_EXPORT DbData::ColProxy
  {
  public:
    ColProxy(DbData& db, const ColID& colid)
      : _db(db)
      , _colid(colid)
    {
      ERole role = _colid.getRole();

      if (role != ERole::UNDEFINED)
      {
        auto it = ERoleAttr.find(role.getKey());
        bool isUnique = (it != ERoleAttr.end()) ? it->second.isUnique : false;

        if (isUnique && _colid.getIndex() > 0)
        {
          message(
            "Warning: Role '%s' is unique. Rank %d ignored, using rank 0.\n",
            std::string(role.getKey()).c_str(), _colid.getIndex());

          _colid = ColID(role, 0);
        }
      }
    }

    // ------------------------------------------------------------------------
    // 1. Scalar Access (Single sample/version -> ValueProxy)
    // ------------------------------------------------------------------------

    ValueProxy operator[](Id isample) { return operator()(isample, 0); }

    ValueProxy operator[](Id isample) const { return operator()(isample, 0); }

    ValueProxy operator()(Id isample)
    {
      ColID copy = _colid;
      copy.setVersion(0);
      return ValueProxy(_db, copy, isample);
    }

    ValueProxy operator()(Id isample) const
    {
      return ValueProxy(_db, ColID(_colid).setVersion(0), isample);
    }

    ValueProxy operator()(Id isample, Id version)
    {
      return ValueProxy(_db, ColID(_colid).setVersion(version), isample);
    }

    ValueProxy operator()(Id isample, Id version) const
    {
      return ValueProxy(_db, ColID(_colid).setVersion(version), isample);
    }

    // ------------------------------------------------------------------------
    // 2. Slice Access via Expanders (ALL / _) returning SliceProxy
    // ------------------------------------------------------------------------

    SliceProxy operator()(AllType /*all*/, Id version)
    {
      return SliceProxy(
        _db, ColID(_colid).setVersion(version),
        SliceProxy::SliceMode::ALL_SAMPLES_SINGLE_VERSION, version);
    }

    SliceProxy operator()(AllType /*all*/, Id version) const
    {
      return SliceProxy(
        const_cast<DbData&>(_db), _colid,
        SliceProxy::SliceMode::ALL_SAMPLES_SINGLE_VERSION, version);
    }

    SliceProxy operator()(Id isample, AllType /*all*/)
    {
      return SliceProxy(
        _db, _colid, SliceProxy::SliceMode::SINGLE_SAMPLE_ALL_VERSIONS,
        isample);
    }

    SliceProxy operator()(Id isample, AllType /*all*/) const
    {
      return SliceProxy(
        const_cast<DbData&>(_db), _colid,
        SliceProxy::SliceMode::SINGLE_SAMPLE_ALL_VERSIONS, isample);
    }

    SliceProxy operator()(AllType /*all1*/ = ALL, AllType /*all2*/ = ALL)
    {
      return SliceProxy(_db, _colid, SliceProxy::SliceMode::FULL_MATRIX, 0);
    }

    SliceProxy operator()(AllType /*all1*/ = ALL, AllType /*all2*/ = ALL) const
    {
      return SliceProxy(
        const_cast<DbData&>(_db), _colid, SliceProxy::SliceMode::FULL_MATRIX,
        0);
    }

    // ------------------------------------------------------------------------
    // 3. Direct Version Setting
    // ------------------------------------------------------------------------

    template<typename VectorType>
    void setVersion(Id version, const VectorType& values)
    {
      _db.setVersion<VectorType>(ColID(_colid), values, version);
    }

    // ------------------------------------------------------------------------
    // 4. Full Column Vector & Scalar Assignment Operators
    // ------------------------------------------------------------------------

    template<typename VectorType>
    std::enable_if_t<is_vector_like_v<VectorType>, ColProxy&>
      operator=(const VectorType& values)
    {
      _db.setColumn<VectorType>(ColID(_colid), values);
      return *this;
    }

    template<typename T>
    std::enable_if_t<!is_vector_like_v<T>, ColProxy&> operator=(const T& value)
    {
      Id nsamples = _db.getNSamples();
      Id nversions = _db.getNVersions(ColID(_colid));
      for (Id v = 0; v < nversions; ++v)
      {
        ColID colIDver = ColID(_colid).setVersion(v);
        for (Id i = 0; i < nsamples; ++i)
        {
          _db.setValue<T>(ColID(colIDver), i, value);
        }
      }
      return *this;
    }

    // ------------------------------------------------------------------------
    // 5. Implicit Conversions
    // ------------------------------------------------------------------------

    template<typename VectorType = VectorDouble>
    operator VectorType&()
    {
      return _db.getColumn<VectorType>(ColID(_colid));
    }

    template<typename VectorType = VectorDouble>
    operator const VectorType&() const
    {
      return _db.getColumn<VectorType>(ColID(_colid));
    }

  private:
    // Private section of SliceProxy
    DbData& _db;
    ColID _colid;
  };

#endif

  template<typename VectorType>
  void DbData::addColumnEmpty(
    const String& name,
    Id nsamples,
    Id nversion,
    const RoleID& roleID,
    std::optional<typename VectorType::value_type> valinit,
    bool forbidNA,
    const Dictionary* dict)
  {
    if (getNColumns() > 0) nsamples = getNSamples();
    if (nsamples <= 0)
    {
      messerr("The number of samples (%d) must be positive.", nsamples);
      return;
    }
    auto array =
      _createEmptyVector<VectorType>(nsamples * nversion, valinit, dict);
    addColumn(name, std::move(array), roleID, nversion, forbidNA);
  }

  template<typename VectorType>
  void DbData::addColumn(
    String&& name,
    VectorType&& array,
    const RoleID& roleID,
    Id nversion,
    bool forbidNA)
  {
    _checkNVersion(nversion);
    if (!_checkNSample(array, nversion)) return;

    auto nameLocal = name;
    for (const auto& col: _cols)
    {
      if (col.getName() == name)
      {
        _updateName(nameLocal);
        break;
      }
    }

    auto roleIDLocal = roleID;
    _updateRoleIDAddition(-1, roleIDLocal);

    if (forbidNA)
    {
      if (!_checkForbidNA(array)) return;
    }

    _cols.emplace_back(
      _uniqueIndexCounter, std::move(nameLocal),
      std::forward<VectorType>(array), nversion, forbidNA);

    _roleIDs.emplace_back(roleIDLocal);

    // Increment the unique index counter for each new column added
    _uniqueIndexCounter++;

    _rebuildRoleMap();
    _rebuildUIDMap();
  }

  template<typename VectorType>
  void DbData::addColumn(
    const String& name,
    VectorType&& array,
    const RoleID& roleID,
    Id nversion,
    bool forbidNA)
  {
    _checkNVersion(nversion);
    if (!_checkNSample(array, nversion)) return;

    auto nameLocal = name;
    for (const auto& col: _cols)
    {
      if (col.getName() == name)
      {
        _updateName(nameLocal);
        break;
      }
    }

    auto roleIDLocal = roleID;
    _updateRoleIDAddition(-1, roleIDLocal);

    if (forbidNA)
    {
      if (!_checkForbidNA(array)) return;
    }

    this->_cols.emplace_back(
      _uniqueIndexCounter, std::move(nameLocal),
      std::forward<VectorType>(array), nversion, forbidNA);

    _roleIDs.emplace_back(roleIDLocal);

    // Increment the unique index counter for each new column added
    _uniqueIndexCounter++;

    _rebuildRoleMap();
    _rebuildUIDMap();
  }

#ifndef SWIG
  template<typename T>
  std::optional<T>
    DbData::getValue(ColID&& colid, Id isample, bool verbose) const
  {
    const auto array = _identifyColumn(std::move(colid), verbose);
    if (!array) return std::nullopt;

    const Id version = colid.getVersion();
    return array->get().getValue<T>(isample, version);
  }

  template<typename T>
  bool DbData::setValue(ColID&& colid, const Id isample, const T& value)
  {
    const auto array = _identifyColumn(std::move(colid));
    if (!array) return false;

    const Id version = colid.getVersion();
    return array->get().setValue<T>(isample, version, value);
  }

  template<typename VectorType>
  VectorType& DbData::getColumn(ColID&& colid)
  {
    static VectorType empty{};

    auto col = _identifyColumn(std::move(colid));
    if (!col) return empty;

    auto vec = col->get().template getVector<VectorType>();
    if (!vec) return empty;

    return vec->get();
  }

  template<typename VectorType>
  const VectorType& DbData::getColumn(ColID&& colid) const
  {
    static const VectorType empty{};

    auto col = _identifyColumn(std::move(colid));
    if (!col) return empty;

    auto vec = col->get().template getVector<VectorType>();
    if (!vec) return empty;

    return vec->get();
  }

  template<typename VectorType>
  span<typename VectorType::value_type>
    DbData::getVersion(ColID&& colid, Id version)
  {
    auto col = _identifyColumn(std::move(colid));
    if (!col) return {};

    auto span = col->get().template getVersion<VectorType>(version);
    if (!span) return {};

    return *span;
  }

  template<typename VectorType>
  span<const typename VectorType::value_type>
    DbData::getVersion(ColID&& colid, Id version) const
  {
    auto col = _identifyColumn(std::move(colid));
    if (!col) return {};

    auto span = col->get().template getVersion<VectorType>(version);
    if (!span) return {};

    return *span;
  }

  template<typename VectorType>
  bool DbData::setColumn(ColID&& colid, const VectorType& values)
  {
    auto col = _identifyColumn(std::move(colid));
    if (!col) return false;

    return col->get().template setVector<VectorType>(values);
  }

  template<typename VectorType>
  bool DbData::setVersion(ColID&& colid, const VectorType& values, Id version)
  {
    auto col = _identifyColumn(std::move(colid));
    if (!col) return false;

    return col->get().template setVersion<VectorType>(values, version);
  }

  template<class VectorType>
  VectorType DbData::_createEmptyVector(
    Id n,
    std::optional<typename VectorType::value_type> value,
    const Dictionary* dict)
  {
    if constexpr (std::is_same_v<VectorType, VectorCategory>)
    {
      if (dict == nullptr)
        throw std::invalid_argument("Dictionary is required");

      VectorCategory vec(n, *dict);
      if (value)
      {
        for (Id i = 0; i < n; i++) vec[i] = *value;
      }
      return vec;
    }
    else
    {
      DECLARE_UNUSED(dict)
      const auto actual =
        value.value_or(getNA<typename VectorType::value_type>());

      return VectorType(n, actual);
    }
  }

  template<typename VectorType>
  bool DbData::_checkForbidNA(const VectorType& tab)
  {
    using ValueType = typename VectorType::value_type;
    for (const auto& val: tab)
    {
      if (isNA<ValueType>(val))
      {
        messerr("Column forbids NA values, but the input tab contains some.");
        return false;
      }
    }
    return true;
  }

  template<typename VectorType>
  bool DbData::_checkNSample(const VectorType& array, Id nversion) const
  {
    const Id nsample = getNSamples();
    if (nsample <= 0) return true;

    const Id size = static_cast<Id>(array.size());
    const Id expected = nsample * nversion;

    if (size == expected) return true;

    messerr(
      "The number of values (%lld) is not compatible with "
      "the number of samples (%lld) and versions (%lld).",
      size, nsample, nversion);

    return false;
  }

#endif
} // namespace gstlrn
