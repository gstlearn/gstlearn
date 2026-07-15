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
#include "Basic/Message.hpp"
#include "Basic/VectorHelper.hpp"
#include "Basic/VectorNumT.hpp"
#include "Basic/VectorT.hpp"
#include "DataBase/Array2D.hpp"
#include "gstlearn_export.hpp"

#include <optional>
#include <span>
#include <type_traits>
#include <utility>
#include <variant>

namespace gstlrn
{
  /**
   * @brief The DbCol class represents a single column in the DbData class:
   *        - it is a heterogeneous container that can hold data of different types
   *        (e.g., double, int, string) and supports multiple versions of the data
   *        - it may contain one or several versions of the data (e.g. outcomes for simulations)²
   *
   *        Each column is identified by its name, and it provides methods to get and set values,
   *        add samples, and retrieve the underlying data as a specific type.
   */
  class GSTLEARN_EXPORT DbCol: public ASerializable
  {
  public:
    DbCol(
      Id uniqueIndex,
      String&& name,
      const Id nsample = 0,
      const Id nversion = 1,
      bool forbidNA = false)
      : _name{std::move(name)}
      , _data{Array2D<VectorDouble>{nsample, nversion}}
      , _forbidNA{forbidNA}
      , _uniqueIndex{uniqueIndex}
    {
    }

    template<typename VectorType>
    DbCol(
      Id uniqueIndex,
      String&& name,
      VectorType&& tab,
      const Id nversion = 1,
      bool forbidNA = false);

    template<typename VectorType>
    DbCol(
      Id uniqueIndex,
      const String& name,
      VectorType&& tab,
      const Id nversion = 1,
      bool forbidNA = false);

    /// ASerializable interface
    String getNFName() const override { return "Column"; }
#ifdef HDF5
    bool deserializeH5(H5::Group& grp) override;
    bool serializeH5(H5::Group& grp) const override;
#endif

    const String& getName() const { return _name; }

    void setName(const String& newName) { this->_name = newName; }

    /**
     * @brief Get the Value of an element of the Column
     *
     * @tparam T
     * @param isample Rank of the sample
     * @param iversion Rank of the version
     * @return std::optional<T>, e.g. the value of the sample if it exists, std::nullopt otherwise
     */
    template<typename T>
    std::optional<T> getValue(const Id isample, const Id iversion = 0) const;

    /**
     * @brief Set the Value of an element of the Column
     *
     * @tparam T
     * @param v Value to be written
     * @param isample Rank of the sample
     * @param iversion Rank of the version
     */
    template<typename T>
    bool setValue(const Id isample, const Id iversion, const T& v);

    template<typename T>
    void addSamples(const Id nnewsample, const T& val);

    template<typename VectorType>
    std::optional<std::reference_wrapper<VectorType>> getVector();

    template<typename VectorType>
    std::optional<std::reference_wrapper<const VectorType>> getVector() const;

    template<typename VectorType>
    std::optional<std::span<typename VectorType::value_type>>
      getVersion(Id iversion = 0);

    template<typename VectorType>
    std::optional<std::span<const typename VectorType::value_type>>
      getVersion(Id iversion = 0) const;

    template<typename VectorType>
    bool setVector(const VectorType& values);

    template<typename VectorType>
    bool setVersion(const VectorType& values, Id iversion = 0);

    template<class F>
    void visitData(F&& f) const;

    Id getNSamples() const;

    Id getNVersions() const;

    bool forbidNA() const { return _forbidNA; }

    String getTypeName() const;

    String getDescr() const;

    void setForbidNA(bool forbidNA) { this->_forbidNA = forbidNA; }

    static DbCol createEmpty() { return DbCol(); }

    void deleteSample(const Id isample);

    Id getUniqueIndex() const { return _uniqueIndex; }

  private:
    // This line allows the deserialization of DbCol from DbData,
    // which needs to create a DbCol without knowing explicitly its type.
    DbCol() = default;

    template<typename VectorType>
    void _reset(String&& name, VectorType&& vec);

    template<typename T, typename VectorType>
    static constexpr bool _isConsistent();

    template<typename T, typename VectorType>
    static constexpr bool _isReadable();

    template<typename T, typename VectorType>
    static constexpr bool _isWritable();

    template<typename VectorType>
    Array2D<VectorType>* _getArray();

    template<typename VectorType>
    const Array2D<VectorType>* _getArray() const;

    template<typename T>
    bool _checkAllowedValue(const T& value, Id isample) const;

    template<typename VectorType>
    bool _checkAllowedValues(const VectorType& values) const;

    template<typename VectorType>
    bool _checkVectorSize(const VectorType& values) const;

    template<typename Container>
    bool _checkVersionSize(const Container& values, Id iversion) const;

    template<typename T>
    void _wrongReadType() const;

    template<typename T>
    void _wrongWriteType() const;

    bool _checkVersion(Id iversion, Id nversion) const;
    bool _checkSample(Id isample, Id nsample) const;

  private:
    String _name;
    std::variant<
      Array2D<VectorDouble>,
      Array2D<VectorFloat>,
      Array2D<VectorInt>,
      Array2D<VectorUChar>,
      Array2D<VectorString>,
      Array2D<VectorBool>
      // Array2D<VectorCategory>
      >
      _data;
    bool _forbidNA = false;
    Id _uniqueIndex = 0;
  };

  template<typename VectorType>
  DbCol::DbCol(
    Id uniqueIndex,
    String&& name,
    VectorType&& tab,
    const Id nversion,
    bool forbidNA)
    : _name{std::move(name)}
    , _data{Array2D<std::decay_t<VectorType>>(
        std::forward<VectorType>(tab),
        nversion)}
    , _forbidNA{forbidNA}
    , _uniqueIndex{uniqueIndex}
  {
  }

  template<typename VectorType>
  DbCol::DbCol(
    Id uniqueIndex,
    const String& name,
    VectorType&& tab,
    const Id nversion,
    bool forbidNA)
    : _name{name}
    , _data{Array2D<std::decay_t<VectorType>>(
        std::forward<VectorType>(tab),
        nversion)}
    , _forbidNA{forbidNA}
    , _uniqueIndex{uniqueIndex}
  {
  }

  /**
   * @brief Get the Value of an element of the Column
   *
   * @tparam T
   * @param isample Rank of the sample
   * @param iversion Rank of the version
   * @return std::optional<T>, e.g. the value of the sample if it exists, std::nullopt otherwise
   */
  template<typename T>
  std::optional<T> DbCol::getValue(const Id isample, const Id iversion) const
  {
    std::optional<T> v{};
    std::visit(
      [this, isample, iversion, &v](auto&& arg)
      {
        if (!_checkSample(isample, arg.inner())) return;

        if (!_checkVersion(iversion, arg.outer())) return;

        using VectorType = std::decay_t<decltype(arg)>::vector_type;
        if constexpr (_isReadable<T, VectorType>())
        {
          auto val = arg.getValue(iversion, isample);
          if (val) v = static_cast<T>(val.value());
        }
        else
        {
          _wrongReadType<T>();
        }
      },
      this->_data);
    return v;
  }

  template<typename T>
  bool DbCol::setValue(const Id isample, const Id iversion, const T& v)
  {
    if (!_checkAllowedValue(v, isample)) return false;

    bool res = false;
    std::visit(
      [this, isample, iversion, &v, &res](auto&& arg)
      {
        if (!_checkSample(isample, arg.inner())) return;

        if (!_checkVersion(iversion, arg.outer())) return;

        using VectorType = std::decay_t<decltype(arg)>::vector_type;
        if constexpr (_isWritable<T, VectorType>())
        {
          res = arg.setValue(iversion, isample, v);
        }
        else
        {
          _wrongWriteType<T>();
        }
      },
      this->_data);
    return res;
  }

  template<typename T>
  void DbCol::addSamples(const Id nnewsample, const T& val)
  {
    std::visit(
      [nnewsample, &val, this](auto&& arg)
      {
        using VectorType = std::decay_t<decltype(arg)>::vector_type;
        if constexpr (_isWritable<T, VectorType>())
        {
          arg.addSamples(nnewsample, val);
        }
        else
        {
          _wrongWriteType<T>();
        }
      },
      this->_data);
  }

  template<typename VectorType>
  std::optional<std::reference_wrapper<VectorType>> DbCol::getVector()
  {
    auto* pv = _getArray<VectorType>();
    if (pv == nullptr) return {};
    return pv->getBuffer();
  }

  template<typename VectorType>
  std::optional<std::reference_wrapper<const VectorType>>
    DbCol::getVector() const
  {
    auto* pv = _getArray<VectorType>();
    if (pv == nullptr) return {};
    return pv->getBuffer();
  }

  template<typename VectorType>
  std::optional<std::span<typename VectorType::value_type>>
    DbCol::getVersion(Id iversion)
  {
    auto* pv = _getArray<VectorType>();
    if (pv == nullptr) return {};

    if (!_checkVersion(iversion, pv->outer())) return {};

    auto& buffer = pv->getBuffer();

    return std::span<typename VectorType::value_type>(
      buffer.data() + iversion * pv->inner(), pv->inner());
  }

  template<typename VectorType>
  std::optional<std::span<const typename VectorType::value_type>>
    DbCol::getVersion(Id iversion) const
  {
    auto* pv = _getArray<VectorType>();
    if (pv == nullptr) return {};

    if (!_checkVersion(iversion, pv->outer())) return {};

    const auto& buffer = pv->getBuffer();

    return std::span<const typename VectorType::value_type>(
      buffer.data() + iversion * pv->inner(), pv->inner());
  }

  template<typename VectorType>
  bool DbCol::setVector(const VectorType& values)
  {
    if (!_checkVectorSize(values)) return false;

    if (!_checkAllowedValues(values)) return false;

    auto vec = getVector<VectorType>();
    if (!vec) return false;

    vec->get() = values;
    return true;
  }

  template<typename VectorType>
  bool DbCol::setVersion(const VectorType& values, Id iversion)
  {
    if (!_checkVersion(iversion, getNVersions())) return false;

    if (!_checkVersionSize(values, iversion)) return false;

    if (!_checkAllowedValues(values)) return false;

    auto span = getVersion<VectorType>(iversion);
    if (!span) return false;

    std::copy(values.begin(), values.end(), span->begin());

    return true;
  }

  template<class F>
  void DbCol::visitData(F&& f) const
  {
    std::visit(std::forward<F>(f), _data);
  }

  template<typename VectorType>
  void DbCol::_reset(String&& name, VectorType&& vec)
  {
    _name = std::move(name);
    _data = Array2D<std::decay_t<VectorType>>(std::forward<VectorType>(vec));
  }

  template<typename T, typename VectorType>
  constexpr bool DbCol::_isConsistent()
  {
    return std::is_same_v<T, typename VectorType::value_type>
        || (std::is_same_v<T, bool> && std::is_same_v<VectorType, VectorBool>)
        || std::is_convertible_v<T, typename VectorType::value_type>;
  }

  template<typename T, typename VectorType>
  constexpr bool DbCol::_isReadable()
  {
    return std::is_convertible_v<typename VectorType::value_type, T>;
  }

  template<typename T, typename VectorType>
  constexpr bool DbCol::_isWritable()
  {
    using ValueType = typename VectorType::value_type;

    return std::is_same_v<T, ValueType>
        || (std::is_arithmetic_v<T> && std::is_arithmetic_v<ValueType>)
        || (std::is_same_v<T, bool> && std::is_same_v<VectorType, VectorBool>);
  }

  template<typename VectorType>
  Array2D<VectorType>* DbCol::_getArray()
  {
    auto* pv = std::get_if<Array2D<VectorType>>(&_data);

    if (pv == nullptr)
    {
      messerr(
        "Returned type is not available as Column '%s' is of type '%s'.",
        getName().c_str(), getTypeName().c_str());
    }
    return pv;
  }

  template<typename VectorType>
  const Array2D<VectorType>* DbCol::_getArray() const
  {
    auto* pv = std::get_if<Array2D<VectorType>>(&_data);

    if (pv == nullptr)
    {
      messerr(
        "Returned type is not available as Column '%s' is of type '%s'.",
        getName().c_str(), getTypeName().c_str());
    }
    return pv;
  }

  template<typename T>
  bool DbCol::_checkAllowedValue(const T& value, Id isample) const
  {
    if (this->forbidNA() && isNA<T>(value))
    {
      messerr(
        "Cannot set NA to sample %d of Column '%s': NA values are forbidden.",
        isample, this->getName().c_str());
      return false;
    }
    return true;
  }

  template<typename VectorType>
  bool DbCol::_checkAllowedValues(const VectorType& values) const
  {
    if (!_forbidNA) return true;

    for (const auto& v: values)
    {
      if (isNA<typename VectorType::value_type>(v))
      {
        messerr(
          "Cannot set NA in Column '%s': NA values are forbidden.",
          getName().c_str());
        return false;
      }
    }
    return true;
  }

  template<typename VectorType>
  bool DbCol::_checkVectorSize(const VectorType& values) const
  {
    const auto expected = static_cast<size_t>(getNSamples()) * getNVersions();
    const auto actual = values.size();

    if (actual != expected)
    {
      messerr(
        "Input vector has size %d while Column '%s' has size %d.",
        static_cast<Id>(actual), getName().c_str(), static_cast<Id>(expected));
      return false;
    }
    return true;
  }

  template<typename Container>
  bool DbCol::_checkVersionSize(const Container& values, Id iversion) const
  {
    const auto expected = static_cast<size_t>(getNSamples());
    const auto actual = values.size();

    if (actual != expected)
    {
      messerr(
        "Input vector has size %d while version %d of Column '%s' has size "
        "%d.",
        static_cast<Id>(actual), static_cast<int>(iversion),
        this->getName().c_str(), static_cast<Id>(expected));
      return false;
    }
    return true;
  }

  template<typename T>
  void DbCol::_wrongReadType() const
  {
    messerr("Cannot read value from Column '%s':", this->getName().c_str());
    messerr(
      "Requested type '%s' is not compatible with Stored type '%s'.",
      getGenericTypeName<T>().c_str(), this->getTypeName().c_str());
  }

  template<typename T>
  void DbCol::_wrongWriteType() const
  {
    messerr("Cannot assign value to Column '%s':", this->getName().c_str());
    messerr(
      "Input type '%s' is not compatible with column type '%s'.",
      getGenericTypeName<T>().c_str(), this->getTypeName().c_str());
  }

} // namespace gstlrn
