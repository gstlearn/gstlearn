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
#include "Basic/AStringable.hpp"
#include "DataBase/Category.hpp"
#include "geoslib_define.h"
#include "gstlearn_export.hpp"

#include <map>
#include <optional>
#include <vector>

namespace gstlrn
{
  /**
   * @brief Dictionary associating category identifiers with their labels.
   *
   * A dictionary stores a set of categories identified by an integer key.
   * Each category is associated with a textual label.
   */
  class GSTLEARN_EXPORT Dictionary: public ASerializable, public AStringable
  {
  public:
    Dictionary() = default;

    /**
     * @brief Constructs a dictionary from an existing map.
     *
     * The map is moved into the dictionary.
     *
     * @param data Map associating category identifiers with their labels.
     */
    Dictionary(std::map<Id, String>&& data)
      : _data{std::move(data)}
    {
    }

    /// AStringable Interface
    String toString(const AStringFormat* strfmt = nullptr) const override;

    /// ASerializable interface
    String getNFName() const override { return "Dictionary"; }
#ifdef HDF5
    bool deserializeH5(H5::Group& grp) override;
    bool serializeH5(H5::Group& grp) const override;
#endif

    bool addCategory(const Id key, const String& val);

    bool hasCategory(const Category& cat) const;

    std::optional<Category> getCategory(const Id key) const;

    std::vector<Category> getCategories() const;

    Id getNCategories() const { return static_cast<Id>(this->_data.size()); }

#ifndef SWIG
    /**
     * @brief Retrieves a category by its identifier.
     *
     * @param key Category identifier.
     *
     * @return The corresponding category, or an empty optional if the
     *         identifier is not present in the dictionary.
     */
    std::optional<Category> operator[](const Id key) const
    {
      return getCategory(key);
    }
#endif

    static Id getCategoryKey(const Category& cat) { return cat.getId(); }

    static String getCategoryName(const Category& cat) { return cat.getName(); }

  private:
    std::map<Id, String> _data;
  };

} // namespace gstlrn
