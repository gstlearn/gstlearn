/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2026) MINES Paris / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://gstlearn.org                                              */
/* License: BSD 3-clause                                                      */
/*                                                                            */
/******************************************************************************/
#pragma once

#include "Basic/ASerializable.hpp"
#include "Basic/AStringable.hpp"
#include "DataBase/Category.hpp"
#include "DataBase/Dictionary.hpp"
#include "gstlearn_export.hpp"

#include <cstddef>
#include <optional>
#include <ostream>
#include <vector>

namespace gstlrn
{
  /**
   * @brief Vector of categorical values associated with a dictionary.
   *
   * A VectorCategory stores one optional category for each sample.
   * Categories are represented by Category and are associated
   * with a Dictionary defining the available category identifiers and
   * their labels.
   *
   * A sample may have no category, in which case its value is represented
   * by an empty optional.
   */
  class GSTLEARN_EXPORT VectorCategory: public AStringable, public ASerializable
  {
  public:
    using value_type = Category;

    /**
     * @brief Constructs a categorical vector.
     *
     * @param nsample Number of samples.
     * @param dict Dictionary associated with the categorical values.
     */
    VectorCategory(
      const size_t nsample = 0,
      const Dictionary& dict = Dictionary())
      : _data(nsample)
      , _dict{dict}
    {
    }

    /// AStringable Interface
    String toString(const AStringFormat* strfmt = nullptr) const override;

    /// ASerializable interface
    String getNFName() const override { return "VectorCategory"; }
#ifdef HDF5
    bool deserializeH5(H5::Group& grp) override;
    bool serializeH5(H5::Group& grp) const override;
#endif

    inline void resize(const size_t count) { this->_data.resize(count); }

    inline size_t size() const { return this->_data.size(); }

    std::optional<Category> getCategory(const size_t isample) const;

    bool setCategory(const size_t isample, const Category& cat);
    bool setCategory(const size_t isample, const Id key);
    bool setCategoriesByKey(const std::vector<Id>& keys);

    Id getNCategories() const { return this->_dict.getNCategories(); }

    const Dictionary& getDictionary() const { return this->_dict; }

#ifndef SWIG
    /**
     * @brief Provides access to a sample category.
     *
     * @param isample Sample index.
     *
     * @return Reference to the optional category associated with the sample.
     *
     * @warning No bounds checking is performed.
     */
    std::optional<Category>& operator[](const size_t isample)
    {
      return this->_data[isample];
    }

    /**
     * @brief Provides read-only access to a sample category.
     *
     * @param isample Sample index.
     *
     * @return Constant reference to the optional category associated with
     *         the sample.
     *
     * @warning No bounds checking is performed.
     */
    const std::optional<Category>& operator[](const size_t isample) const
    {
      return this->_data[isample];
    }
#endif

  private:
    bool _checkSampleIndex(const size_t isample) const;

  private:
    std::vector<std::optional<Category>> _data;
    Dictionary _dict;
  };

#ifndef SWIG
  GSTLEARN_EXPORT std::ostream&
    operator<<(std::ostream& os, const VectorCategory& vec);
#endif

} // namespace gstlrn
