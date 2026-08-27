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
#include "DataBase/VectorCategory.hpp"
#include "Basic/SerializeHDF5.hpp"

namespace gstlrn
{
  /**
   * @brief Sets the category associated with a sample.
   *
   * The operation fails if the sample index is out of range or if the
   * specified category is not defined in the associated dictionary.
   *
   * @param isample Sample index.
   * @param cat Category to assign to the sample.
   *
   * @return @c true if the category was successfully assigned,
   *         @c false otherwise.
   */
  bool VectorCategory::setCategory(const size_t isample, const Category& cat)
  {
    if (!_checkSampleIndex(isample)) return false;

    if (!this->_dict.hasCategory(cat)) return false;

    this->_data[isample] = cat;
    return true;
  }

  /**
   * @brief Sets the category associated with a sample using its identifier.
   *
   * The category is searched in the associated dictionary using its
   * identifier. The operation fails if the sample index is out of range
   * or if no category with the specified identifier is defined in the
   * dictionary.
   *
   * @param isample Sample index.
   * @param key Category identifier.
   *
   * @return @c true if the category was successfully assigned,
   *         @c false otherwise.
   */
  bool VectorCategory::setCategory(const size_t isample, const Id key)
  {
    if (!_checkSampleIndex(isample)) return false;

    const auto cat = this->_dict.getCategory(key);
    if (!cat) return false;

    this->_data[isample] = *cat;
    return true;
  }

  bool VectorCategory::setCategoriesByKey(const std::vector<Id>& keys)
  {
    this->_data.resize(keys.size());
    for (size_t i = 0; i < keys.size(); ++i)
    {
      const Id key = keys[i];

      if (key < 0)
      {
        this->_data[i].reset();
        continue;
      }

      const auto cat = this->_dict.getCategory(key);
      if (!cat) return false;

      this->_data[i] = *cat;
    }

    return true;
  }

  /**
   * @brief Returns the category associated with a sample.
   *
   * @param isample Sample index.
   *
   * @return The category associated with the sample, or an empty optional
   *         if the sample has no category or if the index is out of range.
   */
  std::optional<Category>
    VectorCategory::getCategory(const size_t isample) const
  {
    if (!_checkSampleIndex(isample)) return {};
    return this->_data[isample];
  }

#ifdef HDF5

  bool VectorCategory::serializeH5(H5::Group& grp) const
  {
    /*
     * Dictionary
     */
    H5::Group dictGrp = grp.createGroup("Dictionary");

    if (!this->_dict.serializeH5(dictGrp)) return false;

    /*
     * Values
     *
     * Each value contains the key of the corresponding category
     * in the dictionary.
     *
     * -1 means that the sample has no category (NA).
     */
    std::vector<Id> values(this->_data.size());

    for (size_t i = 0; i < this->_data.size(); i++)
    {
      const auto& cat = this->_data[i];

      if (!cat)
      {
        values[i] = -1;
        continue;
      }

      /*
       * Get the category index in the dictionary.
       *
       * getCategoryKey() returns -1 if the category is not found.
       */
      values[i] = cat->getId();
    }

    /*
     * Store category indices
     */
    hsize_t dims[1] = {static_cast<hsize_t>(values.size())};
    H5::DataSpace space(1, dims);

    const auto& type = SerializeHDF5::getHDF5Type<Id>();

    auto ds = grp.createDataSet("Values", type, space);

    if (!values.empty()) ds.write(values.data(), type);

    return true;
  }

  bool VectorCategory::deserializeH5(H5::Group& grp)
  {
    /*
     * Dictionary
     */
    H5::Group dictGrp = grp.openGroup("Dictionary");

    Dictionary dict;

    if (!dict.deserializeH5(dictGrp)) return false;

    this->_dict = std::move(dict);

    /*
     * Values
     */
    auto ds = grp.openDataSet("Values");

    H5::DataSpace space = ds.getSpace();

    hsize_t dims[1];
    space.getSimpleExtentDims(dims);

    const size_t nsample = dims[0];

    std::vector<Id> values(nsample);

    const auto& type = SerializeHDF5::getHDF5Type<Id>();

    if (!values.empty()) ds.read(values.data(), type);

    /*
     * Rebuild sample values.
     */
    this->_data.resize(nsample);

    for (size_t i = 0; i < nsample; i++)
    {
      const Id categoryKey = values[i];

      /*
       * -1 means NA.
       */
      if (categoryKey < 0)
      {
        this->_data[i].reset();
        continue;
      }

      const auto cat = this->_dict.getCategory(categoryKey);

      if (!cat) return false;

      this->_data[i] = *cat;
    }

    return true;
  }

#endif

  bool VectorCategory::_checkSampleIndex(const size_t isample) const
  {
    if (isample >= this->_data.size())
    {
      messerr("Sample index (%d) is out of bounds.", static_cast<Id>(isample));
      return false;
    }
    return true;
  }

  String VectorCategory::toString(const AStringFormat* strfmt) const
  {
    DECLARE_UNUSED(strfmt);
    std::ostringstream os;
    os << "[";

    for (size_t i = 0; i < size(); ++i)
    {
      const auto cat = getCategory(i);

      if (!cat)
        os << "NA";
      else
        os << cat->getName();

      if (i != size() - 1) os << " ";
    }

    os << "]";
    return os.str();
  }

#ifndef SWIG
  std::ostream& operator<<(std::ostream& os, const VectorCategory& vec)
  {
    os << vec.toString() << std::endl;
    return os;
  }
#endif

} // namespace gstlrn
