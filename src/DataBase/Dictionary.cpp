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
#include "DataBase/Dictionary.hpp"
#include "Basic/SerializeHDF5.hpp"
#include "DataBase/Category.hpp"

namespace gstlrn
{
  String Dictionary::toString(const AStringFormat* strfmt) const
  {
    DECLARE_UNUSED(strfmt);
    std::ostringstream os;

    // Entête facultatif ou titre
    os << "Dictionary (" << getNCategories() << " categories):\n";

    // Parcours de la map interne _data
    for (const auto& [key, label]: _data)
    {
      os << "  " << key << " : " << label << "\n";
    }

    return os.str();
  }

  /**
   * @brief Adds a category to the dictionary.
   *
   * @param key Category identifier.
   * @param val Category label.
   *
   * @return @c true if the category was successfully added,
   *         @c false otherwise.
   */
  bool Dictionary::addCategory(const Id key, const String& val)
  {
    for (const auto& p: this->_data)
    {
      // Check that the same entry does not already exist in the dictionary
      if (p.first == key || p.second == val)
      {
        return false;
      }
    }

    // Add the new entry to the dictionary
    this->_data[key] = val;
    return true;
  }

  /**
   * @brief Checks whether a category is present in the dictionary.
   *
   * Both the category identifier and its label are checked.
   *
   * @param cat Category identifier and label to check.
   *
   * @return @c true if the specified category is present,
   *         @c false otherwise.
   */
  bool Dictionary::hasCategory(const Category& cat) const
  {
    const auto it = this->_data.find(cat.getId());
    if (it == this->_data.end()) return false;
    if (it->second != cat.getName()) return false;
    return true;
  }

  /**
   * @brief Returns all categories stored in the dictionary.
   *
   * Categories are returned in increasing order of their identifiers.
   *
   * @return Vector containing all categories.
   */
  std::vector<Category> Dictionary::getCategories() const
  {
    std::vector<Category> categories;
    categories.reserve(this->_data.size());

    for (const auto& [key, val]: this->_data) categories.emplace_back(key, val);

    return categories;
  }

  /**
   * @brief Retrieves a category by its identifier.
   *
   * @param key Category identifier.
   *
   * @return The corresponding category, or an empty optional if the
   *         identifier is not present in the dictionary.
   */
  std::optional<Category> Dictionary::getCategory(const Id key) const
  {
    const auto it = this->_data.find(key);

    if (it == this->_data.end())
    {
      messerr(
        "Category with Id %d is not present in the Dictionary.",
        static_cast<Id>(key));
      return {};
    }
    return Category{it->first, it->second};
  }

#ifdef HDF5

  bool Dictionary::serializeH5(H5::Group& grp) const
  {
    const auto categories = getCategories();
    const size_t ncategory = categories.size();

    /*
     * Ids
     */
    hsize_t dims[1] = {static_cast<hsize_t>(ncategory)};
    H5::DataSpace space(1, dims);

    const auto& idType = SerializeHDF5::getHDF5Type<Id>();
    auto dsIds = grp.createDataSet("Ids", idType, space);

    std::vector<Id> ids;
    ids.reserve(ncategory);

    for (const auto& cat: categories) ids.push_back(cat.getId());

    if (!ids.empty()) dsIds.write(ids.data(), idType);

    /*
     * Labels
     */
    H5::StrType strType(H5::PredType::C_S1, H5T_VARIABLE);

    auto dsLabels = grp.createDataSet("Labels", strType, space);

    std::vector<const char*> labels;
    labels.reserve(ncategory);

    for (const auto& cat: categories) labels.push_back(cat.getName().c_str());

    if (!labels.empty()) dsLabels.write(labels.data(), strType);

    return true;
  }

  bool Dictionary::deserializeH5(H5::Group& grp)
  {
    auto dsIds = grp.openDataSet("Ids");
    auto dsLabels = grp.openDataSet("Labels");

    H5::DataSpace idSpace = dsIds.getSpace();

    hsize_t dims[1];
    idSpace.getSimpleExtentDims(dims);

    const size_t ncategory = dims[0];

    /*
     * Ids
     */
    VectorInt ids(ncategory);

    dsIds.read(ids.data(), SerializeHDF5::getHDF5Type<Id>());

    /*
     * Labels
     */
    H5::StrType strType = dsLabels.getStrType();

    std::vector<char*> buffer(ncategory);

    if (ncategory > 0) dsLabels.read(buffer.data(), strType);

    /*
     * Rebuild dictionary
     */
    this->_data.clear();

    for (size_t i = 0; i < ncategory; i++) this->_data[ids[i]] = buffer[i];

    /*
     * Release variable-length strings allocated by HDF5
     */
    if (ncategory > 0)
    {
      H5Dvlen_reclaim(
        strType.getId(), idSpace.getId(), H5P_DEFAULT, buffer.data());
    }

    return true;
  }

#endif

} // namespace gstlrn
