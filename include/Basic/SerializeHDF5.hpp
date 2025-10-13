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

#ifdef HDF5

#include "geoslib_define.h"

#include "Basic/ASerializable.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/VectorT.hpp"

#include <H5Cpp.h>

#include <optional>

namespace gstlrn
{
 
namespace SerializeHDF5
{
  /**
   * @brief Map values to corresponding HDF5 C++ types
   */
  inline H5::DataType getHDF5Type([[maybe_unused]] const I32 a)
  {
    return H5::PredType::NATIVE_INT;
  }
  inline H5::DataType getHDF5Type([[maybe_unused]] const double a)
  {
    return H5::PredType::NATIVE_DOUBLE;
  }
  inline H5::DataType getHDF5Type([[maybe_unused]] const long a)
  {
    return H5::PredType::NATIVE_LONG;
  }
  inline H5::DataType getHDF5Type([[maybe_unused]] const bool a)
  {
    return H5::PredType::NATIVE_HBOOL;
  }
  inline H5::DataType getHDF5Type([[maybe_unused]] const std::string &a)
  {
    return H5::StrType {0, H5T_VARIABLE};
  }
  inline H5::DataType getHDF5Type([[maybe_unused]] const char a[])
  {
    return H5::StrType {0, H5T_VARIABLE};
  }

  /**
   * @brief Read individual value (primitive type) from group
   *
   * @param[in] grp HDF5 Group from which to read value
   * @param[in] name Value name
   * @param[out] value Return value
   * @return True if success
   */
  template<typename T>
  bool readValue(const H5::H5Object& grp, const String& name, T& value);

  /**
   * @brief Write individual value (primitive type) to group
   *
   * Here we make use of H5::Attributes to store individual values
   * (e.g.  class members of primitive types). Prefer _createAttribute
   * for strings.
   *
   * @param[in,out] grp HDF5 Group in which to write value
   * @param[in] name Value name
   * @param[out] value Value to write
   * @return True if success
   */
  template<typename T>
  bool writeValue(H5::H5Object& grp, const String& name, const T& value);

  /**
   * @brief Open HDF5 file in read mode, check metadata
   */
  inline H5::H5File fileOpenRead(const String& fname)
  {
    // Build the multi-platform filename
    const auto filepath = ASerializable::buildFileName(1, fname, true);

    H5::H5File file {filepath, H5F_ACC_RDONLY};

    if (!file.nameExists("gstlearn metadata"))
    {
      messerr("File %s doesn't contain Gstlearn metadata…", filepath.c_str());
      return file;
    }

    auto metadata = file.openGroup("gstlearn metadata");
    String version;
    readValue(metadata, "Format version", version);
    if (version != "1.0.0")
    {
      messerr("File %s has format version %s, expected 1.0.0", filepath.c_str(),
              version.c_str());
    }
    return file;
  }

  inline String getFileClass(const String& filepath, bool verbose = false)
  {
    H5::H5File file {filepath, H5F_ACC_RDONLY};

    if (!file.nameExists("gstlearn metadata"))
    {
      messerr("File %s doesn't contain Gstlearn metadata…", filepath.c_str());
      return String();
    }

    auto metadata = file.openGroup("gstlearn metadata");
    if (verbose)
      message("Metadata:\n");

    String version;
    readValue(metadata, "Format version", version);
    if (version != "1.0.0")
    {
      messerr("File %s has format version %s, expected 1.0.0", filepath.c_str(),
              version.c_str());
      return String();
    }
    if (verbose)
      message("- Version = %s\n", version.c_str());

    String classType;
    readValue(metadata, "Class_Type", classType);
    if (verbose)
      message("- Class_Type = %s\n", classType.c_str());

    return classType;
  }

  /**
   * @brief Open HDF5 file in write mode, write metadata
   */
  inline H5::H5File fileOpenWrite(const ASerializable& parent,
                                  const String& fname)
  {
    // Build the multi-platform filename
    const auto filepath = ASerializable::buildFileName(2, fname, true);

    H5::H5File file {filepath, H5F_ACC_TRUNC};
    auto metadata = file.createGroup("gstlearn metadata");
    writeValue(metadata, "Description",
               "This file is used to Serialize gstlearn's internal data structures");
    writeValue(metadata, "Format version", "1.0.0");
    writeValue(metadata, "Class_Type", parent._getNFName());
    return file;
  }

  /**
   * @brief Read a HDF5 DataSet into a generic VectorT
   *
   * @param[in] grp HDF5 group containing the variable to read
   * @param[in] title Name of HDF5 variable to read
   * @param[out] vec Vector to be filled with HDF5 data (will be resized)
   * @return true if success
   */
  template<typename T>
  bool readVec(const H5::Group& grp, const String& title, VectorT<T>& vec);

  /**
   * @brief Read a HDF5 string DataSet into a VectorString
   *
   * @param[in] grp HDF5 group containing the string variable to read
   * @param[in] title Name of HDF5 string variable to read
   * @param[out] vec VectorString to be filled with HDF5 data (will be resized)
   * @return true if success
   */
  template<>
  inline bool readVec(const H5::Group& grp, const String& title, VectorString& vec);

  /**
   * @brief Write a generic VectorT into a HDF5 DataSet
   *
   * @param[in,out] grp HDF5 Group containing the variable to write
   * @param[in] title Name of HDF5 DataSet to write
   * @param[in] vec Vector to write into HDF5
   * @return true if success
   */
  template<typename T>
  bool writeVec(H5::Group& grp, const String& title, const VectorT<T>& vec);

  /**
   * @brief Write a VectorString into a HDF5 string DataSet
   *
   * @param[in,out] grp HDF5 Group containing the string variable to write
   * @param[in] title Name of HDF5 string DataSet to write
   * @param[in] vec VectorString to write to HDF5
   * @return true if success
   */
  template<>
  inline bool writeVec(H5::Group& grp, const String& title, const VectorString& vec);

  /**
   * @brief Extract a group inside a parent group
   *
   * @param[in] parent Parent group
   * @param[in] name Name of the group to find in parent
   * @return Group if found else nullopt
   */
  inline std::optional<H5::Group> getGroup(const H5::Group& parent, const String& name, bool verbose = true)
  {
    if (!parent.nameExists(name) || parent.childObjType(name) != H5O_TYPE_GROUP)
    {
      std::string parent_name;
      parent.getObjName(parent_name);
      if (verbose)
        messerr("Cannot find group %s in parent group %s", name.c_str(),
                parent_name.c_str());
      return std::nullopt;
    }

    auto grp = parent.openGroup(name);
    return grp;
  }
}; // namespace SerializeHDF5

template<typename T>
bool SerializeHDF5::readVec(const H5::Group& grp, const String& title, VectorT<T>& vec)
{

  const auto grp_name = grp.getObjName();

  if (!grp.nameExists(title) || grp.childObjType(title) != H5O_TYPE_DATASET)
  {
    messerr("Cannot read HDF5 Variable of name %s in group %s", title.c_str(),
            grp_name.c_str());
    return false;
  }

  const auto data = grp.openDataSet(title);
  const auto ds   = data.getSpace();

  // Assume variable is of dim 1
  if (ds.getSimpleExtentNdims() != 1)
  {
    messerr("HDF5 Variable of name %s in group %s has %d dims, but we expect only 1",
            title.c_str(), grp_name.c_str(), ds.getSimpleExtentNdims());
    return false;
  }

  hsize_t dim {};
  ds.getSimpleExtentDims(&dim);
  vec.resize(dim);

  data.read(vec.data(), getHDF5Type(T {}));

  return true;
}

template<>
bool SerializeHDF5::readVec(const H5::Group& grp, const String& title, VectorString& vec)
{

  const auto grp_name = grp.getObjName();

  if (!grp.nameExists(title) || grp.childObjType(title) != H5O_TYPE_DATASET)
  {
    messerr("Cannot read HDF5 Variable of name %s in group %s", title.c_str(),
            grp_name.c_str());
    return false;
  }

  const auto data = grp.openDataSet(title);
  const auto ds   = data.getSpace();

  // Assume variable is of dim 1
  if (ds.getSimpleExtentNdims() != 1)
  {
    messerr(
      "HDF5 String Variable of name %s in group %s has %d dims, but we expect only 1",
      title.c_str(), grp_name.c_str(), ds.getSimpleExtentNdims());
    return false;
  }

  const auto dim = ds.getSimpleExtentDims(0);
  vec.resize(dim);

  // Use a vector of char* managed by HDF5 to read string data
  std::vector<char*> data_ptr(dim);
  data.read(static_cast<void *>(data_ptr.data()), H5::StrType {0, H5T_VARIABLE});

  // copy char pointers into gstlearn managed string vector
  for (size_t i = 0; i < data_ptr.size(); ++i)
  {
    vec[i] = data_ptr[i];
  }

  return true;
}

template<typename T>
bool SerializeHDF5::writeVec(H5::Group& grp, const String& title, const VectorT<T>& vec)
{
  // Allow processing empty string: nothing is done
  if (vec.empty()) return true;

  hsize_t dim = vec.size();
  H5::DataSpace ds {1, &dim};

  auto data = grp.createDataSet(title, getHDF5Type(vec[0]), ds);
  data.write(vec.constData(), getHDF5Type(vec[0]));
  return true;
}

template<>
bool SerializeHDF5::writeVec(H5::Group& grp,
                             const String& title,
                             const VectorString& vec)
{
  // Allow processing empty vector: nothing is done
  if (vec.empty()) return true;

  // generate a vector of char * to feed HDF5
  std::vector<const char*> data_ptr(vec.size());
  for (size_t i = 0; i < vec.size(); ++i)
  {
    data_ptr[i] = vec[i].c_str();
  }

  hsize_t dim = vec.size();
  H5::DataSpace ds {1, &dim};

  const auto var = grp.createDataSet(title, H5::StrType {0, H5T_VARIABLE}, ds);
  var.write(static_cast<void *>(data_ptr.data()), H5::StrType {0, H5T_VARIABLE});
  return true;
}

template<typename T>
bool SerializeHDF5::readValue(const H5::H5Object& grp, const String& name, T& value)
{
  const auto grp_name = grp.getObjName();

  if (!grp.attrExists(name))
  {
    messerr("Could not read value '%s' in group '%s': attribute does not exist", name.c_str(),
            grp_name.data());
    return false;
  }

  auto attr = grp.openAttribute(name);
  if (attr.getDataType() != getHDF5Type(value))
  {
    messerr("Could not read value '%s' in group '%s': mismatch in datatypes", name.c_str(),
            grp_name.data());
    return false;
  }

  const auto type = getHDF5Type(value);
  if constexpr (std::is_same<T, std::string>::value)
    attr.read(type, value);
  else if constexpr (std::is_convertible<T, std::string>::value)
    attr.read(type, std::string {value});
  else
    attr.read(type, &value);
  attr.close();
  return true;
}

template<typename T>
bool SerializeHDF5::writeValue(H5::H5Object& grp, const String& name, const T& value)
{
  std::string grp_name;
  grp.getObjName(grp_name);

  const hsize_t dim = 1;
  const H5::DataSpace ds {1, &dim};
  auto attr       = grp.createAttribute(name, getHDF5Type(value), ds);
  const auto type = getHDF5Type(value);
  if constexpr (std::is_same<T, std::string>::value)
    attr.write(type, value);
  else if constexpr (std::is_convertible<T, std::string>::value)
    attr.write(type, std::string {value});
  else
    attr.write(type, &value);

  return true;
}

} // namespace gstlrn

#endif // HDF5
