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
#include "DataBase/DbCol.hpp"
#include "Basic/SerializeHDF5.hpp"
#include "Basic/VectorHelper.hpp"

namespace gstlrn
{
#ifdef HDF5

  bool DbCol::serializeH5(H5::Group& grp) const
  {
    SerializeHDF5::writeValue(grp, "Name", getName());
    SerializeHDF5::writeValue(grp, "NVersion", getNVersions());
    SerializeHDF5::writeValue(grp, "ForbidNA", forbidNA());

    const bool isCategory =
      std::holds_alternative<Array2D<VectorCategory>>(this->_data);

    SerializeHDF5::writeValue(grp, "IsCategory", isCategory);

    bool success = true;

    visitData(
      [&](const auto& array)
      {
        if (!success) return;

        using ArrayType = std::decay_t<decltype(array)>;
        using ValueType = typename ArrayType::value_type;

        hsize_t dims[2] = {static_cast<hsize_t>(array.outer()),
                           static_cast<hsize_t>(array.inner())};
        H5::DataSpace space(2, dims);

        /*
         * Categorical case
         */
        if constexpr (std::is_same_v<ArrayType, Array2D<VectorCategory>>)
        {
          /*
           * VectorCategory case.
           *
           * The dictionary is common to all versions and is therefore
           * serialized only once.
           */
          const auto& vc = array.getBuffer();

          H5::Group dictGrp = grp.createGroup("Dictionary");

          if (!vc.getDictionary().serializeH5(dictGrp))
          {
            success = false;
            return;
          }

          std::vector<Id> values;
          values.reserve(array.outer() * array.inner());

          for (Id iversion = 0; iversion < array.outer(); ++iversion)
          {
            /*
             * The VectorCategory contains all versions consecutively.
             */
            const auto offset = static_cast<size_t>(iversion * array.inner());

            for (Id isample = 0; isample < array.inner(); ++isample)
            {
              const auto cat =
                vc.getCategory(offset + static_cast<size_t>(isample));

              if (!cat)
              {
                values.push_back(-1);
                continue;
              }
              values.push_back(cat->getIdentifier());
            }
          }

          const auto& type = SerializeHDF5::getHDF5Type<Id>();

          auto ds = grp.createDataSet("Values", type, space);

          if (!values.empty()) ds.write(values.data(), type);
        }

        /*
         * String case
         */
        else if constexpr (std::is_same_v<ValueType, String>
                           || std::is_same_v<ValueType, std::string>)
        {
          H5::StrType strType(H5::PredType::C_S1, H5T_VARIABLE);

          auto ds = grp.createDataSet("Values", strType, space);

          std::vector<const char*> buffer;

          buffer.reserve(array.getBuffer().size());

          for (const auto& s: array.getBuffer()) buffer.push_back(s.c_str());

          if (!buffer.empty()) ds.write(buffer.data(), strType);
        }

        else
        {
          /*
           * Numeric case
           */
          const auto& type = SerializeHDF5::getHDF5Type<ValueType>();

          auto ds = grp.createDataSet("Values", type, space);

          ds.write(array.data(), type);
        }
      });

    return success;
  }

  bool DbCol::deserializeH5(H5::Group& grp)
  {
    String name;
    Id nversion = 1;
    bool forbidNA = false;

    SerializeHDF5::readValue(grp, "Name", name);
    SerializeHDF5::readValue(grp, "NVersion", nversion);
    SerializeHDF5::readValue(grp, "ForbidNA", forbidNA);

    setForbidNA(forbidNA);
    setName(name);

    auto ds = grp.openDataSet("Values");

    H5::DataSpace space = ds.getSpace();
    H5::DataType dtype = ds.getDataType();

    hsize_t dims[2];
    space.getSimpleExtentDims(dims);

    const Id outer = static_cast<Id>(dims[0]);
    const Id inner = static_cast<Id>(dims[1]);
    const auto size = static_cast<size_t>(outer * inner);

    /*
     * Check Categorical case
     */
    bool isCategory = false;
    SerializeHDF5::readValue(grp, "IsCategory", isCategory);

    if (isCategory)
    {
      /*
       * Read dictionary.
       */
      H5::Group dictGrp = grp.openGroup("Dictionary");

      Dictionary dict;

      if (!dict.deserializeH5(dictGrp)) return false;

      /*
       * Read category indices.
       */
      hid_t h5type = dtype.getId();

      if (H5Tget_class(h5type) != H5T_INTEGER) return false;

      if (H5Tget_size(h5type) != sizeof(Id) || H5Tget_sign(h5type) != H5T_SGN_2)
        return false;

      std::vector<Id> values(size);

      const auto& type = SerializeHDF5::getHDF5Type<Id>();

      if (!values.empty()) ds.read(values.data(), type);

      /*
       * Rebuild the VectorCategory.
       */
      VectorCategory vc(size, dict);

      if (!vc.setCategoriesByKey(values)) return false;

      /*
       * Build the Array2D.
       *
       * The dictionary is already stored in vc and is therefore
       * automatically common to all versions.
       */
      _data = Array2D<VectorCategory>(std::move(vc), outer);

      return true;
    }

    /*
     * Non-categorical cases
     */
    hid_t type = dtype.getId();

    switch (H5Tget_class(type))
    {
      case H5T_FLOAT:
      {
        if (H5Tget_size(type) == sizeof(double))
        {
          VectorDouble vec(size);
          ds.read(vec.data(), SerializeHDF5::getHDF5Type<double>());
          _data = Array2D<VectorDouble>(std::move(vec), outer);
          return true;
        }

        if (H5Tget_size(type) == sizeof(float))
        {
          VectorFloat vec(size);
          ds.read(vec.data(), SerializeHDF5::getHDF5Type<float>());
          _data = Array2D<VectorFloat>(std::move(vec), outer);
          return true;
        }
        break;
      }

      case H5T_INTEGER:
      {
        if (H5Tget_sign(type) == H5T_SGN_2 && H5Tget_size(type) == sizeof(Id))
        {
          VectorInt vec(size);
          ds.read(vec.data(), SerializeHDF5::getHDF5Type<Id>());
          _data = Array2D<VectorInt>(std::move(vec), outer);
          return true;
        }

        if (H5Tget_sign(type) == H5T_SGN_NONE
            && H5Tget_size(type) == sizeof(unsigned char))
        {
          VectorBool vec(size);
          ds.read(vec.data(), SerializeHDF5::getHDF5Type<unsigned char>());
          _data = Array2D<VectorBool>(std::move(vec), outer);
          return true;
        }
        break;
      }

      case H5T_STRING:
      {
        H5::StrType strType = ds.getStrType();

        std::vector<char*> buffer(size);

        ds.read(buffer.data(), strType);

        VectorString vec(size);

        for (size_t i = 0; i < size; i++) vec[i] = buffer[i];

        H5Dvlen_reclaim(
          strType.getId(), space.getId(), H5P_DEFAULT, buffer.data());

        _data = Array2D<VectorString>(std::move(vec), outer);
        return true;
      }

      default: break;
    }

    return false;
  }
#endif

  void DbCol::deleteSample(const Id isample)
  {
    std::visit(
      [isample](auto&& arg) { arg.deleteElement(isample); }, this->_data);
  }

  Id DbCol::getNSamples() const
  {
    return std::visit(
      [](auto&& arg) -> Id { return arg.inner(); }, this->_data);
  }

  Id DbCol::getNVersions() const
  {
    return std::visit(
      [](auto&& arg) -> Id { return arg.outer(); }, this->_data);
  }

  String DbCol::getTypeName() const
  {
    return std::visit(
      [](auto&& arg) -> String
      {
        using VectorType = std::decay_t<decltype(arg)>::vector_type;
        using ValueType = typename VectorType::value_type;
        return getGenericTypeName<ValueType>();
      },
      this->_data);
  }

  String DbCol::getDescr() const
  {
    String descr = this->getName();

    if (this->getNVersions() > 1)
      descr += " [" + std::to_string(this->getNVersions()) + "]";

    descr += " (Type: " + this->getTypeName() + ")";

    return descr;
  }

  bool DbCol::_checkVersion(Id iversion, Id nversion) const
  {
    if (iversion < 0 || iversion >= nversion)
    {
      messerr(
        "Version %d is invalid for Column '%s' which contains %d "
        "version(s).",
        static_cast<Id>(iversion), getName().c_str(),
        static_cast<Id>(nversion));
      return false;
    }
    return true;
  }

  bool DbCol::_checkSample(Id isample, Id nsample) const
  {
    if (isample < 0 || isample >= nsample)
    {
      messerr(
        "Sample index (%d) is out of bounds for Column '%s' which contains %d "
        "samples.",
        static_cast<Id>(isample), getName().c_str(), static_cast<Id>(nsample));
      return false;
    }
    return true;
  }

} // namespace gstlrn
