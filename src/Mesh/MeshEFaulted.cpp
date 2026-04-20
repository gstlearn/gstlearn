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
#include "Mesh/MeshEFaulted.hpp"
#include "Basic/AException.hpp"
#include "Basic/SerializeHDF5.hpp"
#include "Matrix/MatrixDense.hpp"
#include "Mesh/MeshEStandard.hpp"
#include "Tree/BallFaulted.hpp"

namespace gstlrn
{

  MeshEFaulted::MeshEFaulted()
    : MeshEStandard()
    , _faultIds()
  {
  }

  MeshEFaulted::MeshEFaulted(const MeshEFaulted& m)
    : MeshEStandard(m)
  {
    _recopy(m);
  }

  MeshEFaulted& MeshEFaulted::operator=(const MeshEFaulted& m)
  {
    _recopy(m);
    return *this;
  }

  MeshEFaulted::~MeshEFaulted()
  {
    _deallocate();
  }

  /****************************************************************************/
  /*!x
  ** Print the contents of the meshing
  **
  ** \param[in] strfmt    Format for printout
  **
  *****************************************************************************/
  String MeshEFaulted::toString(const AStringFormat* strfmt) const
  {
    std::stringstream sstr;
    sstr << toStrTitle(0, "Faulted Meshing");
    sstr << "Number of Virtual meshes  = " << getNVirtualFaultMeshes()
         << std::endl;
    sstr << MeshEStandard::toString(strfmt);
    return sstr.str();
  }

  MeshEFaulted* MeshEFaulted::createFromExternal(
    const MatrixDense& apices,
    const MatrixInt& meshes,
    const VectorInt& fault_ids,
    bool verbose)
  {
    MeshEFaulted* mesh = new MeshEFaulted;
    mesh->reset(apices, meshes, fault_ids, verbose);
    return mesh;
  }

  /****************************************************************************/
  /*!
  ** Create the meshing (from mesh information)
  **
  ** \param[in]  apices          Matrix of Apices
  ** \param[in]  meshes          Array of mesh indices
  ** \param[in]  faultIds        Vector of fault ids
  ** \param[in]  verbose         Verbose flag
  **
  *****************************************************************************/
  Id MeshEFaulted::reset(
    const MatrixDense& apices,
    const MatrixInt& meshes,
    const VectorInt& faultIds,
    bool verbose)
  {
    // Core allocation
    _faultIds = faultIds;

    MeshEStandard::reset(apices, meshes, verbose);

    return (0);
  }

  /****************************************************************************/
  /*!
  ** Create the meshing (from mesh information)
  **
  ** \param[in]  ndim            Space Dimension
  ** \param[in]  napexpermesh    Number of apices per mesh
  ** \param[in]  apices          Vector of Apex information
  ** \param[in]  meshes          Vector of mesh indices
  ** \param[in]  faultIds        Vector of fault ids
  ** \param[in]  byCol           true for Column major; false for Row Major
  ** \param[in]  verbose         Verbose flag
  **
  ** \remark The argument 'byCol' concerns 'apices' and 'meshes'
  **
  *****************************************************************************/
  Id MeshEFaulted::resetFromVectors(
    Id ndim,
    Id napexpermesh,
    const VectorDouble& apices,
    const VectorInt& meshes,
    const VectorInt& faultIds,
    bool byCol,
    bool verbose)
  {
    // Core allocation
    _faultIds = faultIds;

    MeshEStandard::resetFromVectors(
      ndim, napexpermesh, apices, meshes, byCol, verbose);

    return (0);
  }

  Id MeshEFaulted::resetFromTurbo(const MeshETurbo& turbo, bool verbose)
  {
    _resetFaultIds();
    MeshEStandard::resetFromTurbo(turbo, verbose);

    return 0;
  }

  /****************************************************************************/
  /*!
  ** Create the meshing (from mesh information)
  **
  ** \param[in]  apices          Matrix of Apices
  ** \param[in]  meshes          Array of mesh indices
  ** \param[in]  verbose         Verbose flag
  **
  *****************************************************************************/
  Id MeshEFaulted::reset(
    const MatrixDense& apices,
    const MatrixInt& meshes,
    bool verbose)
  {
    _resetFaultIds();
    MeshEStandard::reset(apices, meshes, verbose);

    return 0;
  }

  /****************************************************************************/
  /*!
  ** Create the meshing (from mesh information)
  **
  ** \param[in]  ndim            Space Dimension
  ** \param[in]  napexpermesh    Number of apices per mesh
  ** \param[in]  apices          Vector of Apex information
  ** \param[in]  meshes          Vector of mesh indices
  ** \param[in]  byCol           true for Column major; false for Row Major
  ** \param[in]  verbose         Verbose flag
  **
  ** \remark The argument 'byCol' concerns 'apices' and 'meshes'
  **
  *****************************************************************************/
  Id MeshEFaulted::resetFromVectors(
    Id ndim,
    Id napexpermesh,
    const VectorDouble& apices,
    const VectorInt& meshes,
    bool byCol,
    bool verbose)
  {
    _resetFaultIds();
    MeshEStandard::resetFromVectors(
      ndim, napexpermesh, apices, meshes, byCol, verbose);
    return 0;
  }

  Id MeshEFaulted::getNVirtualFaultMeshes() const
  {
    if (_faultIds.empty()) return 0;

    Id nFaultMeshes = 0;
    Id nmesh = getNMeshes();

    for (Id imesh = 0; imesh < nmesh; imesh++)
    {
      nFaultMeshes += _faultIds[imesh] != 0;
    }

    return nFaultMeshes;
  }

  void MeshEFaulted::_createBalltree() const
  {
    _ballTree = std::make_unique<BallFaulted>(this, 10, true);
  }

  void MeshEFaulted::_deallocate() {}

  Id MeshEFaulted::_recopy(const MeshEFaulted& m)
  {
    _deallocate();

    _faultIds = m._faultIds;
    MeshEStandard::_recopy(m);
    return (0);
  }

  void MeshEFaulted::_resetFaultIds()
  {
    _faultIds.resize(getNMeshes());
  }

#ifdef HDF5
  bool MeshEFaulted::deserializeH5(H5::Group& grp)
  {
    auto meshG = SerializeHDF5::getGroup(grp, "MeshEFaulted");
    if (!meshG)
    {
      return false;
    }

    /* Read the grid characteristics */
    bool ret = true;

    ret = ret && MeshEFaulted::deserializeH5(*meshG);
    VectorInt faultIds;

    ret = ret && SerializeHDF5::readVec(*meshG, "FaultIds", faultIds);

    if (ret)
    {
      _faultIds = faultIds;
    }

    return ret;
  }

  bool MeshEFaulted::serializeH5(H5::Group& grp) const
  {
    auto meshG = grp.createGroup("MeshEFaulted");

    bool ret = true;
    ret = ret && MeshEStandard::serializeH5(meshG);
    ret = ret && SerializeHDF5::writeVec(meshG, "FaultIds", _faultIds);

    return ret;
  }
#endif
} // namespace gstlrn
