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
#include "Mesh/VectorMeshes.hpp"
#include "Mesh/MeshETurbo.hpp"

namespace gstlrn
{
  VectorMeshes::VectorMeshes(const std::vector<const AMesh*>& meshes)
    : AStringable()
    , _meshes(meshes)
  {
  }

  VectorMeshes::VectorMeshes(Id nmesh)
    : AStringable()
    , _meshes(std::vector<const AMesh*>(nmesh, nullptr))
  {
  }

  VectorMeshes::VectorMeshes(const VectorMeshes& r)
    : AStringable(r)
    , _meshes(r._meshes)
  {
  }

  VectorMeshes& VectorMeshes::operator=(const VectorMeshes& r)
  {
    if (this != &r)
    {
      AStringable::operator=(r);
      _meshes = r._meshes;
    }
    return *this;
  }

  VectorMeshes::~VectorMeshes() {}

  /****************************************************************************/
  /*!
  ** Print the contents of the meshing
  **
  ** \param[in] strfmt    Format for printout
  **
  *****************************************************************************/
  String VectorMeshes::toString(const AStringFormat* strfmt) const
  {
    DECLARE_UNUSED(strfmt);

    std::stringstream sstr;
    Id nmesh = static_cast<Id>(_meshes.size());
    sstr << toStrTitle(1, "Contents of the VectorMeshes") << std::endl;
    sstr << "It contains " << nmesh << " mesh(es)" << std::endl;

    for (Id imesh = 0; imesh < nmesh; imesh++)
      sstr << _meshes[imesh]->toString();

    return sstr.str();
  }

  /**
   * @brief Check if a series of Meshes (included in 'meshes') are Turbo
   *
   * @return True if ALL meshes are TURBO
   */
  bool VectorMeshes::isTurbo() const
  {
    if (_meshes.empty()) return false;
    for (Id imesh = 0, nmesh = static_cast<Id>(_meshes.size()); imesh < nmesh;
         imesh++)
    {
      const auto* mTurbo = dynamic_cast<const MeshETurbo*>(_meshes[imesh]);
      if (mTurbo == nullptr) return false;
    }
    return true;
  }

  bool VectorMeshes::allDefined() const
  {
    Id nmesh = static_cast<Id>(_meshes.size());
    for (Id imesh = 0; imesh < nmesh; imesh++)
      if (_meshes[imesh] == nullptr) return false;
    return true;
  }

  void VectorMeshes::replace(Id ind, const AMesh* mesh)
  {
    _meshes[ind] = mesh;
  }

} // namespace gstlrn
