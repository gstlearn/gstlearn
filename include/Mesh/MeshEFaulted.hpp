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
#pragma once

#include "Basic/VectorNumT.hpp"
#include "Matrix/MatrixDense.hpp"
#include "Matrix/MatrixInt.hpp"
#include "Mesh/MeshEStandard.hpp"

namespace gstlrn
{
  /**
   * Faulted meshing defined in the Euclidean space
   */
  class GSTLEARN_EXPORT MeshEFaulted: public MeshEStandard
  {
  public:
    MeshEFaulted();
    MeshEFaulted(const MeshEFaulted& m);
    MeshEFaulted& operator=(const MeshEFaulted& m);
    virtual ~MeshEFaulted();

    /// Interface to AStringable
    String toString(const AStringFormat* strfmt = nullptr) const override;

    /// ASeriazable interface
    String getNFName() const override { return "MeshEFaulted"; }

#ifdef HDF5
    bool deserializeH5(H5::Group& grp) override;
    bool serializeH5(H5::Group& grp) const override;
#endif

    /// Interface for AMesh
    static MeshEFaulted* createFromExternal(
      const MatrixDense& apices,
      const MatrixInt& meshes,
      const VectorInt& fault_ids,
      bool verbose = false);

    Id reset(
      const MatrixDense& apices,
      const MatrixInt& meshes,
      const VectorInt& faultIds,
      bool verbose = false);
    Id resetFromVectors(
      Id ndim,
      Id napexpermesh,
      const VectorDouble& apices,
      const VectorInt& meshes,
      const VectorInt& faultIds,
      bool byCol = true,
      bool verbose = false);

    // Override MeshEStandard
    Id reset(
      const MatrixDense& apices,
      const MatrixInt& meshes,
      bool verbose = false) override;
    Id resetFromVectors(
      Id ndim,
      Id napexpermesh,
      const VectorDouble& apices,
      const VectorInt& meshes,
      bool byCol = true,
      bool verbose = false) override;
    Id resetFromTurbo(const MeshETurbo& turbo, bool verbose = false) override;

    Id getNVirtualFaultMeshes() const;

    inline const VectorInt& getFaultIds() const { return _faultIds; }

  private:
    void _deallocate();
    Id _recopy(const MeshEFaulted& m);
    void _createBalltree() const override;
    void _resetFaultIds();

  private:
    VectorInt _faultIds; // Dimension: N=Nmesh;
  };
} // namespace gstlrn
