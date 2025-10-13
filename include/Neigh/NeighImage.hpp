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

#include "Basic/ICloneable.hpp"
#include "gstlearn_export.hpp"
#include "geoslib_define.h"

#include "Space/ASpace.hpp"
#include "Enum/ENeigh.hpp"
#include "Neigh/ANeigh.hpp"

namespace gstlrn
{
class Db;
class DbGrid;

/**
 * \brief
 * Image Neighborhood definition.
 *
 * The Neighborhood is usually meant to select a sub-population from the input Data Base,
 * containing the active samples close to the target.
 *
 * This Neighborhood is only defined in the case when the Data and the Target belong
 * to the same grid.
 * This neighborhood is defined as a rectangular set of pixels, located around the target.
 * This rectangle is given by its half-extension in each space dimension (called 'radius')
 * As the number of pixels grows fast with the space dimension, it is offered to sample
 * them by specifying a skipping factor, so as to retain only 1 / (1 + skip) of them.
 */
class GSTLEARN_EXPORT NeighImage: public ANeigh
{
public:
  NeighImage(const VectorInt &radius = VectorInt(),
             Id skip = 0,
             const ASpaceSharedPtr& space = ASpaceSharedPtr());
  NeighImage(const NeighImage& r);
  NeighImage& operator=(const NeighImage& r);
  virtual ~NeighImage();

  IMPLEMENT_CLONING(NeighImage)
  /// Interface for ANeigh
  void getNeigh(Id iech_out, VectorInt& ranks) override;
  Id getNSampleMax(const Db* db) const override;
  bool hasChanged(Id iech_out) const override;
  ENeigh getType() const override { return ENeigh::fromKey("IMAGE"); }

  /// Interface for AStringable
  String toString(const AStringFormat* strfmt = nullptr) const override;

  static NeighImage* create(const VectorInt& radius, Id skip = 0, const ASpaceSharedPtr& space = ASpaceSharedPtr());
  static NeighImage* createFromNF(const String& NFFilename, bool verbose = true);

  Id getSkip() const { return _skip; }
  const VectorInt& getImageRadius() const { return _imageRadius; }
  Id getImageRadius(Id idim) const { return _imageRadius[idim]; }

  void setImageRadius(const VectorInt& imageRadius) { _imageRadius = imageRadius; }
  void setSkip(Id skip) { _skip = skip; }

  DbGrid* buildImageGrid(const DbGrid* dbgrid, Id seed) const;

protected:
  bool _deserializeAscii(std::istream& is, bool verbose = false) override;
  bool _serializeAscii(std::ostream& os, bool verbose = false) const override;
#ifdef HDF5
  bool _deserializeH5(H5::Group& grp, bool verbose = false) override;
  bool _serializeH5(H5::Group& grp, bool verbose = false) const override;
#endif
  String _getNFName() const override { return "NeighImage"; }

private:
  void _uimage(Id iech_out, VectorInt& ranks);

private:
  Id _skip;                  /* Skipping factor */
  VectorInt _imageRadius;     /* Vector of image neighborhood radius */
};
}