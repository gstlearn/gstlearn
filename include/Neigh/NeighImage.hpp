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

#include "gstlearn_export.hpp"
#include "geoslib_define.h"

#include "Space/ASpace.hpp"
#include "Enum/ENeigh.hpp"
#include "Neigh/ANeigh.hpp"

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
             int skip = 0,
             const ASpaceSharedPtr& space = ASpaceSharedPtr());
  NeighImage(const NeighImage& r);
  NeighImage& operator=(const NeighImage& r);
  virtual ~NeighImage();

  /// Interface for ANeigh
  virtual void getNeigh(int iech_out, VectorInt& ranks) override;
  virtual int getNSampleMax(const Db* db) const override;
  virtual bool hasChanged(int iech_out) const override;
  virtual ENeigh getType() const override { return ENeigh::fromKey("IMAGE"); }

  /// Interface for AStringable
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  static NeighImage* create(const VectorInt& radius, int skip = 0, const ASpaceSharedPtr& space = ASpaceSharedPtr());
  static NeighImage* createFromNF(const String& neutralFilename, bool verbose = true);

  int getSkip() const { return _skip; }
  const VectorInt& getImageRadius() const { return _imageRadius; }
  int getImageRadius(int idim) const { return _imageRadius[idim]; }

  void setImageRadius(const VectorInt& imageRadius) { _imageRadius = imageRadius; }
  void setSkip(int skip) { _skip = skip; }

  DbGrid* buildImageGrid(const DbGrid* dbgrid, int seed) const;

protected:
  /// Interface for ASerializable
  virtual bool _deserialize(std::istream& is, bool verbose = false) override;
  virtual bool _serialize(std::ostream& os, bool verbose = false) const override;
  String _getNFName() const override { return "NeighImage"; }

private:
  void _uimage(int iech_out, VectorInt& ranks);

private:
  int _skip;                  /* Skipping factor */
  VectorInt _imageRadius;     /* Vector of image neighborhood radius */
};
