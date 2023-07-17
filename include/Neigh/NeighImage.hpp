/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://github.com/gstlearn                                       */
/* License: BSD 3 clauses                                                     */
/*                                                                            */
/******************************************************************************/
#pragma once

#include "gstlearn_export.hpp"
#include "geoslib_define.h"

#include "Enum/ENeigh.hpp"

#include "Neigh/ANeigh.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/ASerializable.hpp"

class Db;

class GSTLEARN_EXPORT NeighImage: public ANeigh
{
public:
  NeighImage(const VectorInt &radius = VectorInt(),
             int skip = 0,
             const ASpace *space = nullptr);
  NeighImage(const NeighImage& r);
  NeighImage& operator=(const NeighImage& r);
  virtual ~NeighImage();

  /// Interface for ANeigh
  virtual VectorInt getNeigh(int iech_out) override;
  virtual int getMaxSampleNumber(const Db* db) const override;
  virtual bool hasChanged(int iech_out) const override;
  virtual ENeigh getType() const override { return ENeigh::fromKey("IMAGE"); }

  /// Interface for AStringable
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  static NeighImage* create(const VectorInt& image, int skip = 0, const ASpace* space = nullptr);
  static NeighImage* createFromNF(const String& neutralFilename, bool verbose = true);

  int getSkip() const { return _skip; }
  const VectorInt& getImageRadius() const { return _imageRadius; }
  int getImageRadius(int idim) const { return _imageRadius[idim]; }

  void setImageRadius(const VectorInt& imageRadius) { _imageRadius = imageRadius; }
  void setSkip(int skip) { _skip = skip; }

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
