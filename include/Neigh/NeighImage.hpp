/******************************************************************************/
/* COPYRIGHT ARMINES, ALL RIGHTS RESERVED                                     */
/*                                                                            */
/* THE CONTENT OF THIS WORK CONTAINS CONFIDENTIAL AND PROPRIETARY             */
/* INFORMATION OF ARMINES. ANY DUPLICATION, MODIFICATION,                     */
/* DISTRIBUTION, OR DISCLOSURE IN ANY FORM, IN WHOLE, OR IN PART, IS STRICTLY */
/* PROHIBITED WITHOUT THE PRIOR EXPRESS WRITTEN PERMISSION OF ARMINES         */
/*                                                                            */
/* TAG_SOURCE_CG                                                              */
/******************************************************************************/
#pragma once

#include "gstlearn_export.hpp"
#include "geoslib_define.h"

// Enums
#include "Neigh/ENeigh.hpp"
#include "Neigh/ANeighParam.hpp"

#include "Basic/AStringable.hpp"
#include "Basic/ASerializable.hpp"

class Db;

class GSTLEARN_EXPORT NeighImage: public ANeighParam
{
public:
  NeighImage(int ndim = 2, const VectorInt radius = VectorInt(), int skip = 0);
  NeighImage(const NeighImage& r);
  NeighImage& operator=(const NeighImage& r);
  virtual ~NeighImage();

  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  virtual int getMaxSampleNumber(const Db* db) const override;
  virtual ENeigh getType() const override { return ENeigh::IMAGE; }

  int reset(int ndim, const VectorInt& image, int skip = 0);
  static NeighImage* create(int ndim, const VectorInt& image, int skip = 0);
  static NeighImage* createFromNF(const String& neutralFilename, bool verbose = false);
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
  int _skip;                  /* Skipping factor */
  VectorInt _imageRadius;     /* Vector of image neighborhood radius */
};
