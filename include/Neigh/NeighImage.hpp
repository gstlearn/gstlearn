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
  NeighImage(int ndim = 2, int skip = 0, const VectorInt image = VectorInt());
  NeighImage(const NeighImage& r);
  NeighImage& operator=(const NeighImage& r);
  virtual ~NeighImage();

  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  virtual int getMaxSampleNumber(const Db* db) const;
  virtual ENeigh getType() const override { return ENeigh::IMAGE; }

  int reset(int ndim = 2, int skip = 0, const VectorInt& image = VectorInt());
  static NeighImage* create(int ndim, int skip, const VectorInt& image);
  static NeighImage* createFromNF2(const String& neutralFilename, bool verbose = false);

  int dumpToNF2(const String& neutralFilename, bool verbose = false) const;

  int getSkip() const { return _skip; }
  const VectorInt& getImageRadius() const { return _imageRadius; }
  int getImageRadius(int idim) const { return _imageRadius[idim]; }

  void setImageRadius(const VectorInt& imageRadius) { _imageRadius = imageRadius; }
  void setSkip(int skip) { _skip = skip; }

protected:
  virtual int _deserialize2(std::istream& is, bool verbose = false) override;
  virtual int _serialize2(std::ostream& os, bool verbose = false) const override;

private:
  int _skip;                  /* Skipping factor */
  VectorInt _imageRadius;     /* Vector of image neighborhood radius */
};
