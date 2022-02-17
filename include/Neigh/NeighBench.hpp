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

// TODO : inherits from ASpaceObject (see _init)
class GSTLEARN_EXPORT NeighBench: public ANeighParam
{
public:
  NeighBench(int ndim = 2, bool flag_xvalid = false, double width = 0.);
  NeighBench(const NeighBench& r);
  NeighBench& operator=(const NeighBench& r);
  virtual ~NeighBench();

  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  int reset(int ndim = 2, bool flag_xvalid = false, double width = 0);

  int dumpToNF(const String& neutralFilename, bool verbose = false) const;
  static NeighBench* create(int ndim = 2, bool flag_xvalid = false, double width = 0);
  static NeighBench* createFromNF(const String& neutralFilename, bool verbose = false);
  static NeighBench* createFromNF2(const String& neutralFilename, bool verbose = false);

  virtual int getMaxSampleNumber(const Db* db) const override;
  virtual ENeigh getType() const override { return ENeigh::BENCH; }

  double getWidth() const { return _width; }
  void setWidth(double width) { _width = width; }

protected:
  virtual int _deserialize(FILE* file, bool verbose = false);
  virtual int _serialize(FILE* file, bool verbose = false) const override;

  virtual int _deserialize2(std::istream& is, bool verbose = false) override;

private:
  double _width;                 /* Width of the slice - bench */
};
