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

#include "Enum/ENeigh.hpp"

#include "Neigh/ANeighParam.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/ASerializable.hpp"

class Db;

class GSTLEARN_EXPORT NeighUnique: public ANeighParam
{
public:
  NeighUnique(int ndim = 2, bool flag_xvalid = false);
  NeighUnique(const NeighUnique& r);
  NeighUnique& operator=(const NeighUnique& r);
  virtual ~NeighUnique();

  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  virtual ENeigh getType() const override { return ENeigh::UNIQUE; }
  virtual int getMaxSampleNumber(const Db* db) const override;

  static NeighUnique* create(int ndim, bool flag_xvalid = false);
  static NeighUnique* createFromNF(const String& neutralFilename, bool verbose = true);

  int reset(int ndim, bool flag_xvalid = false);

protected:
  /// Interface for ASerializable
  virtual bool _deserialize(std::istream& is, bool verbose = false) override;
  virtual bool _serialize(std::ostream& os, bool verbose = false) const override;
  String _getNFName() const override { return "NeighUnique"; }
};
