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

#include "Basic/AStringable.hpp"
#include "Basic/ASerializable.hpp"
#include "Space/ASpaceObject.hpp"

class Db;

// TODO : Inherits from ASpaceParam which inherits from ASPaceObject and AParam, which inherits from ASerializable, AStringable, IClonable
class GSTLEARN_EXPORT ANeighParam: public ASpaceObject, public ASerializable
{
public:
  ANeighParam(bool flag_xvalid = false, const ASpace* space = nullptr);
  ANeighParam(const ANeighParam& r);
  ANeighParam& operator=(const ANeighParam& r);
  virtual ~ANeighParam();

  /// ASpaceObject Interface
  virtual bool isConsistent(const ASpace* space) const override;

  // AStringable Interface overriding
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  // ANeighParam Interface definition
  virtual int getMaxSampleNumber(const Db* db) const = 0;
  virtual ENeigh getType() const = 0;
  virtual bool getFlagContinuous() const { return false; }

  bool getFlagXvalid() const { return _flagXvalid; }
  bool getFlagKFold() const { return _flagKFold; }

  void setFlagXvalid(bool flagXvalid) { _flagXvalid = flagXvalid; }
  void setFlagKFold(bool flagKFold) { _flagKFold = flagKFold; }

  VectorInt eval(Db *dbin, Db *dbout, int iech0) const;

protected:
  // Interface for ASerializable
  virtual bool _deserialize(std::istream& is, bool verbose = false) override;
  virtual bool _serialize(std::ostream& os, bool verbose = false) const override;
  String _getNFName() const override { return "ANeighParam"; }

private:
  bool _isDimensionValid(int idim) const;

private:
  bool _flagXvalid;              /* True to suppress the target */
  bool _flagKFold;               /* True to perform a KFold Cross-validation */
};
