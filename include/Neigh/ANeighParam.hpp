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

class Db;

class GSTLEARN_EXPORT ANeighParam: public AStringable, public ASerializable
{
public:
  ANeighParam(int ndim = 2, bool flag_xvalid = false);
  ANeighParam(const ANeighParam& r);
  ANeighParam& operator=(const ANeighParam& r);
  virtual ~ANeighParam();

  // AStringable Interface overriding
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  // ANeighParam Interface definition
  virtual int getMaxSampleNumber(const Db* db) const = 0;
  virtual ENeigh getType() const = 0;
  virtual bool getFlagContinuous() const { return false; }
  virtual bool hasFault() const { return false; }

  int getNDim() const { return _nDim; }
  bool getFlagXvalid() const { return _flagXvalid; }
  bool getFlagKFold() const { return _flagKFold; }

  void setFlagXvalid(bool flagXvalid) { _flagXvalid = flagXvalid; }
  void setNDim(int dim) { _nDim = dim; }
  void setFlagKFold(bool flagKFold) { _flagKFold = flagKFold; }

protected:
  // Interface for ASerializable
  virtual bool _deserialize(std::istream& is, bool verbose = false) override;
  virtual bool _serialize(std::ostream& os, bool verbose = false) const override;
  String _getNFName() const override { return "ANeighParam"; }

private:
  bool _isDimensionValid(int idim) const;

private:
  int  _nDim;                    /* Space dimension */
  bool _flagXvalid;              /* True to suppress the target */
  bool _flagKFold;               /* True to perform a KFold Cross-validation */
};
