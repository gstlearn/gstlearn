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

#include "Neigh/ENeigh.hpp"
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

  int getNDim() const { return _nDim; }
  double getDistCont() const { return _distCont; }
  int getFlagContinuous() const { return _flagContinuous; }
  int getFlagXvalid() const { return _flagXvalid; }

  void setDistCont(double distCont) { _distCont = distCont; }
  void setFlagXvalid(int flagXvalid) { _flagXvalid = flagXvalid; }
  void setNDim(int dim) { _nDim = dim; }
  void setFlagContinuous(int flagContinuous) { _flagContinuous = flagContinuous; }

protected:
  // ASerializable Interface overriding
  virtual int _deserialize(FILE* file, bool verbose = false);
  virtual int _serialize(FILE* file, bool verbose = false) const override;

  virtual int _deserialize2(std::istream& is, bool verbose = false) override;

  bool _isDimensionValid(int idim) const;

private:
  int _nDim;                     /* Space dimension */
  int _flagXvalid;               /* 1 to suppress the target */
  int _flagContinuous;           /* 1 for continuous moving neighborhood */
  double _distCont;              /* Distance for continuous ANeighParamborhood */
};
