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

#include "Basic/VectorNumT.hpp"
#include "Space/ASpaceObject.hpp"
#include "Space/SpacePoint.hpp"

class GSTLEARN_EXPORT SpaceTarget : public SpacePoint
{
public:
  SpaceTarget(const ASpace* space = nullptr,
              bool checkExtend    = true,
              bool checkCode      = true,
              bool checkDate      = true);
  SpaceTarget(const SpaceTarget& r);
  SpaceTarget& operator=(const SpaceTarget& r);
  virtual ~SpaceTarget();

  static SpaceTarget* create(const VectorDouble &center = VectorDouble(),
                             const VectorDouble &extend = VectorDouble(),
                             double code = TEST,
                             double date = TEST,
                             const ASpace *space = nullptr);

//  SpacePoint& getCoordAsSP() { return _center; }
  const SpacePoint& getCoordAsSP() const { return *this; }

  const VectorDouble& getExtend() const { return _extend; }
  double getExtend(int idim) const { return _extend[idim]; }
  void setExtend(const VectorDouble& extend) { _extend = extend; }
  void setExtend(int i, double val){ _extend[i] = val; }
  const double* getExtendP() const {return _extend.data(); }
  double* getExtendPM() {return _extend.data(); }
  void setCode(double code) { _code = code; }
  void setDate(double date) { _date = date; }
  double getCode() const { return _code; }
  double getDate() const { return _date; }
  bool checkExtend() const { return _checkExtend; }
  bool checkCode() const { return _checkCode; }
  bool checkDate() const { return _checkDate; }

  /// Convert space point to string
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

private:
  void _initialize();

protected:
  bool _checkExtend;
  bool _checkCode;
  bool _checkDate;
  VectorDouble _extend;
  double _code;
  double _date;
};
