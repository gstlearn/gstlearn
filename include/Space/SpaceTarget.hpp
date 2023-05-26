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

#include "Basic/VectorNumT.hpp"
#include "Space/ASpaceObject.hpp"
#include "Space/SpacePoint.hpp"

class GSTLEARN_EXPORT SpaceTarget : public ASpaceObject
{
public:
  SpaceTarget(const ASpace* space = nullptr);

  SpaceTarget(const SpaceTarget& r);
  SpaceTarget& operator=(const SpaceTarget& r);
  virtual ~SpaceTarget();

  static SpaceTarget* create(const VectorDouble &center = VectorDouble(),
                             const VectorDouble &extend = VectorDouble(),
                             const ASpace *space = nullptr);

  SpacePoint& getCoordAsSP() { return _center; }
  const SpacePoint& getCoordAsSPP() const { return _center; }
  const VectorDouble& getCoord() const { return _center.getCoord(); }
  double getCoord(int idim) const { return _center.getCoord(idim); }
  const VectorDouble& getExtend() const { return _extend; }
  double getExtend(int idim) const { return _extend[idim]; }

  void setCoord(const VectorDouble& center) { _center = center; };
  void setCoord(int i, double val){_center.setCoord(i, val); }
  const double* getCoordsP() const {return _center.getCoordsP(); }
  double* getCoordsPM() {return _center.getCoordsPM(); }

  void setExtend(const VectorDouble& extend) { _extend = extend; }
  void setExtend(int i, double val){ _extend[i] = val; }
  const double* getExtendP() const {return _extend.data(); }
  double* getExtendPM() {return _extend.data(); }

  /// Return true if the point is consistent with the provided space
  virtual bool isConsistent(const ASpace* space) const override;

  /// Convert space point to string
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

private:
  void _initialize();

protected:
  SpacePoint   _center;
  VectorDouble _extend;
};
