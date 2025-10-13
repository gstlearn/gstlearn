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

#include "Basic/AStringable.hpp"
#include "gstlearn_export.hpp"
#include "LinearOp/IProj.hpp"
#include <vector>

namespace gstlrn
{
class GSTLEARN_EXPORT ProjMulti : public IProj, public AStringable
{
public:
  ProjMulti(const std::vector<std::vector<const IProj*>> &projs, bool silent = false);

  /// AStringable Interface
  String toString(const AStringFormat* strfmt = nullptr) const override;
  
  Id getNApex() const override;
  Id getNPoint() const override;
  Id getNVariable() const { return _nvariable; }
  Id getNLatent() const { return _nlatent; }
  virtual ~ProjMulti();
  bool empty() const { return _projs.empty();}

#ifndef SWIG

protected:
  Id _addPoint2mesh(const constvect inv, vect outv) const override;
  Id _addMesh2point(const constvect inv, vect outv) const override;
#endif

private:
  bool _checkArg(const std::vector<std::vector<const IProj*>>& projs) const;
  void _init();
  virtual void _clear() {};

protected:
  Id findFirstNoNullOnRow(Id j) const;
  Id findFirstNoNullOnCol(Id j) const;
  const std::vector<Id>& getNPoints() const { return _pointNumbers; }
  const std::vector<Id>& getNApexs() const { return _apexNumbers; }

protected:
std::vector<std::vector<const IProj*>> _projs; // NOT TO BE DELETED

private:
Id _pointNumber;
Id _apexNumber;
Id _nlatent;
Id _nvariable;
std::vector<Id> _pointNumbers;
std::vector<Id> _apexNumbers;
bool _silent;
mutable std::vector<double> _work;
mutable std::vector<double> _workmesh;
};
}