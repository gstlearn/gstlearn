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
#include "LinearOp/IProj.hpp"
#include <vector>

class GSTLEARN_EXPORT ProjMulti : public IProj
{
public:
  ProjMulti(const std::vector<std::vector<const IProj*>> &projs,bool silent = false);
  int getNApex() const override;
  int getNPoint() const override;
  int getNVariable() const { return _nvariable; }
  int getNLatent() const { return _nlatent; }
  virtual ~ProjMulti();
  bool empty() const { return _projs.empty();}

#ifndef SWIG           
  protected:
    virtual int _addPoint2mesh(const constvect inv, vect outv) const override;
    virtual int _addMesh2point(const constvect inv, vect outv) const override;
#endif

private : 
  bool _checkArg(const std::vector<std::vector<const IProj*>> &projs) const;
  void _init();
  virtual void _clear(){};
protected:
int findFirstNoNullOnRow(int j) const;
int findFirstNoNullOnCol(int j) const;
const std::vector<int>& getNPoints() const {return _pointNumbers;}
const std::vector<int>& getNApexs()  const {return _apexNumbers;}

protected:
std::vector<std::vector<const IProj*>> _projs; // NOT TO BE DELETED

private:
int _pointNumber;
int _apexNumber;
int _nlatent;
int _nvariable;
std::vector<int> _pointNumbers;
std::vector<int> _apexNumbers;
bool _silent;
mutable std::vector<double> _work;
mutable std::vector<double> _workmesh;
};
