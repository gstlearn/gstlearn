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

#include "Matrix/MatrixSparse.hpp"
#include "LinearOp/ProjMulti.hpp"
#include "Matrix/MatrixT.hpp"

class ProjMatrix;
class AMesh;
class Db;

class GSTLEARN_EXPORT ProjMultiMatrix : public ProjMulti
{
public:
  ProjMultiMatrix(const std::vector<std::vector<const ProjMatrix*>> &proj,bool toClean = false,bool silent = false);
  virtual ~ProjMultiMatrix();
  static std::vector<std::vector<const ProjMatrix*>> create(std::vector<const ProjMatrix*> &vectproj, int nvariable);
  static ProjMultiMatrix createFromDbAndMeshes(const Db* db,const std::vector<const AMesh*> &meshes,bool verbose = false);

  const MatrixSparse* getProj() const { return &_Proj;} 
#ifndef SWIG           
  protected:
    virtual int _addPoint2mesh(const constvect inv, vect outv) const override;
    virtual int _addMesh2point(const constvect inv, vect outv) const override;
#endif
private:
  MatrixSparse  _Proj;
  void _clear() override;
private:
  bool _toClean;
};
