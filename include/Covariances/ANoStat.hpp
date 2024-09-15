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

#include "Basic/VectorNumT.hpp"
#include "Mesh/AMesh.hpp"
#include "gstlearn_export.hpp"
#include "Basic/AStringable.hpp"


class ACov;
class AMesh;

class GSTLEARN_EXPORT ANoStat : public AStringable
{
public:
  ANoStat();
  ANoStat(const ANoStat &m) = delete;
  double getValueOnDbOut(int iech) const;
  double getValueOnDbIn(int iech) const;
  double getValueOnDb(int iech,int icas) const;
  bool   getValuesOnDb(int icas1, int iech1,double* val1, 
                       int icas2, int iech2, double* val2) const;
  double getValueOnMeshByMesh(int imesh) const;
  double getValueOnMeshByApex(int iapex) const;
  double getValueOnMesh(int iapex,bool center = false) const;
  void informField(const VectorVectorDouble & coords, VectorDouble& tab, bool verbose = false);
  void informMeshByMesh(const AMesh* amesh, bool verbose = false);
  void informMeshByApex(const AMesh* amesh, bool verbose = false);
  void informDbIn(const Db* dbin, bool verbose = false);
  void informDbOout(const Db* dbout, bool verbose = false);

  ANoStat& operator= (const ANoStat &m) = delete;
  virtual ~ANoStat();

private:
  void _informDb(const Db* db, VectorDouble &res, bool verbose = false);
  virtual void _informField(const VectorVectorDouble &coords,
                            VectorDouble& res,
                            bool verbose = false) = 0;

protected:
  bool _isValid(int icas, int rank) const;
  mutable VectorDouble _tabdbin;
  mutable VectorDouble _tabdbout;
  mutable VectorDouble _tabmesh; // Dimension: nmeshes
  mutable VectorDouble _tabvertices; // Dimension: nvertex

};
