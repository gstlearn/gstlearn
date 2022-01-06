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
#include "IProjMatrix.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/Vector.hpp"
#include "csparse_d.h"

class AMesh;
class Db;
class SPDE_Mesh;

class GSTLEARN_EXPORT ProjMatrix: public IProjMatrix, public AStringable
{
public:
  ProjMatrix();
  ProjMatrix(const Db* db, AMesh *a_mesh, int verbose = 0);
  ProjMatrix(int npoint, int napices, const cs *aproj);
  ProjMatrix(const ProjMatrix &m);
  ProjMatrix& operator= (const ProjMatrix &m);
  virtual ~ProjMatrix();

  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  int resetFromDb(const Db* db, AMesh *a_mesh, int verbose = 0);
  int resetFromPoints(int npoint, int napices, const cs *aproj);
  int resetFromDbOldStyle(Db* db, SPDE_Mesh* s_mesh, int verbose = 0);
  int resetFromDbByNeighOldStyle(const Db* db,
                                 SPDE_Mesh* s_mesh,
                                 double radius,
                                 int flag_exact = 0,
                                 int verbose = 0);

  int point2mesh(const VectorDouble& in, VectorDouble& out) const override;
  int mesh2point(const VectorDouble& in, VectorDouble& out) const override;
  int getApexNumber() const override { return _nApices; }
  int getPointNumber() const override { return _nPoint; }
  const cs* getAproj() const { return _Aproj; }

private:
  int  _nPoint;
  int  _nApices;
  cs*  _Aproj;
};
