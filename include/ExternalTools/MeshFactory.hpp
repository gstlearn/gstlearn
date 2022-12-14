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

#include "Basic/VectorNumT.hpp"

class AMesh;
class MatrixRectangular;
class MatrixInt;

class GSTLEARN_EXPORT MeshFactory
{

public:
  MeshFactory();
  virtual ~MeshFactory();
  /*! Returns a pointer to the optimally create Meshing */

  static AMesh *createMesh(int variety,
                           const VectorDouble& extendmin,
                           const VectorDouble& extendmax,
                           const VectorDouble& cellsize,
                           const VectorDouble& rotmat,
                           const VectorDouble& extperc,
                           Db *dbin,
                           Db *dbout,
                           const String& triswitch,
                           MatrixRectangular& apices,
                           MatrixInt& meshes,
                           bool flag_polarize,
                           int verbose = 0);
};
