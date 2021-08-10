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

#include "MatrixC/MatrixCRectangular.hpp"
#include "Mesh/AMesh.hpp"
#include "Basic/Vector.hpp"
#include "Basic/AStringable.hpp"
#include "Model/ANoStat.hpp"

class NoStatArray : public ANoStat
{
public:
	NoStatArray();
	NoStatArray(const VectorString& codes, const Db* dbnostat);
  NoStatArray(const NoStatArray &m);
  NoStatArray& operator=(const NoStatArray &m);
  virtual ~NoStatArray();

  virtual String toString(int level = 0) const override;

  double getValue(int igrf, int icov, ENUM_CONS type, int iv1, int iv2, int icas,
                  int rank) const override;
  double getValue(int ipar, int icas, int rank) const override;

  int attachToMesh(const AMesh* mesh, bool verbose = false) const;
  void detachFromMesh(bool verbose) const;
  int attachToDb(Db* db, bool verbose = false) const;
  void detachFromDb(Db* db, bool verbose) const;
  String displayStats(int ipar, int icas) const;
  String displayStats(int icas) const;

  void   updateModel(Model* model, int iech1, int iech2) const;
  void   updateModel(Model* model, int vertex) const;
  bool   isEmpty(int ipar, int icas) const;

private:
  bool _checkValid() const;
  int  _getNpoints() const { return _tab.getNRows(); }
  void _getInfoFromDb(int ipar,
                      int iech1,
                      int iech2,
                      double *val1,
                      double *val2) const;
  double _interpolate(int ipar, int iech1, int iech2) const;
  int _informField(int ipar,
                   int nech,
                   double* coor[3],
                   VectorDouble& tab,
                   bool verbose) const;

private:
  const Db* _dbnostat;
  mutable const Db* _dbin;
  mutable VectorInt _attIn;
  mutable const Db* _dbout;
  mutable VectorInt _attOut;
	mutable MatrixCRectangular _tab; // Dimension: nvertex * npar
};
