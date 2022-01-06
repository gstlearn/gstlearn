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
#include "Matrix/MatrixRectangular.hpp"
#include "Mesh/AMesh.hpp"
#include "Basic/Vector.hpp"
#include "Basic/AStringable.hpp"
#include "Model/ANoStat.hpp"

class GSTLEARN_EXPORT NoStatArray : public ANoStat
{
public:
	NoStatArray();
	NoStatArray(const VectorString& codes, const Db* dbnostat);
  NoStatArray(const NoStatArray &m);
  NoStatArray& operator=(const NoStatArray &m);
  virtual ~NoStatArray();

  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  virtual IClonable* clone() const override { return new NoStatArray(*this); };

  double getValue(int igrf, int icov, const EConsElem& type, int iv1, int iv2,
                  int icas, int rank) const override;
  double getValueByParam(int ipar, int icas, int rank) const override;

  int  attachToMesh(const AMesh* mesh, bool verbose = false) const override;
  void detachFromMesh() const override;
  int  attachToDb(Db* db, int icas, bool verbose = false) const override;
  void detachFromDb(Db* db, int icas) const override;

  bool   isEmpty(int icas) const;

private:
  bool   _checkValid() const;
  int    _getNpoints() const { return _tab.getNRows(); }
  double _interpolate(int ipar, int icas1, int iech1, int icas2, int iech2) const;
  int    _informField(int ipar,
                      int nech,
                      double* coor[3],
                      VectorDouble& tab,
                      bool verbose) const;
  String _displayStats(int ipar, int icas) const;
  String _displayStats(int icas) const;

private:
  const Db* _dbnostat;
	mutable MatrixRectangular _tab; // Dimension: nvertex * npar
};
