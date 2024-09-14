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

#include "Matrix/MatrixRectangular.hpp"
#include "Mesh/AMesh.hpp"
#include "Basic/VectorNumT.hpp"
#include "Covariances/ANoStatCov.hpp"
class GSTLEARN_EXPORT NoStatArrayCov : public ANoStatCov
{
public:
  NoStatArrayCov();
  NoStatArrayCov(const VectorString& codes, const Db* dbnostat);
  NoStatArrayCov(const NoStatArrayCov &m);
  NoStatArrayCov& operator=(const NoStatArrayCov &m);
  virtual ~NoStatArrayCov();

  /// ICloneable interface
  IMPLEMENT_CLONING(NoStatArrayCov)

  /// AStringable Interface
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  double getValue(const EConsElem &type,
                  int icas,
                  int rank,
                  int iv1 = -1,
                  int iv2 = -1) const override;
  double getValueByParam(int ipar, int icas, int rank) const override;

  int  attachToMesh(const AMesh* mesh, bool center = true, bool verbose = false) const override;
  void detachFromMesh() const override;
  int  attachToDb(Db* db, int icas, bool verbose = false) const override;
  void detachFromDb(Db* db, int icas) const override;

  bool   isEmpty(int icas) const;

private:
  bool   _checkValid() const;
  int    _getNpoints() const { return _tab.getNRows(); }
  int    _informField(int ipar,
                      const VectorVectorDouble& coords,
                      VectorDouble& tab,
                      bool verbose) const;
  String _displayStats(int ipar, int icas) const;
  String _displayStats(int icas) const;

private:
  const Db* _dbnostat;
  mutable MatrixRectangular _tab; // Dimension: nvertex * npar
};
