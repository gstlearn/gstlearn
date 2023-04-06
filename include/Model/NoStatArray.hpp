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

#include "Matrix/MatrixRectangular.hpp"
#include "Mesh/AMesh.hpp"
#include "Basic/VectorNumT.hpp"
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

  /// ICloneable interface
  IMPLEMENT_CLONING(NoStatArray)

  /// AStringable Interface
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  double getValue(int igrf,
                  int icov,
                  const EConsElem &type,
                  int iv1,
                  int iv2,
                  int icas,
                  int rank) const override;
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
                      const VectorVectorDouble& coords,
                      VectorDouble& tab,
                      bool verbose) const;
  String _displayStats(int ipar, int icas) const;
  String _displayStats(int icas) const;

private:
  const Db* _dbnostat;
  mutable MatrixRectangular _tab; // Dimension: nvertex * npar
};
