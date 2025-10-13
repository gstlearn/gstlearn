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

#include "Enum/EDirGen.hpp"
#include "Enum/EGaussInv.hpp"

#include "Anamorphosis/AnamHermite.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/ICloneable.hpp"
#include "Basic/NamingConvention.hpp"
#include "Basic/VectorNumT.hpp"

namespace gstlrn
{

class Db;
class MatrixSymmetric;
class MatrixDense;
class AMatrix;

class GSTLEARN_EXPORT PPMT: public AStringable, public ICloneable
{
public:
  PPMT(Id ndir                     = 50,
       bool flagPreprocessing       = false,
       const EDirGen& methodDir     = EDirGen::fromKey("VDC"),
       const EGaussInv& methodTrans = EGaussInv::fromKey("EMP"),
       Id nbpoly                   = 30,
       double alpha                 = 2.);
  PPMT(const PPMT& m);
  PPMT& operator=(const PPMT& m);
  virtual ~PPMT();

  /// AStringable Interface
  String toString(const AStringFormat* strfmt = nullptr) const override;

  /// ICloneable Interface
  IMPLEMENT_CLONING(PPMT)

  static PPMT* create(Id ndir                     = 50,
                      bool flagPreprocessing       = false,
                      const EDirGen& methodDir     = EDirGen::fromKey("VDC"),
                      const EGaussInv& methodTrans = EGaussInv::fromKey("EMP"),
                      Id nbpoly                   = 30,
                      double alpha                 = 2.);

  Id getNiter() const { return _niter; }
  double getAlpha() const { return _alpha; }
  Id getNdir() const { return _ndir; }
  Id getNdim() const { return _ndim; }
  Id getNbpoly() const { return _nbpoly; }
  bool isFitted() const { return _isFitted; }
  const EDirGen& getMethodDir() const { return _methodDir; }
  const EGaussInv& getMethodTrans() const { return _methodTrans; }

  VectorDouble getSerieAngle() const { return _serieAngle; }
  VectorDouble getSerieScore(bool flagLog = false) const;

  Id fit(Db* db,
          const VectorString& names,
          bool flagStoreInDb              = false,
          Id niter                       = 100,
          bool verbose                    = false,
          const NamingConvention& namconv = NamingConvention("Y"));
  Id fitFromMatrix(AMatrix* Y, Id niter, bool verbose = false);
  Id rawToGaussian(Db* db,
                    const VectorString& names,
                    Id niter                       = 0,
                    const NamingConvention& namconv = NamingConvention("Y"));
  Id gaussianToRaw(Db* db,
                    const VectorString& names,
                    Id niter                       = 0,
                    const NamingConvention& namconv = NamingConvention("Z"));

private:
  void _generateAllDirections();
  void _fitInitHermite(AMatrix* Y);
  void _initGaussianizeForward(AMatrix* Y);
  void _initGaussianizeBackward(AMatrix* Y);
  void _iterationFit(AMatrix* Y, const VectorDouble& N0);
  void _iterationForward(AMatrix* Y, const VectorDouble& N0, Id iter = 0);
  void _iterationBackward(AMatrix* Y, const VectorDouble& N0, Id iter = 0);
  static double _gaussianizeForward(double Yi,
                                    Id rank,
                                    const AnamHermite* anam,
                                    const VectorDouble& N0);
  static double _gaussianizeBackward(double Yi, const AnamHermite* anam);
  void _projectOnDirection(const AMatrix* Y, Id id, VectorDouble& Y0);
  double _getGaussianDistance(const VectorDouble& Yi,
                              const VectorInt& Ri,
                              const VectorDouble& N0) const;
  void _shiftForward(AMatrix* Y,
                     Id id,
                     const AnamHermite* anam,
                     const VectorDouble& Y0,
                     const VectorInt& R0,
                     const VectorDouble& N0) const;
  void _shiftBackward(AMatrix* Y,
                      Id id,
                      const AnamHermite* anam,
                      const VectorDouble& Y0) const;

private:
  Id _niter;
  Id _ndir;
  Id _nbpoly;
  double _alpha;
  EDirGen _methodDir;
  EGaussInv _methodTrans;
  bool _flagPreprocessing;

  bool _isFitted;

  mutable Id _ndim;
  mutable VectorDouble _serieAngle;
  mutable VectorDouble _serieScore;
  mutable MatrixDense* _dirmat;
  mutable std::vector<AnamHermite*> _anams;
  mutable std::vector<AnamHermite*> _initAnams;
  mutable MatrixDense* _initSphering;
};

} // namespace gstlrn