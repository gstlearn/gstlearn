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
#include "geoslib_define.h"

#include "Neigh/ENeigh.hpp"
#include "Neigh/ANeighParam.hpp"

#include "Basic/AStringable.hpp"
#include "Basic/ASerializable.hpp"
#include "Basic/IClonable.hpp"

class Db;

class GSTLEARN_EXPORT NeighMoving: public ANeighParam
{
public:
  NeighMoving(int ndim = 2, bool flag_xvalid = false);
  NeighMoving(const NeighMoving& r);
  NeighMoving& operator=(const NeighMoving& r);
  virtual ~NeighMoving();

  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  virtual int getMaxSampleNumber(const Db* db) const override;
  virtual ENeigh getType() const override { return ENeigh::MOVING; }

  int reset(int ndim,
            bool flag_xvalid,
            int nmaxi,
            double radius,
            int nmini = 1,
            int nsect = 1,
            int nsmax = ITEST,
            VectorDouble coeffs = VectorDouble(),
            VectorDouble angles = VectorDouble());

  int dumpToNF(const String& neutralFilename, bool verbose = false) const;
  static NeighMoving* create(int ndim,
                             bool flag_xvalid,
                             int nmaxi,
                             double radius = TEST,
                             int nmini = 1,
                             int nsect = 1,
                             int nsmax = ITEST,
                             VectorDouble coeffs = VectorDouble(),
                             VectorDouble angles = VectorDouble());
  static NeighMoving* createFromNF(const String& neutralFilename, bool verbose = false);

  const VectorDouble& getAnisoCoeffs() const { return _anisoCoeffs; }
  double getAnisoCoeff(int i) const { return _anisoCoeffs[i]; }
  const VectorDouble& getAnisoRotMats() const { return _anisoRotMat; }
  double getAnisoRotMat(int i) const { return _anisoRotMat[i]; }
  int getFlagAniso() const { return _flagAniso; }
  int getFlagRotation() const { return _flagRotation; }
  int getFlagSector() const { return _flagSector; }
  int getNMaxi() const { return _nMaxi; }
  int getNMini() const { return _nMini; }
  int getNSect() const { return _nSect; }
  int getNSMax() const { return _nSMax; }
  double getRadius() const { return _radius; }

  void setAnisoCoeffs(const VectorDouble& anisoCoeffs) { _anisoCoeffs = anisoCoeffs; }
  void setAnisoCoeff(int idim, double value);
  void anisoRescale();
  void setAnisoRotMat(const VectorDouble& anisoRotMat) { _anisoRotMat = anisoRotMat; }
  void setFlagAniso(int flagAniso) { _flagAniso = flagAniso; }
  void setFlagRotation(int flagRotation) { _flagRotation = flagRotation; }
  void setFlagSector(int flagSector) { _flagSector = flagSector; }
  void setNMaxi(int nmaxi) { _nMaxi = nmaxi; }
  void setNMini(int nmini) { _nMini = nmini; }
  void setNSect(int nsect) { _nSect = nsect; }
  void setNSMax(int nsmax) { _nSMax = nsmax; }
  void setRadius(double radius) { _radius = radius; }

protected:
  virtual int _deserialize(std::istream& is, bool verbose = false) override;
  virtual int _serialize(std::ostream& os, bool verbose = false) const override;

private:
  int _flagSector;               /* 1 if MOVING neigh. used sector search */
  int _flagAniso;                /* 1 if the MOVING neigh. is anisotropic */
  int _flagRotation;             /* 1 if the anisotropy is rotated */
  int _nMini;                    /* Minimum number of points in neigh. */
  int _nMaxi;                    /* Maximum number of points in neigh. */
  int _nSect;                    /* Number of 2-D angular sectors */
  int _nSMax;                    /* Maximum number of points per 2-D sector */
  double _radius;                /* Maximum isotropic distance */
  VectorDouble _anisoCoeffs;     /* Anisotropy ratio for MOVING neigh. */
  VectorDouble _anisoRotMat;     /* Anisotropy rotation matrix */
};
