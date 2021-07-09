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

#include "Basic/AStringable.hpp"
#include "Basic/ASerializable.hpp"
#include "Basic/Vector.hpp"

class Neigh: public AStringable , ASerializable
{
public:
  Neigh();
  Neigh(int ndim);
  Neigh(int ndim, int nmaxi, double radius,
        int nmini=1, int nsect=1, int nsmax=ITEST, double width=0, double distcont=0,
        VectorDouble coeffs = VectorDouble(), VectorDouble angles = VectorDouble());
  Neigh(int ndim, int skip, VectorDouble image);
  Neigh(const String& neutralFileName, bool verbose);
  Neigh(const Neigh& r);
  Neigh& operator=(const Neigh& r);
  virtual ~Neigh();

  virtual std::string toString(int level = 0) const override;
  int deSerialize(const String& filename, bool verbose = false) override;
  int serialize(const String& filename, bool verbose = false) override;

  const VectorDouble& getAnisoCoeff() const { return _anisoCoeffs; }
  double getAnisoCoeff(int i) const { return _anisoCoeffs[i]; }
  const VectorDouble& getAnisoRotMat() const { return _anisoRotMat; }
  double getAnisoRotMat(int i) const { return _anisoRotMat[i]; }
  double getDistCont() const { return _distCont; }
  int getFlagAniso() const { return _flagAniso; }
  int getFlagContinuous() const { return _flagContinuous; }
  int getFlagRotation() const { return _flagRotation; }
  int getFlagSector() const { return _flagSector; }
  int getFlagXvalid() const { return _flagXvalid; }
  const VectorDouble& getImageRadius() const { return _imageRadius; }
  double getImageRadius(int idim) const { return _imageRadius[idim]; }
  int getNDim() const { return _nDim; }
  int getNMaxi() const { return _nMaxi; }
  int getNMini() const { return _nMini; }
  int getNSect() const { return _nSect; }
  int getNSMax() const { return _nSMax; }
  double getRadius() const { return _radius; }
  int getSkip() const { return _skip; }
  int getType() const { return _type; }
  double getWidth() const { return _width; }

  void setAnisoCoeff(const VectorDouble& anisoCoeffs) { _anisoCoeffs = anisoCoeffs; }
  void setAnisoCoeff(int idim, double value) { _anisoCoeffs[idim] = value; }
  void setAnisoRotMat(const VectorDouble& anisoRotMat) { _anisoRotMat = anisoRotMat; }
  void setDistCont(double distCont) { _distCont = distCont; }
  void setFlagAniso(int flagAniso) { _flagAniso = flagAniso; }
  void setFlagContinuous(int flagContinuous) { _flagContinuous = flagContinuous; }
  void setFlagRotation(int flagRotation) { _flagRotation = flagRotation; }
  void setFlagSector(int flagSector) { _flagSector = flagSector; }
  void setFlagXvalid(int flagXvalid) { _flagXvalid = flagXvalid; }
  void setImageRadius(const VectorDouble& imageRadius) { _imageRadius = imageRadius; }
  void setNDim(int ndim)   { _nDim = ndim; }
  void setNMaxi(int nmaxi) { _nMaxi = nmaxi; }
  void setNMini(int nmini) { _nMini = nmini; }
  void setNSect(int nsect) { _nSect = nsect; }
  void setNSMax(int nsmax) { _nSMax = nsmax; }
  void setRadius(double radius) { _radius = radius; }
  void setSkip(int skip) { _skip = skip; }
  void setType(int type) { _type = type; }
  void setWidth(double width) { _width = width; }

private:
  bool _isDimensionValid(int idim) const;
  void _init(int ndim,
             int type,
             int flag_xvalid,
             int flag_sector,
             int flag_aniso,
             int flag_rotation,
             int flag_continuous,
             int nmini,
             int nmaxi,
             int nsect,
             int nsmax,
             int skip,
             double width,
             double radius,
             double dist_cont,
             const VectorDouble& nbgh_radius,
             const VectorDouble& nbgh_rotmat,
             const VectorDouble& nbgh_image);

public:
  int _nDim;                     /* Space dimension */
  int _type;                     /* Neighborhood type: NEIGH_* */
  int _flagXvalid;               /* 1 to suppress the target */
  int _flagSector;               /* 1 if MOVING neigh. used sector search */
  int _flagAniso;                /* 1 if the MOVING neigh. is anisotropic */
  int _flagRotation;             /* 1 if the anisotropy is rotated */
  int _flagContinuous;           /* 1 for continuous moving neighborhood */
  int _nMini;                    /* Minimum number of points in neigh. */
  int _nMaxi;                    /* Maximum number of points in neigh. */
  int _nSect;                    /* Number of 2-D angular sectors */
  int _nSMax;                    /* Maximum number of points per 2-D sector */
  int _skip;                     /* Skipping factor */
  double _width;                 /* Width of the slice - bench */
  double _radius;                /* Maximum isotropic distance */
  double _distCont;              /* Distance for continuous neighborhood */
  VectorDouble _anisoCoeffs;     /* Anisotropy ratio for MOVING neigh. */
  VectorDouble _anisoRotMat;     /* Anisotropy rotation matrix */
  VectorDouble _imageRadius;     /* Vector of image neighborhood radius */
};
