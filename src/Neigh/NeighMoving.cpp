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
#include "geoslib_old_f.h"

#include "Neigh/NeighMoving.hpp"
#include "Morpho/Morpho.hpp"
#include "Basic/Utilities.hpp"
#include "Basic/AException.hpp"
#include "Basic/VectorHelper.hpp"
#include "Faults/Faults.hpp"
#include "Db/Db.hpp"

#include <math.h>

NeighMoving::NeighMoving(bool flag_xvalid,
                         int nmaxi,
                         double radius,
                         int nmini,
                         int nsect,
                         int nsmax,
                         VectorDouble coeffs,
                         VectorDouble angles,
                         double distcont,
                         const Faults *faults,
                         const ASpace *space)
    : ANeighParam(flag_xvalid, space),
      _flagAniso(0),
      _flagRotation(0),
      _nMini(nmini),
      _nMaxi(nmaxi),
      _nSect(nsect),
      _nSMax(nsmax),
      _forceWithinBlock(false),
      _radius(radius),
      _distCont(distcont),
      _anisoCoeffs(),
      _anisoRotMat(),
      _faults(faults)
{
  _radius = radius;
  _distCont = distcont;

  if (! coeffs.empty())
  {
    //    _flagAniso = (ut_vector_constant(coeffs)) ? 0 : 1;
    _flagAniso = 1;
    _anisoCoeffs = coeffs;

    if (! angles.empty())
    {
      int ndim = getNDim();
      _flagRotation = (VH::isConstant(angles, 0.)) ? 0 : 1;
      _anisoRotMat.resize(ndim * ndim);
      ut_rotation_matrix(ndim, angles.data(), _anisoRotMat.data());
    }
  }
}

NeighMoving::NeighMoving(const NeighMoving& r)
    : ANeighParam(r),
      _flagAniso(r._flagAniso),
      _flagRotation(r._flagRotation),
      _nMini(r._nMini),
      _nMaxi(r._nMaxi),
      _nSect(r._nSect),
      _nSMax(r._nSMax),
      _forceWithinBlock(r._forceWithinBlock),
      _radius(r._radius),
      _distCont(r._distCont),
      _anisoCoeffs(r._anisoCoeffs),
      _anisoRotMat(r._anisoRotMat),
      _faults(r._faults)
{
}

NeighMoving& NeighMoving::operator=(const NeighMoving& r)
{
  if (this != &r)
  {
    ANeighParam::operator=(r);
    _flagAniso = r._flagAniso;
    _flagRotation = r._flagRotation;
    _nMini = r._nMini;
    _nMaxi = r._nMaxi;
    _nSect = r._nSect;
    _nSMax = r._nSMax;
    _forceWithinBlock = r._forceWithinBlock;
    _radius = r._radius;
    _distCont = r._distCont;
    _anisoCoeffs = r._anisoCoeffs;
    _anisoRotMat = r._anisoRotMat;
    _faults = r._faults;
   }
  return *this;
}

NeighMoving::~NeighMoving()
{
}

void NeighMoving::setForceWithinBlock(bool forceWithinBlock)
{
  _forceWithinBlock = forceWithinBlock;

  // In the case of forcing data selection per cell, cancel the other options
  if (_forceWithinBlock)
  {
    _flagAniso = false;
    _flagRotation = false;
    _nSect = 1;
    _nSMax = ITEST;
  }
}

String NeighMoving::toString(const AStringFormat* strfmt) const
{
  std::stringstream sstr;
  int ndim = getNDim();

  sstr << toTitle(0,"Moving Neighborhood");
  sstr << ANeighParam::toString(strfmt);

  if (_nMini > 0)
    sstr << "Minimum number of samples           = " << _nMini << std::endl;

  if (_forceWithinBlock)
  {
    sstr << "Force Selection of all samples within target Block" << std::endl;
  }
  else
  {
    if (_nMaxi > 0)
      sstr << "Maximum number of samples           = " << _nMaxi << std::endl;

    if (_nSect > 1)
    {
      sstr << "Number of angular sectors           = " << _nSect << std::endl;
      if (_nSMax > 0)
        sstr << "Maximum number of points per sector = " << _nSMax << std::endl;
    }
    if (!FFFF(_radius))
    {
      if (!_flagAniso)
      {
        sstr << "Maximum horizontal distance         = " << _radius
             << std::endl;
      }
      else
      {
        VectorDouble ranges(ndim);
        for (int idim = 0; idim < ndim; idim++)
          ranges[idim] = _radius * _anisoCoeffs[idim];
        sstr
            << toMatrix("Anisotropic Ranges :", VectorString(), VectorString(),
                        true, ndim, 1, ranges);

        if (_flagRotation)
        {
          sstr
              << toMatrix("Anisotropy Rotation :", VectorString(),
                          VectorString(), true, ndim, ndim, _anisoRotMat);
        }
      }
    }
  }
  if (getFlagContinuous())
  {
    sstr << "Norm. dist. for continuous NeighMoving.   = " << getDistCont()
         << std::endl;
  }
  return sstr.str();
}

bool NeighMoving::_deserialize(std::istream& is, bool verbose)
{
  bool ret = true;
  ret = ret && ANeighParam::_deserialize(is, verbose);
  if (! ret) return ret;

  int ndim = getNDim();
  VectorDouble radius(ndim);
  VectorDouble nbgh_coeffs(ndim);
  VectorDouble nbgh_rotmat(ndim * ndim);

  int flag_aniso = 0;
  int flag_rotation = 0;
  int flag_sector = 0;
  double dmax = 0.;

  ret = true;
  ret = ret && _recordRead<int>(is, "NeighMovingborhood sector search", flag_sector);
  ret = ret && _recordRead<int>(is, "Minimum Number of samples", _nMini);
  ret = ret && _recordRead<int>(is, "Maximum Number of samples", _nMaxi);
  ret = ret && _recordRead<int>(is, "Optimum Number of samples per sector", _nSect);
  ret = ret && _recordRead<int>(is, "Maximum Number of samples per sector", _nSMax);
  ret = ret && _recordRead<double>(is, "Maximum Isotropic Radius", dmax);
  ret = ret && _recordRead<int>(is, "Flag for Anisotropy", flag_aniso);
  if (flag_aniso)
  {
    for (int idim = 0; ret && idim < ndim; idim++)
      ret = ret && _recordRead<double>(is, "Anisotropy Coefficient", nbgh_coeffs[idim]);
    ret = ret && _recordRead<int>(is, "Flag for Anisotropy Rotation", flag_rotation);
    if (flag_rotation)
    {
      int lec = 0;
      for (int idim = 0; ret && idim < ndim; idim++)
        for (int jdim = 0; ret && jdim < ndim; jdim++, lec++)
          ret = ret && _recordRead<double>(is, "Anisotropy Rotation Matrix", nbgh_rotmat[lec]);
    }
  }
  if (! ret) return ret;

  if (!nbgh_coeffs.empty())
    for (int idim = 0; idim < ndim; idim++)
      nbgh_coeffs[idim] *= dmax;

  setNSect((getFlagSector()) ? MAX(_nSect, 1) : 1);
  setRadius(dmax);
  setFlagAniso(flag_aniso && !nbgh_coeffs.empty());
  setFlagRotation(flag_rotation && flag_aniso && !nbgh_rotmat.empty());

  if (getFlagAniso() && !nbgh_coeffs.empty())
  {
    setRadius(0.);
    for (int i = 0; i < (int) getNDim(); i++)
      setRadius(MAX(getRadius(), nbgh_coeffs[i]));
    for (int i = 0; i < (int) getNDim(); i++)
      setAnisoCoeff(i, nbgh_coeffs[i] / getRadius());
  }
  if (getFlagRotation() && !nbgh_rotmat.empty())
  {
    setAnisoRotMat(nbgh_rotmat);
  }
  return ret;
}

bool NeighMoving::_serialize(std::ostream& os, bool verbose) const
{
  bool ret = true;
  ret = ret && ANeighParam::_serialize(os, verbose);
  ret = ret && _recordWrite<int>(os, "Use angular sectors", getFlagSector());
  ret = ret && _recordWrite<int>(os, "", getNMini());
  ret = ret && _recordWrite<int>(os, "", getNMaxi());
  ret = ret && _recordWrite<int>(os, "", getNSect());
  ret = ret && _recordWrite<int>(os, "", getNSMax());
  ret = ret && _commentWrite(os, "Parameters (nmini,nmaxi,nsect,nsmax)");
  ret = ret && _recordWrite<double>(os, "Maximum distance radius", getRadius());
  ret = ret && _recordWrite<int>(os, "Anisotropy Flag", getFlagAniso());

  if (getFlagAniso())
  {
    for (int idim = 0; ret && idim < (int) getNDim(); idim++)
      ret = ret && _recordWrite<double>(os, "", getAnisoCoeff(idim));
    ret = ret && _commentWrite(os, "Anisotropy Coefficients");
    ret = ret && _recordWrite<int>(os, "Anisotropy Rotation Flag", getFlagRotation());
    if (getFlagRotation())
    {
      int ecr = 0;
      for (int idim = 0; ret && idim < (int) getNDim(); idim++)
        for (int jdim = 0; ret && jdim < (int) getNDim(); jdim++)
          ret = ret && _recordWrite<double>(os, "", getAnisoRotMat(ecr++));
      ret = ret && _commentWrite(os, "Anisotropy Rotation Matrix");
    }
  }
  return ret;
}

void NeighMoving::setAnisoCoeff(int idim, double value)
{
  if (_anisoCoeffs.size() != getNDim())
    _anisoCoeffs.resize(getNDim(),1.);
  _anisoCoeffs[idim] = value;
}

void NeighMoving::anisoRescale()
{
  if (FFFF(_radius)) return;
  for (int idim = 0; idim < (int) getNDim(); idim++)
    _anisoCoeffs[idim] /= _radius;
}

NeighMoving* NeighMoving::create(bool flag_xvalid,
                                 int nmaxi,
                                 double radius,
                                 int nmini,
                                 int nsect,
                                 int nsmax,
                                 VectorDouble coeffs,
                                 VectorDouble angles,
                                 double distcont,
                                 const Faults* faults,
                                 const ASpace* space)
{
  return new NeighMoving(flag_xvalid, nmaxi, radius, nmini, nsect, nsmax,
                         coeffs, angles, distcont, faults, space);
}

/**
 * Create a NeighMovingborhood by loading the contents of a Neutral File
 * @param neutralFilename Name of the Neutral File
 * @param verbose         Verbose flag
 * @return
 */
NeighMoving* NeighMoving::createFromNF(const String& neutralFilename, bool verbose)
{
  NeighMoving* neigh = nullptr;
  std::ifstream is;
  neigh = new NeighMoving();
  bool success = false;
  if (neigh->_fileOpenRead(neutralFilename, is, verbose))
  {
    success =  neigh->deserialize(is, verbose);
  }
  if (! success)
  {
    delete neigh;
    neigh = nullptr;
  }
  return neigh;
}

/**
 * Given a Db, returns the maximum number of samples per NeighMovingborhood
 * @param db Pointer to the target Db
 * @return
 */
int NeighMoving::getMaxSampleNumber(const Db* /*db*/) const
{
  return (getFlagSector()) ? _nSect * _nSMax : _nMaxi;
}

bool NeighMoving::getFlagSector() const
{
  return (getNDim() > 1 && _nSect > 1);
}

// TODO: add the rotation and possible ellipse
VectorVectorDouble NeighMoving::getEllipsoid(const VectorDouble& target, int count) const
{
  if (_forceWithinBlock) return VectorVectorDouble();
  VectorVectorDouble coords(2);
  coords[0].resize(count+1,0.);
  coords[1].resize(count+1,0.);

  for (int i = 0; i < count; i++)
  {
    double angle = 2. * i * GV_PI / count;
    coords[0][i] = target[0] + _radius * cos(angle);
    coords[1][i] = target[1] + _radius * sin(angle);
  }

  coords[0][count] = target[0] + _radius;
  coords[1][count] = target[1];
  return coords;
}

/**
 * Generate the end-points of the sectors. By default, the extension is set to radius
 * @param target Coordinates of the Target
 * @return
 */
VectorVectorDouble NeighMoving::getSectors(const VectorDouble& target) const
{
  if (_forceWithinBlock) return VectorVectorDouble();
  VectorVectorDouble coords(2);
  coords[0].resize(_nSect,0.);
  coords[1].resize(_nSect,0.);

  for (int i = 0; i < _nSect; i++)
  {
    double angle = 2. * i * GV_PI / _nSect;
    coords[0][i] = target[0] + _radius * cos(angle);
    coords[1][i] = target[1] + _radius * sin(angle);
  }
  return coords;
}
