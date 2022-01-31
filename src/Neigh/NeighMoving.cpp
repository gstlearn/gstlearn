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
#include "Neigh/NeighMoving.hpp"
#include "Morpho/Morpho.hpp"
#include "Basic/Utilities.hpp"
#include "Basic/AException.hpp"
#include "Basic/Vector.hpp"
#include "Db/Db.hpp"
#include "geoslib_f.h"
#include "geoslib_old_f.h"

NeighMoving::NeighMoving(int ndim, bool flag_xvalid)
    : ANeighParam(ndim, flag_xvalid),
      _flagSector(0),
      _flagAniso(0),
      _flagRotation(0),
      _nMini(1),
      _nMaxi(0),
      _nSect(1),
      _nSMax(0),
      _radius(0.),
      _anisoCoeffs(),
      _anisoRotMat()
{
}

NeighMoving::NeighMoving(const NeighMoving& r)
    : ANeighParam(r),
      _flagSector(r._flagSector),
      _flagAniso(r._flagAniso),
      _flagRotation(r._flagRotation),
      _nMini(r._nMini),
      _nMaxi(r._nMaxi),
      _nSect(r._nSect),
      _nSMax(r._nSMax),
      _radius(r._radius),
      _anisoCoeffs(r._anisoCoeffs),
      _anisoRotMat(r._anisoRotMat)
{
}

NeighMoving& NeighMoving::operator=(const NeighMoving& r)
{
  if (this != &r)
  {
    ANeighParam::operator=(r);
    _flagSector = r._flagSector;
    _flagAniso = r._flagAniso;
    _flagRotation = r._flagRotation;
    _nMini = r._nMini;
    _nMaxi = r._nMaxi;
    _nSect = r._nSect;
    _nSMax = r._nSMax;
    _radius = r._radius;
    _anisoCoeffs = r._anisoCoeffs;
    _anisoRotMat = r._anisoRotMat;
   }
  return *this;
}

NeighMoving::~NeighMoving()
{
}

int NeighMoving::reset(int ndim,
                       bool flag_xvalid,
                       int nmaxi,
                       double radius,
                       int nmini,
                       int nsect,
                       int nsmax,
                       VectorDouble coeffs,
                       VectorDouble angles)
{
  setNDim(ndim);
  setFlagXvalid(flag_xvalid);

  _nMini = nmini;
  _nMaxi = nmaxi;
  _nSect = nsect;
  _nSMax = nsmax;
  _radius = radius;

  if (! coeffs.empty())
  {
    _flagAniso = 1;
    _anisoCoeffs = coeffs;

    if (! angles.empty())
    {
      _flagRotation = 1;
      _anisoRotMat.resize(ndim * ndim);
      ut_rotation_matrix(ndim, angles.data(), _anisoRotMat.data());
    }
  }
  return 0;
}

String NeighMoving::toString(const AStringFormat* strfmt) const
{
  std::stringstream sstr;
  int ndim = getNDim();

  sstr << ANeighParam::toString(strfmt);

  if (_nMini > 0)
    sstr << "Minimum number of samples           = " << _nMini << std::endl;
  if (_nMaxi > 0)
    sstr << "Maximum number of samples           = " << _nMaxi << std::endl;
  if (_flagSector)
  {
    sstr << "Number of angular sectors           = " << _nSect << std::endl;
    if (_nSMax > 0)
      sstr << "Maximum number of points per sector = " << _nSMax << std::endl;
  }
  if (getFlagContinuous())
  {
    sstr << "Norm. dist. for continuous NeighMoving.   = " << getDistCont()
         << std::endl;
  }

  if (!_flagAniso)
  {
    if (!FFFF(_radius))
      sstr << "Maximum horizontal distance         = " << _radius << std::endl;
  }
  else
  {
    VectorDouble ranges(ndim);
    for (int idim = 0; idim < ndim; idim++)
      ranges[idim] = _radius * _anisoCoeffs[idim];
    sstr << toMatrix("Anisotropic Ranges :", VectorString(), VectorString(),
                    true, ndim, 1, ranges);

    if (_flagRotation)
    {
      sstr << toMatrix("Anisotropy Rotation :", VectorString(), VectorString(),
                      true, ndim, ndim, _anisoRotMat);
    }
  }

  /* Cross-validation option */

  if (getFlagXvalid() != 0)
    sstr << "The Cross-Validation Option is switched ON" << std::endl;

  return sstr.str();
}

int NeighMoving::_deserialize(FILE* file, bool verbose)
{
  if (ANeighParam::_deserialize(file, verbose))
  {
    if (verbose)
      messerr("Problem reading from the Neutral File.");
    return 1;
  }

  int ndim = getNDim();
  VectorDouble radius(ndim);
  VectorDouble nbgh_coeffs(ndim);
  VectorDouble nbgh_rotmat(ndim * ndim);

  int flag_aniso = 0;
  int flag_rotation = 0;
  int flag_sector = 0;
  double dmax = 0.;

  if (_recordRead(file, "NeighMovingborhood sector search", "%d", &flag_sector))
    return 1;
  if (_recordRead(file, "Minimum Number of samples", "%d", &_nMini)) return 1;
  if (_recordRead(file, "Maximum Number of samples", "%d", &_nMaxi)) return 1;
  if (_recordRead(file, "Optimum Number of samples per sector", "%d", &_nSect))
    return 1;
  if (_recordRead(file, "Maximum Number of samples per sector", "%d", &_nSMax))
    return 1;
  if (_recordRead(file, "Maximum Isotropic Radius", "%lf", &dmax)) return 1;
  if (_recordRead(file, "Flag for Anisotropy", "%d", &flag_aniso)) return 1;
  if (flag_aniso)
  {
    for (int idim = 0; idim < ndim; idim++)
      if (_recordRead(file, "Anisotropy Coefficient", "%lf",
                      &nbgh_coeffs[idim])) return 1;
    if (_recordRead(file, "Flag for Anisotropy Rotation", "%d", &flag_rotation))
      return 1;
    if (flag_rotation)
    {
      int lec = 0;
      for (int idim = 0; idim < ndim; idim++)
        for (int jdim = 0; jdim < ndim; jdim++, lec++)
          if (_recordRead(file, "Anisotropy Rotation Matrix", "%lf",
                          &nbgh_rotmat[lec])) return 1;
    }
  }
  if (!nbgh_coeffs.empty())
    for (int idim = 0; idim < ndim; idim++)
      nbgh_coeffs[idim] *= dmax;

  setNSect((flag_sector) ? MAX(_nSect, 1) : 1);
  setRadius(dmax);
  setFlagSector(flag_sector && ndim >= 2);
  setFlagAniso(flag_aniso && !nbgh_coeffs.empty());
  setFlagRotation(flag_rotation && flag_aniso && !nbgh_rotmat.empty());

  if (getFlagAniso() && !nbgh_coeffs.empty())
  {
    setRadius(0.);
    for (int i = 0; i < getNDim(); i++)
      setRadius(MAX(getRadius(), nbgh_coeffs[i]));
    for (int i = 0; i < getNDim(); i++)
      setAnisoCoeff(i, nbgh_coeffs[i] / getRadius());
  }
  if (getFlagRotation() && !nbgh_rotmat.empty())
  {
    setAnisoRotMat(nbgh_rotmat);
  }

  return 0;
}

int NeighMoving::_serialize(FILE* file, bool verbose) const
{
  if (ANeighParam::_serialize(file, verbose))
    {
      if (verbose) messerr("Problem writing in the Neutral File.");
      return 1;
    }

  _recordWrite(file, "%d", getFlagSector());
  _recordWrite(file, "#", "Use angular sectors");
  _recordWrite(file, "%d", getNMini());
  _recordWrite(file, "%d", getNMaxi());
  _recordWrite(file, "%d", getNSect());
  _recordWrite(file, "%d", getNSMax());
  _recordWrite(file, "#", "Parameters (nmini,nmaxi,nsect,nsmax)");
  _recordWrite(file, "%lf", getRadius());
  _recordWrite(file, "#", "Maximum distance radius");
  _recordWrite(file, "%d", getFlagAniso());
  _recordWrite(file, "#", "Anisotropy Flag");

  if (getFlagAniso())
  {
    for (int idim = 0; idim < getNDim(); idim++)
      _recordWrite(file, "%lf", getAnisoCoeff(idim));
    _recordWrite(file, "#", "Anisotropy Coefficients");
    _recordWrite(file, "%d", getFlagRotation());
    _recordWrite(file, "#", "Anisotropy Rotation Flag");
    if (getFlagRotation())
    {
      int ecr = 0;
      for (int idim = 0; idim < getNDim(); idim++)
        for (int jdim = 0; jdim < getNDim(); jdim++)
          _recordWrite(file, "%lf", getAnisoRotMat(ecr++));
      _recordWrite(file, "#", "Anisotropy Rotation Matrix");
    }
  }
  return 0;
}

void NeighMoving::setAnisoCoeff(int idim, double value)
{
  if ((int) _anisoCoeffs.size() != getNDim())
    _anisoCoeffs.resize(getNDim(),1.);
  _anisoCoeffs[idim] = value;
}

void NeighMoving::anisoRescale()
{
  for (int idim = 0; idim < getNDim(); idim++)
    _anisoCoeffs[idim] /= _radius;
}

NeighMoving* NeighMoving::create(int ndim,
                                 bool flag_xvalid,
                                 int nmaxi,
                                 double radius,
                                 int nmini,
                                 int nsect,
                                 int nsmax,
                                 VectorDouble coeffs,
                                 VectorDouble angles)
{
  NeighMoving* neighM = new NeighMoving;
  if (neighM->reset(ndim, flag_xvalid, nmaxi, radius, nmini, nsect, nsmax,
                    coeffs, angles))
  {
    messerr("Problem when creating Moving NeighMovingborhood");
    delete neighM;
    neighM = nullptr;
  }
  return neighM;
}

int NeighMoving::dumpToNF(const String& neutralFilename, bool verbose) const
{
  FILE* file = _fileOpen(neutralFilename, "NeighMoving", "w", verbose);
  if (file == nullptr) return 1;

  if (_serialize(file, verbose))
  {
    if (verbose) messerr("Problem writing in the Neutral File.");
    _fileClose(file, verbose);
    return 1;
  }
  _fileClose(file, verbose);
  return 0;
}

/**
 * Create a NeighMovingborhood by loading the contents of a Neutral File
 * @param neutralFilename Name of the Neutral File
 * @param verbose         Verbose flag
 * @return
 */
NeighMoving* NeighMoving::createFromNF(const String& neutralFilename, bool verbose)
{
  FILE* file = _fileOpen(neutralFilename, "NeighMoving", "r", verbose);
  if (file == nullptr) return nullptr;

  NeighMoving* neigh = new NeighMoving();
  if (neigh->_deserialize(file, verbose))
  {
    if (verbose) messerr("Problem reading the Neutral File.");
    delete neigh;
    neigh = nullptr;
  }
  _fileClose(file, verbose);
  return neigh;
}

/**
 * Given a Db, returns the maximum number of samples per NeighMovingborhood
 * @param db Pointer to the taregt Db
 * @return
 */
int NeighMoving::getMaxSampleNumber(const Db* /*db*/) const
{
  return (_flagSector) ? _nSect * _nSMax : _nMaxi;
}
