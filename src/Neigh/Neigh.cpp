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
#include "Morpho/Morpho.hpp"
#include "Basic/Utilities.hpp"
#include "geoslib_f.h"

Neigh::Neigh()
    : _nDim(0),
      _type(0),
      _flagXvalid(0),
      _flagSector(0),
      _flagAniso(0),
      _flagRotation(0),
      _flagContinuous(0),
      _nMini(1),
      _nMaxi(0),
      _nSect(1),
      _nSMax(0),
      _skip(0),
      _width(0.),
      _radius(0.),
      _distCont(0.),
      _anisoCoeffs(),
      _anisoRotMat(),
      _imageRadius()
{
}

Neigh::Neigh(int ndim)
    : _nDim(ndim),
      _type(NEIGH_UNIQUE),
      _flagXvalid(0),
      _flagSector(0),
      _flagAniso(0),
      _flagRotation(0),
      _flagContinuous(0),
      _nMini(0),
      _nMaxi(0),
      _nSect(1),
      _nSMax(0),
      _skip(0),
      _width(0.),
      _radius(0.),
      _distCont(0.),
      _anisoCoeffs(),
      _anisoRotMat(),
      _imageRadius()
{
}

Neigh::Neigh(int ndim, int nmaxi, double radius, int nmini, int nsect,
             int nsmax, double width, double distcont, VectorDouble coeffs,
             VectorDouble angles)
    : _nDim(ndim),
      _type(NEIGH_MOVING),
      _flagXvalid(0),
      _flagSector(0),
      _flagAniso(0),
      _flagRotation(0),
      _flagContinuous(0),
      _nMini(nmini),
      _nMaxi(nmaxi),
      _nSect(nsect),
      _nSMax(nsmax),
      _skip(0),
      _width(width),
      _radius(radius),
      _distCont(distcont),
      _anisoCoeffs(),
      _anisoRotMat(),
      _imageRadius()
{
  VectorDouble nbgh_radius;
  VectorDouble nbgh_rotmat;

  nbgh_rotmat.resize(ndim * ndim);

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
}

Neigh::Neigh(int ndim, int skip, VectorDouble image)
    : _nDim(ndim),
      _type(NEIGH_IMAGE),
      _flagXvalid(0),
      _flagSector(0),
      _flagAniso(0),
      _flagRotation(0),
      _flagContinuous(0),
      _nMini(0),
      _nMaxi(0),
      _nSect(0),
      _nSMax(0),
      _skip(skip),
      _width(0.),
      _radius(0.),
      _distCont(0.),
      _anisoCoeffs(),
      _anisoRotMat(),
      _imageRadius(image)
{

}

Neigh::Neigh(const Neigh& r)
    : _nDim(r._nDim),
      _type(r._type),
      _flagXvalid(r._flagXvalid),
      _flagSector(r._flagSector),
      _flagAniso(r._flagAniso),
      _flagRotation(r._flagRotation),
      _flagContinuous(r._flagContinuous),
      _nMini(r._nMini),
      _nMaxi(r._nMaxi),
      _nSect(r._nSect),
      _nSMax(r._nSMax),
      _skip(r._skip),
      _width(r._width),
      _radius(r._radius),
      _distCont(r._distCont),
      _anisoCoeffs(r._anisoCoeffs),
      _anisoRotMat(r._anisoRotMat),
      _imageRadius(r._imageRadius)
{
}

Neigh& Neigh::operator=(const Neigh& r)
{
  if (this != &r)
  {
    _nDim = r._nDim;
    _type = r._type;
    _flagXvalid = r._flagXvalid;
    _flagSector = r._flagSector;
    _flagAniso = r._flagAniso;
    _flagRotation = r._flagRotation;
    _flagContinuous = r._flagContinuous;
    _nMini = r._nMini;
    _nMaxi = r._nMaxi;
    _nSect = r._nSect;
    _nSMax = r._nSMax;
    _skip = r._skip;
    _width = r._width;
    _radius = r._radius;
    _distCont = r._distCont;
    _anisoCoeffs = r._anisoCoeffs;
    _anisoRotMat = r._anisoRotMat;
    _imageRadius = r._imageRadius;
   }
  return *this;
}

Neigh::~Neigh()
{
}

std::string Neigh::toString(int level) const
{
  std::stringstream sstr;

  /* Neighborhood options */

  sstr << toTitle(0,"Neighborhood characteristics");

  switch (_type)
  {
    case NEIGH_UNIQUE:
      sstr << "Unique neighborhood option" << std::endl;
      sstr << "Space dimension = " << _nDim << std::endl;
      break;

    case NEIGH_BENCH:
      sstr << "Bench neighborhood option" << std::endl;
      sstr << "Space dimension = " << _nDim << std::endl;
      sstr << "Bench width     = " << _width << std::endl;
      break;

    case NEIGH_MOVING:
      sstr << "Moving neighborhood option" << std::endl;
      sstr << "Space dimension                     = " << _nDim << std::endl;
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
      if (_flagContinuous)
      {
        sstr << "Norm. dist. for continuous Neigh.   = " << _distCont << std::endl;
      }

      if (! _flagAniso)
      {
        if (! FFFF(_radius))
          sstr << "Maximum horizontal distance         = " << _radius << std::endl;
      }
      else
      {
        VectorDouble ranges(_nDim);
        for (int idim=0; idim<_nDim; idim++)
          ranges[idim] = _radius * _anisoCoeffs[idim];
        sstr << toMatrix("Anisotropic Ranges :",VectorString(), VectorString(),
                         true,_nDim,1,ranges);

        if (_flagRotation)
        {
          sstr << toMatrix("Anisotropy Rotation :",VectorString(),VectorString(),
                           true, _nDim,_nDim,_anisoRotMat);
        }
      }
      break;

    case NEIGH_IMAGE:
      sstr << "Image neighborhood option" << std::endl;
      sstr << "Skipping factor = " << _skip << std::endl;
      sstr << toMatrix("Image radius :",VectorString(),VectorString(),
                       true,_nDim,1,_imageRadius);
      break;
  }

  /* Cross-validation option */

  if (_flagXvalid != 0)
    sstr << "The Cross-Validation Option is switched ON" << std::endl;

  return sstr.str();
}

bool Neigh::_isDimensionValid(int idim) const
{
  if (idim < 0 || idim >= getNDim())
  {
    messerr("Error in 'idim'(%d). It should lie within [0,%d[",idim,getNDim());
    return false;
  }
  return true;
}
