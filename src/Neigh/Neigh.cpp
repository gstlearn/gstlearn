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
#include "Neigh/Neigh.hpp"
#include "Morpho/Morpho.hpp"
#include "Basic/Utilities.hpp"
#include "Basic/AException.hpp"
#include "geoslib_f.h"
#include "geoslib_old_f.h"

Neigh::Neigh()
    : AStringable(),
      ASerializable(),
      _nDim(0),
      _type(),
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

/**
 * Constructor of a Unique Neighborhood
 */
Neigh::Neigh(int ndim)
    : AStringable(),
      ASerializable(),
      _nDim(ndim),
      _type(ENeigh::UNIQUE),
      _flagXvalid(0),
      _flagSector(0),
      _flagAniso(0),
      _flagRotation(0),
      _flagContinuous(0),
      _nMini(1), // Put 1 for consistency with default construtor
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
    : AStringable(),
      ASerializable(),
      _nDim(ndim),
      _type(ENeigh::MOVING),
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

Neigh::Neigh(int ndim, int skip, const VectorInt& image)
    : AStringable(),
      ASerializable(),
      _nDim(ndim),
      _type(ENeigh::IMAGE),
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

Neigh::Neigh(const String& neutralFileName, bool verbose)
    : AStringable(),
      ASerializable(),
      _nDim(0),
      _type(ENeigh::MOVING),
      _flagXvalid(0),
      _flagSector(0),
      _flagAniso(0),
      _flagRotation(0),
      _flagContinuous(0),
      _nMini(0),
      _nMaxi(0),
      _nSect(0),
      _nSMax(0),
      _skip(0),
      _width(0.),
      _radius(0.),
      _distCont(0.),
      _anisoCoeffs(),
      _anisoRotMat(),
      _imageRadius()
{
  if (deSerialize(neutralFileName, verbose))
  {
    messerr("Error when reading a Neutral File");
    messerr("The Neigh is not entirely defined");
  }
}

Neigh::Neigh(const Neigh& r)
    : AStringable(),
      ASerializable(),
      _nDim(r._nDim),
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

String Neigh::toString(int /*level*/) const
{
  std::stringstream sstr;

  /* Neighborhood options */

  sstr << toTitle(0,"Neighborhood characteristics");

  switch (_type.toEnum())
  {
    case ENeigh::E_UNIQUE:
      sstr << "Unique neighborhood option" << std::endl;
      sstr << "Space dimension = " << _nDim << std::endl;
      break;

    case ENeigh::E_BENCH:
      sstr << "Bench neighborhood option" << std::endl;
      sstr << "Space dimension = " << _nDim << std::endl;
      sstr << "Bench width     = " << _width << std::endl;
      break;

    case ENeigh::E_MOVING:
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

    case ENeigh::E_IMAGE:
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

int Neigh::deSerialize(const String& filename, bool verbose)
{
  int type, idim, ndim, flag_sector, flag_xvalid, nmini, nmaxi, nsect, nsmax, skip;
  int flag_aniso, flag_rotation, lec, jdim;
  double width, dmax;
  VectorInt radius;
  VectorDouble nbgh_coeffs;
  VectorDouble nbgh_rotmat;

  /* Opening the Data file */

  if (_fileOpen(filename, "Neigh", "r", verbose)) return 1;

  /* Create the Model structure */

  if (_recordRead("Space Dimension", "%d", &ndim)) return 1;
  if (_recordRead("Neighborhood type", "%d", &type)) return 1;

  /* Core allocation */

  radius.resize(ndim);
  nbgh_coeffs.resize(ndim);
  nbgh_rotmat.resize(ndim * ndim);

  switch (type) // Real integer read from the file
  {
    case ENeigh::E_UNIQUE:
      _init(ndim, ENeigh::UNIQUE, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0., 0., 0.,
            VectorDouble(), VectorDouble(), VectorInt());
      break;

    case ENeigh::E_BENCH:
      if (_recordRead("Flag for Cross-validation", "%d", &flag_xvalid)) return 1;
      if (_recordRead("Bench Width", "%lf", &width)) return 1;
      _init(ndim, ENeigh::BENCH, flag_xvalid, 0, 0, 0, 0, 0, 0, 0, 0, 0, width,
            0., 0., VectorDouble(), VectorDouble(), VectorInt());
      break;

    case ENeigh::E_MOVING:
      flag_aniso = flag_rotation = 0;
      if (_recordRead("Flag for Cross-validaiton", "%d", &flag_xvalid)) return 1;
      if (_recordRead("Neighborhood sector search", "%d", &flag_sector)) return 1;
      if (_recordRead("Neighborhood Width", "%lf", &width)) return 1;
      if (_recordRead("Minimum Number of samples", "%d", &nmini)) return 1;
      if (_recordRead("Maximum Number of samples", "%d", &nmaxi)) return 1;
      if (_recordRead("Optimum Number of samples per sector", "%d", &nsect)) return 1;
      if (_recordRead("Maximum Number of samples per sector", "%d", &nsmax)) return 1;
      if (_recordRead("Maximum Isotropic Radius", "%lf", &dmax)) return 1;
      if (_recordRead("Flag for Anisotropy", "%d", &flag_aniso)) return 1;
      if (flag_aniso)
      {
        for (idim = 0; idim < ndim; idim++)
          if (_recordRead("Anisotropy Coefficient", "%lf", &nbgh_coeffs[idim])) return 1;
        if (_recordRead("Flag for Anisotropy Rotation", "%d", &flag_rotation)) return 1;
        if (flag_rotation)
        {
          for (idim = lec = 0; idim < ndim; idim++)
            for (jdim = 0; jdim < ndim; jdim++, lec++)
              if (_recordRead("Anisotropy Rotation Matrix", "%lf", &nbgh_rotmat[lec])) return 1;
        }
      }
      if (!nbgh_coeffs.empty()) for (idim = 0; idim < ndim; idim++) nbgh_coeffs[idim] *= dmax;

      _init(ndim, ENeigh::MOVING, flag_xvalid, flag_sector, flag_aniso,
            flag_rotation, 0, nmini, nmaxi, nsect, nsmax, 0, 0., dmax, 0.,
            nbgh_coeffs, nbgh_rotmat, VectorInt());
      break;

    case ENeigh::E_IMAGE:
      if (_recordRead("Flag for Cross-Validation", "%d", &flag_xvalid)) return 1;
      if (_recordRead("Skipping factor", "%d", &skip)) return 1;
      for (idim = 0; idim < ndim; idim++)
      {
        double loc_radius;
        if (_recordRead("Image Neighborhood Radius", "%lf", &loc_radius)) return 1;
        radius[idim] = static_cast<int> (loc_radius);
      }
      _init(ndim, ENeigh::IMAGE, flag_xvalid, 0, 0, 0, 0, 0, 0, 0, 0, skip, 0.,
            0., 0., VectorDouble(), VectorDouble(), radius);
  }

  _fileClose(verbose);

  return 0;
}

void Neigh::_init(int ndim,
                  ENeigh type,
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
                  const VectorInt& nbgh_image)
{
  /// TODO : Force SpaceRN creation (deSerialization do not know yet how to manage other space types)
  //SpaceRN space(ndim);
  // To be used when Neigh object will manage space context (and not only ndim)

  setNDim(ndim);
  setType(type);
  setNMini(nmini);
  setNMaxi(nmaxi);
  setNSect((flag_sector) ? MAX(nsect, 1) :
                           1);
  setNSMax(nsmax);
  setWidth(width);
  setRadius(radius);
  setDistCont(dist_cont);
  setSkip(skip);
  setFlagXvalid(flag_xvalid);
  setFlagSector(flag_sector && ndim >= 2);
  setFlagAniso(flag_aniso && !nbgh_radius.empty());
  setFlagRotation(flag_rotation && flag_aniso && !nbgh_rotmat.empty());
  setFlagContinuous((!IFFFF(flag_continuous)) ? flag_continuous : 0);

  /* Core allocation */

  if (getFlagAniso() && !nbgh_radius.empty())
  {
    setRadius(0.);
    for (int i = 0; i < getNDim(); i++)
      setRadius(MAX(getRadius(), nbgh_radius[i]));
    for (int i = 0; i < getNDim(); i++)
      setAnisoCoeff(i, nbgh_radius[i] / getRadius());
  }
  if (getFlagRotation() && !nbgh_rotmat.empty())
  {
    setAnisoRotMat(nbgh_rotmat);
  }
  if (type == ENeigh::IMAGE && !nbgh_image.empty())
  {
    setImageRadius(nbgh_image);
  }
}

int Neigh::serialize(const String& filename, bool verbose) const
{
  int ecr;
  if (_fileOpen(filename, "Neigh", "w", verbose)) return 1;

  /* Create the Model structure */

  _recordWrite("%d", getNDim());
  _recordWrite("#", "Space Dimension");
  _recordWrite("%d", getType().getValue());
  _recordWrite("#", "Neighborhood type");

  switch (getType().toEnum())
  {
    case ENeigh::E_UNIQUE:
      break;

    case ENeigh::E_BENCH:
      _recordWrite("%d", getFlagXvalid());
      _recordWrite("#", "Cross-Validation flag");
      _recordWrite("%lf", getWidth());
      _recordWrite("#", "Bench Width");
      break;

    case ENeigh::E_MOVING:
      _recordWrite("%d", getFlagXvalid());
      _recordWrite("#", "Cross-Validation flag");
      _recordWrite("%d", getFlagSector());
      _recordWrite("#", "Use angular sectors");
      _recordWrite("%lf", getWidth());
      _recordWrite("#", "Bench Width");
      _recordWrite("%d", getNMini());
      _recordWrite("%d", getNMaxi());
      _recordWrite("%d", getNSect());
      _recordWrite("%d", getNSMax());
      _recordWrite("#", "Parameters (nmini,nmaxi,nsect,nsmax)");
      _recordWrite("%lf", getRadius());
      _recordWrite("#", "Maximum distance radius");
      _recordWrite("%d", getFlagAniso());
      _recordWrite("#", "Anisotropy Flag");
      if (!getFlagAniso()) break;

      for (int idim = 0; idim < getNDim(); idim++)
        _recordWrite("%lf", getAnisoCoeff(idim));
      _recordWrite("#", "Anisotropy Coefficients");
      _recordWrite("%d", getFlagRotation());
      _recordWrite("#", "Anisotropy Rotation Flag");
      if (!getFlagRotation()) break;

      ecr = 0;
      for (int idim = 0; idim < getNDim(); idim++)
        for (int jdim = 0; jdim < getNDim(); jdim++)
          _recordWrite("%lf", getAnisoRotMat(ecr++));
      _recordWrite("#", "Anisotropy Rotation Matrix");
      break;

    case ENeigh::E_IMAGE:
      _recordWrite("%d", getFlagXvalid());
      _recordWrite("#", "Cross-Validation flag");
      _recordWrite("%d", getSkip());
      for (int idim = 0; idim < getNDim(); idim++)
        _recordWrite("%lf", (double) getImageRadius(idim));
      _recordWrite("#", "Image neighborhood parameters");
      break;
  }

  _fileClose(verbose);

  return 0;
}

void Neigh::setAnisoCoeff(int idim, double value)
{
  if ((int) _anisoCoeffs.size() != _nDim)
    _anisoCoeffs.resize(_nDim,1.);
  _anisoCoeffs[idim] = value;
}

void Neigh::anisoRescale()
{
  for (int idim = 0; idim < _nDim; idim++)
    _anisoCoeffs[idim] /= _radius;
}
