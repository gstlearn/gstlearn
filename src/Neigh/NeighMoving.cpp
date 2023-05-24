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
#include "geoslib_old_f.h"

#include "Neigh/NeighMoving.hpp"
#include "Morpho/Morpho.hpp"
#include "Basic/OptDbg.hpp"
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
    : ANeigh(nullptr, nullptr, space),
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
      _faults(faults),
      _movingInd(),
      _movingIsect(),
      _movingNsect(),
      _movingX1(),
      _movingX2(),
      _movingDst()
{
  setFlagXvalid(flag_xvalid);
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
      GH::rotationInit(ndim, angles.data(), _anisoRotMat.data());
    }
  }
}

NeighMoving::NeighMoving(const NeighMoving& r)
    : ANeigh(r),
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
      _faults(r._faults),
      _movingInd(r._movingInd),
      _movingIsect(r._movingIsect),
      _movingNsect(r._movingNsect),
      _movingX1(r._movingX1),
      _movingX2(r._movingX2),
      _movingDst(r._movingDst)
{
}

NeighMoving& NeighMoving::operator=(const NeighMoving& r)
{
  if (this != &r)
  {
    ANeigh::operator=(r);
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
    _movingInd = r._movingInd;
    _movingIsect = r._movingIsect;
    _movingNsect = r._movingNsect;
    _movingX1 = r._movingX1;
    _movingX2 = r._movingX2;
    _movingDst = r._movingDst;
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
  sstr << ANeigh::toString(strfmt);

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
  ret = ret && ANeigh::_deserialize(is, verbose);
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
  ret = ret && ANeigh::_serialize(os, verbose);
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

VectorVectorDouble NeighMoving::getZoomLimits(const VectorDouble& target, double percent) const
{
  VectorVectorDouble coords(2);
  coords[0].resize(2, 0.);
  coords[1].resize(2, 0.);

  coords[0][0] = target[0] - _radius * (1. + percent / 100.);
  coords[0][1] = target[0] + _radius * (1. + percent / 100.);
  coords[1][0] = target[1] - _radius * (1. + percent / 100.);
  coords[1][1] = target[1] + _radius * (1. + percent / 100.);
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

/****************************************************************************/
/*!
 **  Initialize the neighborhood search
 **
 ** \param[in]  dbin          input Db structure
 ** \param[in]  dbout         output Db structure (optional)
 **
 *****************************************************************************/
int NeighMoving::attach(const Db *dbin, const Db *dbout)
{
  if (ANeigh::attach(dbin, dbout)) return 1;
  int nech = _dbin->getSampleNumber();
  int ndim = _dbin->getNDim();
  int nsect = getNSect();
  _movingInd = VectorInt(nech);
  _movingDst = VectorDouble(nech);
  _movingIsect = VectorInt(nsect);
  _movingNsect = VectorInt(nsect);
  _movingX1 = VectorDouble(ndim);
  _movingX2 = VectorDouble(ndim);

  return 0;
}

/****************************************************************************/
/*!
 **  Returns the main Neighborhood Parameters for a given target as a vector:
 ** \li                    0 : Number of active samples
 ** \li                    1 : Minimum distance
 ** \li                    2 : Maximum distance
 ** \li                    3 : Number of non-empty sectors
 ** \li                    4 : Number of consecutive empty sectors
 **
 ** \param[in]  iech_out      Valid Rank of the sample in the output Db
 **
 *****************************************************************************/
VectorDouble NeighMoving::summary(int iech_out)
{
  VectorDouble tab(5,0.);

  /* Number of selected samples */

  VectorInt nbgh_ranks = select(iech_out);
  int nsel = (int) nbgh_ranks.size();
  tab[0] = (double) nsel;

  /* Maximum distance */

  double dmax = TEST;
  for (int iech = 0; iech < nsel; iech++)
  {
    double dist = _movingDst[iech];
    if (FFFF(dmax) || dist > dmax) dmax = dist;
  }
  tab[1] = dmax;

  /* Minimum distance */

  double dmin = TEST;
  for (int iech = 0; iech < nsel; iech++)
  {
    double dist = _movingDst[iech];
    if (FFFF(dmin) || dist < dmin) dmin = dist;
  }
  tab[2] = dmin;

  /* Number of sectors containing neighborhood information */

  int number = 0;
  for (int isect = 0; isect < getNSect(); isect++)
  {
    if (_movingNsect[isect] > 0) number++;
  }
  tab[3] = (double) number;

  /* Number of consecutive empty sectors */

  number = 0;
  int n_empty = 0;
  for (int isect = 0; isect < getNSect(); isect++)
  {
    if (_movingNsect[isect] > 0)
      n_empty = 0;
    else
    {
      n_empty++;
      if (n_empty > number) number = n_empty;
    }
  }
  if (getNSect() > 0)
  {
    if(_movingNsect[0] > 0)
      n_empty = 0;
    else
    {
      n_empty++;
      if (n_empty > number) number = n_empty;
    }
  }
  tab[4] = (double) number;

  return tab;
}

bool NeighMoving::hasChanged(int iech_out) const
{
  if (_iechMemo < 0 || _isNbghMemoEmpty()) return true;

  return true;
}

/****************************************************************************/
/*!
 **  Select the neighborhood
 **
 ** \return  Vector of sample ranks in neighborhood (empty when error)
 **
 ** \param[in]  iech_out      Valid Rank of the sample in the output Db
 **
 *****************************************************************************/
VectorInt NeighMoving::getNeigh(int iech_out)
{
  int nech = _dbin->getSampleNumber();
  VectorInt ranks(nech, -1);

  // Select the neighborhood samples as the target sample has changed
  if (_moving(iech_out, ranks)) return VectorInt();

  // In case of debug option, dump out neighborhood characteristics
  if (OptDbg::query(EDbg::NBGH)) _display(ranks);

  // Compress the vector of returned sample ranks
  _neighCompress(ranks);

  return ranks;
}

/****************************************************************************/
/*!
 **  Moving neighborhood search
 **
 ** \return  Error return code
 ** \return  0 : No error
 ** \return  1 : The number of data is smaller than the minimum number of
 ** \return      data in Moving Neighborhood
 **
 ** \param[in]  iech_out  Rank of the target in the output Db structure
 ** \param[in]  eps       Tolerance
 **
 ** \param[out]  ranks    Vector of samples elected in the Neighborhood
 **
 *****************************************************************************/
int NeighMoving::_moving(int iech_out, VectorInt& ranks, double eps)
{
  int nech = _dbin->getSampleNumber();
  int isect = 0;
  if (nech < getNMini()) return 1;

  /* Loop on the data points */

  double distmax = 0.;
  int nsel = 0;

  for (int iech = 0; iech < nech; iech++)
  {

    /* Discard the masked input sample */

    if (! _dbin->isActive(iech)) continue;

    /* Discard samples where all variables are undefined */

    if (_discardUndefined(iech)) continue;

    /* Discard the target sample for the cross-validation option */

    if (getFlagXvalid())
    {
      if (_xvalid(iech, iech_out)) continue;
    }

    // In presence of Faults, check that sample 'iech' is still eligible

    if (hasFaults())
    {
      if (_hiddenByFault(iech, iech_out)) continue;
    }

    // Force selection if the sample belongs to the block

    double dist;
    if (getForceWithinBlock() && _dbout->isGrid())
    {
      if (! _belongsToCell(iech, iech_out)) continue;
      dist = _movingDist(iech, iech_out);
    }
    else
    {
      /* Calculate the distance between data and target */

      dist = _movingDist(iech, iech_out);
      if (! FFFF(getRadius()) && dist > getRadius()) continue;
      if (dist > distmax) distmax = dist;

      /* Calculate the angular sector to which the sample belongs */

      if (getFlagSector())
        isect = _movingSectorDefine(_movingX1[0], _movingX1[1]);
    }

    /* The sample may be selected */

    _movingInd[nsel] = iech;
    _movingDst[nsel] = dist;
    ranks[iech] = isect;
    nsel++;
  }
  if (nsel < getNMini()) return 1;

  /* Slightly modify the distances in order to ensure the sorting results */
  /* In the case of equal distances                                       */

  for (int isel = 0; isel < nsel; isel++)
    _movingDst[isel] += distmax * isel * eps;

  /* Sort the selected samples according to the distance */

  ut_sort_double(0, nsel, _movingInd.data(), _movingDst.data());

  /* For each angular sector, select the first sample up to the maximum */

  if (! getForceWithinBlock() && getFlagSector() && getNSMax() > 0)
  {
    _movingSectorNsmax(nsel, ranks);
    if (nsel < getNMini()) return 1;
  }

    /* Select the first data samples (skipped if forcing all samples in block) */

  if (! getForceWithinBlock())
  {
    _movingSelect(nsel, ranks);
  }

  return 0;
}

bool NeighMoving::_hiddenByFault(int iech, int iech_out) const
{
  double xt1 = _dbout->getCoordinate(iech_out, 0);
  double yt1 = _dbout->getCoordinate(iech_out, 1);
  double xt2 = _dbin->getCoordinate(iech, 0);
  double yt2 = _dbin->getCoordinate(iech, 1);

  // Loop on the Fault polylines

  if (_faults == nullptr) return false;

  return _faults->isSplitByFault(xt1, yt1, xt2, yt2);
}

bool NeighMoving::_belongsToCell(int iech, int iech_out)
{
  if (_dbgrid == nullptr) return false;

  // Get the coordinates of the sample
  VectorDouble coor = _dbin->getSampleCoordinates(iech);

  // Identify the dimensions of the cell
  VectorDouble dxsPerCell = _dbgrid->getBlockExtensions(iech_out);

  // Check if the sample belongs to the cell
  return _dbgrid->getGrid().sampleBelongsToCell(coor, iech_out, dxsPerCell);
}

/****************************************************************************/
/*!
 **  Calculates the distance between an input sample and the target
 **
 ** \return  Distance
 **
 ** \param[in]  iech_in  Rank of the sample in the input Db structure
 ** \param[in]  iech_out Rank of the sample in the output Db structure
 **
 *****************************************************************************/
double NeighMoving::_movingDist(int iech_in, int iech_out)
{
  int ndim = _dbin->getNDim();

  /* Calculate the distance to the target */

  for (int idim = 0; idim < ndim; idim++)
    _movingX1[idim] = _dbout->getCoordinate(iech_out, idim)
        - _dbin->getCoordinate(iech_in, idim);

  /* Anisotropic neighborhood */

  if (getFlagAniso())
  {

    /* Rotated anisotropy ellipsoid */

    if (getFlagRotation())
    {
      matrix_product_safe(1, ndim, ndim, _movingX1.data(),
                          getAnisoRotMats().data(), _movingX2.data());
      _movingX1 = _movingX2;
    }
    for (int idim = 0; idim < ndim; idim++)
      _movingX1[idim] /= getAnisoCoeff(idim);
  }

  /* Calculate the distance */

  double dist;
  matrix_product_safe(1, ndim, 1, _movingX1.data(), _movingX1.data(), &dist);
  dist = sqrt(dist);
  return dist;
}

/****************************************************************************/
/*!
 **  Returns the rank of the sector (using the first two coordinates)
 **
 ** \return  Sector index
 **
 ** \param[in]  dx    increment along X
 ** \param[in]  dy    increment along Y
 **
 *****************************************************************************/
int NeighMoving::_movingSectorDefine(double dx, double dy)
{
  double angle;

  int isect = 0;
  if (getNSect() > 1)
  {
    if (dx == 0.)
    {
      if (dy >= 0.)
        angle = GV_PI / 2.;
      else
        angle = 1.5 * GV_PI;
    }
    else if (dx > 0.)
    {
      if (dy >= 0.)
        angle = atan(dy / dx);
      else
        angle = 2. * GV_PI - atan(-dy / dx);
    }
    else
    {
      if (dy > 0.)
        angle = GV_PI / 2. + atan(-dx / dy);
      else
        angle = GV_PI + atan(dy / dx);
    }
    isect = (int) (getNSect() * angle / (2. * GV_PI));
  }
  return (isect);
}

/****************************************************************************/
/*!
 **  For each angular sector, select the first samples until
 **  the maximum number of samples is reached
 **
 ** \param[in]  nsel   Number of selected samples
 **
 ** \param[out]  ranks Array of active data point ranks
 **
 *****************************************************************************/
void NeighMoving::_movingSectorNsmax(int nsel, VectorInt& ranks)
{
  int n_end = 0;
  for (int isect = 0; isect < getNSect(); isect++)
  {
    int n_ang = 0;
    for (int i = 0; i < nsel; i++)
    {
      int j = _movingInd[i];
      if (ranks[j] != isect) continue;
      if (n_ang < getNSMax())
        n_ang++;
      else
        ranks[j] = -1;
    }
    n_end += n_ang;
  }
  return;
}

/****************************************************************************/
/*!
 **  Select the closest data per sector, until the maximum number
 **  neighbors is reached
 **
 ** \param[in]  nsel  Number of ellegible data points
 ** \param[in]  ranks Rank of the ellegible samples
 **
 ** \remark  The samples beyond the maximum number of neighbors have their
 ** \remark  rank turned into -1
 **
 *****************************************************************************/
void NeighMoving::_movingSelect(int nsel, VectorInt& ranks)
{
  int number;

  if (getNMaxi() <= 0) return;
  for (int isect = 0; isect < getNSect(); isect++)
    _movingNsect[isect] = _movingIsect[isect] = 0;

  /* Count the number of samples per sector */

  number = 0;
  for (int i = 0; i < nsel; i++)
  {
    int j = _movingInd[i];
    int isect = ranks[j];
    if (isect < 0) continue;
    _movingNsect[isect]++;
    number++;
  }
  if (number < getNMaxi()) return;

  /* Find the rank of the admissible data per sector */

  number = 0;
  while (number < getNMaxi())
  {
    for (int isect = 0; isect < getNSect(); isect++)
    {
      if (_movingIsect[isect] >= _movingNsect[isect]) continue;
      _movingIsect[isect]++;
      number++;
      if (number >= getNMaxi()) break;
    }
  }

  /* Discard the data beyond the admissible rank per sector */

  for (int isect = 0; isect < getNSect(); isect++)
  {
    if (_movingIsect[isect] >= _movingNsect[isect]) continue;
    number = 0;
    for (int i = 0; i < nsel; i++)
    {
      int j = _movingInd[i];
      int jsect = ranks[j];
      if (jsect < 0) continue;
      if (isect != jsect) continue;
      number++;
      if (number > _movingIsect[isect]) ranks[j] = -1;
    }
  }
}

