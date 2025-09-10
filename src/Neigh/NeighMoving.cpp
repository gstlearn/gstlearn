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
#include "Neigh/NeighMoving.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/OptCustom.hpp"
#include "Basic/OptDbg.hpp"
#include "Basic/SerializeHDF5.hpp"
#include "Basic/Utilities.hpp"
#include "Basic/VectorHelper.hpp"
#include "Db/Db.hpp"
#include "Db/DbGrid.hpp"
#include "Geometry/BiTargetCheckDistance.hpp"
#include "Geometry/GeometryHelper.hpp"

#include <cmath>

namespace gstlrn
{

NeighMoving::NeighMoving(bool flag_xvalid,
                         Id nmaxi,
                         double radius,
                         Id nmini,
                         Id nsect,
                         Id nsmax,
                         const VectorDouble& coeffs,
                         const VectorDouble& angles,
                         bool useBallTree,
                         Id leaf_size,
                         const ASpaceSharedPtr& space)
  : ANeigh(space)
  , _nMini(nmini)
  , _nMaxi(nmaxi)
  , _nSect(nsect)
  , _nSMax(nsmax)
  , _distCont(TEST)
  , _biPtDist(nullptr)
  , _bipts()
  , _movingInd()
  , _movingIsect()
  , _movingNsect()
  , _movingDst()
  , _T1(space)
  , _T2(space)
{
  setFlagXvalid(flag_xvalid);

  setBallSearch(useBallTree, leaf_size);

  _biPtDist = BiTargetCheckDistance::create(radius, coeffs, angles);
}

NeighMoving::NeighMoving(const NeighMoving& r)
  : ANeigh(r)
  , _nMini(r._nMini)
  , _nMaxi(r._nMaxi)
  , _nSect(r._nSect)
  , _nSMax(r._nSMax)
  , _distCont(r._distCont)
  , _biPtDist(nullptr)
  , _bipts()
  , _movingInd(r._movingInd)
  , _movingIsect(r._movingIsect)
  , _movingNsect(r._movingNsect)
  , _movingDst(r._movingDst)
  , _dbgrid(r._dbgrid)
  , _T1(r._T1)
  , _T2(r._T2)
{
  for (Id ipt = 0, npt = static_cast<Id>(r._bipts.size()); ipt < npt; ipt++)
    _bipts.push_back(r._bipts[ipt]);

  // _biPtDist = r._biPtDist;
  _biPtDist = new BiTargetCheckDistance(*r._biPtDist);
}

NeighMoving& NeighMoving::operator=(const NeighMoving& r)
{
  if (this != &r)
  {
    ANeigh::operator=(r);
    _nMini       = r._nMini;
    _nMaxi       = r._nMaxi;
    _nSect       = r._nSect;
    _nSMax       = r._nSMax;
    _distCont    = r._distCont;
    _movingInd   = r._movingInd;
    _movingIsect = r._movingIsect;
    _movingNsect = r._movingNsect;
    _movingDst   = r._movingDst;
    _T1          = r._T1;
    _T2          = r._T2;

    for (Id ipt = 0, npt = static_cast<Id>(r._bipts.size()); ipt < npt; ipt++)
      //_bipts.push_back(dynamic_cast<ABiTargetCheck*>(r._bipts[ipt]->clone()));
      _bipts.push_back(r._bipts[ipt]);
    delete _biPtDist;
    //  _biPtDist = r._biPtDist;
    _biPtDist = new BiTargetCheckDistance(*r._biPtDist);
  }
  return *this;
}

NeighMoving::~NeighMoving()
{
  _bipts.clear();
  delete _biPtDist;
  _biPtDist = nullptr;
}

String NeighMoving::toString(const AStringFormat* strfmt) const
{
  std::stringstream sstr;

  sstr << toTitle(0, "Moving Neighborhood");

  if (_nMini > 0)
    sstr << "Minimum number of samples           = " << _nMini << std::endl;

  if (_nMaxi > 0)
    sstr << "Maximum number of samples           = " << _nMaxi << std::endl;

  if (_nSect > 1)
  {
    sstr << "Number of angular sectors           = " << _nSect << std::endl;
    if (_nSMax > 0)
      sstr << "Maximum number of points per sector = " << _nSMax << std::endl;
  }

  sstr << _biPtDist->toString(strfmt);

  auto number = _getNBiPts();
  for (Id ipt = 0; ipt < number; ipt++)
  {
    sstr << _bipts[ipt]->toString(strfmt);
  }

  if (getFlagContinuous())
  {
    sstr << "Norm. dist. for continuous NeighMoving.   = " << getDistCont()
         << std::endl;
  }
  return sstr.str();
}

bool NeighMoving::_deserializeAscii(std::istream& is, bool verbose)
{
  bool ret = true;
  ret      = ret && ANeigh::_deserializeAscii(is, verbose);
  if (!ret) return ret;

  auto ndim = getNDim();
  VectorDouble radius(ndim);
  VectorDouble nbgh_coeffs;
  VectorDouble nbgh_rotmat;

  Id flag_aniso    = 0;
  Id flag_rotation = 0;
  Id flag_sector   = 0;
  double dmax      = 0.;

  ret = true;
  ret = ret && _recordRead<Id>(is, "NeighMovingborhood sector search", flag_sector);
  ret = ret && _recordRead<Id>(is, "Minimum Number of samples", _nMini);
  ret = ret && _recordRead<Id>(is, "Maximum Number of samples", _nMaxi);
  ret = ret && _recordRead<Id>(is, "Optimum Number of samples per sector", _nSect);
  ret = ret && _recordRead<Id>(is, "Maximum Number of samples per sector", _nSMax);
  ret = ret && _recordRead<double>(is, "Maximum Isotropic Radius", dmax);
  ret = ret && _recordRead<Id>(is, "Flag for Anisotropy", flag_aniso);
  if (flag_aniso)
  {
    nbgh_coeffs.resize(ndim);
    for (size_t idim = 0; ret && idim < ndim; idim++)
      ret = ret && _recordRead<double>(is, "Anisotropy Coefficient", nbgh_coeffs[idim]);
    ret = ret && _recordRead<Id>(is, "Anisotropy Rotation Flag", flag_rotation);
    if (flag_rotation)
    {
      nbgh_rotmat.resize(ndim * ndim);
      Id lec = 0;
      for (size_t idim = 0; ret && idim < ndim; idim++)
        for (size_t jdim = 0; ret && jdim < ndim; jdim++, lec++)
          ret = ret && _recordRead<double>(is, "Anisotropy Rotation Matrix", nbgh_rotmat[lec]);
    }
  }
  if (!ret) return ret;

  setNSect((getFlagSector()) ? MAX(_nSect, 1) : 1);

  delete _biPtDist;
  _biPtDist = BiTargetCheckDistance::create(dmax, nbgh_coeffs);
  _biPtDist->setFlagRotation(flag_rotation);
  if (!nbgh_rotmat.empty())
    _biPtDist->setAnisoRotMat(nbgh_rotmat);

  return ret;
}

bool NeighMoving::_serializeAscii(std::ostream& os, bool verbose) const
{
  bool ret = true;
  ret      = ret && ANeigh::_serializeAscii(os, verbose);
  ret      = ret && _recordWrite<Id>(os, "Use angular sectors", getFlagSector());
  ret      = ret && _recordWrite<Id>(os, "", getNMini());
  ret      = ret && _recordWrite<Id>(os, "", getNMaxi());
  ret      = ret && _recordWrite<Id>(os, "", getNSect());
  ret      = ret && _recordWrite<Id>(os, "", getNSMax());
  ret      = ret && _commentWrite(os, "Parameters (nmini,nmaxi,nsect,nsmax)");

  // Store information from the Bipoint Checker based on distances
  ret = ret && _recordWrite<double>(os, "Maximum distance radius", _biPtDist->getRadius());
  ret = ret && _recordWrite<Id>(os, "Anisotropy Flag", _biPtDist->getFlagAniso());

  Id ndim = _biPtDist->getNDim();
  if (_biPtDist->getFlagAniso())
  {
    for (Id idim = 0; ret && idim < ndim; idim++)
      ret = ret && _recordWrite<double>(os, "", _biPtDist->getAnisoCoeff(idim));
    ret = ret && _commentWrite(os, "Anisotropy Coefficients");
    ret = ret && _recordWrite<Id>(os, "Anisotropy Rotation Flag",
                                  _biPtDist->getFlagRotation());
    if (_biPtDist->getFlagRotation())
    {
      Id ecr = 0;
      for (Id idim = 0; ret && idim < ndim; idim++)
        for (Id jdim = 0; ret && jdim < ndim; jdim++)
          ret = ret && _recordWrite<double>(os, "", _biPtDist->getAnisoRotMat(ecr++));
      ret = ret && _commentWrite(os, "Anisotropy Rotation Matrix");
    }
  }
  return ret;
}

NeighMoving* NeighMoving::create(bool flag_xvalid,
                                 Id nmaxi,
                                 double radius,
                                 Id nmini,
                                 Id nsect,
                                 Id nsmax,
                                 const VectorDouble& coeffs,
                                 const VectorDouble& angles,
                                 bool useBallTree,
                                 Id leaf_size,
                                 const ASpaceSharedPtr& space)
{
  return new NeighMoving(flag_xvalid, nmaxi, radius, nmini, nsect, nsmax,
                         coeffs, angles, useBallTree, leaf_size, space);
}

/**
 * Create a NeighMovingborhood by loading the contents of a Neutral File
 * @param NFFilename Name of the Neutral File
 * @param verbose    Verbose flag
 * @return
 */
NeighMoving* NeighMoving::createFromNF(const String& NFFilename, bool verbose)
{
  auto* neigh = new NeighMoving();
  if (neigh->_fileOpenAndDeserialize(NFFilename, verbose)) return neigh;
  delete neigh;
  return nullptr;
}

/**
 * Given a Db, returns the maximum number of samples per NeighMovingborhood
 * @param db Pointer to the target Db
 * @return
 */
Id NeighMoving::getNSampleMax(const Db* /*db*/) const
{
  return (getFlagSector()) ? _nSect * _nSMax : _nMaxi;
}

void NeighMoving::addBiTargetCheck(ABiTargetCheck* abpc)
{
  _bipts.push_back(abpc);
}

bool NeighMoving::getFlagSector() const
{
  return (getNDim() > 1 && _nSect > 1);
}

bool NeighMoving::_getAnisotropyElements(double* rx, double* ry, double* theta, double* cosp, double* sinp) const
{
  double radius = _getRadius();
  if (FFFF(radius)) return false;
  VectorDouble anisoRatio = getAnisoCoeffs();
  Id ndim                 = static_cast<Id>(anisoRatio.size());
  if (ndim != 2) return false;
  *rx = radius * anisoRatio[0];
  *ry = radius * anisoRatio[1];
  VectorDouble angrot(ndim);
  GH::rotationGetAnglesInPlace(getAnisoRotMats(), angrot);
  double angref = angrot[0] * GV_PI / 180.;
  *theta        = angrot[0];
  *cosp         = cos(angref);
  *sinp         = sin(angref);
  return true;
}

VectorVectorDouble NeighMoving::getEllipsoid(const VectorDouble& target, Id count) const
{
  double rx, ry, theta, cosp, sinp;
  if (!_getAnisotropyElements(&rx, &ry, &theta, &cosp, &sinp)) return VectorVectorDouble();

  return GH::getEllipse(target, rx, ry, theta, count);
}

VectorVectorDouble NeighMoving::getZoomLimits(const VectorDouble& target, double percent) const
{
  double rx, ry, theta, cosp, sinp;
  if (!_getAnisotropyElements(&rx, &ry, &theta, &cosp, &sinp)) return VectorVectorDouble();
  double radius = MAX(rx, ry);

  VectorVectorDouble coords(2);
  coords[0].resize(2, 0.);
  coords[1].resize(2, 0.);

  coords[0][0] = target[0] - radius * (1. + percent / 100.);
  coords[0][1] = target[0] + radius * (1. + percent / 100.);
  coords[1][0] = target[1] - radius * (1. + percent / 100.);
  coords[1][1] = target[1] + radius * (1. + percent / 100.);
  return coords;
}

/**
 * Generate the end-points of the sectors. By default, the extension is set to radius
 * @param target Coordinates of the Target
 * @return
 */
VectorVectorDouble NeighMoving::getSectors(const VectorDouble& target) const
{
  double rx, ry, theta, cosp, sinp;
  if (!_getAnisotropyElements(&rx, &ry, &theta, &cosp, &sinp)) return VectorVectorDouble();

  VectorVectorDouble coords(2);
  coords[0].resize(_nSect, 0.);
  coords[1].resize(_nSect, 0.);

  for (Id i = 0; i < _nSect; i++)
  {
    double angle = 2. * i * GV_PI / _nSect;
    double cosa  = cos(angle);
    double sina  = sin(angle);
    coords[0][i] = target[0] + rx * cosa * cosp - ry * sina * sinp;
    coords[1][i] = target[1] + rx * cosa * sinp + ry * sina * cosp;
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
Id NeighMoving::attach(const Db* dbin, const Db* dbout)
{
  if (ANeigh::attach(dbin, dbout)) return 1;

  _dbgrid = dynamic_cast<const DbGrid*>(_dbout);

  for (Id ipt = 0, nbt = _getNBiPts(); ipt < nbt; ipt++)
  {
    if (!_bipts[ipt]->isValid(dbin, dbout)) return 1;
  }

  Id nech      = _dbin->getNSample();
  auto nsect   = getNSect();
  _movingInd   = VectorInt(nech);
  _movingDst   = VectorDouble(nech);
  _movingIsect = VectorInt(nsect);
  _movingNsect = VectorInt(nsect);

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
VectorDouble NeighMoving::summary(Id iech_out)
{
  VectorDouble tab(5, 0.);

  /* Number of selected samples */

  VectorInt nbgh_ranks;
  select(iech_out, nbgh_ranks);
  Id nsel = static_cast<Id>(nbgh_ranks.size());
  tab[0]  = static_cast<double>(nsel);

  /* Maximum distance */

  double dmax = TEST;
  for (Id iech = 0; iech < nsel; iech++)
  {
    double dist = _movingDst[iech];
    if (FFFF(dmax) || dist > dmax) dmax = dist;
  }
  tab[1] = dmax;

  /* Minimum distance */

  double dmin = TEST;
  for (Id iech = 0; iech < nsel; iech++)
  {
    double dist = _movingDst[iech];
    if (FFFF(dmin) || dist < dmin) dmin = dist;
  }
  tab[2] = dmin;

  /* Number of sectors containing neighborhood information */

  Id number = 0;
  for (Id isect = 0; isect < getNSect(); isect++)
  {
    if (_movingNsect[isect] > 0) number++;
  }
  tab[3] = static_cast<double>(number);

  /* Number of consecutive empty sectors */

  number     = 0;
  Id n_empty = 0;
  for (Id isect = 0; isect < getNSect(); isect++)
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
    if (_movingNsect[0] > 0)
      n_empty = 0;
    else
    {
      n_empty++;
      if (n_empty > number) number = n_empty;
    }
  }
  tab[4] = static_cast<double>(number);

  return tab;
}

bool NeighMoving::hasChanged(Id iech_out) const
{
  DECLARE_UNUSED(iech_out);

  if (_iechMemo < 0 || _isNbghMemoEmpty()) return true;

  return true;
}

/**
 * Select the neighborhood
 * @param iech_out Valid Rank of the sample in the output Db
 * @param ranks Vector of sample ranks in neighborhood (empty when error)
 */
void NeighMoving::getNeigh(Id iech_out, VectorInt& ranks)
{
  // Select the neighborhood samples as the target sample has changed
  if (_moving(iech_out, ranks))
  {
    ranks.clear();
    return;
  }

  // In case of debug option, dump out neighborhood characteristics
  Id optim = OptCustom::query("Optim", 0);
  if (!optim || (optim && (iech_out + 1 == OptDbg::getReference())))
    if (OptDbg::query(EDbg::NBGH)) _display(ranks);

  // Compress the vector of returned sample ranks
  _neighCompress(ranks);
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
Id NeighMoving::_moving(Id iech_out, VectorInt& ranks, double eps)
{
  Id nech = _dbin->getNSample();
  ranks.resize(nech);
  ranks.fill(-1);
  Id isect = 0;
  if (nech < getNMini()) return 1;

  /* Loop on the data points */

  double distmax = 0.;
  Id nsel        = 0;

  // Load the target sample as a Space Point
  if (_dbgrid != nullptr)
    _dbgrid->getSampleAsSTInPlace(iech_out, _T1);
  else
    _dbout->getSampleAsSTInPlace(iech_out, _T1);

  // Select the elligible points when using Ball Tree search
  VectorInt elligibles;
  if (_useBallSearch)
  {
    elligibles = getBall().getIndices(_T1, _nMaxi);
    nech       = static_cast<Id>(elligibles.size());
  }

  for (Id jech = 0; jech < nech; jech++)
  {
    Id iech;
    if (_useBallSearch)
    {
      iech = elligibles[jech];
    }
    else
    {
      iech = jech;
      if (!_dbin->isActive(iech)) continue;
    }

    /* Discard sample if masked by a selection */

    if (!_dbin->isActive(iech)) continue;

    /* Discard samples where all variables are undefined */

    if (_discardUndefined(iech)) continue;

    /* Discard the target sample for the cross-validation option */

    if (getFlagXvalid())
    {
      if (_xvalid(iech, iech_out)) continue;
    }

    _dbin->getSampleAsSTInPlace(iech, _T2);

    // Reject the point due to BiTargetChecker
    // (other than the one based on distance which must come last)

    bool reject = false;
    for (Id ipt = 0, npt = _getNBiPts(); ipt < npt && !reject; ipt++)
    {
      if (!_bipts[ipt]->isOK(_T1, _T2)) reject = true;
    }
    if (reject) continue;

    // Calculate the distance between data and target
    // The rejection with respect to maximum distance is bypassed if
    // '_forceWithinCell'

    if (!_biPtDist->isOK(_T1, _T2)) continue;
    double dist = _biPtDist->getDistance();
    if (dist > distmax) distmax = dist;

    /* Calculate the angular sector to which the sample belongs */

    if (getFlagSector())
    {
      const VectorDouble& incr = _biPtDist->getIncr();
      isect                    = _movingSectorDefine(incr[0], incr[1]);
    }

    /* The sample may be selected */

    _movingInd[nsel] = iech;
    _movingDst[nsel] = dist;
    ranks[iech]      = isect;
    nsel++;
  }
  if (nsel < getNMini()) return 1;

  /* Slightly modify the distances in order to ensure the sorting results */
  /* In the case of equal distances                                       */

  for (Id isel = 0; isel < nsel; isel++)
    _movingDst[isel] += distmax * isel * eps;

  /* Sort the selected samples according to the distance */

  VH::arrangeInPlace(0, _movingInd, _movingDst, true, nsel);

  /* For each angular sector, select the first sample up to the maximum */

  if (getFlagSector() && getNSMax() > 0)
  {
    _movingSectorNsmax(nsel, ranks);
    if (nsel < getNMini()) return 1;
  }

  /* Select the first data samples (skipped if forcing all samples in block) */

  _movingSelect(nsel, ranks);

  return 0;
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
Id NeighMoving::_movingSectorDefine(double dx, double dy) const
{
  double angle;

  Id isect = 0;
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
    isect = static_cast<Id>(getNSect() * angle / (2. * GV_PI));
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
void NeighMoving::_movingSectorNsmax(Id nsel, VectorInt& ranks)
{
  for (Id isect = 0; isect < getNSect(); isect++)
  {
    Id n_ang = 0;
    for (Id i = 0; i < nsel; i++)
    {
      Id j = _movingInd[i];
      if (ranks[j] != isect) continue;
      if (n_ang < getNSMax())
        n_ang++;
      else
        ranks[j] = -1;
    }
  }
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
void NeighMoving::_movingSelect(Id nsel, VectorInt& ranks)
{
  Id number;

  if (getNMaxi() <= 0) return;
  for (Id isect = 0; isect < getNSect(); isect++)
    _movingNsect[isect] = _movingIsect[isect] = 0;

  /* Count the number of samples per sector */

  number = 0;
  for (Id i = 0; i < nsel; i++)
  {
    Id j     = _movingInd[i];
    Id isect = ranks[j];
    if (isect < 0) continue;
    _movingNsect[isect]++;
    number++;
  }
  if (number < getNMaxi()) return;

  /* Find the rank of the admissible data per sector */

  number = 0;
  while (number < getNMaxi())
  {
    for (Id isect = 0; isect < getNSect(); isect++)
    {
      if (_movingIsect[isect] >= _movingNsect[isect]) continue;
      _movingIsect[isect]++;
      number++;
      if (number >= getNMaxi()) break;
    }
  }

  /* Discard the data beyond the admissible rank per sector */

  for (Id isect = 0; isect < getNSect(); isect++)
  {
    if (_movingIsect[isect] >= _movingNsect[isect]) continue;
    number = 0;
    for (Id i = 0; i < nsel; i++)
    {
      Id j     = _movingInd[i];
      Id jsect = ranks[j];
      if (jsect < 0) continue;
      if (isect != jsect) continue;
      number++;
      if (number > _movingIsect[isect]) ranks[j] = -1;
    }
  }
}

#ifdef HDF5
bool NeighMoving::_deserializeH5(H5::Group& grp, [[maybe_unused]] bool verbose)
{
  auto neighG = SerializeHDF5::getGroup(grp, "NeighMoving");
  if (!neighG)
  {
    return false;
  }

  /* Read the grid characteristics */
  bool ret = true;

  ret = ret && ANeigh::_deserializeH5(*neighG, verbose);

  Id flag_aniso    = 0;
  Id flag_rotation = 0;
  Id flag_sector   = 0;
  double dmax      = 0.;
  VectorDouble nbgh_coeffs;
  VectorDouble nbgh_rotmat;

  ret = ret && SerializeHDF5::readValue(*neighG, "UseSectors", flag_sector);
  ret = ret && SerializeHDF5::readValue(*neighG, "NMini", _nMini);
  ret = ret && SerializeHDF5::readValue(*neighG, "NMaxi", _nMaxi);
  ret = ret && SerializeHDF5::readValue(*neighG, "NSect", _nSect);
  ret = ret && SerializeHDF5::readValue(*neighG, "NSMax", _nSMax);

  // Store information from the Bipoint Checker based on distances
  ret = ret && SerializeHDF5::readValue(*neighG, "Radius", dmax);
  ret = ret && SerializeHDF5::readValue(*neighG, "AnisoFlag", flag_aniso);

  if (flag_aniso)
  {
    ret = ret && SerializeHDF5::readVec(*neighG, "AnisoCoeff", nbgh_coeffs);
    ret = ret && SerializeHDF5::readValue(*neighG, "RotationFlag", flag_rotation);
    if (flag_rotation)
      ret = ret && SerializeHDF5::readVec(*neighG, "AnisoRotmat", nbgh_rotmat);
  }

  delete _biPtDist;
  _biPtDist = BiTargetCheckDistance::create(dmax, nbgh_coeffs);
  _biPtDist->setFlagRotation(flag_rotation);
  if (!nbgh_rotmat.empty())
    _biPtDist->setAnisoRotMat(nbgh_rotmat);

  return ret;
}

bool NeighMoving::_serializeH5(H5::Group& grp, [[maybe_unused]] bool verbose) const
{
  auto neighG = grp.createGroup("NeighMoving");

  bool ret = true;

  /* Writing the tail of the file */

  ret = ret && ANeigh::_serializeH5(neighG, verbose);

  ret = ret && SerializeHDF5::writeValue(neighG, "UseSectors", static_cast<Id>(getFlagSector()));
  ret = ret && SerializeHDF5::writeValue(neighG, "NMini", getNMini());
  ret = ret && SerializeHDF5::writeValue(neighG, "NMaxi", getNMaxi());
  ret = ret && SerializeHDF5::writeValue(neighG, "NSect", getNSect());
  ret = ret && SerializeHDF5::writeValue(neighG, "NSMax", getNSMax());

  // Store information from the Bipoint Checker based on distances
  ret = ret && SerializeHDF5::writeValue(neighG, "Radius", _biPtDist->getRadius());
  ret = ret && SerializeHDF5::writeValue(neighG, "AnisoFlag", _biPtDist->getFlagAniso());

  if (_biPtDist->getFlagAniso())
  {
    ret = ret && SerializeHDF5::writeVec(neighG, "AnisoCoeff", _biPtDist->getAnisoCoeffs());
    ret = ret && SerializeHDF5::writeValue(neighG, "RotationFlag", _biPtDist->getFlagRotation());
    if (_biPtDist->getFlagRotation())
      ret = ret && SerializeHDF5::writeVec(neighG, "AnisoRotmat", _biPtDist->getAnisoRotMats());
  }
  return ret;
}
#endif
} // namespace gstlrn
