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
#include "Variogram/Vario.hpp"
#include "Variogram/VarioParam.hpp"
#include "Variogram/DirParam.hpp"
#include "Db/Db.hpp"
#include "Db/DbGrid.hpp"
#include "Model/Model.hpp"
#include "Basic/Utilities.hpp"
#include "Stats/Classical.hpp"
#include "Space/ASpaceObject.hpp"
#include "Space/ASpace.hpp"

VarioParam::VarioParam(double scale,
                       const VectorDouble& dates,
                       const Faults* faults)
  : AStringable()
  , _scale(scale)
  , _dates(dates)
  , _dirparams()
  , _faults(faults)
{
}

VarioParam::VarioParam(const VarioParam& VarioParam,
                       const VectorInt& dircols,
                       const Faults* faults)
  : AStringable()
  , _scale()
  , _dates()
  , _dirparams()
  , _faults(faults)
{
    _scale = VarioParam.getScale();
    _dates = VarioParam.getDates();

    for (int idir = 0; idir < (int) dircols.size(); idir++)
    {
      _dirparams.push_back(VarioParam.getDirParam(dircols[idir]));
    }
}

VarioParam::VarioParam(const VarioParam& r)
    : AStringable(r),
      _scale(r._scale),
      _dates(r._dates),
      _dirparams(r._dirparams),
      _faults(r._faults)
{
}

VarioParam& VarioParam::operator=(const VarioParam& r)
{
  if (this != &r)
  {
    AStringable::operator=(r);
    _scale = r._scale;
    _dates = r._dates;
    _dirparams  = r._dirparams;
    _faults = r._faults;
  }
  return *this;
}

VarioParam::~VarioParam()
{
}

bool VarioParam::isDefinedForGrid() const
{
  int ndir = getNDir();
  if (ndir <= 0) return false;
  return _dirparams[0].isDefinedForGrid();
}

bool VarioParam::_validDefinedFromGrid(const DirParam& dirparam) const
{
  int ndir = getNDir();
  bool definedFromGrid = dirparam.isDefinedForGrid();
  if (ndir > 0)
  {
    for (int idir = 0; idir < ndir; idir++)
    {
      if (_dirparams[idir].isDefinedForGrid() != definedFromGrid)
      {
        messerr("The current 'dirParam' cannot be added to 'varioParam'");
        if (_dirparams[idir].isDefinedForGrid())
          messerr("Element (%d) is defined using Grid definition",idir+1);
        else
          messerr("Element(%d) is defined NOT using Grid definition",idir+1);

        if (definedFromGrid)
          messerr("Current 'dirparam' is defined using Grid definition");
        else
          messerr("Current 'dirparam' is defined NOT using Grid definition");
        return false;
      }
    }
  }
  return true;
}

/**
 * Create one Calculation Direction corresponding to the Omni-direction calculation
 * For details, see DirParam::createOmniDirection documentation
 * @param nlag Number of lags
 * @param dlag Lag value
 * @param toldis Tolerance on distance
 * @param opt_code Option for usage of the code
 * @param idate Reference date
 * @param bench Bench value
 * @param cylrad Value of radius of the Cylinder search
 * @param tolcode Tolerance on the code
 * @param breaks Definition of the irregular intervals
 * @param scale Scaling factor
 * @param dates Range of dates
 * @param space Pointer to the space definition
 * @return
 */
VarioParam* VarioParam::createOmniDirection(int nlag,
                                            double dlag,
                                            double toldis,
                                            int opt_code,
                                            int idate,
                                            double bench,
                                            double cylrad,
                                            double tolcode,
                                            const VectorDouble& breaks,
                                            double scale,
                                            const VectorDouble& dates,
                                            const ASpaceSharedPtr& space)
{
  DirParam* dir = DirParam::createOmniDirection(nlag, dlag, toldis,
                                                opt_code, idate, bench, cylrad,
                                                tolcode, breaks, space);
  VarioParam* varioparam = new VarioParam(scale, dates);
  varioparam->addDir(*dir);
  delete dir;
  return varioparam;
}

VarioParam* VarioParam::createMultiple(int ndir,
                                       int nlag,
                                       double dlag,
                                       double toldis,
                                       double angref,
                                       double scale,
                                       const VectorDouble &dates,
                                       const ASpaceSharedPtr& space)
{
  std::vector<DirParam> dirs = DirParam::createMultiple(ndir, nlag, dlag,
                                                        toldis, angref, space);
  if (dirs.empty()) return nullptr;
  VarioParam* varioparam = new VarioParam(scale, dates);
  varioparam->addMultiDirs(dirs);
  return varioparam;
}

/**
 * Automatically create several calculation directions from Grid information:
 * For details, see DirParam::createMultipleFromGrid documentation
 * @param dbgrid a DbGrid structure
 * @param nlag Number of lags
 * @param scale Scaling factor
 * @param dates Range of dates
 * @param space Pointer to the Space definition
 * @param ndimax Maximum dimension (see note)
 * @return A pointer to the newly created VarioParam object
 *
 * @note This method creates as many calculation direction as space dimension
 * @note However, this number can be truncated to 'ndimax' (when defined)
 */
VarioParam* VarioParam::createMultipleFromGrid(const DbGrid* dbgrid,
                                               int nlag,
                                               double scale,
                                               const VectorDouble& dates,
                                               const ASpaceSharedPtr& space,
                                               int ndimax)
{
  VarioParam* varioparam = new VarioParam(scale, dates);
  int ndim               = dbgrid->getNDim();
  int ncalc = (ndimax <= 0) ? ndim : ndimax;
  VectorInt grincr(ndim, 0);
  for (int idim = 0; idim < ncalc; idim++)
  {
    VH::fill(grincr,  0.);
    grincr[idim] = 1;
    DirParam* dirparam = DirParam::createFromGrid(dbgrid, nlag, grincr, space);
    varioparam->addDir(*dirparam);
    delete dirparam;
  }
  return varioparam;
}

/**
 * Automatically create a set of calculation directions for a given Space Direction:
 * - one calculation direction per space direction
 * - the same parameters are used for each direction, such as:
 * @param nlag Number of lags
 * @param dlag Value of the lag
 * @param toldis Tolerance on distancecomputeFromDb
 * @param tolang Tolerance on angle
 * @param scale Scaling factor
 * @param dates Range of dates
 * @param space Pointer to the Space definition
 * @return
 */
VarioParam* VarioParam::createFromSpaceDimension(int nlag,
                                                 double dlag,
                                                 double toldis,
                                                 double tolang,
                                                 double scale,
                                                 const VectorDouble &dates,
                                                 const ASpaceSharedPtr& space)
{
  int ndim = getDefaultSpaceDimension();
  if (space != nullptr) ndim = space->getNDim();

  VarioParam* varioparam = new VarioParam(scale, dates);

  for (int idim = 0; idim < ndim; idim++)
  {
    DirParam dirparam(nlag, dlag, toldis, tolang, 0, 0, TEST, TEST, 0.,
                      VectorDouble(), VectorDouble(), TEST, space);

    VectorDouble codir(ndim,0.);
    codir[idim] = 1.;
    dirparam.setCodir(codir);
    varioparam->addDir(dirparam);
  }
  return varioparam;
}

VarioParam* VarioParam::createSeveral2D(const VectorDouble &angles,
                                        int nlag,
                                        double dlag,
                                        double toldis,
                                        double tolang,
                                        double scale,
                                        const VectorDouble& dates,
                                        const ASpaceSharedPtr& space)
{
  std::vector<DirParam> dirs = DirParam::createSeveral2D(angles, nlag, dlag,
                                                         toldis, tolang, space);
  if (dirs.empty()) return nullptr;
  VarioParam* varioparam = new VarioParam(scale, dates);
  varioparam->addMultiDirs(dirs);
  return varioparam;
}

void VarioParam::addDir(const DirParam& dirparam)
{
  if (! _validDefinedFromGrid(dirparam)) return;
  _dirparams.push_back(dirparam);
}

void VarioParam::addMultiDirs(const std::vector<DirParam>& dirparams)
{
  for (int i = 0; i < (int) dirparams.size(); i++)
  {
    if (! _validDefinedFromGrid(dirparams[i])) return;
    _dirparams.push_back(dirparams[i]);
  }
}

void VarioParam::delDir(int rank)
{
  if (rank < 0 || rank >= getNDir()) return;
  _dirparams.erase(_dirparams.begin() + rank);
}

void VarioParam::delAllDirs()
{
  _dirparams.clear();
}

String VarioParam::toString(const AStringFormat* strfmt) const
{
  std::stringstream sstr;
  if (getNDir() <= 0) return sstr.str();

  // Print the Main part

  sstr << toStringMain(strfmt);

  /* Loop on the directions */

  for (int idir=0; idir<getNDir(); idir++)
  {
    sstr << toTitle(1,"Direction #%d",idir+1);
    sstr << _dirparams[idir].toString(strfmt);
  }

  return sstr.str();
}

String VarioParam::toStringMain(const AStringFormat* /*strfmt*/) const
{
  std::stringstream sstr;
  int ndir = getNDir();

  /* General parameters */

  sstr << "Number of direction(s)      = " << ndir << std::endl;
  sstr << "Space dimension             = " << getNDim() << std::endl;

  if (hasDate())
  {
    sstr << "Number of Date Intervals    = " << getNDate() << std::endl;
    sstr << toMatrix("Matrix of Bounds for Data Intervals",VectorString(),VectorString(),
                     false,getNDate(),2,getDates());
  }

  if (hasFaults())
  {
    sstr << "Calculation takes Faults into account" << std::endl;
  }
  return sstr.str();
}

double VarioParam::getDate(int idate, int icas) const
{
  if (!_isDateValid(idate)) return 0.;
  return _dates[2 * idate + icas];
}

int VarioParam::getNLag(int idir) const
{
  if (! _isDirectionValid(idir)) return 0;
  return _dirparams[idir].getNLag();
}

VectorDouble VarioParam::getCodirs(int idir) const
{
  if (! _isDirectionValid(idir)) return VectorDouble();
  return _dirparams[idir].getCodirs();
}

bool VarioParam::_isDirectionValid(int idir) const
{
  return checkArg("Direction Index", idir, getNDir());
}

bool VarioParam::_isDateValid(int idate) const
{
  if (!hasDate()) return false;
  return checkArg("Date Index", idate, getNDate());
}

VectorDouble VarioParam::_getDirectionInterval(int idir) const
{
  VectorDouble bounds(2);
  if (idir < 0 || idir >= getNDim())
  {
    bounds[0] = 0;
    bounds[1] = getNDir();
  }
  else
  {
    bounds[0] = idir;
    bounds[1] = idir + 1;
  }
  return bounds;
}

void VarioParam::setDPas(int idir,const DbGrid* db)
{
  if (! _isDirectionValid(idir)) return;
  _dirparams[idir].setDPas(db);
}

void VarioParam::setGrincr(int idir, const VectorInt& grincr)
{
  if (! _isDirectionValid(idir)) return;
  _dirparams[idir].setGrincr(grincr);
}

int VarioParam::getNDim() const
{
  if (getNDir() <= 0) return 0;
  return _dirparams[0].getNDim();
}

/****************************************************************************/
/*!
 **  Check if dates are involved in the variogram calculation
 **
 ** \return  1 if dates are used; 0 otherwise
 **
 ** \param[in]  db1        First Db structure
 ** \param[in]  db2        Second Db structure
 **
 *****************************************************************************/
bool VarioParam::isDateUsed(const Db *db1, const Db *db2) const
{
  if (getDates().empty()) return false;
  if (!db1->hasLocVariable(ELoc::DATE)) return false;
  if (db2 != nullptr && !db2->hasLocVariable(ELoc::DATE)) return false;
  return true;
}

/****************************************************************************/
/*!
 **  Establish a new Db containing the pairs of the Variogram
 **
 ** \return  Pointer to the newly created Db
 **
 ** \param[in]  db           Db structure
 ** \param[in]  varioparam   VarioParam structure
 **
 *****************************************************************************/
Db* buildDbFromVarioParam(Db *db, const VarioParam& varioparam)
{
  if (db == nullptr) return nullptr;
  if (db->getNDim() != 2 && db->getNDim() != 3)
  {
    messerr("This function can only be calculated in dimension equal to 2 or 3");
    return nullptr;
  }
  if (db->getNSampleActive() <= 0)
  {
    messerr("This function requires your 'Db' to have some active samples");
    return nullptr;
  }
  if (db->getNLoc(ELoc::Z) != 1)
  {
    messerr("This function is restricted to the Monovariate case");
    return nullptr;
  }
  if (varioparam.getNDir() <= 0)
  {
    messerr("This function requires some direction to be defined in 'VarioParam'");
    return nullptr;
  }

  // Creating a local Vario structure (to constitute the BiTargetCheck list
  Vario vario = Vario(varioparam);
  vario.setDb(db);
  if (vario.prepare()) return nullptr;

  // Creating the output Db
  Db* newdb = Db::create();
  int ndim = db->getNDim();
  VectorVectorDouble ranks(2);
  VectorDouble lags;
  VectorDouble dirs;
  VectorDouble dists;
  VectorDouble vect(ndim);

  SpaceTarget T1(varioparam.getSpace());
  SpaceTarget T2(varioparam.getSpace());

  // Calculating the admissible pairs
  VectorInt rindex = db->getSortArray();

  // Local variables to speed up calculations
  bool hasSel = db->hasLocVariable(ELoc::SEL);
  bool hasWeight = db->hasLocVariable(ELoc::W);
  bool hasDate = varioparam.isDateUsed(db);
  double dist = 0.;

  for (int idir = 0; idir < varioparam.getNDir(); idir++)
  {
    const DirParam& dirparam = varioparam.getDirParam(idir);
    int nech = db->getNSample();
    double maxdist = dirparam.getMaximumDistance();

    /* Loop on the first point */

    for (int iiech = 0; iiech < nech - 1; iiech++)
    {
      int iech = rindex[iiech];
      if (hasSel && !db->isActive(iech)) continue;
      if (hasWeight && FFFF(db->getWeight(iech))) continue;
      db->getSampleAsSTInPlace(iech, T1);

      int ideb = (hasDate) ? 0 : iiech + 1;
      for (int jjech = ideb; jjech < nech; jjech++)
      {
        int jech = rindex[jjech];
        if (db->getDistance1D(iech, jech) > maxdist) break;
        if (hasSel && !db->isActive(jech)) continue;
        if (hasWeight && FFFF(db->getWeight(jech))) continue;
        db->getSampleAsSTInPlace(jech, T2);

        // Reject the point as soon as one BiTargetChecker is not correct
        if (! vario.keepPair(idir, T1, T2, &dist)) continue;

        /* Get the rank of the lag */

        int ilag = dirparam.getLagRank(dist);
        if (IFFFF(ilag)) continue;

        // The pair is kept

        ranks[0].push_back((double) iech);
        ranks[1].push_back((double) jech);
        ranks[0].push_back((double) jech);
        ranks[1].push_back((double) iech);
        dirs.push_back((double) idir);
        dirs.push_back((double) idir);
        lags.push_back((double) ilag);
        lags.push_back((double) ilag);
        dists.push_back(dist);
        dists.push_back(dist);
      }
    }
  }

  // Loading the coordinate vectors in the newly created Db

  newdb->addColumnsByVVD(ranks, "Sample", ELoc::UNKNOWN);
  newdb->addColumns(lags,  "Lag",ELoc::UNKNOWN);
  newdb->addColumns(dirs,  "Direction",ELoc::UNKNOWN);
  newdb->addColumns(dists, "Distance",ELoc::UNKNOWN);

  return newdb;
}

