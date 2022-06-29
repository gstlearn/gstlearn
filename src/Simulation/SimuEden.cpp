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
#include "geoslib_d.h"
#include "geoslib_old_f.h"

#include "Simulation/ASimulation.hpp"
#include "Simulation/SimuEden.hpp"
#include "Skin/Skin.hpp"
#include "Skin/ISkinFunctions.hpp"
#include "Basic/Law.hpp"
#include "Matrix/MatrixRectangular.hpp"

#define DIR_UP         4
#define DIR_DOWN       5

#define CORK_FACIES   -1
#define SHALE          0

#define CORK_FLUID    -2
#define NO_FLUID      -1
#define UNDEF_FLUID    0

static int invdir[6] = { 1, 0, 3, 2, 5, 4 };

SimuEden::SimuEden(int nbsimu, int seed)
    : ASimulation(nbsimu, seed),
      AStringable()
{
}

SimuEden::SimuEden(const SimuEden &r)
    : ASimulation(r),
      AStringable(r)
{
}

SimuEden& SimuEden::operator=(const SimuEden &r)
{
  if (this != &r)
  {
    ASimulation::operator =(r);
    AStringable::operator =(r);
  }
  return *this;
}

SimuEden::~SimuEden()
{
}


String SimuEden::toString(const AStringFormat* /*strfmt*/) const
{
  std::stringstream sstr;

  return sstr.str();
}

/*****************************************************************************/
/*!
 **  Multivariate multiphase propagation into a set of components
 **  constrained by initial conditions and fluid densities
 **
 ** \return  Error return code : 1 no fluid to propagate
 **
 ** \param[in]  dbgrid        Db grid structure
 ** \param[in]  ind_facies    Rank of the variable containing the Facies
 ** \param[in]  ind_fluid     Rank of the variable containing the Fluid
 ** \param[in]  ind_perm      Rank of the variable containing the Permeability
 ** \param[in]  ind_poro      Rank of the variable containing the Porosity
 ** \param[in]  nfacies       number of facies (facies 0 excluded)
 ** \param[in]  nfluids       number of fluids
 ** \param[in]  niter         Number of iterations
 ** \param[in]  iptr_fluid      Rank for storing the Fluid
 ** \param[in]  iptr_date       Rank for storing the Date
 ** \param[in]  iptr_stat_fluid Optional rank for storing (nfluids)
 ** \param[in]  iptr_stat_cork  Optional rank for storing Cork (1)
 ** \param[in]  speeds        array containing the travel speeds
 ** \param[in]  verbose       1 for a verbose option
 ** \param[in]  show_fluid    1 for modifying the value of the cells to show
 ** \li                       the initial valid fluid information
 ** \li                       the cork (different from shale)
 ** \param[in]  number_max    Maximum count of cells invaded (or TEST)
 ** \param[in]  volume_max    Maximum volume invaded (or TEST)
 **
 ** \remark  Directions are ordered as follows :
 ** \remark  0: +X; 1: -X; 2: +Y; 3: -Y; 4: +Z(up); 5: -Z(down)
 ** \remark  The coding of the matrix is:
 ** \remark              facies + nfacies * fluid
 ** \remark  Facies: 0 (Shale), 1 to nfacies, -1 (Cork)
 ** \remark  Fluids: 0 (undefined), 1 to nfluids, -1 (No Fluid)
 ** \remark  Fluids should be ordered by increasing weight
 ** \remark  A Permeability variable is a value (>=1) which divides
 ** \remark  the velocities. This variable is optional.
 ** \remark  A Porosity variable is a value (in [0,1]) which multiplies
 ** \remark  the volumes. This variable is optional.
 ** \remark  Volume_max represents the volumic part of the invaded area:
 ** \remark  it is always <= number of cells invaded.
 **
 *****************************************************************************/
int SimuEden::simulate(DbGrid *dbgrid,
                       int ind_facies,
                       int ind_fluid,
                       int ind_perm,
                       int ind_poro,
                       int nfacies,
                       int nfluids,
                       int niter,
                       int iptr_fluid,
                       int iptr_date,
                       int iptr_stat_fluid,
                       int iptr_stat_cork,
                       const VectorInt& speeds,
                       bool verbose,
                       bool show_fluid,
                       double number_max,
                       double volume_max)
{
  /* Define global variables */

  _dbgrid  = dbgrid;
  _nxyz    = _dbgrid->getSampleNumber();
  _nfacies = nfacies;
  _nfluids = nfluids;
  _speeds  = speeds;

  _indFacies = ind_facies;
  _indFluid  = ind_fluid;
  _indPerm   = ind_perm;
  _indPoro   = ind_poro;

  _iptrFluid = iptr_fluid;
  _iptrDate = iptr_date;
  _iptrStatFluid = iptr_stat_fluid;
  _iptrStatCork  = iptr_stat_cork;

  law_set_random_seed(getSeed());

  Skin* skin = new Skin(this, _dbgrid);

  /* Preliminary checks */

  if (! _fluid_check()) return 1;

  /* Printout of the fluid propagation parameters */

  _printParams(verbose);

  /* Core allocation */

  _statsDefine();

  /* Loop on the iterations */

  for (int iter = 0; iter < niter; iter++)
  {
    int seed_memo = law_get_random_seed();
    _statsReset();

    /* Check the consistency */

    _checkInconsistency(verbose);

    /* Initialize the grid with the initial values */

    _statsInit();
    if (skin->init(verbose))
    {
      delete skin;
      return 0;
    }

    /* Modifying the peripheral cells using a random walk */

    int idate = 0;
    while (skin->remains(verbose))
    {

      /* Check that the maximum quantities have not been reached */

      if (_checkMax(number_max, volume_max)) break;
      idate++;

      /* Find the next cell to be processed */

      int rank;
      int ipos;
      skin->getNext(&rank, &ipos);

      /* Find the new value of the target cell according to its neighborhood */

      int ref_fluid;
      if (_fluidModify(skin, ipos, &ref_fluid))
      {
        _ncork++;
        _setFACIES_CORK(ipos);
        _setFLUID(ipos, CORK_FLUID);
        _setDATE(ipos, ITEST);
      }
      else
      {
        _addStatNumber(_getFACIES(ipos)-1,ref_fluid-1,1);
        _addStatVolume(_getFACIES(ipos)-1,ref_fluid-1,_getPORO(ipos));
        _setFLUID(ipos, ref_fluid);
        _setDATE(ipos, idate);
      }

      /* Deduce the initial influence of the central cell */

      if (skin->unstack(rank, ipos))
      {
        delete skin;
        return 1;
      }
    }

    /* Final printout */

    if (verbose)
    {
      mestitle(1, "Final status (iteration %d)", iter + 1);
      message("- Seed Value                     = %d\n", seed_memo);
      _statsPrint("Cells already filled");
      _statsEmpty("Cells not reached");
    }

    /* Calculate statistics on fluids and corks */

    if (niter > 1) _calculateCumul();

    /* Update the data (optional) */

    _updateResults(iter < niter - 1, show_fluid);
  }

  /* Normalize the statistics */

  if (niter > 1) _normalizeCumul(niter);

  /* Set the error return flag */

  if (verbose) skin->skinPrint();

  return 0;
}

/****************************************************************************/
/*!
 **  Check the validity of the array Speed
 **
 ** \return  Error return code
 **
 *****************************************************************************/
bool SimuEden::_fluid_check(void)
{

  /* Check that there is no zero value */

  for (int ifluid = 0; ifluid < _nfluids; ifluid++)
    for (int ifacies = 0; ifacies < _nfacies; ifacies++)
      for (int idir = 0; idir < 6; idir++)
      {
        if (_getWT(ifacies + 1, ifluid + 1, 1, idir) <= 0)
        {
          messerr(
              "The Propagation Directional Speed is zero for: Fluid=%d - Facies=%d - Direction=%d",
              ifluid + 1, ifacies + 1, idir + 1);
          messerr("This may cause artifacts. Change it to a low value instead");
          return false;
        }
      }

  /* Check that at least one speed is defined (for each facies/fluid pair) */

  for (int ifluid = 0; ifluid < _nfluids; ifluid++)
    for (int ifacies = 0; ifacies < _nfacies; ifacies++)
    {
      double total = 0;
      for (int idir = 0; idir < 6; idir++)
        total += _getWT(ifacies + 1, ifluid + 1, 1, idir);
      if (total <= 0.)
      {
        messerr("For Facies (%d) and Fluid (%d), no positive speed is defined",
                ifacies + 1, ifluid + 1);
        return false;
      }
    }

  /* Check for wrong order relationship for velocities */

  for (int ifacies = 0; ifacies < _nfacies; ifacies++)
    for (int ifluid = 1; ifluid < _nfluids; ifluid++)
    {
      int jfluid = ifluid - 1;

      /* Z+ Speed */

      if (_getWT(ifacies + 1, jfluid + 1, 1, DIR_UP) <
          _getWT(ifacies + 1, ifluid + 1, 1, DIR_UP))
      {
        messerr("Error for the Z+ Propagation Directional Speed for Facies=%d:",
                ifacies + 1);
        messerr(
            "Speed for Fluid=%d [%d] must not be smaller than Speed for Fluid=%d [%d]",
            jfluid + 1, _getWT(ifacies + 1, jfluid + 1, 1, DIR_UP),
            ifluid + 1, _getWT(ifacies + 1, ifluid + 1, 1, DIR_UP));
        return false;
      }

      /* Z- Speed */

      if (_getWT(ifacies + 1, jfluid + 1, 1, DIR_DOWN) >
          _getWT(ifacies + 1, ifluid + 1, 1, DIR_DOWN))
      {
        messerr("Error for the Z- Propagation Directional Speed for Facies=%d:",
                ifacies + 1);
        messerr(
            "Speed for Fluid=%d [%d] must not be larger than Speed for Fluid=%d  [%d]",
            jfluid + 1, _getWT(ifacies + 1, jfluid + 1, 1, DIR_DOWN),
            ifluid + 1, _getWT(ifacies + 1, ifluid + 1, 1, DIR_DOWN));
        return false;
      }
    }
  return true;
}

/****************************************************************************/
/*!
 **  Transition speed for a Facies/Fluid pair in a given direction
 **
 ** \return  Transition speed
 **
 ** \param[in]  ifacies Facies value
 ** \param[in]  ifluid  Fluid value
 ** \param[in]  perm    Permeability value
 ** \param[in]  idir    Direction value
 **
 *****************************************************************************/
int SimuEden::_getWT(int ifacies, int ifluid, int perm, int idir)
{
  int ind, value;

  if (_speeds.empty())
    value = 1;
  else
  {
    ind = (idir) + 6 * ((ifacies - 1) * _nfluids + (ifluid - 1));
    value = perm * _speeds[ind];
  }
  return (value);
}

/****************************************************************************/
/*!
 **  Print the parameters of the fluid propagation simulation
 **
 *****************************************************************************/
void SimuEden::_printParams(bool verbose)

{
  if (! verbose) return;

  mestitle(0, "Fluid propagation parameters");
  message("Number of facies = %d\n", _nfacies);
  message("Number of fluids = %d\n", _nfluids);

  for (int ifacies = 0; ifacies < _nfacies; ifacies++)
    for (int ifluid = 0; ifluid < _nfluids; ifluid++)
    {
      message("Facies=%d - Fluid=%d -", ifacies + 1, ifluid + 1);
      for (int idir = 0; idir < 6; idir++)
        message(" Dir #%d=%d", idir+1,
                _getWT(ifacies + 1, ifluid + 1, 1, idir));
      message("\n");
    }
  return;
}

/*****************************************************************************/
/*!
 **  Allocate the Eden_Stats structure
 **
 ** \return  Pointer to the allocated Eden_Stats structure (or NULL)
 **
 *****************************************************************************/
void SimuEden::_statsDefine(void)

{
  _number.resize(_nfacies * _nfluids, 0);
  _volume.resize(_nfacies * _nfluids, 0.);
}

/****************************************************************************/
/*!
 **  Check if the cell is already filled with fluid
 **
 ** \return  1 if the cell (filled with facies) is already filled with Fluid
 **
 ** \param[in]  ipos   Absolute grid index of the input grid node
 **
 *****************************************************************************/
int SimuEden::isAlreadyFilled(int ipos) const
{
  return (_getFACIES(ipos) > 0 &&  _getPERM(ipos) > 0
          && _getFLUID(ipos) != UNDEF_FLUID);
}

/****************************************************************************/
/*!
 **  Check if the cell can be filled with fluid
 **
 ** \return  1 if the cell (filled with facies) can be filled with Fluid
 **
 ** \param[in]  ipos   Absolute grid index of the input grid node
 **
 *****************************************************************************/
int SimuEden::isToBeFilled(int ipos) const
{
  return (_getFACIES(ipos) > 0 && _getPERM(ipos) > 0
          && _getFLUID(ipos) == UNDEF_FLUID);
}

/****************************************************************************/
/*!
 **  Get the Facies value for a grid node
 **
 ** \return  Facies value
 **
 ** \param[in]  iech  Rank of the sample
 **
 *****************************************************************************/
int SimuEden::_getFACIES(int iech) const
{
  int ifacies = (int) _dbgrid->getArray(iech, _indFacies);
  if (ifacies < 0 || ifacies > _nfacies || IFFFF(ifacies)) ifacies = SHALE;
  return (ifacies);
}

/****************************************************************************/
/*!
 **  Get the Permeability value for a grid node
 **
 ** \return  Permeability value (>= 1)
 **
 ** \param[in]  iech  Rank of the sample
 **
 *****************************************************************************/
int SimuEden::_getPERM(int iech) const
{
  if (_indPerm <= 0) return (1);
  double perm = _dbgrid->getArray(iech, _indPerm);
  if (FFFF(perm) || perm < 0.) perm = 0.;
  return ((int) perm);
}

/****************************************************************************/
/*!
 **  Get the Date value for a grid node
 **
 ** \return  Date value
 **
 ** \param[in]  iech  Rank of the sample
 **
 *****************************************************************************/
double SimuEden::_getDATE(int iech)
{
  double date;

  if (_indDate <= 0) return (0);
  date = _dbgrid->getArray(iech, _indDate);
  if (FFFF(date)) return (0);
  date = MAX(1., date);
  return (date);
}

/****************************************************************************/
/*!
 **  Get the Fluid value for a grid node
 **
 ** \return  Fluid value
 **
 ** \param[in]  iech  Rank of the grid node
 **
 *****************************************************************************/
int SimuEden::_getFLUID(int iech) const
{
  int ifluid = (int) _dbgrid->getArray(iech, _indFluid);
  if (ifluid < 0 || ifluid > _nfluids || IFFFF(ifluid)) ifluid = UNDEF_FLUID;
  return (ifluid);
}

/****************************************************************************/
/*!
 **  Get the Fluid value for a grid node
 **  This routine is meant for questionning the old Fluid variable
 **
 ** \return  Fluid value
 **
 ** \param[in]  iech  Rank of the grid node
 **
 *****************************************************************************/
int SimuEden::_getFLUID_OLD(int iech) const
{
  double ifluid = _dbgrid->getArray(iech, _indFluid);
  if (ifluid < 0 || ifluid > _nfluids) ifluid = UNDEF_FLUID;
  return ((int) ifluid);
}

/****************************************************************************/
/*!
 **  Get the Porosity value for a grid node
 **
 ** \return  Porosity value (in [0,1])
 **
 ** \param[in]  iech  Rank of the sample
 **
 *****************************************************************************/
double SimuEden::_getPORO(int iech) const
{
  if (_indPoro <= 0) return (1);
  double poro = _dbgrid->getArray(iech, _indPoro);
  if (FFFF(poro)) return (0);
  poro = MIN(1., MAX(0., poro));
  return (poro);
}

/****************************************************************************/
/*!
 **  Returns the weight of a cell in a given direction
 **
 ** \return  The weight
 **
 ** \param[in]  ipos    Cell location
 ** \param[in]  idir    Direction value
 **
 *****************************************************************************/
double SimuEden::getWeight(int ipos, int idir) const
{
  double value;
  if (_speeds.empty())
    value = 1.;
  else
  {
    int ifacies = _getFACIES(ipos);
    int ifluid = _getFLUID(ipos);
    int perm = _getPERM(ipos);
    int ind = (idir) + 6 * ((ifacies - 1) * _nfluids + (ifluid - 1));
    value = perm * _speeds[ind];
  }
  return (value);
}

/*****************************************************************************/
/*!
 **  Reset the Eden_Stats structure
 **
 *****************************************************************************/
void SimuEden::_statsReset()

{
  _ncork = 0;
  for (int ifluid = 0; ifluid < _nfluids; ifluid++)
    for (int ifacies = 0; ifacies < _nfacies; ifacies++)
    {
      _setStatNumber(ifacies,ifluid,0);
      _setStatVolume(ifacies,ifluid,0.);
    }
}

void SimuEden::_setStatNumber(int ifacies, int ifluid, int value)
{
  _number[ifacies * _nfluids + ifluid] = value;
}
void SimuEden::_setStatVolume(int ifacies, int ifluid, double value)
{
  _volume[ifacies * _nfluids + ifluid] = value;
}
void SimuEden::_addStatNumber(int ifacies, int ifluid, int value)
{
  _number[ifacies * _nfluids + ifluid] += value;
}
void SimuEden::_addStatVolume(int ifacies, int ifluid, double value)
{
  _volume[ifacies * _nfluids + ifluid] = +value;
}
int SimuEden::_getStatNumber(int ifacies, int ifluid) const
{
  return _number[ifacies * _nfluids + ifluid];
}
double SimuEden::_getStatVolume(int ifacies, int ifluid) const
{
  return _volume[ifacies * _nfluids + ifluid];
}

/****************************************************************************/
/*!
 **  Initialize and check Facies & Fluid matrices
 **
 *****************************************************************************/
void SimuEden::_checkInconsistency(bool verbose)
{

  /* Loop on the cells of the matrix */

  int n_shale_fluid = 0;
  for (int iech = 0; iech < _nxyz; iech++)
  {
    int ifluid = _getFLUID_OLD(iech);
    int ifacies = _getFACIES(iech);
    double perm = _getPERM(iech);

    if (ifacies == SHALE || perm <= 0)
    {
      if (ifluid > 0)
      {
        if (verbose)
          messerr(
              "Cell %d: Inconsistent Fluid (%d) with Facies (%d) or Perm (%d) -> set to %d",
              iech + 1, ifluid, ifacies, perm, NO_FLUID);
        n_shale_fluid++;
      }
      _setFLUID(iech, NO_FLUID);
      _setFACIES(iech, SHALE);
      _setDATE(iech, ITEST);
    }
    else
    {
      _setFLUID(iech, ifluid);
      _setFACIES(iech, ifacies);
      _setDATE(iech, (ifluid > 0));
    }
  }

  /* Summary */

  if (n_shale_fluid > 0)
    message("Number of cells with inconsistent facies and fluid = %d\n",
            n_shale_fluid);
  return;
}

/****************************************************************************/
/*!
 **  Set the Fluid value for a grid node
 **
 ** \param[in]  iech   Rank of the grid node
 ** \param[in]  ifluid Fluid value
 **
 *****************************************************************************/
void SimuEden::_setFLUID(int iech, int ifluid)
{
  _dbgrid->setArray(iech, _indFluid, ifluid);
}

/****************************************************************************/
/*!
 **  Set the Facies value for a grid node
 **
 ** \param[in]  iech    Rank of the grid node
 ** \param[in]  ifacies Facies value
 **
 *****************************************************************************/
void SimuEden::_setFACIES(int iech, int ifacies)
{
  _dbgrid->setArray(iech, _indFacies, ifacies);
}

/****************************************************************************/
/*!
 **  Turn the Facies into Cork
 **
 ** \param[in]  iech   Rank of the grid node
 **
 *****************************************************************************/
void SimuEden::_setFACIES_CORK(int iech)
{
  int ifacies = (int) _dbgrid->getArray(iech, _indFacies);
  _dbgrid->setArray(iech, _indFacies, -ifacies);
  return;
}

/****************************************************************************/
/*!
 **  Set the Date for a grid node
 **
 ** \param[in]  iech   Rank of the grid node
 ** \param[in]  idate  Rank of the iteration
 **
 *****************************************************************************/
void SimuEden::_setDATE(int iech, int idate)
{
  double value = (IFFFF(idate)) ? TEST : idate;
  _dbgrid->setArray(iech, _iptrDate, value);
  return;
}

/*****************************************************************************/
/*!
 **  Initialize the Eden_Stats structure
 **
 *****************************************************************************/
void SimuEden::_statsInit()
{
  for (int lec = 0; lec < _nxyz; lec++)
  {
    if (! isAlreadyFilled(lec)) continue;

    /* The cell does not belong to the skin: it is already filled */

    _addStatNumber(_getFACIES(lec)-1,_getFLUID(lec)-1,1);
    _addStatVolume(_getFACIES(lec)-1,_getFLUID(lec)-1,_getPORO(lec));
  }
}

/****************************************************************************/
/*!
 **  Check that the Maximum quantities have been reached
 **
 *****************************************************************************/
int SimuEden::_checkMax(double number_max, double volume_max)
{
  if (FFFF(number_max) && FFFF(volume_max)) return (0);

  /* Print the statistics */

  double totnum = 0;
  double totvol = 0.;
  for (int ifluid = 0; ifluid < _nfluids; ifluid++)
    for (int ifacies = 0; ifacies < _nfacies; ifacies++)
    {
      int number = _getStatNumber(ifacies, ifluid);
      double volume = _getStatVolume(ifacies, ifluid);
      totnum += number;
      totvol += volume;
      if (!FFFF(number_max) && totnum >= number_max) return (1);
      if (!FFFF(volume_max) && totvol >= volume_max) return (1);
    }
  return (0);
}

/****************************************************************************/
/*!
 **  Calculate the new value for the fluid of the target cell
 **
 ** \return  0 for a regular cell and 1 for a newly generated cork
 **
 ** \param[in]  skin       Pointer to the skin
 ** \param[in]  ipos       Cell location
 **
 ** \param[out] ref_fluid_loc Current fluid value for the target cell
 **
 *****************************************************************************/
int SimuEden::_fluidModify(Skin *skin, int ipos, int *ref_fluid_loc)
{
  int ecr;
  int ref_fluid = UNDEF_FLUID;

  /* Loop on the directions */

  for (int dir = 0; dir < 6; dir++)
  {
    ecr = skin->gridShift(ipos, dir);
    if (! IFFFF(ecr))
    {
      if (isAlreadyFilled(ecr))
      {
        int fluid = _getFLUID(ecr);

        if (ref_fluid == UNDEF_FLUID)
        {
          if (dir == DIR_DOWN)
          {
            if (_getWT(_getFACIES(ecr), _getFLUID(ecr), _getPERM(ecr),
                       invdir[dir]) > 0)
              ref_fluid = fluid;
          }
          else if (dir == DIR_UP)
          {
            if (_getWT(_getFACIES(ecr), _getFLUID(ecr), _getPERM(ecr),
                       invdir[dir]) > 0)
              ref_fluid = fluid;
          }
          else
          {
            ref_fluid = fluid;
          }
        }
        else
        {
          if (dir == DIR_DOWN)
          {
            if (ref_fluid > fluid) return (1);
          }
          else if (dir == DIR_UP)
          {
            if (ref_fluid < fluid) return (1);
          }
          else
          {
            if (ref_fluid != fluid) return (1);
          }
        }
      }
    }
  }

  /* Returning argument */

  if (ref_fluid == UNDEF_FLUID) messageAbort("Undefined replacement Fluid");
  *ref_fluid_loc = ref_fluid;
  return (0);
}

/****************************************************************************/
/*!
 **  Print the statistics
 **
 ** \param[in]  title    Title
 **
 *****************************************************************************/
void SimuEden::_statsPrint(const char *title)

{
  /* Print the title */

  message("- %s\n", title);

  /* Print the statistics */

  double totnum = 0;
  double totvol = 0.;
  for (int ifluid = 0; ifluid < _nfluids; ifluid++)
    for (int ifacies = 0; ifacies < _nfacies; ifacies++)
    {
      int number    = _getStatNumber(ifacies, ifluid);
      double volume = _getStatVolume(ifacies, ifluid);
      totnum += number;
      totvol += volume;
      if (number > 0)
        message("  . Facies %d - Fluid %d  : Number = %d - Volume = %lf\n",
                ifacies + 1, ifluid + 1, number, volume);
    }
  if (totnum > 0)
  {
    message("           Total Number = %d\n", totnum);
    message("           Total Volume = %lf\n", totvol);
  }

  if (_ncork > 0) message("  . Cork                = %d\n", _ncork);

  return;
}

/****************************************************************************/
/*!
 **  Print the statistics for the cells not filled
 **
 ** \param[in]  title    Title
 **
 *****************************************************************************/
void SimuEden::_statsEmpty(const char *title)

{
  /* Print the statistics */

  double total = 0;
  int flag_title = 1;
  for (int ifacies = 0; ifacies < _nfacies; ifacies++)
  {
    int number = 0;
    for (int i = 0; i < _nxyz; i++)
    {
      if (! isToBeFilled(i)) continue;
      if (_getFACIES(i) == (ifacies + 1)) number++;
    }
    total += number;
    if (total > 0 && flag_title)
    {
      flag_title = 0;
      message("- %s\n", title);
    }
    if (number > 0)
      message("  . Facies %d not filled = %d\n", ifacies + 1, number);
  }
  if (total > 0) message("                  Total = %d\n", total);

  return;
}

/****************************************************************************/
/*!
 **  Calculate the statistics on Fluids and Corks
 **
 *****************************************************************************/
void SimuEden::_calculateCumul(void)

{
  /* Loop on the cells of the matrix */

  for (int iech = 0; iech < _nxyz; iech++)
  {

    /* Update the Fluid statistics */

    int ifluid = _getFLUID(iech);
    if (ifluid > 0) _dbgrid->updArray(iech, _iptrStatFluid + ifluid - 1, 0, 1);

    /* Update the Cork statistics */

    int ifacies = (int) _dbgrid->getArray(iech, _indFacies);
    if (ifacies < 0) _dbgrid->updArray(iech, _iptrStatCork, 0, 1);
  }
  return;
}

/****************************************************************************/
/*!
 **  Update the fluid at data location
 **
 ** \param[in]  reset_facies  option
 ** \li                       1 to reset the cork facies to initial value
 ** \li                       0 to set the facies to CORK_FACIES
 ** \param[in]  show_fluid    1 for modifying the value of the cells to show:
 ** \li                       the initial valid fluid information
 ** \li                       the cork (different from shale)
 **
 *****************************************************************************/
void SimuEden::_updateResults(int reset_facies, int show_fluid)

{
  /* Loop on the cells of the matrix */

  for (int iech = 0; iech < _nxyz; iech++)
  {
    int ifluid = _getFLUID_OLD(iech);
    int ifacies = (int) _dbgrid->getArray(iech, _indFacies);

    /* Update the Facies information */

    if (ifacies < 0)
    {
      if (reset_facies)
        _setFACIES(iech, -ifacies);
      else
        _setFACIES(iech, CORK_FACIES);
    }

    /* Update the Fluid information */

    if (show_fluid)
    {
      if (ifluid == NO_FLUID || ifluid == UNDEF_FLUID || ifluid == CORK_FLUID)
        continue;
      _setFLUID(iech, ifluid + _nfacies * _nfluids);
    }
    else
    {
      if (ifluid == CORK_FLUID) _setFLUID(iech, UNDEF_FLUID);
    }
  }
  return;
}

/****************************************************************************/
/*!
 **  Normalize the statistics on Fluids and Corks
 **
 ** \param[in]  niter  Number of iterations
 **
 *****************************************************************************/
void SimuEden::_normalizeCumul(int niter)

{
  /* Loop on the cells of the matrix */

  for (int iech = 0; iech < _nxyz; iech++)
  {

    /* Normalize the Fluid statistics */

    for (int ifluid = 0; ifluid < _nfluids; ifluid++)
      _dbgrid->updArray(iech, _iptrStatFluid + ifluid, 3, (double) niter);

    /* Update the Cork statistics */

    _dbgrid->updArray(iech, _iptrStatCork, 3, (double) niter);
  }

  return;
}

/****************************************************************************/
/*!
 **  Check if the sample belongs to the time slice
 **
 ** \return  Rank of the time slice (or -1)
 **
 ** \param[in]  date     Date attached to a sample
 ** \param[in]  ntime    Number of time intervals
 ** \param[in]  time0    Origin of the first time itnerval
 ** \param[in]  dtime    Time interval
 **
 *****************************************************************************/
int SimuEden::getTimeInterval(double date,
                              int ntime,
                              double time0,
                              double dtime)
{
  for (int itime = 0; itime < ntime; itime++)
  {
    double time_deb = time0 + dtime * itime;
    double time_fin = time0 + dtime * (itime + 1);
    if (date >= time_deb && date < time_fin) return (itime);
  }
  return (-1);
}

/*****************************************************************************/
/*!
 **  Extract time charts from the fluid propagation block
 **
 ** \return  Error return code
 **
 ** \param[in]  dbgrid        Pointer to the DbGrid structure
 ** \param[in]  ind_facies    Rank of the variable containing the Facies
 ** \param[in]  ind_fluid     Rank of the variable containing the Fluid
 ** \param[in]  ind_poro      Rank of the variable containing the Porosity
 ** \param[in]  ind_date      Rank of variable containing Date
 ** \param[in]  nfacies       number of facies (facies 0 excluded)
 ** \param[in]  nfluids       number of fluids
 ** \param[in]  facies0       Value of the target facies
 ** \param[in]  fluid0        Value of the target fluid
 ** \param[in]  ntime         Number of Time intervals
 ** \param[in]  time0         Starting time
 ** \param[in]  dtime         Time interval
 ** \param[in]  verbose       1 for a verbose option
 **
 *****************************************************************************/
MatrixRectangular SimuEden::fluidExtract(DbGrid* dbgrid,
                                         int ind_facies,
                                         int ind_fluid,
                                         int ind_poro,
                                         int ind_date,
                                         int nfacies,
                                         int nfluids,
                                         int facies0,
                                         int fluid0,
                                         int ntime,
                                         double time0,
                                         double dtime,
                                         bool verbose)
{
  MatrixRectangular tab;

  /* Preliminary checks */

  /* Define global variables */

  _dbgrid  = dbgrid;
  _nxyz    = _dbgrid->getSampleNumber();
  _nfacies = nfacies;
  _nfluids = nfluids;

  _indFacies = ind_facies;
  _indFluid  = ind_fluid;
  _indPoro   = ind_poro;
  _indDate   = ind_date;

  /* Initialize the array */

  tab = MatrixRectangular(ntime, 4);
  for (int itime = 0; itime < ntime; itime++)
  {
    tab.setValue(itime, 0, time0 + dtime * itime);
    tab.setValue(itime, 1, time0 + dtime * (itime + 1));
    tab.setValue(itime, 2, 0.);
    tab.setValue(itime, 3, 0.);
  }

  /* Loop on the blocks */

  double totnum = 0.;
  double totvol = 0.;
  double locnum = 0.;
  double locvol = 0.;
  double datmax = 0;
  for (int iech = 0; iech < _nxyz; iech++)
  {

    if (_getFACIES(iech) != facies0) continue;
    if (_getFLUID(iech) != fluid0) continue;
    double volume = _getPORO(iech);
    double date = _getDATE(iech);
    if (date > datmax) datmax = date;

    totnum += 1;
    totvol += volume;
    int itime = getTimeInterval(date, ntime, time0, dtime);
    if (itime < 0) continue;
    locnum += 1;
    locvol += volume;

    tab.setValue(itime, 2, tab.getValue(itime, 2) + 1);
    tab.setValue(itime, 3, tab.getValue(itime, 3) + volume);
  }

  /* Final printout */

  if (verbose)
  {
    mestitle(1, "Extraction for Fluid(%d) and Facies(%d)", facies0, fluid0);
    message("Time slices: From %lf to %lf by step of %lf\n",
            time0, time0 + dtime * ntime, dtime);
    message("Total Number of Cells               = %d\n", _nxyz);
    message("Maximum Date                        = %lf\n", datmax);
    message("Total Number of Invaded Cells       = %lf\n", totnum);
    message("Total Volume of Invaded Cells       = %lf\n", totvol);
    message("Total Number of Cells in Time Slice = %lf\n", locnum);
    message("Total Volume of Cells in Time Slice = %lf\n", locvol);
  }

  return tab;
}
