/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://gstlearn.org                                              */
/* License: BSD 3 clauses                                                     */
/*                                                                            */
/******************************************************************************/
#include "geoslib_d.h"
#include "geoslib_old_f.h"

#include "Simulation/CalcSimuEden.hpp"
#include "Simulation/ACalcSimulation.hpp"
#include "Skin/Skin.hpp"
#include "Skin/ISkinFunctions.hpp"
#include "Basic/Law.hpp"
#include "Matrix/MatrixRectangular.hpp"
#include "Db/Db.hpp"

#define DIR_UP         4
#define DIR_DOWN       5

#define CORK_FACIES   -1
#define SHALE          0

#define CORK_FLUID    -2
#define NO_FLUID      -1
#define UNDEF_FLUID    0

static int invdir[6] = { 1, 0, 3, 2, 5, 4 };

CalcSimuEden::CalcSimuEden(int nfacies, int nfluids, int niter, int nbsimu, int seed, bool verbose)
    : ACalcSimulation(nbsimu, seed),
      AStringable(),
      _verbose(verbose),
      _showFluid(false),
      _iptrStatFluid(-1),
      _iptrStatCork(-1),
      _iptrFluid(-1),
      _iptrDate(-1),
      _niter(niter),
      _nfacies(nfacies),
      _nfluids(nfluids),
      _speeds(),
      _numberMax(TEST),
      _volumeMax(TEST),
      _indFacies(-1),
      _indFluid(-1),
      _indPerm(-1),
      _indPoro(-1),
      _indDate(-1),
      _nxyz(0),
      _ncork(0),
      _numbers(),
      _volumes()
{
}

CalcSimuEden::~CalcSimuEden()
{
}

String CalcSimuEden::toString(const AStringFormat* /*strfmt*/) const
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
 *****************************************************************************/
bool CalcSimuEden::_simulate()
{
  DbGrid* dbgrid = dynamic_cast<DbGrid*>(getDbout());
  _nxyz    = dbgrid->getSampleNumber();

  Skin* skin = new Skin(this, dbgrid);

  /* Preliminary checks */

  if (! _fluid_check()) return false;

  /* Printout of the fluid propagation parameters */

  _printParams(_verbose);

  /* Core allocation */

  _statsDefine();

  /* Loop on the iterations */

  for (int iter = 0; iter < _niter; iter++)
  {
    int seed_memo = law_get_random_seed();
    _statsReset();

    /* Check the consistency */

    _checkInconsistency(_verbose);

    /* Initialize the grid with the initial values */

    _statsInit();

    if (skin->init(_verbose))
    {
      delete skin;
      return true;
    }

    /* Modifying the peripheral cells using a random walk */

    int idate = 0;
    while (skin->remains(_verbose))
    {

      /* Check that the maximum quantities have not been reached */

      if (_checkMax(_numberMax, _volumeMax)) break;
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
        return false;
      }
    }

    /* Final printout */

    if (_verbose)
    {
      mestitle(1, "Final status (iteration %d)", iter + 1);
      message("- Seed Value                     = %d\n", seed_memo);
      _statsPrint("Cells already filled");
      _statsEmpty("Cells not reached");
    }

    /* Calculate statistics on fluids and corks */

    if (_niter > 1) _calculateCumul();

    /* Update the data (optional) */

    _updateResults(iter < _niter - 1, _showFluid);
  }

  /* Normalize the statistics */

  if (_niter > 1) _normalizeCumul(_niter);

  /* Print statistics */

  if (_verbose) skin->skinPrint();

  return true;
}

/****************************************************************************/
/*!
 **  Check the validity of the array Speed
 **
 ** \return  Error return code
 **
 *****************************************************************************/
bool CalcSimuEden::_fluid_check(void)
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
int CalcSimuEden::_getWT(int ifacies, int ifluid, int perm, int idir)
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
void CalcSimuEden::_printParams(bool verbose)

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
 *****************************************************************************/
void CalcSimuEden::_statsDefine(void)

{
  _numbers.resize(_nfacies * _nfluids, 0);
  _volumes.resize(_nfacies * _nfluids, 0.);
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
int CalcSimuEden::isAlreadyFilled(int ipos) const
{
  bool answer = _getFACIES(ipos) > 0 &&
                _getPERM(ipos)   > 0 &&
                _getFLUID(ipos) != UNDEF_FLUID;
  return answer;
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
int CalcSimuEden::isToBeFilled(int ipos) const
{
  bool answer = _getFACIES(ipos) > 0 &&
                _getPERM(ipos)   > 0 &&
                _getFLUID(ipos) == UNDEF_FLUID;
  return answer;
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
int CalcSimuEden::_getFACIES(int iech) const
{
  DbGrid* dbgrid = dynamic_cast<DbGrid*>(getDbout());
  int ifacies = (int) dbgrid->getArray(iech, _indFacies);
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
int CalcSimuEden::_getPERM(int iech) const
{
  if (_indPerm <= 0) return (1);
  DbGrid* dbgrid = dynamic_cast<DbGrid*>(getDbout());
  double perm = dbgrid->getArray(iech, _indPerm);
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
double CalcSimuEden::_getDATE(int iech)
{
  double date;

  if (_indDate <= 0) return (0);
  DbGrid* dbgrid = dynamic_cast<DbGrid*>(getDbout());
  date = dbgrid->getArray(iech, _indDate);
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
int CalcSimuEden::_getFLUID(int iech) const
{
  DbGrid* dbgrid = dynamic_cast<DbGrid*>(getDbout());
  int ifluid = (int) dbgrid->getArray(iech, _indFluid);
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
int CalcSimuEden::_getFLUID_OLD(int iech) const
{
  DbGrid* dbgrid = dynamic_cast<DbGrid*>(getDbout());
  double ifluid = dbgrid->getArray(iech, _indFluid);
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
double CalcSimuEden::_getPORO(int iech) const
{
  if (_indPoro <= 0) return (1);
  DbGrid* dbgrid = dynamic_cast<DbGrid*>(getDbout());
  double poro = dbgrid->getArray(iech, _indPoro);
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
double CalcSimuEden::getWeight(int ipos, int idir) const
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
void CalcSimuEden::_statsReset()

{
  _ncork = 0;
  for (int ifluid = 0; ifluid < _nfluids; ifluid++)
    for (int ifacies = 0; ifacies < _nfacies; ifacies++)
    {
      _setStatNumber(ifacies,ifluid,0);
      _setStatVolume(ifacies,ifluid,0.);
    }
}

void CalcSimuEden::_setStatNumber(int ifacies, int ifluid, int value)
{
  _numbers[ifacies * _nfluids + ifluid] = value;
}
void CalcSimuEden::_setStatVolume(int ifacies, int ifluid, double value)
{
  _volumes[ifacies * _nfluids + ifluid] = value;
}
void CalcSimuEden::_addStatNumber(int ifacies, int ifluid, int value)
{
  _numbers[ifacies * _nfluids + ifluid] += value;
}
void CalcSimuEden::_addStatVolume(int ifacies, int ifluid, double value)
{
  _volumes[ifacies * _nfluids + ifluid] += value;
}
int CalcSimuEden::_getStatNumber(int ifacies, int ifluid) const
{
  return _numbers[ifacies * _nfluids + ifluid];
}
double CalcSimuEden::_getStatVolume(int ifacies, int ifluid) const
{
  return _volumes[ifacies * _nfluids + ifluid];
}

/****************************************************************************/
/*!
 **  Initialize and check Facies & Fluid matrices
 **
 *****************************************************************************/
void CalcSimuEden::_checkInconsistency(bool verbose)
{

  /* Loop on the cells of the matrix */

  int n_shale_fluid = 0;
  for (int iech = 0; iech < _nxyz; iech++)
  {
    int ifluid  = _getFLUID_OLD(iech);
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
void CalcSimuEden::_setFLUID(int iech, int ifluid)
{
  DbGrid* dbgrid = dynamic_cast<DbGrid*>(getDbout());
  dbgrid->setArray(iech, _indFluid, ifluid);
}

/****************************************************************************/
/*!
 **  Set the Facies value for a grid node
 **
 ** \param[in]  iech    Rank of the grid node
 ** \param[in]  ifacies Facies value
 **
 *****************************************************************************/
void CalcSimuEden::_setFACIES(int iech, int ifacies)
{
  DbGrid* dbgrid = dynamic_cast<DbGrid*>(getDbout());
  dbgrid->setArray(iech, _indFacies, ifacies);
}

/****************************************************************************/
/*!
 **  Turn the Facies into Cork
 **
 ** \param[in]  iech   Rank of the grid node
 **
 *****************************************************************************/
void CalcSimuEden::_setFACIES_CORK(int iech)
{
  DbGrid* dbgrid = dynamic_cast<DbGrid*>(getDbout());
  int ifacies = (int) dbgrid->getArray(iech, _indFacies);
  dbgrid->setArray(iech, _indFacies, -ifacies);
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
void CalcSimuEden::_setDATE(int iech, int idate)
{
  double value = (IFFFF(idate)) ? TEST : idate;
  DbGrid* dbgrid = dynamic_cast<DbGrid*>(getDbout());
  dbgrid->setArray(iech, _iptrDate, value);
  return;
}

/*****************************************************************************/
/*!
 **  Initialize the Eden_Stats structure
 **
 *****************************************************************************/
void CalcSimuEden::_statsInit()
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
int CalcSimuEden::_checkMax(double number_max, double volume_max)
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
int CalcSimuEden::_fluidModify(Skin *skin, int ipos, int *ref_fluid_loc)
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
void CalcSimuEden::_statsPrint(const char *title)

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
    message("           Total Number = %lf\n", totnum);
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
void CalcSimuEden::_statsEmpty(const char *title)

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
void CalcSimuEden::_calculateCumul(void)

{
  DbGrid* dbgrid = dynamic_cast<DbGrid*>(getDbout());

  /* Loop on the cells of the matrix */

  for (int iech = 0; iech < _nxyz; iech++)
  {

    /* Update the Fluid statistics */

    int ifluid = _getFLUID(iech);
    if (ifluid > 0) dbgrid->updArray(iech, _iptrStatFluid + ifluid - 1, 0, 1);

    /* Update the Cork statistics */

    int ifacies = (int) dbgrid->getArray(iech, _indFacies);
    if (ifacies < 0) dbgrid->updArray(iech, _iptrStatCork, 0, 1);
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
void CalcSimuEden::_updateResults(int reset_facies, int show_fluid)

{
  DbGrid* dbgrid = dynamic_cast<DbGrid*>(getDbout());

  /* Loop on the cells of the matrix */

  for (int iech = 0; iech < _nxyz; iech++)
  {
    int ifluid = _getFLUID_OLD(iech);
    int ifacies = (int) dbgrid->getArray(iech, _indFacies);

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
void CalcSimuEden::_normalizeCumul(int niter)

{
  DbGrid* dbgrid = dynamic_cast<DbGrid*>(getDbout());

  /* Loop on the cells of the matrix */

  for (int iech = 0; iech < _nxyz; iech++)
  {

    /* Normalize the Fluid statistics */

    for (int ifluid = 0; ifluid < _nfluids; ifluid++)
      dbgrid->updArray(iech, _iptrStatFluid + ifluid, 3, (double) niter);

    /* Update the Cork statistics */

    dbgrid->updArray(iech, _iptrStatCork, 3, (double) niter);
  }

  return;
}

int CalcSimuEden::_countAlreadyFilled() const
{
  int count = 0;
  for (int lec = 0; lec < _nxyz; lec++)
    count += isAlreadyFilled(lec);
  return count;
}

int CalcSimuEden::_countIsToBeFilled() const
{
  int count = 0;
  for (int lec = 0; lec < _nxyz; lec++)
    count += isToBeFilled(lec);
  return count;
}

bool CalcSimuEden::_check()
{
  if (! ACalcSimulation::_check()) return false;

  if (! hasDbout()) return false;
  int ndim = _getNDim();
  if (ndim > 3)
  {
    messerr("The Turning Band Method is not a relevant simulation model");
    messerr("for this Space Dimension (%d)", ndim);
    return false;
  }
  if (! getDbout()->isGrid())
  {
    messerr("The argument 'dbout'  should be a grid");
    return false;
  }

  if (_indFacies < 0)
  {
    messerr("Variable 'Facies' must be provided");
    return false;
  }
  if (_indFluid < 0)
  {
    messerr("Variable 'Fluid' must be provided");
    return false;
  }

  return true;
}

bool CalcSimuEden::_preprocess()
{

  /* Add the attributes for storing the results */

  if (_niter > 1)
  {
    _iptrStatFluid = _addVariableDb(2, 1, ELoc::UNKNOWN, 0, _nfluids, 0.);
    if (_iptrStatFluid < 0) return false;
    _iptrStatCork = _addVariableDb(2, 1, ELoc::UNKNOWN, 0, 1, 0.);
    if (_iptrStatCork < 0) return false;
  }

  /* Add the attributes for storing the Fluid and Data informations */

  _iptrFluid = _addVariableDb(2, 1, ELoc::UNKNOWN, 0, 1, 0.);
  if (_iptrFluid < 0) return false;
  _iptrDate = _addVariableDb(2, 1, ELoc::UNKNOWN, 0, 1, TEST);
  if (_iptrDate < 0) return false;

  return true;
}

bool CalcSimuEden::_run()
{
  law_set_random_seed(getSeed());

  return (_simulate());
}

bool CalcSimuEden::_postprocess()
{
  /* Free the temporary variables */
  _cleanVariableDb(2);

  if (_iptrStatFluid >= 0)
    _renameVariable(2, 1, _iptrStatFluid, "Stat_Fluid", _niter);
  if (_iptrStatCork >= 0)
    _renameVariable(2, 1, _iptrStatCork, "Stat_Cork", _niter);
  if (_iptrFluid)
    _renameVariable(2, 1, _iptrFluid, "Fluid", 1);
  if (_iptrDate)
    _renameVariable(2, 1, _iptrDate, "Date", 1);
  return true;
}

void CalcSimuEden::_rollback()
{
  _cleanVariableDb(1);
}

/*****************************************************************************/
/*!
**  Multivariate multiphase propagation into a set of components
**  constrained by initial conditions and fluid densities
**
** \return  Error return code : 1 no fluid to propagate
**
** \param[in]  dbgrid        Db grid structure
** \param[in]  name_facies   Name of the variable containing the Facies
** \param[in]  name_fluid    Name of the variable containing the Fluid
** \param[in]  name_perm     Name of the variable containing the Permeability
** \param[in]  name_poro     Name of the variable containing the Porosity
** \param[in]  nfacies       number of facies (facies 0 excluded)
** \param[in]  nfluids       number of fluids
** \param[in]  niter         Number of iterations
** \param[in]  speeds        Array containing the travel speeds
** \param[in]  show_fluid    1 for modifying the value of the cells to show
** \li                       the initial valid fluid information
** \li                       the cork (different from shale)
** \param[in]  number_max    Maximum count of cells invaded (or TEST)
** \param[in]  volume_max    Maximum volume invaded (or TEST)
** \param[in]  seed          Seed for random number generator (or 0)
** \param[in]  verbose       1 for a verbose option
** \param[in]  namconv       Naming convention
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
** \remark  the velocities. This variable is optional.
**
** \remark  the volumes. This variable is optional.
** \remark  Volume_max represents the volumic part of the invaded area:
** \remark  it is always <= number of cells invaded.
**
*****************************************************************************/
int fluid_propagation(DbGrid *dbgrid,
                      const String& name_facies,
                      const String& name_fluid,
                      const String& name_perm,
                      const String& name_poro,
                      int     nfacies,
                      int     nfluids,
                      int     niter,
                      const VectorInt& speeds,
                      bool    show_fluid,
                      double  number_max,
                      double  volume_max,
                      int seed,
                      bool verbose,
                      const NamingConvention& namconv)
{
  CalcSimuEden seden(nfacies, nfluids, niter, 1, seed, verbose);

  seden.setDbout(dbgrid);
  seden.setNamingConvention(namconv);

  seden.setIndFacies(dbgrid->getUID(name_facies));
  seden.setIndFluid(dbgrid->getUID(name_fluid));
  if (! name_poro.empty()) seden.setIndPoro(dbgrid->getUID(name_poro));
  if (! name_perm.empty()) seden.setIndPerm(dbgrid->getUID(name_perm));

  seden.setSpeeds(speeds);
  seden.setShowFluid(show_fluid);
  seden.setNumberMax(number_max);
  seden.setVolumeMax(volume_max);

  // Run the calculator
  int error = (seden.run()) ? 0 : 1;
  return error;
}
