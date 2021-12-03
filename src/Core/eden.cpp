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
#include "Basic/Utilities.hpp"
#include "Basic/Law.hpp"
#include "Db/Db.hpp"

/*! \cond */
#define DIR_UP         4
#define DIR_DOWN       5

#define CORK_FACIES   -1
#define SHALE          0

#define CORK_FLUID    -2
#define NO_FLUID      -1
#define UNDEF_FLUID    0

#define STAT_NUM(ifacies,ifluid) (stats->number[(ifacies) * NFLUIDS + (ifluid)])
#define STAT_VOL(ifacies,ifluid) (stats->volume[(ifacies) * NFLUIDS + (ifluid)])
#define TAB(i,itime)             (tab[4 * (itime) + (i)])
/*! \endcond */

typedef struct
{
  int ncork;
  int *number;
  double *volume;
} Eden_Stats;

static int invdir[6] = { 1, 0, 3, 2, 5, 4 };
static const char *dirchar[] = { "X+", "X-", "Y+", "Y-", "Z+", "Z-" };
static int NFLUIDS, NXYZ, IND_FACIES, IND_FLUID, IND_PERM, IND_PORO, IND_DATE;
static int NFACIES, IPTR_FLUID, IPTR_DATE, IPTR_STAT_FLUID, IPTR_STAT_CORK;
static int *SPEEDS, VERBOSE;
static Db *DB;

/*****************************************************************************/
/*!
 **  Deallocate the Eden_Stats structure
 **
 ** \return  Pointer to the newly deallocated stats
 **
 ** \param[in]  stats    Pointer to the Eden_Stats to be deallocated
 **
 *****************************************************************************/
static Eden_Stats* st_stats_undefine(Eden_Stats *stats)

{

  /* Deallocation */

  if (stats != nullptr)
  {
    stats->number = (int*) mem_free((char* ) stats->number);
    stats->volume = (double*) mem_free((char* ) stats->volume);
  }
  stats = (Eden_Stats*) mem_free((char* ) stats);

  return (stats);
}

/*****************************************************************************/
/*!
 **  Reset the Eden_Stats structure
 **
 ** \param[in]  stats  Eden_Stats structure to be reset
 **
 *****************************************************************************/
static void st_stats_reset(Eden_Stats *stats)

{
  int ifacies, ifluid;

  stats->ncork = 0;
  for (ifluid = 0; ifluid < NFLUIDS; ifluid++)
    for (ifacies = 0; ifacies < NFACIES; ifacies++)
    {
      STAT_NUM(ifacies,ifluid) = 0;
      STAT_VOL(ifacies,ifluid) = 0.;
    }
}

/*****************************************************************************/
/*!
 **  Allocate the Eden_Stats structure
 **
 ** \return  Pointer to the allocated Eden_Stats structure (or NULL)
 **
 *****************************************************************************/
static Eden_Stats* st_stats_define(void)

{
  Eden_Stats *stats;
  int error;

  /* Initializations */

  error = 1;

  /* Allocation */

  stats = (Eden_Stats*) mem_alloc(sizeof(Eden_Stats), 0);
  if (stats == nullptr) goto label_end;
  stats->number = (int*) mem_alloc(sizeof(int) * NFACIES * NFLUIDS, 0);
  if (stats->number == nullptr) goto label_end;
  stats->volume = (double*) mem_alloc(sizeof(double) * NFACIES * NFLUIDS, 0);
  if (stats->volume == nullptr) goto label_end;
  st_stats_reset(stats);

  /* Set the error return code */

  error = 0;

  /* Returned arguments */

  label_end: if (error) stats = st_stats_undefine(stats);
  return (stats);
}

/****************************************************************************/
/*!
 **  Print the statistics
 **
 ** \param[in]  title    Title
 ** \param[in]  stats    Eden_Stats structure
 **
 *****************************************************************************/
static void st_stats_print(const char *title, Eden_Stats *stats)

{
  int ifacies, ifluid, number, totnum;
  double volume, totvol;

  /* Print the title */

  message("- %s\n", title);

  /* Print the statistics */

  totnum = 0;
  totvol = 0.;
  for (ifluid = 0; ifluid < NFLUIDS; ifluid++)
    for (ifacies = 0; ifacies < NFACIES; ifacies++)
    {
      number = STAT_NUM(ifacies, ifluid);
      volume = STAT_VOL(ifacies, ifluid);
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

  if (stats->ncork > 0) message("  . Cork                = %d\n", stats->ncork);

  return;
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
static int get_FACIES(int iech)
{
  int ifacies;

  ifacies = (int) DB->getArray(iech, IND_FACIES);
  if (ifacies < 0 || ifacies > NFACIES || IFFFF(ifacies)) ifacies = SHALE;
  return (ifacies);
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
static int get_FLUID(int iech)
{
  int ifluid;

  ifluid = (int) DB->getArray(iech, IPTR_FLUID);
  if (ifluid < 0 || ifluid > NFLUIDS || IFFFF(ifluid)) ifluid = UNDEF_FLUID;
  return (ifluid);
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
static int get_PERM(int iech)
{
  double perm;

  if (IND_PERM <= 0) return (1);
  perm = DB->getArray(iech, IND_PERM);
  if (FFFF(perm) || perm < 0.) perm = 0.;
  return ((int) perm);
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
static double get_PORO(int iech)
{
  double poro;

  if (IND_PORO <= 0) return (1);
  poro = DB->getArray(iech, IND_PORO);
  if (FFFF(poro)) return (0);
  poro = MIN(1., MAX(0., poro));
  return (poro);
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
static double get_DATE(int iech)
{
  double date;

  if (IND_DATE <= 0) return (0);
  date = DB->getArray(iech, IND_DATE);
  if (FFFF(date)) return (0);
  date = MAX(1., date);
  return (date);
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
static int ALREADY_FILLED(int ipos)

{
  return (get_FACIES(ipos) > 0 && get_PERM(ipos) > 0
          && get_FLUID(ipos) != UNDEF_FLUID);
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
static int TO_BE_FILLED(int ipos)

{
  return (get_FACIES(ipos) > 0 && get_PERM(ipos) > 0
          && get_FLUID(ipos) == UNDEF_FLUID);
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
static int get_FLUID_OLD(int iech)
{
  double ifluid;

  ifluid = DB->getArray(iech, IND_FLUID);
  if (ifluid < 0 || ifluid > NFLUIDS) ifluid = UNDEF_FLUID;
  return ((int) ifluid);
}

/****************************************************************************/
/*!
 **  Set the Facies value for a grid node
 **
 ** \param[in]  iech    Rank of the grid node
 ** \param[in]  ifacies Facies value
 **
 *****************************************************************************/
static void set_FACIES(int iech, int ifacies)
{
  DB->setArray(iech, IND_FACIES, ifacies);
  return;
}

/****************************************************************************/
/*!
 **  Turn the Facies into Cork
 **
 ** \param[in]  iech   Rank of the grid node
 **
 *****************************************************************************/
static void set_FACIES_CORK(int iech)

{
  int ifacies;

  ifacies = (int) DB->getArray(iech, IND_FACIES);
  DB->setArray(iech, IND_FACIES, -ifacies);
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
static void set_FLUID(int iech, int ifluid)
{
  DB->setArray(iech, IPTR_FLUID, ifluid);
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
static void set_DATE(int iech, int idate)
{
  double value;

  value = (IFFFF(idate)) ? TEST :
                           idate;
  DB->setArray(iech, IPTR_DATE, value);
  return;
}

/*****************************************************************************/
/*!
 **  Initialize the Eden_Stats structure
 **
 ** \param[in]  stats  Eden_Stats structure to be reset
 **
 *****************************************************************************/
static void st_stats_init(Eden_Stats *stats)

{
  int lec;

  for (lec = 0; lec < NXYZ; lec++)
  {
    if (!ALREADY_FILLED(lec)) continue;

    /* The cell does not belong to the skin: it is already filled */

    STAT_NUM(get_FACIES(lec)-1,get_FLUID(lec)-1) += 1;
    STAT_VOL(get_FACIES(lec)-1,get_FLUID(lec)-1) += get_PORO(lec);
  }
}

/****************************************************************************/
/*!
 **  Check that the Maximum quantities have been reached
 **
 *****************************************************************************/
static int st_check_max(Eden_Stats *stats, double number_max, double volume_max)
{
  int totnum, number;
  double totvol, volume;

  if (FFFF(number_max) && FFFF(volume_max)) return (0);

  /* Print the statistics */

  totnum = 0;
  totvol = 0.;
  for (int ifluid = 0; ifluid < NFLUIDS; ifluid++)
    for (int ifacies = 0; ifacies < NFACIES; ifacies++)
    {
      number = STAT_NUM(ifacies, ifluid);
      volume = STAT_VOL(ifacies, ifluid);
      totnum += number;
      totvol += volume;
      if (!FFFF(number_max) && totnum >= number_max) return (1);
      if (!FFFF(volume_max) && totvol >= volume_max) return (1);
    }
  return (0);
}

/****************************************************************************/
/*!
 **  Initialize and check Facies & Fluid matrices
 **
 *****************************************************************************/
static void st_check_inconsistency(void)

{
  int iech, ifacies, ifluid, perm, n_shale_fluid;

  /* Loop on the cells of the matrix */

  n_shale_fluid = 0;
  for (iech = 0; iech < NXYZ; iech++)
  {
    ifluid = get_FLUID_OLD(iech);
    ifacies = get_FACIES(iech);
    perm = get_PERM(iech);

    if (ifacies == SHALE || perm <= 0)
    {
      if (ifluid > 0)
      {
        if (VERBOSE)
          messerr(
              "Cell %d: Inconsistent Fluid (%d) with Facies (%d) or Perm (%d) -> set to %d",
              iech + 1, ifluid, ifacies, perm, NO_FLUID);
        n_shale_fluid++;
      }
      set_FLUID(iech, NO_FLUID);
      set_FACIES(iech, SHALE);
      set_DATE(iech, ITEST);
    }
    else
    {
      set_FLUID(iech, ifluid);
      set_FACIES(iech, ifacies);
      set_DATE(iech, (ifluid > 0));
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
static void st_update_results(int reset_facies, int show_fluid)

{
  int iech, ifluid, ifacies;

  /* Loop on the cells of the matrix */

  for (iech = 0; iech < NXYZ; iech++)
  {
    ifluid = get_FLUID_OLD(iech);
    ifacies = (int) DB->getArray(iech, IND_FACIES);

    /* Update the Facies information */

    if (ifacies < 0)
    {
      if (reset_facies)
        set_FACIES(iech, -ifacies);
      else
        set_FACIES(iech, CORK_FACIES);
    }

    /* Update the Fluid information */

    if (show_fluid)
    {
      if (ifluid == NO_FLUID || ifluid == UNDEF_FLUID || ifluid == CORK_FLUID)
        continue;
      set_FLUID(iech, ifluid + NFACIES * NFLUIDS);
    }
    else
    {
      if (ifluid == CORK_FLUID) set_FLUID(iech, UNDEF_FLUID);
    }
  }

  return;
}

/****************************************************************************/
/*!
 **  Calculate the statistics on Fluids and Corks
 **
 *****************************************************************************/
static void st_calculate_cumul(void)

{
  int iech, ifluid, ifacies;

  /* Loop on the cells of the matrix */

  for (iech = 0; iech < NXYZ; iech++)
  {

    /* Update the Fluid statistics */

    ifluid = get_FLUID(iech);
    if (ifluid > 0) DB->updArray(iech, IPTR_STAT_FLUID + ifluid - 1, 0, 1);

    /* Update the Cork statistics */

    ifacies = (int) DB->getArray(iech, IND_FACIES);
    if (ifacies < 0) DB->updArray(iech, IPTR_STAT_CORK, 0, 1);
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
static void st_normalize_cumul(int niter)

{
  int iech, ifluid;

  /* Loop on the cells of the matrix */

  for (iech = 0; iech < NXYZ; iech++)
  {

    /* Normalize the Fluid statistics */

    for (ifluid = 0; ifluid < NFLUIDS; ifluid++)
      DB->updArray(iech, IPTR_STAT_FLUID + ifluid, 3, (double) niter);

    /* Update the Cork statistics */

    DB->updArray(iech, IPTR_STAT_CORK, 3, (double) niter);
  }

  return;
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
 ** \param[in]  idir    Direciton value
 **
 *****************************************************************************/
static int WT(int ifacies, int ifluid, int perm, int idir)
{
  int ind, value;

  if (SPEEDS == nullptr)
    value = 1;
  else
  {
    ind = (idir) + 6 * ((ifacies - 1) * NFLUIDS + (ifluid - 1));
    value = perm * SPEEDS[ind];
  }
  return (value);
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
static double GET_WEIGHT(int ipos, int idir)
{
  int ind, ifacies, ifluid, perm;
  double value;

  if (SPEEDS == nullptr)
    value = 1.;
  else
  {
    ifacies = get_FACIES(ipos);
    ifluid = get_FLUID(ipos);
    perm = get_PERM(ipos);
    ind = (idir) + 6 * ((ifacies - 1) * NFLUIDS + (ifluid - 1));
    value = perm * SPEEDS[ind];
  }
  return (value);
}

/****************************************************************************/
/*!
 **  Print the parameters of the fluid propagation simulation
 **
 *****************************************************************************/
static void st_print_params(void)

{
  int ifacies, ifluid, idir;

  if (!VERBOSE) return;

  mestitle(0, "Fluid propagation parameters");
  message("Number of facies = %d\n", NFACIES);
  message("Number of fluids = %d\n", NFLUIDS);

  for (ifacies = 0; ifacies < NFACIES; ifacies++)
    for (ifluid = 0; ifluid < NFLUIDS; ifluid++)
    {
      message("Facies=%d - Fluid=%d -", ifacies + 1, ifluid + 1);
      for (idir = 0; idir < 6; idir++)
        message(" %s=%d", dirchar[idir], WT(ifacies + 1, ifluid + 1, 1, idir));
      message("\n");
    }
  return;
}

/****************************************************************************/
/*!
 **  Print the statistics for the cells not filled
 **
 ** \param[in]  title    Title
 **
 *****************************************************************************/
static void st_stats_empty(const char *title)

{
  int ifacies, i, number, total, flag_title;

  /* Print the statistics */

  total = 0;
  flag_title = 1;
  for (ifacies = 0; ifacies < NFACIES; ifacies++)
  {
    for (i = number = 0; i < NXYZ; i++)
    {
      if (!TO_BE_FILLED(i)) continue;
      if (get_FACIES(i) == (ifacies + 1)) number++;
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
static int st_fluid_modify(Skin *skin, int ipos, int *ref_fluid_loc)
{
  int ecr, dir, fluid, ref_fluid;

  /* Initializations */

  ref_fluid = UNDEF_FLUID;

  /* Loop on the directions */

  for (dir = 0; dir < 6; dir++)
  {
    if (skin_grid_shift(skin, ipos, dir, &ecr))
    {
      if (ALREADY_FILLED(ecr))
      {
        fluid = get_FLUID(ecr);

        if (ref_fluid == UNDEF_FLUID)
        {
          if (dir == DIR_DOWN)
          {
            if (WT(get_FACIES(ecr), get_FLUID(ecr), get_PERM(ecr), invdir[dir]) > 0)
              ref_fluid = fluid;
          }
          else if (dir == DIR_UP)
          {
            if (WT(get_FACIES(ecr), get_FLUID(ecr), get_PERM(ecr), invdir[dir]) > 0)
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
 **  Check the validity of the array Speed
 **
 ** \return  Error return code
 **
 *****************************************************************************/
static int st_fluid_check(void)

{
  int ifacies, ifluid, jfluid, idir, error, total;

  /* Initializations */

  error = 0;

  /* Check that there is no zero value */

  for (ifluid = 0; ifluid < NFLUIDS; ifluid++)
    for (ifacies = 0; ifacies < NFACIES; ifacies++)
      for (idir = 0; idir < 6; idir++)
      {
        if (WT(ifacies + 1, ifluid + 1, 1, idir) <= 0)
        {
          messerr(
              "The Propagation Directional Speed is zero for: Fluid=%d - Facies=%d - Direction='%s'",
              ifluid + 1, ifacies + 1, dirchar[idir]);
          messerr("This may cause artifacts. Change it to a low value instead");
          error = 1;
        }
      }

  /* Check that at least one speed is defined (for each facies/fluid pair) */

  for (ifluid = 0; ifluid < NFLUIDS; ifluid++)
    for (ifacies = 0; ifacies < NFACIES; ifacies++)
    {
      total = 0;
      for (idir = 0; idir < 6; idir++)
        total += WT(ifacies + 1, ifluid + 1, 1, idir);
      if (total <= 0.)
      {
        messerr("For Facies (%d) and Fluid (%d), no positive speed is defined",
                ifacies + 1, ifluid + 1);
        error = 1;
      }
    }

  /* Check for wrong order relationship for velocities */

  for (ifacies = 0; ifacies < NFACIES; ifacies++)
    for (ifluid = 1; ifluid < NFLUIDS; ifluid++)
    {
      jfluid = ifluid - 1;

      /* Z+ Speed */

      if (WT(ifacies + 1, jfluid + 1, 1, DIR_UP) < WT(ifacies + 1, ifluid + 1,
                                                      1, DIR_UP))
      {
        messerr("Error for the Z+ Propagation Directional Speed for Facies=%d:",
                ifacies + 1);
        messerr(
            "Speed for Fluid=%d [%d] must not be smaller than Speed for Fluid=%d [%d]",
            jfluid + 1, WT(ifacies + 1, jfluid + 1, 1, DIR_UP), ifluid + 1,
            WT(ifacies + 1, ifluid + 1, 1, DIR_UP));
        error = 1;
      }

      /* Z- Speed */

      if (WT(ifacies + 1, jfluid + 1, 1, DIR_DOWN) > WT(ifacies + 1, ifluid + 1,
                                                        1, DIR_DOWN))
      {
        messerr("Error for the Z- Propagation Directional Speed for Facies=%d:",
                ifacies + 1);
        messerr(
            "Speed for Fluid=%d [%d] must not be larger than Speed for Fluid=%d  [%d]",
            jfluid + 1, WT(ifacies + 1, jfluid + 1, 1, DIR_DOWN), ifluid + 1,
            WT(ifacies + 1, ifluid + 1, 1, DIR_DOWN));
        error = 1;
      }
    }
  return (error);
}

/*****************************************************************************/
/*!
 **  Multivariate multiphase propagation into a set of components
 **  constrained by initial conditions and fluid densities
 **
 ** \return  Error return code : 1 no fluid to propagate
 **
 ** \param[in]  dbgrid        Db grid structure
 ** \param[in]  verbose       1 for a verbose option
 ** \param[in]  seed          Seed for random number generator (or 0)
 ** \param[in]  niter         Number of iterations
 ** \param[in]  ind_facies    Rank of the variable containing the Facies
 ** \param[in]  ind_fluid     Rank of the variable containing the Fluid
 ** \param[in]  ind_perm      Rank of the variable containing the Permeability
 ** \param[in]  ind_poro      Rank of the variable containing the Porosity
 ** \param[in]  nfacies       number of facies (facies 0 excluded)
 ** \param[in]  nfluids       number of fluids
 ** \param[in]  speeds        array containing the travel speeds
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
int fluid_propagation(Db *dbgrid,
                                      int verbose,
                                      int seed,
                                      int niter,
                                      int ind_facies,
                                      int ind_fluid,
                                      int ind_perm,
                                      int ind_poro,
                                      int nfacies,
                                      int nfluids,
                                      int *speeds,
                                      int show_fluid,
                                      double number_max,
                                      double volume_max)
{
  Eden_Stats *stats;
  Skin *skin;
  int ipos, rank, error, seed_memo, ref_fluid, idate, iter;

  /* Initializations */

  error = 1;
  stats = nullptr;

  /* Preliminary checks */

  if (!is_grid(dbgrid))
  {
    messerr("The Fluid Propagation is restricted to regular grid");
    return (1);
  }
  if (dbgrid->getNDim() > 3)
  {
    messerr("Fluid propagation is limited to 3-D space (maximum)");
    return (1);
  }
  if (ind_facies < 0 || ind_facies > dbgrid->getFieldNumber() || ind_fluid < 0
      || ind_fluid > dbgrid->getFieldNumber())
  {
    messerr("Error in the ranks of the facies (%d) and fluid (%d) variables",
            ind_facies, ind_fluid);
    return (1);
  }

  /* Define global variables */

  DB = dbgrid;
  NXYZ = DB->getSampleNumber();
  NFACIES = nfacies;
  NFLUIDS = nfluids;
  SPEEDS = speeds;
  IND_FACIES = ind_facies;
  IND_FLUID = ind_fluid;
  IND_PERM = ind_perm;
  IND_PORO = ind_poro;
  VERBOSE = verbose;
  skin = nullptr;

  /* Preliminary checks */

  if (st_fluid_check()) goto label_end;

  /* Printout of the fluid propagation parameters */

  st_print_params();

  /* Core allocation */

  law_set_random_seed(seed);
  stats = st_stats_define();
  skin = skin_define(dbgrid, ALREADY_FILLED, TO_BE_FILLED, GET_WEIGHT);
  if (skin == nullptr) goto label_end;

  /* Add the attributes for storing the Fluid and Cork statistics */

  if (niter > 1)
  {
    IPTR_STAT_FLUID = DB->addFields(NFLUIDS, 0.);
    if (IPTR_STAT_FLUID < 0) goto label_end;
    IPTR_STAT_CORK = DB->addFields(1, 0.);
    if (IPTR_STAT_CORK < 0) goto label_end;
  }

  /* Add the attributes for storing the Fluid and Data informations */

  IPTR_FLUID = DB->addFields(1, 0.);
  if (IPTR_FLUID < 0) goto label_end;
  IPTR_DATE = DB->addFields(1, TEST);
  if (IPTR_DATE < 0) goto label_end;

  /* Loop on the iterations */

  for (iter = 0; iter < niter; iter++)
  {
    seed_memo = law_get_random_seed();
    st_stats_reset(stats);

    /* Check the consistency */

    st_check_inconsistency();

    /* Initialize the grid with the initial values */

    st_stats_init(stats);
    if (skin_init(skin, VERBOSE))
    {
      error = 0;
      goto label_end;
    }

    /* Modifying the peripheral cells using a random walk */

    idate = 0;
    while (skin_remains(skin))
    {

      /* Check that the maximum quantities have not been reached */

      if (st_check_max(stats, number_max, volume_max)) break;
      idate++;

      /* Find the next cell to be processed */

      skin_next(skin, &rank, &ipos);

      /* Find the new value of the target cell according to its neighborhood */

      if (st_fluid_modify(skin, ipos, &ref_fluid))
      {
        stats->ncork++;
        set_FACIES_CORK(ipos);
        set_FLUID(ipos, CORK_FLUID);
        set_DATE(ipos, ITEST);
      }
      else
      {
        STAT_NUM(get_FACIES(ipos)-1,ref_fluid-1) += 1;
        STAT_VOL(get_FACIES(ipos)-1,ref_fluid-1) += get_PORO(ipos);
        set_FLUID(ipos, ref_fluid);
        set_DATE(ipos, idate);
      }

      /* Deduce the initial influence of the central cell */

      if (skin_unstack(skin, rank, ipos)) goto label_end;
    }

    /* Final printout */

    if (VERBOSE)
    {
      mestitle(1, "Final status (iteration %d)", iter + 1);
      message("- Seed Value                     = %d\n", seed_memo);
      st_stats_print("Cells already filled", stats);
      st_stats_empty("Cells not reached");
    }

    /* Calculate statistics on fluids and corks */

    if (niter > 1) st_calculate_cumul();

    /* Update the data (optional) */

    st_update_results(iter < niter - 1, show_fluid);
  }

  /* Normalize the statistics */

  if (niter > 1) st_normalize_cumul(niter);

  /* Set the error return flag */

  if (VERBOSE) skin_print(skin);
  error = 0;

  label_end: stats = st_stats_undefine(stats);
  skin = skin_undefine(skin);
  if (niter > 1)
  {
    if (IPTR_FLUID > 0) DB->deleteFieldByAttribute(IPTR_FLUID);
    if (IPTR_DATE > 0) DB->deleteFieldByAttribute(IPTR_DATE);
  }
  return (error);
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
static int st_get_time_interval(double date,
                                int ntime,
                                double time0,
                                double dtime)
{
  double time_deb, time_fin;
  int itime;

  /* Loop on the intervals */

  for (itime = 0; itime < ntime; itime++)
  {
    time_deb = time0 + dtime * itime;
    time_fin = time0 + dtime * (itime + 1);
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
 ** \param[in]  dbgrid        Db grid structure
 ** \param[in]  verbose       1 for a verbose option
 ** \param[in]  ind_date      Rank of variable containing Date
 ** \param[in]  ind_facies    Rank of variable containing Facies
 ** \param[in]  ind_fluid     Rank of variable containing Fluid
 ** \param[in]  ind_poro      Rank of variable containing Porosity (optional)
 ** \param[in]  nfacies       number of facies (facies 0 excluded)
 ** \param[in]  nfluids       number of fluids
 ** \param[in]  facies0       Value of the target facies
 ** \param[in]  fluid0        Value of the target fluid
 ** \param[in]  ntime         Number of Time intervals
 ** \param[in]  time0         Starting time
 ** \param[in]  dtime         Time interval
 ** \param[in]  tab           Array of extracted statistics
 **
 *****************************************************************************/
int fluid_extract(Db *dbgrid,
                                  int verbose,
                                  int ind_date,
                                  int ind_facies,
                                  int ind_fluid,
                                  int ind_poro,
                                  int nfacies,
                                  int nfluids,
                                  int facies0,
                                  int fluid0,
                                  int ntime,
                                  double time0,
                                  double dtime,
                                  double *tab)
{
  int itime, iech;
  double totnum, totvol, locnum, locvol, volume, datmax, date;

  /* Preliminary checks */

  if (!is_grid(dbgrid))
  {
    messerr("The Fluid Propagation is restricted to regular grid");
    return (1);
  }
  if (dbgrid->getNDim() > 3)
  {
    messerr("Fluid propagation is limited to 3-D space (maximum)");
    return (1);
  }
  if (ind_facies < 0 || ind_facies > dbgrid->getFieldNumber() || ind_fluid < 0
      || ind_fluid > dbgrid->getFieldNumber() || ind_date < 0
      || ind_fluid > dbgrid->getFieldNumber())
  {
    messerr("Error in ranks of facies (%d), fluid (%d) or date (%d) variables",
            ind_facies, ind_fluid, ind_date);
    return (1);
  }
  if (ntime < 0 || time0 < 0 || dtime <= 0)
  {
    messerr("Error in Time Interval Definition");
    messerr("Origin=%lf - Step=%lf - Number=%d", time0, dtime, ntime);
    return (1);
  }

  /* Define global variables */

  DB = dbgrid;
  NXYZ = DB->getSampleNumber();
  NFACIES = nfacies;
  NFLUIDS = nfluids;
  IND_FACIES = ind_facies;
  IND_FLUID = ind_fluid;
  IND_PORO = ind_poro;
  IND_DATE = ind_date;
  VERBOSE = verbose;

  /* Initialize the array */

  for (itime = 0; itime < ntime; itime++)
  {
    TAB(0,itime) = time0 + dtime * itime;
    TAB(1,itime) = time0 + dtime * (itime + 1);
    TAB(2,itime) = 0.;
    TAB(3,itime) = 0.;
  }

  /* Loop on the blocks */

  totnum = totvol = locnum = locvol = datmax = 0;
  for (iech = 0; iech < NXYZ; iech++)
  {

    if (get_FACIES(iech) != facies0) continue;
    if (get_FLUID(iech) != fluid0) continue;
    volume = get_PORO(iech);
    date = get_DATE(iech);
    if (date > datmax) datmax = date;

    totnum += 1;
    totvol += volume;
    itime = st_get_time_interval(date, ntime, time0, dtime);
    if (itime < 0) continue;
    locnum += 1;
    locvol += volume;

    TAB(2,itime) += 1;
    TAB(3,itime) += volume;
  }

  /* Final printout */

  if (VERBOSE)
  {
    mestitle(1, "Extraction for Fluid(%d) and Facies(%d)", facies0, fluid0);
    message("Time slices: From %lf to %lf by step of %lf\n", time0,
            time0 + dtime * ntime, dtime);
    message("Total Number of Cells               = %d\n", NXYZ);
    message("Maximum Date                        = %lf\n", datmax);
    message("Total Number of Invaded Cells       = %lf\n", totnum);
    message("Total Volume of Invaded Cells       = %lf\n", totvol);
    message("Total Number of Cells in Time Slice = %lf\n", locnum);
    message("Total Volume of Cells in Time Slice = %lf\n", locvol);
  }

  return (0);
}
