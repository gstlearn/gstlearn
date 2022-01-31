/******************************************************************************/
/* COPYRIGHT ARMINES, ALL RIGHTS RESERVED                                     */
/*                                                                            */
/* THE CONTENT OF THIS WORK CONTAINS CONFIDENTIAL AND PROPRIETARY             */
/* INFORMATION OF ARMINES. ANY DUPLICATION, MODIFICATION,                     */
/* DISTRIBUTION, OR DISCLOSURE IN ANY FORM, IN WHOLE, OR IN PART, IS STRICTLY */
/* PROHIBITED WITHOUT THE PRIOR EXPRESS WRITT%EN PERMISSION OF ARMINES        */
/*                                                                            */
/* TAG_SOURCE_CG                                                              */
/******************************************************************************/
#include "geoslib_f.h"
#include "geoslib_old_f.h"
#include "Basic/Law.hpp"
#include "Basic/OptDbg.hpp"
#include "Db/Db.hpp"

#include <math.h>

/*! \cond */
#define SKIN_QUANT  1000
/*! \endcond */

static int (*FUNC_ALREADY_DONE)(int ipos);
static int (*FUNC_TO_BE_DONE)(int ipos);
static double (*FUNC_GET_WEIGHT)(int ipos, int idir);

static int ndir[4] = { 0, 2, 4, 6 };
static int invdir[6] = { 1, 0, 3, 2, 5, 4 };
static int id[6][3] = { { 1, 0, 0 },
                        { -1, 0, 0 },
                        { 0, 1, 0 },
                        { 0, -1, 0 },
                        { 0, 0, 1 },
                        { 0, 0, -1 } };

/****************************************************************************/
/*!
 **  Returns the weight for a given cell and direction
 **
 ** \return  The weight
 **
 ** \param[in]  ipos  Absolute grid index of the input grid node
 ** \param[in]  dir   Rank of the direction
 **
 *****************************************************************************/
static double st_get_weight(int ipos, int dir)
{
  double value;

  if (FUNC_GET_WEIGHT != NULL)
    value = FUNC_GET_WEIGHT(ipos, dir);
  else
    value = 1.;
  return (value);
}

/****************************************************************************/
/*!
 **  Returns the shifted node of a skin
 **
 ** \return  1 if the shifted node belongs to the grid and 0 otherwise
 **
 ** \param[in]  skin   Skin structure
 ** \param[in]  indg0  Array of directional grid indices
 ** \param[in]  dir    Rank of the direction
 **
 ** \param[out]  iad   Absolute sample address
 **
 *****************************************************************************/
static int st_skin_grid_shift(Skin *skin, int indg0[3], int dir, int *iad)
{
  int i, indg[3];

  /* Shift the target grid node and check if it belongs to the grid */

  for (i = 0; i < skin->ndim; i++)
  {
    indg[i] = indg0[i] + id[dir][i];
    if (indg[i] < 0 || indg[i] >= skin->db->getNX(i)) return (0);
  }

  /* Calculate the absolute address of the shifted grid node */

  *iad = db_index_grid_to_sample(skin->db, indg);
  return (1);
}

/****************************************************************************/
/*!
 **  Returns the shifted node of a skin
 **
 ** \return  1 if the shifted node belongs to the grid and 0 otherwise
 **
 ** \param[in]  skin  Skin structure
 ** \param[in]  lec   Absolute grid index of the input grid node
 ** \param[in]  dir   Rank of the direction
 **
 ** \param[out]  iad  Absolute sample address
 **
 *****************************************************************************/
int skin_grid_shift(Skin *skin, int lec, int dir, int *iad)
{
  int indg[3];

  /* Convert an absolute address into the grid indices */

  db_index_sample_to_grid(skin->db, lec, indg);

  /* Shift the target grid node and check if it belongs to the grid */

  if (!st_skin_grid_shift(skin, indg, dir, iad)) return (0);

  return (1);
}

/*****************************************************************************/
/*!
 **  Allocate the skin of the body
 **
 ** \return  Pointer to the allocated skin (or NULL)
 **
 ** \param[in]  db  Db structure
 ** \param[in]  func_already_done  Check if a cell is already processed
 ** \param[in]  func_to_be_done    Check if the cell can be processed
 ** \param[in]  func_get_weight    Returns the weight of a cell in a direction
 **
 ** \remark  If func_get_weight is not defined, the weight is set to 1
 **
 *****************************************************************************/
Skin* skin_define(DbGrid *db,
                  int (*func_already_done)(int ipos),
                  int (*func_to_be_done)(int ipos),
                  double (*func_get_weight)(int ipos, int dir))
{
  Skin *skin;
  int error;

  /* Initializations */

  error = 1;
  FUNC_ALREADY_DONE = func_already_done;
  FUNC_TO_BE_DONE = func_to_be_done;
  FUNC_GET_WEIGHT = func_get_weight;

  /* Allocation */

  skin = (Skin*) mem_alloc(sizeof(Skin), 0);
  if (skin == nullptr) goto label_end;
  skin->db = db;
  skin->nxyz = db->getSampleNumber();
  skin->ndim = db->getNDim();
  skin->nval = 0;
  skin->size = 0;
  skin->date = 0;
  skin->nalloc = 0;
  skin->nval_max = 0;
  skin->total = 0.;
  skin->total_max = 0.;
  skin->quant = MIN(SKIN_QUANT, (int ) (sqrt(skin->nxyz) / 10));
  skin->address = nullptr;
  skin->energy = nullptr;

  /* Set the error return code */

  error = 0;

  /* Returned arguments */

  label_end: if (error) skin = skin_undefine(skin);
  return (skin);
}

/*****************************************************************************/
/*!
 **  Deallocate the Skin
 **
 ** \return  Pointer to the newly deallocated skin
 **
 ** \param[in]  skin    Pointer to the skin to be deallocated
 **
 *****************************************************************************/
Skin* skin_undefine(Skin *skin)

{

  /* Deallocation */

  if (skin != nullptr)
  {
    skin->address = (int*) mem_free((char* ) skin->address);
    skin->energy = (double*) mem_free((char* ) skin->energy);
  }
  skin = (Skin*) mem_free((char* ) skin);

  return (skin);
}

/*****************************************************************************/
/*!
 **  Delete a cell from the skin
 **
 ** \return  Error returned code
 **
 ** \param[in]  skin     Pointer to the skin
 ** \param[in]  rank     Rank of the cell to be deleted
 **
 ** \remark  When deleting a cell from the skin, the last cell is copied
 ** \remark  in the place of the deleted one
 **
 *****************************************************************************/
static int st_skin_cell_del(Skin *skin, int rank)
{
  /* Delete the target cell : move the last cell to the target location */

  skin->nval--;
  skin->address[rank] = skin->address[skin->nval];
  skin->energy[rank] = skin->energy[skin->nval];

  /* Deallocate complementary room in the skin */

  if (skin->size - skin->nval > skin->quant)
  {
    skin->size -= skin->quant;
    if (skin->size > 0)
    {
      skin->address = (int*) mem_realloc((char* ) skin->address,
                                         sizeof(int) * skin->size, 0);
      if (skin->address == nullptr) return (1);
      skin->energy = (double*) mem_realloc((char* ) skin->energy,
                                           sizeof(double) * skin->size, 0);
      if (skin->energy == nullptr) return (1);
    }
    else
    {
      skin->address = (int*) mem_free((char* ) skin->address);
      skin->energy = (double*) mem_free((char* ) skin->energy);
    }
  }
  return (0);
}

/*****************************************************************************/
/*!
 **  Checks if a cell already belongs to the skin
 **
 ** \return  Rank within the skin or -1 if the cell does not belong to the skin
 **
 ** \param[in]  skin     Pointer to the skin
 ** \param[in]  ipos     Cell location
 **
 *****************************************************************************/
static int st_skin_cell_already(Skin *skin, int ipos)
{
  int i;

  /* Check if the cell already belongs to the skin */

  for (i = 0; i < skin->nval; i++)
    if (skin->address[i] == ipos) return (i);

  return (-1);
}

/*****************************************************************************/
/*!
 **  Modify the energy for a cell which already belongs to the skin
 **
 ** \param[in]  skin     Pointer to the skin
 ** \param[out] rank     Location of the cell in the skin
 ** \param[in]  energy   Additional energy for the new cell
 **
 *****************************************************************************/
static void st_skin_cell_mod(Skin *skin, int rank, double energy)
{
  skin->energy[rank] += energy;
  skin->total += energy;
  if (skin->total > skin->total_max) skin->total_max = skin->total;
  return;
}

/*****************************************************************************/
/*!
 **  Add a cell to the skin (if not already in the skin)
 **
 ** \return  Error returned code
 **
 ** \param[in]  skin     Pointer to the skin
 ** \param[in]  ipos     Cell location
 ** \param[in]  energy   Energy for the new cell
 **
 *****************************************************************************/
static int st_skin_cell_add(Skin *skin, int ipos, double energy)
{
  int rank;

  /* Allocate complementary room in the skin */

  if (skin->nval >= skin->size)
  {
    skin->nalloc++;
    if (skin->size > 0)
    {
      skin->size += skin->quant;
      skin->address = (int*) mem_realloc((char* ) skin->address,
                                         sizeof(int) * skin->size, 0);
      if (skin->address == nullptr) return (1);
      skin->energy = (double*) mem_realloc((char* ) skin->energy,
                                           sizeof(double) * skin->size, 0);
      if (skin->energy == nullptr) return (1);
    }
    else
    {
      skin->size = skin->quant;
      skin->address = (int*) mem_alloc(sizeof(int) * skin->size, 0);
      if (skin->address == nullptr) return (1);
      skin->energy = (double*) mem_alloc(sizeof(double) * skin->size, 0);
      if (skin->energy == nullptr) return (1);
    }
  }
  rank = skin->nval;
  skin->address[rank] = ipos;
  skin->energy[rank] = 0.;
  skin->nval++;
  if (skin->nval > skin->nval_max) skin->nval_max = skin->nval;

  /* Upgrade the energy */

  st_skin_cell_mod(skin, rank, energy);

  return (0);
}

/*****************************************************************************/
/*!
 **  Initialize the skin
 **
 ** \return  Error returned code
 **
 ** \param[in] skin     Skin structure
 ** \param[in] verbose  Verbose flag
 **
 *****************************************************************************/
int skin_init(Skin *skin, int verbose)
{
  int lec, ecr, dir, nb_count, nb_done, nb_mask, total, indg[3];
  double local;

  /* Loop on all the cells */

  nb_mask = nb_count = nb_done = 0;
  total = skin->nxyz;
  for (lec = 0; lec < total; lec++)
  {
    if (FUNC_ALREADY_DONE(lec))
    {

      /* The cell does not belong to the skin: it is already filled */

      nb_done++;
      continue;
    }
    else if (!FUNC_TO_BE_DONE(lec))
    {

      /* The cell does not belong to the cell: it is masked off */

      nb_mask++;
      continue;
    }
    else
    {

      /* The cell is elligible */

      nb_count++;
      local = 0.;
      db_index_sample_to_grid(skin->db, lec, indg);
      for (dir = 0; dir < ndir[skin->ndim]; dir++)
      {
        if (!st_skin_grid_shift(skin, indg, dir, &ecr)) continue;
        if (!FUNC_ALREADY_DONE(ecr)) continue;
        local += st_get_weight(ecr, invdir[dir]);
      }
      if (local > 0.)
      {
        if (st_skin_cell_add(skin, lec, local))
        {
          messerr("Core allocation problem in Skin algorithm");
          return (1);
        }
      }
    }
  }

  /* Print the statistics */

  if (verbose)
  {
    mestitle(1, "Skin algorithm: Initial status");
    message("- Total number of cells           = %d\n", total);
    message("- Number of cells already filled  = %d\n", nb_done);
    message("- Number of cells active          = %d\n", total - nb_mask);
    message("- Number of cells to be processed = %d\n", nb_count);
  }

  if (nb_count <= 0 || skin->total <= 0.)
  {
    messerr("There is no cell to be processed");
    return (1);
  }

  return (0);
}

/*****************************************************************************/
/*!
 **  Returns the number of cells still to be processed
 **
 ** \return  Returns the number of cells still to be processed
 **
 ** \param[out] skin    Skin structure
 **
 *****************************************************************************/
int skin_remains(Skin *skin)

{
  skin->date++;
  if (OptDbg::query(EDbg::MORPHO))
    message("Skin iteration:%5d - Length:%4d - Energy:%lf\n", skin->date,
            skin->nval, skin->total);
  return ((int) skin->total);
}

/*****************************************************************************/
/*!
 **  Find the next cell at random within the skin
 **
 ** \param[in]  skin     Skin structure
 **
 ** \param[out] rank     Location of the cell in the skin
 ** \param[out] ipos     Cell location
 **
 *****************************************************************************/
void skin_next(Skin *skin, int *rank, int *ipos)
{
  int i;
  double tirage, total;

  /* Draw a random cell */

  tirage = skin->total * law_uniform(0., 1.);

  /* Find the cell */

  total = 0.;
  for (i = 0; i < skin->nval; i++)
  {
    total += skin->energy[i];
    if (total >= tirage)
    {
      *rank = i;
      *ipos = skin->address[i];
      if (!FUNC_TO_BE_DONE(*ipos))
        messageAbort(
            "Elligible cell (%d ipos=%d) of the skin is already filled", i,
            *ipos);
      return;
    }
  }

  messageAbort("Cannot find a cell for propagation");
  return;
}

/*****************************************************************************/
/*!
 **  Suppress the current cell from the skin
 **
 ** \return  Error return code
 **
 ** \param[in] skin     Skin structure
 ** \param[in] rank0    Rank of the current cell in the skin
 ** \param[in] ipos0    Cell location
 **
 *****************************************************************************/
int skin_unstack(Skin *skin, int rank0, int ipos0)
{
  int dir, ecr, rank, indg[3];
  double local;

  /* Suppress the current cell from the skin */

  skin->total -= skin->energy[rank0];
  if (st_skin_cell_del(skin, rank0)) return (1);

  /* Update the neighboring cells */

  db_index_sample_to_grid(skin->db, ipos0, indg);
  for (dir = 0; dir < ndir[skin->ndim]; dir++)
  {
    if (!st_skin_grid_shift(skin, indg, dir, &ecr)) continue;

    /* Discard the neighboring cell if it cannot filled */

    if (!FUNC_TO_BE_DONE(ecr)) continue;
    local = st_get_weight(ipos0, dir);
    rank = st_skin_cell_already(skin, ecr);
    if (rank < 0)
    {

      /* The cell does not already belong to the skin: add it */

      if (st_skin_cell_add(skin, ecr, local)) return (1);
    }
    else
    {
      /* If the cell already belongs to the skin, upgrade its energy */

      st_skin_cell_mod(skin, rank, local);
    }
  }
  return (0);
}

/*****************************************************************************/
/*!
 **  Print the computing information concerning the skin algorithm
 **
 ** \param[in] skin    Skin structure
 **
 *****************************************************************************/
void skin_print(Skin *skin)

{
  mestitle(1, "Skin algorithm: Final status");
  message("- Number of iterations          = %d\n", skin->date);
  message("- Maximum skin length           = %d\n", skin->nval_max);
  message("- Maximum energy                = %lf\n", skin->total_max);
  return;
}
