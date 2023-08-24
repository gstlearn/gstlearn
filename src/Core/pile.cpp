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
#include "geoslib_old_f.h"

static int DEBUG = 0;
static int MAX_PILE = 10;
static const char *PILE_NAME[] = { "Db",
                                   "Vario",
                                   "Model",
                                   "Neigh",
                                   "Rule",
                                   "AAnam",
                                   "Tokens",
                                   "Polygon",
                                   "Fracture",
                                   "PCA" };
static int MAX_COUNT[] = { 10, 4, 4, 2, 4, 1, 1, 1, 2, 2 };
static char ***piles;

/****************************************************************************/
/*!
 **  Check that the type is valid or abort if invalid
 **
 ** \param[in]  type    Type of the pile
 ** \param[in]  rank    Rank of the slot (-1 if not defined)
 **
 *****************************************************************************/
static void st_valid(int type, int rank)
{
  if (type < 0 || type >= MAX_PILE)
    messageAbort("Type error in Pile management: %d (>= %d)", type, MAX_PILE);
  if (rank < 0) return;
  if (rank >= MAX_COUNT[type])
    messageAbort("Slot error in Pile Management: rank(%d) >= max[type=%d](%d)",
                 rank, type, MAX_COUNT[type]);
}

/****************************************************************************/
/*!
 **  Reset one pile
 **
 ** \param[in]  type    Type of the pile
 **
 ** \remark  This function must be called compulsorily when entering RGeostats
 **
 *****************************************************************************/
void pile_reset(int type)

{
  int i;

  st_valid(type, -1);

  /* Reset the pile */

  for (i = 0; i < MAX_COUNT[type]; i++)
    piles[type][i] = NULL;

  return;
}

/****************************************************************************/
/*!
 **  Reset all the piles
 **
 ** \remark  This function must be called compulsorily when entering RGeostats
 **
 *****************************************************************************/
void piles_reset(void)

{
  int type;

  if (piles != NULL)
  {
    for (type = 0; type < MAX_PILE; type++)
      piles[type] = (char**) mem_free((char* ) piles[type]);
    piles = (char***) mem_free((char* ) piles);
  }

  piles = (char***) mem_alloc(sizeof(char**) * MAX_PILE, 1);

  for (type = 0; type < MAX_PILE; type++)
  {
    piles[type] = (char**) mem_alloc(sizeof(char*) * MAX_COUNT[type], 1);
    pile_reset(type);
  }

  return;
}

/****************************************************************************/
/*!
 **  Returns the rank of the first free element
 **
 ** \return  Rank of the first available slot or -1 if failure
 **
 ** \param[in]  type    Type of the pile
 **
 *****************************************************************************/
int pile_next(int type)

{
  int ival, i;

  st_valid(type, -1);

  /* Look for the first empty slot */

  ival = -1;
  for (i = 0; i < MAX_COUNT[type] && ival < 0; i++)
  {
    if (piles[type][i] == NULL) ival = i;
  }
  if (ival < 0)
  {
    messerr("No %s Slot available: all entries are used (%d)", PILE_NAME[type],
            MAX_COUNT[type]);
    messerr("All Entries are released for future use");
    piles_reset();
    return (-1);
  }

  return (ival);
}

/****************************************************************************/
/*!
 **  Manage the pile
 **
 ** \param[in]  type    Type of the pile
 ** \param[in]  rank    Rank of the slot
 ** \param[in]  mode    1 for allocation and -1 for deallocation
 ** \param[in]  ptr     Pointor to be stored (allocation mode)
 **
 *****************************************************************************/
void pile_manage(int type, int rank, int mode, char *ptr)
{
  st_valid(type, rank);

  /* Dispatch */

  if (mode > 0)
  {

    /* Allocation */

    piles[type][rank] = ptr;
  }
  else
  {

    /* Deallocation */

    piles[type][rank] = NULL;
  }
}

/****************************************************************************/
/*!
 **  Check if the Pile is managed correctly
 **
 ** \return  1 if the Pile Management is incorrect; 0 otherwise
 **
 ** \param[in]  type   Type of the pile
 ** \param[in]  rank   Rank of the slot
 ** \param[in]  mode   Status
 ** \li                  1 : Check that the element must be already allocated
 ** \li                 -1 : Check that the element will be allocated
 **
 *****************************************************************************/
int pile_correct(int type, int rank, int mode)
{
  if (piles == NULL)
  {
    messerr("The Piles have not been initialized");
    return (1);
  }

  st_valid(type, rank);

  /* Dispatch */

  if (rank < 0) return (1);

  switch (mode)
  {
    case 1: /* The element must be already allocated */
      if (piles[type][rank] != NULL) return (0);
      messerr("Error: the element (%d) of the pile (%s) is not allocated", rank,
              PILE_NAME[type]);
      break;

    case -1: /* The element will be allocated */
      if (!DEBUG) return (0);
      if (piles[type][rank] == NULL) return (0);
      messerr("Error: the element (%d) of the pile (%s) is already allocated",
              rank, PILE_NAME[type]);
      break;
  }
  return (1);
}

/****************************************************************************/
/*!
 **  Returns the address managed in the pile
 **
 ** \param[in]  type    Type of the pile
 ** \param[in]  rank    Rank of the slot
 **
 *****************************************************************************/
char* pile_get(int type, int rank)
{
  if (piles == NULL)
  {
    messerr("The Piles have not been initialized");
    return (NULL);
  }

  st_valid(type, rank);

  /* Returns the address */

  return (piles[type][rank]);
}

/****************************************************************************/
/*!
 **  Dump the piles contents
 **
 *****************************************************************************/
void piles_dump(void)

{
  int type, i;

  if (piles == NULL)
  {
    messerr("The Piles have not been initialized");
    return;
  }
  for (type = 0; type < MAX_PILE; type++)
  {
    if (piles[type] == NULL) continue;
    for (i = 0; i < MAX_COUNT[type]; i++)
    {
      if (piles[type][i] != NULL)
        message("Type %9s : Slot %2d (out of %2d) allocated\n", PILE_NAME[type],
                i + 1, MAX_COUNT[type]);
      else
        message("Type %9s : Slot %2d (out of %2d) free\n", PILE_NAME[type],
                i + 1, MAX_COUNT[type]);
    }
  }
}

