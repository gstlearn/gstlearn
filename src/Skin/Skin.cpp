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
#include "Skin/ISkinFunctions.hpp"
#include "Skin/Skin.hpp"
#include "Db/DbGrid.hpp"
#include "Basic/Law.hpp"
#include "Basic/Vector.hpp"

#include "math.h"

#define SKIN_QUANT  1000

static int ndir[4] = { 0, 2, 4, 6 };
static int invdir[6] = { 1, 0, 3, 2, 5, 4 };
static int id[6][3] = { { 1, 0, 0 },
                        { -1, 0, 0 },
                        { 0, 1, 0 },
                        { 0, -1, 0 },
                        { 0, 0, 1 },
                        { 0, 0, -1 } };

Skin::Skin(const ISkinFunctions* skf, DbGrid* dbgrid)
    : _skf(skf),
      _dbgrid(dbgrid),
      _nxyz(0),
      _ndim(0),
      _nval(0),
      _date(0),
      _nvalMax(0),
      _total(0.),
      _totalMax(0.),
      _address(),
      _energy()
{
  if (dbgrid != nullptr)
  {
    _nxyz = _dbgrid->getSampleNumber();
    _ndim = _dbgrid->getNDim();
  }
}

Skin::Skin(const Skin &r)
    : _skf(r._skf),
      _dbgrid(r._dbgrid),
      _nxyz(r._nxyz),
      _ndim(r._ndim),
      _nval(r._nval),
      _date(r._date),
      _nvalMax(r._nvalMax),
      _total(r._total),
      _totalMax(r._totalMax),
      _address(r._address),
      _energy(r._energy)
{
}

Skin& Skin::operator=(const Skin &r)
{
  if (this != &r)
  {
    _skf = r._skf;
    _dbgrid = r._dbgrid;
    _nxyz = r._nxyz;
    _ndim = r._ndim;
    _nval = r._nval;
    _date = r._date;
    _nvalMax = r._nvalMax;
    _total = r._total;
    _totalMax = r._totalMax;
    _address = r._address;
    _energy = r._energy;
  }
  return *this;
}

Skin::~Skin()
{
}

/****************************************************************************/
/*!
 **  Returns the weight for a given cell and direction
 **
 ** \return  The weight
 **
 ** \param[in]  ipos  Absolute grid index of the input grid node
 ** \param[in]  idir  Rank of the direction
 **
 *****************************************************************************/
double Skin::_getWeight(int ipos, int idir)
{
  return _skf->getWeight(ipos, idir);
}

/****************************************************************************/
/*!
 **  Returns the shifted node of a skin
 **
 ** \return  Absolute sample address (or ITEST)
 **
 ** \param[in]  indg0  Array of directional grid indices
 ** \param[in]  dir    Rank of the direction
 **
 *****************************************************************************/
int Skin::_gridShift(const VectorInt& indg0, int dir)
{
  VectorInt indg = indg0;

  /* Shift the target grid node and check if it belongs to the grid */

  for (int i = 0; i < _ndim; i++)
  {
    indg[i] = indg0[i] + id[dir][i];
    if (indg[i] < 0 || indg[i] >= _dbgrid->getNX(i)) return ITEST;
  }
  return _dbgrid->indiceToRank(indg);
}

/****************************************************************************/
/*!
 **  Returns the shifted node of a skin
 **
 ** \return  The absolute sample address
 **
 ** \param[in]  skin  Skin2 structure
 ** \param[in]  lec   Absolute grid index of the input grid node
 ** \param[in]  dir   Rank of the direction
 **
 *****************************************************************************/
int Skin::gridShift(int lec, int dir)
{
  VectorInt indg(_ndim);

  /* Convert an absolute address into the grid indices */

  _dbgrid->rankToIndice(lec, indg);

  /* Shift the target grid node and check if it belongs to the grid */

  return _gridShift(indg, dir);
}

/*****************************************************************************/
/*!
 **  Delete a cell from the skin
 **
 ** \return  Error returned code
 **
 ** \param[in]  rank     Rank of the cell to be deleted
 **
 ** \remark  When deleting a cell from the skin, the last cell is copied
 ** \remark  in the place of the deleted one
 **
 *****************************************************************************/
void Skin::_cellDelete(int rank)
{
  /* Delete the target cell : move the last cell to the target location */

  _nval--;
  _address[rank] = _address[_nval];
  _energy[rank] = _energy[_nval];

  /* Deallocate complementary room in the skin */

  _address.resize(_nval);
  _energy.resize(_nval);
}

/*****************************************************************************/
/*!
 **  Checks if a cell already belongs to the skin
 **
 ** \return  Rank within the skin or -1 if the cell does not belong to the skin
 **
 ** \param[in]  ipos     Cell location
 **
 *****************************************************************************/
int Skin::_cellAlreadyFilled(int ipos)
{
  for (int i = 0; i < _nval; i++)
    if (_address[i] == ipos) return (i);
  return (-1);
}

/*****************************************************************************/
/*!
 **  Modify the energy for a cell which already belongs to the skin
 **
 ** \param[out] rank     Location of the cell in the skin
 ** \param[in]  energy   Additional energy for the new cell
 **
 *****************************************************************************/
void Skin::_cellModify(int rank, double energy)
{
  _energy[rank] += energy;
  _total += energy;
  if (_total > _totalMax) _totalMax = _total;
}

/*****************************************************************************/
/*!
 **  Add a cell to the skin (if not already in the skin)
 **
 ** \return  Error returned code
 **
 ** \param[in]  ipos     Cell location
 ** \param[in]  energy   Energy for the new cell
 **
 *****************************************************************************/
int Skin::_cellAdd(int ipos, double energy)
{
  int rank = _nval;
  _address.resize(_nval + 1);
  _energy.resize(_nval + 1);
  _address[rank] = ipos;
  _energy[rank] = 0.;
  _nval++;
  if (_nval > _nvalMax) _nvalMax = _nval;

  /* Upgrade the energy */

  _cellModify(rank, energy);

  return 0;
}

/*****************************************************************************/
/*!
 **  Initialize the skin
 **
 ** \return  Error returned code
 **
 ** \param[in] verbose  Verbose flag
 **
 *****************************************************************************/
int Skin::init(bool verbose)
{
  if (_skf == nullptr || _ndim <= 0)
  {
    messerr("SKF and DbGrid must be defined beforehand");
    return 1;
  }
  VectorInt indg(_ndim);
  int nb_mask = 0;
  int nb_count = 0;
  int nb_done = 0;
  int total = _nxyz;

  // Loop on all the cells

  for (int lec = 0; lec < total; lec++)
  {
    if (_skf->isAlreadyFilled(lec))
    {

      /* The cell does not belong to the skin: it is already filled */

      nb_done++;
      continue;
    }
    else if (! _skf->isToBeFilled(lec))
    {

      /* The cell does not belong to the cell: it is masked off */

      nb_mask++;
      continue;
    }
    else
    {

      /* The cell is eligible */

      nb_count++;
      int local = 0.;
      _dbgrid->rankToIndice(lec, indg);
      for (int dir = 0; dir < ndir[_ndim]; dir++)
      {
        int ecr = _gridShift(indg, dir);
        if (IFFFF(ecr)) continue;
        if (! _skf->isAlreadyFilled(ecr)) continue;
        local += _skf->getWeight(ecr, invdir[dir]);
      }
      if (local > 0.)
      {
        if (_cellAdd(lec, local))
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

  if (nb_count <= 0 || _total <= 0.)
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
 *****************************************************************************/
int Skin::remains(bool verbose)
{
  _date++;
  if (verbose)
    message("Skin iteration:%5d - Length:%4d - Energy:%lf\n",
            _date, _nval, _total);
  return ((int) _total);
}

/*****************************************************************************/
/*!
 **  Find the next cell at random within the skin
 **
 ** \param[out] rank     Location of the cell in the skin
 ** \param[out] ipos     Cell location
 **
 *****************************************************************************/
void Skin::next(int *rank, int *ipos)
{
  /* Draw a random cell */

  double tirage = _total * law_uniform(0., 1.);

  /* Find the cell */

  double total = 0.;
  for (int i = 0; i < _nval; i++)
  {
    total += _energy[i];
    if (total >= tirage)
    {
      *rank = i;
      *ipos = _address[i];
      if (! _skf->isToBeFilled(*ipos))
        messageAbort(
            "Elligible cell (%d ipos=%d) of the skin is already filled", i,
            *ipos);
      return;
    }
  }
  messageAbort("Cannot find a cell for propagation");
}

/*****************************************************************************/
/*!
 **  Suppress the current cell from the skin
 **
 ** \return  Error return code
 **
 ** \param[in] rank0    Rank of the current cell in the skin
 ** \param[in] ipos0    Cell location
 **
 *****************************************************************************/
int Skin::unstack(int rank0, int ipos0)
{
  VectorInt indg(_ndim);

  /* Suppress the current cell from the skin */

  _total -= _energy[rank0];
  _cellDelete(rank0);

  /* Update the neighboring cells */

  int local = 0.;
  _dbgrid->rankToIndice(ipos0, indg);
  for (int dir = 0; dir < ndir[_ndim]; dir++)
  {
    int ecr = _gridShift(indg, dir);
    if (IFFFF(ecr)) continue;

    /* Discard the neighboring cell if it cannot filled */

    if (! _skf->isToBeFilled(ecr)) continue;
    local = _skf->getWeight(ipos0, dir);
    int rank = _cellAlreadyFilled(ecr);
    if (rank < 0)
    {

      /* The cell does not already belong to the skin: add it */

      if (_cellAdd(ecr, local)) return (1);
    }
    else
    {
      /* If the cell already belongs to the skin, upgrade its energy */

      _cellModify(rank, local);
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
void Skin::skinPrint()
{
  mestitle(1, "Skin algorithm: Final status");
  message("- Number of iterations          = %d\n", _date);
  message("- Maximum skin length           = %d\n", _nvalMax);
  message("- Maximum energy                = %lf\n", _totalMax);
  return;
}
