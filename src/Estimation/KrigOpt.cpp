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
#include "Estimation/KrigOpt.hpp"
#include "Enum/EKrigOpt.hpp"
#include "Matrix/MatrixRectangular.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/VectorHelper.hpp"
#include "Db/DbGrid.hpp"

KrigOpt::KrigOpt(const EKrigOpt& calcul)
  : _calcul(calcul)
  , _mode()
  , _flagPerCell(false)
  , _ndim(0)
  , _ndiscNumber(0)
  , _ndiscs()
  , _disc1()
  , _disc2()
  , _flagDGM(false)
  , _matLC()
  , _dbgrid()
{
  _mode = CovCalcMode(ECalcMember::RHS);
}

KrigOpt::KrigOpt(const KrigOpt& m)
  : _calcul(m._calcul)
  , _mode(m._mode)
  , _flagPerCell(m._flagPerCell)
  , _ndim(m._ndim)
  , _ndiscNumber(m._ndiscNumber)
  , _ndiscs(m._ndiscs)
  , _disc1(m._disc1)
  , _disc2(m._disc2)
  , _matLC(m._matLC)
  , _dbgrid(m._dbgrid)
{
}

KrigOpt& KrigOpt::operator=(const KrigOpt& m)
{
  if (this != &m)
  {
    _calcul      = m._calcul;
    _mode        = m._mode;
    _flagPerCell = m._flagPerCell;
    _ndim        = m._ndim;
    _ndiscNumber = m._ndiscNumber;
    _ndiscs      = m._ndiscs;
    _disc1       = m._disc1;
    _disc2       = m._disc2;
    _matLC       = m._matLC;
    _dbgrid      = m._dbgrid;
  }
  return *this;
}

KrigOpt::~KrigOpt()
{
}

int KrigOpt::getNvarCL() const
{
  if (_matLC == nullptr) return 0;
  return _matLC->getNRows();
}

double KrigOpt::getMatCLValue(int ivarcl, int ivar) const
{
  if (_matLC == nullptr) return TEST;
  return _matLC->getValue(ivarcl, ivar);
}

/**
 * Define the output as Linear Combinations of the Input Variables
 * @param matLC Vector of Vectors of weights (see remarks)
 * @param nvar  Number of Input variables
 * @return
 * @remarks The number of Rows of d'matLC' is the number of Output variables
 * @remarks The number of Columns of 'matLC' is the number of input Variables.
 */
int KrigOpt::setMatLC(const MatrixRectangular* matLC, int nvar)
{
  if (matLC == nullptr) return 0;
  int n1 = (int)matLC->getNRows();
  int n2 = (int)matLC->getNCols();

  if (n1 > nvar)
  {
    messerr("Number of Rows of 'matLC' (%d)", (int)n1);
    messerr("should be smaller than the number of variables (%d)", nvar);
    return 1;
  }
  if (n2 != nvar)
  {
    messerr("Number of Columns of 'matLC' (%d)", (int)n2);
    messerr("should be equal to the number of variables (%d)", nvar);
    return 1;
  }
  _matLC = matLC;
  return 0;
}

int KrigOpt::setKrigingOption(const EKrigOpt& calcul,
                              DbGrid* dbgrid,
                              const VectorInt& ndiscs,
                              bool flag_per_cell)
{
  _calcul = calcul;

  // Clear all parameters
  _dbgrid = nullptr;
  _ndiscNumber = 0;
  _ndiscs.clear();
  _disc1.clear();
  _disc2.clear();

  // Particular options for Block Kriging
  if (calcul == EKrigOpt::BLOCK)
  {
    // Preliminary checks
    if (ndiscs.empty())
    {
      messerr("In case of BLOCK kriging, you must define the discretization parameters");
      messerr("i.e. a vector (dimension: Space Dimension) filled with positive numbers");
      return 1;
    }
    if (dbgrid == nullptr)
    {
      messerr("For Block Kriging, the output must be a DbGrid");
      return 1;
    }
    _ndiscs = ndiscs;
    _dbgrid = dbgrid;
    _flagPerCell = flag_per_cell;

    // Prepare auxiliary storage
    _ndim        = (int)ndiscs.size();
    _ndiscNumber = VH::product(_ndiscs);
    _disc1.resize(_ndiscNumber);
    _disc2.resize(_ndiscNumber);
    for (int i = 0; i < _ndiscNumber; i++)
    {
      _disc1[i].resize(_ndim);
      _disc2[i].resize(_ndim);
    }

    // For constant discretization, calculate discretization coordinates
    if (!_flagPerCell) blockDiscretize(0, true);
  }

  return 0;
}

int KrigOpt::setKrigingDGM(bool flag_dgm)
{
  _flagDGM = flag_dgm;
  return 0;
}

void KrigOpt::blockDiscretize(int iechout, bool flagRandom, int seed) const
{
  _disc1 = _dbgrid->getDiscretizedBlock(_ndiscs, iechout, _flagPerCell, false);

  if (flagRandom)
    _disc2 = _dbgrid->getDiscretizedBlock(_ndiscs, iechout, _flagPerCell, true, seed);
}

double KrigOpt::_getDisc1(int idisc, int idim) const
{
  return _disc1[idisc][idim];
}
VectorDouble KrigOpt::getDisc1VD(int idisc) const
{
  VectorDouble vec(_ndim);
  for (int idim = 0; idim < _ndim; idim++) vec[idim] = _disc1[idisc][idim];
  return vec;
}
VectorVectorDouble KrigOpt::getDisc1VVD() const
{
  VectorVectorDouble vecvec(_ndiscNumber);
  for (int idisc = 0; idisc < _ndiscNumber; idisc++) vecvec[idisc] = getDisc1VD(idisc);
  return vecvec;
}
double KrigOpt::_getDisc2(int idisc, int idim) const
{
  return _disc2[idisc][idim];
}
VectorDouble KrigOpt::getDisc2VD(int idisc) const
{
  VectorDouble vec(_ndim);
  for (int idim = 0; idim < _ndim; idim++) vec[idim] = _disc2[idisc][idim];
  return vec;
}
VectorVectorDouble KrigOpt::getDisc2VVD() const
{
  VectorVectorDouble vecvec(_ndiscNumber);
  for (int idisc = 0; idisc < _ndiscNumber; idisc++) vecvec[idisc] = getDisc2VD(idisc);
  return vecvec;
}
