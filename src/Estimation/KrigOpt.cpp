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
#include "Neigh/ANeigh.hpp"
#include "Model/ModelGeneric.hpp"
#include "Model/Model.hpp"
#include "Neigh/ANeigh.hpp"
#include "Model/ModelGeneric.hpp"
#include "Model/Model.hpp"

KrigOpt::KrigOpt(const EKrigOpt& calcul)
  : _calcul(calcul)
  , _mode()
  , _flagPerCell(false)
  , _nDiscDim(0)
  , _nDiscNumber(0)
  , _ndiscs()
  , _disc1()
  , _disc2()
  , _flagDGM(false)
  , _flagColcok(false)
  , _rankColcok()
  , _matLC()
  , _dbgrid()
{
  _mode = CovCalcMode(ECalcMember::RHS);
}

KrigOpt::KrigOpt(const KrigOpt& m)
  : _calcul(m._calcul)
  , _mode(m._mode)
  , _flagPerCell(m._flagPerCell)
  , _nDiscDim(m._nDiscDim)
  , _nDiscNumber(m._nDiscNumber)
  , _ndiscs(m._ndiscs)
  , _disc1(m._disc1)
  , _disc2(m._disc2)
  , _flagDGM(m._flagDGM)
  , _flagColcok(m._flagColcok)
  , _rankColcok(m._rankColcok)
  , _matLC(m._matLC)
  , _dbgrid(m._dbgrid)
{
}

KrigOpt& KrigOpt::operator=(const KrigOpt& m)
{
  if (this != &m)
  {
    _calcul       = m._calcul;
    _mode         = m._mode;
    _flagPerCell  = m._flagPerCell;
    _nDiscDim     = m._nDiscDim;
    _nDiscNumber  = m._nDiscNumber;
    _ndiscs       = m._ndiscs;
    _disc1        = m._disc1;
    _disc2        = m._disc2;
    _flagDGM      = m._flagDGM;
    _flagColcok   = m._flagColcok;
    _rankColcok   = m._rankColcok;
    _matLC        = m._matLC;
    _dbgrid       = m._dbgrid;
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
                              Db* dbout,
                              const VectorInt& ndiscs,
                              bool flag_per_cell)
{
  _calcul = calcul;

  // Clear all parameters
  _dbgrid = nullptr;
  _nDiscNumber = 0;
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
    DbGrid* dbgrid = dynamic_cast<DbGrid*>(dbout);
    if (dbgrid == nullptr)
    {
      messerr("For Block Kriging, the output must be a DbGrid");
      return 1;
    }
    _ndiscs = ndiscs;
    _dbgrid = dbgrid;
    _flagPerCell = flag_per_cell;

    // Prepare auxiliary storage
    _nDiscDim        = (int)ndiscs.size();
    _nDiscNumber = VH::product(_ndiscs);
    _disc1.resize(_nDiscNumber);
    _disc2.resize(_nDiscNumber);
    for (int i = 0; i < _nDiscNumber; i++)
    {
      _disc1[i].resize(_nDiscDim);
      _disc2[i].resize(_nDiscDim);
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

int KrigOpt::setRankColCok(const VectorInt& rank_colcok)
{
  _rankColcok = rank_colcok;
  _flagColcok = ! rank_colcok.empty();
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
  VectorDouble vec(_nDiscDim);
  for (int idim = 0; idim < _nDiscDim; idim++) vec[idim] = _disc1[idisc][idim];
  return vec;
}
VectorVectorDouble KrigOpt::getDisc1VVD() const
{
  VectorVectorDouble vecvec(_nDiscNumber);
  for (int idisc = 0; idisc < _nDiscNumber; idisc++) vecvec[idisc] = getDisc1VD(idisc);
  return vecvec;
}
double KrigOpt::_getDisc2(int idisc, int idim) const
{
  return _disc2[idisc][idim];
}
VectorDouble KrigOpt::getDisc2VD(int idisc) const
{
  VectorDouble vec(_nDiscDim);
  for (int idim = 0; idim < _nDiscDim; idim++) vec[idim] = _disc2[idisc][idim];
  return vec;
}
VectorVectorDouble KrigOpt::getDisc2VVD() const
{
  VectorVectorDouble vecvec(_nDiscNumber);
  for (int idisc = 0; idisc < _nDiscNumber; idisc++) vecvec[idisc] = getDisc2VD(idisc);
  return vecvec;
}

void KrigOpt::setMode(const CovCalcMode* mode)
{
  if (mode != nullptr)
    _mode = *mode;
  else
    _mode = CovCalcMode(ECalcMember::RHS);
}

bool KrigOpt::_isValidCalcul(const Db* dbout, const ANeigh* neigh) const
{
  // Check the Block calculation
  if (_calcul == EKrigOpt::BLOCK)
  {
    const DbGrid* dbgrid = dynamic_cast<const DbGrid*>(dbout);
    if (dbgrid == nullptr)
    {
      messerr("Block Estimation is only possible for Grid '_dbout'");
      return false;
    }

    // Block support is defined per sample
    if (neigh->getType() == ENeigh::CELL)
    {
      _flagPerCell = true;
    }

    // Check that discretization is defined
    if (_ndiscs.empty())
    {
      messerr("In case of BLOCK kriging, you must define the discretization coefficients");
      messerr("i.e. a vector (dimension equal Space Dimension) filled with positive numbers");
      return false;
    }
  }
  return true;
}

bool KrigOpt::_isValidColcok(const Db* dbout, const ModelGeneric* model) const
{
  if (!_flagColcok) return true;

  int nvar = model->getNVar();

  /* Loop on the ranks of the colocated variables */

  for (int ivar = 0; ivar < nvar; ivar++)
  {
    int jvar = _rankColcok[ivar];
    if (jvar < 0) continue;
    if (jvar > dbout->getNLoc(ELoc::Z))
    {
      messerr("Error in the Colocation array:");
      messerr("Input variable (#%d): rank of the colocated variable is %d",
              ivar + 1, jvar);
      messerr("But the Output file only contains %d attributes(s)",
              dbout->getNColumn());
      return false;
    }
  }
  return true;
}

bool KrigOpt::_isValidMatLC(const ModelGeneric* model) const
{
  if (_matLC == nullptr) return true;
  if (_matLC->empty()) return true;
  int nvar = model->getNVar();
  int n1   = (int)_matLC->getNRows();
  int n2   = (int)_matLC->getNCols();

  if (n1 > nvar)
  {
    messerr("First dimension of 'matLC' (%d)", (int)n1);
    messerr("should be smaller than the number of variables in the model (%d)", nvar);
    return false;
  }
  if (n2 != nvar)
  {
    messerr("Second dimension of 'matLC' (%d)", (int)n2);
    messerr("should be equal to the number of variables in the model (%d)", nvar);
    return false;
  }
  return true;
}

bool KrigOpt::_isValidDGM(const ModelGeneric* model) const
{
  if (!_flagDGM) return false;
  int nvar                = model->getNVar();
  const Model* modelAniso = dynamic_cast<const Model*>(model);
  if (modelAniso == nullptr)
  {
    messerr("The option DGM is limited to model Aniso");
    return false;
  }
  if (modelAniso->getCovMinIRFOrder() != -1)
  {
    messerr("The option DGM is limited to Stationary Covariances");
    return false;
  }
  if (nvar != 1)
  {
    messerr("The DGM option is limited to the Monovariate case");
    return false;
  }
  if (ABS(modelAniso->getTotalSill(0, 0) - 1.) > 1.e-6)
  {
    messerr("The DGM option requires a Model with Total Sill equal to 1.");
    return false;
  }

  if (_calcul == EKrigOpt::BLOCK || _calcul == EKrigOpt::DRIFT)
  {
    messerr("The DGM option is incompatible with 'Block' calculation option");
    return false;
  }
  return true;
}

  bool KrigOpt::isValid(const Db* dbout, const ANeigh* neigh, const ModelGeneric* model) const
{
  // Check against Block calculation options
  if (! _isValidCalcul(dbout, neigh)) return false;

  // Check against Colocated CoKriging options
  if (!_isValidColcok(dbout, model)) return false;

  // Check against the matLC option
  if (!_isValidMatLC(model)) return false;

   // Check the validity for Discrete Gaussian Model
  if (!_isValidDGM(model)) return false;
 
  return true;
}

void KrigOpt::dumpOptions() const
{
  switch (_calcul.toEnum())
  {
    case EKrigOpt::E_POINT:
    {
      message("Punctual Estimation\n");
      break;
    }

    case EKrigOpt::E_BLOCK:
    {
      message("Block Estimation : Discretization = ");
      for (int idim = 0; idim < _nDiscDim; idim++)
      {
        if (idim != 0) message(" x ");
        message("%d", getDisc(idim));
      }
      message("\n");
      break;
    }

    case EKrigOpt::E_DRIFT:
    {
      message("Drift Estimation\n");
      break;
    }

    case EKrigOpt::E_DGM:
    {
      message("Discrete Gaussian Model\n");
      break;
    }
  }
  message("\n");
}
