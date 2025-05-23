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
#include "Db/Db.hpp"
#include "LinearOp/ProjMultiMatrix.hpp"
#include "Basic/AStringable.hpp"
#include "LinearOp/IProj.hpp"
#include "LinearOp/ProjMatrix.hpp"
#include "LinearOp/ProjMulti.hpp"
#include "Matrix/MatrixSparse.hpp"

static std::vector<std::vector<const IProj*>> castToBase(const std::vector<std::vector<const ProjMatrix*>>& vect)
{
  std::vector<std::vector<const IProj*>> casted(vect.size());
  int iv = 0;
  for (const auto& e: vect)
  {
    std::vector<const IProj*> temp(e.size());
    int ie = 0;
    for (const auto& f: e)
    {
      temp[ie++] = static_cast<const IProj*>(f);
    }
    casted[iv++] = temp;
  }
  return casted;
}

/**
 * @brief Construct the Projection Matrix starting from 'db' and 'meshes'
 * 
 * @param db  Target Db structure
 * @param meshes List of target meshes
 * @param ncov Number of covariances (nugget excluded)
 * @param nvar Number of variables (see notes)
 * @param checkOnZVariable Check if a sample should be considered or not
 * @param verbose Verbose flag
 * @return ProjMultiMatrix 
 * @note Argument 'nvar' is provided as it cannot be derived from 'db'
 * (when 'db' refers to the output file for example, where no Z-variable is available)
 * @note When Z-variable is defined, you can still bypass checking the validity of 
 * a sample (its Z-value is not NA) if 'checkOnZVariable' is False. 
 */
ProjMultiMatrix* ProjMultiMatrix::createFromDbAndMeshes(const Db* db,
                                                        const std::vector<const AMesh*>& meshes,
                                                        int ncov,
                                                        int nvar,
                                                        bool checkOnZVariable,
                                                        bool verbose)
{
   if (db == nullptr)
  {
    messerr("db is null");
    return nullptr;
  }
  if (nvar <= 0)
  {
    messerr("nvar should be > 0");
    return nullptr;
  }
  int nmeshes = (int)meshes.size();
  if (nmeshes == 0)
  {
    messerr("You have to provide at least one mesh");
    return nullptr;
  }
  if (nmeshes != 1 && nmeshes != ncov)
  {
    messerr("Inconsistent number of meshes (%d) and structures (%d)",
            nmeshes, ncov);
    return nullptr;
  }

  for (const auto& e: meshes)
  {
    if (e == nullptr)
    {
      messerr("All the meshes have to be defined");
      return nullptr;
    }
  }

  std::vector<std::vector<const ProjMatrix*>> stocker;

  int nmesh = (int)meshes.size();
  bool flagIsVar = checkOnZVariable && db->hasLocator(ELoc::Z);
  for (int ivar = 0; ivar < nvar; ivar++)
  {
    stocker.push_back(std::vector<const ProjMatrix*>());
    for (int imesh = 0; imesh < nmesh; imesh++)
      for (int jvar = 0; jvar < nvar; jvar++)
      {
        int kvar = (flagIsVar) ? jvar : -1;
        if (ivar != jvar)
          stocker[ivar].push_back(nullptr);
        else
          stocker[ivar].push_back(new ProjMatrix(db, meshes[imesh], kvar, verbose));
      }
  }
  return new ProjMultiMatrix(stocker, true);
}

void ProjMultiMatrix::_clear()
{
    if (!_toClean || _projs.empty()) return;

    for (auto &e: _projs)
    {
        for (auto &f : e)
        {
            delete const_cast<IProj*>(f);
            f = nullptr;
        }
        e.clear();
    }
    _projs.clear();
}
ProjMultiMatrix::~ProjMultiMatrix()
{
}

std::vector<std::vector<const ProjMatrix*>> ProjMultiMatrix::create(std::vector<const ProjMatrix*>& vectproj,
                                                                    int nvariable)
{
  int nlatent = (int)vectproj.size();
  std::vector<std::vector<const ProjMatrix*>> result;

  for (int i = 0; i < nlatent; i++)
  {
    if (vectproj[i] == nullptr)
    {
      messerr("Projmatrix shouldn't be nullptr.");
      return result;
    }
  }

  int npoint = vectproj[0]->getNPoint();

  for (int i = 1; i < nlatent; i++)
  {
    if (vectproj[i]->getNPoint() != npoint)
    {
      messerr("All the ProjMatrix should have the same number of Point.");
      messerr("Element %d has %d Point instead of %d.", i, vectproj[i]->getNPoint(), npoint);
      return result;
    }
  }

  result.resize(nvariable);
  for (int i = 0; i < nvariable; i++)
  {
    std::vector<const ProjMatrix*> e(nlatent * nvariable, nullptr);
    for (int j = 0; j < nlatent; j++)
    {
      e[j * nvariable + i] = vectproj[j];
    }
    result[i] = e;
  }
  return result;
}

ProjMultiMatrix::ProjMultiMatrix(const std::vector<std::vector<const ProjMatrix*>>& proj,
                                 bool toClean,
                                 bool silent)
  : ProjMulti(castToBase(proj), silent)
  , _Proj(MatrixSparse(0, 0))
  , _toClean(toClean)
{
  if (ProjMulti::empty()) return;
  const VectorInt& pointNumbers = getNPoints();
  const VectorInt& apexNumbers  = getNApexs();

  for (int i = 0; i < getNVariable(); i++)
  {
    MatrixSparse currentrow;
    for (int j = 0; j < getNLatent(); j++)
    {
      if (_projs[i][j] != nullptr)
      {
        MatrixSparse::glueInPlace(&currentrow, ((MatrixSparse*)proj[i][j]), 0, 1);
      }
      else
      {
        auto tempMat = MatrixSparse(pointNumbers[i], apexNumbers[j]);
        MatrixSparse::glueInPlace(&currentrow, &tempMat, 0, 1);
      }
    }
    MatrixSparse::glueInPlace(&_Proj, &currentrow, 1, 0);
  }
}

int ProjMultiMatrix::_addPoint2mesh(const constvect inv, vect outv) const
{
  _Proj.addProdMatVecInPlaceToDest(inv, outv, true);
  return 0;
}
int ProjMultiMatrix::_addMesh2point(const constvect inv, vect outv) const
{
  _Proj.addProdMatVecInPlaceToDest(inv, outv, false);
  return 0;
}
