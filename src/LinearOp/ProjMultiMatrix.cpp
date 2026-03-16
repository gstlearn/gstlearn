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

#include "LinearOp/ProjMultiMatrix.hpp"
#include "Db/Db.hpp"
#include "LinearOp/ProjMatrix.hpp"
#include "LinearOp/ProjMulti.hpp"
#include "Matrix/MatrixSparse.hpp"
#include "Mesh/VectorMeshes.hpp"

namespace gstlrn
{

template<typename T>
static std::vector<std::vector<const IProj*>> castToBase(const std::vector<std::vector<T>>& vec)
{
  std::vector<std::vector<const IProj*>> casted(vec.size());
  Id iv = 0;
  for (const auto& e: vec)
  {
    auto& dst = casted[iv++];
    dst.resize(e.size());
    Id ie = 0;
    for (const auto& f: e)
    {
      if (!f)
      {
        dst[ie++] = {};
        continue;
      }
      if constexpr (std::is_same<T, const ProjMatrix*>::value)
        dst[ie++] = f;
      else if constexpr (std::is_same<T, std::optional<ProjMatrix>>::value)
        dst[ie++] = &f.value();
    }
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
                                                        const VectorMeshes& meshes,
                                                        Id ncov,
                                                        Id nvar,
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
  Id nmeshes = meshes.size();
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

  if (!meshes.allDefined())
  {
    messerr("All the meshes have to be defined");
    return nullptr;
  }

  ProjsStore stocker(nvar);

  const auto nmesh = meshes.size();
  bool flagIsVar   = checkOnZVariable && db->hasLocator(ELoc::Z);
  for (Id ivar = 0; ivar < nvar; ivar++)
  {
    for (Id imesh = 0; imesh < nmesh; imesh++)
      for (Id jvar = 0; jvar < nvar; jvar++)
      {
        Id kvar = (flagIsVar) ? jvar : -1;
        if (ivar != jvar)
          stocker[ivar].emplace_back();
        else
          stocker[ivar].emplace_back(ProjMatrix(db, meshes(imesh), kvar, verbose));
      }
  }
  return new ProjMultiMatrix(std::move(stocker));
}

std::vector<std::vector<const ProjMatrix*>> ProjMultiMatrix::create(std::vector<const ProjMatrix*>& vectproj,
                                                                    Id nvariable)
{
  Id nlatent = static_cast<Id>(vectproj.size());
  std::vector<std::vector<const ProjMatrix*>> result;

  for (Id i = 0; i < nlatent; i++)
  {
    if (vectproj[i] == nullptr)
    {
      messerr("ProjMatrix shouldn't be nullptr.");
      return result;
    }
  }

  Id npoint = vectproj[0]->getNPoint();

  for (Id i = 1; i < nlatent; i++)
  {
    if (vectproj[i]->getNPoint() != npoint)
    {
      messerr("All the ProjMatrix should have the same number of Point.");
      messerr("Element %d has %d Point instead of %d.", i, vectproj[i]->getNPoint(), npoint);
      return result;
    }
  }

  result.resize(nvariable);
  for (Id i = 0; i < nvariable; i++)
  {
    std::vector<const ProjMatrix*> e(nlatent * nvariable, nullptr);
    for (Id j = 0; j < nlatent; j++)
    {
      e[(j * nvariable) + i] = vectproj[j];
    }
    result[i] = e;
  }
  return result;
}

ProjMultiMatrix::ProjMultiMatrix(const std::vector<std::vector<const ProjMatrix*>>& proj,
                                 bool silent)
  : ProjMulti(castToBase(proj), silent)
  , _Proj(MatrixSparse(0, 0))
{
  this->init();
}

ProjMultiMatrix::ProjMultiMatrix(ProjsStore&& proj,
                                 bool silent)
  : ProjMulti(castToBase(proj), silent)
  , _projsStore(std::move(proj))
  , _Proj(MatrixSparse(0, 0))
{
  this->init();
}

void ProjMultiMatrix::init()
{
  if (ProjMulti::empty()) return;
  const VectorInt& pointNumbers = getNPoints();
  const VectorInt& apexNumbers  = getNApexs();

  for (Id i = 0; i < getNVariable(); i++)
  {
    MatrixSparse currentrow;
    for (Id j = 0; j < getNLatent(); j++)
    {
      if (_projs[i][j])
      {
        MatrixSparse::glueInPlace(&currentrow, &static_cast<const ProjMatrix&>(_projs[i][j].value().get()), 0, 1);
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

Id ProjMultiMatrix::_addPoint2mesh(const constvect inv, vect outv) const
{
  AMatrix::prodMVInPlace(outv, _Proj, inv, true, false);
  // _Proj.addProdMatVecInPlaceC(inv, outv, true);
  return 0;
}
Id ProjMultiMatrix::_addMesh2point(const constvect inv, vect outv) const
{
  AMatrix::prodMVInPlace(outv, _Proj, inv, false, false);
  // _Proj.addProdMatVecInPlaceC(inv, outv, false);
  return 0;
}
} // namespace gstlrn
