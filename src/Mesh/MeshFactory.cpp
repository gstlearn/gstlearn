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
#include "Mesh/AMesh.hpp"
#include "Mesh/MeshEStandard.hpp"
#include "Mesh/MeshETurbo.hpp"
#include "Mesh/MeshFactory.hpp"
#include "Mesh/MeshSpherical.hpp"
#include "Db/Db.hpp"
#include "geoslib_f.h"

MeshFactory::MeshFactory()
{
}

MeshFactory::~MeshFactory()
{
}

AMesh* MeshFactory::createMesh(int variety,
                               const VectorDouble& extendmin,
                               const VectorDouble& extendmax,
                               const VectorDouble& cellsize,
                               const VectorDouble& rotmat,
                               const VectorDouble& dilate,
                               Db *dbin,
                               Db *dbout,
                               const String& triswitch,
                               MatrixRectangular& apices,
                               VectorInt& meshes,
                               bool flag_polarize,
                               int verbose)
{
  int ndim = 0;

  // Creating the relevant class element

  if (variety == 0)
  {
    // Euclidean case

    if (dbin != nullptr)
    {
      if (verbose)
      {
        mestitle(0,"Standard Meshing (from Dbs)");
      }
      if (dbin != nullptr)
      {
        if (ndim > 0 && dbin->getNDim() != ndim)
        {
          messerr("Dbin's space dimension (%d) should match 'ndim' (%d)",
                  dbin->getNDim(), ndim);
          return (NULL);
        }
        ndim = dbin->getNDim();
      }
      if (dbout != nullptr)
      {
        if (ndim > 0 && dbout->getNDim() != ndim)
        {
          messerr("Dbout's space dimension (%d) should match 'ndim' (%d)",
                  dbout->getNDim(), ndim);
          return (NULL);
        }
        ndim = dbout->getNDim();
      }
      MeshEStandard* mesh = new MeshEStandard();
      if (!mesh->resetFromDb(dbin, dbout, dilate, triswitch, verbose))
        return ((AMesh *) mesh);
    }
    else if (! extendmin.empty() &&
             ! extendmax.empty() &&
             ! cellsize.empty())
    {
      if (verbose)
      {
        mestitle(0,"Turbo Meshing");
      }

      // Turbo case
      MeshETurbo* mesh = new MeshETurbo();
      if (!mesh->initFromExtend(extendmin, extendmax, cellsize, rotmat,
                                flag_polarize, verbose))
        return ((AMesh *) mesh);
    }
    else if (! apices.isEmpty() && ! meshes.empty())
    {
      if (verbose)
      {
        mestitle(0,"Standard Meshing (from External Mesh)");
      }

      MeshEStandard* mesh = new MeshEStandard();
      if (! mesh->reset(apices, meshes, verbose))
        return ((AMesh *) mesh);
    }
  }
  else if (variety == 1)
  {
    if (verbose)
    {
      mestitle(0,"Spherical Meshing");
    }
    MeshSpherical* mesh = new MeshSpherical();
    if (!mesh->reset(dbin, dbout, triswitch, verbose))
      return ((AMesh *) mesh);
  }
  else
  {
    messerr("No relevant option found in Mesh Factory");
    return (NULL);
  }
  return (NULL);
}
