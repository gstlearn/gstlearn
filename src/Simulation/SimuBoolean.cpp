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
#include "Boolean/AShape.hpp"
#include "Db/DbGrid.hpp"
#include "Db/Db.hpp"
#include "Simulation/ASimulation.hpp"
#include "Simulation/SimuBoolean.hpp"
#include "Simulation/BooleanObject.hpp"
#include "Simulation/SimuBooleanParam.hpp"
#include "Basic/Law.hpp"

#include <math.h>

SimuBoolean::SimuBoolean(int nbsimu, int seed)
    : ASimulation(nbsimu, seed),
      AStringable(),
      _objlist()
{
}

SimuBoolean::SimuBoolean(const SimuBoolean &r)
    : ASimulation(r),
      AStringable(r),
      _objlist(r._objlist)
{
}

SimuBoolean& SimuBoolean::operator=(const SimuBoolean &r)
{
  if (this != &r)
  {
    ASimulation::operator=(r);
    AStringable::operator =(r);
    _objlist = r._objlist;
  }
  return *this;
}

SimuBoolean::~SimuBoolean()
{
  _clearAllObjects();
}

void SimuBoolean::_clearAllObjects()
{
  if (_objlist.empty()) return;
  for (int iobj = 0; iobj < _getNObjects(); iobj++)
    delete _objlist[iobj];
}

String SimuBoolean::toString(const AStringFormat* strfmt) const
{
  std::stringstream sstr;

  for (int iobj = 0; iobj < _getNObjects(); iobj++)
  {
    sstr << "Characteristics of the Object: " << iobj+1 << std::endl;
    sstr << _objlist[iobj]->toString(strfmt);
  }
  return sstr.str();
}

int SimuBoolean::simulate(Db *dbin,
                          DbGrid *dbout,
                          ModelBoolean* tokens,
                          const SimuBooleanParam& boolparam,
                          int iptr_simu,
                          int iptr_rank,
                          bool verbose)
{
  /* Define the global variables */

  law_set_random_seed(getSeed());

  /* Count the number of conditioning pores and grains */

  if (verbose) mestitle(0,"Boolean simulation");

  // Clear any existing object

  _clearAllObjects();

  // Simulate the Initial grains (optional if dbin == nullptr)

  if (_generatePrimary(dbin, dbout, tokens, boolparam, verbose)) return 1;

  // Simulate the Secondary grains

  if (_generateSecondary(dbin, dbout, tokens, boolparam, verbose)) return 1;

  if (verbose) display();

  // Project the objects on the output grid

  _projectToGrid(dbout, boolparam, iptr_simu, iptr_rank);

  return 0;
}

/*****************************************************************************/
/*!
 **  Project the objects on the output grid
 **
 ** \param[in]  dbout         Ouput grid file
 ** \param[in]  boolparam     SimuBooleanParam structure
 ** \param[in]  iptr_simu     Pointer for storing the Boolean simulation (or <0)
 ** \param[in]  iptr_rank     Pointer for storing the object ranks (or <0)
 **
 *****************************************************************************/
void SimuBoolean::_projectToGrid(DbGrid* dbout,
                                 const SimuBooleanParam& boolparam,
                                 int iptr_simu,
                                 int iptr_rank)
{
  for (int iobj = 0; iobj < _getNObjects(); iobj++)
  {
    _objlist[iobj]->projectToGrid(dbout, iptr_simu, iptr_rank,
                                  boolparam.getFacies(), iobj + 1);
  }
  return;
}

int SimuBoolean::_countConditioningPore(const Db* db)
{
  if (db == nullptr) return 0;

  int nbpore = 0;
  for (int iech = 0; iech < db->getSampleNumber(); iech++)
  {
    if (! db->isActive(iech)) continue;
    double data = db->getVariable(iech, 0);
    if (FFFF(data)) continue;
    if (data) continue;
    nbpore++;
  }
  return nbpore;
}

int SimuBoolean::_countConditioningGrain(const Db* db)
{
  if (db == nullptr) return 0 ;

  int nbgrain = 0;
  for (int iech = 0; iech < db->getSampleNumber(); iech++)
  {
    if (! db->isActive(iech)) continue;
    double data = db->getVariable(iech, 0);
    if (FFFF(data)) continue;
    if (! data) continue;
    nbgrain++;
  }
  return nbgrain;
}

int SimuBoolean::_getRankUncovered(const Db* db, int rank)
{
  int number = 0;
  for (int iech = 0; iech < db->getSampleNumber(); iech++)
  {
    if (! db->isActive(iech)) continue;
    int data = (int) db->getVariable(iech, 0);
    if (data <= 0) continue;
    int cover = (int) db->getVariable(iech, 1);
    if (cover > 0) continue;
    if (number == rank) return iech;
    number++;
  }
  messerr("Error when searching for the rank of Uncovered Grain");
  return -1;
}

int SimuBoolean::_getNObjects(int mode) const
{
  if (mode == 0)
    return (int) _objlist.size();
  else
  {
    int number = 0;
    for (int iobj = 0; iobj < (int) _objlist.size(); iobj++)
    {
      if (_objlist[iobj]->getMode() == mode) number++;
    }
    return number;
  }
}

int SimuBoolean::_generatePrimary(Db* dbin,
                                  DbGrid* dbout,
                                  const ModelBoolean* tokens,
                                  const SimuBooleanParam& boolparam,
                                  bool verbose)
{
  if (dbin == nullptr) return 0;
  int ndim = dbout->getNDim();

  // Count the conditioning information
  int nbpore  = _countConditioningPore(dbin);
  int nbgrain = _countConditioningGrain(dbin);
  VectorDouble cdgrain(ndim);

  if (verbose)
  {
    message("- Conditioning option               = YES\n");
    mestitle(1, "Simulating the initial tokens");
    message("- Number of grains to be covered = %d\n", nbgrain);
    message("- Number of conditioning pores      = %d\n", nbpore);
  }

  // Generate the initial objects
  int draw_more = nbgrain;
  int iter = 0;
  while (draw_more)
  {
    iter++;
    if (iter >= boolparam.getMaxiter())
    {
      messerr("Simulation of the initial objects failed after %d iterations",
              boolparam.getMaxiter());
      messerr("to cover %d of the %d grains", draw_more, nbgrain);
      messerr("Check the Token definition or the Intensity value(s)");
      return 1;
    }

    /* Look for a non-covered grain */

    int rank = (int) (draw_more * law_uniform(0., 1.));
    int iref = _getRankUncovered(dbin, rank);
    if (iref < 0) return 1;
    dbin->getCoordinatesInPlace(iref, cdgrain);

    /* Generate an object covering the grain(x,y,z) */

    BooleanObject* object = BooleanObject::generate(dbout, cdgrain, tokens, boolparam);
    if (object == nullptr) continue;

    /* Check if the object is compatible with the constraining pores */

    if (! object->isCompatiblePore(dbin)) continue;

    /* Check if object is compatible with the constraining grains */

    if (! object->isCompatibleGrainAdd(dbin)) continue;

    /* Add the object to the list */

    object->setMode(1);
    _objlist.push_back(object);

    // Update the coverage

    object->coverageUpdate(dbin, 1);
  }

  if (verbose)
    message("- Number of Initial Objects = %d\n",_getNObjects(1));

  return 0;
}

int SimuBoolean::_generateSecondary(Db* dbin,
                                    DbGrid* dbout,
                                    const ModelBoolean* tokens,
                                    const SimuBooleanParam& boolparam,
                                    bool verbose)
{
  int iter = 0;
  double tabtime = 0.;
  int nb_average = _getAverageCount(dbout, tokens, boolparam);

  if (verbose)
  {
    mestitle(1, "Simulating the secondary tokens");
  }

  while (tabtime < boolparam.getTmax())
  {
    iter++;
    if (iter >= boolparam.getMaxiter()) break;

    int nbObject = _getNObjects();
    // The next line is not correct but is kept for compatibility.
    // The correct version should be implemented on next case study
    // update.
    tabtime += law_exponential() / (nb_average + nbObject);
    // This should be the right version
    //    tabtime += law_exponential();

    double ratio = (double) nb_average / (double) (nb_average + nbObject);

    if (law_uniform(0., 1.) <= ratio)
    {

      /* Add an object */

      BooleanObject* object = BooleanObject::generate(dbout, VectorDouble(), tokens, boolparam);
      if (object == nullptr) continue;

      /* Check if the object is compatible with the constraining pores */

      if (! object->isCompatiblePore(dbin)) continue;

      /* Check if the object is compatible with the constraining grains */

      if (! object->isCompatibleGrainAdd(dbin)) continue;

      /* Add the object to the list */

      object->setMode(2);
      _objlist.push_back(object);

      // Update the coverage

      object->coverageUpdate(dbin, 1);
    }
    else
    {

      /* Delete a primary object */

      if (_deleteObject(1, dbin))
      {

        /* Delete a secondary object */

        if (_deleteObject(2, dbin)) continue;
      }
    }
  }

  if (verbose)
  {
    if (dbin != nullptr)
      message("- Ending number of primary objects  = %d\n",_getNObjects(1));
    message("- Total number of objects           = %d\n",_getNObjects());
  }

  return 0;
}

int SimuBoolean::_getObjectRank(int mode, int rank)
{
  int nb_objects = _getNObjects();

  int number = 0;
  for (int iobj = 0; iobj < nb_objects; iobj++)
  {
    if (_objlist[iobj]->getMode() != mode) continue;
    if (number == rank) return iobj;
    number++;
  }
  return -1;
}

/*****************************************************************************/
/*!
 **  Attempts to delete an object (primary or secondary)
 **
 ** \return  Error return code: 1 if the object cannot be deleted; 0 otherwise
 **
 ** \param[in]  mode     0 for primary tokens; 1 for secondary tokens
 ** \param[in]  dbin     Db structure
 **
 *****************************************************************************/
int SimuBoolean::_deleteObject(int mode, Db* dbin)
{
  /* Search for the object to be deleted */

  int count = _getNObjects(mode);
  if (count <= 0) return 1;
  int rank = (int) (count * law_uniform(0., 1.));
  int iref = _getObjectRank(mode, rank);
  if (iref < 0) return 1;

  /* Check if the object can be deleted */

  BooleanObject* object = _objlist[iref];
  if (! object->isCompatibleGrainDelete(dbin)) return 1;

  /* Erase the object from the list*/

  _objlist.erase(_objlist.begin() + iref);

  // Update the coverage

  object->coverageUpdate(dbin, -1);

  // Delete the object

  delete object;

  return 0;
}

int SimuBoolean::_getAverageCount(const DbGrid* dbout,
                                 const ModelBoolean* tokens,
                                 const SimuBooleanParam& boolparam) const
 {
   double theta;
   if (tokens->isFlagStat())
   {
     theta = tokens->getThetaCst();
   }
   else
   {
     VectorDouble vec = dbout->getColumnByLocator(ELoc::P, 0, true);
     theta = ut_vector_mean(vec);
   }

   VectorDouble field = dbout->getExtends();

   // Dilate the field (optional)

   int ndim = dbout->getNDim();
   double volume = 1.;
   for (int idim = 0; idim < ndim; idim++)
   {
     field[idim] += 2 * boolparam.getDilate(idim);
     volume *= field[idim];
   }
   return (int) (theta * volume);
 }

VectorDouble SimuBoolean::extractObjects() const
{
  VectorDouble tabs;

  for (int iobj = 0; iobj < _getNObjects(); iobj++)
  {
    VectorDouble tab = _objlist[iobj]->getValues();
    tabs.insert(tab.end(), tab.begin(), tab.end());
  }
  return tabs;
}
