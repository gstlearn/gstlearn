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
#include "Boolean/ObjectList.hpp"
#include "Boolean/Object.hpp"
#include "Boolean/AToken.hpp"
#include "Db/DbGrid.hpp"
#include "Db/Db.hpp"
#include "Basic/Law.hpp"

#include <math.h>

ObjectList::ObjectList()
    : AStringable(),
      _objlist()
{
}

ObjectList::ObjectList(const ObjectList &r)
    : AStringable(r),
      _objlist(r._objlist)
{
}

ObjectList& ObjectList::operator=(const ObjectList &r)
{
  if (this != &r)
  {
    AStringable::operator =(r);
    _objlist = r._objlist;
  }
  return *this;
}

ObjectList::~ObjectList()
{
  for (int iobj = 0; iobj < getNObjects(); iobj++)
    delete _objlist[iobj];
}

String ObjectList::toString(const AStringFormat* strfmt) const
{
  std::stringstream sstr;

  for (int iobj = 0; iobj < getNObjects(); iobj++)
  {
    sstr << "Characteristics of the Object: " << iobj+1 << std::endl;
    sstr << _objlist[iobj]->toString(strfmt);
  }
  return sstr.str();
}

/*****************************************************************************/
/*!
 **  Project the objects on the output grid
 **
 ** \param[in]  dbout         Ouput grid file
 ** \param[in]  iptr_simu     Pointer for storing the Boolean simulation (or <0)
 ** \param[in]  iptr_rank     Pointer for storing the object ranks (or <0)
 ** \param[in]  facies        Value of the facies assigned
 **
 *****************************************************************************/
void ObjectList::projectToGrid(DbGrid* dbout,
                               int iptr_simu,
                               int iptr_rank,
                               int facies)
{
  for (int iobj = 0; iobj < getNObjects(); iobj++)
  {
    _objlist[iobj]->projectToGrid(dbout, iptr_simu, iptr_rank, facies, iobj + 1);
  }
  return;
}

void ObjectList::countConditioning(const Db* db,
                                   int *nbgrain_arg,
                                   int *nbpore_arg,
                                   bool verbose)
{
  *nbgrain_arg = 0;
  *nbpore_arg = 0;
  if (db == nullptr) return;

  int nbgrain = 0;
  int nbpore = 0;
  for (int iech = 0; iech < db->getSampleNumber(); iech++)
  {
    if (! db->isActive(iech)) continue;
    double data = db->getVariable(iech, 0);
    if (FFFF(data)) continue;
    if (data)
      nbgrain++;
    else
      nbpore++;
  }
  if (verbose)
    message("Conditioning data: %d grains and %d pores\n",
            nbgrain, nbpore);

  *nbgrain_arg = nbgrain;
  *nbpore_arg = nbpore;
}

int ObjectList::_getRankUncovered(const Db* db, int rank)
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

int ObjectList::getNObjects(int mode) const
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

int ObjectList::generatePrimary(Db* dbin,
                                DbGrid* dbout,
                                const Tokens* tokens,
                                bool flagStat,
                                double thetaCst,
                                const VectorDouble& dilate,
                                int maxiter)
{
  int nbgrain, nbpore;
  int ndim = dbout->getNDim();
  VectorDouble cdgrain(ndim);

  // Count the conditioning information
  countConditioning(dbin, &nbgrain, &nbpore, false);

  // Generate the initial objects
  int draw_more = nbgrain;
  int iter = 0;
  while (draw_more)
  {
    iter++;
    if (iter >= maxiter)
    {
      messerr("Simulation of the initial objects failed after %d iterations",
              maxiter);
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

    Object* object = Object::generate(dbout, cdgrain, tokens,
                                      flagStat, thetaCst,
                                      dilate, maxiter);
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
  return 0;
}

int ObjectList::generateSecondary(Db* dbin,
                                  DbGrid* dbout,
                                  const Tokens* tokens,
                                  bool flagStat,
                                  double thetaCst,
                                  int nb_average,
                                  double tmax,
                                  const VectorDouble& dilate,
                                  int maxiter)
{
  int iter = 0;
  double tabtime = 0.;
  int Nb_object = getNObjects();

  while (tabtime < tmax)
  {
    if (iter > maxiter) break;

    // The next line is not correct but is kept for compatibility.
    // The correct version should be implemented on next case study
    // update.
    tabtime += law_exponential() / (nb_average + Nb_object);
    // This should be the right version
    //    tabtime += law_exponential();

    double ratio = (double) nb_average / (double) (nb_average + Nb_object);
    if (law_uniform(0., 1.) <= ratio)
    {

      /* Add an object */

      Object* object = Object::generate(dbout, VectorDouble(), tokens,
                                         flagStat, thetaCst,
                                         dilate, maxiter);
      if (object == nullptr) continue;

      /* Check if the object is compatible with the constraining pores */

      if (! object->isCompatiblePore(dbin)) continue;

      /* Check if the object is compatible with the constraining grains */

      if (! object->isCompatibleGrainAdd(dbin)) continue;

      /* Add the object to the list */

      _objlist.push_back(object);

      // Update the coverage

      object->coverageUpdate(dbin, 1);
    }
    else
    {

      /* Delete a primary object */

      if (_deleteObject(0, dbin)) continue;

     /* Delete a secondary object */

      if (_deleteObject(1, dbin)) continue;
    }
  }
  return 0;
}

int ObjectList::_getObjectRank(int mode, int rank)
{
  int nb_objects = getNObjects();

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
int ObjectList::_deleteObject(int mode, Db* dbin)
{
  /* Search for the object to be deleted */

  int count = getNObjects(mode);
  int rank = (int) (count * law_uniform(0., 1.));
  int iref = _getObjectRank(mode, rank);
  if (iref < 0) return 1;

  /* Check if the object can be deleted */

  Object* object = _objlist[iref];
  if (! object->isCompatibleGrainDelete(dbin)) return 1;

  /* Erase the object from the list*/

  _objlist.erase(_objlist.begin() + iref);

  // Update the coverage

  object->coverageUpdate(dbin, -1);

  // Delete the object

  delete object;

  return 0;
}
