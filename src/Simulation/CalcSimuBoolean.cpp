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
#include "Simulation/CalcSimuBoolean.hpp"
#include "Basic/Law.hpp"
#include "Boolean/AShape.hpp"
#include "Db/Db.hpp"
#include "Db/DbGrid.hpp"
#include "Simulation/ACalcSimulation.hpp"
#include "Simulation/BooleanObject.hpp"
#include "Simulation/SimuBooleanParam.hpp"

#include <cmath>

namespace gstlrn
{
  CalcSimuBoolean::CalcSimuBoolean(Id nbsimu, Id seed, bool verbose)
    : ACalcSimulation(nbsimu, seed, verbose)
    , AStringable()
    , _flagSimu(true)
    , _flagRank(false)
    , _boolparam()
    , _tokens(nullptr)
    , _objlist()
    , _iptrSimu(-1)
    , _iptrRank(-1)
    , _iptrCover(-1)
  {
  }

  CalcSimuBoolean::~CalcSimuBoolean()
  {
    _clearAllObjects();
  }

  void CalcSimuBoolean::_clearAllObjects()
  {
    if (_objlist.empty()) return;
    for (Id iobj = 0; iobj < _getNObjects(); iobj++) delete _objlist[iobj];
  }

  String CalcSimuBoolean::toString(const AStringFormat* strfmt) const
  {
    std::stringstream sstr;

    for (Id iobj = 0; iobj < _getNObjects(); iobj++)
    {
      sstr << "Characteristics of the Object: " << iobj + 1 << std::endl;
      sstr << _objlist[iobj]->toString(strfmt);
    }
    return sstr.str();
  }

  bool CalcSimuBoolean::_simulate()
  {
    /* Define the global variables */
    law_set_random_seed(getSeed());

    /* Count the number of conditioning pores and grains */
    if (getVerbose()) mestitle(0, "Boolean simulation");

    // Clear any existing object
    _clearAllObjects();

    // Simulate the Initial grains (optional if dbin == nullptr)
    if (!_generatePrimary()) return false;

    // Simulate the Secondary grains
    if (!_generateSecondary()) return false;

    // Display all the objects (can be extensive)
    // if (getVerbose()) display();

    // Project the objects on the output grid
    _projectToGrid();

    return true;
  }

  void CalcSimuBoolean::_projectToGrid()
  {
    auto* dbout = dynamic_cast<DbGrid*>(getDbout());
    if (dbout == nullptr) return;
    for (Id iobj = 0, nobj = _getNObjects(); iobj < nobj; iobj++)
    {
      _objlist[iobj]->projectToGrid(
        dbout, _iptrSimu, _iptrRank, static_cast<Id>(_boolparam.getFacies()),
        iobj + 1);
    }
  }

  Id CalcSimuBoolean::_countConditioningPore() const
  {
    auto* db = getDbin();
    if (db == nullptr) return 0;

    Id nbpore = 0;
    for (Id iech = 0, nech = db->getNSample(); iech < nech; iech++)
    {
      if (!db->isActive(iech)) continue;
      if (BooleanObject::isPore(db, iech)) nbpore++;
    }
    return nbpore;
  }

  Id CalcSimuBoolean::_countConditioningGrain() const
  {
    auto* db = getDbin();
    if (db == nullptr) return 0;

    Id nbgrain = 0;
    for (Id iech = 0, nech = db->getNSample(); iech < nech; iech++)
    {
      if (!db->isActive(iech)) continue;
      if (BooleanObject::isGrain(db, iech)) nbgrain++;
    }
    return nbgrain;
  }

  /**
   * @brief Returns the rank of the grain to be covered in the list of non-covered grains
   *
   * @param db Db containing the conditioning information
   * @param rank Rank of the grain to be covered
   * @return Id Index of the uncovered grain (-1 if not found)
   */
  Id CalcSimuBoolean::_getRankGrainUncovered(const Db* db, Id rank) const
  {
    Id number = 0;
    for (Id iech = 0, nech = db->getNSample(); iech < nech; iech++)
    {
      if (!db->isActive(iech)) continue;
      // Skip if it is not a grain
      if (BooleanObject::isPore(db, iech)) continue;

      // Check if the grain is already covered
      if (db->getArray(iech, _iptrCover) > 0) continue;

      // Return the index of the grain if it is the one we are looking for
      if (number == rank) return iech;
      number++;
    }
    messerr("Error when searching for Uncovered Grain #%d", rank);
    return -1;
  }

  Id CalcSimuBoolean::_getNObjects(Id mode) const
  {
    if (mode == 0) return static_cast<Id>(_objlist.size());
    Id number = 0;
    for (Id iobj = 0, nobj = static_cast<Id>(_objlist.size()); iobj < nobj;
         iobj++)
    {
      if (_objlist[iobj]->getMode() == mode) number++;
    }
    return number;
  }

  bool CalcSimuBoolean::_generatePrimary()
  {
    auto* dbin = getDbin();
    if (dbin == nullptr) return true;
    auto* dbout = dynamic_cast<DbGrid*>(getDbout());
    if (dbout == nullptr) return false;
    Id ndim = dbout->getNDim();

    // Count the conditioning information
    Id nbpore = _countConditioningPore();
    Id nbgrain = _countConditioningGrain();
    VectorDouble cdgrain(ndim);

    if (getVerbose())
    {
      message("- Conditioning option               = YES\n");
      mestitle(1, "Simulating the initial tokens");
      message("- Number of grains to be covered    = %d\n", nbgrain);
      message("- Number of conditioning pores      = %d\n", nbpore);
    }

    // Generate the initial objects
    Id draw_more = nbgrain;
    Id iter = 0;
    while (draw_more)
    {
      iter++;
      if (iter >= _boolparam.getMaxiter())
      {
        messerr(
          "Simulation of the initial objects failed after %d iterations to "
          "cover %d of the %d grains",
          _boolparam.getMaxiter(), draw_more, nbgrain);
        messerr("Check the Token definition or the Intensity value(s)");
        return false;
      }

      // Look for a non-covered grain
      Id rank = static_cast<Id>(draw_more * law_uniform(0., 1.));
      auto iech = _getRankGrainUncovered(dbin, rank);
      if (iech < 0) return false;
      dbin->getCoordinatesInPlace(cdgrain, iech);

      // Generate an object covering the grain(x,y,z)
      BooleanObject* object =
        BooleanObject::generate(dbout, cdgrain, _tokens, _boolparam);
      if (object == nullptr) continue;

      // Check if the object is compatible with the constraining pores
      if (!object->isCompatiblePore(dbin)) continue;

      // Check if the object is compatible with the constraining grains
      if (!object->isCompatibleGrainAdd(dbin)) continue;

      /* Add the object to the list */
      object->setMode(1);
      _objlist.push_back(object);

      // Update the coverage
      draw_more = object->coverageUpdate(dbin, _iptrCover, 1);
    }

    if (getVerbose())
    {
      message("- Number of Initial Objects = %d\n", _getNObjects(1));
      message("- Number of iterations      = %d\n", iter);
    }

    return true;
  }

  bool CalcSimuBoolean::_generateSecondary()
  {
    auto* dbin = getDbin();
    auto* dbout = dynamic_cast<DbGrid*>(getDbout());
    if (dbout == nullptr) return false;
    Id iter = 0;
    double tabtime = 0.;
    Id nb_average = _getAverageCount();

    if (getVerbose())
    {
      mestitle(1, "Simulating the secondary tokens");
    }

    while (tabtime < _boolparam.getTmax())
    {
      iter++;
      if (iter >= _boolparam.getMaxiter()) break;

      auto nbObject = _getNObjects();
      tabtime += law_exponential();

      double ratio = static_cast<double>(nb_average)
                   / static_cast<double>(nb_average + nbObject);

      if (law_uniform(0., 1.) <= ratio)
      {
        // Add an object
        BooleanObject* object =
          BooleanObject::generate(dbout, VectorDouble(), _tokens, _boolparam);
        if (object == nullptr) continue;

        // Check if the object is compatible with the constraining pores
        if (!object->isCompatiblePore(dbin)) continue;

        // Check if the object is compatible with the constraining grains
        if (!object->isCompatibleGrainAdd(dbin)) continue;

        // Add the object to the list
        object->setMode(2);
        _objlist.push_back(object);

        // Update the coverage
        object->coverageUpdate(dbin, _iptrCover, 1);
      }
      else
      {

        // Delete a primary object
        if (_deleteObject(1, dbin))
        {

          // Delete a secondary object
          if (_deleteObject(2, dbin)) continue;
        }
      }
    }

    if (getVerbose())
    {
      if (dbin != nullptr)
        message("- Ending number of primary objects  = %d\n", _getNObjects(1));
      message("- Total number of objects           = %d\n", _getNObjects());
    }

    return true;
  }

  Id CalcSimuBoolean::_getObjectRank(Id mode, Id rank)
  {
    Id nb_objects = _getNObjects();

    Id number = 0;
    for (Id iobj = 0, nobj = nb_objects; iobj < nobj; iobj++)
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
   ** \return  Error return code: False if the object cannot be deleted; True otherwise
   **
   ** \param[in]  mode     0 for primary tokens; 1 for secondary tokens
   ** \param[in]  dbin     Db structure
   **
   *****************************************************************************/
  bool CalcSimuBoolean::_deleteObject(Id mode, Db* dbin)
  {
    /* Search for the object to be deleted */

    auto count = _getNObjects(mode);
    if (count <= 0) return false;
    Id rank = static_cast<Id>(count * law_uniform(0., 1.));
    auto iref = _getObjectRank(mode, rank);
    if (iref < 0) return false;

    /* Check if the object can be deleted */

    BooleanObject* object = _objlist[iref];
    if (!object->isCompatibleGrainDelete(dbin, _iptrCover)) return false;

    /* Erase the object from the list*/
    _objlist.erase(_objlist.begin() + iref);

    // Update the coverage
    object->coverageUpdate(dbin, _iptrCover, -1);

    // Delete the object
    delete object;

    return true;
  }

  Id CalcSimuBoolean::_getAverageCount() const
  {
    auto* dbout = getDbout();
    double theta;
    if (_tokens->isFlagStat())
    {
      theta = _tokens->getThetaCst();
    }
    else
    {
      VectorDouble vec = dbout->getColumnByLocator(ELoc::P, 0, true);
      theta = vec.mean();
    }

    VectorDouble field = dbout->getExtends();

    // Dilate the field (optional)

    Id ndim = dbout->getNDim();
    double volume = 1.;
    for (Id idim = 0; idim < ndim; idim++)
    {
      field[idim] += 2 * _boolparam.getDilate(idim);
      volume *= field[idim];
    }
    return static_cast<Id>(theta * volume);
  }

  VectorDouble CalcSimuBoolean::extractObjects() const
  {
    VectorDouble tabs;

    for (Id iobj = 0, nobj = _getNObjects(); iobj < nobj; iobj++)
    {
      VectorDouble tab = _objlist[iobj]->getValues();
      tabs.insert(tab.end(), tab.begin(), tab.end());
    }
    return tabs;
  }

  bool CalcSimuBoolean::_check()
  {
    if (!ACalcSimulation::_check()) return false;

    if (!hasDbout()) return false;
    if (!getDbout()->isGrid())
    {
      messerr("The argument 'dbout' should be a grid");
      return false;
    }
    auto ndim = _getNDim();
    if (ndim > 3)
    {
      messerr("The Boolean Method is not a relevant simulation model");
      messerr("for this Space Dimension (%d)", ndim);
      return false;
    }
    if (getDbin() != nullptr)
    {
      if (getDbin()->getNLoc(ELoc::Z) != 1)
      {
        messerr("Conditional Boolean simulation needs 1 variable");
        return false;
      }
    }

    return true;
  }

  bool CalcSimuBoolean::_preprocess()
  {
    if (!ACalcSimulation::_preprocess()) return false;

    _iptrCover = -1;
    if (getDbin() != nullptr)
    {
      _iptrCover = _addVariableDb(1, 1, ELoc::SIMU, 0, getNbSimu(), 0.);
      if (_iptrCover < 0) return false;
    }
    _iptrSimu = -1;
    if (_flagSimu)
    {
      _iptrSimu = _addVariableDb(
        2, 1, ELoc::SIMU, 0, getNbSimu(), _boolparam.getBackground());
      if (_iptrSimu < 0) return false;
    }
    _iptrRank = -1;
    if (_flagRank)
    {
      _iptrRank = _addVariableDb(
        2, 1, ELoc::SIMU, 0, getNbSimu(), _boolparam.getBackground());
      if (_iptrRank < 0) return false;
    }
    return true;
  }

  bool CalcSimuBoolean::_run()
  {
    return (_simulate());
  }

  bool CalcSimuBoolean::_postprocess()
  {
    /* Free the temporary variables */
    _cleanVariableDb(2);

    _renameVariable(
      2, VectorString(), ELoc::Z, 1, _iptrSimu, "Facies", getNbSimu());

    _renameVariable(
      2, VectorString(), ELoc::Z, 1, _iptrRank, "Rank", getNbSimu());
    return true;
  }

  void CalcSimuBoolean::_rollback()
  {
    _cleanVariableDb(1);
  }

} // namespace gstlrn
