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
#include "Db/DbLine.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/Law.hpp"
#include "Basic/SerializeHDF5.hpp"
#include "Basic/String.hpp"
#include "Basic/VectorNumT.hpp"
#include "Db/Db.hpp"
#include "Db/DbGrid.hpp"
#include "Db/DbStringFormat.hpp"
#include "Enum/ELoadBy.hpp"
#include "Enum/ELoc.hpp"
#include "Polygon/Polygons.hpp"
#include "Space/SpacePoint.hpp"
#include "Stats/Classical.hpp"
#include "geoslib_define.h"

#include <cmath>

namespace gstlrn
{
  DbLine::DbLine()
    : Db()
    , _lineAdds()
  {
    _clear();
  }

  DbLine::DbLine(const DbLine& r)
    : Db(r)
    , _lineAdds(r._lineAdds)
  {
  }

  DbLine& DbLine::operator=(const DbLine& r)
  {
    if (this != &r)
    {
      Db::operator=(r);
      _lineAdds = r._lineAdds;
    }
    return *this;
  }

  DbLine::~DbLine()
  {
  }

  /**
   * @brief Check if the target Line number 'iline' (0-based) is valid or not
   *
   * @param iline Target Line number
   * @return true If the Line rank is valid
   * @return false otherwise
   */
  bool DbLine::_isLineNumberValid(Id iline) const
  {
    if (iline < 0)
    {
      messerr("Argument 'iline' should be non negative");
      return false;
    }
    if (iline >= getNLine())
    {
      messerr("ilin' (%d) should be smaller than Number of Lines (%d)", iline,
              getNLine());
      return false;
    }
    return true;
  }

  Id DbLine::getNLine() const
  {
    if (_lineAdds.empty()) return 0;
    return static_cast<Id>(_lineAdds.size());
  }

  Id DbLine::getNSamplePerLine(Id iline) const
  {
    if (!_isLineNumberValid(iline)) return -1;
    return static_cast<Id>(_lineAdds[iline].size());
  }

  Id DbLine::getNTotal() const
  {
    Id ntotal = 0;
    for (Id iline = 0, nbline = getNLine(); iline < nbline; iline++)
      ntotal += getNSamplePerLine(iline);
    return ntotal;
  }

  double DbLine::getLineLength(Id iline) const
  {
    if (!_isLineNumberValid(iline)) return TEST;
    double total = 0.;
    SpacePoint P1;
    SpacePoint P2;
    getSampleAsSPInPlace(P1, _lineAdds[iline][0]);
    for (Id iech = 1, nech = getNSamplePerLine(iline); iech < nech; iech++)
    {
      getSampleAsSPInPlace(P2, _lineAdds[iline][iech]);
      total += P2.getDistance(P1);
      P1 = P2;
    }
    return total;
  }

  VectorDouble DbLine::getLineLengths() const
  {
    auto nline = getNLine();
    VectorDouble lengths(nline);
    for (Id iline = 0; iline < nline; iline++)
      lengths[iline] = getLineLength(iline);
    return lengths;
  }

  String DbLine::toString(const AStringFormat* strfmt) const
  {
    std::stringstream sstr;

    const auto* dbfmt = dynamic_cast<const DbStringFormat*>(strfmt);
    DbStringFormat dsf;
    if (dbfmt != nullptr) dsf = *dbfmt;

    sstr << toTitle(0, "Data Base Line Characteristics");

    sstr << "Number of Lines = " << getNLine() << std::endl;
    sstr << "Number of samples = " << getNSample() << std::endl;
    sstr << "Line length = ";
    for (Id iline = 0, nbline = getNLine(); iline < nbline; iline++)
    {
      if (iline > 0) sstr << " / ";
      sstr << getNSamplePerLine(iline);
    }
    sstr << std::endl;

    sstr << _toStringCommon(&dsf);

    return sstr.str();
  }

  DbLine* DbLine::createFromSamples(Id nech,
                                    const ELoadBy& order,
                                    const VectorDouble& tab,
                                    const VectorInt& lineCounts,
                                    const VectorString& names,
                                    const VectorString& locatorNames,
                                    bool flagAddSampleRank)
  {
    DbLine* dbline = new DbLine;
    if (dbline->resetFromSamples(nech, order, tab, lineCounts, names, locatorNames,
                                 flagAddSampleRank))
    {
      messerr("Error when creating DbLine from Samples");
      delete dbline;
      return nullptr;
    }
    return dbline;
  }

  DbLine* DbLine::createFromSamplesById(Id nech,
                                        const ELoadBy& order,
                                        const VectorDouble& tab,
                                        const VectorInt& lineIds,
                                        const VectorInt& ranksPerId,
                                        const VectorString& names,
                                        const VectorString& locatorNames,
                                        bool flagAddSampleRank)
  {
    DbLine* dbline = new DbLine;
    if (dbline->resetFromSamplesById(nech, order, tab, lineIds, ranksPerId, names,
                                     locatorNames, flagAddSampleRank))
    {
      messerr("Error when creating DbLine from Samples By Ids");
      delete dbline;
      return nullptr;
    }
    return dbline;
  }

  Id DbLine::_lineLinkage(const VectorInt& lineCounts)
  {
    // Prelimnary check
    Id nech = VH::cumul(lineCounts);
    if (nech != getNSample())
    {
      messerr("Cumulated number of samples given by 'lineCounts' (%d) should "
              "match the number of samples (%d)",
              nech, getNSample());
      return 1;
    }

    // Count the number of lines
    Id nbline = static_cast<Id>(lineCounts.size());

    // Create the Linkage
    _lineAdds.resize(nbline, 0);

    // Loop over the lines
    Id start = 0;
    for (Id iline = 0; iline < nbline; iline++)
    {
      _lineAdds[iline] = VH::sequence(lineCounts[iline], start);
      start += lineCounts[iline];
    }
    return 0;
  }

  Id DbLine::_lineLinkageById(const VectorInt& linesId,
                               const VectorInt& ranksPerId)
  {
    auto nech = getNSample();

    // Preliminary checks by dimensions
    if (static_cast<Id>(linesId.size()) != nech)
    {
      messerr("Dimension of 'linesId' (%d) should match Number of samples (%d)",
              static_cast<Id>(linesId.size()), nech);
      return 1;
    }
    if (static_cast<Id>(ranksPerId.size()) != nech)
    {
      messerr("Dimension of 'ranksPerId' (%d) should match Number of samples (%d)",
              static_cast<Id>(ranksPerId.size()), nech);
      return 1;
    }

    // Find the number of lines
    VectorInt allLines = VH::unique(linesId);
    Id nbline         = static_cast<Id>(allLines.size());

    // Create the Linkage
    _lineAdds.resize(nbline, 0);

    for (Id iline = 0; iline < nbline; iline++)
    {
      Id refLineId = allLines[iline];

      VectorInt ranks;
      VectorInt iadds;
      for (Id iech = 0; iech < nech; iech++)
      {
        if (linesId[iech] != refLineId) continue;
        ranks.push_back(ranksPerId[iech]);
        iadds.push_back(iech);
      }

      VectorInt sortedRanks = VH::orderRanks(ranks);
      _lineAdds[iline]      = VH::reorder(iadds, sortedRanks);
    }
    return static_cast<Id>(!isConsistent());
  }

  /**
   * @brief Reset the contents of a DbLine from arguments (previous contents is
   * cleared beforehand). The line contents is provided in 'lineCounts'.
   *
   * @param nech Number of samples to be loaded
   * @param order Ordering mode used for storing in 'tab' (by column or by sample)
   * @param tab Vector containing the values to be imported
   * @param lineCounts Vector giving the number of samples per Line (see details)
   * @param names Names given to the output variables
   * @param locatorNames Name of the locators given to the output variables
   * @param flagAddSampleRank When TRUE, the 'rank' variable is added
   * @return Id Error returned code
   *
   * @details: Argument 'lineCounts' give the number of samples per Line.
   * @details: This assumes that samples of per line are ordered sequentially
   * @details and that samples of Line 'j' are followed by those of Line 'j+1'.
   */
  Id DbLine::resetFromSamples(Id nech,
                               const ELoadBy& order,
                               const VectorDouble& tab,
                               const VectorInt& lineCounts,
                               const VectorString& names,
                               const VectorString& locatorNames,
                               bool flagAddSampleRank)
  {
    if (Db::resetFromSamples(nech, order, tab, names, locatorNames,
                             flagAddSampleRank) != 0)
      return 1;

    // Create the Line Linkage

    if (_lineLinkage(lineCounts) != 0) return 1;

    return 0;
  }

  /**
   * @brief Reset the contents of a DbLine from arguments (previous contents is
   * cleared beforehand). The line contents is provided in 'lineIds' and
  'ranksPerId'
   *
   * @param nech Number of samples to be loaded
   * @param order Ordering mode used for storing in 'tab' (by column or by sample)
   * @param tab Vector containing the values to be imported
   * @param lineIds Vector giving the LineId to which each sample belongs (see details)
   * @param ranksPerId Vector giving the ordering of samples within Line (see details)
   * @param names Names given to the output variables
   * @param locatorNames Name of the locators given to the output variables
   * @param flagAddSampleRank When TRUE, the 'rank' variable is added
   * @return Id Error returned code
   *
   * @details: Argument 'lineIds' is dimensioned to the total number of samples.
   * @details: For each sample, it gives Id of Line to which the sample belongs.
   * @details: LineId must be numeric: equzl for samples of the same line and
   * @details: different for samples of different lines
   *
   * @details: Argument 'ranksPerId' is dimensionned to total number of samples.
   * @details: Along the samples belonging to one line (sharing the same LineId)
   * @details: it should provide the ordering of the samples.
   * @details: For one line, the values of 'ranksPerId' must be numeric:
   * @details: they do not need to be consecutive ... simply ordered.
   */
  Id DbLine::resetFromSamplesById(Id nech,
                                   const ELoadBy& order,
                                   const VectorDouble& tab,
                                   const VectorInt& lineIds,
                                   const VectorInt& ranksPerId,
                                   const VectorString& names,
                                   const VectorString& locatorNames,
                                   bool flagAddSampleRank)
  {
    if (Db::resetFromSamples(nech, order, tab, names, locatorNames,
                             flagAddSampleRank) != 0)
      return 1;

    // Create the Line Linkage

    if (_lineLinkageById(lineIds, ranksPerId) != 0) return 1;

    return 0;
  }

  bool DbLine::_deserializeAscii(std::istream& is, bool verbose)
  {
    Id ndim   = 0;
    Id nbline = 0;
    Id number = 0;
    VectorString locators;
    VectorString names;
    VectorDouble values;
    VectorDouble allvalues;

    /* Initializations */

    bool ret = true;
    ret      = ret && _recordRead<Id>(is, "Space Dimension", ndim);

    // Writing the set of addresses for Line organization

    ret = ret && _recordRead<Id>(is, "Number of Lines", nbline);
    _lineAdds.resize(nbline);
    for (Id iline = 0; iline < nbline; iline++)
    {
      ret = ret && _recordRead<Id>(is, "Number of Samples", number);
      ret = ret && _recordReadVec<Id>(is, "", _lineAdds[iline], number);
    }
    ret = ret && Db::_deserializeAscii(is, verbose);

    return ret;
  }

  bool DbLine::_serializeAscii(std::ostream& os, bool verbose) const
  {
    bool ret = true;

    /* Writing the header */

    ret = ret && _recordWrite<Id>(os, "Space Dimension", getNDim());

    // Writing the set of addresses for Line organization

    ret = ret && _recordWrite<Id>(os, "Number of Lines", getNLine());
    for (Id iline = 0, nbline = getNLine(); iline < nbline; iline++)
    {
      ret = ret && _recordWrite<Id>(os, "Number of Samples", getNSamplePerLine(iline));
      ret = ret && _recordWriteVec<Id>(os, "", _lineAdds[iline]);
    }

    /* Writing the tail of the file */

    ret&& Db::_serializeAscii(os, verbose);

    return ret;
  }

  /**
   * Create a Db by loading the contents of a Neutral File
   *
   * @param NFFilename Name of the Neutral File (Db format)
   * @param verbose    Verbose
   *
   * @remarks The name does not need to be completed in particular when defined by absolute path
   * @remarks or read from the Data Directory (in the gstlearn distribution)
   */
  DbLine* DbLine::createFromNF(const String& NFFilename, bool verbose)
  {
    DbLine* dbline = new DbLine;
    if (dbline->_fileOpenAndDeserialize(NFFilename, verbose)) return dbline;
    delete dbline;
    return nullptr;
  }

  /**
   * @brief Create a DbLine from the following information provided as input
  arguments
   *
   * @param ndim  Space dimension
   * @param nbline Number of Lines
   * @param nperline Average number of samples per line
   * @param deltaX Average distance between Lines along first space dimension
   * @param delta Average distances between samples along each line (in all directions)
   * @param unifDelta 5half-) width of uniform distribution
   * @param seed Seed used for the random number generator
   * @return DbLine* Pointer to the newly created DbLine structure
   */
  DbLine* DbLine::createFillRandom(Id ndim,
                                   Id nbline,
                                   Id nperline,
                                   double deltaX,
                                   const VectorDouble& delta,
                                   double unifDelta,
                                   Id seed)
  {
    law_set_random_seed(seed);

    // Origin of the lines
    VectorDouble d = delta;
    if (d.empty()) d.resize(ndim, 1.);
    VectorDouble shift = d;
    shift[0]           = 0.;
    VectorVectorDouble coor0(nbline, 0.);
    VectorVectorDouble incr0(nbline, 0.);
    for (Id iline = 0; iline < nbline; iline++)
    {
      coor0[iline].resize(ndim);
      for (Id idim = 0; idim < ndim; idim++)
      {
        coor0[iline][idim] = (idim == 0)
                             ? deltaX * iline + deltaX * law_uniform(1. - unifDelta, 1. + unifDelta)
                             : 0.;
      }
    }

    // Creating the coordinates
    Id nech = 0;
    VectorDouble tab;
    VectorInt lineCounts;
    for (Id iline = 0; iline < nbline; iline++)
    {
      Id nsample = nperline * law_uniform(1. - unifDelta, 1. + unifDelta);
      nech += nsample;
      lineCounts.push_back(nsample);

      // Generate the coordinates along the line
      for (Id is = 0; is < nsample; is++)
        for (Id idim = 0; idim < ndim; idim++)
        {
          double value = coor0[iline][idim] + is * shift[idim] +
                         d[idim] * law_uniform(1 - unifDelta, 1. + unifDelta);
          tab.push_back(value);
        }
    }

    VectorString names    = generateMultipleNames("x", ndim);
    VectorString locnames = generateMultipleNames(String {ELoc::X.getKey()}, ndim, "");
    DbLine* dbline        = createFromSamples(nech, ELoadBy::SAMPLE, tab, lineCounts, names, locnames);

    return dbline;
  }

  /**
   * @brief Check if the contents of private member of this class is compatible
   * with the number of samples stored in the Db
   * @return true if everything is OK; false if a problem occurs
   */
  bool DbLine::isConsistent() const
  {
    // Check on the count of addresses
    auto nech = getNSample();
    if (nech != getNTotal())
    {
      messerr("The number of samples contained in the Db (%d)",
              getNSample());
      messerr("is not equal to the number of addresses referenced in DbLine (%d)",
              getNTotal());
      return false;
    }

    // Check that all addresses are reached
    VectorBool isReached(nech, false);
    for (Id iline = 0, nbline = getNLine(); iline < nbline; iline++)
    {
      for (Id i = 0, number = getNSamplePerLine(iline); i < number; i++)
      {
        Id iadd = _lineAdds[iline][i];
        if (isReached[iadd])
        {
          messerr("Sample %d is reached twice:", iadd);
          messerr("- Line %d:", iline);
          VH::dump("Adds_1", _lineAdds[iline]);
          auto jline = getLineBySample(iadd);
          messerr("- Line %d:", jline);
          VH::dump("Adds_1", _lineAdds[jline]);
          return false;
        }
      }
    }
    return true;
  }

  /**
   * @brief Returns the rank of the line containing the target address
   *
   * @param iech Target address
   * @return Id Returne line number
   */
  Id DbLine::getLineBySample(Id iech) const
  {
    for (Id iline = 0, nbline = getNLine(); iline < nbline; iline++)
    {
      Id rank = VH::whereElement(_lineAdds[iline], iech);
      if (rank >= 0) return iline;
    }
    return -1;
  }

  VectorDouble DbLine::_getHeaderCoordinate(Id idim) const
  {
    auto nbline = getNLine();
    VectorDouble vec(nbline);
    for (Id iline = 0; iline < nbline; iline++)
    {
      Id iech   = _lineAdds[iline][0];
      vec[iline] = getCoordinate(iech, idim);
    }
    return vec;
  }

  VectorDouble DbLine::getCoordinatesPerLine(Id iline, Id idim) const
  {
    VectorDouble vec;
    if (!_isLineNumberValid(iline)) return vec;

    auto number = getNSamplePerLine(iline);
    vec.resize(number);
    for (Id i = 0; i < number; i++)
      vec[i] = getCoordinate(_lineAdds[iline][i], idim);

    return vec;
  }

  /**
   * @brief This is an example for a future more sophisticated method
   * which will collect statistics calculated per line, and store them into a newly
   * created Db.
   * In the current version, the statistics only concerns the number of samples per Line
   *
   * @return Db* Resulting Db
   */
  Db* DbLine::createStatToHeader() const
  {
    // Create the resulting output Db
    auto* db = new Db();

    // Glue the coordinates
    for (Id idim = 0, ndim = getNDim(); idim < ndim; idim++)
    {
      VectorDouble tab = _getHeaderCoordinate(idim);
      String name      = concatenateString("x", idim + 1);
      db->addColumns(tab, name, ELoc::X, idim);
    }

    // Add the line length as variable
    auto nbline = getNLine();
    VectorDouble tab(nbline);
    for (Id iline = 0; iline < nbline; iline++)
      tab[iline] = getNSamplePerLine(iline);
    db->addColumns(tab, "Count");

    return db;
  }

  /**
   * @brief Returns the absolute rank of the sample 'isample' or the line 'iline'
   * within the Db structure (ir -1 if an error occurs)
   *
   * @param iline Target line number
   * @param isample Target sample number within line
   * @return Rank of the sample
   */
  Id DbLine::getLineSampleRank(Id iline, Id isample) const
  {
    if (iline < 0 || iline >= getNLine())
    {
      messerr("Error in Line number (%d): it must lie within [0, %d]\n",
              iline, getNLine());
      return -1;
    }
    auto nsample = getNSamplePerLine(iline);
    if (isample < 0 || isample >= nsample)
    {
      messerr(
        "Error in Sample number (%d) in line (%d): it must lie within [0, %d]\n",
        isample, iline, nsample);
      return -1;
    }
    return _lineAdds[iline][isample];
  }

  DbLine* DbLine::createVerticalFromGrid(const DbGrid& grid,
                                         const VectorString& names,
                                         const VectorInt& xranks,
                                         const VectorInt& yranks,
                                         Id byZ)
  {
    // Preliminary checks
    auto ndim = grid.getNDim();
    if (ndim != 3)
    {
      messerr("This method is coded to extract wells from a 3-D Grid only");
      return nullptr;
    }
    if (static_cast<Id>(xranks.size()) != static_cast<Id>(yranks.size()))
    {
      messerr("Arguments 'xranks' and 'yranks' should have same dimensions");
      return nullptr;
    }
    Id nvar    = static_cast<Id>(names.size());
    Id nwells  = static_cast<Id>(xranks.size());
    auto nz     = grid.getNX(2);
    Id nbywell = nz / byZ;
    Id nsample = nwells * nbywell;
    VectorDouble tab(nsample * (3 + nvar));
    VectorInt lineCounts(nwells);

    VectorDouble coor(3);
    VectorInt indg(3);

    // Loop on the wells
    Id nech = 0;
    Id ecr  = 0;
    for (Id iwell = 0; iwell < nwells; iwell++)
    {
      indg[0] = xranks[iwell];
      indg[1] = yranks[iwell];

      // Loop on the samples
      for (Id iz = 0; iz < nbywell; iz++)
      {
        indg[2] = iz * byZ;

        // Assign the coordinates
        grid.indicesToCoordinateInPlace(indg, coor);
        for (Id idim = 0; idim < ndim; idim++) tab[ecr++] = coor[idim];

        // Assign the variable values
        Id rank = grid.indiceToRank(indg);
        for (Id ivar = 0; ivar < nvar; ivar++)
          tab[ecr++] = grid.getValue(names[ivar], rank);
        nech++;
      }
      lineCounts[iwell] = nbywell;
    }

    // Constitute the list of names
    VectorString locnames = generateMultipleNames("x", ndim);
    for (Id ivar = 0; ivar < nvar; ivar++)
      locnames.push_back((names[ivar]));

    DbLine* dbline = new DbLine;
    if (dbline->resetFromSamples(nech, ELoadBy::SAMPLE, tab, lineCounts, locnames))
      return nullptr;

    return dbline;
  }

  DbLine* DbLine::createMarkersFromGrid(const DbGrid& grid,
                                        const String& name,
                                        const VectorInt& xranks,
                                        const VectorInt& yranks,
                                        const VectorDouble& cuts)
  {
    // Preliminary checks
    auto ndim = grid.getNDim();
    if (ndim != 3)
    {
      messerr("This method is coded to extract wells from a 3-D Grid only");
      return nullptr;
    }
    if (static_cast<Id>(xranks.size()) != static_cast<Id>(yranks.size()))
    {
      messerr("Arguments 'xranks' and 'yranks' should have same dimensions");
      return nullptr;
    }
    Id ncuts  = static_cast<Id>(cuts.size());
    Id nwells = static_cast<Id>(xranks.size());
    auto nz    = grid.getNX(2);
    VectorDouble tab;
    VectorInt lineCounts(nwells);
    VectorDouble coor(3);
    VectorDouble cooriz(3);
    VectorDouble coorjz(3);
    VectorDouble well(nz);
    VectorInt indg(3);

    // Loop on the wells
    Id nech = 0;
    for (Id iwell = 0; iwell < nwells; iwell++)
    {
      indg[0] = xranks[iwell];
      indg[1] = yranks[iwell];

      // Loop on the samples
      for (Id iz = 0; iz < nz; iz++)
      {
        indg[2]  = iz;
        Id rank = grid.indiceToRank(indg);
        well[iz] = grid.getValue(name, rank);
      }

      // Find the markers
      Id nmark = 0;
      for (Id iz = 1; iz < nz; iz++)
      {
        Id jz = iz - 1;

        // Loop on the cuts
        for (Id icut = 0; icut < ncuts; icut++)
        {
          double zcut  = cuts[icut];
          double delta = zcut - well[jz];
          if ((zcut - well[iz]) * delta > 0) continue;

          // Define the marker by interpolation
          double dist  = well[iz] - well[jz];
          double ratio = (dist > 0) ? delta / dist : 0.;

          // Interpolate the coordinates
          indg[2] = jz;
          grid.indicesToCoordinateInPlace(indg, coorjz);
          indg[2] = iz;
          grid.indicesToCoordinateInPlace(indg, cooriz);
          for (Id idim = 0; idim < ndim; idim++)
            coor[idim] = coorjz[idim] * (1. - ratio) + cooriz[idim] * ratio;

          // Add the sample
          for (Id idim = 0; idim < ndim; idim++) tab.push_back(coor[idim]);
          tab.push_back(zcut);
          nech++;
          nmark++;
        }
      }
      lineCounts[iwell] = nmark;
    }

    // Constitute the list of names
    VectorString locnames = generateMultipleNames("x", ndim);
    VectorString auxnames = generateMultipleNames("cut", ncuts);
    for (Id icut = 0; icut < ncuts; icut++) locnames.push_back(auxnames[icut]);

    DbLine* dbline = new DbLine;
    if (dbline->resetFromSamples(nech, ELoadBy::SAMPLE, tab, lineCounts,
                                 locnames)) return nullptr;

    return dbline;
  }
#ifdef HDF5
  bool DbLine::_deserializeH5(H5::Group& grp, [[maybe_unused]] bool verbose)
  {
    auto dbG = SerializeHDF5::getGroup(grp, "DbLine");
    if (!dbG) return false;

    /* Read the grid characteristics */
    bool ret   = true;
    Id ndim   = 0;
    Id nbline = 0;

    ret = ret && SerializeHDF5::readValue(*dbG, "NDim", ndim);
    ret = ret && SerializeHDF5::readValue(*dbG, "NLines", nbline);

    auto linesG = SerializeHDF5::getGroup(*dbG, "Lines");
    if (!linesG) return false;
    _lineAdds.resize(nbline);
    for (Id iline = 0; iline < nbline; iline++)
    {
      String locName = "Line" + std::to_string(iline);
      auto lineg     = SerializeHDF5::getGroup(*linesG, locName);
      if (!lineg) return false;

      Id nsample = 0;
      ret         = ret && SerializeHDF5::readValue(*lineg, "NSamples", nsample);
      ret         = ret && SerializeHDF5::readVec(*lineg, "Samples", _lineAdds[iline]);
    }

    /* Writing the tail of the file */

    ret = ret && Db::_deserializeH5(*dbG, verbose);

    return ret;
  }

  bool DbLine::_serializeH5(H5::Group& grp, [[maybe_unused]] bool verbose) const
  {
    auto dbG = grp.createGroup("DbLine");

    bool ret = true;

    ret = ret && SerializeHDF5::writeValue(dbG, "NDim", getNDim());
    ret = ret && SerializeHDF5::writeValue(dbG, "NLines", getNLine());

    auto linesG = dbG.createGroup("Lines");
    for (Id iline = 0, nbline = getNLine(); iline < nbline; iline++)
    {
      String locName = "Line" + std::to_string(iline);
      auto lineG     = linesG.createGroup(locName);

      ret = ret && SerializeHDF5::writeValue(lineG, "NSamples", getNSamplePerLine(iline));
      ret = ret && SerializeHDF5::writeVec(lineG, "Samples", _lineAdds[iline]);
    }

    /* Writing the tail of the file */

    ret = ret && Db::_serializeH5(dbG, verbose);

    return ret;
  }
#endif
} // namespace gstlrn