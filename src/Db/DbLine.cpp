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
#include "Basic/Law.hpp"
#include "Basic/String.hpp"
#include "Enum/ELoadBy.hpp"
#include "geoslib_define.h"

#include "Db/Db.hpp"
#include "Db/DbLine.hpp"
#include "Db/DbGrid.hpp"
#include "Db/DbStringFormat.hpp"
#include "Polygon/Polygons.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/VectorNumT.hpp"
#include "Stats/Classical.hpp"
#include "Space/SpacePoint.hpp"
#include "Enum/ELoc.hpp"

#include <math.h>

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
bool DbLine::_isLineNumberValid(int iline) const
{
  if (iline < 0)
  {
    messerr("Argument 'iline' should be non negative");
    return false;
  }
  if (iline >= getLineNumber())
  {
    messerr("ilin' (%d) should be smaller than Number of Lines (%d)", iline,
            getLineNumber());
    return false;
  }
  return true;
}

int DbLine::getLineNumber() const
{
  if (_lineAdds.empty()) return 0;
  return (int) _lineAdds.size();
}

int DbLine::getLineSampleNumber(int iline) const
{
  if (! _isLineNumberValid(iline)) return -1;
  return (int) _lineAdds[iline].size();
}

int DbLine::getNTotal() const
{
  int ntotal = 0;
  for (int iline = 0, nbline = getLineNumber(); iline < nbline; iline++)
    ntotal += getLineSampleNumber(iline);
  return ntotal;
}

double DbLine::getLineLength(int iline) const
{
  if (!_isLineNumberValid(iline)) return TEST;
  double total = 0.;
  SpacePoint P1;
  SpacePoint P2;
  P1.setIech(_lineAdds[iline][0]);
  getSampleAsSPInPlace(P1);
  for (int iech = 1, nech = getLineSampleNumber(iline); iech < nech; iech++)
  {
    P2.setIech(_lineAdds[iline][iech]);
    getSampleAsSPInPlace(P2);
    total += P2.getDistance(P1);
    P1 = P2;
  }
  return total;
}

VectorDouble DbLine::getLineLengths() const
{
  int nline = getLineNumber();
  VectorDouble lengths(nline);
  for (int iline = 0; iline < nline; iline++)
    lengths[iline] = getLineLength(iline);
  return lengths;
}

String DbLine::toString(const AStringFormat* strfmt) const
{
  std::stringstream sstr;

  const DbStringFormat* dbfmt = dynamic_cast<const DbStringFormat*>(strfmt);
  DbStringFormat dsf;
  if (dbfmt != nullptr) dsf = *dbfmt;

  sstr << toTitle(0, "Data Base Line Characteristics");

  sstr << "Number of Lines = " << getLineNumber() << std::endl;
  sstr << "Line length = ";
  for (int iline = 0, nbline = getLineNumber(); iline < nbline; iline++)
  {
    if (iline > 0) sstr << " / ";
    sstr << getLineSampleNumber(iline);
  }
  sstr << std::endl;

  sstr << _toStringCommon(&dsf);

  return sstr.str();
}

DbLine* DbLine::createFromSamples(int nech,
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

DbLine* DbLine::createFromSamplesById(int nech,
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

int DbLine::_lineLinkage(const VectorInt& lineCounts)
{
  // Prelimnary check
  int nech = VH::cumul(lineCounts);
  if (nech != getSampleNumber())
  {
    messerr("Cumulated number of samples given by 'lineCounts' (%d) should "
            "match the number of samples (%d)",
            nech, getSampleNumber());
    return 1;
  }

  // Count the number of lines
  int nbline = (int)lineCounts.size();

  // Create the Linkage
  _lineAdds.resize(nbline, 0);

  // Loop over the lines
  int start = 0;
  for (int iline = 0; iline < nbline; iline++)
  {
    _lineAdds[iline] = VH::sequence(lineCounts[iline], start);
    start += lineCounts[iline];
  }
  return 0;
}

int DbLine::_lineLinkageById(const VectorInt& linesId,
                             const VectorInt& ranksPerId)
{
  int nech = getSampleNumber();

  // Preliminary checks by dimensions
  if ((int)linesId.size() != nech)
  {
    messerr("Dimension of 'linesId' (%d) should match Number of samples (%d)",
            (int)linesId.size(), nech);
    return 1;
  }
  if ((int)ranksPerId.size() != nech)
  {
    messerr("Dimension of 'ranksPerId' (%d) should match Number of samples (%d)",
      (int)ranksPerId.size(), nech);
    return 1;
  }

  // Find the number of lines
  VectorInt allLines = VH::unique(linesId);
  int nbline         = (int)allLines.size();

  // Create the Linkage
  _lineAdds.resize(nbline, 0);

  for (int iline = 0; iline < nbline; iline++)
  {
    int refLineId = allLines[iline];

    VectorInt ranks;
    VectorInt iadds;
    for (int iech = 0; iech < nech; iech++)
    {
      if (linesId[iech] != refLineId) continue;
      ranks.push_back(ranksPerId[iech]);
      iadds.push_back(iech);
    }

    VectorInt sortedRanks  = VH::orderRanks(ranks);
    _lineAdds[iline]       = VH::reorder(iadds, sortedRanks);
  }
  return (int) !isConsistent();
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
 * @return int Error returned code
 *
 * @details: Argument 'lineCounts' give the number of samples per Line.
 * @details: This assumes that samples of per line are ordered sequentially
 * @details and that samples of Line 'j' are followed by those of Line 'j+1'.
 */
int DbLine::resetFromSamples(int nech,
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
 * @return int Error returned code
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
int DbLine::resetFromSamplesById(int nech,
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

bool DbLine::_deserialize(std::istream& is, bool verbose)
{
  int ndim = 0;
  int nbline = 0;
  int number = 0;
  VectorString locators;
  VectorString names;
  VectorDouble values;
  VectorDouble allvalues;

  /* Initializations */

  bool ret = true;
  ret = ret && _recordRead<int>(is, "Space Dimension", ndim);

  // Writing the set of addresses for Line organization

  ret      = ret && _recordRead<int>(is, "Number of Lines", nbline);
  _lineAdds.resize(nbline);
  for (int iline = 0; iline < nbline; iline++)
  {
    ret = ret && _recordRead<int>(is, "Number of Samples", number);
    ret = ret && _recordReadVec<int>(is, "", _lineAdds[iline], number);
  }
  ret = ret && Db::_deserialize(is, verbose);

  return ret;
}

bool DbLine::_serialize(std::ostream& os, bool verbose) const
{
  bool ret = true;

  /* Writing the header */

  ret = ret && _recordWrite<int>(os, "Space Dimension", getNDim());

  // Writing the set of addresses for Line organization

  ret      = ret && _recordWrite<int>(os, "Number of Lines", getLineNumber());
  for (int iline = 0, nbline = getLineNumber(); iline < nbline; iline++)
  {
    ret = ret && _recordWrite<int>(os, "Number of Samples", getLineSampleNumber(iline));
    ret = ret && _recordWriteVec<int>(os, "", _lineAdds[iline]);
  }

  /* Writing the tail of the file */

  ret && Db::_serialize(os, verbose);

  return ret;
}

/**
 * Create a Db by loading the contents of a Neutral File
 *
 * @param neutralFilename Name of the Neutral File (Db format)
 * @param verbose         Verbose
 *
 * @remarks The name does not need to be completed in particular when defined by absolute path
 * @remarks or read from the Data Directory (in the gstlearn distribution)
 */
DbLine* DbLine::createFromNF(const String& neutralFilename, bool verbose)
{
  DbLine* dbLine = nullptr;
  std::ifstream is;
  dbLine = new DbLine;
  bool success = false;
  if (dbLine->_fileOpenRead(neutralFilename, is, verbose))
  {
    success = dbLine->deserialize(is, verbose);
  }
  if (! success)
  {
    delete dbLine;
    dbLine = nullptr;
  }
  return dbLine;
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
DbLine* DbLine::createFillRandom(int ndim,
                                 int nbline,
                                 int nperline,
                                 double deltaX,
                                 const VectorDouble& delta,
                                 double unifDelta,
                                 int seed)
{
  law_set_random_seed(seed);

  // Origin of the lines
  VectorDouble d = delta;
  if (d.empty()) d.resize(ndim, 1.);
  VectorDouble shift = d;
  shift[0] = 0.;
  VectorVectorDouble coor0(nbline, 0.);
  VectorVectorDouble incr0(nbline, 0.);
  for (int iline = 0; iline < nbline; iline++)
  {
    coor0[iline].resize(ndim);
    for (int idim = 0; idim < ndim; idim++)
    {
      coor0[iline][idim] = (idim == 0)
          ? deltaX * iline + deltaX * law_uniform(1. - unifDelta, 1. + unifDelta) : 0.;
    }
  }

  // Creating the coordinates
  int nech = 0;
  VectorDouble tab;
  VectorInt lineCounts;
  for (int iline = 0; iline < nbline; iline++)
  {
    int nsample = nperline * law_uniform(1. - unifDelta, 1. + unifDelta);
    nech += nsample;
    lineCounts.push_back(nsample);

    // Generate the coordinates along the line
    for (int is = 0; is < nsample; is++)
      for (int idim = 0; idim < ndim; idim++)
      {
        double value = coor0[iline][idim] + is * shift[idim] +
                       d[idim] * law_uniform(1 - unifDelta, 1. + unifDelta);
        tab.push_back(value);
      }
  }

  VectorString names = generateMultipleNames("x", ndim);
  VectorString locnames = generateMultipleNames(String{ELoc::X.getKey()}, ndim, "");
  DbLine* dbline = createFromSamples(nech, ELoadBy::SAMPLE, tab, lineCounts, names, locnames);

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
  int nech = getSampleNumber();
  if (nech != getNTotal())
  {
    messerr("The number of samples contained in the Db (%d)",
            getSampleNumber());
    messerr("is not equal to the number of addresses referenced in DbLine (%d)",
            getNTotal());
    return false;
  }

  // Check that all addresses are reached
  VectorBool isReached(nech, false);
  for (int iline = 0, nbline = getLineNumber(); iline < nbline; iline++)
  {
    for (int i = 0, number = getLineSampleNumber(iline); i < number; i++)
    {
      int iadd = _lineAdds[iline][i];
      if (isReached[iadd])
      {
        messerr("Sample %d is reached twice:", iadd);
        messerr("- Line %d:", iline);
        VH::display("Adds_1", _lineAdds[iline]);
        int jline = getLineBySample(iadd);
        messerr("- Line %d:", jline);
        VH::display("Adds_1", _lineAdds[jline]);
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
 * @return int Returne line number
 */
int DbLine::getLineBySample(int iech) const
{
  for (int iline = 0, nbline = getLineNumber(); iline < nbline; iline++)
  {
    int rank = VH::whereElement(_lineAdds[iline], iech);
    if (rank >= 0) return iline;
  }
  return -1;
}

VectorDouble DbLine::_getHeaderCoordinate(int idim) const
{
  int nbline = getLineNumber();
  VectorDouble vec(nbline);
  for (int iline = 0; iline < nbline; iline++)
  {
    int iech = _lineAdds[iline][0];
    vec[iline] = getCoordinate(iech, idim);
  }
  return vec;
}

VectorDouble DbLine::getCoordinates(int iline, int idim) const
{
  VectorDouble vec;
  if (!_isLineNumberValid(iline)) return vec;

  int number = getLineSampleNumber(iline);
  vec.resize(number);
  for (int i = 0; i < number; i++)
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
  Db* db     = new Db();

  // Glue the coordinates
  for (int idim = 0, ndim = getNDim(); idim < ndim; idim++)
  {
    VectorDouble tab = _getHeaderCoordinate(idim);
    String name = concatenateString("x", idim+1);
    db->addColumns(tab, name, ELoc::X, idim);
  }

  // Add the line length as variable
  int nbline = getLineNumber();
  VectorDouble tab(nbline);
  for (int iline = 0; iline < nbline; iline++)
    tab[iline] = getLineSampleNumber(iline);
  db->addColumns(tab, "Count");

  return db;
}

/**
 * @brief Returns the absolute rank of the sample 'isample' or the line 'iline'
 * within the Db structure (ir -1 if an error occurs)
 * 
 * @param iline 
 * @param isample 
 * @return int 
 */
int DbLine::getLineSampleRank(int iline, int isample) const
{
  if (iline < 0 || iline >= getLineNumber())
  {
    messerr("Error in Line number (%d): it must lie within [0, %d]\n",
            iline, getLineNumber());
    return -1;
  }
  int nsample = getLineSampleNumber(iline);
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
                                       int byZ)
{
  // Preliminary checks
  int ndim = grid.getNDim();
  if (ndim != 3)
  {
    messerr("This method is coded to extract wells from a 3-D Grid only");
    return nullptr;
    }
  if ((int)xranks.size() != (int)yranks.size())
  {
    messerr("Arguments 'xranks' and 'yranks' should have same dimensions");
    return nullptr;
  }
  int nvar    = (int)names.size();
  int nwells  = (int)xranks.size();
  int nz      = grid.getNX(2);
  int nbywell = nz / byZ;
  int nsample = nwells * nbywell;
  VectorDouble tab(nsample * (3 + nvar));
  VectorInt lineCounts(nsample);

  VectorDouble coor(3);
  VectorInt indg(3);

  // Loop on the wells
  int nech = 0;
  int ecr = 0;
  for (int iwell = 0; iwell < nwells; iwell++)
  {
    indg[0] = xranks[iwell];
    indg[1] = yranks[iwell];

    // Loop on the samples
    for (int iz = 0; iz < nbywell; iz++) indg[2] = iz * byZ;

    // Assign the coordinates
    grid.indicesToCoordinateInPlace(indg, coor);
    for (int idim = 0; idim < ndim; idim++) tab[ecr++] = coor[idim];

    // Assign the variable values
    int rank = grid.indiceToRank(indg);
    for (int ivar = 0; ivar < nvar; ivar++)
      tab[ecr++] = grid.getValue(names[ivar], rank);

    // Concatenate to the array
    lineCounts[nech] = nech;
    nech++;
  }

  DbLine* dbline = new DbLine;
  if (dbline->resetFromSamples(nech, ELoadBy::SAMPLE, tab, lineCounts, names))
    return nullptr;

  return dbline;
}
