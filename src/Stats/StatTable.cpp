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
#include "Stats/StatTable.hpp"
#include "Basic/String.hpp"

StatTable::StatTable()
  : ASerializable(),
    _stats()
{
}

StatTable::StatTable(const StatTable &m)
    : ASerializable(m),
      _stats(m._stats)
{

}

StatTable& StatTable::operator=(const StatTable &m)
{
  if (this != &m)
  {
    _stats = m._stats;
  }
  return *this;
}

StatTable::~StatTable()
{

}

int StatTable::getRowNumber() const
{
  if (isEmpty()) return 0;
  return _stats[0].size();
}

int StatTable::getColNumber() const
{
  if (isEmpty()) return 0;
  return _stats.size();
}

void StatTable::resize(int irow, int ncols)
{
  int nrows = 0;
  if (isEmpty())
    _stats.resize(ncols);
  else
    nrows = getRowNumber();

  if (irow >= nrows)
  {
    int ncols = getColNumber();
    for (int icol = 0; icol < ncols; icol++)
      _stats[icol].resize(irow);
  }
}

double StatTable::getValue(int irow, int icol) const
{
  if (irow <= getRowNumber()) return 0.;
  return _stats[icol][irow];
}

int StatTable::serialize(const String& filename, bool verbose) const
{

  /* Opening the Data file */

  if (_fileOpen(filename, "Table", "w", verbose)) return 1;

  /* Writing the header */

  _recordWrite("%d", getColNumber());
  _recordWrite("#", "Number of Columns");
  _recordWrite("%d", getRowNumber());
  _recordWrite("#", "Number of Rows");

  /* Writing the tail of the file */

  for (int irow = 0; irow < getRowNumber(); irow++)
  {
    for (int icol = 0; icol < getColNumber(); icol++)
    {
      _recordWrite("%lf", _stats[icol][irow]);
    }
   _recordWrite("\n");
  }

  // Close the Neutral file
  _fileClose(verbose);

  return 0;
}

int StatTable::deSerialize(const String& filename, bool verbose)
{
  int ncols, nrows;
  double value;

  /* Opening the Data file */

  if (_fileOpen(filename, "Table", "r", verbose)) return 1;
  int error = 1;

  /* Decoding the header */

  if (_recordRead("Number of Columns", "%d", &ncols)) goto label_end;
  if (_recordRead("Number of Rows", "%d", &nrows)) goto label_end;

  _stats.clear();
  _stats.resize(ncols);

  /* Loop on the lines */

  for (int irow = 0; irow < nrows; irow++)
  {
    for (int icol = 0; icol < ncols; icol++)
    {
      if (_recordRead("Numerical value", "%lf", &value)) goto label_end;
      _stats[icol].push_back(value);
    }
  }

  error = 0;

  label_end:
  if (error) _stats.clear();
  return error;
}

/**
 * Print the contents of the statistics
 */
void StatTable::display(int isimu) const
{
  if (_stats.empty()) return;

  mestitle(1,"Statistics printout for Simulation %d",isimu);

  int ncols = getColNumber();
  int nrows = getRowNumber();

  for (int irow = 0; irow < nrows; irow++)
  {
    for (int icol = 0; icol < ncols; icol++)
    {
      message(" %10.3lf",_stats[icol][irow]);
    }
    message("\n");
  }
}

/**
 * Plot the contents of the statistics
 */
void StatTable::plot(int isimu) const
{
  if (_stats.empty()) return;
  String filename = incrementStringVersion("GibbsStats",isimu+1);
  serialize(filename,false);
}
