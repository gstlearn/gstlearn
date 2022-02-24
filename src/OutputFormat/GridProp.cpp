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
#include "OutputFormat/GridProp.hpp"
#include "OutputFormat/AOF.hpp"
#include "Db/Db.hpp"
#include "Db/DbGrid.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/String.hpp"

#include <string.h>

GridProp::GridProp(const char* filename, const Db* db)
  : AOF(filename, db)
{
}

GridProp::GridProp(const GridProp& r)
    : AOF(r)
{
}

GridProp& GridProp::operator=(const GridProp& r)
{
  if (this != &r)
  {
    AOF::operator=(r);
  }
  return *this;
}

GridProp::~GridProp()
{
}

int GridProp::dumpFile()
{
  static double valnull = 3.0;

  /* Open the file */

  if (_fileOpen()) return 1;

  // Preliminary calculations

  int ncol = _cols.size();
  int ndim = _dbgrid->getNDim();
  VectorInt nx = _dbgrid->getNXs();
  VectorDouble angles = _dbgrid->getAngles();
  int ntot = 1;
  for (int idim = 0; idim < 3; idim++)
  {
    nx[idim] = (idim < ndim) ? nx[idim] : 1;
    ntot *= nx[idim];
  }

  /* Write the header */

  _writeLine( 0, "##########################", 0, 0., NULL);
  _writeLine( 0, "FILE_DESCRIPTION         # PROP", 0, 0., NULL);
  _writeLine( 0, "APPLICATION              #", 0, 0., "# CobraFlow");
  _writeLine( 0, "SURVEY_NAME              #", 0, 0., NULL);
  _writeLine( 0, "MATRIX_NAME              # VPCMatrix_test_export", 0,
                 0., NULL);
  _writeLine( 0, "METHOD                   # BY_CPV", 0, 0., NULL);
  _writeLine( 2, "FLOAT_NULL_VALUE         #", 0, valnull, NULL);
  _writeLine( 0, "ROW_COLUMN_ORIENTATION   # ROW", 0, 0., NULL);
  _writeLine( 0, "REPRESENTATION_CODE      # ASCII", 0, 0., NULL);
  _writeLine( 0, "##########################", 0, 0., NULL);
  _writeLine( 2, "ANGLE                    #", 0, angles[0], "# DEG");
  _writeLine( 1, "ROW_COUNT                #", nx[1], 0., NULL);
  _writeLine( 1, "COLUMN_COUNT             #", nx[0], 0., NULL);
  _writeLine( 2, "ROW_DISTANCE             #", 0, _dbgrid->getDX(1), "# m");
  _writeLine( 2, "COLUMN_DISTANCE          #", 0, _dbgrid->getDX(0), "# m");
  _writeLine( 1, "LAYER_COUNT              #", nx[2], 0., NULL);
  _writeLine( 2, "X_ORIGIN                 #", 0, _dbgrid->getX0(0), "# m");
  _writeLine( 2, "Y_ORIGIN                 #", 0, _dbgrid->getX0(1), "# m");
  _writeLine( 1, "FACIES_COUNT             #", ncol, 0., NULL);
  _writeLine( 0, "DATA_PROP                # CHANNEL1", 0, 0.,
                 "# Facies proportion");
  _writeLine( 0, "##########################", 0, 0., NULL);

  /* Grid description */

  for (int j = 0; j < ncol; j++)
    for (int i = 0; i < ntot; i++)
    {
      double value = _dbgrid->getArray(i, _cols[j]);
      _writeLine( 2, NULL, 0, value, NULL);
    }

  _fileClose();
  return 0;
}

/****************************************************************************/
/*!
 **   Encode a line for IFPEN file
 **
 ** \param[in]  mode       Type of encoding
 ** \li                     0 : Comment
 ** \li                     1 : Integer value
 ** \li                     2 : Real value
 ** \param[in]  comment    Comment string (or NULL)
 ** \param[in]  valint     Integer value
 ** \param[in]  valrel     Float value
 ** \param[in]  combis     Second comment (or NULL)
 **
 *****************************************************************************/
void GridProp::_writeLine(int mode,
                          const char *comment,
                          int valint,
                          double valrel,
                          const char *combis)
{
  char line[100];

  /* Initialize the string */

  (void) gslStrcpy(line, "");

  /* Comment */

  if (comment != NULL) (void) gslSPrintf(&line[strlen(line)], "%s", comment);

  /* Encoding the value */

  if (mode == 1)
  {
    (void) gslSPrintf(&line[strlen(line)], " %d", valint);
  }
  else if (mode == 2)
  {
    (void) gslSPrintf(&line[strlen(line)], " %lf", valrel);
  }

  /* Secondary comment */

  if (combis != NULL) (void) gslSPrintf(&line[strlen(line)], " %s", combis);

  /* Print the line */

  fprintf(_file, "%s\n", line);
}

