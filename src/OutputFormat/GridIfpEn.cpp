/******************************************************************************/
/* gstlearn C++ Library                                                       */
/*                                                                            */
/* Authors: <authors>                                                         */
/*                                                                            */
/* License: BSD 3 Clause                                                      */
/******************************************************************************/
#include "geoslib_f_private.h"

#include "OutputFormat/AOF.hpp"
#include "Db/Db.hpp"
#include "Db/DbGrid.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/String.hpp"

#include <string.h>
#include "OutputFormat/GridIfpEn.hpp"

GridIfpEn::GridIfpEn(const char* filename, const Db* db)
  : AOF(filename, db)
{
}

GridIfpEn::GridIfpEn(const GridIfpEn& r)
    : AOF(r)
{
}

GridIfpEn& GridIfpEn::operator=(const GridIfpEn& r)
{
  if (this != &r)
  {
    AOF::operator=(r);
  }
  return *this;
}

GridIfpEn::~GridIfpEn()
{
}

int GridIfpEn::writeInFile()
{
  static double valnull = 3.0;

  /* Open the file */

  if (_fileWriteOpen()) return 1;

  // Preliminary calculations
  int ncol = (int) _cols.size();
  VectorInt nx = _dbgrid->getNXsExt(3);
  VectorDouble angles = _dbgrid->getAngles();
  int ntot = 1;
  for (int idim = 0; idim < 3; idim++)
  {
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
void GridIfpEn::_writeLine(int mode,
                          const char *comment,
                          int valint,
                          double valrel,
                          const char *combis)
{
  std::stringstream sstr;

  //char line[1000];

  /* Initialize the string */

  //(void) gslStrcpy(line, "");

  /* Comment */

  if (comment != NULL)
    //(void) gslSPrintf(&line[strlen(line)], "%s", comment);
    sstr << comment;

  /* Encoding the value */

  if (mode == 1)
  {
    // (void) gslSPrintf(&line[strlen(line)], " %d", valint);
    sstr << " " << valint;
  }
  else if (mode == 2)
  {
    // (void) gslSPrintf(&line[strlen(line)], " %lf", valrel);
    sstr << " " << valrel;
  }

  /* Secondary comment */

  if (combis != NULL)
    //(void) gslSPrintf(&line[strlen(line)], " %s", combis);
    sstr << " " << combis;

  /* Print the line */

  //fprintf(_file, "%s\n", line);
  fprintf(_file, "%s\n", sstr.str().c_str());
}

DbGrid* GridIfpEn::readGridFromFile()
{
  DbGrid* dbgrid = nullptr;
  int dumint, ncol;
  double dumrel, test, value;
  VectorDouble x0(3);
  VectorDouble dx(3);
  VectorDouble angles(3);
  VectorInt nx(3);

  /* Open the file */

  if (_fileReadOpen()) return dbgrid;

  /* Read the grid characteristics */

  for (int i = 0; i < 3; i++)
  {
    nx[i] = 1;
    dx[i] = 1.;
    x0[i] = 0.;
    angles[i] = 0.;
  }

  /* Read the header */

  if (_readLine( 0, "##########################", &dumint, &dumrel))    return dbgrid;
  if (_readLine( 0, "FILE_DESCRIPTION         #", &dumint, &dumrel))    return dbgrid;
  if (_readLine( 0, "APPLICATION              #", &dumint, &dumrel))    return dbgrid;
  if (_readLine( 0, "SURVEY_NAME              #", &dumint, &dumrel))    return dbgrid;
  if (_readLine( 0, "MATRIX_NAME              #", &dumint, &dumrel))    return dbgrid;
  if (_readLine( 0, "METHOD                   #", &dumint, &dumrel))    return dbgrid;
  if (_readLine( 2, "FLOAT_NULL_VALUE         #", &dumint, &test))      return dbgrid;
  if (_readLine( 0, "ROW_COLUMN_ORIENTATION   #", &dumint, &dumrel))    return dbgrid;
  if (_readLine( 0, "REPRESENTATION_CODE      #", &dumint, &dumrel))    return dbgrid;
  if (_readLine( 0, "##########################", &dumint, &dumrel))    return dbgrid;
  if (_readLine( 2, "ANGLE                    #", &dumint, &angles[0])) return dbgrid;
  if (_readLine( 1, "ROW_COUNT                #", &nx[1], &dumrel))     return dbgrid;
  if (_readLine( 1, "COLUMN_COUNT             #", &nx[0], &dumrel))     return dbgrid;
  if (_readLine( 2, "ROW_DISTANCE             #", &dumint, &dx[1]))     return dbgrid;
  if (_readLine( 2, "COLUMN_DISTANCE          #", &dumint, &dx[0]))     return dbgrid;
  if (_readLine( 1, "LAYER_COUNT              #", &nx[2], &dumrel))     return dbgrid;
  if (_readLine( 2, "X_ORIGIN                 #", &dumint, &x0[0]))     return dbgrid;
  if (_readLine( 2, "Y_ORIGIN                 #", &dumint, &x0[1]))     return dbgrid;
  if (_readLine( 1, "FACIES_COUNT             #", &ncol, &dumrel))      return dbgrid;
  if (_readLine( 0, "DATA_PROP                #", &dumint, &dumrel))    return dbgrid;
  if (_readLine( 0, "##########################", &dumint, &dumrel))    return dbgrid;

  /* Read the array of real values */

  int lec = 0;
  int nech = nx[0] * nx[1] * nx[2];
  VectorDouble tab(nech * ncol);

  for (int icol = 0; icol < ncol; icol++)
    for (int ix = 0; ix < nx[0]; ix++)
      for (int iy = 0; iy < nx[1]; iy++)
        for (int iz = 0; iz < nx[2]; iz++)
        {
          if (_record_read(_file, "%lf", &value)) break;
          if (value == test) value = TEST;
          tab[lec] = value;
          lec++;
        }

  if (lec != nech * ncol)
  {
    messerr("Number of decoded values (%d) is not equal to the number of grid nodes (%d) x number of attributes (%d)",
            lec, nech, ncol);
    return dbgrid;
  }

  /* Set the error return code */

  dbgrid = new DbGrid();
  VectorString names = generateMultipleNames("IfpEn", ncol);
  dbgrid->reset(nx,dx,x0,angles,ELoadBy::SAMPLE,tab,names);

  // Close the file

  _fileClose();

  return dbgrid;
}

/****************************************************************************/
/*!
 **   Decode a line for IFPEN file
 **
 ** \param[in]  mode       Type of encoding
 ** \li                     0 : Comment
 ** \li                     1 : Integer value
 ** \li                     2 : Real value
 ** \param[in]  comment    Comment string (or NULL)
 **
 ** \param[out]  valint     Integer value
 ** \param[out]  valrel     Float value
 **
 *****************************************************************************/
int GridIfpEn::_readLine(int mode,
                         const char *comment,
                         int *valint,
                         double *valrel)
{
  char line[100];
  int start;

  /* Reading the line */

  if (fgets(line, 100, _file) == NULL) return (1);
  line[strlen(line) - 1] = '\0';

  /* Check the comment */

  start = 0;
  if (comment != NULL)
  {
    if (strcmp(line, comment) < 0) return (1);
    start = static_cast<int>(strlen(comment));
  }

  /* Decoding the value */

  if (mode == 1)
  {
    if (gslSScanf(&line[start], "%d", valint) != 1) return (1);
  }
  else if (mode == 2)
  {
    if (gslSScanf(&line[start], "%lf", valrel) != 1) return (1);
  }
  return (0);
}
