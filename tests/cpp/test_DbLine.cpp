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
#include "Enum/ESpaceType.hpp"
#include "Space/ASpaceObject.hpp"

#include "Basic/File.hpp"
#include "Db/Db.hpp"
#include "Db/DbLine.hpp"
#include "Db/DbStringFormat.hpp"

/****************************************************************************/
/*!
 ** Main Program
 **
 ** This program is meant to check the manipulation of the Db
 **
 *****************************************************************************/
int main(int argc, char *argv[])
{
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str(), argc, argv);

  ASerializable::setContainerName(true);

  // Creating the Lines (in a 1-D space)
  defineDefaultSpace(ESpaceType::RN, 1);

  // tab contains the
  VectorInt lineCounts = { 3, 4, 2, 3};
  VectorDouble x1      = { 1., 1., 1., 2., 2., 2., 2., 3., 3., 4., 4., 4.};
  VectorDouble x2      = { 1., 2., 3., 1., 2., 3., 4., 3., 4., 2., 3., 4.};
  VectorDouble z1      = {1.2, 2.5, 3.6, 1.4, 0.3, 0.2, 8.2, 0.3, 3.2, 1.2, 0.4, 0.1};
  VectorDouble z2      = {3.1, 2.7, 3.2, 8.1, 5.3, 7.2, 9.2, 1.1, 0.3, 0.6, 1.5, 5.2};
  VectorDouble tab;
  copy(x1.begin(), x1.end(), std::back_inserter(tab));
  copy(x2.begin(), x2.end(), std::back_inserter(tab));
  copy(z1.begin(), z1.end(), std::back_inserter(tab));
  copy(z2.begin(), z2.end(), std::back_inserter(tab));
  DbLine* dbline =
    DbLine::createFromSamples((int)x1.size(), ELoadBy::COLUMN, tab, lineCounts,
                              {"x1", "x2", "z1", "z2"},
                              {"x1", "x2", "z1", "z2"});
  if (dbline == nullptr) return 1;
  mestitle(1, "Reference DbLine File");
  dbline->display();

  // Dump into a Neutral File
  dbline->dumpToNF("Line.ascii");

  // Read from the Neutral File and print contents again
  DbLine* dbline2 = DbLine::createFromNF("Line.ascii");
  mestitle(1, "DbLine after Serialization / Deserialization");
  dbline2->display();
    
   // Checking the second way to initiate the DbLine (using lineIds)
  VectorInt lineIds    = {2, 2, 2, 5, 5, 5, 5, 1, 1, 6, 6, 6};
  VectorInt ranksPerId = {1, 2, 3, 10, 11, 12, 13, 1, 2, 1, 2, 3};
  DbLine* dbline3      = DbLine::createFromSamplesById(
         (int)x1.size(), ELoadBy::COLUMN, tab, lineIds, ranksPerId,
         {"x1", "x2", "z1", "z2"}, {"x1", "x2", "z1", "z2"});
  if (dbline3 == nullptr) return 1;
  mestitle(1, "DbLine created using alternative solution");
  dbline3->display();

  // Create the corresponding Header file
  // This new file contains one variable which gives the number of samples per
  // Line
  Db* db = dbline->createStatToHeader();

  DbStringFormat* dbfmt =
    DbStringFormat::createFromFlags(true, true, false, false, true);
  db->display(dbfmt);

  delete dbline;
  delete dbline2;
  delete dbline3;
  delete dbfmt;
  delete db;

  return 0;
}

