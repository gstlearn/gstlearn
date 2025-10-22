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
#include "Db/Db.hpp"
#include "geoslib_define.h"

using namespace gstlrn;

int main(int argc, char* argv[])
{
  DECLARE_UNUSED(argc);
  DECLARE_UNUSED(argv);

  auto* db = Db::createFillRandom();

  db->addColumnsByConstant(1, 2., "MyFirst");
  db->addColumnsByConstant(1, 2., "MyFirst");
  db->addColumnsByConstant(1, 2., "MySecond");
  db->addColumnsByConstant(1, 2., "MyFirst");
  db->addColumnsByConstant(1, 2., "MySecond");
  db->addColumnsByConstant(1, 2., "MySecond");
  db->addColumnsByConstant(1, 2., "MyFirst");

  db->display();
  return 0;
}
