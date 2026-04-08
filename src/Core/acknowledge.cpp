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
#include "Core/Acknowledge.hpp"
#include "Basic/Message.hpp"
#include "version.h"

namespace gstlrn
{
  /****************************************************************************/
  /*!
   *  Acknowledgment of the authors for gstlearn Library
   *
   ****************************************************************************/
  void acknowledge_gstlearn(void)
  {
    // Print the header

    message(
      "gstlearn Library (Version: %s - Date: %s - Commit: %s)",
      GSTLEARN_FULL_VERSION, GSTLEARN_DATE, GSTLEARN_COMMIT);

    // Print the list of authors

    message("\n");
    message("Authors:\n");
    message("Didier RENARD    (didier.renard@minesparis.psl.eu)\n");
    message("Fabien ORS       (fabien.ors@minesparis.psl.eu)\n");
    message("Nicolas DESASSIS (nicolas.desassis@minesparis.psl.eu)\n");
    message("Pierre GUILLOU   (pierre.guillou@minesparis.psl.eu)\n");
    message("Xavier FREULON   (xavier.freulon@minesparis.psl.eu)\n");
    message("Mike PEREIRA     (mike.pereira@minesparis.psl.eu)\n");
  }

} // namespace gstlrn
