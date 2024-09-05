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
#include "geoslib_old_f.h"
#include "Basic/String.hpp"
#include "version.h"

#include <string.h>

/****************************************************************************/
/*!
**  Inquiry about the Version release and Date
**
** \param[out]  release : Array of characters with Release Name
** \param[out]  date    : Array of characters with Data
**
** \remarks The output arrays should be freed by the calling program
**
****************************************************************************/
void inquire_gstlearn(char **release, char **date)
{
  char *buffer;

  int size  = static_cast<int> (strlen(GSTLEARN_VERSION));
  buffer = (char *) mem_alloc(sizeof(char) * (size+1),1);
  (void) gslStrcpy(buffer,GSTLEARN_VERSION);
  buffer[size] = '\0';
  *release = buffer;

  size  = static_cast<int> (strlen(GSTLEARN_VERSION));
  buffer = (char *) mem_alloc(sizeof(char) * (size+1),1);
  (void) gslStrcpy(buffer,GSTLEARN_VERSION);
  buffer[size] = '\0';
  *date = buffer;
}
  
/****************************************************************************/
/*!
 *  Acknowledgment of the authors for gstlearn Library
 *
 ****************************************************************************/
void acknowledge_gstlearn(void)

{
  // Print the header 

  message("gstlearn Library (Version: %s - Date: %s - Commit: %s)",
          GSTLEARN_VERSION, GSTLEARN_DATE, GSTLEARN_COMMIT);

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

