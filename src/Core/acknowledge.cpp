/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://gstlearn.org                                              */
/* License: BSD 3-clause                                                      */
/*                                                                            */
/******************************************************************************/
#include "geoslib_old_f.h"
#include "geoslib_f.h"
#include "Basic/File.hpp"
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
void inquire_gstlearn(char **release,
                                      char **date)
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

  message("gstlearn Library (Version:%s - Date:%s)",
          GSTLEARN_VERSION,GSTLEARN_DATE);

  // Print the list of authors

  message("\n");
  message("Authors:\n");
  message("Didier RENARD    (didier.renard@mines-paristech.fr)\n");
  message("Nicolas BEZ      (nicolas.bez@ird.fr)\n");
  message("Nicolas DESASSIS (nicolas.desassis@mines-paristech.fr)\n");
  message("Helene BEUCHER   (helene.beucher@mines-paristech.fr)\n");
  message("Fabien ORS       (fabien.ors@mines-paristech.fr)\n");
  message("Xavier FREULON   (xavier.freulon@mines-paristech.fr)\n");
  return;
}

