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
#include "Basic/File.hpp"
#include "geoslib_e.h"
#include "geoslib_old_f.h"
#include "version.h"

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
GEOSLIB_API void inquire_Geoslib(char **release,
                                 char **date)
{
  char *buffer;

  int size  = static_cast<int> (strlen(GEOSLIB_RELEASE));
  buffer = (char *) mem_alloc(sizeof(char) * (size+1),1);
  (void) gslStrcpy(buffer,size,GEOSLIB_RELEASE);
  buffer[size] = '\0';
  *release = buffer;

  size  = static_cast<int> (strlen(GEOSLIB_DATE));
  buffer = (char *) mem_alloc(sizeof(char) * (size+1),1);
  (void) gslStrcpy(buffer,size,GEOSLIB_DATE);
  buffer[size] = '\0';
  *date = buffer;
}
  
/****************************************************************************/
/*!
 *  Acknowledgment of the authors for Geoslib Library
 *
 ****************************************************************************/
GEOSLIB_API void acknowledge_Geoslib(void)

{
  // Print the header 

  message("Geoslib Library (Version:%s - Date:%s)",
          GEOSLIB_RELEASE,GEOSLIB_DATE);

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

